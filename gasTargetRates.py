#!/usr/bin/env python3
"""
Solar-neutrino nu-e interaction rates for gas targets listed in
DetectorNumbers/GasDensities.csv.

This script reuses the flavor-separated solar-neutrino fluxes built by
plotFluxes.py and the single-electron scattering helpers defined in
scatteringPlots.py. The only new ingredient is the conversion from gas
density to target-electron density.

Outputs for each gas entry:
  - one subfolder per gas target
  - dR/dE_nu tables
  - dR/dT_e tables
  - dR/dcos(theta_e) tables
  - dR/dE_nu plots
  - dR/dT_e plots
  - dR/dcos(theta_e) plots
  - an overall summary CSV with integrated rates
  - a summary plot comparing total rates across gases

All rates are scaled to the requested fiducial volume.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from detector_model import (
    DEFAULT_DIFFUSION_CSV,
    DEFAULT_FLUX_CSV,
    DEFAULT_GAS_CSV,
    DEFAULT_GEOMETRY_CONFIG,
    DEFAULT_RANGE_CSV,
    PLOT_DPI,
    REPO_ROOT,
    ProgressBar,
    electron_density_cm3,
    gas_entry_label,
    gas_entry_slug,
    read_detector_geometry_config,
    read_diffusion_summary_table,
    read_gas_density_table,
    read_recoil_window_table,
    resolve_recoil_window_keV,
)
import scatteringPlots as sp
from cevns_pipeline import (
    add_cevns_cli_args,
    compute_cevns_spectra_for_gas,
    load_cevns_config,
    save_cevns_tables,
    solar_active_flux,
    write_cevns_spectrum_plots,
    write_cevns_summary,
)


SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0

DEFAULT_OUTDIR = REPO_ROOT / "solar_nu_gas_target_rates"
DEFAULT_RESPONSE_CONFIG = REPO_ROOT / "reco_response_config.json"

TRAPEZOID = getattr(np, "trapezoid", np.trapz)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute solar-neutrino nu-e differential interaction rates for the "
            "gas targets listed in DetectorNumbers/GasDensities.csv."
        )
    )
    parser.add_argument(
        "--flux-csv",
        default=str(DEFAULT_FLUX_CSV),
        help="Flavor-separated solar-neutrino flux CSV produced by plotFluxes.py",
    )
    parser.add_argument(
        "--gas-csv",
        default=str(DEFAULT_GAS_CSV),
        help="Gas density table",
    )
    parser.add_argument(
        "--range-csv",
        default=str(DEFAULT_RANGE_CSV),
        help="Recoil-energy thresholds derived from electron ranges",
    )
    parser.add_argument(
        "--diffusion-csv",
        default=str(DEFAULT_DIFFUSION_CSV),
        help=(
            "Diffusion summary CSV used to raise the low-energy threshold when "
            "DL_2sigma exceeds 1 mm for electric fields <= 2 kV/cm"
        ),
    )
    parser.add_argument(
        "--geometry-config",
        default=str(DEFAULT_GEOMETRY_CONFIG),
        help="Detector geometry JSON config",
    )
    parser.add_argument(
        "--response-config",
        default=str(DEFAULT_RESPONSE_CONFIG),
        help="Response JSON config containing the optional CEvNS block",
    )
    parser.add_argument(
        "--outdir",
        default=str(DEFAULT_OUTDIR),
        help="Output directory",
    )
    parser.add_argument(
        "--volume-cm3",
        type=float,
        default=None,
        help="Optional fiducial detector volume override in cm^3. Default: use geometry config.",
    )
    parser.add_argument(
        "--t-grid-points",
        type=int,
        default=1600,
        help="Number of recoil-energy grid points",
    )
    parser.add_argument(
        "--costh-grid-points",
        type=int,
        default=1200,
        help="Number of recoil-direction grid points",
    )
    add_cevns_cli_args(parser)
    return parser.parse_args()


def sigma_total_sm_from_dT(E_MeV: np.ndarray, channel: str) -> np.ndarray:
    E = np.asarray(E_MeV, dtype=float)
    gL, gR = sp.couplings(channel)

    prefactor_gev = 2.0 * (sp.GF_GEV ** 2) * (sp.M_E_MEV * 1e-3) / np.pi
    prefactor = prefactor_gev * sp.GEV2_TO_CM2 / 1000.0

    with np.errstate(divide="ignore", invalid="ignore"):
        Tmax_vals = 2.0 * E**2 / (sp.M_E_MEV + 2.0 * E)
        term = (
            gL**2 * Tmax_vals
            + gR**2 * (
                Tmax_vals
                - (Tmax_vals**2) / E
                + (Tmax_vals**3) / (3.0 * E**2)
            )
            - gL * gR * sp.M_E_MEV * (Tmax_vals**2) / (2.0 * E**2)
        )

    out = prefactor * term
    out[~np.isfinite(out)] = 0.0
    out[E <= 0.0] = 0.0
    return out


def sigma_sm_integral_range(
    E_MeV: np.ndarray,
    T_low_MeV: float,
    T_high_MeV: float,
    channel: str,
) -> np.ndarray:
    E = np.asarray(E_MeV, dtype=float)
    gL, gR = sp.couplings(channel)

    prefactor_gev = 2.0 * (sp.GF_GEV ** 2) * (sp.M_E_MEV * 1e-3) / np.pi
    prefactor = prefactor_gev * sp.GEV2_TO_CM2 / 1000.0

    physical_tmax = 2.0 * E**2 / (sp.M_E_MEV + 2.0 * E)
    t1 = np.full_like(E, T_low_MeV, dtype=float)
    t2 = np.minimum(np.full_like(E, T_high_MeV, dtype=float), physical_tmax)

    valid = (E > 0.0) & (t2 > t1)
    out = np.zeros_like(E, dtype=float)
    if not np.any(valid):
        return out

    ev = E[valid]
    a = t1[valid]
    b = t2[valid]

    term = (
        (gL**2 + gR**2) * (b - a)
        - gR**2 * (b**2 - a**2) / ev
        + gR**2 * (b**3 - a**3) / (3.0 * ev**2)
        - gL * gR * sp.M_E_MEV * (b**2 - a**2) / (2.0 * ev**2)
    )
    out[valid] = prefactor * term
    out[~np.isfinite(out)] = 0.0
    out[out < 0.0] = 0.0
    return out


def sigma_accepted_rescaled(
    E_MeV: np.ndarray,
    channel: str,
    T_low_MeV: float,
    T_high_MeV: float,
) -> np.ndarray:
    sigma_sm_total = sigma_total_sm_from_dT(E_MeV, channel=channel)
    sigma_sm_window = sigma_sm_integral_range(
        E_MeV,
        T_low_MeV=T_low_MeV,
        T_high_MeV=T_high_MeV,
        channel=channel,
    )
    sigma_target = sp.sigma_total_approx(E_MeV, channel=channel)

    out = np.zeros_like(np.asarray(E_MeV, dtype=float))
    valid = sigma_sm_total > 0.0
    out[valid] = sigma_target[valid] * (sigma_sm_window[valid] / sigma_sm_total[valid])
    out[~np.isfinite(out)] = 0.0
    out[out < 0.0] = 0.0
    return out


def dsigma_dT_rescaled(E_MeV: float, T_MeV: np.ndarray, channel: str) -> np.ndarray:
    shape = sp.dsigma_dT_sm(E_MeV, T_MeV, channel=channel)
    sigma_sm = sigma_total_sm_from_dT(np.array([E_MeV]), channel=channel)[0]
    sigma_target = sp.sigma_total_approx(np.array([E_MeV]), channel=channel)[0]
    if sigma_sm <= 0.0:
        return np.zeros_like(T_MeV, dtype=float)
    return shape * (sigma_target / sigma_sm)


def convolved_recoil_energy_rate(
    E_MeV: np.ndarray,
    flux: np.ndarray,
    T_MeV: np.ndarray,
    channel: str,
) -> np.ndarray:
    dE = sp.bin_widths_from_centers(E_MeV)
    rate = np.zeros_like(T_MeV, dtype=float)
    for E_i, flux_i, dE_i in zip(E_MeV, flux, dE):
        rate += flux_i * dsigma_dT_rescaled(E_i, T_MeV, channel=channel) * dE_i
    return rate


def apply_recoil_window_to_rate(
    T_MeV: np.ndarray,
    rate: np.ndarray,
    T_low_MeV: float,
    T_high_MeV: float,
) -> np.ndarray:
    out = np.asarray(rate, dtype=float).copy()
    mask = (T_MeV >= T_low_MeV) & (T_MeV <= T_high_MeV)
    out[~mask] = 0.0
    return out


def differential_rate_per_neutrino_energy_window(
    E_MeV: np.ndarray,
    flux: np.ndarray,
    channel: str,
    T_low_MeV: float,
    T_high_MeV: float,
) -> np.ndarray:
    sigma = sigma_accepted_rescaled(
        E_MeV,
        channel=channel,
        T_low_MeV=T_low_MeV,
        T_high_MeV=T_high_MeV,
    )
    return np.asarray(flux, dtype=float) * sigma


def build_recoil_energy_grid(E_MeV: np.ndarray, n_points: int) -> np.ndarray:
    Tmax_global = sp.Tmax(float(np.max(E_MeV)))
    if n_points < 2:
        raise ValueError("Need at least 2 recoil-energy grid points.")
    positive = np.logspace(-6, np.log10(Tmax_global), n_points - 1)
    return np.concatenate([[0.0], positive])


def convolved_angular_rate_window(
    E_MeV: np.ndarray,
    flux: np.ndarray,
    costh: np.ndarray,
    channel: str,
    T_low_MeV: float,
    T_high_MeV: float,
) -> np.ndarray:
    dE = sp.bin_widths_from_centers(E_MeV)
    rate = np.zeros_like(costh, dtype=float)
    E_min_for_window = float(sp.Emin_from_T(np.array([T_low_MeV]))[0]) if T_low_MeV > 0.0 else 0.0

    for E_i, flux_i, dE_i in zip(E_MeV, flux, dE):
        if flux_i <= 0.0 or E_i < E_min_for_window:
            continue
        dsdc = sp.dsigma_dcosth_rescaled(E_i, costh, channel=channel)
        recoil_t = sp.T_from_costh(E_i, costh)
        dsdc[(recoil_t < T_low_MeV) | (recoil_t > T_high_MeV)] = 0.0
        rate += flux_i * dsdc * dE_i
    return rate


def prepare_base_flux_inputs(
    flux_df: pd.DataFrame,
    t_grid_points: int,
    costh_grid_points: int,
) -> dict[str, np.ndarray]:
    E = flux_df["E_MeV"].to_numpy()
    flux_nue = flux_df["nu_e_dphi_dE"].to_numpy()
    flux_numu = flux_df["nu_mu_dphi_dE"].to_numpy()
    flux_nutau = flux_df["nu_tau_dphi_dE"].to_numpy()

    dE = sp.bin_widths_from_centers(E)

    T_grid = build_recoil_energy_grid(E, n_points=t_grid_points)
    costh_grid = np.linspace(0.0, 1.0, costh_grid_points)

    return {
        "E_MeV": E,
        "dE_MeV": dE,
        "T_MeV": T_grid,
        "costh": costh_grid,
        "theta_deg": np.degrees(np.arccos(np.clip(costh_grid, 0.0, 1.0))),
        "flux_nue": flux_nue,
        "flux_numu": flux_numu,
        "flux_nutau": flux_nutau,
        "dR_dT_nue_full": convolved_recoil_energy_rate(E, flux_nue, T_grid, channel="nu_e"),
        "dR_dT_numu_full": convolved_recoil_energy_rate(E, flux_numu, T_grid, channel="nu_x"),
        "dR_dT_nutau_full": convolved_recoil_energy_rate(E, flux_nutau, T_grid, channel="nu_x"),
    }


def compute_per_electron_spectra(
    base_inputs: dict[str, np.ndarray],
    T_low_MeV: float,
    T_high_MeV: float,
) -> dict[str, np.ndarray]:
    E = base_inputs["E_MeV"]
    dE = base_inputs["dE_MeV"]
    T_grid = base_inputs["T_MeV"]
    costh_grid = base_inputs["costh"]

    dR_dEnu_nue = differential_rate_per_neutrino_energy_window(
        E,
        base_inputs["flux_nue"],
        channel="nu_e",
        T_low_MeV=T_low_MeV,
        T_high_MeV=T_high_MeV,
    )
    dR_dEnu_numu = differential_rate_per_neutrino_energy_window(
        E,
        base_inputs["flux_numu"],
        channel="nu_x",
        T_low_MeV=T_low_MeV,
        T_high_MeV=T_high_MeV,
    )
    dR_dEnu_nutau = differential_rate_per_neutrino_energy_window(
        E,
        base_inputs["flux_nutau"],
        channel="nu_x",
        T_low_MeV=T_low_MeV,
        T_high_MeV=T_high_MeV,
    )
    dR_dEnu_total = dR_dEnu_nue + dR_dEnu_numu + dR_dEnu_nutau

    rate_bin_nue = dR_dEnu_nue * dE
    rate_bin_numu = dR_dEnu_numu * dE
    rate_bin_nutau = dR_dEnu_nutau * dE
    rate_bin_total = dR_dEnu_total * dE

    dR_dT_nue = apply_recoil_window_to_rate(
        T_grid,
        base_inputs["dR_dT_nue_full"],
        T_low_MeV=T_low_MeV,
        T_high_MeV=T_high_MeV,
    )
    dR_dT_numu = apply_recoil_window_to_rate(
        T_grid,
        base_inputs["dR_dT_numu_full"],
        T_low_MeV=T_low_MeV,
        T_high_MeV=T_high_MeV,
    )
    dR_dT_nutau = apply_recoil_window_to_rate(
        T_grid,
        base_inputs["dR_dT_nutau_full"],
        T_low_MeV=T_low_MeV,
        T_high_MeV=T_high_MeV,
    )
    dR_dT_total = dR_dT_nue + dR_dT_numu + dR_dT_nutau

    dR_dcosth_nue = convolved_angular_rate_window(
        E,
        base_inputs["flux_nue"],
        costh_grid,
        channel="nu_e",
        T_low_MeV=T_low_MeV,
        T_high_MeV=T_high_MeV,
    )
    dR_dcosth_numu = convolved_angular_rate_window(
        E,
        base_inputs["flux_numu"],
        costh_grid,
        channel="nu_x",
        T_low_MeV=T_low_MeV,
        T_high_MeV=T_high_MeV,
    )
    dR_dcosth_nutau = convolved_angular_rate_window(
        E,
        base_inputs["flux_nutau"],
        costh_grid,
        channel="nu_x",
        T_low_MeV=T_low_MeV,
        T_high_MeV=T_high_MeV,
    )
    dR_dcosth_total = dR_dcosth_nue + dR_dcosth_numu + dR_dcosth_nutau

    return {
        "E_MeV": E,
        "dE_MeV": dE,
        "T_MeV": T_grid,
        "costh": costh_grid,
        "theta_deg": base_inputs["theta_deg"],
        "dR_dEnu_nue": dR_dEnu_nue,
        "dR_dEnu_numu": dR_dEnu_numu,
        "dR_dEnu_nutau": dR_dEnu_nutau,
        "dR_dEnu_total": dR_dEnu_total,
        "rate_bin_nue": rate_bin_nue,
        "rate_bin_numu": rate_bin_numu,
        "rate_bin_nutau": rate_bin_nutau,
        "rate_bin_total": rate_bin_total,
        "dR_dT_nue": dR_dT_nue,
        "dR_dT_numu": dR_dT_numu,
        "dR_dT_nutau": dR_dT_nutau,
        "dR_dT_total": dR_dT_total,
        "dR_dcosth_nue": dR_dcosth_nue,
        "dR_dcosth_numu": dR_dcosth_numu,
        "dR_dcosth_nutau": dR_dcosth_nutau,
        "dR_dcosth_total": dR_dcosth_total,
    }


def scale_spectra(per_electron: dict[str, np.ndarray], target_electrons: float) -> dict[str, np.ndarray]:
    scaled = {}
    for key, value in per_electron.items():
        if key in {"E_MeV", "dE_MeV", "T_MeV", "costh", "theta_deg"}:
            scaled[key] = value
        else:
            scaled[key] = value * target_electrons
    return scaled


def write_rate_tables(
    gas_outdir: Path,
    scaled: dict[str, np.ndarray],
) -> None:
    dEnu_df = pd.DataFrame(
        {
            "E_MeV": scaled["E_MeV"],
            "dE_MeV": scaled["dE_MeV"],
            "dR_dEnu_nu_e_s^-1_MeV^-1": scaled["dR_dEnu_nue"],
            "dR_dEnu_nu_mu_s^-1_MeV^-1": scaled["dR_dEnu_numu"],
            "dR_dEnu_nu_tau_s^-1_MeV^-1": scaled["dR_dEnu_nutau"],
            "dR_dEnu_total_s^-1_MeV^-1": scaled["dR_dEnu_total"],
            "rate_bin_total_s^-1": scaled["rate_bin_total"],
        }
    )
    dEnu_df.to_csv(gas_outdir / "dR_dEnu.csv", index=False)

    dTe_df = pd.DataFrame(
        {
            "T_MeV": scaled["T_MeV"],
            "T_keV": 1.0e3 * scaled["T_MeV"],
            "dR_dTe_nu_e_s^-1_MeV^-1": scaled["dR_dT_nue"],
            "dR_dTe_nu_mu_s^-1_MeV^-1": scaled["dR_dT_numu"],
            "dR_dTe_nu_tau_s^-1_MeV^-1": scaled["dR_dT_nutau"],
            "dR_dTe_total_s^-1_MeV^-1": scaled["dR_dT_total"],
            "dR_dTe_nu_e_s^-1_keV^-1": scaled["dR_dT_nue"] / 1.0e3,
            "dR_dTe_nu_mu_s^-1_keV^-1": scaled["dR_dT_numu"] / 1.0e3,
            "dR_dTe_nu_tau_s^-1_keV^-1": scaled["dR_dT_nutau"] / 1.0e3,
            "dR_dTe_total_s^-1_keV^-1": scaled["dR_dT_total"] / 1.0e3,
        }
    )
    dTe_df.to_csv(gas_outdir / "dR_dTe.csv", index=False)

    direction_df = pd.DataFrame(
        {
            "cos_theta_e": scaled["costh"],
            "theta_deg": scaled["theta_deg"],
            "dR_dcosth_nu_e_s^-1": scaled["dR_dcosth_nue"],
            "dR_dcosth_nu_mu_s^-1": scaled["dR_dcosth_numu"],
            "dR_dcosth_nu_tau_s^-1": scaled["dR_dcosth_nutau"],
            "dR_dcosth_total_s^-1": scaled["dR_dcosth_total"],
        }
    )
    direction_df.to_csv(gas_outdir / "dR_dcosth.csv", index=False)


def plot_positive_series(ax, x: np.ndarray, y: np.ndarray, label: str, **kwargs) -> None:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0.0) & (y > 0.0)
    if np.any(mask):
        ax.plot(x[mask], y[mask], label=label, **kwargs)


def write_rate_plots(
    gas_outdir: Path,
    gas_label: str,
    T_low_keV: float,
    T_high_keV: float,
    scaled: dict[str, np.ndarray],
) -> None:
    E = scaled["E_MeV"]
    T_keV = 1.0e3 * scaled["T_MeV"]
    costh = scaled["costh"]
    window_text = f"Accepted recoil window: {T_low_keV:.2f} to {T_high_keV:.2f} keV"

    # dR / dE_nu
    fig, ax = plt.subplots(figsize=(8, 5.5))
    plot_positive_series(ax, E, scaled["dR_dEnu_nue"], label=r"$\nu_e$", linewidth=2.0)
    plot_positive_series(ax, E, scaled["dR_dEnu_numu"], label=r"$\nu_\mu$", linewidth=2.0)
    plot_positive_series(ax, E, scaled["dR_dEnu_nutau"], label=r"$\nu_\tau$", linewidth=2.0)
    plot_positive_series(
        ax,
        E,
        scaled["dR_dEnu_total"],
        label="total",
        linewidth=2.5,
        linestyle="--",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Neutrino energy $E_\nu$ [MeV]")
    ax.set_ylabel(r"$dR/dE_\nu$ [s$^{-1}$ MeV$^{-1}$]")
    ax.set_title(f"Gas-target differential rate vs neutrino energy\n{gas_label}\n{window_text}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(gas_outdir / "dR_dEnu.png", dpi=PLOT_DPI)
    plt.close(fig)

    # dR / dT_e
    fig, ax = plt.subplots(figsize=(8, 5.5))
    plot_positive_series(
        ax,
        T_keV,
        scaled["dR_dT_nue"] / 1.0e3,
        label=r"$\nu_e$",
        linewidth=2.0,
    )
    plot_positive_series(
        ax,
        T_keV,
        scaled["dR_dT_numu"] / 1.0e3,
        label=r"$\nu_\mu$",
        linewidth=2.0,
    )
    plot_positive_series(
        ax,
        T_keV,
        scaled["dR_dT_nutau"] / 1.0e3,
        label=r"$\nu_\tau$",
        linewidth=2.0,
    )
    plot_positive_series(
        ax,
        T_keV,
        scaled["dR_dT_total"] / 1.0e3,
        label="total",
        linewidth=2.5,
        linestyle="--",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Recoil kinetic energy $T_e$ [keV]")
    ax.set_ylabel(r"$dR/dT_e$ [s$^{-1}$ keV$^{-1}$]")
    ax.axvline(T_low_keV, color="black", linestyle=":", linewidth=1.2)
    ax.axvline(T_high_keV, color="black", linestyle=":", linewidth=1.2)
    ax.set_title(f"Gas-target differential rate vs recoil energy\n{gas_label}\n{window_text}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(gas_outdir / "dR_dTe.png", dpi=PLOT_DPI)
    plt.close(fig)

    # dR / dcos(theta_e)
    fig, ax = plt.subplots(figsize=(8, 5.5))
    y_series = [
        (scaled["dR_dcosth_nue"], r"$\nu_e$", {"linewidth": 2.0}),
        (scaled["dR_dcosth_numu"], r"$\nu_\mu$", {"linewidth": 2.0}),
        (scaled["dR_dcosth_nutau"], r"$\nu_\tau$", {"linewidth": 2.0}),
        (
            scaled["dR_dcosth_total"],
            "total",
            {"linewidth": 2.5, "linestyle": "--"},
        ),
    ]
    for y, label, style in y_series:
        mask = np.isfinite(costh) & np.isfinite(y) & (y > 0.0)
        if np.any(mask):
            ax.plot(costh[mask], y[mask], label=label, **style)
    ax.set_yscale("log")
    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel(r"$\cos\theta_e$")
    ax.set_ylabel(r"$dR/d\cos\theta_e$ [s$^{-1}$]")
    ax.set_title(f"Gas-target differential rate vs recoil direction\n{gas_label}\n{window_text}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(gas_outdir / "dR_dcosth.png", dpi=PLOT_DPI)
    plt.close(fig)


def write_summary_plot(outdir: Path, summary_df: pd.DataFrame) -> None:
    ordered = summary_df.sort_values("rate_total_year^-1", ascending=True).copy()

    labels = []
    for gas_name, pressure_atm, tmin, tmax, e_nu_min, e_nu_max in zip(
        ordered["gas"],
        ordered["pressure_atm"],
        ordered["recoil_window_min_keV"],
        ordered["recoil_window_max_keV"],
        ordered["E_nu_min_MeV"],
        ordered["E_nu_max_MeV"],
    ):
        gas_lbl = gas_entry_label(gas_name, pressure_atm)
        labels.append(
            f"{gas_lbl}\n$T_e$: [{tmin:.1e}, {tmax:.1e}] keV"
            f"\n$E_\\nu$: [{e_nu_min*1000:.1e}, {e_nu_max*1000:.1e}] keV"
        )

    fig, ax = plt.subplots(figsize=(11, 6))
    ax.barh(labels, ordered["rate_total_year^-1"], color="#3b82f6")
    ax.set_xscale("log")
    ax.set_xlabel(r"Total interaction rate [year$^{-1}$]")
    ax.set_title("Solar neutrino interaction rate by gas target\n(with investigated electron-energy window)")
    ax.grid(True, axis="x", which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / "gas_target_total_rate_comparison.png", dpi=PLOT_DPI)
    plt.close(fig)


def main() -> int:
    args = parse_args()

    flux_csv = Path(args.flux_csv)
    gas_csv = Path(args.gas_csv)
    range_csv = Path(args.range_csv)
    diffusion_csv = Path(args.diffusion_csv)
    geometry_config = Path(args.geometry_config)
    response_config = Path(args.response_config)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not flux_csv.exists():
        raise FileNotFoundError(f"Flux CSV not found: {flux_csv}")
    if not gas_csv.exists():
        raise FileNotFoundError(f"Gas-density CSV not found: {gas_csv}")
    if not range_csv.exists():
        raise FileNotFoundError(f"Recoil-window CSV not found: {range_csv}")
    if not diffusion_csv.exists():
        raise FileNotFoundError(f"Diffusion summary CSV not found: {diffusion_csv}")
    if not geometry_config.exists():
        raise FileNotFoundError(f"Detector geometry config not found: {geometry_config}")

    flux_df = sp.read_flux_csv(flux_csv)
    gas_df = read_gas_density_table(gas_csv)
    recoil_window_df = read_recoil_window_table(range_csv)
    diffusion_df = read_diffusion_summary_table(diffusion_csv)
    geometry = read_detector_geometry_config(geometry_config)
    volume_cm3 = float(args.volume_cm3) if args.volume_cm3 is not None else geometry.volume_cm3
    volume_source = "cli_override" if args.volume_cm3 is not None else "geometry_config"
    cevns_config = load_cevns_config(response_config, args, default_enabled=False)
    cevns_flux = solar_active_flux(flux_df) if cevns_config.enabled else None
    cevns_outdir = outdir / "cevns"
    cevns_summary_rows = []

    base_inputs = prepare_base_flux_inputs(
        flux_df=flux_df,
        t_grid_points=args.t_grid_points,
        costh_grid_points=args.costh_grid_points,
    )

    summary_rows = []
    progress = ProgressBar(len(gas_df), prefix="Gas targets")
    progress.update(0, "starting")
    for gas_index, (_, row) in enumerate(gas_df.iterrows(), start=1):
        gas_name = row["Gas"]
        density_g_cm3 = float(row["density @ 293 K (g/cm3)"])
        pressure_atm = float(row["Pressure (atm)"])
        gas_label = gas_entry_label(gas_name, pressure_atm)
        progress.update(gas_index - 1, gas_label)
        recoil_window = resolve_recoil_window_keV(
            gas_name,
            density_g_cm3,
            pressure_atm,
            recoil_window_df=recoil_window_df,
            diffusion_df=diffusion_df,
        )
        T_low_keV = recoil_window.low_keV
        T_high_keV = recoil_window.high_keV
        T_low_MeV = T_low_keV / 1.0e3
        T_high_MeV = T_high_keV / 1.0e3

        electron_density, molar_mass, electrons_per_mixture = electron_density_cm3(
            gas_name,
            density_g_cm3,
        )
        target_electrons = electron_density * volume_cm3
        gas_slug = gas_entry_slug(gas_name, pressure_atm)
        gas_outdir = outdir / gas_slug
        gas_outdir.mkdir(parents=True, exist_ok=True)
        per_electron = compute_per_electron_spectra(
            base_inputs,
            T_low_MeV=T_low_MeV,
            T_high_MeV=T_high_MeV,
        )
        scaled = scale_spectra(per_electron, target_electrons=target_electrons)
        write_rate_tables(gas_outdir, scaled=scaled)
        write_rate_plots(
            gas_outdir,
            gas_label=gas_label,
            T_low_keV=T_low_keV,
            T_high_keV=T_high_keV,
            scaled=scaled,
        )

        if cevns_config.enabled:
            recoil_df, enu_df, cevns_summary_df = compute_cevns_spectra_for_gas(
                gas_name=gas_name,
                gas_label=gas_label,
                density_g_cm3=density_g_cm3,
                volume_cm3=volume_cm3,
                E_MeV=base_inputs["E_MeV"],
                active_flux_or_fluence=cevns_flux,
                config=cevns_config,
                t_grid_points=args.t_grid_points,
                quantity="rate",
            )
            cevns_gas_outdir = cevns_outdir / gas_slug
            save_cevns_tables(
                cevns_gas_outdir,
                recoil_df=recoil_df,
                enu_df=enu_df,
                reco_df=None,
                skip_save=args.skip_cevns_save,
            )
            write_cevns_spectrum_plots(
                cevns_gas_outdir,
                gas_label=gas_label,
                recoil_df=recoil_df,
                enu_df=enu_df,
                quantity="rate",
                skip_plots=args.skip_cevns_plots,
            )
            cevns_summary_df["output_subdir"] = str(Path("cevns") / gas_slug)
            cevns_summary_rows.append(cevns_summary_df)

        total_rate_from_Enu = float(np.sum(scaled["rate_bin_total"]))
        total_rate_from_T = float(TRAPEZOID(scaled["dR_dT_total"], scaled["T_MeV"]))
        total_rate_from_costh = float(TRAPEZOID(scaled["dR_dcosth_total"], scaled["costh"]))

        E_nu_min_MeV = float(sp.Emin_from_T(np.array([T_low_MeV]))[0]) if T_low_MeV > 0.0 else 0.0
        E_nu_max_MeV = float(sp.Emin_from_T(np.array([T_high_MeV]))[0]) if T_high_MeV > 0.0 else 0.0

        summary_rows.append(
            {
                "gas": gas_name,
                "pressure_atm": pressure_atm,
                "density_g_cm3": density_g_cm3,
                "recoil_window_min_keV": T_low_keV,
                "recoil_window_max_keV": T_high_keV,
                "recoil_window_min_range_mm": recoil_window.low_range_mm,
                "base_recoil_window_min_keV_at_1mm": recoil_window.base_low_keV,
                "threshold_source": recoil_window.threshold_source,
                "diffusion_DL_2sigma_mm": recoil_window.diffusion_dl_2sigma_mm,
                "diffusion_field_V_cm": recoil_window.diffusion_field_v_cm,
                "diffusion_point": recoil_window.diffusion_point,
                "E_nu_min_MeV": E_nu_min_MeV,
                "E_nu_max_MeV": E_nu_max_MeV,
                "detector_name": geometry.name,
                "detector_length_m": geometry.length_m,
                "detector_width_m": geometry.width_m,
                "detector_height_m": geometry.height_m,
                "configured_volume_m3": geometry.volume_m3,
                "volume_source": volume_source,
                "volume_cm3": volume_cm3,
                "volume_m3": volume_cm3 / 1.0e6,
                "molar_mass_g_mol": molar_mass,
                "electrons_per_mixture": electrons_per_mixture,
                "electron_density_cm^-3": electron_density,
                "target_electrons": target_electrons,
                "rate_nu_e_s^-1": float(np.sum(scaled["rate_bin_nue"])),
                "rate_nu_mu_s^-1": float(np.sum(scaled["rate_bin_numu"])),
                "rate_nu_tau_s^-1": float(np.sum(scaled["rate_bin_nutau"])),
                "rate_total_s^-1": total_rate_from_Enu,
                "rate_total_year^-1": total_rate_from_Enu * SECONDS_PER_YEAR,
                "rate_total_from_dTe_s^-1": total_rate_from_T,
                "rate_total_from_dcosth_s^-1": total_rate_from_costh,
                "gas_label": gas_label,
                "output_subdir": gas_slug,
            }
        )
        progress.update(gas_index, gas_label)
    progress.close()

    summary_df = pd.DataFrame(summary_rows)
    summary_path = outdir / "gas_target_rate_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    write_summary_plot(outdir, summary_df)
    cevns_summary_df = write_cevns_summary(
        cevns_outdir,
        cevns_summary_rows,
        quantity="rate",
        skip_save=args.skip_cevns_save,
        skip_plots=args.skip_cevns_plots,
        filename="cevns_rate_summary.csv",
    ) if cevns_config.enabled else pd.DataFrame()

    print(f"Read {len(flux_df)} flux points from: {flux_csv}")
    print(f"Read {len(gas_df)} gas entries from: {gas_csv}")
    print(f"Read {len(diffusion_df)} diffusion rows from: {diffusion_csv}")
    print(
        f"Using detector geometry '{geometry.name}' = "
        f"{geometry.length_m:g} x {geometry.width_m:g} x {geometry.height_m:g} m^3 "
        f"(volume = {volume_cm3 / 1.0e6:g} m^3)"
    )
    print(f"Saved gas-target rate tables and plots in: {outdir.resolve()}")
    print(f"Saved summary: {summary_path.resolve()}")
    if cevns_config.enabled:
        print(
            "CEvNS enabled: "
            f"threshold={cevns_config.nr_threshold_keV:g} keV, "
            f"max={cevns_config.nr_max_keV}, "
            f"form_factor={cevns_config.form_factor}, "
            f"axial_model={cevns_config.axial_model}"
        )
        print(f"CEvNS outputs directory: {cevns_outdir.resolve()}")
        if not cevns_summary_df.empty:
            print(
                "CEvNS total solar rate across all gas-aggregated summary rows: "
                f"{cevns_summary_df['total_rate_s-1'].sum():.6e} s^-1"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
