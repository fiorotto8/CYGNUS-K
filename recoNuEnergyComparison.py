#!/usr/bin/env python3
"""
Estimate reconstructed solar-neutrino energy spectra for:

1. energy-only reconstruction
2. energy + direction reconstruction

The script reuses:
  - the flavor-separated solar-neutrino fluxes
  - the gas tables and recoil-energy acceptance windows
  - the detector geometry

For each gas target, it builds a gas-dependent detector-response model with:

    sigma_E / E = sqrt((res_E_thr * sqrt(E_thr / T_e))^2 + res_E_const^2)
    sigma_theta = sqrt((res_theta_thr * sqrt(E_thr / T_e))^2 + res_theta_const^2)

where E_thr is the recoil threshold energy corresponding to 1 mm range in the
gas, taken from DetectorNumbers/electron_range_energy_table.csv.

Unlike the earlier Monte Carlo study, this script now estimates only the
reconstructed dR/dE_nu spectra. It uses propagated neutrino-energy
uncertainties to smear accepted event rates into coarse reconstructed
neutrino-energy bins for the two reconstruction strategies.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
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
from EnergywithPerformance import neutrino_energy_from_Te_theta, neutrino_energy_min
from gasTargetRates import (
    sigma_accepted_rescaled,
    sigma_total_sm_from_dT,
)
from cevns_pipeline import (
    add_cevns_cli_args,
    build_cevns_lower_bound_reco,
    compute_cevns_spectra_for_gas,
    load_cevns_config,
    save_cevns_tables,
    solar_active_flux,
    write_cevns_reco_plot,
    write_cevns_summary,
)


REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_RESPONSE_CONFIG = REPO_ROOT / "reco_response_config.json"
DEFAULT_OUTDIR = REPO_ROOT / "solar_nu_reco_energy_comparison"
M_E_KEV = 511.0
SQRT_TWO = np.sqrt(2.0)
SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0

OBSOLETE_ROOT_OUTPUTS = [
    "reco_directional_improvement_summary.csv",
    "directional_improvement_summary.png",
]

OBSOLETE_GAS_OUTPUTS = [
    "reco_neutrino_energy_resolution.csv",
    "reco_spectrum_ratio_to_truth.png",
    "neutrino_energy_resolution_vs_true.png",
    "directional_improvement_vs_true.png",
    "reco_neutrino_energy_population.png",
]


@dataclass(frozen=True)
class ReconstructionResponse:
    name: str
    energy_resolution_at_threshold_frac: float
    energy_resolution_constant_frac: float
    angular_resolution_at_threshold_deg: float
    angular_resolution_constant_deg: float
    reco_energy_bins: int

    def sigma_E_over_E(self, T_keV: np.ndarray, threshold_keV: float) -> np.ndarray:
        T = np.asarray(T_keV, dtype=float)
        T = np.clip(T, 1.0e-9, None)
        threshold_term = self.energy_resolution_at_threshold_frac * np.sqrt(threshold_keV / T)
        return np.sqrt(threshold_term**2 + self.energy_resolution_constant_frac**2)

    def sigma_theta_deg(self, T_keV: np.ndarray, threshold_keV: float) -> np.ndarray:
        T = np.asarray(T_keV, dtype=float)
        T = np.clip(T, 1.0e-9, None)
        threshold_term = self.angular_resolution_at_threshold_deg * np.sqrt(threshold_keV / T)
        return np.sqrt(threshold_term**2 + self.angular_resolution_constant_deg**2)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Estimate reconstructed neutrino-energy spectra for energy-only and "
            "energy+direction measurements."
        )
    )
    parser.add_argument(
        "--flux-csv",
        default=str(DEFAULT_FLUX_CSV),
        help="Flavor-separated solar-neutrino flux CSV",
    )
    parser.add_argument(
        "--gas-csv",
        default=str(DEFAULT_GAS_CSV),
        help="Gas density table",
    )
    parser.add_argument(
        "--range-csv",
        default=str(DEFAULT_RANGE_CSV),
        help="Recoil-energy threshold table",
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
        help="Reconstruction response JSON config",
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
        help="Optional detector-volume override in cm^3",
    )
    parser.add_argument(
        "--t-grid-points",
        type=int,
        default=220,
        help="Number of recoil-energy grid points used per true neutrino-energy sample",
    )
    parser.add_argument(
        "--max-true-energy-points",
        type=int,
        default=600,
        help=(
            "Compress the accepted true-energy support before smearing. "
            "Use 0 to disable compression and use the full flux grid."
        ),
    )
    add_cevns_cli_args(parser)
    return parser.parse_args()


def read_response_config(path: Path) -> ReconstructionResponse:
    with open(path, "r", encoding="utf-8") as handle:
        data = json.load(handle)

    required = [
        "energy_resolution_at_threshold_frac",
        "angular_resolution_at_threshold_deg",
        "reco_energy_bins",
    ]
    missing = [key for key in required if key not in data]
    if missing:
        raise ValueError(
            "Missing required keys in reconstruction response config: "
            + ", ".join(missing)
        )

    response = ReconstructionResponse(
        name=str(data.get("name", path.stem)),
        energy_resolution_at_threshold_frac=float(data["energy_resolution_at_threshold_frac"]),
        energy_resolution_constant_frac=float(data.get("energy_resolution_constant_frac", 0.0)),
        angular_resolution_at_threshold_deg=float(data["angular_resolution_at_threshold_deg"]),
        angular_resolution_constant_deg=float(data.get("angular_resolution_constant_deg", 0.0)),
        reco_energy_bins=int(data["reco_energy_bins"]),
    )
    if response.reco_energy_bins < 20:
        raise ValueError("Response config must use at least 20 reconstructed-energy bins.")
    if (
        response.energy_resolution_at_threshold_frac < 0.0
        or response.energy_resolution_constant_frac < 0.0
        or response.angular_resolution_at_threshold_deg < 0.0
        or response.angular_resolution_constant_deg < 0.0
    ):
        raise ValueError("All reconstruction response terms must be non-negative.")
    return response


def dsigma_dT_sm_vectorized(E_MeV: np.ndarray, T_MeV: np.ndarray, channel: str) -> np.ndarray:
    E = np.asarray(E_MeV, dtype=float)
    T = np.asarray(T_MeV, dtype=float)
    gL, gR = sp.couplings(channel)

    prefactor_gev = 2.0 * (sp.GF_GEV**2) * (sp.M_E_MEV * 1e-3) / np.pi
    prefactor = prefactor_gev * sp.GEV2_TO_CM2 / 1000.0

    with np.errstate(divide="ignore", invalid="ignore"):
        bracket = (
            gL**2
            + gR**2 * (1.0 - T / E) ** 2
            - gL * gR * (sp.M_E_MEV * T) / (E**2)
        )

    out = prefactor * bracket
    tmax = 2.0 * E**2 / (sp.M_E_MEV + 2.0 * E)
    out[(T < 0.0) | (T > tmax) | (~np.isfinite(out))] = 0.0
    out[out < 0.0] = 0.0
    return out


def dsigma_dT_rescaled_vectorized(E_MeV: np.ndarray, T_MeV: np.ndarray, channel: str) -> np.ndarray:
    shape = dsigma_dT_sm_vectorized(E_MeV, T_MeV, channel=channel)
    sigma_sm = sigma_total_sm_from_dT(np.asarray(E_MeV, dtype=float), channel=channel)
    sigma_target = sp.sigma_total_approx(np.asarray(E_MeV, dtype=float), channel=channel)
    out = np.zeros_like(shape, dtype=float)
    valid = sigma_sm > 0.0
    out[valid] = shape[valid] * (sigma_target[valid] / sigma_sm[valid])
    out[~np.isfinite(out)] = 0.0
    out[out < 0.0] = 0.0
    return out


def costh_from_E_T(E_MeV: np.ndarray, T_MeV: np.ndarray) -> np.ndarray:
    E = np.asarray(E_MeV, dtype=float)
    T = np.asarray(T_MeV, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        c2 = T * (E + sp.M_E_MEV) ** 2 / (E**2 * (2.0 * sp.M_E_MEV + T))
    c2 = np.clip(c2, 0.0, 1.0)
    out = np.sqrt(c2)
    out[~np.isfinite(out)] = 0.0
    return out


def build_reco_energy_edges(E_min_MeV: float, E_max_MeV: float, n_bins: int) -> np.ndarray:
    return np.logspace(np.log10(E_min_MeV), np.log10(E_max_MeV), n_bins + 1)


def histogram_to_density(values: np.ndarray, bin_edges: np.ndarray) -> np.ndarray:
    widths = np.diff(bin_edges)
    out = np.zeros_like(values, dtype=float)
    valid = widths > 0.0
    out[valid] = values[valid] / widths[valid]
    return out


def build_true_rate_arrays(
    flux_df: pd.DataFrame,
    target_electrons: float,
    T_low_MeV: float,
    T_high_MeV: float,
) -> dict[str, np.ndarray]:
    E = flux_df["E_MeV"].to_numpy(dtype=float)
    dE = sp.bin_widths_from_centers(E)
    flux_nue = flux_df["nu_e_dphi_dE"].to_numpy(dtype=float)
    flux_numu = flux_df["nu_mu_dphi_dE"].to_numpy(dtype=float)
    flux_nutau = flux_df["nu_tau_dphi_dE"].to_numpy(dtype=float)

    dR_dEnu_nue = target_electrons * flux_nue * sigma_accepted_rescaled(
        E, "nu_e", T_low_MeV=T_low_MeV, T_high_MeV=T_high_MeV
    )
    dR_dEnu_numu = target_electrons * flux_numu * sigma_accepted_rescaled(
        E, "nu_x", T_low_MeV=T_low_MeV, T_high_MeV=T_high_MeV
    )
    dR_dEnu_nutau = target_electrons * flux_nutau * sigma_accepted_rescaled(
        E, "nu_x", T_low_MeV=T_low_MeV, T_high_MeV=T_high_MeV
    )

    return {
        "E_MeV": E,
        "dE_MeV": dE,
        "dR_dEnu_nue": dR_dEnu_nue,
        "dR_dEnu_numu": dR_dEnu_numu,
        "dR_dEnu_nutau": dR_dEnu_nutau,
        "rate_bin_nue": dR_dEnu_nue * dE,
        "rate_bin_numu": dR_dEnu_numu * dE,
        "rate_bin_nutau": dR_dEnu_nutau * dE,
    }


def erf_approx(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    sign = np.sign(x)
    ax = np.abs(x)
    t = 1.0 / (1.0 + 0.3275911 * ax)
    poly = (
        (((((1.061405429 * t) - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t
        + 0.254829592
    ) * t
    return sign * (1.0 - poly * np.exp(-(ax**2)))


def normal_cdf(x: np.ndarray) -> np.ndarray:
    return 0.5 * (1.0 + erf_approx(np.asarray(x, dtype=float) / SQRT_TWO))


def gaussian_bin_probabilities(
    means: np.ndarray,
    sigmas: np.ndarray,
    bin_edges: np.ndarray,
) -> np.ndarray:
    mu = np.asarray(means, dtype=float).reshape(-1, 1)
    sigma = np.asarray(sigmas, dtype=float).reshape(-1, 1)
    sigma = np.clip(sigma, 1.0e-12, None)
    edges = np.asarray(bin_edges, dtype=float)
    cdf_hi = normal_cdf((edges[1:].reshape(1, -1) - mu) / sigma)
    cdf_lo = normal_cdf((edges[:-1].reshape(1, -1) - mu) / sigma)
    probs = cdf_hi - cdf_lo
    probs = np.clip(probs, 0.0, 1.0)
    return probs


def compress_true_energy_support(
    E_MeV: np.ndarray,
    accepted_rate_bins: np.ndarray,
    max_points: int,
) -> tuple[np.ndarray, np.ndarray]:
    E = np.asarray(E_MeV, dtype=float)
    rates = np.asarray(accepted_rate_bins, dtype=float)
    valid = np.isfinite(E) & np.isfinite(rates) & (E > 0.0) & (rates > 0.0)
    E = E[valid]
    rates = rates[valid]

    if len(E) == 0:
        return np.array([], dtype=float), np.array([], dtype=float)
    if max_points <= 0 or len(E) <= max_points:
        return E, rates

    edges = np.logspace(np.log10(E.min()), np.log10(E.max()), max_points + 1)
    bin_index = np.clip(np.digitize(E, edges) - 1, 0, max_points - 1)
    rate_sums = np.bincount(bin_index, weights=rates, minlength=max_points)
    energy_weighted = np.bincount(bin_index, weights=rates * E, minlength=max_points)
    keep = rate_sums > 0.0
    return energy_weighted[keep] / rate_sums[keep], rate_sums[keep]


def d_enu_min_dT_keV(T_keV: np.ndarray) -> np.ndarray:
    T = np.asarray(T_keV, dtype=float)
    T = np.clip(T, 1.0e-12, None)
    root = np.sqrt(T**2 + 2.0 * M_E_KEV * T)
    return 0.5 * (1.0 + (T + M_E_KEV) / root)


def d_enu_dir_dT_keV(T_keV: np.ndarray, theta_deg: np.ndarray) -> np.ndarray:
    T = np.asarray(T_keV, dtype=float)
    theta = np.asarray(theta_deg, dtype=float)
    theta = np.clip(theta, 0.0, 89.999999)
    T = np.clip(T, 1.0e-12, None)

    p = np.sqrt(T**2 + 2.0 * M_E_KEV * T)
    cos_theta = np.cos(np.radians(theta))
    denom = p * cos_theta - T
    dp_dT = (T + M_E_KEV) / p
    ddenom_dT = cos_theta * dp_dT - 1.0

    out = np.full_like(T, np.nan, dtype=float)
    valid = denom > 1.0e-12
    out[valid] = M_E_KEV * (denom[valid] - T[valid] * ddenom_dT[valid]) / (denom[valid] ** 2)
    return out


def d_enu_dir_dtheta_rad_keV(T_keV: np.ndarray, theta_deg: np.ndarray) -> np.ndarray:
    T = np.asarray(T_keV, dtype=float)
    theta = np.asarray(theta_deg, dtype=float)
    theta = np.clip(theta, 0.0, 89.999999)
    T = np.clip(T, 1.0e-12, None)

    theta_rad = np.radians(theta)
    p = np.sqrt(T**2 + 2.0 * M_E_KEV * T)
    cos_theta = np.cos(theta_rad)
    sin_theta = np.sin(theta_rad)
    denom = p * cos_theta - T

    out = np.full_like(T, np.nan, dtype=float)
    valid = denom > 1.0e-12
    out[valid] = M_E_KEV * T[valid] * p[valid] * sin_theta[valid] / (denom[valid] ** 2)
    return out


def smear_rates_into_reco_bins(
    bin_edges_MeV: np.ndarray,
    weights_s: np.ndarray,
    means_MeV: np.ndarray,
    sigmas_MeV: np.ndarray,
) -> np.ndarray:
    weights = np.asarray(weights_s, dtype=float)
    means = np.asarray(means_MeV, dtype=float)
    sigmas = np.asarray(sigmas_MeV, dtype=float)

    valid = (
        np.isfinite(weights)
        & np.isfinite(means)
        & np.isfinite(sigmas)
        & (weights > 0.0)
        & (means > 0.0)
        & (sigmas >= 0.0)
    )
    if not np.any(valid):
        return np.zeros(len(bin_edges_MeV) - 1, dtype=float)

    probs = gaussian_bin_probabilities(means[valid], sigmas[valid], bin_edges_MeV)
    row_sums = np.sum(probs, axis=1, keepdims=True)
    overshoot = row_sums > 1.0
    if np.any(overshoot):
        probs[overshoot[:, 0]] /= row_sums[overshoot[:, 0]]
    return np.sum(weights[valid].reshape(-1, 1) * probs, axis=0)


def build_estimated_spectra_for_gas(
    flux_df: pd.DataFrame,
    target_electrons: float,
    T_low_keV: float,
    T_high_keV: float,
    response: ReconstructionResponse,
    t_grid_points: int,
    max_true_energy_points: int,
) -> tuple[pd.DataFrame, dict[str, float]]:
    T_low_MeV = T_low_keV / 1.0e3
    T_high_MeV = T_high_keV / 1.0e3
    truth = build_true_rate_arrays(
        flux_df=flux_df,
        target_electrons=target_electrons,
        T_low_MeV=T_low_MeV,
        T_high_MeV=T_high_MeV,
    )

    E_true = truth["E_MeV"]
    E_true_valid = E_true[E_true > 0.0]
    E_min_truth = max(float(neutrino_energy_min(np.array([T_low_keV]))[0] / 1.0e3), float(E_true_valid.min()))
    E_max_truth = float(E_true.max())
    E_min_reco = max(0.25 * E_min_truth, 1.0e-4)
    E_max_reco = max(3.0 * E_max_truth, E_max_truth + 3.0)
    reco_edges = build_reco_energy_edges(E_min_reco, E_max_reco, response.reco_energy_bins)
    reco_centers = np.sqrt(reco_edges[:-1] * reco_edges[1:])
    reco_widths = np.diff(reco_edges)
    u_grid = np.linspace(0.0, 1.0, t_grid_points)
    du_grid = sp.bin_widths_from_centers(u_grid)

    true_rate_total = truth["rate_bin_nue"] + truth["rate_bin_numu"] + truth["rate_bin_nutau"]
    true_hist = np.histogram(E_true, bins=reco_edges, weights=true_rate_total)[0]

    reco_hist_energy_only = np.zeros(len(reco_centers), dtype=float)
    reco_hist_directional = np.zeros(len(reco_centers), dtype=float)

    flavor_payloads = [
        ("rate_bin_nue", "nu_e"),
        ("rate_bin_numu", "nu_x"),
    ]
    truth["rate_bin_nux"] = truth["rate_bin_numu"] + truth["rate_bin_nutau"]
    flavor_payloads[1] = ("rate_bin_nux", "nu_x")

    for rate_key, channel in flavor_payloads:
        E_support, rate_support = compress_true_energy_support(
            E_true,
            truth[rate_key],
            max_points=max_true_energy_points,
        )

        for E_MeV, accepted_rate_bin in zip(E_support, rate_support):
            if accepted_rate_bin <= 0.0:
                continue

            T_max = 2.0 * E_MeV**2 / (sp.M_E_MEV + 2.0 * E_MeV)
            T_upper = min(T_high_MeV, T_max)
            if T_upper <= T_low_MeV:
                continue

            T_span_MeV = T_upper - T_low_MeV
            T_grid_MeV = T_low_MeV + T_span_MeV * u_grid
            dT_MeV = T_span_MeV * du_grid
            E_grid_MeV = np.full_like(T_grid_MeV, E_MeV)
            dsigma_dT = dsigma_dT_rescaled_vectorized(E_grid_MeV, T_grid_MeV, channel=channel)
            shape_weights = dsigma_dT * dT_MeV
            total_shape = float(np.sum(shape_weights))
            if total_shape <= 0.0:
                continue
            cell_rates = accepted_rate_bin * shape_weights / total_shape

            T_grid_keV = 1.0e3 * T_grid_MeV
            sigma_T_keV = response.sigma_E_over_E(T_grid_keV, threshold_keV=T_low_keV) * T_grid_keV

            mu_energy_only_MeV = neutrino_energy_min(T_grid_keV) / 1.0e3
            sigma_energy_only_MeV = np.abs(d_enu_min_dT_keV(T_grid_keV)) * sigma_T_keV / 1.0e3
            reco_hist_energy_only += smear_rates_into_reco_bins(
                reco_edges,
                cell_rates,
                mu_energy_only_MeV,
                sigma_energy_only_MeV,
            )

            costh_true = costh_from_E_T(E_grid_MeV, T_grid_MeV)
            theta_true_deg = np.degrees(np.arccos(np.clip(costh_true, 0.0, 1.0)))
            mu_directional_MeV = neutrino_energy_from_Te_theta(T_grid_keV, theta_true_deg) / 1.0e3
            sigma_theta_rad = np.deg2rad(
                response.sigma_theta_deg(T_grid_keV, threshold_keV=T_low_keV)
            )
            dE_dT_keV = d_enu_dir_dT_keV(T_grid_keV, theta_true_deg)
            dE_dtheta_keV = d_enu_dir_dtheta_rad_keV(T_grid_keV, theta_true_deg)
            sigma_directional_MeV = (
                np.sqrt((dE_dT_keV * sigma_T_keV) ** 2 + (dE_dtheta_keV * sigma_theta_rad) ** 2)
                / 1.0e3
            )
            reco_hist_directional += smear_rates_into_reco_bins(
                reco_edges,
                cell_rates,
                mu_directional_MeV,
                sigma_directional_MeV,
            )

    reco_df = pd.DataFrame(
        {
            "E_nu_bin_low_MeV": reco_edges[:-1],
            "E_nu_bin_high_MeV": reco_edges[1:],
            "E_nu_bin_center_MeV": reco_centers,
            "bin_width_MeV": reco_widths,
            "rate_true_bin_s^-1": true_hist,
            "rate_reco_energy_only_bin_s^-1": reco_hist_energy_only,
            "rate_reco_energy_direction_bin_s^-1": reco_hist_directional,
            "dR_dEnu_true_s^-1_MeV^-1": histogram_to_density(true_hist, reco_edges),
            "dR_dEnu_reco_energy_only_s^-1_MeV^-1": histogram_to_density(reco_hist_energy_only, reco_edges),
            "dR_dEnu_reco_energy_direction_s^-1_MeV^-1": histogram_to_density(reco_hist_directional, reco_edges),
        }
    )

    metrics = {
        "total_rate_true_s^-1": float(np.sum(true_hist)),
        "total_rate_true_year^-1": float(np.sum(true_hist) * SECONDS_PER_YEAR),
        "total_rate_reco_energy_only_s^-1": float(np.sum(reco_hist_energy_only)),
        "total_rate_reco_energy_only_year^-1": float(np.sum(reco_hist_energy_only) * SECONDS_PER_YEAR),
        "total_rate_reco_energy_direction_s^-1": float(np.sum(reco_hist_directional)),
        "total_rate_reco_energy_direction_year^-1": float(np.sum(reco_hist_directional) * SECONDS_PER_YEAR),
        "reco_energy_min_MeV": float(reco_edges[0]),
        "reco_energy_max_MeV": float(reco_edges[-1]),
    }
    return reco_df, metrics


def plot_positive(ax, x: np.ndarray, y: np.ndarray, label: str, **kwargs) -> None:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0.0) & (y > 0.0)
    if np.any(mask):
        ax.plot(x[mask], y[mask], label=label, **kwargs)


def write_reco_plot(
    gas_outdir: Path,
    gas_label: str,
    response: ReconstructionResponse,
    threshold_keV: float,
    reco_df: pd.DataFrame,
) -> None:
    response_text = (
        f"$\\sigma_E/E = \\sqrt{{({response.energy_resolution_at_threshold_frac:.3g}"
        f"\\sqrt{{E_{{thr}}/T_e}})^2 + {response.energy_resolution_constant_frac:.3g}^2}}$, "
        f"$\\sigma_\\theta = \\sqrt{{({response.angular_resolution_at_threshold_deg:.3g}^\\circ"
        f"\\sqrt{{E_{{thr}}/T_e}})^2 + ({response.angular_resolution_constant_deg:.3g}^\\circ)^2}}$, "
        f"$E_{{thr}}={threshold_keV:.2f}$ keV"
    )

    fig, ax = plt.subplots(figsize=(8.7, 5.8))
    plot_positive(
        ax,
        reco_df["E_nu_bin_center_MeV"],
        reco_df["dR_dEnu_true_s^-1_MeV^-1"],
        "true accepted",
        linewidth=2.2,
    )
    plot_positive(
        ax,
        reco_df["E_nu_bin_center_MeV"],
        reco_df["dR_dEnu_reco_energy_only_s^-1_MeV^-1"],
        "estimated: energy only",
        linewidth=2.0,
    )
    plot_positive(
        ax,
        reco_df["E_nu_bin_center_MeV"],
        reco_df["dR_dEnu_reco_energy_direction_s^-1_MeV^-1"],
        "estimated: energy + direction",
        linewidth=2.0,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Neutrino energy [MeV]")
    ax.set_ylabel(r"$dR/dE_\nu$ [s$^{-1}$ MeV$^{-1}$]")
    ax.set_title(f"Estimated reconstructed neutrino-energy spectra\n{gas_label}\n{response_text}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(gas_outdir / "reco_neutrino_energy_spectra.png", dpi=PLOT_DPI)
    plt.close(fig)


def remove_if_exists(path: Path) -> None:
    if path.exists() and path.is_file():
        path.unlink()


def main() -> int:
    args = parse_args()
    if args.t_grid_points < 20:
        raise ValueError("--t-grid-points must be at least 20.")
    if args.max_true_energy_points < 0:
        raise ValueError("--max-true-energy-points must be non-negative.")

    flux_csv = Path(args.flux_csv)
    gas_csv = Path(args.gas_csv)
    range_csv = Path(args.range_csv)
    diffusion_csv = Path(args.diffusion_csv)
    geometry_config = Path(args.geometry_config)
    response_config = Path(args.response_config)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for filename in OBSOLETE_ROOT_OUTPUTS:
        remove_if_exists(outdir / filename)

    for path, label in [
        (flux_csv, "Flux CSV"),
        (gas_csv, "Gas-density CSV"),
        (range_csv, "Recoil-window CSV"),
        (diffusion_csv, "Diffusion summary CSV"),
        (geometry_config, "Geometry config"),
        (response_config, "Response config"),
    ]:
        if not path.exists():
            raise FileNotFoundError(f"{label} not found: {path}")

    flux_df = sp.read_flux_csv(flux_csv)
    gas_df = read_gas_density_table(gas_csv)
    recoil_window_df = read_recoil_window_table(range_csv)
    diffusion_df = read_diffusion_summary_table(diffusion_csv)
    geometry = read_detector_geometry_config(geometry_config)
    response = read_response_config(response_config)
    volume_cm3 = float(args.volume_cm3) if args.volume_cm3 is not None else geometry.volume_cm3
    volume_source = "cli_override" if args.volume_cm3 is not None else "geometry_config"
    cevns_config = load_cevns_config(response_config, args, default_enabled=False)
    cevns_flux = solar_active_flux(flux_df) if cevns_config.enabled else None
    cevns_outdir = outdir / "cevns"
    cevns_summary_rows = []

    summary_rows = []
    progress = ProgressBar(len(gas_df), prefix="Reconstructed spectra")
    progress.update(0, "starting")
    for gas_index, (_, row) in enumerate(gas_df.iterrows(), start=1):
        gas_name = row["Gas"]
        density_g_cm3 = float(row["density @ 293 K (g/cm3)"])
        pressure_atm = float(row["Pressure (atm)"])
        gas_label = gas_entry_label(gas_name, pressure_atm)
        gas_slug = gas_entry_slug(gas_name, pressure_atm)
        gas_outdir = outdir / gas_slug
        gas_outdir.mkdir(parents=True, exist_ok=True)
        progress.update(gas_index - 1, gas_label)

        for filename in OBSOLETE_GAS_OUTPUTS:
            remove_if_exists(gas_outdir / filename)

        recoil_window = resolve_recoil_window_keV(
            gas_name,
            density_g_cm3,
            pressure_atm,
            recoil_window_df=recoil_window_df,
            diffusion_df=diffusion_df,
        )
        T_low_keV = recoil_window.low_keV
        T_high_keV = recoil_window.high_keV
        electron_density, molar_mass, electrons_per_mixture = electron_density_cm3(
            gas_name,
            density_g_cm3,
        )
        target_electrons = electron_density * volume_cm3

        reco_df, metrics = build_estimated_spectra_for_gas(
            flux_df=flux_df,
            target_electrons=target_electrons,
            T_low_keV=T_low_keV,
            T_high_keV=T_high_keV,
            response=response,
            t_grid_points=args.t_grid_points,
            max_true_energy_points=args.max_true_energy_points,
        )
        reco_df.to_csv(gas_outdir / "reco_neutrino_energy_spectra.csv", index=False)
        write_reco_plot(
            gas_outdir=gas_outdir,
            gas_label=gas_label,
            response=response,
            threshold_keV=T_low_keV,
            reco_df=reco_df,
        )

        if cevns_config.enabled:
            recoil_df, enu_df, cevns_summary_df = compute_cevns_spectra_for_gas(
                gas_name=gas_name,
                gas_label=gas_label,
                density_g_cm3=density_g_cm3,
                volume_cm3=volume_cm3,
                E_MeV=flux_df["E_MeV"].to_numpy(dtype=float),
                active_flux_or_fluence=cevns_flux,
                config=cevns_config,
                t_grid_points=args.t_grid_points,
                quantity="rate",
            )
            cevns_reco_df = build_cevns_lower_bound_reco(
                recoil_df,
                config=cevns_config,
                quantity="rate",
            )
            cevns_gas_outdir = cevns_outdir / gas_slug
            save_cevns_tables(
                cevns_gas_outdir,
                recoil_df=recoil_df,
                enu_df=enu_df,
                reco_df=cevns_reco_df,
                skip_save=args.skip_cevns_save,
            )
            write_cevns_reco_plot(
                cevns_gas_outdir,
                gas_label=gas_label,
                reco_df=cevns_reco_df,
                quantity="rate",
                skip_plots=args.skip_cevns_plots,
            )
            cevns_summary_df["output_subdir"] = str(Path("cevns") / gas_slug)
            cevns_summary_df["reconstruction_note"] = cevns_config.reconstruction_note
            cevns_summary_rows.append(cevns_summary_df)

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
                "detector_name": geometry.name,
                "configured_volume_m3": geometry.volume_m3,
                "volume_source": volume_source,
                "volume_m3": volume_cm3 / 1.0e6,
                "response_name": response.name,
                "energy_resolution_at_threshold_frac": response.energy_resolution_at_threshold_frac,
                "energy_resolution_constant_frac": response.energy_resolution_constant_frac,
                "angular_resolution_at_threshold_deg": response.angular_resolution_at_threshold_deg,
                "angular_resolution_constant_deg": response.angular_resolution_constant_deg,
                "reco_energy_bins": response.reco_energy_bins,
                "t_grid_points": args.t_grid_points,
                "max_true_energy_points": args.max_true_energy_points,
                "molar_mass_g_mol": molar_mass,
                "electrons_per_mixture": electrons_per_mixture,
                "electron_density_cm^-3": electron_density,
                "target_electrons": target_electrons,
                "total_rate_true_s^-1": metrics["total_rate_true_s^-1"],
                "total_rate_true_year^-1": metrics["total_rate_true_year^-1"],
                "total_rate_reco_energy_only_s^-1": metrics["total_rate_reco_energy_only_s^-1"],
                "total_rate_reco_energy_only_year^-1": metrics["total_rate_reco_energy_only_year^-1"],
                "total_rate_reco_energy_direction_s^-1": metrics["total_rate_reco_energy_direction_s^-1"],
                "total_rate_reco_energy_direction_year^-1": metrics["total_rate_reco_energy_direction_year^-1"],
                "reco_energy_min_MeV": metrics["reco_energy_min_MeV"],
                "reco_energy_max_MeV": metrics["reco_energy_max_MeV"],
                "gas_label": gas_label,
                "output_subdir": gas_slug,
            }
        )
        progress.update(gas_index, gas_label)
    progress.close()

    summary_df = pd.DataFrame(summary_rows)
    summary_path = outdir / "reco_neutrino_energy_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    cevns_summary_df = write_cevns_summary(
        cevns_outdir,
        cevns_summary_rows,
        quantity="rate",
        skip_save=args.skip_cevns_save,
        skip_plots=args.skip_cevns_plots,
        filename="cevns_reco_energy_min_summary.csv",
    ) if cevns_config.enabled else pd.DataFrame()

    print(f"Read {len(flux_df)} flux points from: {flux_csv}")
    print(f"Read {len(gas_df)} gas entries from: {gas_csv}")
    print(f"Read {len(diffusion_df)} diffusion rows from: {diffusion_csv}")
    print(
        f"Using detector geometry '{geometry.name}' = "
        f"{geometry.length_m:g} x {geometry.width_m:g} x {geometry.height_m:g} m^3 "
        f"(volume = {volume_cm3 / 1.0e6:g} m^3)"
    )
    print(
        f"Using response '{response.name}' with "
        f"sigma_E/E(thr)={response.energy_resolution_at_threshold_frac:g}, "
        f"sigma_E/E(const)={response.energy_resolution_constant_frac:g}, "
        f"sigma_theta(thr)={response.angular_resolution_at_threshold_deg:g} deg, "
        f"sigma_theta(const)={response.angular_resolution_constant_deg:g} deg, "
        f"reco bins={response.reco_energy_bins}, "
        f"max true-energy points={args.max_true_energy_points}"
    )
    print(f"Saved reconstructed-spectrum outputs in: {outdir.resolve()}")
    print(f"Saved summary: {summary_path.resolve()}")
    if cevns_config.enabled:
        print(f"CEvNS lower-bound reconstruction outputs directory: {cevns_outdir.resolve()}")
        print(cevns_config.reconstruction_note)
        if not cevns_summary_df.empty:
            print(
                "CEvNS total solar rate across all gas-aggregated summary rows: "
                f"{cevns_summary_df['total_rate_s-1'].sum():.6e} s^-1"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
