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
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import scatteringPlots as sp


AVOGADRO = 6.02214076e23
SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0
PLOT_DPI = 180

REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_FLUX_CSV = REPO_ROOT / "solar_neutrino_fluxes" / "solar_neutrino_total_flux.csv"
DEFAULT_GAS_CSV = REPO_ROOT / "DetectorNumbers" / "GasDensities.csv"
DEFAULT_RANGE_CSV = REPO_ROOT / "DetectorNumbers" / "electron_range_energy_table.csv"
DEFAULT_DIFFUSION_CSV = REPO_ROOT / "DetectorNumbers" / "diffusion_2sigma_50cm_summary.csv"
DEFAULT_GEOMETRY_CONFIG = REPO_ROOT / "detector_geometry.json"
DEFAULT_OUTDIR = REPO_ROOT / "solar_nu_gas_target_rates"

TRAPEZOID = getattr(np, "trapezoid", np.trapz)


@dataclass(frozen=True)
class GasSpecies:
    molar_mass_g_mol: float
    electrons_per_molecule: float


@dataclass(frozen=True)
class DetectorGeometry:
    name: str
    length_m: float
    width_m: float
    height_m: float

    @property
    def volume_m3(self) -> float:
        return self.length_m * self.width_m * self.height_m

    @property
    def volume_cm3(self) -> float:
        return self.volume_m3 * 1.0e6


@dataclass(frozen=True)
class RecoilWindow:
    low_keV: float
    high_keV: float
    low_range_mm: float
    base_low_keV: float
    base_low_range_mm: float
    diffusion_dl_2sigma_mm: float
    diffusion_field_v_cm: float
    diffusion_point: str
    threshold_source: str


class ProgressBar:
    def __init__(self, total: int, prefix: str, width: int = 32, enabled: bool = True):
        self.total = max(int(total), 0)
        self.prefix = prefix
        self.width = max(int(width), 8)
        self.enabled = enabled
        self._last_len = 0

    def update(self, current: int, label: str = "") -> None:
        if not self.enabled:
            return

        current = min(max(int(current), 0), self.total)
        frac = 1.0 if self.total == 0 else current / self.total
        filled = min(self.width, int(round(frac * self.width)))
        bar = "#" * filled + "-" * (self.width - filled)
        suffix = f" {label}" if label else ""
        message = f"\r{self.prefix} [{bar}] {current}/{self.total} ({frac * 100:5.1f}%){suffix}"
        padding = " " * max(self._last_len - len(message), 0)
        print(message + padding, end="", file=sys.stderr, flush=True)
        self._last_len = len(message)

    def close(self) -> None:
        if self.enabled and self._last_len:
            print(file=sys.stderr, flush=True)


SPECIES = {
    "He": GasSpecies(molar_mass_g_mol=4.002602, electrons_per_molecule=2.0),
    "CF4": GasSpecies(
        molar_mass_g_mol=12.011 + 4.0 * 18.998403163,
        electrons_per_molecule=6.0 + 4.0 * 9.0,
    ),
    "CH4": GasSpecies(
        molar_mass_g_mol=12.011 + 4.0 * 1.00794,
        electrons_per_molecule=6.0 + 4.0 * 1.0,
    ),
}

GAS_MIXTURES = {
    "CF4": {"CF4": 1.0},
    "HeCF4(60% He 40% CF4)": {"He": 0.60, "CF4": 0.40},
    "HeCF4CH4(58% He 37% CF4 5% CH4)": {"He": 0.58, "CF4": 0.37, "CH4": 0.05},
}

STOPPING_POWER_FILES = {
    "CF4": REPO_ROOT / "DetectorNumbers" / "CF4_Stopping_power_electrons.csv",
    "HeCF4(60% He 40% CF4)": REPO_ROOT / "DetectorNumbers" / "HeCF4_Stopping_power_electrons.csv",
    "HeCF4CH4(58% He 37% CF4 5% CH4)": REPO_ROOT / "DetectorNumbers" / "HeCF4CH4_Stopping_power_electrons.csv",
}

DIFFUSION_GAS_PATTERNS = {
    "CF4": "cf4_100_{pressure_mbar}mbar293K",
    "HeCF4(60% He 40% CF4)": "he-cf4_60-40_{pressure_mbar}mbar293K",
    "HeCF4CH4(58% He 37% CF4 5% CH4)": "he-cf4-ch4_58-37-5_{pressure_mbar}mbar293K",
}


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
    return parser.parse_args()


def read_gas_density_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()

    required = ["Gas", "density @ 293 K (g/cm3)", "Pressure (atm)"]
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(
            "Missing required columns in gas-density table: "
            + ", ".join(missing)
        )

    df = df.copy()
    df["Gas"] = df["Gas"].astype(str).str.strip()
    df["density @ 293 K (g/cm3)"] = pd.to_numeric(
        df["density @ 293 K (g/cm3)"], errors="coerce"
    )
    df["Pressure (atm)"] = pd.to_numeric(df["Pressure (atm)"], errors="coerce")
    df = df.dropna(subset=["Gas", "density @ 293 K (g/cm3)", "Pressure (atm)"])
    return df.reset_index(drop=True)


def read_recoil_window_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()

    required = [
        "Gas",
        "density_g_cm3",
        "energy_min_at_1mm_keV",
        "energy_max_at_1m_keV",
    ]
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(
            "Missing required columns in recoil-window table: "
            + ", ".join(missing)
        )

    df = df.copy()
    df["Gas"] = df["Gas"].astype(str).str.strip()
    for column in required[1:]:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    df = df.dropna(subset=required)
    return df.reset_index(drop=True)


def read_diffusion_summary_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()

    required = [
        "gas",
        "point",
        "electric_field(V/cm)",
        "DL_2sigma(mm)",
    ]
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(
            "Missing required columns in diffusion summary table: "
            + ", ".join(missing)
        )

    df = df.copy()
    df["gas"] = df["gas"].astype(str).str.strip()
    df["point"] = df["point"].astype(str).str.strip()
    df["electric_field(V/cm)"] = pd.to_numeric(df["electric_field(V/cm)"], errors="coerce")
    df["DL_2sigma(mm)"] = pd.to_numeric(df["DL_2sigma(mm)"], errors="coerce")
    df = df.dropna(subset=required)
    return df.reset_index(drop=True)


def read_detector_geometry_config(path: Path) -> DetectorGeometry:
    with open(path, "r", encoding="utf-8") as handle:
        data = json.load(handle)

    required = ["length_m", "width_m", "height_m"]
    missing = [key for key in required if key not in data]
    if missing:
        raise ValueError(
            "Missing required keys in detector geometry config: "
            + ", ".join(missing)
        )

    geometry = DetectorGeometry(
        name=str(data.get("name", path.stem)),
        length_m=float(data["length_m"]),
        width_m=float(data["width_m"]),
        height_m=float(data["height_m"]),
    )
    if geometry.length_m <= 0.0 or geometry.width_m <= 0.0 or geometry.height_m <= 0.0:
        raise ValueError("Detector geometry dimensions must all be positive.")
    return geometry


def electron_density_cm3(gas_name: str, density_g_cm3: float) -> tuple[float, float, float]:
    if gas_name not in GAS_MIXTURES:
        raise KeyError(
            f"Gas mixture '{gas_name}' is not configured in GAS_MIXTURES."
        )

    mixture = GAS_MIXTURES[gas_name]
    molar_mass = 0.0
    electrons_per_mixture = 0.0

    for species_name, fraction in mixture.items():
        species = SPECIES[species_name]
        molar_mass += fraction * species.molar_mass_g_mol
        electrons_per_mixture += fraction * species.electrons_per_molecule

    molecules_per_cm3 = density_g_cm3 * AVOGADRO / molar_mass
    electrons_per_cm3 = molecules_per_cm3 * electrons_per_mixture
    return electrons_per_cm3, molar_mass, electrons_per_mixture


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


def safe_slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def gas_entry_slug(gas_name: str, pressure_atm: float) -> str:
    return safe_slug(f"{gas_name}_{pressure_atm:g}atm")


def gas_entry_label(gas_name: str, pressure_atm: float) -> str:
    return f"{gas_name} @ {pressure_atm:g} atm"


def diffusion_gas_key(gas_name: str, pressure_atm: float) -> str:
    if gas_name not in DIFFUSION_GAS_PATTERNS:
        raise KeyError(f"Gas mixture '{gas_name}' is not configured for diffusion matching.")
    pressure_mbar = int(round(pressure_atm * 1000.0))
    return DIFFUSION_GAS_PATTERNS[gas_name].format(pressure_mbar=pressure_mbar)


def cumulative_trapezoid(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    if len(y) != len(x):
        raise ValueError("x and y must have the same length.")
    if len(y) == 0:
        return np.array([], dtype=float)
    dx = np.diff(x)
    area = 0.5 * (y[1:] + y[:-1]) * dx
    return np.concatenate([[0.0], np.cumsum(area)])


def read_stopping_power_curve(gas_name: str, density_g_cm3: float) -> tuple[np.ndarray, np.ndarray]:
    if gas_name not in STOPPING_POWER_FILES:
        raise KeyError(f"Gas mixture '{gas_name}' is not configured for stopping-power lookup.")

    path = STOPPING_POWER_FILES[gas_name]
    if not path.exists():
        raise FileNotFoundError(f"Stopping-power CSV not found: {path}")

    sp_data = pd.read_csv(path, skiprows=4)
    energy_mev = pd.to_numeric(sp_data.iloc[:, 0], errors="coerce").to_numpy(dtype=float)
    collision_sp = pd.to_numeric(sp_data.iloc[:, 1], errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(energy_mev) & np.isfinite(collision_sp) & (energy_mev > 0.0) & (collision_sp > 0.0)
    energy_mev = energy_mev[valid]
    collision_sp = collision_sp[valid]
    if len(energy_mev) < 2:
        raise ValueError(f"Stopping-power CSV has too few valid rows: {path}")

    order = np.argsort(energy_mev)
    energy_mev = energy_mev[order]
    collision_sp = collision_sp[order]

    e_min = energy_mev[0]
    sp_min = collision_sp[0]
    energy_mev = np.concatenate([[e_min / 10.0], energy_mev])
    collision_sp = np.concatenate([[sp_min * 10.0], collision_sp])

    linear_sp = collision_sp * density_g_cm3
    inverse_sp = 1.0 / linear_sp
    range_mm = cumulative_trapezoid(inverse_sp, energy_mev) * 10.0
    return range_mm, energy_mev


def energy_at_electron_range_keV(
    gas_name: str,
    density_g_cm3: float,
    target_range_mm: float,
) -> float:
    range_mm, energy_mev = read_stopping_power_curve(gas_name, density_g_cm3)
    order = np.argsort(range_mm)
    r_sorted = range_mm[order]
    e_sorted = energy_mev[order]
    r_unique, idx = np.unique(r_sorted, return_index=True)
    e_unique = e_sorted[idx]

    if target_range_mm < r_unique[0] or target_range_mm > r_unique[-1]:
        raise ValueError(
            f"Cannot interpolate electron energy for range {target_range_mm:.4g} mm "
            f"in gas '{gas_name}' at density {density_g_cm3:.8g} g/cm^3. "
            f"Available range is {r_unique[0]:.4g} to {r_unique[-1]:.4g} mm."
        )
    return float(np.interp(target_range_mm, r_unique, e_unique) * 1.0e3)


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


def find_recoil_window_keV(
    gas_name: str,
    density_g_cm3: float,
    recoil_window_df: pd.DataFrame,
) -> tuple[float, float]:
    gas_mask = recoil_window_df["Gas"].astype(str).str.strip() == gas_name
    density_mask = np.isclose(
        recoil_window_df["density_g_cm3"].to_numpy(dtype=float),
        density_g_cm3,
        rtol=1.0e-6,
        atol=1.0e-9,
    )
    matches = recoil_window_df[gas_mask & density_mask]

    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one recoil-window match for gas '{gas_name}' at density "
            f"{density_g_cm3:.8g} g/cm^3, found {len(matches)}."
        )

    row = matches.iloc[0]
    t_low_keV = float(row["energy_min_at_1mm_keV"])
    t_high_keV = float(row["energy_max_at_1m_keV"])
    if t_low_keV < 0.0 or t_high_keV <= t_low_keV:
        raise ValueError(
            f"Invalid recoil window for gas '{gas_name}' at density {density_g_cm3:.8g} g/cm^3."
        )
    return t_low_keV, t_high_keV


def resolve_recoil_window_keV(
    gas_name: str,
    density_g_cm3: float,
    pressure_atm: float,
    recoil_window_df: pd.DataFrame,
    diffusion_df: pd.DataFrame,
) -> RecoilWindow:
    base_low_keV, high_keV = find_recoil_window_keV(
        gas_name,
        density_g_cm3,
        recoil_window_df=recoil_window_df,
    )

    diffusion_key = diffusion_gas_key(gas_name, pressure_atm)
    matches = diffusion_df[
        (diffusion_df["gas"].astype(str).str.strip() == diffusion_key)
        & (diffusion_df["electric_field(V/cm)"].to_numpy(dtype=float) <= 2000.0)
    ].copy()

    if matches.empty:
        raise ValueError(
            f"No diffusion rows found for '{diffusion_key}' with electric field <= 2 kV/cm."
        )

    max_idx = matches["DL_2sigma(mm)"].astype(float).idxmax()
    diffusion_row = matches.loc[max_idx]
    diffusion_dl_mm = float(diffusion_row["DL_2sigma(mm)"])
    diffusion_field_v_cm = float(diffusion_row["electric_field(V/cm)"])
    diffusion_point = str(diffusion_row["point"])

    low_range_mm = max(1.0, diffusion_dl_mm)
    if low_range_mm > 1.0:
        low_keV = energy_at_electron_range_keV(gas_name, density_g_cm3, low_range_mm)
        threshold_source = "diffusion_DL_2sigma"
    else:
        low_keV = base_low_keV
        threshold_source = "electron_range_1mm"

    if low_keV < 0.0 or high_keV <= low_keV:
        raise ValueError(
            f"Invalid recoil window for gas '{gas_name}' at density {density_g_cm3:.8g} g/cm^3: "
            f"{low_keV:.6g} to {high_keV:.6g} keV."
        )

    return RecoilWindow(
        low_keV=float(low_keV),
        high_keV=float(high_keV),
        low_range_mm=float(low_range_mm),
        base_low_keV=float(base_low_keV),
        base_low_range_mm=1.0,
        diffusion_dl_2sigma_mm=diffusion_dl_mm,
        diffusion_field_v_cm=diffusion_field_v_cm,
        diffusion_point=diffusion_point,
        threshold_source=threshold_source,
    )


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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
