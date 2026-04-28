#!/usr/bin/env python3
"""
Shared detector and gas-target helpers for the CYGNUS-K research scripts.

The repository keeps its runnable analysis scripts at the root for convenience.
This module holds the reusable detector-side pieces that are needed by the
solar and supernova workflows: gas densities, target-electron counting,
range/diffusion recoil windows, detector geometry, and simple output labels.
"""

from __future__ import annotations

import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


AVOGADRO = 6.02214076e23
PLOT_DPI = 180

REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_FLUX_CSV = REPO_ROOT / "solar_neutrino_fluxes" / "solar_neutrino_total_flux.csv"
DEFAULT_GAS_CSV = REPO_ROOT / "DetectorNumbers" / "GasDensities.csv"
DEFAULT_RANGE_CSV = REPO_ROOT / "DetectorNumbers" / "electron_range_energy_table.csv"
DEFAULT_DIFFUSION_CSV = REPO_ROOT / "DetectorNumbers" / "diffusion_2sigma_50cm_summary.csv"
DEFAULT_GEOMETRY_CONFIG = REPO_ROOT / "detector_geometry.json"


@dataclass(frozen=True)
class GasSpecies:
    """Molecular properties used to convert gas density into electron density."""

    molar_mass_g_mol: float
    electrons_per_molecule: float


@dataclass(frozen=True)
class DetectorGeometry:
    """Rectangular fiducial detector geometry read from JSON."""

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
    """
    Gas-dependent accepted recoil-energy window.

    Energies are stored in keV. The lower threshold begins at the electron
    energy whose CSDA range is 1 mm and is raised when the selected diffusion
    summary implies a larger 2-sigma longitudinal spread over 50 cm.
    """

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
    """Small stderr progress bar for long per-gas loops."""

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

# Fractions are treated as molecular/volume fractions for ideal-gas mixtures.
GAS_MIXTURES = {
    "CF4": {"CF4": 1.0},
    "HeCF4(60% He 40% CF4)": {"He": 0.60, "CF4": 0.40},
    "HeCF4CH4(58% He 37% CF4 5% CH4)": {"He": 0.58, "CF4": 0.37, "CH4": 0.05},
}

STOPPING_POWER_FILES = {
    "CF4": REPO_ROOT / "DetectorNumbers" / "CF4_Stopping_power_electrons.csv",
    "HeCF4(60% He 40% CF4)": REPO_ROOT / "DetectorNumbers" / "HeCF4_Stopping_power_electrons.csv",
    "HeCF4CH4(58% He 37% CF4 5% CH4)": REPO_ROOT
    / "DetectorNumbers"
    / "HeCF4CH4_Stopping_power_electrons.csv",
}

DIFFUSION_GAS_PATTERNS = {
    "CF4": "cf4_100_{pressure_mbar}mbar293K",
    "HeCF4(60% He 40% CF4)": "he-cf4_60-40_{pressure_mbar}mbar293K",
    "HeCF4CH4(58% He 37% CF4 5% CH4)": "he-cf4-ch4_58-37-5_{pressure_mbar}mbar293K",
}


def read_gas_density_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()

    required = ["Gas", "density @ 293 K (g/cm3)", "Pressure (atm)"]
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError("Missing required columns in gas-density table: " + ", ".join(missing))

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
        raise ValueError("Missing required columns in recoil-window table: " + ", ".join(missing))

    df = df.copy()
    df["Gas"] = df["Gas"].astype(str).str.strip()
    for column in required[1:]:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    df = df.dropna(subset=required)
    return df.reset_index(drop=True)


def read_diffusion_summary_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()

    required = ["gas", "point", "electric_field(V/cm)", "DL_2sigma(mm)"]
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError("Missing required columns in diffusion summary table: " + ", ".join(missing))

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
        raise ValueError("Missing required keys in detector geometry config: " + ", ".join(missing))

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
    """
    Return electron density, mean molar mass, and electrons per mixture molecule.

    The gas-mixture fractions are interpreted as molecular fractions. For an
    ideal gas this is equivalent to volume fraction.
    """

    if gas_name not in GAS_MIXTURES:
        raise KeyError(f"Gas mixture '{gas_name}' is not configured in GAS_MIXTURES.")

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
    """
    Return CSDA-like electron range in mm and kinetic energy in MeV.

    The stopping-power tables store collisional mass stopping power in
    MeV cm^2/g. The code multiplies by gas density to obtain linear stopping
    power and integrates 1/(dE/dx). A low-energy 1/E extrapolation is appended
    before integration to keep the tabulated threshold interpolation stable.
    """

    if gas_name not in STOPPING_POWER_FILES:
        raise KeyError(f"Gas mixture '{gas_name}' is not configured for stopping-power lookup.")

    path = STOPPING_POWER_FILES[gas_name]
    if not path.exists():
        raise FileNotFoundError(f"Stopping-power CSV not found: {path}")

    sp_data = pd.read_csv(path, skiprows=4)
    energy_mev = pd.to_numeric(sp_data.iloc[:, 0], errors="coerce").to_numpy(dtype=float)
    collision_sp = pd.to_numeric(sp_data.iloc[:, 1], errors="coerce").to_numpy(dtype=float)
    valid = (
        np.isfinite(energy_mev)
        & np.isfinite(collision_sp)
        & (energy_mev > 0.0)
        & (collision_sp > 0.0)
    )
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
