#!/usr/bin/env python3
"""
Shared CEvNS pipeline utilities.

This module keeps CEvNS integration out of the existing nu-e driver logic. It
handles configuration, flavor-independent neutral-current flux sums, true
nuclear-recoil spectra, accepted neutrino-energy spectra, and the energy-only
lower-bound estimator E_nu_min(T_N).
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from math import erf
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import cevns
from detector_model import PLOT_DPI, isotope_target_counts


TRAPEZOID = getattr(np, "trapezoid", np.trapz)
CEVNS_AXIAL_MODELS = ("none", "hoferichter_19f_fast", "hoferichter_19f_central", "toy")
CEVNS_FORM_FACTORS = ("default", "helm", "pointlike")


@dataclass(frozen=True)
class CevnsConfig:
    enabled: bool = False
    nr_threshold_keV: float = 1.0
    nr_max_keV: float | None = None
    form_factor: str = "default"
    axial_model: str = "hoferichter_19f_fast"
    reco_energy_bins: int = 100
    energy_resolution_at_threshold_frac: float | None = None
    energy_resolution_constant_frac: float | None = None

    @property
    def has_energy_resolution(self) -> bool:
        return (
            self.energy_resolution_at_threshold_frac is not None
            and self.energy_resolution_constant_frac is not None
        )

    @property
    def reconstruction_note(self) -> str:
        if self.has_energy_resolution:
            return (
                "CEvNS energy-only lower-bound estimator E_nu_min(T_N), "
                "smeared with the CEvNS nuclear-recoil energy-resolution model."
            )
        return (
            "CEvNS energy-only lower-bound estimator E_nu_min(T_N); no CEvNS "
            "nuclear-recoil energy resolution is configured, so this spectrum is unsmeared."
        )


def parse_optional_float(value: str) -> float | None:
    if str(value).strip().lower() in {"", "none", "null"}:
        return None
    return float(value)


def add_cevns_cli_args(parser) -> None:
    group = parser.add_argument_group("CEvNS nuclear-recoil options")
    enabled_group = group.add_mutually_exclusive_group()
    enabled_group.add_argument(
        "--enable-cevns",
        action="store_true",
        help="Enable optional CEvNS nuclear-recoil outputs.",
    )
    enabled_group.add_argument(
        "--disable-cevns",
        action="store_true",
        help="Disable optional CEvNS nuclear-recoil outputs.",
    )
    group.add_argument(
        "--cevns-threshold-kev",
        type=float,
        default=None,
        help="True nuclear-recoil CEvNS threshold in keV.",
    )
    group.add_argument(
        "--cevns-max-kev",
        type=parse_optional_float,
        default=None,
        help="Optional true nuclear-recoil CEvNS maximum in keV; use 'none' for kinematic maximum.",
    )
    group.add_argument(
        "--cevns-form-factor",
        choices=CEVNS_FORM_FACTORS,
        default=None,
        help="CEvNS vector form-factor mode.",
    )
    group.add_argument(
        "--cevns-axial-model",
        choices=CEVNS_AXIAL_MODELS,
        default=None,
        help="CEvNS axial model. Hydrogen axial remains unimplemented for all choices except none.",
    )
    group.add_argument(
        "--skip-cevns-plots",
        action="store_true",
        help="Compute CEvNS but do not write CEvNS plots.",
    )
    group.add_argument(
        "--skip-cevns-save",
        action="store_true",
        help="Compute CEvNS but do not write CEvNS CSV outputs.",
    )


def _coerce_optional_float(value) -> float | None:
    if value is None:
        return None
    return float(value)


def load_cevns_config(path: Path | None, args, default_enabled: bool = False) -> CevnsConfig:
    block = {}
    if path is not None and path.exists():
        with open(path, "r", encoding="utf-8") as handle:
            data = json.load(handle)
        block = data.get("cevns", {}) or {}

    enabled = bool(block.get("enabled", default_enabled))
    if getattr(args, "enable_cevns", False):
        enabled = True
    if getattr(args, "disable_cevns", False):
        enabled = False

    nr_threshold_keV = float(block.get("nr_threshold_keV", 1.0))
    if getattr(args, "cevns_threshold_kev", None) is not None:
        nr_threshold_keV = float(args.cevns_threshold_kev)

    nr_max_keV = _coerce_optional_float(block.get("nr_max_keV", None))
    if hasattr(args, "cevns_max_kev") and args.cevns_max_kev is not None:
        nr_max_keV = args.cevns_max_kev

    form_factor = str(block.get("form_factor", "default"))
    if getattr(args, "cevns_form_factor", None) is not None:
        form_factor = str(args.cevns_form_factor)

    axial_model = str(block.get("axial_model", "hoferichter_19f_fast"))
    if getattr(args, "cevns_axial_model", None) is not None:
        axial_model = str(args.cevns_axial_model)

    reco_energy_bins = int(block.get("reco_energy_bins", 100))
    energy_resolution_at_threshold_frac = _coerce_optional_float(
        block.get("energy_resolution_at_threshold_frac", None)
    )
    energy_resolution_constant_frac = _coerce_optional_float(
        block.get("energy_resolution_constant_frac", None)
    )

    config = CevnsConfig(
        enabled=enabled,
        nr_threshold_keV=nr_threshold_keV,
        nr_max_keV=nr_max_keV,
        form_factor=form_factor,
        axial_model=axial_model,
        reco_energy_bins=reco_energy_bins,
        energy_resolution_at_threshold_frac=energy_resolution_at_threshold_frac,
        energy_resolution_constant_frac=energy_resolution_constant_frac,
    )
    validate_cevns_config(config)
    return config


def validate_cevns_config(config: CevnsConfig) -> None:
    if config.nr_threshold_keV < 0.0:
        raise ValueError("CEvNS nuclear-recoil threshold must be non-negative.")
    if config.nr_max_keV is not None and config.nr_max_keV <= config.nr_threshold_keV:
        raise ValueError("CEvNS nuclear-recoil max must be greater than the threshold.")
    if config.form_factor not in CEVNS_FORM_FACTORS:
        raise ValueError(f"Unsupported CEvNS form factor: {config.form_factor}")
    if config.axial_model not in CEVNS_AXIAL_MODELS:
        raise ValueError(f"Unsupported CEvNS axial model: {config.axial_model}")
    if config.reco_energy_bins < 20:
        raise ValueError("CEvNS reco_energy_bins must be at least 20.")
    if (
        config.energy_resolution_at_threshold_frac is not None
        and config.energy_resolution_at_threshold_frac < 0.0
    ):
        raise ValueError("CEvNS energy-resolution threshold term must be non-negative.")
    if (
        config.energy_resolution_constant_frac is not None
        and config.energy_resolution_constant_frac < 0.0
    ):
        raise ValueError("CEvNS energy-resolution constant term must be non-negative.")


def solar_active_flux(flux_df: pd.DataFrame) -> np.ndarray:
    """Return the SM neutral-current active solar flux sum."""

    columns = ["nu_e_dphi_dE", "nu_mu_dphi_dE", "nu_tau_dphi_dE"]
    missing = [column for column in columns if column not in flux_df.columns]
    if missing:
        raise ValueError("Missing solar flux columns for CEvNS: " + ", ".join(missing))
    return sum(flux_df[column].to_numpy(dtype=float) for column in columns)


def supernova_active_fluence(fluence_df: pd.DataFrame) -> np.ndarray:
    """Return the SM neutral-current active supernova neutrino plus antineutrino sum."""

    columns = ["nue", "anue", "numu", "anumu", "nutau", "anutau"]
    missing = [column for column in columns if column not in fluence_df.columns]
    if missing:
        raise ValueError("Missing supernova fluence columns for CEvNS: " + ", ".join(missing))
    return sum(fluence_df[column].to_numpy(dtype=float) for column in columns)


def bin_widths_from_centers(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    if len(x) < 2:
        return np.ones_like(x)
    dx = np.empty_like(x)
    dx[1:-1] = 0.5 * (x[2:] - x[:-2])
    dx[0] = x[1] - x[0]
    dx[-1] = x[-1] - x[-2]
    return dx


def recoil_grid_keV(
    E_MeV: np.ndarray,
    target: cevns.NuclearTarget,
    config: CevnsConfig,
    n_points: int,
) -> np.ndarray:
    if n_points < 2:
        raise ValueError("CEvNS recoil grid needs at least 2 points.")
    t_kin_max = float(cevns.tmax_nuclear_keV(float(np.max(E_MeV)), target))
    t_upper = t_kin_max if config.nr_max_keV is None else min(config.nr_max_keV, t_kin_max)
    t_low = config.nr_threshold_keV
    if t_upper <= t_low:
        return np.array([t_low, max(t_low, t_upper)], dtype=float)
    if t_low > 0.0:
        return np.logspace(np.log10(t_low), np.log10(t_upper), n_points)
    positive = np.logspace(-6, np.log10(t_upper), n_points - 1)
    return np.concatenate([[0.0], positive])


def sigma_accepted_components(
    E_MeV: np.ndarray,
    target: cevns.NuclearTarget,
    config: CevnsConfig,
    n_t_points: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    E = np.asarray(E_MeV, dtype=float)
    sigma_vector = np.zeros_like(E)
    sigma_axial = np.zeros_like(E)

    for idx, energy in enumerate(E):
        t_upper = float(cevns.tmax_nuclear_keV(energy, target))
        if config.nr_max_keV is not None:
            t_upper = min(t_upper, config.nr_max_keV)
        if t_upper <= config.nr_threshold_keV:
            continue

        T = np.linspace(config.nr_threshold_keV, t_upper, max(n_t_points, 2))
        sigma_vector[idx] = float(
            TRAPEZOID(
                cevns.dsigma_dT_vector(
                    energy,
                    T,
                    target,
                    form_factor=config.form_factor,
                ),
                T,
            )
        )
        sigma_axial[idx] = float(
            TRAPEZOID(
                cevns.dsigma_dT_axial(
                    energy,
                    T,
                    target,
                    axial_model=config.axial_model,
                ),
                T,
            )
        )

    sigma_total = sigma_vector + sigma_axial
    return sigma_vector, sigma_axial, sigma_total


def recoil_spectrum_components(
    E_MeV: np.ndarray,
    active_flux_or_fluence: np.ndarray,
    target: cevns.NuclearTarget,
    target_count: float,
    config: CevnsConfig,
    T_grid_keV: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    E = np.asarray(E_MeV, dtype=float)
    weights = np.asarray(active_flux_or_fluence, dtype=float) * bin_widths_from_centers(E)

    vector = np.zeros_like(T_grid_keV, dtype=float)
    axial = np.zeros_like(T_grid_keV, dtype=float)
    for idx, T_keV in enumerate(T_grid_keV):
        vector[idx] = target_count * float(
            np.sum(
                weights
                * cevns.dsigma_dT_vector(
                    E,
                    T_keV,
                    target,
                    form_factor=config.form_factor,
                )
            )
        )
        axial[idx] = target_count * float(
            np.sum(
                weights
                * cevns.dsigma_dT_axial(
                    E,
                    T_keV,
                    target,
                    axial_model=config.axial_model,
                )
            )
        )
    total = vector + axial
    return vector, axial, total


def compute_cevns_spectra_for_gas(
    gas_name: str,
    gas_label: str,
    density_g_cm3: float,
    volume_cm3: float,
    E_MeV: np.ndarray,
    active_flux_or_fluence: np.ndarray,
    config: CevnsConfig,
    t_grid_points: int,
    quantity: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compute CEvNS recoil, accepted-E_nu, and summary tables for one gas."""

    if quantity not in {"rate", "count"}:
        raise ValueError("quantity must be 'rate' or 'count'.")

    target_counts, _molar_mass, _isotopes_per_mixture = isotope_target_counts(
        gas_name,
        density_g_cm3,
        volume_cm3,
    )
    dE = bin_widths_from_centers(E_MeV)
    recoil_rows = []
    enu_rows = []
    summary_rows = []

    for isotope_key in sorted(target_counts, key=lambda key: cevns.NUCLEAR_TARGETS[key].A):
        target = cevns.NUCLEAR_TARGETS[isotope_key]
        number_of_targets = target_counts[isotope_key]
        T_grid = recoil_grid_keV(E_MeV, target, config, t_grid_points)

        dT_vector, dT_axial, dT_total = recoil_spectrum_components(
            E_MeV,
            active_flux_or_fluence,
            target,
            number_of_targets,
            config,
            T_grid,
        )
        sigma_vector, sigma_axial, sigma_total = sigma_accepted_components(
            E_MeV,
            target,
            config,
            n_t_points=t_grid_points,
        )
        accepted_vector = number_of_targets * active_flux_or_fluence * sigma_vector
        accepted_axial = number_of_targets * active_flux_or_fluence * sigma_axial
        accepted_total = number_of_targets * active_flux_or_fluence * sigma_total

        for T_keV, vector, axial, total in zip(T_grid, dT_vector, dT_axial, dT_total):
            recoil_rows.append(
                {
                    "gas_label": gas_label,
                    "gas": gas_name,
                    "isotope": target.label,
                    "isotope_key": isotope_key,
                    "A": target.A,
                    "Z": target.Z,
                    "N": target.N,
                    "T_N_keV": T_keV,
                    component_column(quantity, "recoil", "vector"): vector,
                    component_column(quantity, "recoil", "axial"): axial,
                    component_column(quantity, "recoil", "total"): total,
                }
            )

        for energy, vector, axial, total in zip(
            E_MeV,
            accepted_vector,
            accepted_axial,
            accepted_total,
        ):
            enu_rows.append(
                {
                    "gas_label": gas_label,
                    "gas": gas_name,
                    "isotope": target.label,
                    "isotope_key": isotope_key,
                    "E_nu_MeV": energy,
                    component_column(quantity, "enu", "vector"): vector,
                    component_column(quantity, "enu", "axial"): axial,
                    component_column(quantity, "enu", "total"): total,
                }
            )

        vector_integral = float(np.sum(accepted_vector * dE))
        axial_integral = float(np.sum(accepted_axial * dE))
        total_integral = vector_integral + axial_integral
        summary_rows.append(
            {
                "gas_label": gas_label,
                "gas": gas_name,
                "isotope": target.label,
                "isotope_key": isotope_key,
                "A": target.A,
                "Z": target.Z,
                "N": target.N,
                "J": target.J,
                "number_of_targets": number_of_targets,
                "nr_threshold_keV": config.nr_threshold_keV,
                "nr_max_keV": config.nr_max_keV,
                "form_factor": config.form_factor,
                "axial_model": config.axial_model,
                summary_column(quantity, "vector"): vector_integral,
                summary_column(quantity, "axial"): axial_integral,
                summary_column(quantity, "total"): total_integral,
                "axial_fraction": axial_integral / total_integral if total_integral > 0.0 else 0.0,
                "hydrogen_axial_warning": (
                    "hydrogen axial neutral-current elastic scattering is not yet implemented"
                    if isotope_key == "H1"
                    else ""
                ),
                "threshold_note": (
                    "true nuclear recoil threshold; quenching/electron-equivalent threshold "
                    "not implemented"
                ),
            }
        )

    return pd.DataFrame(recoil_rows), pd.DataFrame(enu_rows), pd.DataFrame(summary_rows)


def component_column(quantity: str, spectrum: str, component: str) -> str:
    if spectrum == "recoil":
        if quantity == "rate":
            return f"dR_dT_s-1_keV-1_{component}"
        return f"dN_dT_keV-1_{component}"
    if quantity == "rate":
        return f"accepted_rate_s-1_MeV-1_{component}"
    return f"accepted_counts_MeV-1_{component}"


def summary_column(quantity: str, component: str) -> str:
    if quantity == "rate":
        return f"{component}_rate_s-1"
    return f"{component}_counts"


def aggregate_summary_column(quantity: str, component: str) -> str:
    if component == "total":
        return summary_column(quantity, component)
    if quantity == "rate":
        return f"{component}_total_rate_s-1"
    return f"{component}_total_counts"


def reco_component_column(quantity: str, component: str, density: bool = False) -> str:
    if quantity == "rate":
        return (
            f"{component}_rate_s-1_MeV-1"
            if density
            else f"{component}_rate_s-1"
        )
    return f"{component}_count_MeV-1" if density else f"{component}_count"


def aggregate_cevns_recoil_spectrum(recoil_df: pd.DataFrame, quantity: str) -> pd.DataFrame:
    """Return gas-total dR/dT_N or dN/dT_N by interpolating isotope grids."""

    if recoil_df.empty:
        return pd.DataFrame()

    component_cols = [component_column(quantity, "recoil", c) for c in ("vector", "axial", "total")]
    rows = []
    for (gas_label, gas_name), gas_df in recoil_df.groupby(["gas_label", "gas"], sort=False):
        grid = np.unique(gas_df["T_N_keV"].to_numpy(dtype=float))
        grid = grid[np.isfinite(grid)]
        if len(grid) == 0:
            continue
        summed = {col: np.zeros_like(grid, dtype=float) for col in component_cols}
        for _isotope, isotope_df in gas_df.groupby("isotope_key", sort=False):
            order = np.argsort(isotope_df["T_N_keV"].to_numpy(dtype=float))
            t_values = isotope_df["T_N_keV"].to_numpy(dtype=float)[order]
            if len(t_values) == 0:
                continue
            for col in component_cols:
                y_values = isotope_df[col].to_numpy(dtype=float)[order]
                summed[col] += np.interp(grid, t_values, y_values, left=0.0, right=0.0)
        for idx, T_keV in enumerate(grid):
            row = {
                "gas_label": gas_label,
                "gas": gas_name,
                "T_N_keV": T_keV,
            }
            for col in component_cols:
                row[col] = summed[col][idx]
            rows.append(row)
    return pd.DataFrame(rows)


def aggregate_cevns_enu_spectrum(enu_df: pd.DataFrame, quantity: str) -> pd.DataFrame:
    """Return gas-total accepted E_nu spectrum."""

    if enu_df.empty:
        return pd.DataFrame()

    component_cols = [component_column(quantity, "enu", c) for c in ("vector", "axial", "total")]
    group_cols = ["gas_label", "gas", "E_nu_MeV"]
    return (
        enu_df[group_cols + component_cols]
        .groupby(group_cols, as_index=False, sort=False)
        .sum()
    )


def aggregate_cevns_reco_spectrum(reco_df: pd.DataFrame, quantity: str) -> pd.DataFrame:
    """Return gas-total energy-only E_nu_min spectrum."""

    if reco_df is None or reco_df.empty:
        return pd.DataFrame()

    value_cols = [reco_component_column(quantity, c) for c in ("vector", "axial", "total")]
    density_cols = [reco_component_column(quantity, c, density=True) for c in ("vector", "axial", "total")]
    group_cols = [
        "gas_label",
        "gas",
        "E_nu_reco_min_MeV",
        "E_nu_bin_low_MeV",
        "E_nu_bin_high_MeV",
        "bin_width_MeV",
        "reconstruction_note",
    ]
    aggregate = (
        reco_df[group_cols + value_cols]
        .groupby(group_cols, as_index=False, sort=False, dropna=False)
        .sum()
    )
    widths = aggregate["bin_width_MeV"].to_numpy(dtype=float)
    for value_col, density_col in zip(value_cols, density_cols):
        values = aggregate[value_col].to_numpy(dtype=float)
        aggregate[density_col] = np.divide(
            values,
            widths,
            out=np.zeros_like(values),
            where=widths > 0.0,
        )
    return aggregate


def aggregate_cevns_summary(summary_df: pd.DataFrame, quantity: str) -> pd.DataFrame:
    """Return the primary one-row-per-gas CEvNS summary."""

    if summary_df.empty:
        return pd.DataFrame()

    vector_col = summary_column(quantity, "vector")
    axial_col = summary_column(quantity, "axial")
    total_col = summary_column(quantity, "total")
    group_cols = [
        col
        for col in (
            "model_name",
            "fluence_csv",
            "gas_label",
            "gas",
            "nr_threshold_keV",
            "nr_max_keV",
            "form_factor",
            "axial_model",
            "threshold_note",
            "reconstruction_note",
            "output_subdir",
        )
        if col in summary_df.columns
    ]

    rows = []
    for keys, group in summary_df.groupby(group_cols, sort=False, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = dict(zip(group_cols, keys))
        isotopes = [str(value) for value in group["isotope"].tolist()] if "isotope" in group else []
        vector_total = float(group[vector_col].sum())
        axial_total = float(group[axial_col].sum())
        total = float(group[total_col].sum())
        row.update(
            {
                "isotopes": ",".join(isotopes),
                "number_of_targets_total": float(group["number_of_targets"].sum())
                if "number_of_targets" in group
                else np.nan,
                aggregate_summary_column(quantity, "vector"): vector_total,
                aggregate_summary_column(quantity, "axial"): axial_total,
                total_col: total,
                "axial_fraction": axial_total / total if total > 0.0 else 0.0,
                "hydrogen_axial_warning": "; ".join(
                    sorted(
                        {
                            str(value)
                            for value in group.get("hydrogen_axial_warning", pd.Series(dtype=str))
                            if str(value)
                        }
                    )
                ),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def normal_cdf(x: np.ndarray) -> np.ndarray:
    vec_erf = np.vectorize(erf)
    return 0.5 * (1.0 + vec_erf(np.asarray(x, dtype=float) / np.sqrt(2.0)))


def gaussian_bin_probabilities(
    means: np.ndarray,
    sigmas: np.ndarray,
    bin_edges: np.ndarray,
) -> np.ndarray:
    mu = np.asarray(means, dtype=float).reshape(-1, 1)
    sigma = np.clip(np.asarray(sigmas, dtype=float).reshape(-1, 1), 1.0e-12, None)
    edges = np.asarray(bin_edges, dtype=float)
    return np.clip(
        normal_cdf((edges[1:].reshape(1, -1) - mu) / sigma)
        - normal_cdf((edges[:-1].reshape(1, -1) - mu) / sigma),
        0.0,
        1.0,
    )


def cevns_sigma_T_keV(T_keV: np.ndarray, config: CevnsConfig) -> np.ndarray:
    if not config.has_energy_resolution:
        return np.zeros_like(T_keV, dtype=float)
    T = np.clip(np.asarray(T_keV, dtype=float), 1.0e-9, None)
    threshold_term = config.energy_resolution_at_threshold_frac * np.sqrt(
        config.nr_threshold_keV / T
    )
    frac = np.sqrt(threshold_term**2 + config.energy_resolution_constant_frac**2)
    return frac * T


def d_enu_min_dT_nuclear(T_keV: np.ndarray, target: cevns.NuclearTarget) -> np.ndarray:
    T = np.clip(np.asarray(T_keV, dtype=float), 1.0e-12, None)
    m_keV = target.mass_MeV * 1.0e3
    root = np.sqrt(T**2 + 2.0 * m_keV * T)
    return 0.5 * (1.0 + (T + m_keV) / root)


def build_cevns_lower_bound_reco(
    recoil_df: pd.DataFrame,
    config: CevnsConfig,
    quantity: str,
) -> pd.DataFrame:
    """Map dR/dT_N or dN/dT_N into E_nu_min(T_N) bins."""

    if recoil_df.empty:
        return pd.DataFrame()

    positive_emin = []
    for isotope_key, group in recoil_df.groupby("isotope_key"):
        target = cevns.NUCLEAR_TARGETS[str(isotope_key)]
        T = group["T_N_keV"].to_numpy(dtype=float)
        positive_emin.append(cevns.emin_from_recoil_keV(T, target))
    all_emin = np.concatenate(positive_emin)
    all_emin = all_emin[np.isfinite(all_emin) & (all_emin > 0.0)]
    if len(all_emin) == 0:
        return pd.DataFrame()

    e_min = max(float(all_emin.min()) * 0.8, 1.0e-4)
    e_max = max(float(all_emin.max()) * 1.2, e_min * 1.01)
    edges = np.logspace(np.log10(e_min), np.log10(e_max), config.reco_energy_bins + 1)
    centers = np.sqrt(edges[:-1] * edges[1:])
    widths = np.diff(edges)

    rows = []
    gas_label = str(recoil_df["gas_label"].iloc[0])
    gas_name = str(recoil_df["gas"].iloc[0])
    for isotope_key, group in recoil_df.groupby("isotope_key", sort=False):
        target = cevns.NUCLEAR_TARGETS[str(isotope_key)]
        T = group["T_N_keV"].to_numpy(dtype=float)
        dT = bin_widths_from_centers(T)
        means = cevns.emin_from_recoil_keV(T, target)
        sigma_mev = (
            d_enu_min_dT_nuclear(T, target) * cevns_sigma_T_keV(T, config) / 1.0e3
            if config.has_energy_resolution
            else np.zeros_like(T)
        )

        histograms = {}
        for component in ("vector", "axial", "total"):
            density_col = component_column(quantity, "recoil", component)
            weights = group[density_col].to_numpy(dtype=float) * dT
            if config.has_energy_resolution:
                probs = gaussian_bin_probabilities(means, sigma_mev, edges)
                hist = np.sum(weights.reshape(-1, 1) * probs, axis=0)
            else:
                hist = np.histogram(means, bins=edges, weights=weights)[0]
            histograms[component] = hist

        for idx, center in enumerate(centers):
            row = {
                "gas_label": gas_label,
                "gas": gas_name,
                "isotope": target.label,
                "isotope_key": isotope_key,
                "A": target.A,
                "Z": target.Z,
                "N": target.N,
                "E_nu_reco_min_MeV": center,
                "E_nu_bin_low_MeV": edges[idx],
                "E_nu_bin_high_MeV": edges[idx + 1],
                "bin_width_MeV": widths[idx],
                "reconstruction_note": config.reconstruction_note,
            }
            for component in ("vector", "axial", "total"):
                row[reco_component_column(quantity, component)] = histograms[component][idx]
                row[reco_component_column(quantity, component, density=True)] = (
                    histograms[component][idx] / widths[idx] if widths[idx] > 0.0 else 0.0
                )
            rows.append(row)

    return pd.DataFrame(rows)


def ensure_output_dir(path: Path, skip_save: bool, skip_plots: bool) -> None:
    if not (skip_save and skip_plots):
        path.mkdir(parents=True, exist_ok=True)


def save_cevns_tables(
    gas_outdir: Path,
    recoil_df: pd.DataFrame,
    enu_df: pd.DataFrame,
    reco_df: pd.DataFrame | None,
    skip_save: bool,
) -> None:
    if skip_save:
        return
    gas_outdir.mkdir(parents=True, exist_ok=True)
    quantity = "rate" if any(col.startswith("dR_") for col in recoil_df.columns) else "count"

    aggregate_cevns_recoil_spectrum(recoil_df, quantity).to_csv(
        gas_outdir / "cevns_recoil_spectrum.csv",
        index=False,
    )
    aggregate_cevns_enu_spectrum(enu_df, quantity).to_csv(
        gas_outdir / "cevns_accepted_Enu_spectrum.csv",
        index=False,
    )
    recoil_df.to_csv(gas_outdir / "cevns_recoil_spectrum_by_isotope.csv", index=False)
    enu_df.to_csv(gas_outdir / "cevns_accepted_Enu_spectrum_by_isotope.csv", index=False)
    if reco_df is not None:
        aggregate_cevns_reco_spectrum(reco_df, quantity).to_csv(
            gas_outdir / "cevns_reco_energy_min_spectrum.csv",
            index=False,
        )
        reco_df.to_csv(gas_outdir / "cevns_reco_energy_min_spectrum_by_isotope.csv", index=False)


def plot_positive_series(ax, x: Iterable[float], y: Iterable[float], label: str, **kwargs) -> None:
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    mask = np.isfinite(x_arr) & np.isfinite(y_arr) & (x_arr > 0.0) & (y_arr > 0.0)
    if np.any(mask):
        ax.plot(x_arr[mask], y_arr[mask], label=label, **kwargs)


def write_cevns_spectrum_plots(
    gas_outdir: Path,
    gas_label: str,
    recoil_df: pd.DataFrame,
    enu_df: pd.DataFrame,
    quantity: str,
    skip_plots: bool,
) -> None:
    if skip_plots:
        return
    gas_outdir.mkdir(parents=True, exist_ok=True)
    aggregate_recoil_df = aggregate_cevns_recoil_spectrum(recoil_df, quantity)
    aggregate_enu_df = aggregate_cevns_enu_spectrum(enu_df, quantity)
    recoil_total_col = component_column(quantity, "recoil", "total")
    enu_total_col = component_column(quantity, "enu", "total")
    recoil_ylabel = (
        r"$dR/dT_N$ [s$^{-1}$ keV$^{-1}$]"
        if quantity == "rate"
        else r"$dN/dT_N$ [counts keV$^{-1}$]"
    )
    enu_ylabel = (
        r"accepted $dR/dE_\nu$ [s$^{-1}$ MeV$^{-1}$]"
        if quantity == "rate"
        else r"accepted $dN/dE_\nu$ [counts MeV$^{-1}$]"
    )

    fig, ax = plt.subplots(figsize=(8.3, 5.5))
    for component, style in (("total", "-"), ("vector", "--"), ("axial", ":")):
        col = component_column(quantity, "recoil", component)
        plot_positive_series(
            ax,
            aggregate_recoil_df["T_N_keV"],
            aggregate_recoil_df[col],
            component,
            linewidth=2.0,
            linestyle=style,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"True nuclear recoil $T_N$ [keV]")
    ax.set_ylabel(recoil_ylabel)
    ax.set_title(f"Gas-total CEvNS nuclear-recoil spectrum\n{gas_label}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(gas_outdir / "cevns_recoil_spectrum.png", dpi=PLOT_DPI)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.3, 5.5))
    for component, style in (("total", "-"), ("vector", "--"), ("axial", ":")):
        col = component_column(quantity, "enu", component)
        plot_positive_series(
            ax,
            aggregate_enu_df["E_nu_MeV"],
            aggregate_enu_df[col],
            component,
            linewidth=2.0,
            linestyle=style,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Neutrino energy $E_\nu$ [MeV]")
    ax.set_ylabel(enu_ylabel)
    ax.set_title(f"Gas-total CEvNS accepted neutrino-energy spectrum\n{gas_label}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(gas_outdir / "cevns_accepted_Enu_spectrum.png", dpi=PLOT_DPI)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.3, 5.5))
    for isotope, group in recoil_df.groupby("isotope", sort=False):
        plot_positive_series(ax, group["T_N_keV"], group[recoil_total_col], isotope, linewidth=2.0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"True nuclear recoil $T_N$ [keV]")
    ax.set_ylabel(recoil_ylabel)
    ax.set_title(f"CEvNS nuclear-recoil spectrum by isotope\n{gas_label}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(gas_outdir / "cevns_recoil_spectrum_by_isotope.png", dpi=PLOT_DPI)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.3, 5.5))
    for isotope, group in enu_df.groupby("isotope", sort=False):
        plot_positive_series(ax, group["E_nu_MeV"], group[enu_total_col], isotope, linewidth=2.0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Neutrino energy $E_\nu$ [MeV]")
    ax.set_ylabel(enu_ylabel)
    ax.set_title(f"CEvNS accepted neutrino-energy spectrum by isotope\n{gas_label}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(gas_outdir / "cevns_accepted_Enu_spectrum_by_isotope.png", dpi=PLOT_DPI)
    plt.close(fig)


def write_cevns_reco_plot(
    gas_outdir: Path,
    gas_label: str,
    reco_df: pd.DataFrame,
    quantity: str,
    skip_plots: bool,
) -> None:
    if skip_plots or reco_df is None or reco_df.empty:
        return
    gas_outdir.mkdir(parents=True, exist_ok=True)
    aggregate_reco_df = aggregate_cevns_reco_spectrum(reco_df, quantity)
    y_col = reco_component_column(quantity, "total", density=True)
    ylabel = (
        r"CEvNS lower-bound $dR/dE_{\nu,\min}$ [s$^{-1}$ MeV$^{-1}$]"
        if quantity == "rate"
        else r"CEvNS lower-bound $dN/dE_{\nu,\min}$ [counts MeV$^{-1}$]"
    )
    fig, ax = plt.subplots(figsize=(8.3, 5.5))
    for component, style in (("total", "-"), ("vector", "--"), ("axial", ":")):
        col = reco_component_column(quantity, component, density=True)
        plot_positive_series(
            ax,
            aggregate_reco_df["E_nu_reco_min_MeV"],
            aggregate_reco_df[col],
            component,
            linewidth=2.0,
            linestyle=style,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$E_{\nu,\min}(T_N)$ [MeV]")
    ax.set_ylabel(ylabel)
    ax.set_title(f"Gas-total CEvNS energy-only lower-bound spectrum\n{gas_label}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(gas_outdir / "cevns_reco_energy_min_spectrum.png", dpi=PLOT_DPI)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.3, 5.5))
    for isotope, group in reco_df.groupby("isotope", sort=False):
        plot_positive_series(ax, group["E_nu_reco_min_MeV"], group[y_col], isotope, linewidth=2.0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$E_{\nu,\min}(T_N)$ [MeV]")
    ax.set_ylabel(ylabel)
    ax.set_title(f"CEvNS energy-only lower-bound spectrum by isotope\n{gas_label}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(gas_outdir / "cevns_reco_energy_min_spectrum_by_isotope.png", dpi=PLOT_DPI)
    plt.close(fig)


def write_cevns_summary(
    cevns_outdir: Path,
    summary_rows: list[pd.DataFrame],
    quantity: str,
    skip_save: bool,
    skip_plots: bool,
    filename: str,
) -> pd.DataFrame:
    if summary_rows:
        summary_df = pd.concat(summary_rows, ignore_index=True)
    else:
        summary_df = pd.DataFrame()
    if summary_df.empty:
        return summary_df

    aggregate_df = aggregate_cevns_summary(summary_df, quantity)
    detail_filename = filename.replace(".csv", "_by_isotope.csv")
    if not aggregate_df.empty:
        aggregate_df["isotope_detail_csv"] = detail_filename if not skip_save else ""

    if not skip_save:
        cevns_outdir.mkdir(parents=True, exist_ok=True)
        aggregate_df.to_csv(cevns_outdir / filename, index=False)
        summary_df.to_csv(cevns_outdir / detail_filename, index=False)

    if not skip_plots:
        cevns_outdir.mkdir(parents=True, exist_ok=True)
        value_col = summary_column(quantity, "total")
        labels = aggregate_df["gas_label"].copy()
        if "model_name" in aggregate_df.columns:
            labels = aggregate_df["model_name"].astype(str) + "\n" + labels
        values = aggregate_df[value_col].to_numpy(dtype=float)
        order = np.argsort(values)
        fig, ax = plt.subplots(figsize=(11, max(5.5, 0.32 * len(aggregate_df))))
        ax.barh(labels.iloc[order], values[order], color="#2563eb")
        positive = values[values > 0.0]
        if len(positive) > 0 and positive.max() / max(positive.min(), 1.0e-300) > 30.0:
            ax.set_xscale("log")
        ax.set_xlabel("Total CEvNS rate [s$^{-1}$]" if quantity == "rate" else "Total CEvNS counts")
        ax.set_title("Gas-total CEvNS summary")
        ax.grid(True, axis="x", which="both", alpha=0.3)
        fig.tight_layout()
        fig.savefig(cevns_outdir / filename.replace(".csv", ".png"), dpi=PLOT_DPI)
        plt.close(fig)
    return aggregate_df
