#!/usr/bin/env python3
"""
Distance scan for supernova gas targets.

The input supernova spectra are time-integrated fluences, so the primary
outputs are expected counts. Average rates are reported as counts divided by
the requested burst duration. By default every gas row in the gas-density table
is scanned, with one output subfolder per gas row. Distance dependence is
applied explicitly as (fluence_reference_distance / distance)^2 after computing
the spectra once per model and gas.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import scatteringPlots as sp
from cevns_pipeline import (
    CEVNS_AXIAL_MODELS,
    CEVNS_FORM_FACTORS,
    CevnsConfig,
    aggregate_cevns_enu_spectrum,
    aggregate_cevns_recoil_spectrum,
    component_column,
    compute_cevns_spectra_for_gas,
    supernova_active_fluence,
    validate_cevns_config,
)
from detector_model import (
    DEFAULT_DIFFUSION_CSV,
    DEFAULT_GAS_CSV,
    DEFAULT_GEOMETRY_CONFIG,
    DEFAULT_RANGE_CSV,
    PLOT_DPI,
    REPO_ROOT,
    electron_density_cm3,
    gas_entry_label,
    gas_entry_slug,
    isotope_target_counts,
    read_detector_geometry_config,
    read_diffusion_summary_table,
    read_gas_density_table,
    read_recoil_window_table,
    resolve_recoil_window_keV,
    safe_slug,
)
from gasTargetRates import sigma_accepted_rescaled
from supernovaGasTargetEvents import (
    DEFAULT_SUPERNOVA_MODEL,
    FLAVORS,
    OUTPUT_SUFFIX,
    choose_models_for_args,
    convolved_recoil_energy_counts_per_electron,
    default_fluence_root,
    discover_fluence_models,
    normalize_grid_spectrum,
    read_supernova_fluence_csv,
)


DEFAULT_OUTPUT_ROOT = REPO_ROOT / "supernova_distance_gas_scan"
DEFAULT_FLUENCE_REFERENCE_DISTANCE_KPC = 10.0
TRAPEZOID = getattr(np, "trapezoid", np.trapz)


def parse_optional_float(value: str) -> float | None:
    text = str(value).strip().lower()
    if text in {"", "none", "null"}:
        return None
    return float(value)


def parse_distance_list(text: str) -> list[float]:
    values = [float(item.strip()) for item in text.split(",") if item.strip()]
    if not values:
        raise argparse.ArgumentTypeError("distance list cannot be empty")
    if any(value <= 0.0 for value in values):
        raise argparse.ArgumentTypeError("plot distances must be positive")
    return values


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute supernova nu-e and CEvNS counts and spectra from 0.5 to "
            "10 kpc by default. All configured gas rows are scanned unless "
            "--gas is provided."
        )
    )
    parser.add_argument(
        "--fluence-root",
        default=str(default_fluence_root()),
        help="Root folder or CSV containing supernova fluence inputs",
    )
    parser.add_argument(
        "--models",
        default=DEFAULT_SUPERNOVA_MODEL,
        help="Comma-separated model names to process, or 'all'",
    )
    parser.add_argument(
        "--gas-density-csv",
        "--gas-csv",
        dest="gas_csv",
        default=str(DEFAULT_GAS_CSV),
        help="Gas density table",
    )
    parser.add_argument(
        "--gas",
        default=None,
        help=(
            "Optional gas-row filter. Matches the Gas column, full label, or "
            "slug. Omit to process every row in the gas table."
        ),
    )
    parser.add_argument(
        "--pressure-atm",
        type=float,
        default=None,
        help=(
            "Pressure selector when --gas matches multiple pressure rows. "
            "Without --gas, this filters the all-gas scan to this pressure."
        ),
    )
    parser.add_argument("--list-gases", action="store_true", help="Print available gas rows and exit")
    parser.add_argument(
        "--range-csv",
        default=str(DEFAULT_RANGE_CSV),
        help="Electron recoil-window table used for nu-e scattering",
    )
    parser.add_argument(
        "--diffusion-csv",
        default=str(DEFAULT_DIFFUSION_CSV),
        help="Diffusion summary used for nu-e lower threshold adjustment",
    )
    parser.add_argument(
        "--geometry",
        default=str(DEFAULT_GEOMETRY_CONFIG),
        help="Detector geometry JSON",
    )
    parser.add_argument(
        "--volume-cm3",
        type=float,
        default=None,
        help="Optional detector volume override in cm^3",
    )
    parser.add_argument(
        "--output-root",
        default=str(DEFAULT_OUTPUT_ROOT),
        help="Output directory",
    )
    parser.add_argument(
        "--distance-min-kpc",
        type=float,
        default=0.5,
        help="Minimum supernova distance in kpc",
    )
    parser.add_argument(
        "--distance-max-kpc",
        type=float,
        default=10.0,
        help="Maximum supernova distance in kpc",
    )
    parser.add_argument(
        "--distance-step-kpc",
        type=float,
        default=0.5,
        help="Distance step in kpc",
    )
    parser.add_argument(
        "--fluence-reference-distance-kpc",
        type=float,
        default=DEFAULT_FLUENCE_REFERENCE_DISTANCE_KPC,
        help="Distance assumed by the input fluence CSV",
    )
    parser.add_argument(
        "--burst-duration-s",
        type=float,
        default=10.0,
        help="Duration used only to convert counts to average rates",
    )
    parser.add_argument(
        "--electron-grid-points",
        type=int,
        default=600,
        help="Electron-recoil grid points for nu-e spectra",
    )
    parser.add_argument(
        "--er-grid-points",
        type=int,
        default=600,
        help="Nuclear-recoil grid points for CEvNS spectra",
    )
    parser.add_argument(
        "--max-energy-points",
        type=int,
        default=None,
        help="Optional reduced neutrino-energy grid for smoke tests",
    )
    parser.add_argument(
        "--nr-threshold-kev",
        type=float,
        default=1.0,
        help="True nuclear-recoil CEvNS threshold in keV",
    )
    parser.add_argument(
        "--nr-max-kev",
        type=parse_optional_float,
        default=None,
        help="Optional CEvNS true nuclear-recoil upper cut in keV; use 'none' for kinematic max",
    )
    parser.add_argument(
        "--form-factor",
        choices=CEVNS_FORM_FACTORS,
        default="default",
        help="CEvNS vector form-factor mode",
    )
    parser.add_argument(
        "--axial-model",
        choices=CEVNS_AXIAL_MODELS,
        default="hoferichter_19f_fast",
        help="CEvNS axial model",
    )
    parser.add_argument(
        "--plot-distances-kpc",
        type=parse_distance_list,
        default=parse_distance_list("0.5,1,2,5,10"),
        help="Comma-separated distances to overlay in spectrum plots",
    )
    parser.add_argument("--skip-plots", action="store_true", help="Write CSVs only")
    return parser.parse_args()


def distance_grid_kpc(start: float, stop: float, step: float) -> np.ndarray:
    if start <= 0.0 or stop <= 0.0 or step <= 0.0:
        raise ValueError("Distance min, max, and step must be positive.")
    if stop < start:
        raise ValueError("--distance-max-kpc must be >= --distance-min-kpc.")
    n_steps = int(np.floor((stop - start) / step + 1.0e-9))
    grid = start + step * np.arange(n_steps + 1, dtype=float)
    if grid[-1] < stop - 1.0e-9:
        grid = np.append(grid, stop)
    return np.round(grid, 10)


def print_available_gases(gas_df: pd.DataFrame) -> None:
    print("Available gas rows:")
    for idx, row in gas_df.reset_index(drop=True).iterrows():
        gas_name = str(row["Gas"])
        pressure = float(row["Pressure (atm)"])
        density = float(row["density @ 293 K (g/cm3)"])
        label = gas_entry_label(gas_name, pressure)
        slug = gas_entry_slug(gas_name, pressure)
        print(f"  [{idx}] {label} | density={density:.8g} g/cm^3 | slug={slug}")


def select_gas_rows(
    gas_df: pd.DataFrame,
    gas_query: str | None,
    pressure_atm: float | None,
) -> list[pd.Series]:
    if not gas_query:
        matches = []
        for _, row in gas_df.iterrows():
            pressure = float(row["Pressure (atm)"])
            if pressure_atm is None or np.isclose(pressure, pressure_atm, rtol=0.0, atol=1.0e-9):
                matches.append(row)
        if not matches:
            pressure_text = "any" if pressure_atm is None else f"{pressure_atm:g}"
            raise ValueError(
                f"No gas rows matched --pressure-atm {pressure_text}. "
                "Run --list-gases for choices."
            )
        return matches

    query_slug = safe_slug(gas_query)
    matches = []
    for _, row in gas_df.iterrows():
        gas_name = str(row["Gas"])
        pressure = float(row["Pressure (atm)"])
        candidates = {
            gas_name.strip().lower(),
            gas_entry_label(gas_name, pressure).strip().lower(),
            gas_entry_slug(gas_name, pressure),
            safe_slug(gas_name),
        }
        if gas_query.strip().lower() in candidates or query_slug in candidates:
            if pressure_atm is None or np.isclose(pressure, pressure_atm, rtol=0.0, atol=1.0e-9):
                matches.append(row)

    if not matches:
        raise ValueError(f"No gas row matched --gas {gas_query!r}. Run --list-gases for choices.")
    if len(matches) > 1:
        labels = [
            gas_entry_label(str(row["Gas"]), float(row["Pressure (atm)"]))
            for row in matches
        ]
        raise ValueError(
            "Gas selection matched multiple pressure rows: "
            + "; ".join(labels)
            + ". Add --pressure-atm."
        )
    return matches


def thin_fluence_grid(fluence_df: pd.DataFrame, max_energy_points: int | None) -> pd.DataFrame:
    if max_energy_points is None or len(fluence_df) <= max_energy_points:
        return fluence_df.reset_index(drop=True)
    if max_energy_points < 2:
        raise ValueError("--max-energy-points must be at least 2 when provided.")
    indices = np.unique(np.linspace(0, len(fluence_df) - 1, max_energy_points).astype(int))
    return fluence_df.iloc[indices].reset_index(drop=True)


def scale_factor(reference_distance_kpc: float, distance_kpc: float) -> float:
    return (reference_distance_kpc / distance_kpc) ** 2


def compute_electron_spectra_for_gas(
    fluence_df: pd.DataFrame,
    target_electrons: float,
    T_low_MeV: float,
    T_high_MeV: float,
    t_grid_points: int,
) -> dict[str, np.ndarray]:
    E = fluence_df["E_MeV"].to_numpy(dtype=float)
    dE = sp.bin_widths_from_centers(E)
    T_grid = np.linspace(T_low_MeV, T_high_MeV, t_grid_points)
    spectra: dict[str, np.ndarray] = {
        "E_MeV": E,
        "dE_MeV": dE,
        "Te_MeV": T_grid,
    }

    for flavor in FLAVORS:
        fluence = fluence_df[flavor.key].to_numpy(dtype=float)
        sigma_acc = sigma_accepted_rescaled(
            E,
            channel=flavor.channel,
            T_low_MeV=T_low_MeV,
            T_high_MeV=T_high_MeV,
        )
        dN_dEnu = target_electrons * fluence * sigma_acc
        count_bins = dN_dEnu * dE
        target_total = float(np.sum(count_bins))
        dN_dTe = target_electrons * convolved_recoil_energy_counts_per_electron(
            E,
            fluence,
            T_grid,
            channel=flavor.channel,
        )
        dN_dTe = normalize_grid_spectrum(T_grid, dN_dTe, target_total)

        suffix = OUTPUT_SUFFIX[flavor.key]
        spectra[f"dN_dEnu_{suffix}"] = dN_dEnu
        spectra[f"count_bin_{suffix}"] = count_bins
        spectra[f"dN_dTe_{suffix}"] = dN_dTe

    spectra["dN_dEnu_total"] = sum(spectra[f"dN_dEnu_{name}"] for name in OUTPUT_SUFFIX.values())
    spectra["count_bin_total"] = sum(
        spectra[f"count_bin_{name}"] for name in OUTPUT_SUFFIX.values()
    )
    spectra["dN_dTe_total"] = sum(spectra[f"dN_dTe_{name}"] for name in OUTPUT_SUFFIX.values())
    return spectra


def electron_totals(spectra: dict[str, np.ndarray]) -> dict[str, float]:
    out = {
        f"electron_counts_{name}": float(np.sum(spectra[f"count_bin_{name}"]))
        for name in OUTPUT_SUFFIX.values()
    }
    out["electron_counts_total"] = float(np.sum(spectra["count_bin_total"]))
    out["electron_counts_from_recoil_integral"] = float(
        TRAPEZOID(spectra["dN_dTe_total"], spectra["Te_MeV"])
    )
    return out


def scale_electron_recoil_spectra(
    spectra: dict[str, np.ndarray],
    distances_kpc: np.ndarray,
    reference_distance_kpc: float,
    burst_duration_s: float,
) -> pd.DataFrame:
    rows = []
    flavor_names = list(OUTPUT_SUFFIX.values())
    for distance in distances_kpc:
        scale = scale_factor(reference_distance_kpc, float(distance))
        for idx, T_MeV in enumerate(spectra["Te_MeV"]):
            row = {
                "distance_kpc": float(distance),
                "fluence_scale": scale,
                "Te_MeV": T_MeV,
                "Te_keV": T_MeV * 1.0e3,
            }
            for name in flavor_names + ["total"]:
                value = spectra[f"dN_dTe_{name}"][idx] * scale
                row[f"dN_dTe_{name}_MeV-1"] = value
                row[f"dN_dTe_{name}_keV-1"] = value / 1.0e3
                row[f"avg_dR_dTe_{name}_s-1_MeV-1"] = value / burst_duration_s
                row[f"avg_dR_dTe_{name}_s-1_keV-1"] = value / (burst_duration_s * 1.0e3)
            rows.append(row)
    return pd.DataFrame(rows)


def scale_electron_enu_spectra(
    spectra: dict[str, np.ndarray],
    distances_kpc: np.ndarray,
    reference_distance_kpc: float,
    burst_duration_s: float,
) -> pd.DataFrame:
    rows = []
    flavor_names = list(OUTPUT_SUFFIX.values())
    for distance in distances_kpc:
        scale = scale_factor(reference_distance_kpc, float(distance))
        for idx, energy in enumerate(spectra["E_MeV"]):
            row = {
                "distance_kpc": float(distance),
                "fluence_scale": scale,
                "E_nu_MeV": energy,
            }
            for name in flavor_names + ["total"]:
                value = spectra[f"dN_dEnu_{name}"][idx] * scale
                row[f"dN_dEnu_{name}_MeV-1"] = value
                row[f"avg_dR_dEnu_{name}_s-1_MeV-1"] = value / burst_duration_s
            rows.append(row)
    return pd.DataFrame(rows)


def scale_cevns_spectrum(
    base_df: pd.DataFrame,
    distances_kpc: np.ndarray,
    reference_distance_kpc: float,
    burst_duration_s: float,
    spectrum: str,
) -> pd.DataFrame:
    if base_df.empty:
        return base_df.copy()

    out_frames = []
    count_cols = [component_column("count", spectrum, comp) for comp in ("vector", "axial", "total")]
    rate_cols = [
        (
            f"avg_dR_dT_s-1_keV-1_{comp}"
            if spectrum == "recoil"
            else f"avg_accepted_rate_s-1_MeV-1_{comp}"
        )
        for comp in ("vector", "axial", "total")
    ]
    for distance in distances_kpc:
        scale = scale_factor(reference_distance_kpc, float(distance))
        frame = base_df.copy()
        frame.insert(0, "fluence_scale", scale)
        frame.insert(0, "distance_kpc", float(distance))
        for count_col, rate_col in zip(count_cols, rate_cols):
            frame[count_col] = frame[count_col] * scale
            frame[rate_col] = frame[count_col] / burst_duration_s
        out_frames.append(frame)
    return pd.concat(out_frames, ignore_index=True)


def scale_cevns_summary(
    summary_df: pd.DataFrame,
    distances_kpc: np.ndarray,
    reference_distance_kpc: float,
    burst_duration_s: float,
) -> pd.DataFrame:
    rows = []
    for distance in distances_kpc:
        scale = scale_factor(reference_distance_kpc, float(distance))
        for _, base in summary_df.iterrows():
            row = base.to_dict()
            row["distance_kpc"] = float(distance)
            row["fluence_scale"] = scale
            for component in ("vector", "axial", "total"):
                count_col = f"{component}_counts"
                rate_col = f"{component}_average_rate_Hz"
                row[count_col] = float(base[count_col]) * scale
                row[rate_col] = row[count_col] / burst_duration_s
            total = row["total_counts"]
            row["axial_fraction"] = row["axial_counts"] / total if total > 0.0 else 0.0
            rows.append(row)
    return pd.DataFrame(rows)


def nearest_distances(grid: np.ndarray, requested: list[float]) -> list[float]:
    selected = []
    for value in requested:
        nearest = float(grid[np.argmin(np.abs(grid - value))])
        if nearest not in selected:
            selected.append(nearest)
    return selected


def plot_positive(ax, x, y, label: str, **kwargs) -> None:
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    mask = np.isfinite(x_arr) & np.isfinite(y_arr) & (x_arr > 0.0) & (y_arr > 0.0)
    if np.any(mask):
        ax.plot(x_arr[mask], y_arr[mask], label=label, **kwargs)


def write_plots(
    outdir: Path,
    gas_label: str,
    summary_df: pd.DataFrame,
    electron_recoil_df: pd.DataFrame,
    electron_enu_df: pd.DataFrame,
    cevns_recoil_df: pd.DataFrame,
    cevns_enu_df: pd.DataFrame,
    plot_distances: list[float],
) -> None:
    fig, ax = plt.subplots(figsize=(8, 5.4))
    plot_positive(
        ax,
        summary_df["distance_kpc"],
        summary_df["electron_counts_total"],
        r"$\nu$-e",
        linewidth=2.2,
    )
    plot_positive(
        ax,
        summary_df["distance_kpc"],
        summary_df["cevns_counts_total"],
        "CEvNS",
        linewidth=2.2,
    )
    ax.set_yscale("log")
    ax.set_xlabel("Supernova distance [kpc]")
    ax.set_ylabel("Expected counts")
    ax.set_title(f"Supernova counts vs distance\n{gas_label}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "counts_vs_distance.png", dpi=PLOT_DPI)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5.4))
    plot_positive(
        ax,
        summary_df["distance_kpc"],
        summary_df["electron_average_rate_Hz_total"],
        r"$\nu$-e",
        linewidth=2.2,
    )
    plot_positive(
        ax,
        summary_df["distance_kpc"],
        summary_df["cevns_average_rate_Hz_total"],
        "CEvNS",
        linewidth=2.2,
    )
    ax.set_yscale("log")
    ax.set_xlabel("Supernova distance [kpc]")
    ax.set_ylabel("Average rate [Hz]")
    ax.set_title(f"Average burst rate vs distance\n{gas_label}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "average_rate_vs_distance.png", dpi=PLOT_DPI)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5.4))
    for distance in plot_distances:
        group = electron_recoil_df[np.isclose(electron_recoil_df["distance_kpc"], distance)]
        plot_positive(
            ax,
            group["Te_keV"],
            group["dN_dTe_total_keV-1"],
            f"{distance:g} kpc",
            linewidth=1.8,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Electron recoil $T_e$ [keV]")
    ax.set_ylabel(r"$dN/dT_e$ [counts keV$^{-1}$]")
    ax.set_title(f"Supernova neutrino-electron recoil spectra\n{gas_label}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "electron_recoil_spectra_selected_distances.png", dpi=PLOT_DPI)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5.4))
    total_col = component_column("count", "recoil", "total")
    for distance in plot_distances:
        group = cevns_recoil_df[np.isclose(cevns_recoil_df["distance_kpc"], distance)]
        plot_positive(
            ax,
            group["T_N_keV"],
            group[total_col],
            f"{distance:g} kpc",
            linewidth=1.8,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Nuclear recoil $T_N$ [keV]")
    ax.set_ylabel(r"$dN/dT_N$ [counts keV$^{-1}$]")
    ax.set_title(f"Supernova CEvNS recoil spectra\n{gas_label}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "cevns_recoil_spectra_selected_distances.png", dpi=PLOT_DPI)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5.4))
    for distance in plot_distances:
        e_group = electron_enu_df[np.isclose(electron_enu_df["distance_kpc"], distance)]
        c_group = cevns_enu_df[np.isclose(cevns_enu_df["distance_kpc"], distance)]
        plot_positive(
            ax,
            e_group["E_nu_MeV"],
            e_group["dN_dEnu_total_MeV-1"],
            rf"$\nu$-e {distance:g} kpc",
            linewidth=1.4,
            linestyle="-",
        )
        plot_positive(
            ax,
            c_group["E_nu_MeV"],
            c_group[component_column("count", "enu", "total")],
            f"CEvNS {distance:g} kpc",
            linewidth=1.4,
            linestyle="--",
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Neutrino energy $E_\nu$ [MeV]")
    ax.set_ylabel(r"Accepted $dN/dE_\nu$ [counts MeV$^{-1}$]")
    ax.set_title(f"Accepted neutrino-energy spectra\n{gas_label}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "accepted_enu_spectra_selected_distances.png", dpi=PLOT_DPI)
    plt.close(fig)


def write_metadata(
    outdir: Path,
    model_name: str,
    model_csv: Path,
    gas_name: str,
    pressure_atm: float,
    density_g_cm3: float,
    volume_cm3: float,
    detector_name: str,
    electron_threshold_low_keV: float,
    electron_threshold_high_keV: float,
    target_electrons: float,
    cevns_config: CevnsConfig,
    target_counts: dict[str, float],
) -> None:
    rows = [
        ("model_name", model_name),
        ("fluence_csv", str(model_csv)),
        ("gas", gas_name),
        ("pressure_atm", pressure_atm),
        ("density_g_cm3", density_g_cm3),
        ("detector_name", detector_name),
        ("volume_cm3", volume_cm3),
        ("target_electrons", target_electrons),
        ("electron_recoil_threshold_low_keV", electron_threshold_low_keV),
        ("electron_recoil_threshold_high_keV", electron_threshold_high_keV),
        ("cevns_nr_threshold_keV", cevns_config.nr_threshold_keV),
        ("cevns_nr_max_keV", cevns_config.nr_max_keV),
        ("cevns_form_factor", cevns_config.form_factor),
        ("cevns_axial_model", cevns_config.axial_model),
    ]
    for isotope, count in sorted(target_counts.items()):
        rows.append((f"target_count_{isotope}", count))
    pd.DataFrame(rows, columns=["key", "value"]).to_csv(outdir / "run_metadata.csv", index=False)


def process_model(
    model,
    fluence_df: pd.DataFrame,
    gas_row: pd.Series,
    volume_cm3: float,
    geometry_name: str,
    recoil_window,
    electron_grid_points: int,
    cevns_grid_points: int,
    cevns_config: CevnsConfig,
    distances_kpc: np.ndarray,
    reference_distance_kpc: float,
    burst_duration_s: float,
    output_root: Path,
    plot_distances: list[float],
    skip_plots: bool,
) -> pd.DataFrame:
    gas_name = str(gas_row["Gas"])
    density_g_cm3 = float(gas_row["density @ 293 K (g/cm3)"])
    pressure_atm = float(gas_row["Pressure (atm)"])
    gas_label = gas_entry_label(gas_name, pressure_atm)
    gas_slug = gas_entry_slug(gas_name, pressure_atm)
    model_outdir = output_root / model.name / gas_slug
    model_outdir.mkdir(parents=True, exist_ok=True)

    electron_density, molar_mass, electrons_per_mixture = electron_density_cm3(
        gas_name,
        density_g_cm3,
    )
    target_electrons = electron_density * volume_cm3
    target_counts, cevns_molar_mass, isotopes_per_mixture = isotope_target_counts(
        gas_name,
        density_g_cm3,
        volume_cm3,
    )

    electron_spectra = compute_electron_spectra_for_gas(
        fluence_df,
        target_electrons=target_electrons,
        T_low_MeV=recoil_window.low_keV / 1.0e3,
        T_high_MeV=recoil_window.high_keV / 1.0e3,
        t_grid_points=electron_grid_points,
    )
    active_fluence = supernova_active_fluence(fluence_df)
    cevns_recoil_by_isotope, cevns_enu_by_isotope, cevns_summary_by_isotope = (
        compute_cevns_spectra_for_gas(
            gas_name=gas_name,
            gas_label=gas_label,
            density_g_cm3=density_g_cm3,
            volume_cm3=volume_cm3,
            E_MeV=fluence_df["E_MeV"].to_numpy(dtype=float),
            active_flux_or_fluence=active_fluence,
            config=cevns_config,
            t_grid_points=cevns_grid_points,
            quantity="count",
        )
    )
    cevns_recoil_total = aggregate_cevns_recoil_spectrum(cevns_recoil_by_isotope, "count")
    cevns_enu_total = aggregate_cevns_enu_spectrum(cevns_enu_by_isotope, "count")

    electron_recoil_distance = scale_electron_recoil_spectra(
        electron_spectra,
        distances_kpc,
        reference_distance_kpc,
        burst_duration_s,
    )
    electron_enu_distance = scale_electron_enu_spectra(
        electron_spectra,
        distances_kpc,
        reference_distance_kpc,
        burst_duration_s,
    )
    cevns_recoil_distance = scale_cevns_spectrum(
        cevns_recoil_total,
        distances_kpc,
        reference_distance_kpc,
        burst_duration_s,
        spectrum="recoil",
    )
    cevns_enu_distance = scale_cevns_spectrum(
        cevns_enu_total,
        distances_kpc,
        reference_distance_kpc,
        burst_duration_s,
        spectrum="enu",
    )
    cevns_recoil_isotope_distance = scale_cevns_spectrum(
        cevns_recoil_by_isotope,
        distances_kpc,
        reference_distance_kpc,
        burst_duration_s,
        spectrum="recoil",
    )
    cevns_enu_isotope_distance = scale_cevns_spectrum(
        cevns_enu_by_isotope,
        distances_kpc,
        reference_distance_kpc,
        burst_duration_s,
        spectrum="enu",
    )
    cevns_summary_distance = scale_cevns_summary(
        cevns_summary_by_isotope,
        distances_kpc,
        reference_distance_kpc,
        burst_duration_s,
    )

    base_electron_totals = electron_totals(electron_spectra)
    cevns_base_vector = float(cevns_summary_by_isotope["vector_counts"].sum())
    cevns_base_axial = float(cevns_summary_by_isotope["axial_counts"].sum())
    cevns_base_total = float(cevns_summary_by_isotope["total_counts"].sum())

    summary_rows = []
    for distance in distances_kpc:
        scale = scale_factor(reference_distance_kpc, float(distance))
        electron_total = base_electron_totals["electron_counts_total"] * scale
        cevns_total = cevns_base_total * scale
        row = {
            "model_name": model.name,
            "fluence_csv": str(model.csv_path),
            "gas_label": gas_label,
            "gas": gas_name,
            "pressure_atm": pressure_atm,
            "density_g_cm3": density_g_cm3,
            "distance_kpc": float(distance),
            "fluence_reference_distance_kpc": reference_distance_kpc,
            "fluence_scale": scale,
            "burst_duration_s": burst_duration_s,
            "detector_volume_cm3": volume_cm3,
            "target_electrons": target_electrons,
            "electron_threshold_low_keV": recoil_window.low_keV,
            "electron_threshold_high_keV": recoil_window.high_keV,
            "electron_threshold_source": recoil_window.threshold_source,
            "electron_counts_total": electron_total,
            "electron_average_rate_Hz_total": electron_total / burst_duration_s,
            "electron_counts_from_recoil_integral": (
                base_electron_totals["electron_counts_from_recoil_integral"] * scale
            ),
            "cevns_counts_vector": cevns_base_vector * scale,
            "cevns_counts_axial": cevns_base_axial * scale,
            "cevns_counts_total": cevns_total,
            "cevns_average_rate_Hz_vector": cevns_base_vector * scale / burst_duration_s,
            "cevns_average_rate_Hz_axial": cevns_base_axial * scale / burst_duration_s,
            "cevns_average_rate_Hz_total": cevns_total / burst_duration_s,
            "cevns_to_electron_count_ratio": (
                cevns_total / electron_total if electron_total > 0.0 else np.nan
            ),
            "all_channels_counts_total": electron_total + cevns_total,
            "all_channels_average_rate_Hz": (electron_total + cevns_total) / burst_duration_s,
            "cevns_nr_threshold_keV": cevns_config.nr_threshold_keV,
            "cevns_nr_max_keV": cevns_config.nr_max_keV,
            "cevns_form_factor": cevns_config.form_factor,
            "cevns_axial_model": cevns_config.axial_model,
            "molar_mass_g_mol": molar_mass,
            "cevns_molar_mass_g_mol": cevns_molar_mass,
            "electrons_per_mixture": electrons_per_mixture,
        }
        for name in OUTPUT_SUFFIX.values():
            count = base_electron_totals[f"electron_counts_{name}"] * scale
            row[f"electron_counts_{name}"] = count
            row[f"electron_average_rate_Hz_{name}"] = count / burst_duration_s
        for isotope, count in sorted(target_counts.items()):
            row[f"number_of_targets_{isotope}"] = count
            row[f"isotopes_per_mixture_{isotope}"] = isotopes_per_mixture[isotope]
        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)

    summary_df.to_csv(model_outdir / "distance_summary.csv", index=False)
    electron_recoil_distance.to_csv(model_outdir / "electron_recoil_spectra_by_distance.csv", index=False)
    electron_enu_distance.to_csv(
        model_outdir / "electron_accepted_enu_spectra_by_distance.csv",
        index=False,
    )
    cevns_recoil_distance.to_csv(model_outdir / "cevns_recoil_spectra_by_distance.csv", index=False)
    cevns_enu_distance.to_csv(
        model_outdir / "cevns_accepted_enu_spectra_by_distance.csv",
        index=False,
    )
    cevns_recoil_isotope_distance.to_csv(
        model_outdir / "cevns_recoil_spectra_by_isotope_by_distance.csv",
        index=False,
    )
    cevns_enu_isotope_distance.to_csv(
        model_outdir / "cevns_accepted_enu_spectra_by_isotope_by_distance.csv",
        index=False,
    )
    cevns_summary_distance.to_csv(
        model_outdir / "cevns_summary_by_isotope_by_distance.csv",
        index=False,
    )
    write_metadata(
        model_outdir,
        model_name=model.name,
        model_csv=model.csv_path,
        gas_name=gas_name,
        pressure_atm=pressure_atm,
        density_g_cm3=density_g_cm3,
        volume_cm3=volume_cm3,
        detector_name=geometry_name,
        electron_threshold_low_keV=recoil_window.low_keV,
        electron_threshold_high_keV=recoil_window.high_keV,
        target_electrons=target_electrons,
        cevns_config=cevns_config,
        target_counts=target_counts,
    )

    if not skip_plots:
        plot_distances = nearest_distances(distances_kpc, plot_distances)
        write_plots(
            model_outdir,
            gas_label=gas_label,
            summary_df=summary_df,
            electron_recoil_df=electron_recoil_distance,
            electron_enu_df=electron_enu_distance,
            cevns_recoil_df=cevns_recoil_distance,
            cevns_enu_df=cevns_enu_distance,
            plot_distances=plot_distances,
        )

    print(f"Saved distance scan for {model.name}, {gas_label}: {model_outdir.resolve()}")
    return summary_df


def main() -> int:
    args = parse_args()
    if args.burst_duration_s <= 0.0:
        raise ValueError("--burst-duration-s must be positive.")
    if args.electron_grid_points < 20:
        raise ValueError("--electron-grid-points must be at least 20.")
    if args.er_grid_points < 20:
        raise ValueError("--er-grid-points must be at least 20.")

    gas_csv = Path(args.gas_csv)
    range_csv = Path(args.range_csv)
    diffusion_csv = Path(args.diffusion_csv)
    geometry_path = Path(args.geometry)
    fluence_root = Path(args.fluence_root)
    output_root = Path(args.output_root)

    gas_df = read_gas_density_table(gas_csv)
    if args.list_gases:
        print_available_gases(gas_df)
        return 0

    gas_rows = select_gas_rows(gas_df, args.gas, args.pressure_atm)

    recoil_window_df = read_recoil_window_table(range_csv)
    diffusion_df = read_diffusion_summary_table(diffusion_csv)
    geometry = read_detector_geometry_config(geometry_path)
    volume_cm3 = float(args.volume_cm3) if args.volume_cm3 is not None else geometry.volume_cm3

    cevns_config = CevnsConfig(
        enabled=True,
        nr_threshold_keV=float(args.nr_threshold_kev),
        nr_max_keV=args.nr_max_kev,
        form_factor=args.form_factor,
        axial_model=args.axial_model,
    )
    validate_cevns_config(cevns_config)

    distances_kpc = distance_grid_kpc(
        args.distance_min_kpc,
        args.distance_max_kpc,
        args.distance_step_kpc,
    )
    if args.fluence_reference_distance_kpc <= 0.0:
        raise ValueError("--fluence-reference-distance-kpc must be positive.")

    models = choose_models_for_args(
        discover_fluence_models(fluence_root),
        fluence_root=fluence_root,
        model_arg=args.models,
    )
    output_root.mkdir(parents=True, exist_ok=True)

    all_summary = []
    for model in models:
        fluence_df = read_supernova_fluence_csv(model.csv_path)
        fluence_df = thin_fluence_grid(fluence_df, args.max_energy_points)
        for gas_row in gas_rows:
            gas_name = str(gas_row["Gas"])
            density_g_cm3 = float(gas_row["density @ 293 K (g/cm3)"])
            pressure_atm = float(gas_row["Pressure (atm)"])
            recoil_window = resolve_recoil_window_keV(
                gas_name,
                density_g_cm3,
                pressure_atm,
                recoil_window_df=recoil_window_df,
                diffusion_df=diffusion_df,
            )
            summary_df = process_model(
                model=model,
                fluence_df=fluence_df,
                gas_row=gas_row,
                volume_cm3=volume_cm3,
                geometry_name=geometry.name,
                recoil_window=recoil_window,
                electron_grid_points=args.electron_grid_points,
                cevns_grid_points=args.er_grid_points,
                cevns_config=cevns_config,
                distances_kpc=distances_kpc,
                reference_distance_kpc=float(args.fluence_reference_distance_kpc),
                burst_duration_s=float(args.burst_duration_s),
                output_root=output_root,
                plot_distances=args.plot_distances_kpc,
                skip_plots=args.skip_plots,
            )
            all_summary.append(summary_df)

    all_summary_df = pd.concat(all_summary, ignore_index=True)
    all_summary_df.to_csv(output_root / "all_distance_summaries.csv", index=False)

    if len(gas_rows) == 1:
        only = gas_rows[0]
        gas_label = gas_entry_label(str(only["Gas"]), float(only["Pressure (atm)"]))
        print(f"Processed {len(models)} model(s) for one gas row: {gas_label}")
    else:
        print(f"Processed {len(models)} model(s) for {len(gas_rows)} gas rows")
    print(
        f"Distance grid: {distances_kpc[0]:g} to {distances_kpc[-1]:g} kpc "
        f"({len(distances_kpc)} points)"
    )
    print("Primary output is counts; average rates use burst_duration_s.")
    print(f"Root summary: {(output_root / 'all_distance_summaries.csv').resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
