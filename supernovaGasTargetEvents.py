#!/usr/bin/env python3
"""
Gas-target event-count spectra for supernova neutrino fluence inputs.

This is the supernova-fluence analogue of gasTargetRates.py.  The key unit
difference is that the input is a time-integrated fluence,

    Phi(E) [cm^-2 MeV^-1],

so the convolution with the nu-e cross section and the number of target
electrons gives expected counts, not rates.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from detector_model import (
    DEFAULT_DIFFUSION_CSV,
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
from gasTargetRates import (
    dsigma_dT_rescaled,
    sigma_accepted_rescaled,
)
from cevns_pipeline import (
    add_cevns_cli_args,
    compute_cevns_spectra_for_gas,
    load_cevns_config,
    save_cevns_tables,
    supernova_active_fluence,
    write_cevns_spectrum_plots,
    write_cevns_summary,
)


REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_SUPERNOVA_MODEL = "SN1987A_like"
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "supernova_nu_gas_target_events"
DEFAULT_RESPONSE_CONFIG = REPO_ROOT / "reco_response_config.json"
TRAPEZOID = getattr(np, "trapezoid", np.trapz)


@dataclass(frozen=True)
class FluenceModel:
    name: str
    csv_path: Path


@dataclass(frozen=True)
class FlavorSpec:
    key: str
    label: str
    channel: str


FLAVORS = (
    FlavorSpec("nue", r"$\nu_e$", "nu_e"),
    FlavorSpec("anue", r"$\bar{\nu}_e$", "anti_nu_e"),
    FlavorSpec("numu", r"$\nu_\mu$", "nu_x"),
    FlavorSpec("anumu", r"$\bar{\nu}_\mu$", "anti_nu_x"),
    FlavorSpec("nutau", r"$\nu_\tau$", "nu_x"),
    FlavorSpec("anutau", r"$\bar{\nu}_\tau$", "anti_nu_x"),
)

FLAVOR_COLUMNS = {
    "nue": "fluence_nue_cm-2_MeV-1",
    "anue": "fluence_anue_cm-2_MeV-1",
    "numu": "fluence_nux_single_cm-2_MeV-1",
    "anumu": "fluence_anux_single_cm-2_MeV-1",
    "nutau": "fluence_nux_single_cm-2_MeV-1",
    "anutau": "fluence_anux_single_cm-2_MeV-1",
}

GROUP_FALLBACK_COLUMNS = {
    "numu": "fluence_nux_group_numu_nutau_cm-2_MeV-1",
    "nutau": "fluence_nux_group_numu_nutau_cm-2_MeV-1",
    "anumu": "fluence_anux_group_anumu_anutau_cm-2_MeV-1",
    "anutau": "fluence_anux_group_anumu_anutau_cm-2_MeV-1",
}

OUTPUT_SUFFIX = {
    "nue": "nue",
    "anue": "anue",
    "numu": "numu",
    "anumu": "anumu",
    "nutau": "nutau",
    "anutau": "anutau",
}


def default_fluence_root() -> Path:
    candidates = [
        REPO_ROOT / "supernova_fluence",
        REPO_ROOT / "outputs" / "supernova_fluence",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute accepted supernova-neutrino nu-e event counts for the gas "
            "targets listed in DetectorNumbers/GasDensities.csv."
        )
    )
    parser.add_argument(
        "--fluence-root",
        default=str(default_fluence_root()),
        help="Root folder containing one fluence CSV per supernova model",
    )
    parser.add_argument(
        "--output-root",
        default=str(DEFAULT_OUTPUT_ROOT),
        help="Output directory for supernova gas-target event counts",
    )
    parser.add_argument(
        "--models",
        default=DEFAULT_SUPERNOVA_MODEL,
        help=(
            "Comma-separated model names to process, or 'all'. "
            f"Default: {DEFAULT_SUPERNOVA_MODEL}."
        ),
    )
    parser.add_argument("--gas-csv", default=str(DEFAULT_GAS_CSV), help="Gas density table")
    parser.add_argument(
        "--range-csv",
        default=str(DEFAULT_RANGE_CSV),
        help="Recoil-energy thresholds derived from electron ranges",
    )
    parser.add_argument(
        "--diffusion-csv",
        default=str(DEFAULT_DIFFUSION_CSV),
        help="Diffusion summary used to adjust the lower recoil threshold",
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
        "--volume-cm3",
        type=float,
        default=None,
        help="Optional fiducial detector-volume override in cm^3",
    )
    parser.add_argument(
        "--burst-duration-s",
        type=float,
        default=10.0,
        help="Optional burst duration used only for average-rate diagnostics",
    )
    parser.add_argument(
        "--t-grid-points",
        type=int,
        default=600,
        help="Number of recoil-energy grid points in the accepted window",
    )
    parser.add_argument(
        "--costh-grid-points",
        type=int,
        default=600,
        help="Number of recoil-direction grid points",
    )
    add_cevns_cli_args(parser)
    return parser.parse_args()


def csv_has_fluence_columns(path: Path) -> bool:
    try:
        header = pd.read_csv(path, nrows=0)
    except Exception:
        return False
    columns = {str(column).strip() for column in header.columns}
    required_base = {
        "E_MeV",
        "fluence_nue_cm-2_MeV-1",
        "fluence_anue_cm-2_MeV-1",
    }
    has_heavy_single = {
        "fluence_nux_single_cm-2_MeV-1",
        "fluence_anux_single_cm-2_MeV-1",
    }.issubset(columns)
    has_heavy_group = {
        "fluence_nux_group_numu_nutau_cm-2_MeV-1",
        "fluence_anux_group_anumu_anutau_cm-2_MeV-1",
    }.issubset(columns)
    return required_base.issubset(columns) and (has_heavy_single or has_heavy_group)


def discover_fluence_models(root: Path) -> list[FluenceModel]:
    if root.is_file():
        if not csv_has_fluence_columns(root):
            raise ValueError(f"Fluence CSV does not contain the expected columns: {root}")
        model_name = root.stem.removesuffix("_fluence")
        return [FluenceModel(model_name, root)]

    if not root.exists():
        raise FileNotFoundError(f"Supernova fluence root not found: {root}")

    models = []
    for path in sorted(root.rglob("*.csv")):
        if csv_has_fluence_columns(path):
            model_name = path.parent.name
            models.append(FluenceModel(model_name, path))

    if not models:
        raise FileNotFoundError(
            f"No supernova fluence CSVs with flavor-separated columns were found under: {root}"
        )
    return models


def select_models(models: list[FluenceModel], model_arg: str) -> list[FluenceModel]:
    if model_arg.strip().lower() == "all":
        return models

    wanted = {name.strip() for name in model_arg.split(",") if name.strip()}
    selected = [
        model
        for model in models
        if model.name in wanted or model.csv_path.stem.removesuffix("_fluence") in wanted
    ]
    found = {model.name for model in selected}
    missing = sorted(wanted - found)
    if missing:
        available = ", ".join(model.name for model in models)
        raise ValueError(
            "Requested supernova model(s) not found: "
            + ", ".join(missing)
            + f". Available models: {available}"
        )
    return selected


def choose_models_for_args(
    discovered_models: list[FluenceModel],
    fluence_root: Path,
    model_arg: str,
) -> list[FluenceModel]:
    if fluence_root.is_file() and model_arg == DEFAULT_SUPERNOVA_MODEL:
        return discovered_models
    return select_models(discovered_models, model_arg)


def read_supernova_fluence_csv(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path)
    raw.columns = raw.columns.str.strip()
    if "E_MeV" not in raw.columns:
        raise ValueError(f"Missing E_MeV in fluence CSV: {path}")

    df = pd.DataFrame()
    df["E_MeV"] = pd.to_numeric(raw["E_MeV"], errors="coerce")

    for flavor, column in FLAVOR_COLUMNS.items():
        if column in raw.columns:
            values = pd.to_numeric(raw[column], errors="coerce")
        elif flavor in GROUP_FALLBACK_COLUMNS and GROUP_FALLBACK_COLUMNS[flavor] in raw.columns:
            values = 0.5 * pd.to_numeric(raw[GROUP_FALLBACK_COLUMNS[flavor]], errors="coerce")
        else:
            raise ValueError(
                f"Missing fluence column for {flavor} in {path}. Expected '{column}'."
            )
        df[flavor] = values

    df = df.dropna().copy()
    df = df[df["E_MeV"] > 0.0]
    for flavor in OUTPUT_SUFFIX:
        df = df[df[flavor] >= 0.0]
    df = df.sort_values("E_MeV").reset_index(drop=True)

    if len(df) < 2:
        raise ValueError(f"Not enough valid energy points in fluence CSV: {path}")
    return df


def build_recoil_energy_grid(T_low_MeV: float, T_high_MeV: float, n_points: int) -> np.ndarray:
    if n_points < 2:
        raise ValueError("Need at least 2 recoil-energy grid points.")
    if T_low_MeV <= 0.0 or T_high_MeV <= T_low_MeV:
        raise ValueError("Invalid recoil-energy window.")
    return np.linspace(T_low_MeV, T_high_MeV, n_points)


def convolved_recoil_energy_counts_per_electron(
    E_MeV: np.ndarray,
    fluence: np.ndarray,
    T_MeV: np.ndarray,
    channel: str,
) -> np.ndarray:
    dE = sp.bin_widths_from_centers(E_MeV)
    counts = np.zeros_like(T_MeV, dtype=float)
    for E_i, fluence_i, dE_i in zip(E_MeV, fluence, dE):
        if fluence_i <= 0.0:
            continue
        counts += fluence_i * dsigma_dT_rescaled(E_i, T_MeV, channel) * dE_i
    return counts


def convolved_angular_counts_window_per_electron(
    E_MeV: np.ndarray,
    fluence: np.ndarray,
    costh: np.ndarray,
    channel: str,
    T_low_MeV: float,
    T_high_MeV: float,
) -> np.ndarray:
    dE = sp.bin_widths_from_centers(E_MeV)
    counts = np.zeros_like(costh, dtype=float)
    E_min_for_window = float(sp.Emin_from_T(np.array([T_low_MeV]))[0])

    for E_i, fluence_i, dE_i in zip(E_MeV, fluence, dE):
        if fluence_i <= 0.0 or E_i < E_min_for_window:
            continue
        dsdc = sp.dsigma_dcosth_rescaled(E_i, costh, channel=channel)
        recoil_t = sp.T_from_costh(E_i, costh)
        dsdc[(recoil_t < T_low_MeV) | (recoil_t > T_high_MeV)] = 0.0
        counts += fluence_i * dsdc * dE_i
    return counts


def normalize_grid_spectrum(
    x: np.ndarray,
    y: np.ndarray,
    target_integral: float,
) -> np.ndarray:
    out = np.asarray(y, dtype=float).copy()
    integral = float(TRAPEZOID(out, x))
    if target_integral > 0.0 and integral > 0.0:
        out *= target_integral / integral
    return out


def compute_event_spectra_for_gas(
    fluence_df: pd.DataFrame,
    target_electrons: float,
    T_low_MeV: float,
    T_high_MeV: float,
    t_grid_points: int,
    costh_grid_points: int,
) -> dict[str, np.ndarray]:
    E = fluence_df["E_MeV"].to_numpy(dtype=float)
    dE = sp.bin_widths_from_centers(E)
    T_grid = build_recoil_energy_grid(T_low_MeV, T_high_MeV, n_points=t_grid_points)
    costh_grid = np.linspace(0.0, 1.0, costh_grid_points)

    spectra: dict[str, np.ndarray] = {
        "E_MeV": E,
        "dE_MeV": dE,
        "Te_MeV": T_grid,
        "costh": costh_grid,
        "theta_deg": np.degrees(np.arccos(np.clip(costh_grid, 0.0, 1.0))),
    }

    unique_flavors = [FLAVORS[0], FLAVORS[1], FLAVORS[2], FLAVORS[3]]
    for flavor in unique_flavors:
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
        dN_dcosth = target_electrons * convolved_angular_counts_window_per_electron(
            E,
            fluence,
            costh_grid,
            channel=flavor.channel,
            T_low_MeV=T_low_MeV,
            T_high_MeV=T_high_MeV,
        )

        dN_dTe = normalize_grid_spectrum(T_grid, dN_dTe, target_total)
        dN_dcosth = normalize_grid_spectrum(costh_grid, dN_dcosth, target_total)

        suffix = OUTPUT_SUFFIX[flavor.key]
        spectra[f"dN_dEnu_{suffix}"] = dN_dEnu
        spectra[f"count_bin_{suffix}"] = count_bins
        spectra[f"dN_dTe_{suffix}"] = dN_dTe
        spectra[f"dN_dcosth_{suffix}"] = dN_dcosth

    for source, target in [("numu", "nutau"), ("anumu", "anutau")]:
        spectra[f"dN_dEnu_{target}"] = spectra[f"dN_dEnu_{source}"].copy()
        spectra[f"count_bin_{target}"] = spectra[f"count_bin_{source}"].copy()
        spectra[f"dN_dTe_{target}"] = spectra[f"dN_dTe_{source}"].copy()
        spectra[f"dN_dcosth_{target}"] = spectra[f"dN_dcosth_{source}"].copy()

    spectra["dN_dEnu_total"] = sum(spectra[f"dN_dEnu_{name}"] for name in OUTPUT_SUFFIX)
    spectra["count_bin_total"] = sum(spectra[f"count_bin_{name}"] for name in OUTPUT_SUFFIX)
    spectra["dN_dTe_total"] = sum(spectra[f"dN_dTe_{name}"] for name in OUTPUT_SUFFIX)
    spectra["dN_dcosth_total"] = sum(spectra[f"dN_dcosth_{name}"] for name in OUTPUT_SUFFIX)
    return spectra


def plot_positive_series(ax, x: np.ndarray, y: np.ndarray, label: str, **kwargs) -> None:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0.0) & (y > 0.0)
    if np.any(mask):
        ax.plot(x[mask], y[mask], label=label, **kwargs)


def write_event_tables(gas_outdir: Path, spectra: dict[str, np.ndarray]) -> None:
    dEnu_payload = {
        "E_MeV": spectra["E_MeV"],
        "dE_MeV": spectra["dE_MeV"],
    }
    for name in OUTPUT_SUFFIX.values():
        dEnu_payload[f"dN_dEnu_{name}_per_MeV"] = spectra[f"dN_dEnu_{name}"]
    dEnu_payload["dN_dEnu_total_per_MeV"] = spectra["dN_dEnu_total"]
    dEnu_payload["dN_bin_total"] = spectra["count_bin_total"]
    pd.DataFrame(dEnu_payload).to_csv(gas_outdir / "dN_dEnu.csv", index=False)

    accepted_mask = np.ones_like(spectra["Te_MeV"], dtype=bool)
    dTe_payload = {"Te_MeV": spectra["Te_MeV"]}
    for name in OUTPUT_SUFFIX.values():
        dTe_payload[f"dN_dTe_{name}_per_MeV"] = spectra[f"dN_dTe_{name}"]
    dTe_payload["dN_dTe_total_per_MeV"] = spectra["dN_dTe_total"]
    dTe_payload["accepted_window_mask"] = accepted_mask
    pd.DataFrame(dTe_payload).to_csv(gas_outdir / "dN_dTe.csv", index=False)

    total_costh_integral = float(TRAPEZOID(spectra["dN_dcosth_total"], spectra["costh"]))
    angular_pdf = np.zeros_like(spectra["dN_dcosth_total"])
    if total_costh_integral > 0.0:
        angular_pdf = spectra["dN_dcosth_total"] / total_costh_integral

    dcosth_payload = {
        "cos_theta_e": spectra["costh"],
        "theta_deg": spectra["theta_deg"],
    }
    for name in OUTPUT_SUFFIX.values():
        dcosth_payload[f"dN_dcosth_{name}"] = spectra[f"dN_dcosth_{name}"]
    dcosth_payload["dN_dcosth_total"] = spectra["dN_dcosth_total"]
    dcosth_payload["angular_pdf_total"] = angular_pdf
    pd.DataFrame(dcosth_payload).to_csv(gas_outdir / "dN_dcosth.csv", index=False)


def write_event_plots(
    gas_outdir: Path,
    model_name: str,
    gas_label_text: str,
    T_low_MeV: float,
    T_high_MeV: float,
    spectra: dict[str, np.ndarray],
) -> None:
    window_text = f"Accepted recoil window: {T_low_MeV:.4g} to {T_high_MeV:.4g} MeV"

    fig, ax = plt.subplots(figsize=(8, 5.5))
    for flavor in FLAVORS:
        suffix = OUTPUT_SUFFIX[flavor.key]
        plot_positive_series(
            ax,
            spectra["E_MeV"],
            spectra[f"dN_dEnu_{suffix}"],
            flavor.label,
            linewidth=1.6,
        )
    plot_positive_series(
        ax,
        spectra["E_MeV"],
        spectra["dN_dEnu_total"],
        "total",
        linewidth=2.4,
        linestyle="--",
        color="black",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Neutrino energy $E_\nu$ [MeV]")
    ax.set_ylabel(r"$dN/dE_\nu$ [counts MeV$^{-1}$]")
    ax.set_title(f"{model_name}: accepted event counts vs neutrino energy\n{gas_label_text}\n{window_text}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(ncol=2, fontsize=8.5)
    fig.tight_layout()
    fig.savefig(gas_outdir / "dN_dEnu.png", dpi=PLOT_DPI)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5.5))
    for flavor in FLAVORS:
        suffix = OUTPUT_SUFFIX[flavor.key]
        plot_positive_series(
            ax,
            spectra["Te_MeV"],
            spectra[f"dN_dTe_{suffix}"],
            flavor.label,
            linewidth=1.6,
        )
    plot_positive_series(
        ax,
        spectra["Te_MeV"],
        spectra["dN_dTe_total"],
        "total",
        linewidth=2.4,
        linestyle="--",
        color="black",
    )
    ax.axvline(T_low_MeV, color="black", linestyle=":", linewidth=1.1)
    ax.axvline(T_high_MeV, color="black", linestyle=":", linewidth=1.1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Recoil kinetic energy $T_e$ [MeV]")
    ax.set_ylabel(r"$dN/dT_e$ [counts MeV$^{-1}$]")
    ax.set_title(f"{model_name}: accepted recoil-electron spectrum\n{gas_label_text}\n{window_text}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(ncol=2, fontsize=8.5)
    fig.tight_layout()
    fig.savefig(gas_outdir / "dN_dTe.png", dpi=PLOT_DPI)
    plt.close(fig)

    total_integral = float(TRAPEZOID(spectra["dN_dcosth_total"], spectra["costh"]))
    pdf = (
        spectra["dN_dcosth_total"] / total_integral
        if total_integral > 0.0
        else np.zeros_like(spectra["dN_dcosth_total"])
    )
    fig, ax = plt.subplots(figsize=(8, 5.5))
    plot_positive_series(
        ax,
        spectra["costh"],
        spectra["dN_dcosth_total"],
        r"$dN/d\cos\theta_e$",
        linewidth=2.2,
        color="#2563eb",
    )
    ax.set_yscale("log")
    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel(r"$\cos\theta_e$")
    ax.set_ylabel(r"$dN/d\cos\theta_e$ [counts]")
    ax.grid(True, which="both", alpha=0.3)
    ax2 = ax.twinx()
    ax2.plot(spectra["costh"], pdf, color="#dc2626", linewidth=1.6, alpha=0.8, label="PDF")
    ax2.set_ylabel("normalized angular PDF")
    ax.set_title(f"{model_name}: accepted recoil-electron angular spectrum\n{gas_label_text}\n{window_text}")
    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines + lines2, labels + labels2, loc="best")
    fig.tight_layout()
    fig.savefig(gas_outdir / "dN_dcosth.png", dpi=PLOT_DPI)
    plt.close(fig)


def relative_difference(value: float, reference: float) -> float:
    scale = max(abs(reference), 1.0e-300)
    return abs(value - reference) / scale


def write_model_summary_plot(model_outdir: Path, summary_df: pd.DataFrame, model_name: str) -> None:
    ordered = summary_df.sort_values("total_events_all_flavors", ascending=True)
    labels = ordered["gas_label"].to_list()
    values = ordered["total_events_all_flavors"].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    ax.barh(labels, values, color="#2563eb")
    positive = values[values > 0.0]
    if len(positive) > 0 and positive.max() / positive.min() > 30.0:
        ax.set_xscale("log")
    ax.set_xlabel("Total expected events [counts]")
    ax.set_title(f"{model_name}: accepted supernova nu-e events by gas target")
    ax.grid(True, axis="x", which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(model_outdir / "gas_target_total_event_comparison.png", dpi=PLOT_DPI)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    if args.t_grid_points < 20:
        raise ValueError("--t-grid-points must be at least 20.")
    if args.costh_grid_points < 20:
        raise ValueError("--costh-grid-points must be at least 20.")
    if args.burst_duration_s <= 0.0:
        raise ValueError("--burst-duration-s must be positive.")

    fluence_root = Path(args.fluence_root)
    output_root = Path(args.output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    gas_csv = Path(args.gas_csv)
    range_csv = Path(args.range_csv)
    diffusion_csv = Path(args.diffusion_csv)
    geometry_config = Path(args.geometry_config)
    response_config = Path(args.response_config)
    for path, label in [
        (gas_csv, "Gas-density CSV"),
        (range_csv, "Recoil-window CSV"),
        (diffusion_csv, "Diffusion summary CSV"),
        (geometry_config, "Geometry config"),
    ]:
        if not path.exists():
            raise FileNotFoundError(f"{label} not found: {path}")

    models = choose_models_for_args(
        discover_fluence_models(fluence_root),
        fluence_root=fluence_root,
        model_arg=args.models,
    )
    gas_df = read_gas_density_table(gas_csv)
    recoil_window_df = read_recoil_window_table(range_csv)
    diffusion_df = read_diffusion_summary_table(diffusion_csv)
    geometry = read_detector_geometry_config(geometry_config)
    volume_cm3 = float(args.volume_cm3) if args.volume_cm3 is not None else geometry.volume_cm3
    cevns_config = load_cevns_config(response_config, args, default_enabled=False)

    all_summary_rows = []
    all_cevns_summary_rows = []
    progress = ProgressBar(len(models) * len(gas_df), prefix="Supernova gas targets")
    progress.update(0, "starting")
    progress_index = 0

    for model in models:
        fluence_df = read_supernova_fluence_csv(model.csv_path)
        cevns_fluence = supernova_active_fluence(fluence_df) if cevns_config.enabled else None
        model_outdir = output_root / model.name
        model_outdir.mkdir(parents=True, exist_ok=True)
        model_cevns_outdir = model_outdir / "cevns"
        model_summary_rows = []
        model_cevns_summary_rows = []

        for _, row in gas_df.iterrows():
            gas_name = row["Gas"]
            density_g_cm3 = float(row["density @ 293 K (g/cm3)"])
            pressure_atm = float(row["Pressure (atm)"])
            gas_label_text = gas_entry_label(gas_name, pressure_atm)
            gas_slug = gas_entry_slug(gas_name, pressure_atm)
            gas_outdir = model_outdir / gas_slug
            gas_outdir.mkdir(parents=True, exist_ok=True)
            progress.update(progress_index, f"{model.name}: {gas_label_text}")

            recoil_window = resolve_recoil_window_keV(
                gas_name,
                density_g_cm3,
                pressure_atm,
                recoil_window_df=recoil_window_df,
                diffusion_df=diffusion_df,
            )
            T_low_MeV = recoil_window.low_keV / 1.0e3
            T_high_MeV = recoil_window.high_keV / 1.0e3

            electron_density, molar_mass, electrons_per_mixture = electron_density_cm3(
                gas_name,
                density_g_cm3,
            )
            target_electrons = electron_density * volume_cm3

            spectra = compute_event_spectra_for_gas(
                fluence_df,
                target_electrons=target_electrons,
                T_low_MeV=T_low_MeV,
                T_high_MeV=T_high_MeV,
                t_grid_points=args.t_grid_points,
                costh_grid_points=args.costh_grid_points,
            )
            write_event_tables(gas_outdir, spectra)
            write_event_plots(
                gas_outdir,
                model_name=model.name,
                gas_label_text=gas_label_text,
                T_low_MeV=T_low_MeV,
                T_high_MeV=T_high_MeV,
                spectra=spectra,
            )

            if cevns_config.enabled:
                recoil_df, enu_df, cevns_summary_df = compute_cevns_spectra_for_gas(
                    gas_name=gas_name,
                    gas_label=gas_label_text,
                    density_g_cm3=density_g_cm3,
                    volume_cm3=volume_cm3,
                    E_MeV=fluence_df["E_MeV"].to_numpy(dtype=float),
                    active_flux_or_fluence=cevns_fluence,
                    config=cevns_config,
                    t_grid_points=args.t_grid_points,
                    quantity="count",
                )
                cevns_gas_outdir = model_cevns_outdir / gas_slug
                save_cevns_tables(
                    cevns_gas_outdir,
                    recoil_df=recoil_df,
                    enu_df=enu_df,
                    reco_df=None,
                    skip_save=args.skip_cevns_save,
                )
                write_cevns_spectrum_plots(
                    cevns_gas_outdir,
                    gas_label=gas_label_text,
                    recoil_df=recoil_df,
                    enu_df=enu_df,
                    quantity="count",
                    skip_plots=args.skip_cevns_plots,
                )
                cevns_summary_df["model_name"] = model.name
                cevns_summary_df["fluence_csv"] = str(model.csv_path)
                cevns_summary_df["output_subdir"] = str(Path("cevns") / gas_slug)
                model_cevns_summary_rows.append(cevns_summary_df)
                all_cevns_summary_rows.append(cevns_summary_df)

            total_by_flavor = {
                flavor.key: float(np.sum(spectra[f"count_bin_{OUTPUT_SUFFIX[flavor.key]}"]))
                for flavor in FLAVORS
            }
            total_events = float(np.sum(spectra["count_bin_total"]))
            total_from_T = float(TRAPEZOID(spectra["dN_dTe_total"], spectra["Te_MeV"]))
            total_from_costh = float(TRAPEZOID(spectra["dN_dcosth_total"], spectra["costh"]))
            rel_T = relative_difference(total_from_T, total_events)
            rel_costh = relative_difference(total_from_costh, total_events)

            if rel_T > 1.0e-3 or rel_costh > 1.0e-3:
                print(
                    f"Warning: consistency check for {model.name}/{gas_slug}: "
                    f"dTe rel diff={rel_T:.3g}, dcosth rel diff={rel_costh:.3g}",
                    file=sys.stderr,
                )

            summary_row = {
                "model_name": model.name,
                "fluence_csv": str(model.csv_path),
                "gas_label": gas_label_text,
                "gas_output_dir": str(gas_outdir),
                "gas": gas_name,
                "pressure_atm": pressure_atm,
                "density_g_cm3": density_g_cm3,
                "detector_volume_cm3": volume_cm3,
                "electron_density_cm-3": electron_density,
                "target_electrons": target_electrons,
                "Te_min_MeV": T_low_MeV,
                "Te_max_MeV": T_high_MeV,
                "lower_threshold_source": recoil_window.threshold_source,
                "effective_lower_range_mm": recoil_window.low_range_mm,
                "diffusion_field_kV_cm_used": recoil_window.diffusion_field_v_cm / 1000.0,
                "total_events_nue": total_by_flavor["nue"],
                "total_events_anue": total_by_flavor["anue"],
                "total_events_numu": total_by_flavor["numu"],
                "total_events_anumu": total_by_flavor["anumu"],
                "total_events_nutau": total_by_flavor["nutau"],
                "total_events_anutau": total_by_flavor["anutau"],
                "total_events_all_flavors": total_events,
                "burst_duration_s": float(args.burst_duration_s),
                "average_rate_Hz_all_flavors": total_events / float(args.burst_duration_s),
                "total_events_from_dTe_integral": total_from_T,
                "total_events_from_dcosth_integral": total_from_costh,
                "consistency_rel_diff_dTe": rel_T,
                "consistency_rel_diff_dcosth": rel_costh,
                "molar_mass_g_mol": molar_mass,
                "electrons_per_mixture": electrons_per_mixture,
                "t_grid_points": args.t_grid_points,
                "costh_grid_points": args.costh_grid_points,
            }
            model_summary_rows.append(summary_row)
            all_summary_rows.append(summary_row)

            progress_index += 1
            progress.update(progress_index, f"{model.name}: {gas_label_text}")

        model_summary_df = pd.DataFrame(model_summary_rows)
        model_summary_path = model_outdir / "gas_target_event_summary.csv"
        model_summary_df.to_csv(model_summary_path, index=False)
        write_model_summary_plot(model_outdir, model_summary_df, model.name)
        if cevns_config.enabled:
            write_cevns_summary(
                model_cevns_outdir,
                model_cevns_summary_rows,
                quantity="count",
                skip_save=args.skip_cevns_save,
                skip_plots=args.skip_cevns_plots,
                filename="cevns_event_summary.csv",
            )

    progress.close()

    all_summary_df = pd.DataFrame(all_summary_rows)
    all_summary_path = output_root / "all_models_event_summary.csv"
    all_summary_df.to_csv(all_summary_path, index=False)
    all_cevns_summary_df = write_cevns_summary(
        output_root / "cevns",
        all_cevns_summary_rows,
        quantity="count",
        skip_save=args.skip_cevns_save,
        skip_plots=args.skip_cevns_plots,
        filename="all_models_cevns_event_summary.csv",
    ) if cevns_config.enabled else pd.DataFrame()

    print(f"Processed {len(models)} supernova fluence model(s) from: {fluence_root}")
    print(f"Read {len(gas_df)} gas entries from: {gas_csv}")
    print(
        f"Using detector geometry '{geometry.name}' = "
        f"{geometry.length_m:g} x {geometry.width_m:g} x {geometry.height_m:g} m^3 "
        f"(volume = {volume_cm3 / 1.0e6:g} m^3)"
    )
    print("Primary supernova outputs are expected counts, not rates.")
    print(f"Saved supernova gas-target event outputs in: {output_root.resolve()}")
    print(f"Saved root summary: {all_summary_path.resolve()}")
    if cevns_config.enabled:
        print(
            "CEvNS enabled: "
            f"threshold={cevns_config.nr_threshold_keV:g} keV, "
            f"max={cevns_config.nr_max_keV}, "
            f"form_factor={cevns_config.form_factor}, "
            f"axial_model={cevns_config.axial_model}"
        )
        print(f"CEvNS outputs directory: {(output_root / 'cevns').resolve()}")
        if not all_cevns_summary_df.empty:
            print(
                "CEvNS total supernova counts across all model/gas summary rows: "
                f"{all_cevns_summary_df['total_counts'].sum():.6e}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
