#!/usr/bin/env python3
"""
Reconstructed neutrino-energy spectra for supernova fluence event counts.

This mirrors recoNuEnergyComparison.py, but the input is a time-integrated
supernova fluence, so the accepted truth and reconstructed spectra are counts
per burst rather than rates.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import scatteringPlots as sp
from EnergywithPerformance import neutrino_energy_from_Te_theta, neutrino_energy_min
from gasTargetRates import (
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
    sigma_accepted_rescaled,
)
from recoNuEnergyComparison import (
    ReconstructionResponse,
    build_reco_energy_edges,
    compress_true_energy_support,
    costh_from_E_T,
    d_enu_dir_dT_keV,
    d_enu_dir_dtheta_rad_keV,
    d_enu_min_dT_keV,
    dsigma_dT_rescaled_vectorized,
    histogram_to_density,
    read_response_config,
    smear_rates_into_reco_bins,
)
from supernovaGasTargetEvents import (
    DEFAULT_SUPERNOVA_MODEL,
    FLAVORS,
    OUTPUT_SUFFIX,
    FluenceModel,
    choose_models_for_args,
    default_fluence_root,
    discover_fluence_models,
    read_supernova_fluence_csv,
)


REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_RESPONSE_CONFIG = REPO_ROOT / "reco_response_config.json"
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "supernova_nu_reco_energy_comparison"
M_E_KEV = 511.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Estimate reconstructed supernova-neutrino energy spectra for "
            "energy-only and energy+direction measurements."
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
        help="Output directory for reconstructed supernova spectra",
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
        help="Recoil-energy threshold table",
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
        help="Reconstruction response JSON config",
    )
    parser.add_argument(
        "--volume-cm3",
        type=float,
        default=None,
        help="Optional fiducial detector-volume override in cm^3",
    )
    parser.add_argument(
        "--t-grid-points",
        type=int,
        default=220,
        help="Number of accepted recoil-energy grid points per true-energy sample",
    )
    parser.add_argument(
        "--max-true-energy-points",
        type=int,
        default=600,
        help=(
            "Compress each accepted true-energy support before smearing. "
            "Use 0 to disable compression."
        ),
    )
    return parser.parse_args()


def build_true_count_arrays(
    fluence_df: pd.DataFrame,
    target_electrons: float,
    T_low_MeV: float,
    T_high_MeV: float,
) -> dict[str, np.ndarray]:
    E = fluence_df["E_MeV"].to_numpy(dtype=float)
    dE = sp.bin_widths_from_centers(E)

    truth: dict[str, np.ndarray] = {
        "E_MeV": E,
        "dE_MeV": dE,
    }

    for flavor in FLAVORS:
        suffix = OUTPUT_SUFFIX[flavor.key]
        fluence = fluence_df[flavor.key].to_numpy(dtype=float)
        dN_dEnu = target_electrons * fluence * sigma_accepted_rescaled(
            E,
            channel=flavor.channel,
            T_low_MeV=T_low_MeV,
            T_high_MeV=T_high_MeV,
        )
        truth[f"dN_dEnu_{suffix}"] = dN_dEnu
        truth[f"count_bin_{suffix}"] = dN_dEnu * dE

    truth["count_bin_total"] = sum(truth[f"count_bin_{name}"] for name in OUTPUT_SUFFIX)
    truth["dN_dEnu_total"] = sum(truth[f"dN_dEnu_{name}"] for name in OUTPUT_SUFFIX)
    return truth


def build_estimated_spectra_for_gas(
    fluence_df: pd.DataFrame,
    target_electrons: float,
    T_low_keV: float,
    T_high_keV: float,
    response: ReconstructionResponse,
    t_grid_points: int,
    max_true_energy_points: int,
) -> tuple[pd.DataFrame, dict[str, float]]:
    T_low_MeV = T_low_keV / 1.0e3
    T_high_MeV = T_high_keV / 1.0e3
    truth = build_true_count_arrays(
        fluence_df=fluence_df,
        target_electrons=target_electrons,
        T_low_MeV=T_low_MeV,
        T_high_MeV=T_high_MeV,
    )

    E_true = truth["E_MeV"]
    true_counts_total = truth["count_bin_total"]
    valid_truth = E_true[(E_true > 0.0) & (true_counts_total > 0.0)]
    if len(valid_truth) == 0:
        raise ValueError("No accepted supernova events for this gas target.")

    E_min_truth = max(
        float(neutrino_energy_min(np.array([T_low_keV]))[0] / 1.0e3),
        float(valid_truth.min()),
    )
    E_max_truth = float(E_true.max())
    E_min_reco = max(0.25 * E_min_truth, 1.0e-4)
    E_max_reco = max(3.0 * E_max_truth, E_max_truth + 3.0)
    reco_edges = build_reco_energy_edges(E_min_reco, E_max_reco, response.reco_energy_bins)
    reco_centers = np.sqrt(reco_edges[:-1] * reco_edges[1:])
    reco_widths = np.diff(reco_edges)

    true_hist = np.histogram(E_true, bins=reco_edges, weights=true_counts_total)[0]
    reco_hist_energy_only = np.zeros(len(reco_centers), dtype=float)
    reco_hist_directional = np.zeros(len(reco_centers), dtype=float)

    u_grid = np.linspace(0.0, 1.0, t_grid_points)
    du_grid = sp.bin_widths_from_centers(u_grid)

    truth["count_bin_nux"] = truth["count_bin_numu"] + truth["count_bin_nutau"]
    truth["count_bin_anux"] = truth["count_bin_anumu"] + truth["count_bin_anutau"]
    flavor_payloads = [
        ("count_bin_nue", "nu_e"),
        ("count_bin_anue", "anti_nu_e"),
        ("count_bin_nux", "nu_x"),
        ("count_bin_anux", "anti_nu_x"),
    ]

    for count_key, channel in flavor_payloads:
        E_support, count_support = compress_true_energy_support(
            E_true,
            truth[count_key],
            max_points=max_true_energy_points,
        )

        for E_MeV, accepted_count_bin in zip(E_support, count_support):
            if accepted_count_bin <= 0.0:
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
            cell_counts = accepted_count_bin * shape_weights / total_shape

            T_grid_keV = 1.0e3 * T_grid_MeV
            sigma_T_keV = response.sigma_E_over_E(T_grid_keV, threshold_keV=T_low_keV) * T_grid_keV

            mu_energy_only_MeV = neutrino_energy_min(T_grid_keV) / 1.0e3
            sigma_energy_only_MeV = np.abs(d_enu_min_dT_keV(T_grid_keV)) * sigma_T_keV / 1.0e3
            reco_hist_energy_only += smear_rates_into_reco_bins(
                reco_edges,
                cell_counts,
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
                cell_counts,
                mu_directional_MeV,
                sigma_directional_MeV,
            )

    reco_df = pd.DataFrame(
        {
            "reco_E_nu_bin_low_MeV": reco_edges[:-1],
            "reco_E_nu_bin_high_MeV": reco_edges[1:],
            "reco_E_nu_bin_center_MeV": reco_centers,
            "bin_width_MeV": reco_widths,
            "true_accepted_counts_per_bin": true_hist,
            "reco_energy_only_counts_per_bin": reco_hist_energy_only,
            "reco_energy_plus_direction_counts_per_bin": reco_hist_directional,
            "true_accepted_dN_dE_per_MeV": histogram_to_density(true_hist, reco_edges),
            "reco_energy_only_dN_dE_per_MeV": histogram_to_density(reco_hist_energy_only, reco_edges),
            "reco_energy_plus_direction_dN_dE_per_MeV": histogram_to_density(
                reco_hist_directional,
                reco_edges,
            ),
        }
    )

    metrics = {
        "total_true_accepted_events": float(np.sum(true_hist)),
        "total_reco_energy_only_events": float(np.sum(reco_hist_energy_only)),
        "total_reco_energy_plus_direction_events": float(np.sum(reco_hist_directional)),
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
    model_name: str,
    gas_label_text: str,
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
        reco_df["reco_E_nu_bin_center_MeV"],
        reco_df["true_accepted_dN_dE_per_MeV"],
        "true accepted",
        linewidth=2.2,
    )
    plot_positive(
        ax,
        reco_df["reco_E_nu_bin_center_MeV"],
        reco_df["reco_energy_only_dN_dE_per_MeV"],
        "estimated: energy only",
        linewidth=2.0,
    )
    plot_positive(
        ax,
        reco_df["reco_E_nu_bin_center_MeV"],
        reco_df["reco_energy_plus_direction_dN_dE_per_MeV"],
        "estimated: energy + direction",
        linewidth=2.0,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Reconstructed neutrino energy [MeV]")
    ax.set_ylabel(r"$dN/dE_\nu$ [counts MeV$^{-1}$]")
    ax.set_title(f"{model_name}: reconstructed supernova neutrino-energy spectra\n{gas_label_text}\n{response_text}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(gas_outdir / "reco_neutrino_energy_spectra.png", dpi=PLOT_DPI)
    plt.close(fig)


def process_model(
    model: FluenceModel,
    output_root: Path,
    gas_df: pd.DataFrame,
    recoil_window_df: pd.DataFrame,
    diffusion_df: pd.DataFrame,
    volume_cm3: float,
    response: ReconstructionResponse,
    t_grid_points: int,
    max_true_energy_points: int,
    progress: ProgressBar,
    progress_index: int,
) -> tuple[list[dict[str, float | str]], int]:
    fluence_df = read_supernova_fluence_csv(model.csv_path)
    model_outdir = output_root / model.name
    model_outdir.mkdir(parents=True, exist_ok=True)
    model_rows = []

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
        electron_density, molar_mass, electrons_per_mixture = electron_density_cm3(
            gas_name,
            density_g_cm3,
        )
        target_electrons = electron_density * volume_cm3

        reco_df, metrics = build_estimated_spectra_for_gas(
            fluence_df=fluence_df,
            target_electrons=target_electrons,
            T_low_keV=recoil_window.low_keV,
            T_high_keV=recoil_window.high_keV,
            response=response,
            t_grid_points=t_grid_points,
            max_true_energy_points=max_true_energy_points,
        )
        reco_df.to_csv(gas_outdir / "reco_neutrino_energy_spectra.csv", index=False)
        write_reco_plot(
            gas_outdir,
            model_name=model.name,
            gas_label_text=gas_label_text,
            response=response,
            threshold_keV=recoil_window.low_keV,
            reco_df=reco_df,
        )

        row_payload = {
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
            "Te_min_MeV": recoil_window.low_keV / 1.0e3,
            "Te_max_MeV": recoil_window.high_keV / 1.0e3,
            "lower_threshold_source": recoil_window.threshold_source,
            "effective_lower_range_mm": recoil_window.low_range_mm,
            "diffusion_field_kV_cm_used": recoil_window.diffusion_field_v_cm / 1000.0,
            "response_name": response.name,
            "energy_resolution_at_threshold_frac": response.energy_resolution_at_threshold_frac,
            "energy_resolution_constant_frac": response.energy_resolution_constant_frac,
            "angular_resolution_at_threshold_deg": response.angular_resolution_at_threshold_deg,
            "angular_resolution_constant_deg": response.angular_resolution_constant_deg,
            "reco_energy_bins": response.reco_energy_bins,
            "t_grid_points": t_grid_points,
            "max_true_energy_points": max_true_energy_points,
            "molar_mass_g_mol": molar_mass,
            "electrons_per_mixture": electrons_per_mixture,
            "total_true_accepted_events": metrics["total_true_accepted_events"],
            "total_reco_energy_only_events": metrics["total_reco_energy_only_events"],
            "total_reco_energy_plus_direction_events": metrics[
                "total_reco_energy_plus_direction_events"
            ],
            "reco_energy_min_MeV": metrics["reco_energy_min_MeV"],
            "reco_energy_max_MeV": metrics["reco_energy_max_MeV"],
        }
        model_rows.append(row_payload)
        progress_index += 1
        progress.update(progress_index, f"{model.name}: {gas_label_text}")

    model_summary_path = model_outdir / "reco_neutrino_energy_summary.csv"
    pd.DataFrame(model_rows).to_csv(model_summary_path, index=False)
    return model_rows, progress_index


def main() -> int:
    args = parse_args()
    if args.t_grid_points < 20:
        raise ValueError("--t-grid-points must be at least 20.")
    if args.max_true_energy_points < 0:
        raise ValueError("--max-true-energy-points must be non-negative.")

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
        (response_config, "Response config"),
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
    response = read_response_config(response_config)
    volume_cm3 = float(args.volume_cm3) if args.volume_cm3 is not None else geometry.volume_cm3

    all_rows = []
    progress = ProgressBar(len(models) * len(gas_df), prefix="Supernova reco spectra")
    progress.update(0, "starting")
    progress_index = 0
    for model in models:
        model_rows, progress_index = process_model(
            model,
            output_root=output_root,
            gas_df=gas_df,
            recoil_window_df=recoil_window_df,
            diffusion_df=diffusion_df,
            volume_cm3=volume_cm3,
            response=response,
            t_grid_points=args.t_grid_points,
            max_true_energy_points=args.max_true_energy_points,
            progress=progress,
            progress_index=progress_index,
        )
        all_rows.extend(model_rows)
    progress.close()

    summary_path = output_root / "all_models_reco_energy_summary.csv"
    pd.DataFrame(all_rows).to_csv(summary_path, index=False)

    print(f"Processed {len(models)} supernova fluence model(s) from: {fluence_root}")
    print(f"Read {len(gas_df)} gas entries from: {gas_csv}")
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
        f"reco bins={response.reco_energy_bins}"
    )
    print("Primary supernova reconstructed spectra are counts, not rates.")
    print(f"Saved supernova reconstructed-energy outputs in: {output_root.resolve()}")
    print(f"Saved root summary: {summary_path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
