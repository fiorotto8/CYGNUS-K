#!/usr/bin/env python3
"""
Compute electron CSDA-like ranges for the configured detector gases.

The script reads DetectorNumbers/GasDensities.csv and the stopping-power tables
for the supported mixtures. It writes:

  - DetectorNumbers/electron_ranges.png
  - DetectorNumbers/electron_range_energy_table.csv

The range estimate is a simple detector-study input, not a detailed transport
simulation. It integrates 1/(rho * S_col) over electron kinetic energy, where
S_col is the tabulated collisional mass stopping power in MeV cm^2/g.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from detector_model import read_gas_density_table, read_stopping_power_curve  # noqa: E402


def energy_at_range(target_mm: float, range_mm: np.ndarray, energy_mev: np.ndarray) -> float:
    """Interpolate kinetic energy in MeV at a target range in mm."""

    order = np.argsort(range_mm)
    r_sorted = range_mm[order]
    e_sorted = energy_mev[order]
    r_unique, idx = np.unique(r_sorted, return_index=True)
    e_unique = e_sorted[idx]

    if target_mm < r_unique[0] or target_mm > r_unique[-1]:
        return float("nan")
    return float(np.interp(target_mm, r_unique, e_unique))


def main() -> None:
    gas_csv = SCRIPT_DIR / "GasDensities.csv"
    output_png = SCRIPT_DIR / "electron_ranges.png"
    output_csv = SCRIPT_DIR / "electron_range_energy_table.csv"

    gas_data = read_gas_density_table(gas_csv)
    fig, ax = plt.subplots(figsize=(10, 6))
    summary_rows = []

    for _, row in gas_data.iterrows():
        gas_name = row["Gas"]
        density_g_cm3 = float(row["density @ 293 K (g/cm3)"])

        try:
            range_mm, energy_mev = read_stopping_power_curve(gas_name, density_g_cm3)
        except Exception as exc:
            print(f"Warning: skipping {gas_name} at {density_g_cm3:g} g/cm^3: {exc}")
            continue

        e_at_1mm_mev = energy_at_range(1.0, range_mm, energy_mev)
        e_at_1m_mev = energy_at_range(1000.0, range_mm, energy_mev)
        summary_rows.append(
            {
                "Gas": gas_name,
                "density_g_cm3": density_g_cm3,
                "energy_min_at_1mm_keV": e_at_1mm_mev * 1.0e3
                if np.isfinite(e_at_1mm_mev)
                else np.nan,
                "energy_max_at_1m_keV": e_at_1m_mev * 1.0e3
                if np.isfinite(e_at_1m_mev)
                else np.nan,
            }
        )

        label = f"{gas_name} ({density_g_cm3:g} g/cm^3)"
        ax.plot(energy_mev * 1.0e3, range_mm, label=label)

    ax.set_xlabel("Kinetic energy [keV]")
    ax.set_ylabel("Range [mm]")
    ax.set_title("Electron range vs energy")
    ax.set_xlim(left=1.0, right=10000.0)
    ax.set_ylim(bottom=0.01, top=10000.0)
    ax.axhline(y=1.0, color="red", linestyle="--", alpha=0.5, linewidth=1.0)
    ax.text(1.5, 1.1, "1 mm", fontsize=10, color="red")
    ax.axhline(y=3.0, color="green", linestyle="--", alpha=0.5, linewidth=1.0)
    ax.text(1.5, 3.1, "3 mm", fontsize=10, color="green")
    ax.axhline(y=1000.0, color="blue", linestyle="--", alpha=0.5, linewidth=1.0)
    ax.text(1.5, 1000.5, "1 m", fontsize=10, color="blue")
    ax.legend(title="Gas mixture", fontsize=8, title_fontsize=9, loc="lower right")
    ax.grid(True, alpha=0.3)
    ax.set_xscale("log")
    ax.set_yscale("log")

    fig.tight_layout()
    fig.savefig(output_png, dpi=150)
    plt.close(fig)

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(output_csv, index=False)

    print("Energy thresholds by gas:")
    print(summary_df.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print(f"Saved plot: {output_png}")
    print(f"Saved table: {output_csv}")


if __name__ == "__main__":
    main()
