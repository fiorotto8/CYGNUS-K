#!/usr/bin/env python3
"""
Generate analytic supernova neutrino fluence spectra at Earth.

Based on the quasi-thermal fluence parameterization used in:
G. Angloher et al., JCAP 03 (2025) 037,
"Neutrino flux sensitivity to the next galactic core-collapse supernova in COSINUS".

The fluence model is:

    dPhi_i/dE = eps_i / (4*pi*d^2) *
                E^alpha * exp(-E/T_i) / (T_i^(alpha+2) * Gamma(alpha+2))

    T_i = <E_i> / (alpha + 1)

where:
    E        = neutrino energy [MeV]
    eps_i    = emitted energy in flavor i [MeV]
    d        = distance to supernova [cm]
    dPhi/dE  = time-integrated fluence at Earth [cm^-2 MeV^-1]

Notes:
    - SN1987A_like uses the parameters explicitly quoted by Angloher/COSINUS.
    - The 27 M_sun model shown in Angloher/COSINUS is simulation-based and is
      intentionally not reproduced here.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -----------------------------
# Constants
# -----------------------------

ERG_TO_MEV = 6.241509074e5       # 1 erg in MeV
KPC_TO_CM = 3.085677581491367e21 # 1 kpc in cm


# -----------------------------
# Quasi-thermal fluence model
# -----------------------------

def quasi_thermal_fluence(
    energy_mev: np.ndarray,
    distance_kpc: float,
    emitted_energy_erg: float,
    average_energy_mev: float,
    alpha: float,
) -> np.ndarray:
    """
    Time-integrated differential neutrino fluence at Earth.

    Parameters
    ----------
    energy_mev:
        Energy grid [MeV].
    distance_kpc:
        Supernova distance [kpc].
    emitted_energy_erg:
        Total emitted energy in this single flavor/species [erg].
    average_energy_mev:
        Average neutrino energy for this flavor/species [MeV].
    alpha:
        Pinching parameter.

    Returns
    -------
    fluence:
        dPhi/dE [cm^-2 MeV^-1].
    """
    E = np.asarray(energy_mev, dtype=float)
    d_cm = distance_kpc * KPC_TO_CM
    eps_mev = emitted_energy_erg * ERG_TO_MEV
    T = average_energy_mev / (alpha + 1.0)

    spectral_shape = (
        E**alpha
        * np.exp(-E / T)
        / (T ** (alpha + 2.0) * math.gamma(alpha + 2.0))
    )

    fluence = eps_mev / (4.0 * math.pi * d_cm**2) * spectral_shape
    return fluence


def build_model_spectra(
    model_name: str,
    model: dict,
    energy_mev: np.ndarray,
    output_root: Path,
) -> pd.DataFrame:
    """
    Build spectra, save CSV and plot.
    """
    distance_kpc = model["distance_kpc"]
    alpha = model["alpha"]

    # Per physical species.
    nue = quasi_thermal_fluence(
        energy_mev,
        distance_kpc,
        model["eps_erg"]["nue"],
        model["avg_energy_mev"]["nue"],
        alpha,
    )

    anue = quasi_thermal_fluence(
        energy_mev,
        distance_kpc,
        model["eps_erg"]["anue"],
        model["avg_energy_mev"]["anue"],
        alpha,
    )

    # nux and anux are per single heavy-flavor species:
    # nu_mu = nu_tau = nux_single
    # anti_nu_mu = anti_nu_tau = anux_single
    nux_single = quasi_thermal_fluence(
        energy_mev,
        distance_kpc,
        model["eps_erg"]["nux_single"],
        model["avg_energy_mev"]["nux_single"],
        alpha,
    )

    anux_single = quasi_thermal_fluence(
        energy_mev,
        distance_kpc,
        model["eps_erg"]["anux_single"],
        model["avg_energy_mev"]["anux_single"],
        alpha,
    )

    # Grouped heavy flavors.
    nux_group = 2.0 * nux_single       # nu_mu + nu_tau
    anux_group = 2.0 * anux_single     # anti_nu_mu + anti_nu_tau
    heavy_total = nux_group + anux_group

    total_all_flavors = nue + anue + heavy_total

    df = pd.DataFrame(
        {
            "E_MeV": energy_mev,
            "fluence_nue_cm-2_MeV-1": nue,
            "fluence_anue_cm-2_MeV-1": anue,
            "fluence_nux_single_cm-2_MeV-1": nux_single,
            "fluence_anux_single_cm-2_MeV-1": anux_single,
            "fluence_nux_group_numu_nutau_cm-2_MeV-1": nux_group,
            "fluence_anux_group_anumu_anutau_cm-2_MeV-1": anux_group,
            "fluence_heavy_total_4species_cm-2_MeV-1": heavy_total,
            "fluence_total_6species_cm-2_MeV-1": total_all_flavors,
        }
    )

    outdir = output_root / model_name
    outdir.mkdir(parents=True, exist_ok=True)

    csv_path = outdir / f"{model_name}_fluence.csv"
    png_path = outdir / f"{model_name}_fluence.png"

    df.to_csv(csv_path, index=False)

    # Plot in the same style/use-case as Angloher/COSINUS Figure 2:
    # fluence vs E at 10 kpc, separate flavor groups and total.
    plt.figure(figsize=(8, 5.5))
    plt.plot(energy_mev, nue, label=r"$\nu_e$")
    plt.plot(energy_mev, anue, label=r"$\bar{\nu}_e$")
    plt.plot(energy_mev, nux_group, label=r"$\nu_x$ group: $\nu_\mu+\nu_\tau$")
    plt.plot(energy_mev, anux_group, label=r"$\bar{\nu}_x$ group: $\bar{\nu}_\mu+\bar{\nu}_\tau$")
    plt.plot(energy_mev, total_all_flavors, linewidth=2.5, label="Total, 6 species")

    plt.xlabel(r"$E_\nu$ [MeV]")
    plt.ylabel(r"$d\Phi/dE$ [cm$^{-2}$ MeV$^{-1}$]")
    plt.title(f"{model_name}: supernova neutrino fluence at {distance_kpc:g} kpc")
    plt.xlim(0, model["plot_energy_max_mev"])
    plt.ylim(bottom=0)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(png_path, dpi=200)
    plt.close()

    print(f"[OK] {model_name}")
    print(f"     CSV : {csv_path}")
    print(f"     Plot: {png_path}")

    return df


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate analytic, time-integrated supernova neutrino fluence spectra."
    )
    parser.add_argument(
        "--output-root",
        default="supernova_fluence",
        help="Output directory for generated model fluence CSVs and plots.",
    )
    parser.add_argument(
        "--energy-max-mev",
        type=float,
        default=80.0,
        help="Maximum neutrino energy in MeV for the generated grid.",
    )
    parser.add_argument(
        "--energy-points",
        type=int,
        default=4000,
        help="Number of points in the generated energy grid.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.energy_max_mev <= 0.01:
        raise ValueError("--energy-max-mev must be greater than 0.01 MeV.")
    if args.energy_points < 100:
        raise ValueError("--energy-points must be at least 100.")

    output_root = Path(args.output_root)

    # Energy grid.
    # Angloher/COSINUS Figure 2 shows 0--50 MeV.
    energy_mev = np.linspace(0.01, args.energy_max_mev, args.energy_points)

    # -----------------------------
    # Model definitions
    # -----------------------------

    # Angloher/COSINUS SN1987A-like analytic benchmark:
    # <E_nue> = 9 MeV
    # <E_anue> = 12 MeV
    # <E_nux> = <E_anux> = 16 MeV
    # alpha = 3
    # total energy = 3e53 erg, equally shared among 6 species.
    eps_per_species_1987a = 3.0e53 / 6.0

    sn1987a_like = {
        "distance_kpc": 10.0,
        "alpha": 3.0,
        "eps_erg": {
            "nue": eps_per_species_1987a,
            "anue": eps_per_species_1987a,
            "nux_single": eps_per_species_1987a,
            "anux_single": eps_per_species_1987a,
        },
        "avg_energy_mev": {
            "nue": 9.0,
            "anue": 12.0,
            "nux_single": 16.0,
            "anux_single": 16.0,
        },
        "plot_energy_max_mev": 50.0,
    }

    # The 27 M_sun model shown in Angloher/COSINUS is simulation-based and is
    # intentionally not reproduced here. This script only provides the
    # SN1987A-like analytic benchmark.
    models = {
        "SN1987A_like": sn1987a_like,
    }

    for model_name, model in models.items():
        build_model_spectra(model_name, model, energy_mev, output_root)


if __name__ == "__main__":
    main()
