#!/usr/bin/env python3
"""
CEvNS nuclear target helpers for the C/F/H/He gas constituents in this repo.

The default isotope map is intentionally explicit:

    C  -> 12C,  Z=6,  N=6,  A=12, J=0
    F  -> 19F,  Z=9,  N=10, A=19, J=1/2
    H  -> 1H,   Z=1,  N=0,  A=1,  J=1/2
    He -> 4He,  Z=2,  N=2,  A=4,  J=0

Vector CEvNS is implemented for all four isotopes. By default H and He are
treated as pointlike weak targets, because the Helm model is only a placeholder
for such light systems at the very low momentum transfers relevant here. Passing
form_factor="helm" applies the same Helm expression to every isotope for
comparison.

Axial CEvNS is implemented only where the current validation basis is clear:
12C and 4He are zero because J=0, 19F uses the Hoferichter-Menendez-Schwenk
transverse Sigma-prime response, and 1H axial neutral-current elastic scattering
is set to zero as an explicit validation gap. Hydrogen therefore currently
contributes only through the suppressed vector term.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from detector_model import (
    DEFAULT_GAS_CSV,
    DEFAULT_GEOMETRY_CONFIG,
    isotope_target_counts,
    read_detector_geometry_config,
    read_gas_density_table,
)


GF_GEV = 1.1663787e-5
GEV2_TO_CM2 = 0.389379338e-27
SIN2_THETA_W = 0.23126
AMU_MEV = 931.49410242
HBARC_MEV_FM = 197.3269804
HELM_C_COEFF_FM = 1.23
HELM_C_OFFSET_FM = 0.60
HELM_A_FM = 0.52
HELM_SKIN_FM = 0.9
G_A = 1.27641


@dataclass(frozen=True)
class NuclearTarget:
    """CEvNS target isotope metadata."""

    key: str
    label: str
    element: str
    Z: int
    N: int
    A: int
    J: float
    default_form_factor: str
    axial_model: str

    @property
    def mass_MeV(self) -> float:
        # Atomic binding/electron-mass corrections are far below the precision
        # targeted by the current detector-design calculations.
        return self.A * AMU_MEV

    @property
    def mass_GeV(self) -> float:
        return self.mass_MeV * 1.0e-3

    @property
    def weak_charge(self) -> float:
        return weak_charge(self)


NUCLEAR_TARGETS = {
    "C12": NuclearTarget(
        key="C12",
        label="12C",
        element="C",
        Z=6,
        N=6,
        A=12,
        J=0.0,
        default_form_factor="helm",
        axial_model="spin_zero",
    ),
    "F19": NuclearTarget(
        key="F19",
        label="19F",
        element="F",
        Z=9,
        N=10,
        A=19,
        J=0.5,
        default_form_factor="helm",
        axial_model="hms_sigma_prime",
    ),
    "H1": NuclearTarget(
        key="H1",
        label="1H",
        element="H",
        Z=1,
        N=0,
        A=1,
        J=0.5,
        default_form_factor="pointlike",
        axial_model="not_implemented",
    ),
    "He4": NuclearTarget(
        key="He4",
        label="4He",
        element="He",
        Z=2,
        N=2,
        A=4,
        J=0.0,
        default_form_factor="pointlike",
        axial_model="spin_zero",
    ),
}

TARGET_ALIASES = {
    "C": "C12",
    "12C": "C12",
    "C12": "C12",
    "C-12": "C12",
    "F": "F19",
    "19F": "F19",
    "F19": "F19",
    "F-19": "F19",
    "H": "H1",
    "1H": "H1",
    "H1": "H1",
    "H-1": "H1",
    "p": "H1",
    "proton": "H1",
    "He": "He4",
    "4He": "He4",
    "He4": "He4",
    "He-4": "He4",
}
TARGET_ALIASES.update(
    {alias.lower(): key for alias, key in list(TARGET_ALIASES.items())}
)

# Table VIII of Papers/CEvENS/HMS_axialFitparams.pdf. For 19F only the L=1
# transverse Sigma-prime amplitudes are nonzero in the fit table.
F19_HMS_B_FM = 1.7623
F19_SIGMA_PRIME_COEFFS = {
    "p": np.array([0.269513, -0.18098, 0.0296873], dtype=float),
    "n": np.array([-0.00113172, 0.00038188, 0.000744991], dtype=float),
}


def resolve_target(target: str | NuclearTarget) -> NuclearTarget:
    if isinstance(target, NuclearTarget):
        return target

    key = TARGET_ALIASES.get(str(target).strip())
    if key is None:
        raise KeyError(
            f"Unknown CEvNS isotope '{target}'. Valid entries are: "
            + ", ".join(sorted(NUCLEAR_TARGETS))
        )
    return NUCLEAR_TARGETS[key]


def weak_charge(target: str | NuclearTarget, sin2_theta_w: float = SIN2_THETA_W) -> float:
    """Return Q_W = Z(1 - 4 sin^2(theta_W)) - N."""

    isotope = resolve_target(target)
    return isotope.Z * (1.0 - 4.0 * sin2_theta_w) - isotope.N


def tmax_nuclear_keV(E_MeV: np.ndarray | float, target: str | NuclearTarget) -> np.ndarray:
    """Maximum nuclear recoil energy in keV for a neutrino energy in MeV."""

    isotope = resolve_target(target)
    E = np.asarray(E_MeV, dtype=float)
    tmax_MeV = 2.0 * E**2 / (isotope.mass_MeV + 2.0 * E)
    return tmax_MeV * 1.0e3


def emin_from_recoil_keV(T_keV: np.ndarray | float, target: str | NuclearTarget) -> np.ndarray:
    """Minimum neutrino energy in MeV needed to produce recoil T_keV."""

    isotope = resolve_target(target)
    T_MeV = np.asarray(T_keV, dtype=float) * 1.0e-3
    return 0.5 * (T_MeV + np.sqrt(T_MeV**2 + 2.0 * isotope.mass_MeV * T_MeV))


def momentum_transfer_MeV(T_keV: np.ndarray | float, target: str | NuclearTarget) -> np.ndarray:
    """Nonrelativistic momentum transfer q = sqrt(2 m_A T)."""

    isotope = resolve_target(target)
    T_MeV = np.asarray(T_keV, dtype=float) * 1.0e-3
    return np.sqrt(np.maximum(2.0 * isotope.mass_MeV * T_MeV, 0.0))


def helm_form_factor(q_MeV: np.ndarray | float, target: str | NuclearTarget) -> np.ndarray:
    """
    Helm weak form factor.

    Uses F(q)=3 j1(qR_n)/(qR_n) exp[-(qs)^2/2] with
    R_n^2 = c^2 + (7/3) pi^2 a^2 - 5s^2,
    c = 1.23 A^(1/3) - 0.60 fm, a = 0.52 fm, and s = 0.90 fm.
    For H and He this is retained only as a comparison placeholder.
    """

    isotope = resolve_target(target)
    q_fm_inv = np.asarray(q_MeV, dtype=float) / HBARC_MEV_FM
    c_fm = HELM_C_COEFF_FM * isotope.A ** (1.0 / 3.0) - HELM_C_OFFSET_FM
    rn2_fm2 = c_fm**2 + (7.0 / 3.0) * np.pi**2 * HELM_A_FM**2 - 5.0 * HELM_SKIN_FM**2
    if rn2_fm2 <= 0.0:
        raise ValueError(f"Helm R_n^2 is non-positive for {isotope.label}: {rn2_fm2:g} fm^2")
    x = q_fm_inv * np.sqrt(rn2_fm2)
    qs = q_fm_inv * HELM_SKIN_FM

    sphere = np.ones_like(x, dtype=float)
    mask = np.abs(x) > 1.0e-6
    xm = x[mask]
    sphere[mask] = 3.0 * (np.sin(xm) - xm * np.cos(xm)) / xm**3
    if np.any(~mask):
        xs = x[~mask]
        sphere[~mask] = 1.0 - xs**2 / 10.0 + xs**4 / 280.0

    return sphere * np.exp(-0.5 * qs**2)


def weak_form_factor(
    q_MeV: np.ndarray | float,
    target: str | NuclearTarget,
    form_factor: str = "default",
) -> np.ndarray:
    """
    Return the vector weak form factor.

    form_factor="default" uses Helm for 12C and 19F, and pointlike F_W=1 for
    1H and 4He. form_factor="helm" applies Helm to every isotope.
    """

    isotope = resolve_target(target)
    mode = form_factor.lower().strip()
    if mode == "default":
        mode = isotope.default_form_factor

    if mode in {"point", "pointlike", "none", "unity"}:
        return np.ones_like(np.asarray(q_MeV, dtype=float), dtype=float)
    if mode == "helm":
        return helm_form_factor(q_MeV, isotope)
    raise ValueError("form_factor must be 'default', 'helm', or 'pointlike'.")


def hms_sigma_prime_amplitudes_19f(T_keV: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
    """
    Return HMS transverse Sigma-prime proton/neutron amplitudes for 19F.

    The fit variable is u=(qb)^2/2 with q in fm^-1 and b=1.7623 fm:
    F(u)=exp(-u/2) sum_i c_i u^i.
    """

    q_fm_inv = momentum_transfer_MeV(T_keV, "F19") / HBARC_MEV_FM
    u = 0.5 * (q_fm_inv * F19_HMS_B_FM) ** 2

    amplitudes = []
    for coeffs in (F19_SIGMA_PRIME_COEFFS["p"], F19_SIGMA_PRIME_COEFFS["n"]):
        poly = np.zeros_like(u, dtype=float)
        for power, coeff in enumerate(coeffs):
            poly += coeff * u**power
        amplitudes.append(np.exp(-0.5 * u) * poly)
    return amplitudes[0], amplitudes[1]


def hms_transverse_sigma_prime_structure_19f(
    T_keV: np.ndarray | float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return transverse S00, S01, S11 structure factors for 19F."""

    fp, fn = hms_sigma_prime_amplitudes_19f(T_keV)
    f_plus = fp + fn
    f_minus = fp - fn
    s00 = f_plus**2
    s01 = 2.0 * f_plus * f_minus
    s11 = f_minus**2
    return s00, s01, s11


def axial_form_factor(
    T_keV: np.ndarray | float,
    target: str | NuclearTarget,
    g_a: float = G_A,
    strange_axial: float = 0.0,
    axial_model: str = "hoferichter_19f_fast",
) -> np.ndarray:
    """
    Pure axial CEvNS form factor F_A(q^2) in the Eq. (66)/(67) convention.

    Hydrogen axial neutral-current elastic scattering is intentionally not
    implemented yet and returns zero.
    """

    isotope = resolve_target(target)
    T = np.asarray(T_keV, dtype=float)

    model = axial_model.lower().strip()
    if model == "none":
        return np.zeros_like(T, dtype=float)

    if isotope.axial_model in {"spin_zero", "not_implemented"}:
        return np.zeros_like(T, dtype=float)
    if isotope.axial_model != "hms_sigma_prime":
        raise ValueError(f"Unsupported axial model for {isotope.label}: {isotope.axial_model}")
    if model not in {"hoferichter_19f_fast", "hoferichter_19f_central", "toy"}:
        raise ValueError(f"Unsupported CEvNS axial model: {axial_model}")
    if model == "hoferichter_19f_central":
        raise NotImplementedError(
            "hoferichter_19f_central needs the HMS two-body delta_0(q^2) "
            "constants/convention pinned down from the local reference before use."
        )

    if model == "toy":
        T_for_response = np.zeros_like(T, dtype=float)
    else:
        T_for_response = T

    s00, s01, s11 = hms_transverse_sigma_prime_structure_19f(T_for_response)
    return (
        8.0
        * np.pi
        / (2.0 * isotope.J + 1.0)
        * (strange_axial**2 * s00 - g_a * strange_axial * s01 + g_a**2 * s11)
    )


def dsigma_dT_vector(
    E_MeV: np.ndarray | float,
    T_keV: np.ndarray | float,
    target: str | NuclearTarget,
    form_factor: str = "default",
    sin2_theta_w: float = SIN2_THETA_W,
) -> np.ndarray:
    """Vector CEvNS d sigma / dT in cm^2 / keV."""

    isotope = resolve_target(target)
    E = np.asarray(E_MeV, dtype=float)
    T = np.asarray(T_keV, dtype=float)
    T_MeV = T * 1.0e-3
    q = momentum_transfer_MeV(T, isotope)
    fw = weak_form_factor(q, isotope, form_factor=form_factor)
    qw = weak_charge(isotope, sin2_theta_w=sin2_theta_w)

    bracket = 1.0 - isotope.mass_MeV * T_MeV / (2.0 * E**2) - T_MeV / E
    prefactor = GF_GEV**2 * isotope.mass_GeV / (4.0 * np.pi)
    out = prefactor * bracket * qw**2 * fw**2 * GEV2_TO_CM2 * 1.0e-6
    valid = (E > 0.0) & (T >= 0.0) & (T <= tmax_nuclear_keV(E, isotope)) & (bracket > 0.0)
    return np.where(valid, out, 0.0)


def dsigma_dT_axial(
    E_MeV: np.ndarray | float,
    T_keV: np.ndarray | float,
    target: str | NuclearTarget,
    g_a: float = G_A,
    strange_axial: float = 0.0,
    axial_model: str = "hoferichter_19f_fast",
) -> np.ndarray:
    """Pure axial CEvNS d sigma / dT in cm^2 / keV."""

    isotope = resolve_target(target)
    E = np.asarray(E_MeV, dtype=float)
    T = np.asarray(T_keV, dtype=float)
    T_MeV = T * 1.0e-3
    fa = axial_form_factor(
        T,
        isotope,
        g_a=g_a,
        strange_axial=strange_axial,
        axial_model=axial_model,
    )

    bracket = 1.0 + isotope.mass_MeV * T_MeV / (2.0 * E**2) - T_MeV / E
    prefactor = GF_GEV**2 * isotope.mass_GeV / (4.0 * np.pi)
    out = prefactor * bracket * fa * GEV2_TO_CM2 * 1.0e-6
    valid = (E > 0.0) & (T >= 0.0) & (T <= tmax_nuclear_keV(E, isotope)) & (bracket > 0.0)
    return np.where(valid, out, 0.0)


def dsigma_dT_total(
    E_MeV: np.ndarray | float,
    T_keV: np.ndarray | float,
    target: str | NuclearTarget,
    form_factor: str = "default",
    include_axial: bool = True,
    axial_model: str = "hoferichter_19f_fast",
) -> np.ndarray:
    """Vector plus pure-axial CEvNS d sigma / dT in cm^2 / keV."""

    vector = dsigma_dT_vector(E_MeV, T_keV, target, form_factor=form_factor)
    if not include_axial:
        return vector
    return vector + dsigma_dT_axial(E_MeV, T_keV, target, axial_model=axial_model)


def gas_cevns_target_counts(
    gas_name: str,
    density_g_cm3: float,
    volume_cm3: float,
) -> tuple[dict[str, float], float, dict[str, float]]:
    """Return CEvNS isotope target counts for a gas entry."""

    return isotope_target_counts(gas_name, density_g_cm3, volume_cm3)


def build_gas_summary(
    gas_df: pd.DataFrame,
    volume_cm3: float,
    E_MeV: float,
    T_keV: float,
    form_factor: str,
    include_axial: bool,
) -> pd.DataFrame:
    rows: list[dict[str, float | str]] = []
    for _, row in gas_df.iterrows():
        gas_name = str(row["Gas"])
        density = float(row["density @ 293 K (g/cm3)"])
        pressure_atm = float(row["Pressure (atm)"])
        counts, molar_mass, isotopes_per_mixture = gas_cevns_target_counts(
            gas_name,
            density,
            volume_cm3,
        )

        for isotope_key in sorted(counts, key=lambda key: NUCLEAR_TARGETS[key].A):
            target = NUCLEAR_TARGETS[isotope_key]
            vector = float(dsigma_dT_vector(E_MeV, T_keV, target, form_factor=form_factor))
            axial = float(dsigma_dT_axial(E_MeV, T_keV, target)) if include_axial else 0.0
            rows.append(
                {
                    "gas": gas_name,
                    "pressure_atm": pressure_atm,
                    "density_g_cm3": density,
                    "molar_mass_g_mol": molar_mass,
                    "isotope": target.label,
                    "isotope_key": isotope_key,
                    "isotopes_per_mixture": isotopes_per_mixture[isotope_key],
                    "target_count": counts[isotope_key],
                    "Q_W": target.weak_charge,
                    "J": target.J,
                    "axial_model": target.axial_model,
                    "Tmax_keV_at_E": float(tmax_nuclear_keV(E_MeV, target)),
                    "dsigma_vector_cm2_per_keV": vector,
                    "dsigma_axial_cm2_per_keV": axial,
                    "dsigma_total_cm2_per_keV": vector + axial,
                }
            )
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize CEvNS isotope targets and sample differential cross sections."
    )
    parser.add_argument("--gas-csv", default=str(DEFAULT_GAS_CSV), help="Gas density table")
    parser.add_argument(
        "--geometry-config",
        default=str(DEFAULT_GEOMETRY_CONFIG),
        help="Detector geometry JSON",
    )
    parser.add_argument(
        "--energy-mev",
        type=float,
        default=10.0,
        help="Sample neutrino energy for d sigma / dT [MeV]",
    )
    parser.add_argument(
        "--recoil-kev",
        type=float,
        default=0.1,
        help="Sample nuclear recoil energy for d sigma / dT [keV]",
    )
    parser.add_argument(
        "--form-factor",
        choices=["default", "helm", "pointlike"],
        default="default",
        help="Vector form-factor mode",
    )
    parser.add_argument(
        "--no-axial",
        action="store_true",
        help="Report vector-only sample cross sections",
    )
    parser.add_argument(
        "--out",
        default="",
        help="Optional CSV output path; stdout is used when omitted",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    gas_csv = Path(args.gas_csv)
    geometry_config = Path(args.geometry_config)
    geometry = read_detector_geometry_config(geometry_config)
    gas_df = read_gas_density_table(gas_csv)

    summary = build_gas_summary(
        gas_df=gas_df,
        volume_cm3=geometry.volume_cm3,
        E_MeV=args.energy_mev,
        T_keV=args.recoil_kev,
        form_factor=args.form_factor,
        include_axial=not args.no_axial,
    )

    if args.out:
        out_path = Path(args.out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(out_path, index=False)
        print(f"Saved CEvNS target summary to: {out_path.resolve()}")
    else:
        print(summary.to_csv(index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
