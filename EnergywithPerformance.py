import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Callable, Optional


# ============================================================
# Constants
# ============================================================

M_E_KEV = 511.0  # electron mass in keV
EPS = 1e-12


# ============================================================
# Kinematics
# ============================================================

def neutrino_energy_min(T_e_keV, m_e_keV=M_E_KEV):
    """
    Minimum neutrino energy needed to produce an electron recoil
    kinetic energy T_e in nu-e elastic scattering.

    Parameters
    ----------
    T_e_keV : float or array
        Electron kinetic energy in keV.
    m_e_keV : float
        Electron mass in keV.

    Returns
    -------
    float or array
        Minimum neutrino energy in keV.
    """
    T = np.asarray(T_e_keV, dtype=float)
    T = np.clip(T, 0.0, None)
    return 0.5 * (T + np.sqrt(T**2 + 2.0 * m_e_keV * T))


def electron_momentum_from_Te(T_e_keV, m_e_keV=M_E_KEV):
    """
    Electron momentum from kinetic energy.

    p_e = sqrt(T^2 + 2 m_e T)
    """
    T = np.asarray(T_e_keV, dtype=float)
    T = np.clip(T, 0.0, None)
    return np.sqrt(T**2 + 2.0 * m_e_keV * T)


def neutrino_energy_from_Te_theta(T_e_keV, theta_deg, m_e_keV=M_E_KEV):
    """
    Reconstruct neutrino energy from electron kinetic energy T_e
    and recoil angle theta (deg), where theta is the angle between
    incoming neutrino and outgoing electron.

    Uses:
        E_nu = m_e T / (p_e cos(theta) - T)
    with
        p_e = sqrt(T^2 + 2 m_e T)

    Returns NaN for non-physical combinations.
    """
    T = np.asarray(T_e_keV, dtype=float)
    theta = np.asarray(theta_deg, dtype=float)

    # Physical recoil angles are forward; clip to just below 90 deg
    theta = np.clip(theta, 0.0, 89.999999)
    T = np.clip(T, 0.0, None)

    p = electron_momentum_from_Te(T, m_e_keV=m_e_keV)
    cos_theta = np.cos(np.radians(theta))
    denom = p * cos_theta - T

    out = np.full(np.broadcast(T, theta).shape, np.nan, dtype=float)
    valid = (T > 0.0) & (denom > EPS)
    out[valid] = (m_e_keV * T[valid]) / denom[valid]
    return out


# ============================================================
# Detector performance model
# ============================================================

@dataclass
class DetectorPerformance:
    """
    Generic detector performance container.

    The two performance functions must accept T_e in keV and return:
      - fractional energy resolution sigma_E / E
      - angular resolution in degrees
    """
    name: str
    threshold_keV: float
    energy_resolution_frac: Callable[[np.ndarray], np.ndarray]
    angular_resolution_deg: Optional[Callable[[np.ndarray], np.ndarray]] = None

    def sigma_E_over_E(self, T_e_keV):
        T = np.asarray(T_e_keV, dtype=float)
        return np.asarray(self.energy_resolution_frac(T), dtype=float)

    def sigma_theta_deg(self, T_e_keV):
        T = np.asarray(T_e_keV, dtype=float)
        if self.angular_resolution_deg is None:
            return np.full_like(T, np.nan, dtype=float)
        return np.asarray(self.angular_resolution_deg(T), dtype=float)


# ============================================================
# Example performance functions
# ============================================================

def borexino_resolution_frac(T_e_keV):
    """
    Example Borexino-like fractional energy resolution:
        sigma_E / E = 0.05 / sqrt(E[MeV])

    Input T_e in keV.
    """
    T = np.asarray(T_e_keV, dtype=float)
    T = np.clip(T, 1e-9, None)
    return 0.05 / np.sqrt(T / 1000.0)


def constant_resolution_frac(value):
    """
    Returns a function sigma_E/E = constant.
    """
    def f(T_e_keV):
        T = np.asarray(T_e_keV, dtype=float)
        return np.full_like(T, value, dtype=float)
    return f


def constant_angular_resolution_deg(value):
    """
    Returns a function sigma_theta = constant in deg.
    """
    def f(T_e_keV):
        T = np.asarray(T_e_keV, dtype=float)
        return np.full_like(T, value, dtype=float)
    return f


def powerlaw_angular_resolution_deg(a_deg, T0_keV=100.0, alpha=0.5, floor_deg=0.0):
    """
    Example tunable angular model:
        sigma_theta(T) = max(floor_deg, a_deg * (T/T0)^(-alpha))

    Useful for future studies.
    """
    def f(T_e_keV):
        T = np.asarray(T_e_keV, dtype=float)
        T = np.clip(T, 1e-9, None)
        sigma = a_deg * (T / T0_keV) ** (-alpha)
        sigma = np.maximum(sigma, floor_deg)
        return sigma
    return f


# ============================================================
# Directional energy band utility
# ============================================================

def directional_energy_band(T_meas_keV, theta_meas_deg, perf: DetectorPerformance):
    """
    Compute a simple compatible neutrino-energy band for a directional detector.

    Assumption:
      T_true in [T_meas*(1-res), T_meas*(1+res)]
      theta_true in [theta_meas-sigma_theta, theta_meas+sigma_theta]

    Because E_nu(T, theta) is monotonic increasing in both T and theta
    in the physical region, we approximate:
      E_low  = E(T_low,  theta_low)
      E_high = E(T_high, theta_high)

    Returns
    -------
    E_low, E_high, valid : arrays
    """
    T = np.asarray(T_meas_keV, dtype=float)
    theta = np.asarray(theta_meas_deg, dtype=float)

    res = perf.sigma_E_over_E(T)
    sigma_theta = perf.sigma_theta_deg(T)

    T_low = np.clip(T * (1.0 - res), 0.0, None)
    T_high = np.clip(T * (1.0 + res), 0.0, None)

    theta_low = np.clip(theta - sigma_theta, 0.0, 89.999999)
    theta_high = np.clip(theta + sigma_theta, 0.0, 89.999999)

    E1 = neutrino_energy_from_Te_theta(T_low, theta_low)
    E2 = neutrino_energy_from_Te_theta(T_high, theta_high)

    E_low = np.minimum(E1, E2)
    E_high = np.maximum(E1, E2)

    valid = np.isfinite(E_low) & np.isfinite(E_high) & (E_low > 0.0) & (E_high > 0.0)
    return E_low, E_high, valid


# ============================================================
# Plotting helpers
# ============================================================

def performance_text(perf_bx: DetectorPerformance, perf_dir: DetectorPerformance):
    """
    Small text summary for plot boxes.
    """
    txt = []
    txt.append(f"{perf_bx.name}:")
    txt.append(f"  threshold = {perf_bx.threshold_keV:.1f} keV")
    txt.append("  energy resolution =  5% / sqrt(E[MeV])")
    txt.append("")
    txt.append(f"{perf_dir.name}:")
    txt.append(f"  threshold = {perf_dir.threshold_keV:.1f} keV")
    txt.append("  energy resolution = 20%")
    txt.append("  angular resolution = 20 degrees")
    return "\n".join(txt)


# ============================================================
# Main analysis / plots
# ============================================================

def main():
    # ----------------------------
    # User configuration
    # ----------------------------
    T_e_min_plot = 10.0
    T_e_max_plot = 10000.0
    E_nu_min_plot = 30.0
    E_nu_max_plot = T_e_max_plot * 10.0
    n_points = 1000

    theta_values_deg = [0, 5, 10, 20, 30, 45]
    fixed_Te_values = [10, 30, 100, 300, 1000, 3000]
    theta_scan = np.linspace(0.0, 89.0, 1000)

    # Example detector definitions
    borexino_like = DetectorPerformance(
        name="Borexino-like (non-directional)",
        threshold_keV=80.0,
        energy_resolution_frac=borexino_resolution_frac,
        angular_resolution_deg=None,
    )

    directional = DetectorPerformance(
        name="Directional detector",
        threshold_keV=10.0,
        energy_resolution_frac=constant_resolution_frac(0.20),
        angular_resolution_deg=constant_angular_resolution_deg(20.0),
        # Example alternative:
        # angular_resolution_deg=powerlaw_angular_resolution_deg(
        #     a_deg=35.0, T0_keV=50.0, alpha=0.4, floor_deg=8.0
        # ),
    )

    # ----------------------------
    # Common x-grid
    # ----------------------------
    T_e_meas = np.logspace(np.log10(T_e_min_plot), np.log10(T_e_max_plot), n_points)

    # ============================================================
    # Borexino-like band
    # ============================================================
    mask_bx = T_e_meas >= borexino_like.threshold_keV
    T_bx = T_e_meas[mask_bx]

    res_bx = borexino_like.sigma_E_over_E(T_bx)
    T_bx_low = np.clip(T_bx * (1.0 - res_bx), 0.0, None)

    E_nu_bx_low = neutrino_energy_min(T_bx_low)
    E_nu_bx_high = np.full_like(E_nu_bx_low, E_nu_max_plot)

    # ============================================================
    # Plot 1: E_nu vs T_e for different recoil angles
    # ============================================================
    fig1, ax1 = plt.subplots(figsize=(11, 8))

    # pick one color once
    borex_color = f"C0"

    ax1.fill_between(
        T_bx,
        E_nu_bx_low,
        E_nu_bx_high,
        color=borex_color,
        alpha=0.12,
        label=borexino_like.name
    )
    ax1.plot(T_bx, E_nu_bx_low, color=borex_color, lw=2.0)

    mask_dir = T_e_meas >= directional.threshold_keV
    T_dir = T_e_meas[mask_dir]

    for i, theta_meas_deg in enumerate(theta_values_deg, start=1):
        E_low, E_high, valid = directional_energy_band(T_dir, theta_meas_deg, directional)

        if not np.any(valid):
            continue

        band_color = f"C{i}"

        ax1.fill_between(
            T_dir[valid],
            E_low[valid],
            E_high[valid],
            color=band_color,
            alpha=0.18,
            label=fr"{directional.name}, $\theta_{{\rm meas}}={theta_meas_deg}^\circ$"
        )
        ax1.plot(T_dir[valid], E_low[valid], color=band_color, lw=1.6)
        ax1.plot(T_dir[valid], E_high[valid], color=band_color, lw=1.6)

    ax1.text(
        0.02, 0.98,
        performance_text(borexino_like, directional),
        transform=ax1.transAxes,
        fontsize=9,
        va="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.2),
        family="monospace"
    )

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim(T_e_min_plot, T_e_max_plot)
    ax1.set_ylim(E_nu_min_plot, E_nu_max_plot)

    ax1.set_xlabel(r"Measured electron kinetic energy $T_e$ [keV]", fontsize=16)
    ax1.set_ylabel(r"Compatible neutrino energy $E_\nu$ [keV]", fontsize=16)
    ax1.set_title(r"$\nu_e$-$e$ elastic scattering: $E_\nu$ vs $T_e$", fontsize=14)

    ax1.grid(True, which="both", alpha=0.25)
    ax1.legend(loc="lower right", fontsize=10, frameon=True)
    fig1.tight_layout()
    fig1.savefig("enu_vs_te_different_theta.png", dpi=300)
    plt.close(fig1)

    # ============================================================
    # Plot 2: E_nu vs theta for fixed electron energies
    # ============================================================
    fig2, ax2 = plt.subplots(figsize=(11, 8))

    for i, T_meas in enumerate(fixed_Te_values):
        if T_meas < directional.threshold_keV:
            continue

        E_low, E_high, valid = directional_energy_band(
            np.full_like(theta_scan, T_meas, dtype=float),
            theta_scan,
            directional
        )

        if not np.any(valid):
            continue

        band_color = f"C{i}"

        ax2.fill_between(
            theta_scan[valid],
            E_low[valid],
            E_high[valid],
            color=band_color,
            alpha=0.20,
            label=fr"$T_e={T_meas}$ keV"
        )
        ax2.plot(theta_scan[valid], E_low[valid], color=band_color, lw=1.6)
        ax2.plot(theta_scan[valid], E_high[valid], color=band_color, lw=1.6)

    ax2.text(
        0.70, 0.98,
        performance_text(borexino_like, directional),
        transform=ax2.transAxes,
        fontsize=9,
        va="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.2),
        family="monospace"
    )

    ax2.set_yscale("log")
    ax2.set_xlim(0.0, 69.0)
    ax2.set_ylim(E_nu_min_plot, E_nu_max_plot)

    ax2.set_xlabel(r"Measured recoil angle $\theta_{\rm meas}$ [deg]", fontsize=16)
    ax2.set_ylabel(r"Compatible neutrino energy $E_\nu$ [keV]", fontsize=16)
    ax2.set_title(r"Directional detector: $E_\nu$ vs recoil angle for fixed $T_e$", fontsize=14)

    ax2.grid(True, which="both", alpha=0.25)
    ax2.legend(loc="lower right", fontsize=10, frameon=True)
    fig2.tight_layout()
    fig2.savefig("enu_vs_theta_fixed_te.png", dpi=300)
    plt.close(fig2)

    print("Done. Saved:")
    print("  - enu_vs_te_different_theta.png")
    print("  - enu_vs_theta_fixed_te.png")


if __name__ == "__main__":
    main()