import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Kinematic comparison:
# Borexino-like (non-directional) vs directional detector
# for nu-e elastic scattering
# ============================================================

# ----------------------------
# Physical constant
# ----------------------------
m_e = 511.0  # electron mass in keV

# ----------------------------
# Kinematic relations
# ----------------------------
def neutrino_energy_min(T_e):
    """
    Minimum neutrino energy required to produce an electron
    recoil kinetic energy T_e (nu-e elastic scattering).
    """
    T_e = np.asarray(T_e, dtype=float)
    return 0.5 * (T_e + np.sqrt(T_e**2 + 2.0 * m_e * T_e))


def neutrino_energy_from_Te_theta(T_e, theta_deg):
    """
    Reconstruct neutrino energy from electron kinetic energy T_e
    and recoil angle theta (angle between incoming neutrino and
    outgoing electron).

    cos(theta) = (1 + m_e / E_nu) * sqrt(T_e / (T_e + 2 m_e))

    Solved for E_nu:
        E_nu = m_e / ( cos(theta)/A - 1 )
    where
        A = sqrt(T_e / (T_e + 2 m_e))

    Returns NaN for non-physical combinations.
    """
    T_e = np.asarray(T_e, dtype=float)
    theta_deg = np.asarray(theta_deg, dtype=float)

    theta_rad = np.radians(theta_deg)
    cos_theta = np.cos(theta_rad)
    A = np.sqrt(T_e / (T_e + 2.0 * m_e))
    denom = (cos_theta / A) - 1.0

    return np.where(denom > 0.0, m_e / denom, np.nan)


# ----------------------------
# Detector response models
# ----------------------------
def borexino_energy_resolution_frac(T_e):
    """
    Borexino fractional energy resolution:
        sigma_E / E = 0.05 / sqrt(E / MeV)
    with T_e in keV.
    """
    T_e = np.asarray(T_e, dtype=float)
    return 0.05 / np.sqrt(T_e / 1000.0)


def directional_energy_resolution_frac(T_e):
    """
    Directional detector fractional energy resolution:
    fixed 20% at all energies.
    """
    T_e = np.asarray(T_e, dtype=float)
    return np.full_like(T_e, 0.20, dtype=float)


def directional_angular_resolution_deg(T_e):
    """
    Directional detector RMS angular resolution:
    fixed 20 degrees at all energies.
    """
    T_e = np.asarray(T_e, dtype=float)
    return np.full_like(T_e, 20.0, dtype=float)


# ----------------------------
# Utility: directional band
# ----------------------------
def directional_energy_band(T_meas, theta_meas_deg):
    """
    Compatible neutrino-energy band for a directional detector
    given measured electron kinetic energy T_meas and measured
    recoil angle theta_meas_deg.

    Band is approximated by propagating:
      T_true in [T_meas*(1-res), T_meas*(1+res)]
      theta_true in [theta_meas-sigma_theta, theta_meas+sigma_theta]

    Returns:
      E_low, E_high, valid
    """
    T_meas = np.asarray(T_meas, dtype=float)

    res = directional_energy_resolution_frac(T_meas)
    sigma_theta = directional_angular_resolution_deg(T_meas)

    T_low = T_meas * (1.0 - res)
    T_high = T_meas * (1.0 + res)

    theta_low = np.clip(theta_meas_deg - sigma_theta, 0.0, None)
    theta_high = theta_meas_deg + sigma_theta

    E_low = neutrino_energy_from_Te_theta(T_low, theta_low)
    E_high = neutrino_energy_from_Te_theta(T_high, theta_high)

    valid = np.isfinite(E_low) & np.isfinite(E_high) & (E_low > 0.0) & (E_high > 0.0)
    return E_low, E_high, valid


# ----------------------------
# Helper: get next automatic matplotlib color
# ----------------------------
def get_next_color(ax):
    return next(ax._get_lines.prop_cycler)["color"]


# ============================================================
# Analysis setup
# ============================================================
T_e_min_plot = 10.0       # keV
T_e_max_plot = 10000.0    # keV
n_points = 1000

T_e_meas = np.logspace(np.log10(T_e_min_plot), np.log10(T_e_max_plot), n_points)

threshold_borexino = 80.0
threshold_directional = 10.0

E_nu_max_plot = T_e_max_plot * 10.0
E_nu_min_plot = 30.0  # for log-scale plotting

# Angles for plot 1
theta_values_deg = [0, 5, 10, 20, 30, 45]

# Fixed electron energies for plot 2
fixed_Te_values = [10, 30, 100, 300, 1000, 3000]  # keV

# Angle scan for plot 2
theta_scan = np.linspace(0.0, 89.0, 1000)


# ============================================================
# Borexino-like non-directional band
# ============================================================
mask_bx = T_e_meas >= threshold_borexino
T_bx = T_e_meas[mask_bx]

res_bx = borexino_energy_resolution_frac(T_bx)
T_bx_low = T_bx * (1.0 - res_bx)

# Non-directional: only lower bound is physical, upper edge is plotting ceiling
E_nu_bx_low = neutrino_energy_min(T_bx_low)
E_nu_bx_high = np.full_like(E_nu_bx_low, E_nu_max_plot)

# Useful values for Borexino box sides
E_bx_left = E_nu_bx_low[0]
E_bx_right = E_nu_bx_low[-1]


# ============================================================
# Plot 1: E_nu vs T_e for different measured recoil angles
# ============================================================
fig1, ax1 = plt.subplots(figsize=(11, 8))

# Borexino region: automatic single color for fill + all borders
borex_color = get_next_color(ax1)

ax1.fill_between(
    T_bx,
    E_nu_bx_low,
    E_nu_bx_high,
    color=borex_color,
    alpha=0.12,
    label="Borexino-like (non-directional)"
)

# Curved lower boundary
ax1.plot(T_bx, E_nu_bx_low, color=borex_color, lw=2.0)

# Top boundary
ax1.plot(
    [threshold_borexino, T_e_max_plot],
    [E_nu_max_plot, E_nu_max_plot],
    color=borex_color,
    lw=2.0,
    linestyle=":"
)

# Left vertical side
ax1.plot(
    [threshold_borexino, threshold_borexino],
    [E_bx_left, E_nu_max_plot],
    color=borex_color,
    lw=2.0,
    linestyle=":"
)

# Right vertical side
ax1.plot(
    [T_e_max_plot, T_e_max_plot],
    [E_bx_right, E_nu_max_plot],
    color=borex_color,
    lw=2.0,
    linestyle=":"
)

# Directional bands for different theta_meas, automatic consistent color per set
mask_dir = T_e_meas >= threshold_directional
T_dir = T_e_meas[mask_dir]

for theta_meas_deg in theta_values_deg:
    E_low, E_high, valid = directional_energy_band(T_dir, theta_meas_deg)

    T_plot = T_dir[valid]
    E_low_plot = E_low[valid]
    E_high_plot = E_high[valid]

    if len(T_plot) == 0:
        continue

    this_color = get_next_color(ax1)

    ax1.fill_between(
        T_plot,
        E_low_plot,
        E_high_plot,
        color=this_color,
        alpha=0.18,
        label=fr"Directional, $\theta_{{\rm meas}}={theta_meas_deg}^\circ$"
    )
    ax1.plot(T_plot, E_low_plot, color=this_color, lw=1.6)
    ax1.plot(T_plot, E_high_plot, color=this_color, lw=1.6)

# Add text box with detector performance
textstr = f"Borexino-like detector:\n" \
          f"• Energy threshold: {threshold_borexino} keV\n" \
          f"• Energy resolution: 5%/√(E[MeV])\n\n" \
          f"Directional detector:\n" \
          f"• Energy threshold: {threshold_directional} keV\n" \
          f"• Angular resolution: 20° RMS\n" \
          f"• Energy resolution: 20% RMS"
props = dict(boxstyle="round", facecolor="wheat", alpha=0.2)
plt.text(0.02, 0.98, textstr, transform=plt.gca().transAxes, fontsize=9,
         verticalalignment="top", bbox=props, family="monospace")


ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlim(T_e_min_plot, T_e_max_plot)
ax1.set_ylim(E_nu_min_plot, E_nu_max_plot)

ax1.set_xlabel(r"Measured electron kinetic energy $T_e$ [keV]", fontsize=18)
ax1.set_ylabel(r"Compatible neutrino energy $E_\nu$ [keV]", fontsize=18)
ax1.set_title(r"$\nu_e$-$e$ elastic scattering: $E_\nu$ vs $T_e$ for different recoil angles", fontsize=15)

ax1.grid(True, which="both", alpha=0.25)
ax1.legend(loc="lower right", fontsize=11, frameon=True)
fig1.tight_layout()
fig1.savefig("enu_vs_te_different_theta_fixedres.png", dpi=500)
plt.close(fig1)


# ============================================================
# Plot 2: E_nu vs theta for fixed electron energies
# ============================================================
fig2, ax2 = plt.subplots(figsize=(11, 8))

for T_meas in fixed_Te_values:
    if T_meas < threshold_directional:
        continue

    res = directional_energy_resolution_frac(np.array([T_meas]))[0]
    sigma_theta = directional_angular_resolution_deg(np.array([T_meas]))[0]

    T_low = T_meas * (1.0 - res)
    T_high = T_meas * (1.0 + res)

    theta_low = np.clip(theta_scan - sigma_theta, 0.0, None)
    theta_high = theta_scan + sigma_theta

    E_low = neutrino_energy_from_Te_theta(T_low, theta_low)
    E_high = neutrino_energy_from_Te_theta(T_high, theta_high)

    valid = np.isfinite(E_low) & np.isfinite(E_high) & (E_low > 0.0) & (E_high > 0.0)
    theta_plot = theta_scan[valid]
    E_low_plot = E_low[valid]
    E_high_plot = E_high[valid]

    if len(theta_plot) == 0:
        continue

    this_color = get_next_color(ax2)

    ax2.fill_between(
        theta_plot,
        E_low_plot,
        E_high_plot,
        color=this_color,
        alpha=0.20,
        label=fr"$T_e={T_meas}$ keV"
    )
    ax2.plot(theta_plot, E_low_plot, color=this_color, lw=1.6)
    ax2.plot(theta_plot, E_high_plot, color=this_color, lw=1.6)

# Add text box with detector performance
textstr = f"Borexino-like detector:\n" \
          f"• Energy threshold: {threshold_borexino} keV\n" \
          f"• Energy resolution: 5%/√(E[MeV])\n\n" \
          f"Directional detector:\n" \
          f"• Energy threshold: {threshold_directional} keV\n" \
          f"• Angular resolution: 20° RMS\n" \
          f"• Energy resolution: 20% RMS"
props = dict(boxstyle="round", facecolor="wheat", alpha=0.2)
plt.text(0.73, 0.98, textstr, transform=plt.gca().transAxes, fontsize=9,
         verticalalignment="top", bbox=props, family="monospace")



ax2.set_yscale("log")
ax2.set_xlim(0.0, 69.0)
ax2.set_ylim(E_nu_min_plot, E_nu_max_plot)

ax2.set_xlabel(r"Measured recoil angle $\theta_{\rm meas}$ [deg]", fontsize=18)
ax2.set_ylabel(r"Compatible neutrino energy $E_\nu$ [keV]", fontsize=18)
ax2.set_title(r"Directional detector: $E_\nu$ vs recoil angle for fixed $T_e$", fontsize=15)

ax2.grid(True, which="both", alpha=0.25)
ax2.legend(loc="lower right", fontsize=11, frameon=True)
fig2.tight_layout()
fig2.savefig("enu_vs_theta_fixed_te_fixedres.png", dpi=500)
plt.close(fig2)