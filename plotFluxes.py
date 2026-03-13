import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


# ============================================================
# User settings
# ============================================================
DATA_DIR = Path("./solar_neutrino_tables")

OUT_COMPONENTS = "solar_neutrino_fluxes_loglog_physical.png"
OUT_TOTAL = "solar_neutrino_total_flux_loglog.png"
OUT_FLAVORS = "solar_neutrino_flavor_fluxes_loglog.png"
OUT_PEE = "solar_neutrino_pee_logx.png"
OUT_CSV = "solar_neutrino_total_flux.csv"

LINE_BIN_WIDTH_MEV = 0.002   # 2 keV
YMIN = 1e-2

EGRID_MIN_MEV = 1e-3         # 1 keV
EGRID_MAX_MEV = 20.0
EGRID_NPTS = 20000


# ============================================================
# Total solar-neutrino source fluxes at Earth (1 AU)
# Units: cm^-2 s^-1
# ============================================================
TOTAL_FLUX = {
    "pp": 5.98e10,
    "pep": 1.44e8,
    "hep": 8.04e3,
    "7Be_total": 4.93e9,
    "8B": 5.46e6,
    "13N": 2.78e8,
    "15O": 2.05e8,
    "17F": 5.29e6,
}

BE7_LINES = [
    (0.3843, 0.1044),   # energy [MeV], branching fraction
    (0.8613, 0.8956),
]

PEP_LINE_ENERGY_MEV = 1.442


# ============================================================
# Minimal MSW-LMA parameters
# ============================================================
# Consistent with PDG 2025 central values for a minimal implementation.
SIN2_THETA12 = 0.307
SIN2_THETA13 = 0.02215
DM21_EV2 = 7.49e-5

# Very simple effective production-region electron densities in mol/cm^3
# Used only for a minimal adiabatic MSW-LMA approximation.
NE_MOL_CM3 = {
    "pp": 58.0,
    "pep": 61.0,
    "7Be": 66.0,
    "8B": 93.0,
    "hep": 93.0,
    "13N": 81.0,
    "15O": 86.0,
    "17F": 86.0,
}


# ============================================================
# Generic helpers
# ============================================================
def is_float(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def compute_bin_edges(x):
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or len(x) < 2:
        raise ValueError("Need at least 2 x points to build bin edges.")

    edges = np.empty(len(x) + 1, dtype=float)
    edges[1:-1] = 0.5 * (x[:-1] + x[1:])
    edges[0] = x[0] - 0.5 * (x[1] - x[0])
    edges[-1] = x[-1] + 0.5 * (x[-1] - x[-2])
    return edges


def sort_unique_xy(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]

    order = np.argsort(x)
    x = x[order]
    y = y[order]

    x_unique, idx = np.unique(x, return_index=True)
    y_unique = y[idx]
    return x_unique, y_unique


def build_global_energy_grid(e_min=EGRID_MIN_MEV, e_max=EGRID_MAX_MEV, n_points=EGRID_NPTS):
    return np.logspace(np.log10(e_min), np.log10(e_max), n_points)


# ============================================================
# Flux normalization / extrapolation
# ============================================================
def normalize_shape_to_physical_flux(E, shape, total_flux):
    """
    Convert a tabulated source shape into physical differential flux:
        dPhi/dE [cm^-2 s^-1 MeV^-1]
    """
    E = np.asarray(E, dtype=float)
    shape = np.asarray(shape, dtype=float)

    valid = np.isfinite(E) & np.isfinite(shape) & (shape >= 0.0)
    E = E[valid]
    shape = shape[valid]

    E, shape = sort_unique_xy(E, shape)

    edges = compute_bin_edges(E)
    dE = np.diff(edges)

    norm = np.sum(shape * dE)
    if norm <= 0:
        raise ValueError("Non-positive normalization encountered.")

    dphi_dE = total_flux * shape / norm
    return E, dphi_dE


def extend_continuum_to_low_energy(E, dphi_dE, E_min_new=1e-3, n_points=400):
    """
    Extend a continuum spectrum down to E_min_new [MeV] with
    a matched E^2 law:
        dPhi/dE ~ A E^2

    Matching is done to the first strictly positive-energy point
    with positive flux, so files containing an explicit (0,0) row
    are handled correctly.
    """
    E = np.asarray(E, dtype=float)
    dphi_dE = np.asarray(dphi_dE, dtype=float)

    E, dphi_dE = sort_unique_xy(E, dphi_dE)

    mask_pos = (E > 0.0) & (dphi_dE > 0.0)
    if not np.any(mask_pos):
        return E, dphi_dE

    first_idx = np.where(mask_pos)[0][0]
    E0 = E[first_idx]
    y0 = dphi_dE[first_idx]

    if E_min_new >= E0:
        return E, dphi_dE

    A = y0 / (E0 ** 2)

    E_low = np.logspace(np.log10(E_min_new), np.log10(E0), n_points, endpoint=False)
    y_low = A * E_low**2

    E_ext = np.concatenate([E_low, E[first_idx:]])
    y_ext = np.concatenate([y_low, dphi_dE[first_idx:]])
    return E_ext, y_ext


def renormalize_differential_flux(E, dphi_dE, total_flux):
    edges = compute_bin_edges(E)
    dE = np.diff(edges)
    integral = np.sum(dphi_dE * dE)

    if integral <= 0:
        raise ValueError("Non-positive integral in renormalize_differential_flux.")

    return dphi_dE * (total_flux / integral)


# ============================================================
# Interpolation / line handling on common grid
# ============================================================
def add_continuum_to_grid(E_src, dphi_dE_src, E_grid):
    return np.interp(E_grid, E_src, dphi_dE_src, left=0.0, right=0.0)


def add_line_to_grid(E_grid, E0, total_flux, bin_width):
    """
    Represent a line as a finite-width top-hat bin, preserving
    the integrated line flux.
    """
    y = np.zeros_like(E_grid)
    half = 0.5 * bin_width
    mask = (E_grid >= E0 - half) & (E_grid < E0 + half)

    if not np.any(mask):
        idx = np.argmin(np.abs(E_grid - E0))
        y[idx] = total_flux / bin_width
        return y

    y[mask] = total_flux / bin_width
    return y


# ============================================================
# Parsers
# ============================================================
def parse_pp_file(path):
    E, Y = [], []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if "P(q)" in s or s.startswith("q"):
                continue
            nums = [float(tok) for tok in s.split() if is_float(tok)]
            for i in range(0, len(nums) - 1, 2):
                E.append(nums[i])
                Y.append(nums[i + 1])
    return np.array(E), np.array(Y)


def parse_two_column_file(path):
    E, Y = [], []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            parts = s.split()
            if len(parts) >= 2 and is_float(parts[0]) and is_float(parts[1]):
                E.append(float(parts[0]))
                Y.append(float(parts[1]))
    return np.array(E), np.array(Y)


def parse_b8_file(path):
    E, Y = [], []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            parts = s.split()
            if len(parts) >= 4 and all(is_float(p) for p in parts[:4]):
                E.append(float(parts[0]))
                Y.append(float(parts[1]))
    return np.array(E), np.array(Y)


def parse_be7_lineshape(path):
    """
    Assumes be7_lineshape.dat contains:
      energy_offset_keV  shape_value

    converted to absolute energy around the 0.8613 MeV line.
    """
    E, Y = [], []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            parts = s.split()
            if len(parts) >= 2 and is_float(parts[0]) and is_float(parts[1]):
                offset_keV = float(parts[0])
                val = float(parts[1])
                E.append(0.8613 + offset_keV / 1000.0)
                Y.append(val)
    return np.array(E), np.array(Y)


def load_spectra(data_dir):
    spectra = {}

    files = {
        "pp": data_dir / "pp.dat",
        "hep": data_dir / "hep.dat",
        "13N": data_dir / "n13.dat",
        "15O": data_dir / "o15.dat",
        "17F": data_dir / "f17.dat",
        "8B": data_dir / "b8.txt",
        "7Be_lineshape": data_dir / "be7_lineshape.dat",
    }

    if files["pp"].exists():
        spectra["pp"] = parse_pp_file(files["pp"])
    if files["hep"].exists():
        spectra["hep"] = parse_two_column_file(files["hep"])
    if files["13N"].exists():
        spectra["13N"] = parse_two_column_file(files["13N"])
    if files["15O"].exists():
        spectra["15O"] = parse_two_column_file(files["15O"])
    if files["17F"].exists():
        spectra["17F"] = parse_two_column_file(files["17F"])
    if files["8B"].exists():
        spectra["8B"] = parse_b8_file(files["8B"])
    if files["7Be_lineshape"].exists():
        spectra["7Be_lineshape"] = parse_be7_lineshape(files["7Be_lineshape"])

    return spectra


# ============================================================
# Minimal adiabatic MSW-LMA
# ============================================================
def pee_msw_lma(E_MeV, source, sin2_theta12=SIN2_THETA12,
                sin2_theta13=SIN2_THETA13, dm21_ev2=DM21_EV2):
    """
    Compute the electron-neutrino survival probability P_ee(E)
    in a minimal adiabatic MSW-LMA approximation.

    Physics model
    -------------
    The function uses the standard approximate 3-flavor reduction

        P_ee^(3nu) ≈ sin^4(theta13) + cos^4(theta13) * P_ee^(2nu),

    where the effective 2-flavor survival probability is evaluated
    in matter at production under the adiabatic approximation:

        P_ee^(2nu) = 1/2 * [1 + cos(2 theta12) * cos(2 theta12^m)].

    The matter mixing angle theta12^m is encoded through

        cos(2 theta12^m) =
            (cos(2 theta12) - beta)
            / sqrt[(cos(2 theta12) - beta)^2 + sin^2(2 theta12)],

    with

        beta = A / Δm21^2

    and matter potential

        A [eV^2] = 1.52e-7 * n_e[mol/cm^3] * E[MeV].

    Inputs
    ------
    E_MeV : float or array-like
        Neutrino energy in MeV.

    source : str
        Source label used to select an effective production-region
        electron density from NE_MOL_CM3. Supported values are:
        "pp", "pep", "7Be", "8B", "hep", "13N", "15O", "17F".

    sin2_theta12 : float
        sin^2(theta12).

    sin2_theta13 : float
        sin^2(theta13).

    dm21_ev2 : float
        Solar mass splitting Δm^2_21 in eV^2.

    Returns
    -------
    pee_3nu : ndarray
        Electron-neutrino survival probability at Earth.

    Important approximation
    -----------------------
    This function does NOT solve full neutrino propagation through the
    Sun. Instead, each source is assigned a single effective electron
    density n_e^eff meant to represent its production region.

    Therefore this is a compact phenomenological MSW-LMA model suitable
    for plotting and detector studies, but not a precision oscillation
    calculation.
    """
    E = np.asarray(E_MeV, dtype=float)

    s13_2 = sin2_theta13
    c13_2 = 1.0 - sin2_theta13
    s13_4 = s13_2**2
    c13_4 = c13_2**2

    sin2_2theta12 = 4.0 * sin2_theta12 * (1.0 - sin2_theta12)
    cos2_theta12 = 1.0 - 2.0 * sin2_theta12

    # A[eV^2] = 1.52e-7 * (n_e in mol/cm^3) * E[MeV]
    ne_eff = NE_MOL_CM3[source] * c13_2
    A = 1.52e-7 * ne_eff * E

    beta = A / dm21_ev2

    cos2_theta12_m = (cos2_theta12 - beta) / np.sqrt(
        (cos2_theta12 - beta) ** 2 + sin2_2theta12
    )

    pee_2nu = 0.5 * (1.0 + cos2_theta12 * cos2_theta12_m)
    pee_3nu = s13_4 + c13_4 * pee_2nu
    return pee_3nu


def split_flux_by_flavor(E_MeV, phi_source, source_name, equal_mu_tau=True, sin2_theta23=0.5):
    pee = pee_msw_lma(E_MeV, source_name)

    phi_e = pee * phi_source
    phi_x = (1.0 - pee) * phi_source

    if equal_mu_tau:
        phi_mu = 0.5 * phi_x
        phi_tau = 0.5 * phi_x
    else:
        phi_tau = sin2_theta23 * phi_x
        phi_mu = (1.0 - sin2_theta23) * phi_x

    return phi_e, phi_mu, phi_tau, pee


# ============================================================
# Saving
# ============================================================
def save_total_flux_file(path, E_grid, phi_total, phi_e, phi_mu, phi_tau):
    header = (
        "E_MeV,total_dphi_dE,nu_e_dphi_dE,nu_mu_dphi_dE,nu_tau_dphi_dE\n"
        "Units: cm^-2 s^-1 MeV^-1\n"
        "Flavor fluxes at Earth after MSW-LMA oscillation."
    )
    data = np.column_stack([E_grid, phi_total, phi_e, phi_mu, phi_tau])
    np.savetxt(path, data, delimiter=",", header=header, comments="")


# ============================================================
# Main
# ============================================================
def main():
    spectra = load_spectra(DATA_DIR)

    print("Loaded sources:", list(spectra.keys()))
    for name, arrs in spectra.items():
        print(name, len(arrs[0]), "points", "Emin=", np.min(arrs[0]), "Emax=", np.max(arrs[0]))

    E_grid = build_global_energy_grid()

    # Source fluxes before oscillation
    total_source_grid = np.zeros_like(E_grid)

    # Flavor fluxes at Earth after oscillation
    total_nue_grid = np.zeros_like(E_grid)
    total_numu_grid = np.zeros_like(E_grid)
    total_nutau_grid = np.zeros_like(E_grid)

    # For plotting components
    component_fluxes = {}
    component_pee = {}

    # --------------------------------------------------------
    # Continuum sources
    # --------------------------------------------------------
    for src in ["pp", "13N", "15O", "17F", "hep", "8B"]:
        if src not in spectra:
            continue

        E_raw, shape_raw = spectra[src]

        E_src, dphi_dE_src = normalize_shape_to_physical_flux(
            E_raw, shape_raw, TOTAL_FLUX[src]
        )
        E_src, dphi_dE_src = extend_continuum_to_low_energy(
            E_src, dphi_dE_src, E_min_new=EGRID_MIN_MEV
        )
        dphi_dE_src = renormalize_differential_flux(
            E_src, dphi_dE_src, TOTAL_FLUX[src]
        )

        phi_e_src, phi_mu_src, phi_tau_src, pee_src = split_flux_by_flavor(
            E_src, dphi_dE_src, src, equal_mu_tau=True
        )

        component_fluxes[src] = (E_src, dphi_dE_src)
        component_pee[src] = (E_src, pee_src)

        total_source_grid += add_continuum_to_grid(E_src, dphi_dE_src, E_grid)
        total_nue_grid += add_continuum_to_grid(E_src, phi_e_src, E_grid)
        total_numu_grid += add_continuum_to_grid(E_src, phi_mu_src, E_grid)
        total_nutau_grid += add_continuum_to_grid(E_src, phi_tau_src, E_grid)

    # --------------------------------------------------------
    # 7Be
    # --------------------------------------------------------
    use_be7_lineshape = "7Be_lineshape" in spectra

    if use_be7_lineshape:
        # thermal lineshape for main 0.8613 MeV line
        E_raw, shape_raw = spectra["7Be_lineshape"]
        E_src, dphi_dE_src = normalize_shape_to_physical_flux(
            E_raw, shape_raw, TOTAL_FLUX["7Be_total"] * 0.8956
        )

        phi_e_src, phi_mu_src, phi_tau_src, pee_src = split_flux_by_flavor(
            E_src, dphi_dE_src, "7Be", equal_mu_tau=True
        )

        component_fluxes["7Be_lineshape"] = (E_src, dphi_dE_src)
        component_pee["7Be"] = (E_src, pee_src)

        total_source_grid += add_continuum_to_grid(E_src, dphi_dE_src, E_grid)
        total_nue_grid += add_continuum_to_grid(E_src, phi_e_src, E_grid)
        total_numu_grid += add_continuum_to_grid(E_src, phi_mu_src, E_grid)
        total_nutau_grid += add_continuum_to_grid(E_src, phi_tau_src, E_grid)

        # low-energy 7Be line
        E0 = 0.3843
        frac = 0.1044
        line_src = add_line_to_grid(E_grid, E0, TOTAL_FLUX["7Be_total"] * frac, LINE_BIN_WIDTH_MEV)
        pee_line = pee_msw_lma(np.array([E0]), "7Be")[0]
        total_source_grid += line_src
        total_nue_grid += pee_line * line_src
        total_numu_grid += 0.5 * (1.0 - pee_line) * line_src
        total_nutau_grid += 0.5 * (1.0 - pee_line) * line_src

    else:
        # both 7Be lines as finite-width bins
        pee_plot_E = []
        pee_plot_y = []
        for E0, frac in BE7_LINES:
            line_src = add_line_to_grid(E_grid, E0, TOTAL_FLUX["7Be_total"] * frac, LINE_BIN_WIDTH_MEV)
            pee_line = pee_msw_lma(np.array([E0]), "7Be")[0]

            total_source_grid += line_src
            total_nue_grid += pee_line * line_src
            total_numu_grid += 0.5 * (1.0 - pee_line) * line_src
            total_nutau_grid += 0.5 * (1.0 - pee_line) * line_src

            pee_plot_E.append(E0)
            pee_plot_y.append(pee_line)

        component_pee["7Be"] = (np.array(pee_plot_E), np.array(pee_plot_y))

    # --------------------------------------------------------
    # pep line
    # --------------------------------------------------------
    pep_src = add_line_to_grid(E_grid, PEP_LINE_ENERGY_MEV, TOTAL_FLUX["pep"], LINE_BIN_WIDTH_MEV)
    pee_pep = pee_msw_lma(np.array([PEP_LINE_ENERGY_MEV]), "pep")[0]

    total_source_grid += pep_src
    total_nue_grid += pee_pep * pep_src
    total_numu_grid += 0.5 * (1.0 - pee_pep) * pep_src
    total_nutau_grid += 0.5 * (1.0 - pee_pep) * pep_src

    component_pee["pep"] = (np.array([PEP_LINE_ENERGY_MEV]), np.array([pee_pep]))

    # --------------------------------------------------------
    # Final checks / save CSV
    # --------------------------------------------------------
    print("source total max =", np.max(total_source_grid))
    print("nu_e   total max =", np.max(total_nue_grid))
    print("nu_mu  total max =", np.max(total_numu_grid))
    print("nu_tau total max =", np.max(total_nutau_grid))

    save_total_flux_file(
        OUT_CSV,
        E_grid,
        total_source_grid,
        total_nue_grid,
        total_numu_grid,
        total_nutau_grid
    )

    # ========================================================
    # Plot 1: source components
    # ========================================================
    fig, ax = plt.subplots(figsize=(12, 8))

    for src in ["pp", "13N", "15O", "17F", "hep", "8B"]:
        if src in component_fluxes:
            E_src, dphi_dE_src = component_fluxes[src]
            positive = dphi_dE_src > 0
            style = {"lw": 2.0}
            if src == "8B":
                style["lw"] = 3.0
                style["linestyle"] = "--"
            ax.plot(E_src[positive], dphi_dE_src[positive], label=src, **style)

    if use_be7_lineshape and "7Be_lineshape" in component_fluxes:
        E_src, dphi_dE_src = component_fluxes["7Be_lineshape"]
        positive = dphi_dE_src > 0
        ax.plot(E_src[positive], dphi_dE_src[positive], lw=2.0, label=r"$^7$Be lineshape (861 keV)")

        be7_low_height = TOTAL_FLUX["7Be_total"] * 0.1044 / LINE_BIN_WIDTH_MEV
        ax.vlines(0.3843, YMIN, be7_low_height, linestyles=":", linewidth=2.2, label=r"$^7$Be (0.3843 MeV)")
    else:
        for E0, frac in BE7_LINES:
            height = TOTAL_FLUX["7Be_total"] * frac / LINE_BIN_WIDTH_MEV
            ax.vlines(E0, YMIN, height, linestyles=":", linewidth=2.2, label=fr"$^7$Be ({E0:.4f} MeV)")

    pep_height = TOTAL_FLUX["pep"] / LINE_BIN_WIDTH_MEV
    ax.vlines(PEP_LINE_ENERGY_MEV, YMIN, pep_height, linestyles="--", linewidth=2.2, label="pep")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(EGRID_MIN_MEV, EGRID_MAX_MEV)
    ax.set_ylim(YMIN, 2e13)
    ax.set_xlabel(r"Neutrino energy $E_\nu$ [MeV]", fontsize=14)
    ax.set_ylabel(r"$d\Phi/dE_\nu$ [cm$^{-2}$ s$^{-1}$ MeV$^{-1}$]", fontsize=14)
    ax.set_title("Solar neutrino source fluxes at Earth (before oscillation)", fontsize=15)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=10, ncol=2)
    fig.tight_layout()
    fig.savefig(OUT_COMPONENTS, dpi=300)
    plt.close(fig)

    # ========================================================
    # Plot 2: total source flux
    # ========================================================
    fig2, ax2 = plt.subplots(figsize=(12, 8))
    positive = total_source_grid > 0
    ax2.plot(E_grid[positive], total_source_grid[positive], lw=2.2, label="Total source flux")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(EGRID_MIN_MEV, EGRID_MAX_MEV)
    ax2.set_ylim(YMIN, 2e13)
    ax2.set_xlabel(r"Neutrino energy $E_\nu$ [MeV]", fontsize=14)
    ax2.set_ylabel(r"$d\Phi/dE_\nu$ [cm$^{-2}$ s$^{-1}$ MeV$^{-1}$]", fontsize=14)
    ax2.set_title("Total solar neutrino source flux at Earth", fontsize=15)
    ax2.grid(True, which="both", alpha=0.25)
    ax2.legend(fontsize=10)
    fig2.tight_layout()
    fig2.savefig(OUT_TOTAL, dpi=300)
    plt.close(fig2)

    # ========================================================
    # Plot 3: flavor-separated fluxes at Earth
    # ========================================================
    fig3, ax3 = plt.subplots(figsize=(12, 8))
    positive = total_source_grid > 0
    ax3.plot(E_grid[positive], total_nue_grid[positive], lw=2.2, label=r"$\nu_e$")
    ax3.plot(E_grid[positive], total_numu_grid[positive], lw=2.2, label=r"$\nu_\mu$")
    ax3.plot(E_grid[positive], total_nutau_grid[positive], lw=2.2, label=r"$\nu_\tau$")
    ax3.plot(E_grid[positive], total_source_grid[positive], lw=2.5, linestyle="--", label="source total")
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.set_xlim(EGRID_MIN_MEV, EGRID_MAX_MEV)
    ax3.set_ylim(YMIN, 2e13)
    ax3.set_xlabel(r"Neutrino energy $E_\nu$ [MeV]", fontsize=14)
    ax3.set_ylabel(r"$d\Phi/dE_\nu$ [cm$^{-2}$ s$^{-1}$ MeV$^{-1}$]", fontsize=14)
    ax3.set_title("Solar neutrino flavor fluxes at Earth (after MSW-LMA)", fontsize=15)
    ax3.grid(True, which="both", alpha=0.25)
    ax3.legend(fontsize=10)
    fig3.tight_layout()
    fig3.savefig(OUT_FLAVORS, dpi=300)
    plt.close(fig3)

    # ========================================================
    # Plot 4: Pee(E)
    # ========================================================
    fig4, ax4 = plt.subplots(figsize=(12, 8))

    pee_order = ["pp", "13N", "15O", "17F", "pep", "7Be", "hep", "8B"]
    for src in pee_order:
        if src not in component_pee:
            continue
        E_src, pee_src = component_pee[src]
        positive = (E_src > 0) & np.isfinite(pee_src)
        style = {"lw": 2.0}
        if src == "8B":
            style["lw"] = 3.0
            style["linestyle"] = "--"
        ax4.plot(E_src[positive], pee_src[positive], label=src, **style)

    ax4.set_xscale("log")
    ax4.set_xlim(EGRID_MIN_MEV, EGRID_MAX_MEV)
    ax4.set_ylim(0.2, 0.6)
    ax4.set_xlabel(r"Neutrino energy $E_\nu$ [MeV]", fontsize=14)
    ax4.set_ylabel(r"$P_{ee}(E)$", fontsize=14)
    ax4.set_title(r"Minimal adiabatic MSW-LMA survival probability $P_{ee}(E)$", fontsize=15)
    ax4.grid(True, which="both", alpha=0.25)
    ax4.legend(fontsize=10, ncol=2)
    fig4.tight_layout()
    fig4.savefig(OUT_PEE, dpi=300)
    plt.close(fig4)

    print(f"Saved {OUT_COMPONENTS}")
    print(f"Saved {OUT_TOTAL}")
    print(f"Saved {OUT_FLAVORS}")
    print(f"Saved {OUT_PEE}")
    print(f"Saved {OUT_CSV}")


if __name__ == "__main__":
    main()