import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# User settings
# ============================================================
DATA_DIR = Path("./solar_neutrino_tables")   # folder containing the downloaded files

OUTFILE = "solar_neutrino_fluxes_loglog_physical.png"

# Finite bin width used to display monochromatic lines
# (pep and 7Be) as cm^-2 s^-1 MeV^-1 on the plot
LINE_BIN_WIDTH_MEV = 0.002  # 2 keV

# Minimum y shown on log axis
YMIN = 1e-2

# ============================================================
# Total solar neutrino fluxes at Earth (B16-GS98 benchmark)
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

# 7Be branching fractions
BE7_LINES = [
    (0.3843, 0.1044),   # energy in MeV, branching fraction
    (0.8613, 0.8956),
]

PEP_LINE_ENERGY_MEV = 1.442


# ============================================================
# Helpers
# ============================================================
def is_float(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def compute_bin_edges(x):
    """
    Build bin edges from tabulated bin centers x.
    Works for nonuniform spacing.
    """
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or len(x) < 2:
        raise ValueError("Need at least 2 x points to build bin edges.")

    edges = np.empty(len(x) + 1, dtype=float)
    edges[1:-1] = 0.5 * (x[:-1] + x[1:])
    edges[0] = x[0] - 0.5 * (x[1] - x[0])
    edges[-1] = x[-1] + 0.5 * (x[-1] - x[-2])
    return edges


def normalize_shape_to_physical_flux(E, shape, total_flux):
    """
    Convert a tabulated shape into physical differential flux:
        dPhi/dE [cm^-2 s^-1 MeV^-1]

    using bin-aware normalization:
        weight_i ~ shape_i * dE_i
        sum(weight_i) = 1
        dPhi/dE_i = total_flux * shape_i / sum(shape_j dE_j)
    """
    E = np.asarray(E, dtype=float)
    shape = np.asarray(shape, dtype=float)

    valid = np.isfinite(E) & np.isfinite(shape) & (shape >= 0.0)
    E = E[valid]
    shape = shape[valid]

    order = np.argsort(E)
    E = E[order]
    shape = shape[order]

    # Remove duplicate energies if present
    E_unique, idx = np.unique(E, return_index=True)
    shape = shape[idx]
    E = E_unique

    edges = compute_bin_edges(E)
    dE = np.diff(edges)

    norm = np.sum(shape * dE)
    if norm <= 0:
        raise ValueError("Non-positive normalization encountered.")

    dphi_dE = total_flux * shape / norm
    return E, dphi_dE


# ============================================================
# Parsers for the downloaded files
# ============================================================
def parse_pp_file(path):
    """
    pp.dat format:
    repeated (q, P(q)) pairs on each line
    """
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
    """
    hep.dat, n13.dat, o15.dat, f17.dat
    """
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
    """
    b8.txt format:
    E, best, +3sigma, -3sigma
    We use the best column.
    """
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
    be7_lineshape.dat:
    assumed to contain energy offset from 0.8613 MeV line in keV
    and line-shape values.
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


# ============================================================
# Load all available spectra
# ============================================================
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


def build_global_energy_grid(spectra, line_bin_width_mev, e_min=1e-3, e_max=20.0, n_points=20000):
    """
    Build a common energy grid for the total flux file.
    A dense log grid is usually convenient for plotting and saving.
    """
    return np.logspace(np.log10(e_min), np.log10(e_max), n_points)


def add_continuum_to_grid(E_src, dphi_dE_src, E_grid):
    """
    Interpolate a continuum spectrum onto the common grid.
    Outside the tabulated range, set flux to zero.
    """
    return np.interp(E_grid, E_src, dphi_dE_src, left=0.0, right=0.0)


def add_line_to_grid(E_grid, E0, total_flux, bin_width):
    """
    Represent a line source as a finite-width top-hat bin:
        height = total_flux / bin_width
    so that the integrated flux is preserved.

    The bin is centered at E0.
    """
    y = np.zeros_like(E_grid)
    half = 0.5 * bin_width
    mask = (E_grid >= E0 - half) & (E_grid < E0 + half)
    y[mask] = total_flux / bin_width
    return y


def save_total_flux_file(path, E_grid, total_dphi_dE):
    """
    Save the total differential flux to a text file.
    """
    header = (
        "E_MeV,total_dphi_dE_cm^-2_s^-1_MeV^-1\n"
        "Global solar neutrino flux on a common grid.\n"
        "Line sources (pep, 7Be) are represented with finite bin width."
    )
    data = np.column_stack([E_grid, total_dphi_dE])
    np.savetxt(path, data, delimiter=",", header=header, comments="")

def extend_to_low_energy(E, dphi_dE, E_min_new=1e-3, n_points=200):
    """
    Extend a continuum spectrum down to E_min_new [MeV]
    using a low-energy E^2 extrapolation matched to the
    first tabulated point.

    Parameters
    ----------
    E : array
        Energy grid in MeV, ascending.
    dphi_dE : array
        Differential flux in cm^-2 s^-1 MeV^-1.
    E_min_new : float
        New minimum energy in MeV.
    n_points : int
        Number of added points below the original minimum.

    Returns
    -------
    E_ext, dphi_ext
    """
    E = np.asarray(E, dtype=float)
    dphi_dE = np.asarray(dphi_dE, dtype=float)

    order = np.argsort(E)
    E = E[order]
    dphi_dE = dphi_dE[order]

    E0 = E[0]
    y0 = dphi_dE[0]

    if E_min_new >= E0:
        return E, dphi_dE

    # match y = A E^2 at the first tabulated point
    A = y0 / (E0**2)

    E_low = np.logspace(np.log10(E_min_new), np.log10(E0), n_points, endpoint=False)
    y_low = A * E_low**2

    E_ext = np.concatenate([E_low, E])
    dphi_ext = np.concatenate([y_low, dphi_dE])

    return E_ext, dphi_ext

def extend_continuum_to_low_energy(E, dphi_dE, E_min_new=1e-3, n_points=400):
    """
    Extend a continuum spectrum down to E_min_new [MeV] with a matched E^2 law:
        dPhi/dE ~ A * E^2
    matched to the first strictly positive tabulated point.

    This avoids problems when the file contains an explicit (0,0) row.
    """
    E = np.asarray(E, dtype=float)
    dphi_dE = np.asarray(dphi_dE, dtype=float)

    order = np.argsort(E)
    E = E[order]
    dphi_dE = dphi_dE[order]

    # Find first strictly positive-energy point with positive flux
    mask_pos = (E > 0.0) & (dphi_dE > 0.0)
    if not np.any(mask_pos):
        return E, dphi_dE

    first_idx = np.where(mask_pos)[0][0]
    E0 = E[first_idx]
    y0 = dphi_dE[first_idx]

    # If already below requested threshold, do nothing
    if E_min_new >= E0:
        return E, dphi_dE

    # Match A * E^2 to first positive tabulated point
    A = y0 / (E0 ** 2)

    E_low = np.logspace(np.log10(E_min_new), np.log10(E0), n_points, endpoint=False)
    y_low = A * E_low**2

    # Keep only original points at/above E0 to avoid overlap with the added low-energy part
    E_ext = np.concatenate([E_low, E[first_idx:]])
    y_ext = np.concatenate([y_low, dphi_dE[first_idx:]])

    return E_ext, y_ext


def renormalize_differential_flux(E, dphi_dE, total_flux):
    """
    Renormalize tabulated dPhi/dE so that its integral equals total_flux.
    """
    edges = compute_bin_edges(E)
    dE = np.diff(edges)
    integral = np.sum(dphi_dE * dE)

    if integral <= 0:
        raise ValueError("Non-positive integral in renormalize_differential_flux.")

    return dphi_dE * (total_flux / integral)

# ============================================================
# Plot
# ============================================================
def main():
    spectra = load_spectra(DATA_DIR)

    print("Loaded sources:", list(spectra.keys()))
    for name, arrs in spectra.items():
        print(name, len(arrs[0]), "points", "Emin=", np.min(arrs[0]), "Emax=", np.max(arrs[0]))

    # ============================================================
    # Build common grid from 1 keV to 20 MeV
    # ============================================================
    E_grid = build_global_energy_grid(
        spectra,
        line_bin_width_mev=LINE_BIN_WIDTH_MEV,
        e_min=1e-3,   # 1 keV
        e_max=20.0,
        n_points=20000
    )

    total_flux_grid = np.zeros_like(E_grid)

    # Store component spectra after extension, for plotting
    component_fluxes = {}

    # ============================================================
    # Add continuous sources
    # ============================================================
    for src in ["pp", "13N", "15O", "17F", "hep", "8B"]:
        if src not in spectra:
            continue

        E_raw, shape_raw = spectra[src]

        # Convert tabulated shape to physical dPhi/dE
        E_src, dphi_dE_src = normalize_shape_to_physical_flux(
            E_raw, shape_raw, TOTAL_FLUX[src]
        )

        # Extend down to 1 keV
        E_src, dphi_dE_src = extend_continuum_to_low_energy(
            E_src, dphi_dE_src, E_min_new=1e-3
        )

        # Renormalize so total integrated flux stays exact
        dphi_dE_src = renormalize_differential_flux(
            E_src, dphi_dE_src, TOTAL_FLUX[src]
        )

        component_fluxes[src] = (E_src, dphi_dE_src)
        total_flux_grid += add_continuum_to_grid(E_src, dphi_dE_src, E_grid)

    # ============================================================
    # Add 7Be main line shape if available
    # ============================================================
    use_be7_lineshape = "7Be_lineshape" in spectra

    if use_be7_lineshape:
        E_raw, shape_raw = spectra["7Be_lineshape"]

        E_src, dphi_dE_src = normalize_shape_to_physical_flux(
            E_raw, shape_raw, TOTAL_FLUX["7Be_total"] * 0.8956
        )

        # Put it on the same common grid down to 1 keV.
        # Below the first tabulated point it should stay zero,
        # because this is not a continuum source extending to zero.
        component_fluxes["7Be_lineshape"] = (E_src, dphi_dE_src)
        total_flux_grid += add_continuum_to_grid(E_src, dphi_dE_src, E_grid)

        # add only the low-energy 7Be line as top-hat
        be7_low = add_line_to_grid(
            E_grid,
            E0=0.3843,
            total_flux=TOTAL_FLUX["7Be_total"] * 0.1044,
            bin_width=LINE_BIN_WIDTH_MEV
        )
        total_flux_grid += be7_low
    else:
        for Eline, frac in BE7_LINES:
            total_flux_grid += add_line_to_grid(
                E_grid,
                E0=Eline,
                total_flux=TOTAL_FLUX["7Be_total"] * frac,
                bin_width=LINE_BIN_WIDTH_MEV
            )

    # ============================================================
    # Add pep line
    # ============================================================
    total_flux_grid += add_line_to_grid(
        E_grid,
        E0=PEP_LINE_ENERGY_MEV,
        total_flux=TOTAL_FLUX["pep"],
        bin_width=LINE_BIN_WIDTH_MEV
    )

    # ============================================================
    # Save total flux CSV
    # ============================================================
    save_total_flux_file(
        "solar_neutrino_total_flux.csv",
        E_grid,
        total_flux_grid
    )

    # ============================================================
    # Plot 1: all components
    # ============================================================
    fig, ax = plt.subplots(figsize=(12, 8))

    for src in ["pp", "13N", "15O", "17F", "hep", "8B"]:
        if src in component_fluxes:
            E_src, dphi_dE_src = component_fluxes[src]
            positive = dphi_dE_src > 0
            ax.plot(E_src[positive], dphi_dE_src[positive], lw=2.0, label=src)

    if use_be7_lineshape:
        E_src, dphi_dE_src = component_fluxes["7Be_lineshape"]
        positive = dphi_dE_src > 0
        ax.plot(E_src[positive], dphi_dE_src[positive], lw=2.0, label=r"$^7$Be lineshape (861 keV)")

        be7_low_height = TOTAL_FLUX["7Be_total"] * 0.1044 / LINE_BIN_WIDTH_MEV
        ax.vlines(
            0.3843,
            YMIN,
            be7_low_height,
            linestyles=":",
            linewidth=2.2,
            label=r"$^7$Be (0.3843 MeV)"
        )
    else:
        for Eline, frac in BE7_LINES:
            height = TOTAL_FLUX["7Be_total"] * frac / LINE_BIN_WIDTH_MEV
            ax.vlines(
                Eline,
                YMIN,
                height,
                linestyles=":",
                linewidth=2.2,
                label=fr"$^7$Be ({Eline:.4f} MeV)"
            )

    pep_height = TOTAL_FLUX["pep"] / LINE_BIN_WIDTH_MEV
    ax.vlines(
        PEP_LINE_ENERGY_MEV,
        YMIN,
        pep_height,
        linestyles="--",
        linewidth=2.2,
        label="pep"
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-3, 20.0)
    ax.set_ylim(YMIN, 2e13)

    ax.set_xlabel(r"Neutrino energy $E_\nu$ [MeV]", fontsize=14)
    ax.set_ylabel(r"$d\Phi/dE_\nu$ [cm$^{-2}$ s$^{-1}$ MeV$^{-1}$]", fontsize=14)
    ax.set_title("Solar neutrino fluxes at Earth", fontsize=15)

    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=10, ncol=2)
    fig.tight_layout()
    fig.savefig(OUTFILE, dpi=300)
    plt.close(fig)

    # ============================================================
    # Plot 2: total flux only
    # ============================================================
    fig2, ax2 = plt.subplots(figsize=(12, 8))

    positive = total_flux_grid > 0
    ax2.plot(E_grid[positive], total_flux_grid[positive], lw=2.2, label="Total solar neutrino flux")

    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(1e-3, 20.0)
    ax2.set_ylim(YMIN, 2e13)

    ax2.set_xlabel(r"Neutrino energy $E_\nu$ [MeV]", fontsize=14)
    ax2.set_ylabel(r"$d\Phi/dE_\nu$ [cm$^{-2}$ s$^{-1}$ MeV$^{-1}$]", fontsize=14)
    ax2.set_title("Total solar neutrino flux at Earth", fontsize=15)

    ax2.grid(True, which="both", alpha=0.25)
    ax2.legend(fontsize=10)
    fig2.tight_layout()
    fig2.savefig("solar_neutrino_total_flux_loglog.png", dpi=300)
    plt.close(fig2)

    print(f"Saved {OUTFILE}")
    print("Saved solar_neutrino_total_flux_loglog.png")
    print("Saved solar_neutrino_total_flux.csv")
    print(f"All continuum sources extended down to 1 keV.")
    print(f"Line bin width used for pep and 7Be: {LINE_BIN_WIDTH_MEV} MeV")

if __name__ == "__main__":    main()