# Tools for Directional Neutrinos

## `EnergywithPerformance.py`

This script compares directional and non-directional detectors for elastic scattering:

\[
\nu_e + e^- \rightarrow \nu_e + e^-
\]

It makes 2 simple plots:

1. `enu_vs_te_different_theta.png`  
   - neutrino energy vs measured electron energy  
   - for different measured recoil angles

2. `enu_vs_theta_fixed_te.png`  
   - neutrino energy vs measured recoil angle  
   - for fixed measured electron energies

### Detector models

The code uses a small detector-performance object.

For each detector you can choose:
- threshold
- energy resolution function
- angular resolution function

Example:

- **Borexino-like**
  - threshold = 80 keV
  - energy resolution = 5% / sqrt(E[MeV])

- **Directional**
  - threshold = 10 keV
  - energy resolution = 20%
  - angular resolution = 20 degrees

### How to change detector performance

Inside `main()` you can modify:

```python
borexino_like = DetectorPerformance(...)
directional = DetectorPerformance(...)
```

---

## `downladFluxes.py`

This script downloads tabulated solar-neutrino source spectra from the Bahcall SNdata archive and saves them into:

`solar_neutrino_tables/`

Downloaded continuum sources:
- pp
- 8B
- hep
- 13N
- 15O
- 17F

It also writes a small `manifest.csv` file with:
- source name
- filename
- file size
- download status
- source URL

### Notes

- `pep` is not downloaded because it is a line source, not a continuum spectrum.
- `7Be` is also a line-like source and is handled separately from the continuum files.

---

## `plotFluxes.py`

This script reads the downloaded solar-neutrino tables from `solar_neutrino_tables/`, converts the tabulated spectral shapes into physical differential fluxes in:

`cm^-2 s^-1 MeV^-1`

and makes log-log plots of the solar neutrino fluxes at Earth.

Included components:
- continuum sources: pp, 8B, hep, 13N, 15O, 17F
- pep as a finite-width line
- 7Be as line components
- optional `be7_lineshape.dat` if available

### What the script does

The script:
- normalizes each tabulated shape to the corresponding total solar flux
- extends continuum spectra down to 1 keV
- builds a common energy grid
- sums all components into a total solar-neutrino flux
- saves plots and a CSV file with the total flux

### Main outputs

- `solar_neutrino_fluxes_loglog_physical.png`
- `solar_neutrino_total_flux_loglog.png`
- `solar_neutrino_total_flux.csv`

### Units

The final differential fluxes are saved and plotted in:

\[
\mathrm{cm^{-2}\,s^{-1}\,MeV^{-1}}
\]

For line sources like pep and 7Be, the script uses a finite bin width to represent them on the plot and in the saved total-flux file.

