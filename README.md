# CYGNUS-K Directional Neutrino Research Sandbox

This repository is a daily-use research sandbox for directional neutrino studies
in gas detectors. It is not a polished public package. The current code supports
first-pass calculations for solar-neutrino and supernova-neutrino elastic
scattering on electrons, then scales those single-electron calculations to
simple gas-target detector configurations.

The repository is meant to stay simple: root-level scripts are the runnable
entry points, while small shared modules hold repeated detector logic. Generated
tables and plots may exist locally but are mostly ignored by git.

## Purpose

The scientific workflow addressed here is:

1. Build solar-neutrino source spectra and flavor-separated fluxes at Earth.
2. Compute Standard Model shaped neutrino-electron scattering observables.
3. Convert single-electron rates into accepted detector rates for gas targets.
4. Estimate reconstructed neutrino-energy spectra with configurable recoil
   energy and angular resolution.
5. Repeat the detector-rate and reconstruction logic for time-integrated
   supernova fluence spectra.
6. Use Garfield/Magboltz gas transport outputs and stopping-power tables to set
   detector-motivated recoil windows.

The calculations are useful for detector-design comparisons and research
iteration. They are not a precision oscillation or detector-response framework.
Unclear or under-referenced physics choices are intentionally documented as
validation gaps rather than hidden.

## Quick Start

The user workflow assumed by this cleanup is WSL:

```bash
cd /mnt/c/Users/david/MyDrive/WORK/CYGNUS-K
python3 -m pip install -r requirements.txt
```

The main solar chain is:

```bash
python3 plotFluxes.py
python3 scatteringPlots.py
python3 gasTargetRates.py
python3 recoNuEnergyComparison.py
```

The main supernova chain is:

```bash
python3 generate_sn_fluence.py
python3 supernovaGasTargetEvents.py
python3 supernovaRecoNuEnergyComparison.py
```

To refresh the raw solar-neutrino input tables from the Bahcall SNdata archive:

```bash
python3 downloadFluxes.py
```

Typical output directories are:

- `solar_neutrino_fluxes/`
- `solar_nu_e_plots_by_flavor/`
- `solar_nu_gas_target_rates/`
- `solar_nu_reco_energy_comparison/`
- `supernova_fluence/`
- `supernova_nu_gas_target_events/`
- `supernova_nu_reco_energy_comparison/`

Most generated CSV and PNG products are ignored by `.gitignore`. Rerun the
workflow after changing configs, gas tables, or physics helpers.

## Environment

Required for the Python analysis scripts:

- Python 3
- `numpy`
- `pandas`
- `matplotlib`
- `requests`, only for downloading solar tables

Installable Python requirements are listed in `requirements.txt`.

Optional external dependencies:

- ROOT and Garfield++/Magboltz are required only for the scripts in
  `GarfieldSim/`.
- The main solar and supernova Python analysis chain does not require ROOT or
  Garfield at runtime if the derived detector tables already exist.

## Repository Layout

- `README.md`: this practical repository guide.
- `requirements.txt`: minimal Python dependencies for the non-Garfield scripts.
- `detector_model.py`: shared detector/gas helper module. It contains detector
  geometry loading, gas electron-density conversion, recoil-window lookup,
  diffusion-threshold adjustment, stopping-power range interpolation, labels,
  slugs, and the small progress bar.
- `plotFluxes.py`: builds solar-neutrino source and flavor fluxes.
- `scatteringPlots.py`: computes single-electron neutrino-electron scattering
  plots and is also used as the scattering helper module by later scripts.
- `gasTargetRates.py`: scales solar single-electron scattering rates to gas
  detector targets and writes accepted detector spectra.
- `recoNuEnergyComparison.py`: estimates reconstructed solar neutrino-energy
  spectra for energy-only and energy-plus-direction strategies.
- `generate_sn_fluence.py`: builds analytic, time-integrated supernova fluence
  spectra.
- `supernovaGasTargetEvents.py`: computes accepted supernova event-count
  spectra for gas targets.
- `supernovaRecoNuEnergyComparison.py`: estimates reconstructed supernova
  neutrino-energy count spectra.
- `EnergywithPerformance.py`: standalone kinematic/performance plotting utility.
- `downloadFluxes.py`: downloader entry point for solar-neutrino source tables.
- `detector_geometry.json`: default rectangular detector geometry.
- `reco_response_config.json`: default reconstruction response model.
- `DetectorNumbers/`: gas densities, stopping powers, derived range and
  diffusion tables, and small detector-number post-processing scripts.
- `GarfieldSim/`: Garfield/ROOT macros, Python helpers, gas tables, and local
  transport plots/tables.
- `Papers/`: local reference material already present in the repository.
- `CAD/`: detector/CAD side material, not used by the Python analysis chain.

## Current Default Configuration

`detector_geometry.json` defines a rectangular fiducial volume:

```text
length = 2 m
width  = 10 m
height = 10 m
V      = 200 m^3 = 2.0e8 cm^3
```

`DetectorNumbers/GasDensities.csv` currently lists seven gas targets:

- CF4 at 1, 3, and 10 atm
- HeCF4(60% He 40% CF4) at 1 atm
- HeCF4CH4(58% He 37% CF4 5% CH4) at 1, 0.5, and 0.1 atm

`reco_response_config.json` currently uses:

```text
energy_resolution_at_threshold_frac = 0.3
energy_resolution_constant_frac     = 0.1
angular_resolution_at_threshold_deg = 20.0
angular_resolution_constant_deg     = 5.0
reco_energy_bins                    = 100
```

The response model combines a threshold-scaled term with a constant floor.

## Workflow Overview

### Solar Workflow

1. `downloadFluxes.py` optionally refreshes raw source spectral tables in
   `solar_neutrino_tables/`.
2. `plotFluxes.py` reads those source tables, normalizes them to hard-coded
   total solar-neutrino fluxes, applies a minimal MSW-LMA survival model, and
   writes `solar_neutrino_fluxes/solar_neutrino_total_flux.csv`.
3. `scatteringPlots.py` reads the flavor flux CSV and computes one-electron
   scattering cross sections, rate-per-energy-bin plots, and angular spectra.
4. `gasTargetRates.py` reads the solar flux CSV, gas densities, stopping-power
   thresholds, diffusion summary, and detector geometry. It computes accepted
   detector-level spectra in `E_nu`, `T_e`, and `cos(theta_e)`.
5. `recoNuEnergyComparison.py` uses the same accepted recoil windows and gas
   targets, then smears accepted true event rates into reconstructed neutrino
   energy bins for energy-only and energy-plus-direction assumptions.

### Supernova Workflow

1. `generate_sn_fluence.py` writes analytic time-integrated fluence spectra in
   `supernova_fluence/`.
2. `supernovaGasTargetEvents.py` applies the gas-target recoil-window workflow
   to the fluence inputs. Since the fluence is already time-integrated, outputs
   are counts, not rates.
3. `supernovaRecoNuEnergyComparison.py` mirrors the solar reconstruction script,
   but smears accepted burst counts instead of steady-state rates.

### Detector-Number Workflow

1. `GarfieldSim/read.py` loads Garfield `.gas` files and writes selected
   transport points to `GarfieldSim/plots/*_transport_points.csv`.
2. `DetectorNumbers/computeDiffusion.py` converts Garfield diffusion
   coefficients into 2-sigma widths over 50 cm drift.
3. `DetectorNumbers/computePlotRanges.py` integrates stopping-power tables to
   derive electron energies for 1 mm and 1 m ranges.
4. `detector_model.py` combines the 1 mm range threshold with diffusion
   information to set the accepted lower recoil threshold for each gas.

## Physics And Mathematical Model

### Solar Source Fluxes

For each continuum source, `plotFluxes.py` treats the input table as a spectral
shape `s_i(E)` and normalizes it to a total source flux `Phi_i`:

```math
\frac{d\Phi_i}{dE}(E)
=
\Phi_i
\frac{s_i(E)}{\int s_i(E')\,dE'} .
```

Units:

- `E` is in MeV.
- `Phi_i` is in cm^-2 s^-1.
- `dPhi_i/dE` is in cm^-2 s^-1 MeV^-1.

Line sources are represented as finite top-hat bins with width
`LINE_BIN_WIDTH_MEV = 0.002`:

```math
\frac{d\Phi_{\rm line}}{dE} =
\frac{\Phi_{\rm line}}{\Delta E}
```

inside the line bin and zero outside. The current checked-in table set does not
include `be7_lineshape.dat`, so the default code treats both 7Be branches as
finite-width line bins.

For continuum tables whose first positive point lies above the global grid
minimum, the code extends to `1e-3 MeV` with a matched low-energy law:

```math
\frac{d\Phi}{dE} \propto E^2 .
```

This is a practical interpolation/extrapolation choice and should be validated
for precision work.

### Solar Oscillation Approximation

`plotFluxes.py` uses a compact adiabatic MSW-LMA approximation:

```math
P_{ee}^{3\nu}
=
\sin^4\theta_{13}
+
\cos^4\theta_{13} P_{ee}^{2\nu},
```

with

```math
P_{ee}^{2\nu}
=
\frac{1}{2}
\left[
1 + \cos(2\theta_{12})\cos(2\theta_{12}^m)
\right],
```

and

```math
\cos(2\theta_{12}^m)
=
\frac{\cos(2\theta_{12})-\beta}
{\sqrt{(\cos(2\theta_{12})-\beta)^2+\sin^2(2\theta_{12})}},
\qquad
\beta = \frac{A}{\Delta m_{21}^2}.
```

The matter potential is implemented as

```math
A[{\rm eV}^2] = 1.52\times 10^{-7}
n_e[{\rm mol\,cm^{-3}}] E[{\rm MeV}].
```

Each solar source is assigned one effective production-region electron density
from `NE_MOL_CM3` in `plotFluxes.py`. The code does not solve propagation
through a full solar model. After computing `P_ee`, flavor splitting is:

```math
\phi_e = P_{ee}\phi_{\rm source},
\qquad
\phi_\mu = \phi_\tau =
\frac{1-P_{ee}}{2}\phi_{\rm source}.
```

### Neutrino-Electron Cross Sections

`scatteringPlots.py` defines simple linear total cross-section approximations:

```math
\sigma(E_\nu) = c_\alpha E_\nu ,
```

with `E_nu` in MeV and `sigma` in cm^2 [Supernova Neutrino Detection in a liquid Argon TPC](https://arxiv.org/pdf/hep-ph/0307222) and [Supernova Relic Neutrinos in Liquid Argon detectors](https://arxiv.org/pdf/hep-ph/0408031):

| channel | coefficient `c_alpha` [cm^2 MeV^-1] |
| --- | ---: |
| `nu_e` | `9.20e-45` |
| `anti_nu_e` | `3.83e-45` |
| `nu_x` (`nu_mu`, `nu_tau`) | `1.57e-45` |
| `anti_nu_x` | `1.29e-45` |

The differential recoil-energy shape is the tree-level Standard Model form:

```math
\frac{d\sigma}{dT}
=
\frac{2G_F^2m_e}{\pi}
\left[
g_L^2
+ g_R^2\left(1-\frac{T}{E_\nu}\right)^2
- g_Lg_R\frac{m_eT}{E_\nu^2}
\right],
```

with the code applying the GeV-to-cm^2 and GeV-to-MeV unit conversions.

The couplings are:

```text
nu_e      : g_L =  1/2 + sin^2(theta_W), g_R =  sin^2(theta_W)
anti_nu_e : g_L =  sin^2(theta_W),       g_R =  1/2 + sin^2(theta_W)
nu_x      : g_L = -1/2 + sin^2(theta_W), g_R =  sin^2(theta_W)
anti_nu_x : g_L =  sin^2(theta_W),       g_R = -1/2 + sin^2(theta_W)
```

Kinematic limits and conversions use:

```math
T_{\max}(E_\nu)
=
\frac{2E_\nu^2}{m_e+2E_\nu},
```

```math
E_{\nu,\min}(T)
=
\frac{1}{2}
\left[
T+\sqrt{T^2+2m_eT}
\right],
```

and for recoil angle `theta_e`,

```math
T(E_\nu,\cos\theta_e)
=
\frac{2m_eE_\nu^2\cos^2\theta_e}
{(E_\nu+m_e)^2-E_\nu^2\cos^2\theta_e}.
```

The angular spectra are formed from the recoil-energy shape and the Jacobian
`dT/dcos(theta_e)`, then rescaled so the integral matches the approximate
linear total cross section. This hybrid normalization is an explicit repository
convention.

For accepted gas-target rates, `gasTargetRates.py` integrates the SM
differential shape inside the accepted recoil window and applies that fraction
to the approximate total cross section:

```math
\sigma_\alpha^{\rm acc}(E_\nu)
=
\sigma_\alpha^{\rm approx}(E_\nu)
\frac{
\int_{T_{\min}}^{\min(T_{\max},T_{\max}^{\rm gas})}
(d\sigma_\alpha^{\rm SM}/dT)\,dT
}{
\int_0^{T_{\max}}(d\sigma_\alpha^{\rm SM}/dT)\,dT
}.
```

### Gas Target Normalization

`detector_model.py` treats gas-mixture fractions as molecular fractions. For
species `i` with fraction `x_i`, molecular electron count `Z_i`, and molar mass
`M_i`,

```math
\langle Z\rangle = \sum_i x_i Z_i,
\qquad
\langle M\rangle = \sum_i x_i M_i .
```

Given gas density `rho` in g/cm^3,

```math
n_e =
\rho
\frac{N_A}{\langle M\rangle}
\langle Z\rangle ,
```

and for detector volume `V` in cm^3,

```math
N_e = n_e V .
```

The solar gas-target differential spectra are:

```math
\frac{dR}{dE_\nu}
=
N_e\sum_\alpha
\phi_\alpha(E_\nu)\sigma_\alpha^{\rm acc}(E_\nu),
```

```math
\frac{dR}{dT}
=
N_e\sum_\alpha\int dE_\nu\,
\phi_\alpha(E_\nu)
\frac{d\sigma_\alpha}{dT}(E_\nu,T),
```

and

```math
\frac{dR}{d\cos\theta_e}
=
N_e\sum_\alpha\int dE_\nu\,
\phi_\alpha(E_\nu)
\frac{d\sigma_\alpha}{d\cos\theta_e}.
```

Solar flux units make `R` a rate in s^-1. Supernova fluence units make the
corresponding output a count.

### Recoil Windows

The accepted recoil window for each gas is derived from
`DetectorNumbers/electron_range_energy_table.csv` and
`DetectorNumbers/diffusion_2sigma_50cm_summary.csv` source [NIST](https://www.nist.gov/pml/stopping-power-range-tables-electrons-protons-and-helium-ions).

The baseline lower threshold is the electron kinetic energy whose CSDA-like
range is 1 mm. The upper threshold is the energy whose range is 1 m. The lower
threshold is then raised if the largest longitudinal 2-sigma diffusion value at
electric fields up to 2 kV/cm is larger than 1 mm:

```math
T_{\min}^{\rm gas}
=
T(R=\max[1\,{\rm mm},D_L^{2\sigma}]),
\qquad
T_{\max}^{\rm gas}=T(R=1\,{\rm m}).
```

The stopping-power range estimate used by `DetectorNumbers/computePlotRanges.py`
is:

```math
R(E) = \int_0^E
\frac{dE'}{\rho S_{\rm col}(E')},
```

where `S_col` is the collisional mass stopping power in MeV cm^2/g. The code
uses a simple low-energy `1/E` stopping-power extrapolation below the first
tabulated point.

The diffusion summary uses Garfield diffusion coefficients `D` in sqrt(cm):

```math
2\sigma = 2D\sqrt{L},
```

with `L = 50 cm`. Outputs are also converted to mm.

### Reconstruction Model

The energy-only estimator uses the minimum neutrino energy compatible with a
recoil kinetic energy:

```math
E_{\nu,\min}(T)
=
\frac{1}{2}
\left[
T+\sqrt{T^2+2m_eT}
\right].
```

The directional estimator uses measured recoil kinetic energy and angle:

```math
E_\nu(T,\theta)
=
\frac{m_eT}{p_e\cos\theta - T},
\qquad
p_e = \sqrt{T^2+2m_eT}.
```

The response model in `reco_response_config.json` is:

```math
\frac{\sigma_E}{E}
=
\sqrt{
\left(r_E\sqrt{\frac{T_{\rm thr}}{T}}\right)^2
+ r_{E,0}^2
},
```

```math
\sigma_\theta
=
\sqrt{
\left(r_\theta\sqrt{\frac{T_{\rm thr}}{T}}\right)^2
+ r_{\theta,0}^2
}.
```

`T_thr` is the gas-dependent lower accepted recoil threshold. The scripts
propagate these uncertainties through the estimators and smear accepted true
rates or counts into logarithmic reconstructed-energy bins. This is a fast
Gaussian propagated-resolution model, not a full event-level detector
simulation.

### Supernova Fluence

`generate_sn_fluence.py` implements the quasi-thermal fluence parameterization
documented in its comments as coming from the [Neutrino flux sensitivity to the next galactic core-collapse supernova in COSINUS](https://iopscience.iop.org/article/10.1088/1475-7516/2025/03/037)
sensitivity study:

```math
\frac{d\Phi_i}{dE}
=
\frac{\epsilon_i}{4\pi d^2}
\frac{E^\alpha e^{-E/T_i}}
{T_i^{\alpha+2}\Gamma(\alpha+2)},
\qquad
T_i = \frac{\langle E_i\rangle}{\alpha+1}.
```

Units:

- `epsilon_i` is emitted energy per species, converted from erg to MeV.
- `d` is distance in cm.
- `dPhi_i/dE` is in cm^-2 MeV^-1.

The current implementation contains only `SN1987A_like`:

- distance `10 kpc`
- `alpha = 3`
- total emitted energy `3e53 erg`, equally shared among six species
- average energies: `nu_e = 9 MeV`, `anti_nu_e = 12 MeV`, heavy flavors
  `16 MeV`

The numerical 27 solar-mass model mentioned in the code comments is not
implemented because the repository does not contain its simulation table or a
SNEWPY/Nakazato/Garching input. Adding such a model should be table-driven.

Because the supernova input is time-integrated fluence, the detector outputs are
counts:

```math
\Phi(E_\nu)\ [{\rm cm}^{-2}{\rm MeV}^{-1}]
\longrightarrow
N\ [{\rm counts}].
```

The optional burst duration in `supernovaGasTargetEvents.py` is used only for an
average-rate diagnostic `N / Delta t`. It is not multiplied into the event
count.

## Code And Tool Reference

| File | Purpose | Main inputs | Main outputs | Example |
| --- | --- | --- | --- | --- |
| `downloadFluxes.py` | Download solar source tables | remote Bahcall SNdata URLs | `solar_neutrino_tables/`, `manifest.csv` | `python3 downloadFluxes.py` |
| `plotFluxes.py` | Build solar source and flavor fluxes | `solar_neutrino_tables/` | `solar_neutrino_fluxes/solar_neutrino_total_flux.csv`, flux plots | `python3 plotFluxes.py` |
| `scatteringPlots.py` | One-electron nu-e scattering plots | solar flux CSV | `solar_nu_e_plots_by_flavor/` | `python3 scatteringPlots.py` |
| `detector_model.py` | Shared detector/gas utilities | JSON configs, gas tables, range and diffusion summaries | imported helpers only | imported by driver scripts |
| `gasTargetRates.py` | Solar gas-target rates | flux CSV, gas density, range, diffusion, geometry | `solar_nu_gas_target_rates/` | `python3 gasTargetRates.py` |
| `recoNuEnergyComparison.py` | Solar reconstructed spectra | gas workflow inputs plus response config | `solar_nu_reco_energy_comparison/` | `python3 recoNuEnergyComparison.py` |
| `generate_sn_fluence.py` | Analytic supernova fluence generation | built-in model parameters | `supernova_fluence/` | `python3 generate_sn_fluence.py` |
| `supernovaGasTargetEvents.py` | Supernova gas-target event counts | fluence CSVs, gas workflow inputs | `supernova_nu_gas_target_events/` | `python3 supernovaGasTargetEvents.py` |
| `supernovaRecoNuEnergyComparison.py` | Supernova reconstructed count spectra | fluence CSVs, gas workflow inputs, response config | `supernova_nu_reco_energy_comparison/` | `python3 supernovaRecoNuEnergyComparison.py` |
| `EnergywithPerformance.py` | Standalone kinematic performance comparison | built-in example detector models | `enu_vs_te_different_theta.png`, `enu_vs_theta_fixed_te.png` | `python3 EnergywithPerformance.py` |
| `DetectorNumbers/computePlotRanges.py` | Electron range and threshold table | gas densities, stopping powers | `electron_range_energy_table.csv`, `electron_ranges.png` | `python3 DetectorNumbers/computePlotRanges.py` |
| `DetectorNumbers/computeDiffusion.py` | Diffusion summary over 50 cm | Garfield transport-point CSVs | `diffusion_2sigma_50cm_summary.csv`, comparison plot | `python3 DetectorNumbers/computeDiffusion.py` |
| `GarfieldSim/read.py` | Read `.gas` files and export transport summaries | `GarfieldSim/tables/*.gas`, ROOT/Garfield | `GarfieldSim/plots/` | run inside a Garfield environment |

## Validation And Sanity Checks

There is not yet a formal test suite. Useful checks are:

```bash
python3 -m py_compile detector_model.py *.py DetectorNumbers/*.py GarfieldSim/*.py
python3 DetectorNumbers/computePlotRanges.py
python3 DetectorNumbers/computeDiffusion.py
python3 plotFluxes.py
python3 scatteringPlots.py
python3 gasTargetRates.py --t-grid-points 80 --costh-grid-points 80 --outdir /tmp/cygnus_gas_smoke
python3 recoNuEnergyComparison.py --t-grid-points 40 --max-true-energy-points 80 --outdir /tmp/cygnus_reco_smoke
python3 generate_sn_fluence.py --output-root /tmp/cygnus_sn_fluence
python3 supernovaGasTargetEvents.py --fluence-root /tmp/cygnus_sn_fluence --output-root /tmp/cygnus_sn_events_smoke --t-grid-points 80 --costh-grid-points 80
python3 supernovaRecoNuEnergyComparison.py --fluence-root /tmp/cygnus_sn_fluence --output-root /tmp/cygnus_sn_reco_smoke --t-grid-points 40 --max-true-energy-points 80
python3 EnergywithPerformance.py
```

The reduced grid commands are intended as smoke tests. Production plots should
use the script defaults unless there is a reason to trade precision for speed.

Internal consistency checks currently include:

- input-column validation for the gas, recoil-window, diffusion, geometry, and
  response config files;
- positivity checks for detector geometry and response parameters;
- rate/count comparisons between the accepted `E_nu` sum and the numerical
  `T_e` or `cos(theta_e)` integrals in the gas-target summaries;
- explicit tracking of the threshold source for each gas.

## Extending The Repository

Add a new solar source by:

1. Adding or downloading its source spectrum into `solar_neutrino_tables/`.
2. Adding a parser or using an existing parser in `plotFluxes.py`.
3. Adding its total flux to `TOTAL_FLUX`.
4. Adding an effective production density to `NE_MOL_CM3` if the MSW-LMA helper
   should process it.
5. Rerunning `plotFluxes.py` and downstream scripts.

Add a new gas target by:

1. Adding a row to `DetectorNumbers/GasDensities.csv`.
2. Adding molecular composition data to `GAS_MIXTURES` in `detector_model.py`.
3. Adding a stopping-power table and `STOPPING_POWER_FILES` entry.
4. Producing or adding Garfield transport-point CSVs for diffusion matching.
5. Rerunning `DetectorNumbers/computePlotRanges.py`,
   `DetectorNumbers/computeDiffusion.py`, and the gas-target scripts.

Add a new detector response by editing or copying `reco_response_config.json`.
Keep all resolution terms non-negative and choose enough reconstructed-energy
bins for the intended comparison.

Add a new supernova model by writing a fluence CSV with:

- `E_MeV`
- `fluence_nue_cm-2_MeV-1`
- `fluence_anue_cm-2_MeV-1`
- either single-heavy-species columns
  `fluence_nux_single_cm-2_MeV-1` and
  `fluence_anux_single_cm-2_MeV-1`, or grouped heavy-flavor columns
  `fluence_nux_group_numu_nutau_cm-2_MeV-1` and
  `fluence_anux_group_anumu_anutau_cm-2_MeV-1`

Then run the supernova scripts with:

```bash
python3 supernovaGasTargetEvents.py --fluence-root PATH --models all
python3 supernovaRecoNuEnergyComparison.py --fluence-root PATH --models all
```

When adding a new cross-section or detector-response model, prefer introducing
one small shared helper rather than adding another copy of the same formula to a
driver script.

## References Present In The Repository

Only references already present or explicit in the code/comments are listed
here:

- `Papers/SuperK-2003.pdf` is stored locally.
- `downloadFluxes.py` points to the Bahcall SNdata archive URLs used for solar
  source tables.
- `generate_sn_fluence.py` cites [Neutrino flux sensitivity to the next galactic core-collapse supernova in COSINUS](https://iopscience.iop.org/article/10.1088/1475-7516/2025/03/037), for the quasi-thermal fluence
  parameterization and SN1987A-like benchmark.
- `GarfieldSim/` assumes Garfield++/Magboltz gas tables and transport APIs.

Documentation gaps that require expert validation:

- The hard-coded total solar source flux values in `plotFluxes.py` need an
  explicit local reference.
- The effective production-region electron densities `NE_MOL_CM3` need an
  explicit source or derivation.
- The approximate total nu-e cross-section coefficients in `scatteringPlots.py`
  need an explicit source.
- The stopping-power input provenance should be documented beyond the file
  names.
- The detector response parameter choices in `reco_response_config.json` are
  assumptions, not validated detector performance results.

## Known Limitations And WIP Notes

- The solar oscillation treatment is a compact effective-density MSW-LMA model,
  not full solar propagation.
- The current 7Be default uses finite-width top-hat lines unless a
  `be7_lineshape.dat` file is supplied.
- The cross-section convention mixes a Standard Model differential shape with
  approximate linear total cross sections. This is intentional for continuity
  with earlier scripts but should be reviewed before precision use.
- The gas target rates include only recoil-window acceptance. They do not include
  detector efficiency, backgrounds, topology cuts, trigger behavior, or
  time-dependent detector conditions.
- Reconstruction uses propagated Gaussian uncertainties and coarse bins, not a
  full detector simulation or event-level inference.
- Gas mixture fractions are treated as molecular/volume fractions.
- Range thresholds are CSDA-like estimates with a simple low-energy
  extrapolation and no straggling model.
- The supernova model set currently contains only one analytic SN1987A-like
  benchmark.
- Garfield scripts require an external ROOT/Garfield environment and are not
  covered by the ordinary Python smoke tests.
- Generated outputs in local ignored folders may be stale relative to code or
  config changes. Rerun the workflow when in doubt.
