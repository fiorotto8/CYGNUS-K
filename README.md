# CYGNUS-K Neutrino Gas-Detector Sandbox

This repository is a research sandbox for first-pass CYGNUS-K neutrino studies
in gas detectors. It connects solar neutrino fluxes, neutrino-electron
scattering, gas-target normalization, supernova fluences, reconstructed
neutrino-energy spectra, and optional CEvNS nuclear-recoil calculations.

The code is intended to make detector and target-scaling questions easy to
inspect. It is not a full detector simulation, and several physics ingredients
are deliberately simple or pending validation.

## What This Repository Computes

The supported workflows are:

- Solar neutrino flux generation and flavor splitting at Earth.
- Standard Model neutrino-electron scattering on gas targets.
- Gas-target accepted rates using gas densities, detector volume, electron
  ranges, and diffusion-derived recoil-window adjustments.
- Supernova neutrino fluence generation with analytic quasi-thermal spectra.
- Supernova gas-target event counts for neutrino-electron scattering.
- Reconstructed neutrino-energy spectra for neutrino-electron scattering.
- CEvNS nuclear-recoil rates and supernova counts for the explicit gas
  constituents `12C`, `19F`, `1H`, and `4He`.
- Optional CEvNS energy-only lower-bound spectra using `E_nu_min(T_N)`.

The main workflows remain simple root-level commands such as:

```bash
python3 plotFluxes.py
python3 gasTargetRates.py --enable-cevns
python3 supernovaGasTargetEvents.py --enable-cevns
```

The cleanup keeps those commands as the stable user interface and centralizes
shared CEvNS logic in helper modules instead of duplicating it in each driver.

## Repository Structure

Root workflow scripts:

- `plotFluxes.py`: builds solar source spectra, applies the compact solar
  oscillation/flavor split, and writes flavor-separated flux CSVs.
- `scatteringPlots.py`: computes and plots one-electron neutrino-electron
  scattering spectra from the solar flux CSV.
- `gasTargetRates.py`: computes solar gas-target neutrino-electron rates and,
  when enabled, gas-total CEvNS rates.
- `recoNuEnergyComparison.py`: computes solar reconstructed neutrino-energy
  spectra for neutrino-electron scattering and optional CEvNS lower-bound
  spectra.
- `generate_sn_fluence.py`: generates analytic time-integrated supernova
  fluence spectra.
- `supernovaGasTargetEvents.py`: computes supernova gas-target
  neutrino-electron event counts and optional gas-total CEvNS counts.
- `supernovaDistanceGasScan.py`: scans supernova neutrino-electron plus CEvNS
  counts/spectra from 0.5 to 10 kpc for every gas row by default, or for one
  selected gas row with `--gas`.
- `supernovaRecoNuEnergyComparison.py`: computes supernova reconstructed
  neutrino-energy count spectra and optional CEvNS lower-bound count spectra.

Shared modules:

- `detector_model.py`: detector geometry, gas density parsing, mixture
  normalization, electron target counts, CEvNS isotope target counts, electron
  range windows, diffusion adjustments, progress bars, and plotting constants.
- `cevns.py`: CEvNS isotope definitions, weak charges, Helm/pointlike vector
  form factors, the `19F` HMS axial response, and vector plus pure-axial
  differential cross sections.
- `cevns_pipeline.py`: CEvNS config loading, CLI options, active
  flux/fluence summing, recoil spectra, accepted `E_nu` spectra, gas-total
  aggregation, isotope diagnostic outputs, lower-bound reconstruction, and
  CEvNS plotting/saving helpers.
- `EnergywithPerformance.py`: neutrino-energy estimator and response helper
  functions used by the reconstructed-energy workflows.

Inputs and local references:

- `detector_geometry.json`: fiducial detector dimensions.
- `reco_response_config.json`: reconstructed-energy response settings and the
  CEvNS response block.
- `DetectorNumbers/GasDensities.csv`: gas composition, pressure, and density
  table.
- `DetectorNumbers/electron_range_energy_table.csv`: electron range thresholds
  produced by `DetectorNumbers/computePlotRanges.py`.
- `DetectorNumbers/diffusion_2sigma_50cm_summary.csv`: diffusion summary
  produced by `DetectorNumbers/computeDiffusion.py`.
- `solar_neutrino_tables/`: local solar source spectra used by `plotFluxes.py`.
- `Papers/CEvENS/`: local CEvNS references used for the implemented CEvNS
  conventions.
- `GarfieldSim/`: Garfield-related transport utilities and generated tables.

Generated output folders are ignored by git and can be regenerated:

- `solar_neutrino_fluxes/`
- `solar_nu_e_plots_by_flavor/`
- `solar_nu_gas_target_rates/`
- `solar_nu_reco_energy_comparison/`
- `supernova_fluence/`
- `supernova_nu_gas_target_events/`
- `supernova_nu_reco_energy_comparison/`

## Working Principle

### Solar Fluxes

`plotFluxes.py` loads local solar source spectra, normalizes each source to the
hard-coded source totals in the script, and interpolates them onto a common
energy grid. It writes:

```text
solar_neutrino_fluxes/solar_neutrino_total_flux.csv
```

The output columns include:

- `E_MeV`
- `phi_total_cm-2_s-1_MeV-1`
- `phi_nu_e_cm-2_s-1_MeV-1`
- `phi_nu_mu_cm-2_s-1_MeV-1`
- `phi_nu_tau_cm-2_s-1_MeV-1`

The flavor split is a compact MSW-LMA style survival-probability model. It is
adequate for this sandbox, but it is not a precision solar global-fit
oscillation calculation.

### Neutrino-Electron Scattering

`scatteringPlots.py`, `gasTargetRates.py`, and the reconstructed-energy scripts
use Standard Model neutrino-electron scattering cross sections. The gas-target
scripts integrate the flavor-separated flux or fluence over neutrino energy and
scale per-electron rates/counts by the number of target electrons in the gas
volume.

The electron-recoil workflow uses a gas-dependent accepted electron-recoil
window. The lower edge is based on the electron kinetic energy whose range is
1 mm, then raised when the diffusion summary indicates that the largest
longitudinal `2 sigma` diffusion value at low drift field exceeds 1 mm. The
upper edge is based on a 1 m range. This recoil window is specific to the
neutrino-electron analysis and is not reused for CEvNS.

### Gas Normalization

Gas normalization is handled in `detector_model.py`. The detector volume comes
from `detector_geometry.json` unless a workflow receives `--volume-cm3`.

For a gas entry with mass density `rho`, average molar mass `M`, and detector
volume `V`, the total number of gas mixture units is:

```text
N_mixture = rho * V * N_A / M
```

Electron targets are obtained from the electron count per mixture unit. CEvNS
isotope targets are obtained from the isotope count per mixture unit.

The gas mixtures currently supported in the target-count logic are:

- `CF4`: `1 x 12C + 4 x 19F`
- `CH4`: `1 x 12C + 4 x 1H`
- `He`: `1 x 4He`
- `HeCF4(60% He 40% CF4)`: `0.60 He + 0.40 CF4`
- `HeCF4CH4(58% He 37% CF4 5% CH4)`: `0.58 He + 0.37 CF4 + 0.05 CH4`

### Supernova Fluences

`generate_sn_fluence.py` writes time-integrated fluence spectra. The supernova
gas-target outputs are counts, not rates, because the input spectra have already
been integrated over the burst.

The supernova CEvNS convention sums all active neutrino and antineutrino
fluence columns produced by the current fluence CSV convention:

- `fluence_nue_cm-2_MeV-1`
- `fluence_anue_cm-2_MeV-1`
- either single heavy-flavor columns or grouped heavy-flavor columns, depending
  on the CSV.

### Reconstructed and Lower-Bound Spectra

For neutrino-electron scattering, `recoNuEnergyComparison.py` and
`supernovaRecoNuEnergyComparison.py` compare reconstructed neutrino-energy
spectra using the electron recoil energy and, where configured, angular
information.

For CEvNS there is no directionality in this repository. The CEvNS reconstructed
product is an energy-only lower-bound spectrum. It maps a nuclear recoil energy
`T_N` to the minimum neutrino energy compatible with that recoil:

```text
E_nu_min(T_N) = 1/2 * [T_N + sqrt(T_N^2 + 2 m_N T_N)]
```

with consistent units. This is not a true event-by-event neutrino-energy
reconstruction. It is a lower-bound estimator.

No CEvNS nuclear-recoil energy-resolution model is configured by default, so
the CEvNS lower-bound spectra are currently unsmeared. The code can use a CEvNS
energy-resolution block if one is added to `reco_response_config.json`, but no
validated detector model is supplied yet.

## CEvNS Model

CEvNS is implemented for the gas constituents explicitly:

| Element | Isotope | Z | N | A | J |
| --- | --- | ---: | ---: | ---: | ---: |
| C | `12C` | 6 | 6 | 12 | 0 |
| F | `19F` | 9 | 10 | 19 | 1/2 |
| H | `1H` | 1 | 0 | 1 | 1/2 |
| He | `4He` | 2 | 2 | 4 | 0 |

### Vector Term

The vector weak charge is:

```text
Q_W = Z * (1 - 4 sin^2(theta_W)) - N
```

Vector CEvNS is computed for all four isotopes. Hydrogen is intentionally kept:
because `N = 0`, its vector weak charge is small, but it is not silently
ignored.

The default vector form-factor mode is:

- `12C`: Helm
- `19F`: Helm
- `1H`: pointlike `F_W = 1`
- `4He`: pointlike `F_W = 1`

The global CLI option `--cevns-form-factor helm` applies the same Helm
placeholder expression to all isotopes for comparison. For `1H` and `4He`, Helm
is not treated as a precision model.

### Axial Term

The cross section is implemented as vector plus pure axial. Vector-axial
interference is not included; the local axial-vector CEvNS reference states
that these terms either vanish or are recoil-suppressed for the current
convention.

Axial status by isotope:

- `12C`: zero because `J = 0`.
- `4He`: zero because `J = 0`.
- `19F`: implemented with the Hoferichter-Menendez-Schwenk transverse
  Sigma-prime response and local fit coefficients.
- `1H`: not implemented for the axial/proton neutral-current elastic term.

Supported axial switches are:

- `none`: disable pure axial terms.
- `hoferichter_19f_fast`: current default `19F` HMS response.
- `hoferichter_19f_central`: currently aliases the same implemented `19F`
  central response convention.
- `toy`: `19F` only, using the zero-momentum response as a diagnostic; never
  used for hydrogen.

### Hydrogen Axial Review

Only the local files under `Papers/CEvENS/` were used for this review:

- `Papers/CEvENS/AccoppiamentoAssiale_Paper_2603.05281v1 (1).pdf`
- `Papers/CEvENS/HMS_axialFitparams.pdf`
- `Papers/CEvENS/HELM_vectorModel.pdf`

The axial-vector paper gives the vector plus pure-axial CEvNS convention used
here, discusses light nonzero-spin nuclei, and highlights fluorine compounds.
The HMS paper provides nuclear response formalism and fit coefficients for
selected nuclei, including `19F`. It does not provide a ready, validated
free-proton neutral-current elastic implementation in the same normalization as
the current `A > 1` CEvNS helper.

Therefore hydrogen axial/proton neutral-current elastic scattering is not
implemented because the current local references do not provide a complete,
validated implementation convention for the free-proton axial term in this
code. Hydrogen therefore contributes only through the suppressed vector term,
and CH4-rich mixtures are incomplete for proton-recoil studies.

The code must not approximate hydrogen axial scattering with:

- the `19F` HMS response;
- the `19F` toy model;
- any nuclear shell-model response table for `A > 1`.

### CEvNS Threshold

CEvNS uses a detector-level threshold interpreted for now as a true nuclear
recoil threshold:

```text
cevns.nr_threshold_keV
```

It can be overridden with:

```bash
--cevns-threshold-kev 1.0
```

An optional true nuclear-recoil maximum can be supplied with:

```bash
--cevns-max-kev VALUE
```

or left unset to integrate to the kinematic maximum. CEvNS does not reuse the
electron-recoil range/diffusion thresholds. Quenching and electron-equivalent
threshold conversion are not implemented. When quenching is added later, the
detector threshold must be mapped between electron-equivalent energy and true
nuclear-recoil energy; that mapping is intentionally outside the current code.

## Configuration

`detector_geometry.json` currently contains the default fiducial detector:

```json
{
  "name": "default_10x10x2m3_detector",
  "length_m": 2.0,
  "width_m": 10.0,
  "height_m": 10.0
}
```

`reco_response_config.json` contains electron-reconstruction response settings
and the optional CEvNS block:

```json
{
  "name": "default_threshold_scaled_response",
  "energy_resolution_at_threshold_frac": 0.3,
  "energy_resolution_constant_frac": 0.1,
  "angular_resolution_at_threshold_deg": 20.0,
  "angular_resolution_constant_deg": 5.0,
  "reco_energy_bins": 100,
  "cevns": {
    "enabled": true,
    "nr_threshold_keV": 1.0,
    "nr_max_keV": null,
    "form_factor": "default",
    "axial_model": "hoferichter_19f_fast",
    "reco_energy_bins": 100
  }
}
```

If the `cevns` block is absent, scripts use safe defaults:

- `enabled = false`
- `nr_threshold_keV = 1.0`
- `nr_max_keV = null`
- `form_factor = default`
- `axial_model = hoferichter_19f_fast`
- `reco_energy_bins = 100`

CLI options override config values:

```bash
--enable-cevns
--disable-cevns
--cevns-threshold-kev 1.0
--cevns-max-kev 10.0
--cevns-form-factor default|helm|pointlike
--cevns-axial-model none|hoferichter_19f_fast|hoferichter_19f_central|toy
--skip-cevns-plots
--skip-cevns-save
```

## Outputs

### Neutrino-Electron Outputs

Solar gas-target rates are written under `solar_nu_gas_target_rates/` by
default. Each gas subdirectory contains:

- `dR_dEnu.csv`: accepted rate density versus neutrino energy.
- `dR_dTe.csv`: accepted rate density versus electron recoil energy.
- `dR_dcosth.csv`: accepted angular spectrum.
- Matching PNG plots.

The root summary is:

```text
solar_nu_gas_target_rates/gas_target_rate_summary.csv
```

Supernova gas-target event counts are written under
`supernova_nu_gas_target_events/`. Each model/gas directory contains:

- `dN_dEnu.csv`
- `dN_dTe.csv`
- `dN_dcosth.csv`
- Matching PNG plots.

The root supernova summary is:

```text
supernova_nu_gas_target_events/all_models_event_summary.csv
```

Reconstructed neutrino-energy spectra are written as:

```text
reco_neutrino_energy_spectra.csv
reco_neutrino_energy_spectra.png
```

inside each gas subdirectory.

### CEvNS Outputs

CEvNS outputs live in a `cevns/` subdirectory of the corresponding workflow
output root.

The primary CEvNS summary is now gas-aggregated. There is one row per gas for
solar, or one row per model/gas for supernova. Isotope-level calculations are
still performed, but isotope summaries are diagnostic files.

Primary solar CEvNS summary:

```text
solar_nu_gas_target_rates/cevns/cevns_rate_summary.csv
```

Primary supernova CEvNS summary:

```text
supernova_nu_gas_target_events/cevns/all_models_cevns_event_summary.csv
```

Primary summary columns include:

- `gas_label`
- `gas`
- `isotopes`
- `number_of_targets_total`
- `nr_threshold_keV`
- `nr_max_keV`
- `form_factor`
- `axial_model`
- `vector_total_rate_s-1` or `vector_total_counts`
- `axial_total_rate_s-1` or `axial_total_counts`
- `total_rate_s-1` or `total_counts`
- `axial_fraction`
- `hydrogen_axial_warning`
- `isotope_detail_csv`

The isotope-resolved diagnostic summaries are:

```text
cevns_rate_summary_by_isotope.csv
all_models_cevns_event_summary_by_isotope.csv
```

Each CEvNS gas subdirectory contains primary gas-total spectra:

```text
cevns_recoil_spectrum.csv
cevns_accepted_Enu_spectrum.csv
cevns_reco_energy_min_spectrum.csv
```

and isotope-resolved diagnostic spectra:

```text
cevns_recoil_spectrum_by_isotope.csv
cevns_accepted_Enu_spectrum_by_isotope.csv
cevns_reco_energy_min_spectrum_by_isotope.csv
```

Solar CEvNS recoil spectrum columns:

- `gas_label`
- `gas`
- `T_N_keV`
- `dR_dT_s-1_keV-1_vector`
- `dR_dT_s-1_keV-1_axial`
- `dR_dT_s-1_keV-1_total`

Supernova CEvNS recoil spectrum columns use counts:

- `dN_dT_keV-1_vector`
- `dN_dT_keV-1_axial`
- `dN_dT_keV-1_total`

CEvNS lower-bound spectra use:

- `E_nu_reco_min_MeV`
- component totals per bin;
- component densities per MeV;
- `reconstruction_note`.

## Practical Commands

### Environment

Use `python3` in WSL:

```bash
python3 --version
python3 -m pip install -r requirements.txt
```

Compile-check the code:

```bash
python3 -m py_compile detector_model.py cevns.py *.py DetectorNumbers/*.py
```

### Detector-Number Preprocessing

Run these when gas densities, stopping-power inputs, or Garfield transport
summaries change:

```bash
python3 DetectorNumbers/computePlotRanges.py
python3 DetectorNumbers/computeDiffusion.py
```

### Solar, Electron Scattering Only

```bash
python3 plotFluxes.py
python3 scatteringPlots.py
python3 gasTargetRates.py --disable-cevns
python3 recoNuEnergyComparison.py --disable-cevns
```

### Solar With CEvNS

```bash
python3 plotFluxes.py
python3 gasTargetRates.py --enable-cevns --cevns-threshold-kev 1.0
python3 recoNuEnergyComparison.py --enable-cevns --cevns-threshold-kev 1.0
```

### Supernova, Electron Scattering Only

```bash
python3 generate_sn_fluence.py
python3 supernovaGasTargetEvents.py --disable-cevns
python3 supernovaRecoNuEnergyComparison.py --disable-cevns
```

### Supernova With CEvNS

```bash
python3 generate_sn_fluence.py
python3 supernovaGasTargetEvents.py --enable-cevns --cevns-threshold-kev 1.0
python3 supernovaRecoNuEnergyComparison.py --enable-cevns --cevns-threshold-kev 1.0
```

### Supernova Distance Scan

Run every gas row in the gas-density table:

```bash
python3 supernovaDistanceGasScan.py --nr-threshold-kev 1.0
```

To inspect or restrict the gas rows:

```bash
python3 supernovaDistanceGasScan.py --list-gases
python3 supernovaDistanceGasScan.py --gas cf4_1atm --nr-threshold-kev 1.0
```

The distance scan writes one subfolder per gas row under
`supernova_distance_gas_scan/<model>/<gas_slug>/`. It reports expected counts
and average burst rates for both neutrino-electron scattering and CEvNS. The
average rates are `counts / --burst-duration-s`; the input supernova fluence is
time-integrated.

### CEvNS Plot/Save Switches

Compute CEvNS but do not write CEvNS plots:

```bash
python3 gasTargetRates.py --enable-cevns --skip-cevns-plots
```

Compute CEvNS but do not write CEvNS CSV tables:

```bash
python3 gasTargetRates.py --enable-cevns --skip-cevns-save
```

### Standalone CEvNS Diagnostic

The standalone helper is useful for checking isotope target counts and sample
cross sections:

```bash
python3 cevns.py --energy-mev 10 --recoil-kev 0.1 --form-factor default
```

## Smoke Tests

The following smoke sequence exercises the main workflows without writing into
the repository output directories:

```bash
python3 -m py_compile detector_model.py cevns.py *.py DetectorNumbers/*.py
python3 plotFluxes.py
python3 scatteringPlots.py --outdir /tmp/cygnus_review_scattering

python3 gasTargetRates.py \
  --disable-cevns \
  --t-grid-points 40 \
  --costh-grid-points 40 \
  --outdir /tmp/cygnus_review_solar_no_cevns

python3 gasTargetRates.py \
  --enable-cevns \
  --cevns-threshold-kev 1.0 \
  --t-grid-points 80 \
  --costh-grid-points 80 \
  --outdir /tmp/cygnus_review_solar_with_cevns

python3 recoNuEnergyComparison.py \
  --enable-cevns \
  --cevns-threshold-kev 1.0 \
  --t-grid-points 40 \
  --max-true-energy-points 80 \
  --outdir /tmp/cygnus_review_solar_reco_with_cevns

python3 generate_sn_fluence.py --output-root /tmp/cygnus_review_sn_fluence

python3 supernovaGasTargetEvents.py \
  --disable-cevns \
  --fluence-root /tmp/cygnus_review_sn_fluence \
  --output-root /tmp/cygnus_review_sn_events_no_cevns \
  --t-grid-points 40 \
  --costh-grid-points 40

python3 supernovaGasTargetEvents.py \
  --enable-cevns \
  --fluence-root /tmp/cygnus_review_sn_fluence \
  --output-root /tmp/cygnus_review_sn_events_with_cevns \
  --t-grid-points 80 \
  --costh-grid-points 80

python3 supernovaRecoNuEnergyComparison.py \
  --enable-cevns \
  --fluence-root /tmp/cygnus_review_sn_fluence \
  --output-root /tmp/cygnus_review_sn_reco_with_cevns \
  --t-grid-points 40 \
  --max-true-energy-points 80
```

Additional switch checks:

```bash
python3 gasTargetRates.py \
  --enable-cevns \
  --skip-cevns-plots \
  --t-grid-points 20 \
  --costh-grid-points 20 \
  --outdir /tmp/cygnus_review_skip_plots

python3 gasTargetRates.py \
  --enable-cevns \
  --skip-cevns-save \
  --t-grid-points 20 \
  --costh-grid-points 20 \
  --outdir /tmp/cygnus_review_skip_save
```

Validation record from April 29, 2026:

- Compile check passed.
- Solar flux generation passed.
- One-electron scattering plots passed.
- Solar gas rates passed with CEvNS disabled.
- Solar gas rates passed with CEvNS enabled.
- Solar CEvNS lower-bound spectra passed.
- Supernova fluence generation passed.
- Supernova gas counts passed with CEvNS disabled.
- Supernova gas counts passed with CEvNS enabled.
- Supernova CEvNS lower-bound spectra passed.
- `--skip-cevns-plots` wrote CEvNS CSV files and no CEvNS PNG summary.
- `--skip-cevns-save` wrote CEvNS PNG files and no CEvNS CSV summary.
- Primary CEvNS summaries are gas-aggregated: 7 gas rows for the current gas
  table.
- Isotope-resolved CEvNS summaries are still available: 21 isotope rows for the
  current gas table.
- CF4 aggregation check passed: gas-total `CF4 @ 1 atm` equals `12C + 19F`
  from the isotope diagnostic file within numerical precision.

Representative CEvNS values for the default 200 m3 detector, threshold
`1.0 keV`, `form_factor=default`, and
`axial_model=hoferichter_19f_fast`:

| Quantity | Value |
| --- | ---: |
| Solar CEvNS total, all configured gases | `2.6119615624430152e-05 s^-1` |
| Solar CEvNS total, `CF4 @ 1 atm` | `1.7377926123377398e-06 s^-1` |
| SN1987A-like CEvNS total, all configured gases | `30.19570376170965 counts` |
| SN1987A-like CEvNS total, `CF4 @ 1 atm` | `2.0093971116078406 counts` |

The final smoke script printed:

```text
aggregation_checks_ok
FULL_SMOKE_OK
```

## Limitations

- The solar oscillation/flavor split is compact and approximate.
- The neutrino-electron gas-target rates include range/diffusion acceptance,
  but not a full detector simulation.
- CEvNS has no directionality.
- CEvNS uses a true nuclear-recoil threshold. Quenching and
  electron-equivalent threshold conversion are not implemented.
- Hydrogen axial/proton neutral-current elastic scattering is not implemented;
  hydrogen contributes only through the suppressed vector term.
- The default `1H` and `4He` vector form factors are pointlike. Helm for those
  isotopes is only a comparison placeholder.
- Isotope masses use `A * atomic_mass_unit`; atomic binding and electron-mass
  corrections are below the intended precision of this sandbox.
- The `19F` HMS axial response is used as a central first-pass model. The HMS
  uncertainty band is not propagated through the rates.
- CEvNS output aggregation changes only the reporting layer: isotope-level
  calculations are still performed and written as diagnostic files.
- Generated output folders can be stale. Rerun the workflows after changing
  inputs, detector settings, or physics helpers.
