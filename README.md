# Tools for Directional Neutrinos

## Short recap of the codes

This repository contains small analysis tools for solar-neutrino fluxes, oscillations, and neutrino-electron scattering, plus Garfield-based gas simulation utilities.

- `downladFluxes.py`: downloads solar-neutrino spectrum tables and writes a local manifest.
- `plotFluxes.py`: builds source and flavor-separated solar-neutrino fluxes at Earth and exports plots/CSV.
- `scatteringPlots.py`: computes and visualizes neutrino-electron scattering observables from the flux inputs.
- `EnergywithPerformance.py`: compares directional vs non-directional detector response assumptions.
- `GarfieldSim/`: Garfield/C++ and Python helpers for gas-mixture generation, transport studies, and table/plot production. 
- `DetectorNumbers/`: tabulated gas densities and stopping powers for relevant mixtures, plus utilities for computing diffusion and range.

Detailed descriptions for each code are reported below.

## `EnergywithPerformance.py`

This script compares directional and non-directional detectors for elastic scattering:

$$
\nu_e + e^- \rightarrow \nu_e + e^-
$$

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

This script downloads tabulated solar-neutrino source spectra from the Bahcall SNdata archive [Ref](https://www.sns.ias.edu/%E2%88%BCjnb) and saves them into:

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

## Solar neutrino flux tools (plotFluxes.py)

### What the code does

This script:

- reads tabulated solar-neutrino source spectra from local files
- converts each tabulated spectrum into a physical differential flux at Earth
- extends continuum spectra down to 1 keV with a simple matched low-energy extrapolation
- adds line-like sources (`pep`, `^7$Be`)
- builds the total solar-neutrino source flux
- applies a minimal adiabatic MSW-LMA oscillation model
- produces flavor-separated fluxes at Earth
- saves plots and a CSV table

The goal is a practical, detector-oriented flux model suitable for rate estimates and first-order phenomenology, not a precision solar-neutrino oscillation calculation.

---

### Input spectra and normalization

The continuum source spectra are read from tabulated files for:

- `pp`
- `8B`
- `hep`
- `13N`
- `15O`
- `17F`

These tables provide source spectral **shapes** as functions of neutrino energy. The code does **not** assume they are already in physical units of flux. Instead, each spectrum is rescaled to the desired total source flux at Earth.

For each continuum component, the code constructs

$$
\frac{d\Phi}{dE}(E)
=
\Phi_{\rm tot}
\frac{s(E)}{\int s(E')\,dE'}
$$

where:

- $s(E)$ is the tabulated spectral shape
- $\Phi_{\rm tot}$ is the adopted total flux at Earth [Ref](https://iopscience.iop.org/article/10.1088/1742-6596/1056/1/012058/pdf)
- $d\Phi/dE$ is the physical differential flux in

$$
\mathrm{cm^{-2}\,s^{-1}\,MeV^{-1}}.
$$

Numerically, the integral is approximated from the tabulated points by constructing bin edges from midpoints between adjacent energy samples.

---

### Treatment of line-like sources

Two neutrino sources are effectively monochromatic:

- `pep` at 1.442 MeV
- `^7$Be` at 0.3843 MeV and 0.8613 MeV

Since a true delta-function line cannot be displayed directly on a finite plotting grid, the code represents each line as a narrow top-hat bin of width

$$
\Delta E_{\rm line} = 0.002\ {\rm MeV}.
$$

The plotted line height is therefore

$$
\left(\frac{d\Phi}{dE}\right)_{\rm line}
=
\frac{\Phi_{\rm line}}{\Delta E_{\rm line}},
$$

so that the integrated flux is preserved:

$$
\int_{\rm line} \frac{d\Phi}{dE}\, dE = \Phi_{\rm line}.
$$

For `^7$Be`, the code can either:

- use a tabulated thermal lineshape for the 0.8613 MeV branch, plus a narrow bin for the 0.3843 MeV branch, or
- represent both branches as finite-width bins.

The adopted branching fractions are [Reference](http://www.lnhb.fr/nuclides/Be-7_com.pdf):

- 0.1044 at 0.3843 MeV
- 0.8956 at 0.8613 MeV

---

### Low-energy extrapolation

Some continuum tables do not start at very low energy. For plotting and rate studies, the code extends each continuum spectrum down to

$$
E_{\min} = 1\ {\rm keV} = 10^{-3}\ {\rm MeV}
$$

using a matched low-energy law

$$
\frac{d\Phi}{dE} \propto E^2.
$$

More explicitly, if the first valid tabulated point is $(E_0, y_0)$, the low-energy continuation is taken as

$$
\frac{d\Phi}{dE}(E) = A E^2,
\qquad
A = \frac{y_0}{E_0^2},
\qquad E < E_0.
$$

After this extension, the spectrum is renormalized so that the total integrated flux remains equal to the target value $\Phi_{\rm tot}$.

This is only a practical extrapolation. It is not intended as a precision model of the extreme low-energy endpoint behavior.

---

### Oscillations: minimal adiabatic MSW-LMA model

The code applies a simplified MSW-LMA treatment for the solar electron-neutrino survival probability.

The 3-flavor survival probability is approximated as [Ref](https://arxiv.org/pdf/1001.4524.pdf):

$$
P_{ee}^{3\nu}(E)
\simeq
\sin^4\theta_{13}
+
\cos^4\theta_{13}\,
P_{ee}^{2\nu}(E),
$$

which is the standard reduction used in solar-neutrino phenomenology when the third mixing angle is included perturbatively.

The effective 2-flavor survival probability in matter is evaluated in the adiabatic approximation:

$$
P_{ee}^{2\nu}(E)
=
\frac{1}{2}
\left[
1 + \cos 2\theta_{12}\,\cos 2\theta_{12}^m(E)
\right].
$$

Here $\theta_{12}^m$ is the matter mixing angle at production, defined through

$$
\cos 2\theta_{12}^m
=
\frac{\cos 2\theta_{12} - \beta}
{\sqrt{(\cos 2\theta_{12}-\beta)^2 + \sin^2 2\theta_{12}}},
$$

with

$$
\beta = \frac{A}{\Delta m_{21}^2}.
$$

The matter potential is written in the code as

$$
A\,[{\rm eV}^2]
=
1.52\times 10^{-7}\,
\left(n_e \,[{\rm mol/cm^3}]\right)\,
\left(E_\nu\,[{\rm MeV}]\right).
$$

This is the standard charged-current matter term written in convenient solar-neutrino units.

Because the calculation is meant to stay simple, each neutrino source is assigned a single effective electron density $n_e^{\rm eff}$ representing its production region in the Sun. The code therefore does **not** solve the full radial propagation problem and does **not** average over a detailed source-production profile.

The flavor-separated fluxes at Earth are then built as

$$
\phi_{\nu_e}(E) = P_{ee}(E)\,\phi_{\rm source}(E),
$$

$$
\phi_{\nu_\mu}(E) + \phi_{\nu_\tau}(E)
=
\left[1-P_{ee}(E)\right]\phi_{\rm source}(E).
$$

In the default implementation, the non-electron component is split equally:

$$
\phi_{\nu_\mu}(E)=\phi_{\nu_\tau}(E)
=
\frac{1-P_{ee}(E)}{2}\,\phi_{\rm source}(E).
$$

This corresponds to the choice `equal_mu_tau=True`.

---

### Oscillation parameters used

The script uses fixed benchmark oscillation parameters:

- $\sin^2\theta_{12} = 0.307$
- $\sin^2\theta_{13} = 0.02215$
- $\Delta m^2_{21} = 7.49\times10^{-5}\ {\rm eV}^2$

These are intended as representative central values for a minimal implementation.

---

### Effective production-region electron densities

The script uses the following source-dependent effective electron densities:

- `pp`: 58 mol/cm$^3$
- `pep`: 61 mol/cm$^3$
- `7Be`: 66 mol/cm$^3$
- `8B`: 93 mol/cm$^3$
- `hep`: 93 mol/cm$^3$
- `13N`: 81 mol/cm$^3$
- `15O`: 86 mol/cm$^3$
- `17F`: 86 mol/cm$^3$

These numbers should be described carefully.

They are **not** taken from a single published table and they were **not** obtained here from a full solar-model propagation calculation. They were infered by a LLM and introduced as a practical source-dependent approximation to capture the fact that different neutrino components are produced at different depths in the solar core.
The resulted probabilities are consistent with the Figure 14.3 of the PDG 2024 review. See `solar_neutrino_pee_logx.png`  

---

### Output files

The script produces:

- `solar_neutrino_fluxes_loglog_physical.png`  
  individual source flux components before oscillation

- `solar_neutrino_total_flux_loglog.png`  
  total source flux before oscillation

- `solar_neutrino_flavor_fluxes_loglog.png`  
  flavor-separated fluxes at Earth after the minimal MSW-LMA treatment

- `solar_neutrino_pee_logx.png`  
  electron-neutrino survival probability $P_{ee}(E)$

- `solar_neutrino_total_flux.csv`  
  table containing:
  - total source flux
  - $\nu_e$ flux
  - $\nu_\mu$ flux
  - $\nu_\tau$ flux

All differential fluxes are saved in units of

$$
\mathrm{cm^{-2}\,s^{-1}\,MeV^{-1}}.
$$

---

### Limitations

This is a useful first-level implementation, but it remains simplified.

Main limitations:

- the low-energy continuation of continuum spectra is an extrapolation
- monochromatic lines are represented with finite plotting width
- the oscillation treatment assumes adiabatic MSW-LMA propagation
- each source is assigned a single effective production density
- the code does not integrate over the full solar density profile
- Earth matter effects and day/night asymmetry are neglected
- no seasonal variation is included
- no dedicated uncertainty propagation is performed

Therefore, the code is appropriate for detector studies, plotting, and first-order rate estimates, but not for precision solar-neutrino oscillation analyses.

# plotFluxes.py

## Solar neutrino $\nu$-e scattering script

This script takes as input the **flavor-separated solar-neutrino fluxes at Earth** produced by the previous flux/oscillation code, and computes the expected elastic-scattering signal on a **single electron target**.

The input CSV is assumed to contain:

- $E_\nu$ in MeV
- total differential flux
- $\nu_e$ differential flux
- $\nu_\mu$ differential flux
- $\nu_\tau$ differential flux

with flux units

$$
\frac{d\Phi}{dE_\nu}
\quad [\mathrm{cm^{-2}\ s^{-1}\ MeV^{-1}}].
$$

The script then combines these fluxes with approximate $\nu e^-$ elastic-scattering cross sections and with the Standard Model angular shape of the recoil electron distribution.

---

## Physics implemented

### 1. Total $\nu$-e cross sections

The script uses simple linear approximations for the total cross section:

$$
\sigma(E_\nu) = k\,E_\nu,
$$

with $E_\nu$ in MeV and $\sigma$ in cm$^2$.

The coefficients used are [Ref](https://arxiv.org/pdf/hep-ph/0307222) and [Ref](https://arxiv.org/pdf/hep-ph/0408031):

$$
\sigma(\nu_e e^-)\simeq 9.20\times 10^{-45}\,E_\nu(\mathrm{MeV})\ \mathrm{cm^2},
$$

$$
\sigma(\bar{\nu}_e e^-)\simeq 3.83\times 10^{-45}\,E_\nu(\mathrm{MeV})\ \mathrm{cm^2},
$$

$$
\sigma(\nu_{\mu,\tau} e^-)\simeq 1.57\times 10^{-45}\,E_\nu(\mathrm{MeV})\ \mathrm{cm^2},
$$

$$
\sigma(\bar{\nu}_{\mu,\tau} e^-)\simeq 1.29\times 10^{-45}\,E_\nu(\mathrm{MeV})\ \mathrm{cm^2}.
$$

For solar neutrinos, only neutrino channels are relevant, so in practice:

- $\nu_e$ uses $\sigma_{\nu_e e}$,
- $\nu_\mu$ and $\nu_\tau$ use $\sigma_{\nu_x e}$, with $x=\mu,\tau$.

---

### 2. Differential cross section in recoil energy

To describe the angular shape of the scattered electron, the script starts from the Standard Model tree-level differential cross section with respect to the electron recoil kinetic energy $T$:

$$
\frac{d\sigma}{dT}
=
\frac{2G_F^2 m_e}{\pi}
\left[
g_L^2
+
g_R^2\left(1-\frac{T}{E_\nu}\right)^2
-
g_L g_R \frac{m_e T}{E_\nu^2}
\right].
$$

Here:

- $G_F$ is the Fermi constant,
- $m_e$ is the electron mass,
- $g_L$ and $g_R$ are the weak couplings,
- $E_\nu$ is the neutrino energy,
- $T$ is the recoil kinetic energy of the electron.

The couplings depend on the neutrino flavor:

For $\nu_e$ scattering:

$$
g_L = \frac{1}{2} + \sin^2\theta_W,
\qquad
g_R = \sin^2\theta_W,
$$

while for $\nu_{\mu,\tau}$ scattering:

$$
g_L = -\frac{1}{2} + \sin^2\theta_W,
\qquad
g_R = \sin^2\theta_W.
$$

This difference reflects the fact that $\nu_e e^-$ scattering receives both charged-current and neutral-current contributions, while $\nu_\mu e^-$ and $\nu_\tau e^-$ receive only neutral-current contributions.

---

### 3. Recoil-energy kinematics

For elastic scattering of a neutrino on an electron initially at rest, the maximum allowed recoil energy is

$$
T_{\max}(E_\nu)
=
\frac{2E_\nu^2}{m_e + 2E_\nu}.
$$

This defines the physical range of the recoil electron kinetic energy.

The script also uses the inverse relation giving the minimum neutrino energy required to produce a recoil $T$:

$$
E_\nu^{\min}(T)
=
\frac{1}{2}
\left(
T + \sqrt{T^2 + 2m_e T}
\right).
$$

---

### 4. Relation between recoil angle and recoil energy

To obtain the angular distribution of the outgoing electron, the script uses the relation between the recoil angle $\theta_e$ and the recoil kinetic energy $T$.

Writing $c=\cos\theta_e$, the implemented formula is

$$
T(E_\nu,c)
=
\frac{2m_e E_\nu^2 c^2}
{(E_\nu+m_e)^2 - E_\nu^2 c^2}.
$$

The Jacobian for the change of variable from $T$ to $\cos\theta_e$ is

$$
\frac{dT}{d\cos\theta_e}
=
\frac{4m_e E_\nu^2 (E_\nu+m_e)^2 \cos\theta_e}
{\left[(E_\nu+m_e)^2 - E_\nu^2 \cos^2\theta_e\right]^2}.
$$

The differential cross section in angle is then obtained as

$$
\frac{d\sigma}{d\cos\theta_e}
=
\frac{d\sigma}{dT}
\frac{dT}{d\cos\theta_e}.
$$

This gives the **shape** of the electron recoil angular distribution.

---

### 5. Rescaling to the chosen total cross section

The script uses the Standard Model formula above only for the **shape** of the angular distribution. Since the total rate is normalized to the simpler approximate cross sections, the angular differential cross section is rescaled so that its integral matches the chosen total cross section:

$$
\left(\frac{d\sigma}{d\cos\theta_e}\right)_{\mathrm{rescaled}}
=
\left(\frac{d\sigma}{d\cos\theta_e}\right)_{\mathrm{SM}}
\,
\frac{\sigma_{\mathrm{approx}}(E_\nu)}
{\int d(\cos\theta_e)\,
\left(\frac{d\sigma}{d\cos\theta_e}\right)_{\mathrm{SM}} }.
$$

In this way:

- the **shape** comes from the Standard Model tree-level calculation,
- the **normalization** comes from the chosen approximate total cross section.

---

## Flux convolution

### 6. Rate per neutrino-energy bin

For a single electron target, the expected contribution from one neutrino-energy bin is

$$
R_i
=
\phi(E_i)\,\sigma(E_i)\,\Delta E_i,
$$

where:

- $\phi(E_i)$ is the differential neutrino flux in that bin,
- $\sigma(E_i)$ is the total scattering cross section,
- $\Delta E_i$ is the bin width.

The result has units:

$$
[\mathrm{s^{-1}\ per\ electron\ per\ bin}].
$$

This is computed separately for each flavor:

$$
R_i^{(\nu_e)} = \phi_{\nu_e}(E_i)\,\sigma_{\nu_e e}(E_i)\,\Delta E_i,
$$

$$
R_i^{(\nu_\mu)} = \phi_{\nu_\mu}(E_i)\,\sigma_{\nu_x e}(E_i)\,\Delta E_i,
$$

$$
R_i^{(\nu_\tau)} = \phi_{\nu_\tau}(E_i)\,\sigma_{\nu_x e}(E_i)\,\Delta E_i.
$$

The total rate per bin is

$$
R_i^{\mathrm{tot}}
=
R_i^{(\nu_e)}
+
R_i^{(\nu_\mu)}
+
R_i^{(\nu_\tau)}.
$$

---

### 7. Flux-convolved angular rate

The script also computes the angular recoil rate by integrating over the full neutrino spectrum:

$$
\frac{dR}{d\cos\theta_e}
=
\int dE_\nu\,
\phi(E_\nu)\,
\frac{d\sigma}{d\cos\theta_e}(E_\nu).
$$

Numerically, this is approximated as a sum over the input energy bins:

$$
\frac{dR}{d\cos\theta_e}
\simeq
\sum_i
\phi(E_i)\,
\frac{d\sigma}{d\cos\theta_e}(E_i)\,
\Delta E_i.
$$

Again, this is done separately for $\nu_e$, $\nu_\mu$, and $\nu_\tau$, and then summed.

The result has units:

$$
[\mathrm{s^{-1}\ per\ electron}].
$$

---

### 8. Normalized angular probability density

Finally, the script constructs a normalized angular distribution, representing the probability density for the recoil direction given that one interaction has occurred:

$$
P(\cos\theta_e)
=
\frac{\dfrac{dR}{d\cos\theta_e}}
{\int d(\cos\theta_e)\,\dfrac{dR}{d\cos\theta_e}}.
$$

This quantity is dimensionless and integrates to unity:

$$
\int_0^1 P(\cos\theta_e)\,d(\cos\theta_e)=1.
$$

This is useful when one wants to study only the **shape** of the recoil direction distribution, independently of the absolute interaction rate.

---

## What the plots mean

The script produces several plots:

- `cross_sections_plain.png`  
  shows the approximate total $\nu e$ cross sections as functions of neutrino energy.

- `input_fluxes_by_flavor.png`  
  shows the input solar-neutrino fluxes at Earth for $\nu_e$, $\nu_\mu$, and $\nu_\tau$.

- `angular_differential_cross_sections_nu_e.png`  
  shows $d\sigma/d\cos\theta_e$ for $\nu_e e^-$ scattering at selected neutrino energies.

- `angular_differential_cross_sections_nu_x.png`  
  shows $d\sigma/d\cos\theta_e$ for $\nu_{\mu,\tau} e^-$ scattering at selected neutrino energies.

- `convolved_rate_per_energy_bin_by_flavor.png`  
  shows the rate contribution from each neutrino-energy bin for each flavor and for the total.

- `convolved_angular_rate_by_flavor.png`  
  shows the absolute recoil-angle rate distribution for each flavor and for the total.

- `convolved_angular_pdf_by_flavor.png`  
  shows the normalized recoil-angle distributions, useful to compare the directional shapes.

---

## Main approximations

This script is intended for detector studies and first-order estimates. The main simplifications are:

- the total cross sections are taken from linear approximations in $E_\nu$,
- the angular shape is taken from the tree-level Standard Model formula,
- the angular differential cross section is rescaled to match the approximate total cross section,
- detector effects such as threshold, energy resolution, angular resolution, and efficiency are not included,
- the calculation is done for a **single electron target**.

To obtain rates for a real detector, the result must be multiplied by the total number of target electrons:

$$
R_{\mathrm{detector}}
=
N_e \times R_{\mathrm{single\ electron}}.
$$

---

## `GarfieldSim/`

This folder contains the Garfield-based transport and gas-processing tools used to study electron transport in the detector gases. It is the part of the repository that turns a gas mixture definition into concrete transport outputs such as drift velocity, diffusion, Townsend coefficient, and attachment coefficient, and then stores those results in plots or tabulated files for later use.

The main goal of this folder is to provide a consistent workflow around Garfield/Magboltz gas files:

- define or load a gas mixture,
- compute the electron transport properties on the electric-field grid available in the gas table,
- inspect the transport coefficients with plots,
- export summary values into CSV files,
- and keep the resulting plots and tables organized in the local `plots/` and `tables/` subfolders.

### What is inside

The folder includes both C++ ROOT macros and Python scripts. The C++ files are used for gas generation, table merging, Penning handling, and Garfield-side plotting or table production. The Python files are used for lighter post-processing and for reading the gas tables in a more flexible way.

Typical files include:

- `generate.C` and `generate_impurity.C`: scripts used to define or prepare gas mixtures.
- `merge.C`: combines or merges gas-table related outputs.
- `penning.C`: handles Penning transfer assumptions when needed.
- `plotcs.C`, `plotExtrap.C`, `plotvsEN.C`: ROOT macros for cross-section and extrapolation studies.
- `printTable.C` and `printTable.py`: utilities to inspect or print transport tables.
- `read.C` and `read.py`: scripts that load a `.gas` file and produce transport plots and summary values.
- `microscopic.py`: a Python helper used in the more detailed Garfield workflow.
- `tables/`: input gas tables used by the scripts.
- `plots/`: generated figures and transport-summary CSV files.

### What the scripts do

The core transport workflow is centered on a Garfield gas file. A gas file encodes the mixture composition, pressure, temperature, and the tabulated transport data as a function of field. From that table, the scripts can extract and compare:

- electron drift velocity,
- longitudinal and transverse diffusion,
- Townsend coefficient,
- attachment coefficient,
- and related transport observables.

The plotting scripts are used to turn those values into a readable visual summary. The output plots are meant for quick detector studies and comparisons across gas mixtures, not for a full Garfield analysis chain.

### Output philosophy

The emphasis is on practical detector-oriented outputs. In particular, the scripts usually save:

- one plot per transport observable,
- a CSV summary for selected field points,
- and, when needed, comparison plots that place several gas mixtures on the same axes.

This makes it easy to compare gases such as CF4 and He-based mixtures at a glance and to pass the values into later detector calculations.

### Notes on usage

The folder is meant to be used with Garfield/Magboltz data already prepared for the chosen gas mixtures. The scripts do not attempt a full microscopic simulation from scratch every time; instead, they load the prepared gas tables and extract the relevant transport information from them.

In the current workflow, this is also where the transport-point CSV files are created, which are then consumed by the detector-level post-processing scripts in `DetectorNumbers/`.

---

## `DetectorNumbers/`

This folder collects the numerical inputs and post-processing tools used to translate detector gas properties into useful macroscopic quantities. It contains the tabulated gas densities, stopping powers, and the scripts that turn those tables into diffusion and range estimates for the detector studies.

The role of this folder is to bridge the gap between raw gas-property tables and the quantities used in detector design:

- gas density tables provide the material density for a given mixture and pressure,
- stopping-power tables provide the energy-loss input for electron transport calculations,
- and the scripts combine those inputs with Garfield transport outputs to produce range and diffusion summaries.

### Main inputs

The folder contains the gas-property tables used by the analysis, including:

- `GasDensities.csv`: the density table for the gas mixtures considered in the study,
- `CF4_Stopping_power_electrons.csv`, `HeCF4_Stopping_power_electrons.csv`, and `HeCF4CH4_Stopping_power_electrons.csv`: collisional stopping-power tables for the corresponding base gases or mixtures,
- and the derived transport CSVs and plots generated from the Garfield-side workflow.

These files provide the numerical basis for computing how far electrons travel and how much they spread in the gas.

### Main scripts

The two main scripts in this folder are:

- `computePlotRanges.py`: reads the gas-density table and the collisional stopping-power table, then computes electron ranges under the constant slowing-down approximation. The result is plotted as range versus energy for the available gas entries.
- `computeDiffusion.py`: reads the Garfield transport-point CSV files from the plots folder, computes the 2σ diffusion expected over a 50 cm drift length, and writes a summary table plus comparison plots of diffusion versus electric field.

### What the calculations mean

For range calculations, the stopping power is used together with the gas density to estimate the electron range as a function of energy. The constant slowing-down approximation treats the inverse stopping power as the local contribution to the path length. This is a practical detector-level estimate that is good for comparing gases and energy scales.

For diffusion calculations, the Garfield transport CSVs store the diffusion coefficients in units of $\sqrt{\rm cm}$. The script converts those coefficients into a physical spread over a chosen drift length using

$$
\sigma = D\,\sqrt{L},
$$

and then reports the two-sigma width as

$$
2\sigma = 2D\sqrt{L}.
$$

In the current setup, the drift length is 50 cm, so the script reports the expected 2σ spread over that distance for each gas and each field point.

### Outputs

This folder is where the derived detector-number summaries are written. Typical outputs include:

- range-versus-energy plots,
- a CSV summary of diffusion values over 50 cm,
- and comparison figures showing how the diffusion changes with electric field for each gas.

These outputs are intended for quick comparison of mixtures and operating points, especially when deciding which gas is better suited for a directional detector or for a given diffusion target.

### How it fits into the workflow

`DetectorNumbers/` is the numerical post-processing layer after the Garfield transport study. The Garfield folder produces the transport-point tables, and this folder consumes those tables to make higher-level detector metrics. In other words, `GarfieldSim/` answers “what does the gas do at a given field?”, while `DetectorNumbers/` answers “what does that mean for range and drift over a detector-scale distance?”.

---
