# Tools for Directional Neutrinos

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
