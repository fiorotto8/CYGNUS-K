# Tools for Directional Neutrinos

## EnergywithPerformance.py

This script comapres directional and non directional detectors and makes 2 simple plots for elastic scattering: nu_e + e -> nu_e + e

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

- Borexino-like:
  - threshold = 80 keV
  - energy resolution = 5% / sqrt(E[MeV])

- Directional:
  - threshold = 10 keV
  - energy resolution = 20%
  - angular resolution = 20 degrees

### How to change detector performance

Inside `main()` you can modify:

```python
borexino_like = DetectorPerformance(...)
directional = DetectorPerformance(...)