#include <iostream>
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/FundamentalConstants.hh"

using namespace Garfield;
using namespace std;

int main(int argc, char* argv[]) {

  // ===== User‐adjustable parameters =====
  const double baseHe      = 60.0;    // % He before impurities
  const double baseCF4     = 40.0;    // % CF4 before impurities
  const double impurityO2  = 0.04;    // 400 ppm → 0.04%
  //const double impurityH2O = 0.0442;  // 1.7% RH 
  const double impurityH2O = 2.6;  
  // ======================================

  // Rescale He/CF4 so that He+CF4+O2+H2O = 100%
  const double totalImp = impurityO2 + impurityH2O;
  if (totalImp >= 100.0) {
    cerr << "Error: Total impurities ≥ 100%\n";
    return 1;
  }
  const double scale   = (100.0 - totalImp) / (baseHe + baseCF4);
  const double fracHe  = baseHe  * scale;
  const double fracCF4 = baseCF4 * scale;

  // Operating conditions
  const double pressure    = 900.0 * 0.750062; // 900 mbar → Torr
  const double temperature = 20.0  + 273.15;   // 20 °C → K

  cout << "Detector P = " << pressure    << " Torr\n"
       << "Detector T = " << temperature << " K\n\n";

  // === Build the gas mixture in one call ===
  // Note: the constructor supports up to 6 pairs (gas, fraction)
  MediumMagboltz gas(
    "He",  fracHe,
    "CF4", fracCF4,
    "O2",  impurityO2,
    "H2O", impurityH2O
    // the remaining two slots default to "" and 0.
  );

  gas.SetPressure(pressure);
  gas.SetTemperature(temperature);

  // Field grid (logarithmic in E)
  const size_t nE      = 20;
  const double emin    = 100.0;    // V/cm
  const double emax    = 10000.0;  // V/cm
  constexpr bool useLog = true;
  gas.SetFieldGrid(emin, emax, nE, useLog);

  // Run Magboltz
  const int ncoll = 20;
  gas.GenerateGasTable(ncoll);

  // Write out the .gas file
  gas.WriteGasFile("hecf4_o2_h2o_impure_2.6percent_h20.gas");

  // Summary
  cout << "\n=== Final Mixture ===\n"
       << "  He   = " << fracHe       << " %\n"
       << "  CF4  = " << fracCF4      << " %\n"
       << "  O2   = " << impurityO2   << " %\n"
       << "  H2O  = " << impurityH2O  << " %\n";

  return 0;
}
