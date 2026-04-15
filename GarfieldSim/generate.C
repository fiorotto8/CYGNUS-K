#include <iostream>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/FundamentalConstants.hh"

using namespace Garfield;
using namespace std;

int main(int argc, char * argv[]) {

  // Define the pressure in Torr (500 mbar converted to Torr).
  const double pressure = 500 * 0.750062;
  // Define the temperature in Kelvin.
  const double temperature = 20 + 273.15;

  // Output the atmospheric pressure in some units (conversion factor applied).
  cout <<"Atmopheric pressure (mbar): " << AtmosphericPressure * 1.33322 << endl;
  cout <<"Atmopheric pressure (torr): " << AtmosphericPressure<< endl;
  cout << "Detector Pressure: " << pressure << " Torr" << endl;

  // Setup the gas mixture with 60% Helium and 40% CF4.
  MediumMagboltz gas("He", 58., "CF4", 37., "CH4", 5.);
  //MediumMagboltz gas("CF4", 100.);
  // Set the temperature of the gas.
  gas.SetTemperature(temperature);
  // Set the pressure of the gas.
  gas.SetPressure(pressure);

  // Define the field range to be covered by the gas table.
  const size_t nE = 20; // Number of field points.
  const double emin = 100.; // Minimum electric field in V/cm.
  const double emax = 10000.; // Maximum electric field in V/cm.
  // Flag to request logarithmic spacing of the field points.
  constexpr bool useLog = true;
  gas.SetFieldGrid(emin, emax, nE, useLog); 

  // Define the number of collisions for Magboltz simulation.
  const int ncoll = 20;
  // Run Magboltz to generate the gas table.
  gas.GenerateGasTable(ncoll);
  // Save the generated gas table to a file.
  gas.WriteGasFile("hecf4ch4_58-37-5_500mbar20C.gas");

}
