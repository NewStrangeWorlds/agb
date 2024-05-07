

#include "species_definition.h"
#include "opacity_species.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <omp.h>


#include "../config/config.h"
#include "../spectral_grid/spectral_grid.h"
#include "../chemistry/chem_species.h"
#include "../additional/physical_const.h"


namespace agb{


bool GasH2m::calcContinuumAbsorption(const double temperature, const std::vector<double>& number_densities, std::vector<double>& absorption_coeff)
{
  const size_t nb_wavelengths = spectral_grid->nbSpectralPoints();

  std::vector<double> kappa_ff =  freeFreeAbsorption(temperature);

  
  const double electron_pressure = number_densities[_e_m] * constants::boltzmann_k * temperature; //electron pressure in dyne/cm2
  const double h2_density = number_densities[_H2];  //H2 number density in cm-3

  //and now finally the absorption_coeff (in cm-1)
  #pragma omp parallel for
  for (size_t i=0; i<nb_wavelengths; ++i)
      absorption_coeff[i] = kappa_ff[i] * electron_pressure * h2_density;

  return true;
}



std::vector<double> GasH2m::freeFreeAbsorption(const double temperature)
{
  double theta = 5040.4/temperature;
  
  if (theta > max_theta) theta = max_theta;
  if (theta < min_theta) theta = min_theta;

  
  const size_t nb_wavelengths = spectral_grid->nbSpectralPoints();
  std::vector<double> kappa_ff(nb_wavelengths, 0.0);

  //fit function for the absorption coefficients
  //x is log10(lambda), y is theta
  //the result is log10(kappa_ff)
  auto kappaFit = [&](const double x, const double y)
                    {return fit_coeff[0] + 
                            fit_coeff[1]*x + fit_coeff[2]*y + 
                            fit_coeff[3]*x*x + fit_coeff[4]*x*y + fit_coeff[5]*y*y + 
                            fit_coeff[6]*std::pow(x,3) + fit_coeff[7]*x*x*y + fit_coeff[8]*x*y*y + fit_coeff[9]*std::pow(y,3) + 
                            fit_coeff[10]*std::pow(x,4) + fit_coeff[11]*std::pow(x,3)*y + fit_coeff[12]*x*x*y*y + fit_coeff[13]*x*std::pow(y,3) + fit_coeff[14]*std::pow(y,4) +
                            fit_coeff[15]*std::pow(x,5) + fit_coeff[16]*std::pow(x,4)*y + fit_coeff[17]*std::pow(x,3)*y*y + fit_coeff[18]*x*x*std::pow(y,3) + fit_coeff[19]*x*std::pow(y,4) + fit_coeff[20]*std::pow(y,5);};
  
  
  //compute free-free kappa (in cm4/dyne)
  #pragma omp parallel for
  for (size_t i=0; i<nb_wavelengths; ++i)
  {
    if (spectral_grid->wavelength_list[i] >= min_lambda && spectral_grid->wavelength_list[i] <= max_lambda)
      kappa_ff[i] = 1e-26 * std::pow(10, kappaFit(std::log10(spectral_grid->wavelength_list[i]), theta));
  }
  


  return kappa_ff;
}

}