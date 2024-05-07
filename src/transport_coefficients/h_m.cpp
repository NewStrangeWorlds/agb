
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


bool GasHm::calcContinuumAbsorption(const double temperature, const std::vector<double>& number_densities, std::vector<double>& absorption_coeff)
{
  const size_t nb_wavelengths = spectral_grid->nbSpectralPoints();

  std::vector<double> kappa_bf =  boundFreeAbsorption(temperature);
  std::vector<double> kappa_ff =  freeFreeAbsorption(temperature);

  
  const double electron_pressure = number_densities[_e_m] * constants::boltzmann_k * temperature; //electron pressure in dyne/cm2
  const double hydrogen_density = number_densities[_H];  //hydrogen number density in cm-3

  //and now finally the absorption_coeff (in cm-1)
  #pragma omp parallel for
  for (size_t i=0; i<nb_wavelengths; ++i)
      absorption_coeff[i] = (kappa_bf[i] + kappa_ff[i]) * electron_pressure * hydrogen_density;

  return true;
}


std::vector<double> GasHm::boundFreeAbsorption(const double temperature)
{
  const double alpha = constants::light_c*constants::planck_h/constants::boltzmann_k*10000.0; // 1.439e4;
  const double lambda_0 = 1.6419; //photo-detachment threshold
  const double C_n[7] = {0.0, 152.519, 49.534, -118.858, 92.536, -34.194, 4.982};


  auto f = [&] (const double &lambda)
    {
      double x = 0;

      for (unsigned int i=1; i<7; ++i)
        x += C_n[i] * std::pow((1.0/lambda - 1.0/lambda_0), (i-1)/2.0);

      return x;
    };


  const size_t nb_wavelengths = spectral_grid->nbSpectralPoints();
  std::vector<double> kappa_bf(nb_wavelengths, 0.0);
 
  //first, we calculate the photo-detachment cross-section (in cm2)
  #pragma omp parallel for
  for (size_t i=0; i<nb_wavelengths; ++i)
    kappa_bf[i] = 1e-18 * spectral_grid->wavelength_list[i] * spectral_grid->wavelength_list[i] * spectral_grid->wavelength_list[i]
                  * std::pow(1.0/spectral_grid->wavelength_list[i] - 1.0/lambda_0, 1.5) * f(spectral_grid->wavelength_list[i]);

  //and now the kappa (in cm4/dyne)
  #pragma omp parallel for
  for (size_t i=0; i<nb_wavelengths; ++i)
  {
    if (spectral_grid->wavelength_list[i] <= lambda_0 && spectral_grid->wavelength_list[i] >=0.125)
      kappa_bf[i] = 0.750 * std::pow(temperature, -2.5) * std::exp(alpha/lambda_0 / temperature)
                    * (1.0 - std::exp( -alpha/spectral_grid->wavelength_list[i] / temperature)) * kappa_bf[i];
      //kappa_bf[i] = 0.750 * std::pow(temperature, -2.5) * std::exp(1.73830)
      //              * (1.0 - std::exp( -alpha/spectral_grid->wavelength_list[i] / temperature)) * kappa_bf[i];  
    else
      kappa_bf[i] = 0.0;
  }
    

  return kappa_bf;
}



std::vector<double> GasHm::freeFreeAbsorption(const double temperature)
{
  //for wavelengths larger than 0.3645 micron
  const double A_n1[7] = {0.0, 0.0, 2483.3460, -3449.8890, 2200.0400, -696.2710, 88.2830};
  const double B_n1[7] = {0.0, 0.0, 285.8270, -1158.3820, 2427.7190, -1841.4000, 444.5170};
  const double C_n1[7] = {0.0, 0.0, -2054.2910, 8746.5230, -13651.1050, 8624.9700, -1863.8650};
  const double D_n1[7] = {0.0, 0.0, 2827.7760, -11485.6320, 16755.5240, -10051.5300, 2095.2880};
  const double E_n1[7] = {0.0, 0.0, -1341.5370, 5303.6090, -7510.4940, 4400.0670, -901.7880};
  const double F_n1[7] = {0.0, 0.0, 208.9520, -812.9390, 1132.7380, -655.0200, 132.9850};

  //for wavelengths between 0.1823 micron and 0.3645 micron
  const double A_n2[7] = {0.0, 518.1021, 473.2636, -482.2089, 115.5291, 0.0, 0.0};
  const double B_n2[7] = {0.0, -734.8666, 1443.4137, -737.1616, 169.6374, 0.0, 0.0};
  const double C_n2[7] = {0.0, 1021.1775, -1977.3395, 1096.8827, -245.6490, 0.0, 0.0};
  const double D_n2[7] = {0.0, -479.0721, 922.3575, -521.1341, 114.2430, 0.0, 0.0};
  const double E_n2[7] = {0.0, 93.1373, -178.9275, 101.7963, -21.9972, 0.0, 0.0};
  const double F_n2[7] = {0.0, -6.4285, 12.3600, -7.0571, 1.5097, 0.0, 0.0};


  auto ff = [&] (const double &lambda, 
                 const double A_n[7], const double B_n[7], const double C_n[7], const double D_n[7], const double E_n[7], const double F_n[7])
    {
      double x = 0;

      for (unsigned int i=1; i<7; ++i)
        x += std::pow(5040.0/temperature, (i+1)/2.0) * 
             (lambda * lambda * A_n[i] + B_n[i] + C_n[i]/lambda + D_n[i]/lambda/lambda 
               + E_n[i]/lambda/lambda/lambda + F_n[i]/lambda/lambda/lambda/lambda);

      return x*1e-29;
    };


  const size_t nb_wavelengths = spectral_grid->nbSpectralPoints();
  std::vector<double> kappa_ff(nb_wavelengths, 0.0);


  //compute free-free kappa (in cm4/dyne)
  #pragma omp parallel for
  for (size_t i=0; i<nb_wavelengths; ++i)
  {
    if (spectral_grid->wavelength_list[i] > 0.3645)
      kappa_ff[i] = ff(spectral_grid->wavelength_list[i], A_n1, B_n1, C_n1, D_n1, E_n1, F_n1);
    else if (spectral_grid->wavelength_list[i] >= 0.1823)
      kappa_ff[i] = ff(spectral_grid->wavelength_list[i], A_n2, B_n2, C_n2, D_n2, E_n2, F_n2);
  }
  


  return kappa_ff;
}

}