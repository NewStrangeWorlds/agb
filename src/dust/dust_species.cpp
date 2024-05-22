
#include "dust_species.h" 

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <omp.h>
#include <sstream>
#include <algorithm> 
#include <assert.h> 
#include <cmath>

#include "../additional/aux_functions.h"
#include "../spectral_grid/spectral_grid.h"
#include "../config/config.h"
#include "../additional/exceptions.h"
#include "../atmosphere/atmosphere.h"
#include "../additional/quadrature.h"
#include "../additional/physical_const.h"
#include "../../_deps/lx_mie-src/mie/mie.h"



namespace agb{


DustSpecies::DustSpecies(
  ModelConfig* config_,
  SpectralGrid* spectral_grid_,
  Atmosphere* atmosphere_)
  : config(config_)
  , spectral_grid(spectral_grid_)
  , atmosphere(atmosphere_)
  , nb_grid_points(atmosphere->nb_grid_points)
  , nb_spectral_points(spectral_grid->nbSpectralPoints())
{
  number_density.assign(nb_grid_points, 0);
  size_distribution.assign(nb_grid_points, std::vector<double>(1, 0.));
  particle_radius.assign(nb_grid_points, std::vector<double>(1, 0.));

  readRefractiveIndexFile(config->model_folder + config->refractive_index_file);
}



void DustSpecies::readRefractiveIndexFile(const std::string file_path)
{
  std::fstream file;
  file.open(file_path.c_str(), std::ios::in);

  if (file.fail()) 
  {
    if (file.fail())
    throw FileNotFound(std::string ("Refractive index file "), file_path);
  }

  std::string line;
  std::string input;

  //header
  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  
  std::vector<double> wavelength;
  std::vector<double> real_part;
  std::vector<double> imag_part;
  
  wavelength.reserve(10000);
  real_part.reserve(10000);
  real_part.reserve(10000);

  while (std::getline(file, line))
  {
    std::istringstream input(line);

    double mu, a, b;

    input >> mu >> a >> b;

    wavelength.push_back(mu);
    real_part.push_back(a);
    imag_part.push_back(b);
  }

  file.close();


  std::vector<double> real_part_interpol = spectral_grid->interpolateToWavelengthGrid(
    wavelength, 
    real_part, 
    false);
  std::vector<double> imag_part_interpol = spectral_grid->interpolateToWavelengthGrid(
    wavelength, 
    imag_part, 
    true);

  refractive_index.reserve(nb_spectral_points);

  for (size_t i=0; i<nb_spectral_points; ++i)
    refractive_index.push_back(std::complex<double>(real_part_interpol[i], std::abs(imag_part_interpol[i])*-1.));
}



void DustSpecies::calcTransportCoefficients(
  const size_t radius_idx,
  std::vector<double>& absorption_coeff,
  std::vector<double>& scattering_coeff)
{
  std::vector<double> absorption_efficiency;
  std::vector<double> scattering_efficiency;
  std::vector<double> asymmetry_parameter;

  opticalProperties(
    particle_radius[radius_idx][0], 
    absorption_efficiency, 
    scattering_efficiency, 
    asymmetry_parameter);

  absorption_coeff.assign(nb_spectral_points, 0.);
  scattering_coeff.assign(nb_spectral_points, 0.);


  for (size_t i=0; i<nb_spectral_points; ++i)
  {
    absorption_coeff[i] = 2. * constants::pi * particle_radius[radius_idx][0]*particle_radius[radius_idx][0] 
                          * absorption_efficiency[i]
                          * number_density[radius_idx];
    scattering_coeff[i] = 2. * constants::pi * particle_radius[radius_idx][0]*particle_radius[radius_idx][0] 
                          * scattering_efficiency[i]
                          * number_density[radius_idx];
  }

}




void DustSpecies::opticalProperties(
  double radius, 
  std::vector<double>& absorption_efficiency, 
  std::vector<double>& scattering_efficiency, 
  std::vector<double>& asymmetry_parameter)
{
  absorption_efficiency.assign(nb_spectral_points, 0.);
  scattering_efficiency.assign(nb_spectral_points, 0.);
  asymmetry_parameter.assign(nb_spectral_points, 0.);


  for (size_t i=0; i<nb_spectral_points; ++i)
  {
    const double size_parameter = 2. * constants::pi * radius/spectral_grid->wavelength_list_cm[i];
    
    double qext;

    lxmie::Mie(
      refractive_index[i], 
      size_parameter, 
      qext, 
      scattering_efficiency[i], 
      absorption_efficiency[i], 
      asymmetry_parameter[i]);
  }


}


}