
#include <iostream>
#include <omp.h>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "transport_coeff.h"

#include "species_definition.h"
#include "../spectral_grid/spectral_grid.h"
#include "../chemistry/chem_species.h"
#include "../config/config.h"


namespace agb{


TransportCoefficients::TransportCoefficients(
  ModelConfig* config_ptr,
  SpectralGrid* grid_ptr, 
  const std::vector<std::string>& opacity_species_symbol,
  const std::vector<std::string>& opacity_species_folder)
{
  config = config_ptr;
  spectral_grid = grid_ptr;


  std::vector<size_t> spectral_indices = spectral_grid->spectralIndexList();
  
  gas_species.reserve(opacity_species_symbol.size());

  bool all_species_added = true;

  for (size_t i=0; i<opacity_species_symbol.size(); ++i)
  {
    bool added = addOpacitySpecies(opacity_species_symbol[i], opacity_species_folder[i]);

    if (!added) all_species_added = false;
  }

  
  std::cout << "\nOpacity specied added:\n";
  for (auto & i : gas_species)
  {
    std::cout << i->species_name << "\t" << i->species_folder << "\tfound: ";
    if (i->dataAvailable()) 
      std::cout << "yes";
    else
      std::cout << "no";
    
    std::cout << "\n\n";
  }
    
    
  if (!all_species_added) 
    std::cout << "Warning, not all opacities species from the model config file could be added!\n\n";
}




bool TransportCoefficients::addOpacitySpecies(
  const std::string& species_symbol, const std::string& species_folder)
{
  //first, the species cases for which separate classes are available
  if (species_symbol == "CIA-H2-H2")
  {
    gas_species.push_back(
      new GasGeneric(config, spectral_grid, _H2, "CIA H2-H2", species_folder, std::vector<size_t>{_H2}));

    return true;
  }

  //H2 Rayleigh scattering
  if (species_symbol == "H2")
  {
    gas_species.push_back(new GasH2(config, spectral_grid, ""));
   
    return true;
  }

  //He Rayleigh scattering
  if (species_symbol == "He")
  {
    gas_species.push_back(new GasHe(config, spectral_grid, ""));
   
    return true;
  }


  //H Rayleigh scattering
  if (species_symbol == "H")
  {
    gas_species.push_back(new GasH(config, spectral_grid, ""));

    return true;
  }


  if (species_symbol == "CIA-H2-He")
  {
    gas_species.push_back(
      new GasGeneric(config, spectral_grid, _H2, "CIA H2-He", species_folder, std::vector<size_t>{_He}));
   
    return true;
  }


  if (species_symbol == "CIA-H-He")
  {
    gas_species.push_back(
      new GasGeneric(config, spectral_grid, _H, "CIA H-He", species_folder, std::vector<size_t>{_He}));

    return true;
  }

  //H- free-free and bound-free continuum
  if (species_symbol == "H-")
  {
    gas_species.push_back(new GasHm(config, spectral_grid));

    return true;
  }


  //H2- free-free  continuum
  if (species_symbol == "H2-")
  {
    gas_species.push_back(new GasH2m(config, spectral_grid));

    return true;
  }


  //now we try the generic ones
  for (size_t i=0; i<constants::species_data.size(); ++i)
  {
    if (constants::species_data[i].symbol == species_symbol)
    {
      gas_species.push_back(
        new GasGeneric(
          config, spectral_grid, constants::species_data[i].id, constants::species_data[i].symbol, species_folder));

      return true;
    } 
  }


  //now we try the generic ones
  for (size_t i=0; i<constants::species_data.size(); ++i)
  {
    if (constants::species_data[i].symbol == species_symbol)
    {
      gas_species.push_back(
        new GasGeneric(
          config, spectral_grid, constants::species_data[i].id, constants::species_data[i].symbol, species_folder));

      return true;
    } 
  }


  //we haven't found the corresponding species
  std::cout << "Opacity species " 
    << species_symbol 
    << " has not been found in the internal list located in chem_species.h!\n";


  return false;
}


//calculates the transport coefficients on the CPU
//calls the calculation method of the individual opacity species
void TransportCoefficients::calculate(
  const double temperature,
  const double pressure,
  const std::vector<double>& number_densities,
  std::vector<double>& absorption_coeff,
  std::vector<double>& scattering_coeff)
{
  absorption_coeff.assign(spectral_grid->nbSpectralPoints(), 0);
  scattering_coeff.assign(spectral_grid->nbSpectralPoints(), 0);


  for (unsigned int i=0; i<gas_species.size(); i++)
    gas_species[i]->calcTransportCoefficients(
      temperature,
      pressure,
      number_densities,
      absorption_coeff,
      scattering_coeff);
}



TransportCoefficients::~TransportCoefficients()
{
  for (unsigned int i=0; i<gas_species.size(); ++i)
    delete gas_species[i];
}



}
