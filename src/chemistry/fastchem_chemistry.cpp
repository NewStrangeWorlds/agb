
#include "fastchem_chemistry.h"

#include "chem_species.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"
#include "../../_deps/fastchem-src/fastchem_src/fastchem.h"


#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h> 
#include <iostream>


namespace agb {


FastChemChemistry::FastChemChemistry(
  const std::string& fastchen_parameter_file,
  const double metallicity_factor,
  const double co_ratio)
  : fastchem(fastchen_parameter_file, 1) 
{
  std::cout << "- Chemistry model: " << "equilibrium/FastChem" << "\n";
  std::cout << "  - Parameter file: " << fastchen_parameter_file << "\n\n";
  
  //the original, unscaled abundances from the FastChem input file
  reference_element_abundances = fastchem.getElementAbundances();

  fastchem_species_indices.assign(constants::species_data.size(), fastchem::FASTCHEM_UNKNOWN_SPECIES);

  for (size_t i=0; i<constants::species_data.size(); ++i)
    fastchem_species_indices[i] = fastchem.getSpeciesIndex(constants::species_data[i].fastchem_symbol);
    
  
  //check if C, O, and H are present in FastChem
  if (fastchem_species_indices[_H] == fastchem::FASTCHEM_UNKNOWN_SPECIES 
      || fastchem_species_indices[_O] == fastchem::FASTCHEM_UNKNOWN_SPECIES 
      || fastchem_species_indices[_C] == fastchem::FASTCHEM_UNKNOWN_SPECIES)
  {
    std::string error_message = "Critical elements (H, C, or O) not found in FastChem\n";
    throw InvalidInput(std::string ("FastChemChemistry::FastChemChemistry"), error_message);
  }

  //set metallicity and element abundances
  element_abundances = reference_element_abundances;
  
  for (size_t i=0; i<element_abundances.size(); ++i)
    if (i != fastchem_species_indices[_H] && (i != fastchem_species_indices[_He] && fastchem_species_indices[_He] != fastchem::FASTCHEM_UNKNOWN_SPECIES) )
      element_abundances[i] *= metallicity_factor;

  element_abundances[fastchem_species_indices[_C]] = element_abundances[fastchem_species_indices[_O]] * co_ratio;


  fastchem.setElementAbundances(element_abundances);

  //get back the new, normalised element abundances
  element_abundances = fastchem.getElementAbundances();
}


void FastChemChemistry::calcChemicalComposition(
  const std::vector<double>& parameters,
  const std::vector<double>& temperature,
  const std::vector<double>& pressure,
  const std::vector<double>& degree_of_condensation_c,
  std::vector<std::vector<double>>& number_densities,
  std::vector<double>& mean_molecular_weight,
  std::vector<double>& total_element_density,
  std::vector<double>& total_h_density)
{
  //in case we want to change element abundances on-the-fly
  if (parameters.size() == 2)
  {
    const double metallicity_factor = parameters[0];
    const double co_ratio = parameters[1];

    //set metallicity and element abundances
    element_abundances = reference_element_abundances;
  
    for (size_t i=0; i<element_abundances.size(); ++i)
      if (i != fastchem_species_indices[_H] && (i != fastchem_species_indices[_He] && fastchem_species_indices[_He] != fastchem::FASTCHEM_UNKNOWN_SPECIES) )
        element_abundances[i] *= metallicity_factor;

    element_abundances[fastchem_species_indices[_C]] = element_abundances[fastchem_species_indices[_O]] * co_ratio;

    fastchem.setElementAbundances(element_abundances);

    //get back the new, normalised element abundances
    element_abundances = fastchem.getElementAbundances();
  }


  mean_molecular_weight.assign(temperature.size(), 0);
  number_densities.assign(temperature.size(), std::vector<double>(constants::species_data.size(), 0));
  total_element_density.assign(temperature.size(), 0);
  total_h_density.assign(temperature.size(), 0);


  //run FastChem point-by-point since the element abundances might change
  for (size_t r=0; r<temperature.size(); ++r)
  {
    //adjust carbon element abundances
    std::vector<double> element_abundances_cond = element_abundances;

    element_abundances_cond[fastchem_species_indices[_C]] *= (1 - degree_of_condensation_c[r]);

    fastchem.setElementAbundances(element_abundances_cond);

    //set up the input & output structures and run the chemistry
    fastchem::FastChemInput input;
    fastchem::FastChemOutput output;

    input.temperature = std::vector<double>(1, temperature[r]);
    input.pressure = std::vector<double>(1, pressure[r]);

    size_t status = fastchem.calcDensities(input, output);

    if (status == fastchem::FASTCHEM_INITIALIZATION_FAILED)
    {
      std::string error_message = "FastChem initialisation failed!\n";
      throw InvalidInput(std::string ("FastChemChemistry::calcChemicalComposition"), error_message);
    }

    if (status != fastchem::FASTCHEM_SUCCESS)
      std::cout << "Something went wrong in FastChem :O!\n";
    
    mean_molecular_weight[r] = output.mean_molecular_weight[0];

    for (size_t i=0; i<constants::species_data.size(); ++i)
    {
      if (fastchem_species_indices[i] != fastchem::FASTCHEM_UNKNOWN_SPECIES)
        number_densities[r][i] = output.number_densities[0][fastchem_species_indices[i]];

      if (i == _TOTAL)
        number_densities[r][_TOTAL] = pressure[r] * 1.e6 / constants::boltzmann_k / temperature[r];
    }

    total_element_density[r] = output.total_element_density[0];

    total_h_density[r] = element_abundances_cond[fastchem_species_indices[_H]] * total_element_density[r];
  }

  fastchem.setElementAbundances(element_abundances);
}


}

