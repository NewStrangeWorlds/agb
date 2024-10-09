
#ifndef _fastchem_chemistry_h
#define _fastchem_chemistry_h


#include "../../_deps/fastchem-src/fastchem_src/fastchem.h"
#include "chem_species.h"

#include <string>
#include <vector>


namespace agb {


class FastChemChemistry{
  public:
    FastChemChemistry(
      const std::string& fastchen_parameter_file,
      const double metallicity,
      const double c_o_ratio);
    ~FastChemChemistry() {}

    std::vector<double> element_abundances;
    std::vector<size_t> fastchem_species_indices;

    void calcChemicalComposition(
      const std::vector<double>& parameters,
      const std::vector<double>& temperature,
      const std::vector<double>& pressure,
      const std::vector<double>& degree_of_condensation_c,
      std::vector<std::vector<double>>& number_densities,
      std::vector<double>& mean_molecular_weight,
      std::vector<double>& total_element_density,
      std::vector<double>& total_h_density);
  private:
    fastchem::FastChem<long double> fastchem;

    std::vector<chemical_species_id> species;
    std::vector<double> reference_element_abundances;
};


}
#endif 
