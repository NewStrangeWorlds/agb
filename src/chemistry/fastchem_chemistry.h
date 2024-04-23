
#ifndef _fastchem_chemistry_h
#define _fastchem_chemistry_h


#include "../../_deps/fastchem-src/fastchem_src/fastchem.h"
#include "chem_species.h"

#include <string>
#include <vector>


namespace agb {


class FastChemChemistry{
  public:
    FastChemChemistry(const std::string& fastchen_parameter_file);
    ~FastChemChemistry() {}
    
    void calcChemicalComposition(
      const std::vector<double>& parameters,
      const std::vector<double>& temperature,
      const std::vector<double>& pressure,
      std::vector<std::vector<double>>& number_densities,
      std::vector<double>& mean_molecular_weight);
  private:
    fastchem::FastChem<long double> fastchem;

    std::vector<chemical_species_id> species;
    std::vector<double> reference_element_abundances;
    std::vector<size_t> fastchem_species_indices;
};


}
#endif 
