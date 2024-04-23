
#ifndef TRANSPORT_COEFF_H
#define TRANSPORT_COEFF_H

#include <vector>

#include "opacity_species.h"
#include "../spectral_grid/spectral_grid.h"
#include "../config/config.h"


namespace agb {


class TransportCoefficients {
  public:
    TransportCoefficients(
      ModelConfig* config_ptr,
      SpectralGrid* grid_ptr, 
      const std::vector<std::string>& opacity_species_symbol,
      const std::vector<std::string>& opacity_species_folder);
    ~TransportCoefficients();

    void calculate(
      const double temperature,
      const double pressure,
      const std::vector<double>& number_densities,
      std::vector<double>& absorption_coeff,
      std::vector<double>& scattering_coeff);

  private:
    ModelConfig* config = nullptr;
    SpectralGrid* spectral_grid = nullptr;

    std::vector<OpacitySpecies*> gas_species;

    bool addOpacitySpecies(
      const std::string& species_symbol, const std::string& species_folder);
};


}

#endif
