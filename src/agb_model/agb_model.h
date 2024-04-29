
#ifndef _agb_model_h
#define _agb_model_h

#include <iostream>
#include <vector>
#include <string>

#include "../config/config.h"
#include "../spectral_grid/spectral_grid.h"
#include "../atmosphere/atmosphere.h"
#include "../chemistry/fastchem_chemistry.h"
#include "../transport_coefficients/transport_coeff.h"
#include "../radiative_transfer/radiative_transfer.h"
#include "../temperature_correction/temperature_correction.h"
#include "../dust/dust_species.h"

namespace agb {


class AGBStarModel{
  public:
    AGBStarModel(const std::string folder);
    ~AGBStarModel() {}

    void calcModel();

    ModelConfig config;
    SpectralGrid spectral_grid;
    
    Atmosphere atmosphere;
    FastChemChemistry chemistry;
    DustSpecies dust_species;
    TransportCoefficients transport_coeff;
    RadiativeTransfer radiative_transfer;
    TemperatureCorrection temperature_correction;
  protected:
};



}

#endif