
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
    TransportCoefficients transport_coeff;
  protected:
};



}

#endif