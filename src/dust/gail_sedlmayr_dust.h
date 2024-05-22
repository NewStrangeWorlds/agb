 
#ifndef _gail_sedlmayr_dust_h
#define _gail_sedlmayr_dust__h

#include <vector>
#include <string>
#include <iostream>
#include <complex>
#include <cmath>

#include "dust_species.h"


namespace agb {

//forward declaration
class ModelConfig;
class SpectralGrid;
class Atmosphere;


class GailSedlmayrDust: public DustSpecies{
  public:
    GailSedlmayrDust(
      ModelConfig* config_,
      SpectralGrid* spectral_grid_,
      Atmosphere* atmosphere_);
    ~GailSedlmayrDust() {}

    void calcDistribution();
    void calcTransportCoefficients(
      const size_t radius_idx,
      std::vector<double>& absorption_coeff,
      std::vector<double>& scattering_coeff);
  protected:
    const double monomer_radius = 1.28e-8;
};



}

#endif