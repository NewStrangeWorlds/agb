 
#ifndef _analytic_dust_h
#define _analytic_dust__h

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


class AnalyticDust: public DustSpecies{
  public:
    AnalyticDust(
      ModelConfig* config_,
      SpectralGrid* spectral_grid_,
      Atmosphere* atmosphere_,
      const double particle_radius);
    ~AnalyticDust() {}

    void calcDistribution();
    void calcTransportCoefficients(
      const size_t radius_idx,
      std::vector<double>& absorption_coeff,
      std::vector<double>& scattering_coeff);
    std::vector<double> degreeOfCondensation(
      const double carbon_abundance);
    void saveOutput(const std::string file_path);
  protected:
    const double monomer_radius = 1.28e-8;
    const double const_radius = 0;
};



}

#endif