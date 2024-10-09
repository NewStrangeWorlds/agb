
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
#include "../dust/analytic_dust.h"
#include "../dust/gail_sedlmayr_dust.h"
#include "../hydrodynamics/hydrodynamics.h"



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
    DustSpecies* dust_species;
    TransportCoefficients transport_coeff;
    RadiativeTransfer radiative_transfer;
    TemperatureCorrection temperature_correction;
    Hydrodynamics hydrodynamics;
  protected:
    bool temperatureIteration();
    bool chemistryDustIteration();
    bool chemistryHydroIteration();
    void radiativeTransfer();
    std::pair<double, size_t> checkFluxConvergence();
    std::pair<double, size_t> checkConvergence(
      const std::vector<double>& old_data,
      const std::vector<double>& new_data);
    std::pair<double, size_t> checkEnergyBalance(
      std::vector<double>& temperature,
      std::vector<std::vector<double>>& absorption__coeff,
      std::vector<double>& deviation);
    std::pair<double, size_t> checkTemperatureConvergence(
      const std::vector<double>& temperature,
      const std::vector<double>& temperature_old,
      std::vector<double>& change);
    void forceMonotonicProfile(std::vector<double>& data);
    void smoothProfile(std::vector<double>& data);
};



}

#endif