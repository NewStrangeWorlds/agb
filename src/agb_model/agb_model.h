
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
    //persistent state for the radiative-equilibrium temperature iteration
    unsigned int temp_iter_count = 0;
    double prev_max_rel_change = 1e30;                         //last applied max |dT|/T
    bool anderson_was_active = false;                          //activation-edge tracker
    //two-phase corrector: begin with Unsoeld-Lucy, latch on the linearisation once it is
    //safe (see config.linearisation_start_unsoeld_lucy and the switch thresholds)
    bool linearisation_active = false;
    unsigned int linearisation_ready_count = 0;
    std::vector<double> relaxation_gas, relaxation_dust;       //per-layer omega
    std::vector<double> prev_delta_b_gas, prev_delta_b_dust;   //previous delta_B
    std::vector<std::vector<double>> anderson_x_gas, anderson_f_gas;   //T and residual history
    std::vector<std::vector<double>> anderson_x_dust, anderson_f_dust;

    bool temperatureIteration();
    bool chemistryDustIteration();
    bool chemistryHydroIteration();
    void radiativeTransfer();

    //movable-grid step: redistribute the radial nodes by equidistribution of a
    //monitor over flux-mean opacity, gas temperature and nucleation rate, remap the
    //structure, rebuild the RT geometry and reset stale per-node iteration state
    void applyMovableGrid();

    //Anderson mixing on a temperature profile. x_k is the current profile, f_k the
    //relaxed correction (so G(x_k) = x_k + f_k). Pushes the pair onto the history,
    //trims to the window, and returns the accelerated next profile.
    std::vector<double> andersonStep(
      std::vector<std::vector<double>>& x_history,
      std::vector<std::vector<double>>& f_history,
      const std::vector<double>& x_k,
      const std::vector<double>& f_k);

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