 
#ifndef _temperature_correction_h
#define _temperature_correction_h

#include <vector>
#include <string>
#include <iostream>

namespace agb {

//forward declaration
class ModelConfig;
class SpectralGrid;
class RadiationField;
class RadiativeTransfer;


class TemperatureCorrection{
  public:
    TemperatureCorrection(
      ModelConfig* config_,
      SpectralGrid* spectral_grid_,
      std::vector<RadiationField>& radiation_field_)
      : config(config_)
      , spectral_grid(spectral_grid_)
      , radiation_field(radiation_field_) {}
    ~TemperatureCorrection() {}

    //Computes one relaxed temperature correction (delta T) for the given profile.
    //relaxation and prev_delta_b carry the per-layer adaptive under-relaxation
    //factor omega and the previous (unrelaxed) Planck correction across calls;
    //they are resized/initialised on first use and updated in place.
    std::vector<double> calculate(
      std::vector<double>& temperature,
      std::vector<double>& radius,
      std::vector<std::vector<double>>& extinction_coeff,
      std::vector<std::vector<double>>& absorption_coeff,
      std::vector<double>& relaxation,
      std::vector<double>& prev_delta_b);

    //Full-linearisation (Newton) temperature correction (thesis 3.2.3 / App. B.1): the
    //radiative-equilibrium peer of calculate(). Uses RadiativeTransfer for the per-frequency
    //linearised moment operator (buildLinearisedMomentSystem) and owns the RE constraints,
    //the Rybicki elimination and the dense 2D x 2D Newton solve; returns the (relaxed but
    //un-capped) temperature changes. Requires the Taylor moment discretisation.
    void linearisedCorrection(
      RadiativeTransfer& radiative_transfer,
      const std::vector<double>& temperature_gas,
      const std::vector<double>& temperature_dust,
      const std::vector<double>& radius,
      const std::vector<std::vector<double>>& extinction_coeff,
      const std::vector<std::vector<double>>& absorption_coeff_gas,
      const std::vector<std::vector<double>>& absorption_coeff_dust,
      std::vector<double>& delta_temperature_gas,
      std::vector<double>& delta_temperature_dust);
  protected:
    ModelConfig* config;
    SpectralGrid* spectral_grid;
    std::vector<RadiationField>& radiation_field;

    void calcIntegratedQuantities(
      const std::vector<double>& radius,
      const std::vector<double>& temperature,
      const std::vector<std::vector<double>>& extinction_coeff,
      const std::vector<std::vector<double>>& absorption_coeff,
      std::vector<double>& fq_J,
      std::vector<double>& qx_H,
      std::vector<double>& kappa_j,
      std::vector<double>& kappa_b,
      std::vector<double>& planck_function_int);

    //Returns the Planck-integrated correction delta_B per layer (not yet converted
    //to a temperature change).
    std::vector<double> unsoeldLucyCorrection(
      const std::vector<double>& temperature,
      const std::vector<double>& radius,
      const std::vector<double>& fq_j,
      const std::vector<double>& qx_h,
      const std::vector<double>& kappa_j,
      const std::vector<double>& kappa_b,
      const std::vector<double>& planck_function_int);

    //1-2-1 filter applied to the correction (not the profile), so it vanishes at
    //convergence and cannot diffuse the converged solution.
    void smoothCorrection(std::vector<double>& correction);

    std::vector<double> lambdaIteration(
      const std::vector<double>& temperature,
      const std::vector<std::vector<double>>& absorption_coeff);
};



}

#endif