 
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

    std::vector<double> calculate(
      std::vector<double>& temperature,
      std::vector<double>& radius,
      std::vector<std::vector<double>>& extinction_coeff,
      std::vector<std::vector<double>>& absorption_coeff);
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

    std::vector<double> unsoeldLucyCorrection(
      const std::vector<double>& temperature,
      const std::vector<double>& radius,
      const std::vector<double>& fq_j,
      const std::vector<double>& qx_h,
      const std::vector<double>& kappa_j,
      const std::vector<double>& kappa_b,
      const std::vector<double>& planck_function_int);

    std::vector<double> lambdaIteration(
      const std::vector<double>& temperature,
      const std::vector<std::vector<double>>& absorption_coeff);
};



}

#endif