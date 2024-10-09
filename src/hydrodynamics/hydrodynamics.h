
#ifndef _hydrodynamics_h
#define _hydrodynamics_h


#include <string>
#include <vector>


namespace agb {

//forward declaration
class ModelConfig;
class SpectralGrid;
class Atmosphere;
class RadiationField;


class Hydrodynamics{
  public:
    Hydrodynamics(
      ModelConfig* config_,
      SpectralGrid* spectral_grid_,
      Atmosphere* atmosphere_,
      std::vector<RadiationField>& radiation_field_);
    ~Hydrodynamics() {}

    void calcWindVelocity();
    void saveOutput(const std::string file_path);

    std::vector<double> isothermal_sound_speed;
    std::vector<double> alpha;
    double mass_loss_rate = 0;
  private:
    ModelConfig* config;
    SpectralGrid* spectral_grid;
    Atmosphere* atmosphere;
    std::vector<RadiationField>& radiation_field;

    const size_t nb_grid_points = 0;
    std::vector<double> sound_speed_derivative;
    std::vector<double> phi;

    void speedOfSound();
    void calcAlpha(
      const std::vector<double>& flux_weighted_extinction);
    int findCriticalPoint();
    double findCriticalPoint(
      unsigned int& lower_index,
      unsigned int& upper_index);
    std::vector<double> soundSpeedDerivative();
    std::vector<double> calcPhi(
      const int critical_point);
    std::vector<double> windVelocity(
      const std::vector<double>& phi,
      const int critical_point);
    double calcMassLossRate(
      const int critical_point,
      const double sound_speed_derivative,
      const double flux_weighted_extinction);
    std::vector<double> massDensity(
      const double mass_loss_rate);

    void integratePhiOutward(
      const int critical_point,
      const double start_velocity,
      std::vector<double>& phi);
    void integratePhiFromCriticalPoint(
      const int critical_point,
      std::vector<double>& phi);
    void integratePhiInward(
      const int critical_point,
      const double end_velocity,
      std::vector<double>& phi);

    void calcPhiInward(
      const int critical_point, 
      const double critical_point_radius,
      std::vector<double>& phi);
    void calcPhiOutward(
      const int critical_point, 
      const double critical_point_radius,
      std::vector<double>& phi);
};

}
#endif 
