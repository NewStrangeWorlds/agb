
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

    //provide the (frozen) dust kernels so the Henyey solver can couple the
    //dust-moment equations: J* [1/(cm^3 s)] and growth timescale tau [s] per point
    void setDustState(
      const std::vector<double>& nucleation_rate,
      const std::vector<double>& growth_timescale,
      const std::vector<double>& second_moment);

    std::vector<double> isothermal_sound_speed;
    std::vector<double> alpha;
    double mass_loss_rate = 0;

    //true if the last calcWindVelocity could not update the structure (no interior
    //critical point, or a rejected Henyey solve). The structure is then frozen, so
    //a near-zero alpha change must NOT be mistaken for convergence by the caller.
    bool last_solve_rejected = false;
  private:
    //adaptive under-relaxation state for the velocity update (see calcWindVelocity)
    double velocity_relaxation = 0.25;
    double prev_velocity_residual = -1.;

    ModelConfig* config;
    SpectralGrid* spectral_grid;
    Atmosphere* atmosphere;
    std::vector<RadiationField>& radiation_field;

    const size_t nb_grid_points = 0;
    std::vector<double> sound_speed_derivative;
    std::vector<double> phi;

    //warm-start state for the Henyey structure solver (HENYEY_SOLVE path)
    std::vector<double> henyey_x_prev;
    std::vector<double> henyey_chi_prev;  //flux-mean extinction used last solve (for alpha damping)
    std::vector<double> dust_nucleation_rate;  //J* [1/(cm^3 s)] (frozen, from the dust module)
    std::vector<double> dust_growth_timescale; //tau [s]      (frozen, from the dust module)
    std::vector<double> dust_second_moment;    //K2 [cm^-3]   (frozen, dust feedback reference)

    //solve the wind structure + mass-loss eigenvalue with the CppAD Henyey
    //solver (dust decoupled; dust opacity enters via the frozen chi_F) and write
    //velocity/density/mass_loss_rate back into the atmosphere. Returns false if
    //the solve was rejected (non-finite / no interior critical point), leaving
    //the structure unchanged.
    bool solveStructureHenyey(
      const std::vector<double>& flux_weighted_extinction,
      const double seed_mass_loss_rate);

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
