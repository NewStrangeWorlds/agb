 
#ifndef _radiative_transfer_h
#define _radiative_transfer_h

#include <vector>
#include <string>
#include <iostream>

#include "../additional/tri_diagonal_matrix.h"

namespace agb {

//forward declaration
class ModelConfig;
class SpectralGrid;
class Atmosphere;



struct AnglePoint
{
  AnglePoint(const size_t nb_spectral_points)
  {
    u.assign(nb_spectral_points, 0.0);
    v.assign(nb_spectral_points, 0.0);
    intensity_in.assign(nb_spectral_points, 0.0);
    intensity_out.assign(nb_spectral_points, 0.0);
  }

  std::vector<double> u;
  std::vector<double> v;
  double angle;
  double quadrature_weight;
  std::vector<double> intensity_out;
  std::vector<double> intensity_in;
};


struct zPoint
{
  double z = 0;
  double radius = 0;
  size_t radius_index = 0;
  AnglePoint* angle_point = nullptr;
};



class ImpactParam{
  public:
    std::vector<zPoint> z_grid;
    double p = 0;
    size_t nb_z_points = 0;

    void solveRadiativeTransfer(
      const size_t nu,
      const bool use_spline_discretisation,
      const double boundary_planck_derivative,
      const double boundary_flux_correction,
      const std::vector<double>& extinction_coeff,
      const std::vector<double>& source_function);
  private:
    void opticalDepth(
      const std::vector<double>& extinction_coeff,
      std::vector<double>& optical_depth);
    void assembleSystemTaylor(
      const std::vector<double>& optical_depth,
      const std::vector<double>& source_function,
      const double boundary_planck_derivative,
      const double boundary_flux_correction,
      const double boundary_exctinction_coeff,
      aux::TriDiagonalMatrix& M,
      std::vector<double>& rhs);
    void assembleSystemSpline(
      const std::vector<double>& optical_depth,
      const std::vector<double>& source_function,
      const double boundary_planck_derivative,
      const double boundary_flux_correction,
      const double boundary_exctinction_coeff,
      aux::TriDiagonalMatrix& M,
      std::vector<double>& rhs);
    void calcV(
      const size_t nu,
      const std::vector<double>& u,
      const std::vector<double>& optical_depth,
      const std::vector<double>& source_function);
};



struct RadiationField
{
  RadiationField(SpectralGrid* spectral_grid_);

  std::vector<double> eddington_flux;
  std::vector<double> flux;
  std::vector<double> mean_intensity;
  std::vector<double> eddington_k;

  double eddington_flux_int = 0;
  double flux_int = 0;
  double mean_intensity_int = 0;
  double eddington_k_int = 0;
  //Conservative (eq. 2.59) frequency-integrated flux, set by conservativeFluxIntegral().
  //Used where the flux VALUE vs the target matters (UL flux term, flux-convergence check);
  //eddington_flux_int itself stays = wavelengthIntegration(eddington_flux) so the flux-mean
  //ratios that divide by it (fluxWeightedExtinction, qx_h) stay consistent with their
  //eq. 2.58 numerators. Defaults to eddington_flux_int when the conservative form is off.
  double eddington_flux_int_conservative = 0;

  std::vector<double> eddington_factor;
  std::vector<double> sphericality_factor;

  std::vector<AnglePoint> angle_grid;
  std::vector<double> angles;

  std::vector<double> mean_intensity_impact;
  std::vector<double> eddington_flux_impact;
  std::vector<double> eddington_k_impact;

  size_t radius_index = 0;
  size_t nb_angles = 0;

  SpectralGrid* spectral_grid;

  void createAngleGrid(
    std::vector<ImpactParam>& impact_parameter_grid,
    const double radius,
    const size_t radius_index_);
  
  void angularIntegration();
  void wavelengthIntegration();
  double wavelengthIntegration(const std::vector<double>& data);
  double fluxWeightedExtinction(
    const std::vector<double>& extinction_coeff);
};



class RadiativeTransfer{
  public:
    RadiativeTransfer(
      ModelConfig* config_,
      SpectralGrid* spectral_grid_,
      Atmosphere* atmosphere_);
    ~RadiativeTransfer() {}

    void solveRadiativeTransfer();

    //Rebuild the radius-dependent ray geometry (impact parameters, z-grids and the
    //per-shell angle grids) after the radial grid has moved. Node count is
    //unchanged, so only positions are refreshed.
    void rebuildGeometry();

    //Full-linearisation temperature correction: one Newton step on (T_gas, T_dust)
    //with the converged RT operator frozen. Eliminates the per-frequency mean
    //intensities by a Rybicki-type reduction and solves a dense 2D x 2D temperature
    //system, returning the (unrelaxed) temperature changes. Requires the Taylor
    //moment discretisation. Must be called after solveRadiativeTransfer().
    void linearisedTemperatureCorrection(
      std::vector<double>& delta_temperature_gas,
      std::vector<double>& delta_temperature_dust);

    //Diagnostic (env FLUX_CONSIST): compare the node-centred flux from calcFlux with a
    //face-centred (conservative) flux r^2 H_{i+1/2} = d(a J)/dX whose divergence equals
    //the solveMomentSystem second-order stencil exactly. Prints the deviation of each
    //r^2 H_int from constancy (the flux-conservation error of the two discretisations).
    void fluxConsistencyDiagnostic();

    std::vector<RadiationField> radiation_field;
    void saveSpectrum(const std::string file_path);
  protected:
    ModelConfig* config;
    SpectralGrid* spectral_grid;
    Atmosphere* atmosphere;

    const size_t nb_core_impact_param = 0;
    const size_t nb_impact_param = 0;

    const size_t nb_grid_points = 0;
    const size_t nb_spectral_points = 0;

    std::vector<ImpactParam> impact_parameter_grid;

    std::vector<std::vector<double>> eddington_factor_f;
    std::vector<double> boundary_eddington_factor_h;
    std::vector<std::vector<double>> sphericality_factor;

    std::vector<std::vector<double>> extinction_coeff;
    std::vector<std::vector<double>> scattering_coeff;

    //Iteration-invariant thermal emission (absorption-weighted Planck) per
    //spectral point and grid point: absorption_gas*B(T_gas) + absorption_dust*B(T_dust).
    //Temperature and opacities do not change during a single RT solve, so this is
    //computed once per solveRadiativeTransfer instead of every iteration.
    //Indexed [nu][grid_point].
    std::vector<std::vector<double>> planck_emission;

    void createImpactParameterGrid();
    void createZGrids();

    void precomputeEmission();

    void calcEddingtonFactors();
    void calcSphericalityFactor();
    
    double boundaryFluxCorrection();
    const std::vector<double>& sourceFunction(const int nu);

    void solveMomentSystem(
      const size_t nu,
      const std::vector<double>& radius,
      const std::vector<double>& radius2,
      const double boundary_planck_derivative,
      const double boundary_flux_correction);
    const std::vector<double>& generateXGrid(
      const std::vector<double>& extinction_coeff,
      const std::vector<double>& radius,
      const std::vector<double>& sphericality_factor);
    void assembleMomentSystemSpline(
      const std::vector<double>& x_grid,
      const std::vector<double>& radius,
      const std::vector<double>& radius2,
      const std::vector<double>& emission_coeff,
      const std::vector<double>& extinction_coeff,
      const std::vector<double>& scattering_coeff,
      const std::vector<double>& eddington_factor,
      const double boundary_eddington_h,
      const std::vector<double>& sphericality_factor,
      const double boundary_planck_derivative,
      const double boundary_flux_correction,
      aux::TriDiagonalMatrix& m,
      std::vector<double>& rhs);
    void assembleMomentSystemTaylor(
      const std::vector<double>& x_grid,
      const std::vector<double>& radius,
      const std::vector<double>& radius2,
      const std::vector<double>& emission_coeff,
      const std::vector<double>& extinction_coeff,
      const std::vector<double>& scattering_coeff,
      const std::vector<double>& eddington_factor,
      const double boundary_eddington_h,
      const std::vector<double>& sphericality_factor,
      const double boundary_planck_derivative,
      const double boundary_flux_correction,
      aux::TriDiagonalMatrix& m,
      std::vector<double>& rhs);
    void calcFlux(
      const size_t nu,
      const std::vector<double>& radius,
      const std::vector<double>& radius2,
      const std::vector<double>& source_function);

    //Fill eddington_flux_int_conservative with the conservative (thesis eq. 2.59) flux:
    //r^2 H_int = L/(16 pi^2) + integral r^2 [int kappa_abs(B-J) dnu] dr. Consistent with the
    //moment-equation J-solve, so r^2 H_int is conserved at radiative equilibrium. Kept
    //separate from eddington_flux_int (and the per-frequency eddington_flux), which stay on
    //the eq. 2.58 transport form for the spectrum and the flux-mean ratios.
    void conservativeFluxIntegral();
    void assembleMomentSystemFluxSpline(
      const std::vector<double>& x_grid,
      const std::vector<double>& radius2,
      const std::vector<double>& source_function,
      const std::vector<double>& mean_intensity,
      const std::vector<double>& eddington_factor,
      const std::vector<double>& sphericality_factor,
      aux::TriDiagonalMatrix& m,
      std::vector<double>& rhs);
    void assembleMomentSystemFluxTaylor(
      const std::vector<double>& x_grid,
      const std::vector<double>& radius2,
      const std::vector<double>& source_function,
      const std::vector<double>& mean_intensity,
      const std::vector<double>& eddington_factor,
      const std::vector<double>& sphericality_factor,
      aux::TriDiagonalMatrix& m,
      std::vector<double>& rhs);

    std::pair<double, std::pair<size_t, size_t>> checkConvergence(
      const std::vector<std::vector<double>>& old_values,
      const std::vector<std::vector<double>>& new_values);
};



}

#endif