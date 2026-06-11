 
#include "radiative_transfer.h" 

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <omp.h>
#include <sstream>
#include <algorithm> 
#include <assert.h>
#include <cmath>

#include "../spectral_grid/spectral_grid.h"
#include "../config/config.h"
#include "../additional/physical_const.h"
#include "../atmosphere/atmosphere.h"
#include "../additional/tri_diagonal_matrix.h"
#include "../additional/aux_functions.h"


namespace agb{


void RadiativeTransfer::solveMomentSystem(
  const size_t nu,
  const std::vector<double>& radius,
  const std::vector<double>& radius2,
  const double boundary_planck_derivative,
  const double boundary_flux_correction)
{
  const std::vector<double>& x_grid = generateXGrid(
    extinction_coeff[nu],
    radius,
    sphericality_factor[nu]);

  //iteration-invariant thermal emission, precomputed once per RT solve in
  //precomputeEmission() (absorption_gas*B(T_gas) + absorption_dust*B(T_dust));
  //bit-identical to the previous inline evaluation since products commute
  const std::vector<double>& emission_coeff = planck_emission[nu];

  //per-thread reusable matrix / scratch (called inside omp-parallel-over-nu)
  static thread_local aux::TriDiagonalMatrix m(0);
  static thread_local std::vector<double> rhs;
  static thread_local std::vector<double> result;
  m.resize(nb_grid_points);
  rhs.assign(nb_grid_points, 0.);

  if (config->use_spline_discretisation)
    assembleMomentSystemSpline(
      x_grid,
      radius,
      radius2,
      emission_coeff,
      extinction_coeff[nu],
      scattering_coeff[nu],
      eddington_factor_f[nu],
      boundary_eddington_factor_h[nu],
      sphericality_factor[nu],
      boundary_planck_derivative,
      boundary_flux_correction,
      m,
      rhs);
  else
    assembleMomentSystemTaylor(
      x_grid,
      radius,
      radius2,
      emission_coeff,
      extinction_coeff[nu],
      scattering_coeff[nu],
      eddington_factor_f[nu],
      boundary_eddington_factor_h[nu],
      sphericality_factor[nu],
      boundary_planck_derivative,
      boundary_flux_correction,
      m,
      rhs);


  m.solveInto(rhs, result);

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    if (result[i] < 0)
    {
      for (size_t i=0; i<nb_grid_points; ++i)
        std::cout << "ms " << nu << "  " << i << "\t" << result[i] << "\t" << x_grid[i] << "\t" << scattering_coeff[nu][i] << "\t" << extinction_coeff[nu][i] << "\t" << emission_coeff[i] << "\t" << sphericality_factor[nu][i] << "\t" << m.a[i] << "\t" << m.b[i] << "\t" << m.c[i] << "\t" << rhs[i] << "\n";
      exit(0);
    }
  }

  for (size_t i=0; i<nb_grid_points; ++i)
    radiation_field[i].mean_intensity[nu] = result[i];
}


const std::vector<double>& RadiativeTransfer::generateXGrid(
  const std::vector<double>& extinction_coeff,
  const std::vector<double>& radius,
  const std::vector<double>& sphericality_factor)
{
  //per-thread reusable buffer (callers are inside omp-parallel-over-nu regions)
  static thread_local std::vector<double> x_grid;
  x_grid.resize(nb_grid_points);
  x_grid[nb_grid_points-1] = 0.;

  for (int i=nb_grid_points-2; i>-1; --i)
  {
    double integral = - (radius[i] - radius[i+1]) * (extinction_coeff[i]*sphericality_factor[i] + extinction_coeff[i+1]*sphericality_factor[i+1]) / 2.0;

    x_grid[i] = x_grid[i+1] + integral;
  }

  return x_grid;
}



//Spline discretisation
void RadiativeTransfer::assembleMomentSystemSpline(
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
  std::vector<double>& rhs)
{
  static thread_local std::vector<double> a;
  a.assign(nb_grid_points, 0.);

  for (size_t i=0; i<nb_grid_points; ++i)
    a[i] = eddington_factor[i] * sphericality_factor[i] * radius2[i];

  auto b = [&](size_t i) {
    return radius2[i]/sphericality_factor[i] * (1. - scattering_coeff[i]/extinction_coeff[i]);};

  auto c = [&](size_t i) {
    return radius2[i]/sphericality_factor[i] / extinction_coeff[i];};

  const double hr = x_grid[1] - x_grid.front();
  
  m.b[0] = hr/3. * b(0) + 1/hr * a[0];
  m.c[0] = hr/6. * b(1) - 1/hr * a[1];

  rhs[0] = hr/3. * c(0) * emission_coeff[0]
         + hr/6. * c(1) * emission_coeff[1]
         - radius2[0]/extinction_coeff[0] * boundary_flux_correction * boundary_planck_derivative;
  
  for (size_t i=1 ; i<nb_grid_points-1; ++i)
  {
    const double hr = x_grid[i+1] - x_grid[i];
	  const double hl = x_grid[i] - x_grid[i-1];

    m.a[i] =  hl/6. * b(i-1) - 1./hl * a[i-1];
	  m.b[i] = (hr+hl)/3. * b(i) + 1./hr * a[i] + 1./hl * a[i];
	  m.c[i] = hr/6. * b(i+1) - 1./hr * a[i+1];

	  rhs[i] = hl/6. * c(i-1) * emission_coeff[i-1]
	         + (hr+hl)/3. * c(i) * emission_coeff[i]
		       + hr/6. * c(i+1) * emission_coeff[i+1];
  }

  const double hl = x_grid.back() - x_grid[nb_grid_points-2];

  m.a.back() = hl/6. * b(nb_grid_points-2) - 1./hl * a[nb_grid_points-2];
  m.b.back() = hl/3. * b(nb_grid_points-1) + 1./hl * a[nb_grid_points-1] - radius2.back() * boundary_eddington_h;
  rhs.back() = hl/6. * c(nb_grid_points-2) * emission_coeff[nb_grid_points-2] + hl/3. * c(nb_grid_points-1) * emission_coeff[nb_grid_points-1];
}



void RadiativeTransfer::calcFlux(
  const size_t nu,
  const std::vector<double>& radius,
  const std::vector<double>& radius2,
  const std::vector<double>& source_function)
{
  //Monochromatic flux H_nu from the first-moment relation (thesis eq. 2.58),
  //d(f q r^2 J)/dX. This transport form keeps each H_nu physical (>= 0); the
  //frequency-integrated flux is made conservative separately in conservativeFluxIntegral()
  //(eq. 2.59) - the divergence form would conserve the integral but produce unphysical
  //negative monochromatic fluxes (e.g. where the cool outer dust absorbs starlight).
  const std::vector<double>& x_grid = generateXGrid(
    extinction_coeff[nu],
    radius,
    sphericality_factor[nu]);

  //per-thread reusable scratch (called inside omp-parallel-over-nu)
  static thread_local std::vector<double> mean_intensity;
  static thread_local aux::TriDiagonalMatrix m(0);
  static thread_local std::vector<double> rhs;
  static thread_local std::vector<double> result;
  mean_intensity.resize(nb_grid_points);
  m.resize(nb_grid_points);
  rhs.assign(nb_grid_points, 0.);

  for (size_t i=0; i<nb_grid_points; ++i)
    mean_intensity[i] = radiation_field[i].mean_intensity[nu];

  if (config->use_spline_discretisation)
    assembleMomentSystemFluxSpline(
      x_grid,
      radius2,
      source_function,
      mean_intensity,
      eddington_factor_f[nu],
      sphericality_factor[nu],
      m,
      rhs);
  else
    assembleMomentSystemFluxTaylor(
      x_grid,
      radius2,
      source_function,
      mean_intensity,
      eddington_factor_f[nu],
      sphericality_factor[nu],
      m,
      rhs);

  m.solveInto(rhs, result);

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    radiation_field[i].eddington_flux[nu] = result[i]/radius2[i];
    radiation_field[i].flux[nu] = 4. * constants::pi * radiation_field[i].eddington_flux[nu];
  }

}


//Conservative frequency-integrated flux (thesis eq. 2.59), used to OVERWRITE
//eddington_flux_int after the per-frequency calcFlux. The flux divergence equals the
//frequency-integrated local balance, so r^2 H_int is conserved at radiative equilibrium:
//  r^2 H_int(r) = L/(16 pi^2) + integral_{r1}^{r} r'^2 [int kappa_abs (B - J) dnu] dr'.
//Same sign convention as the linearisation's S_bal (which drives this to the target).
void RadiativeTransfer::conservativeFluxIntegral()
{
  const double target = config->stellar_luminosity / (16. * constants::pi * constants::pi);

  //per-node radial source  S_i = r^2 * wavelengthIntegration( kappa_g(J-B_g) + kappa_d(J-B_d) )
  std::vector<double> source(nb_grid_points, 0.);

  #pragma omp parallel for
  for (size_t i=0; i<nb_grid_points; ++i)
  {
    std::vector<double> y(nb_spectral_points, 0.);
    for (size_t nu=0; nu<nb_spectral_points; ++nu)
    {
      const double J  = radiation_field[i].mean_intensity[nu];
      const double wl = spectral_grid->wavelength_list[nu];
      const double Bg = aux::planckFunctionWavelength(atmosphere->temperature_gas[i],  wl);
      const double Bd = aux::planckFunctionWavelength(atmosphere->temperature_dust[i], wl);
      y[nu] = atmosphere->absorption_coeff_gas[i][nu]  * (J - Bg)
            + atmosphere->absorption_coeff_dust[i][nu] * (J - Bd);
    }
    const double r2 = atmosphere->radius[i]*atmosphere->radius[i];
    source[i] = r2 * radiation_field[i].wavelengthIntegration(y);
  }

  //cumulative trapezoid outward from the inner boundary r^2 H_int(r1) = L/(16 pi^2).
  //Stored in the separate eddington_flux_int_conservative so the (eq. 2.58) eddington_flux_int
  //- and the flux-mean ratios that divide by it - stay self-consistent.
  double r2h = target;
  radiation_field[0].eddington_flux_int_conservative =
    r2h / (atmosphere->radius[0]*atmosphere->radius[0]);

  for (size_t i=1; i<nb_grid_points; ++i)
  {
    r2h += 0.5*(source[i] + source[i-1])*(atmosphere->radius[i] - atmosphere->radius[i-1]);
    radiation_field[i].eddington_flux_int_conservative =
      r2h / (atmosphere->radius[i]*atmosphere->radius[i]);
  }
}


//spline discretisation
void RadiativeTransfer::assembleMomentSystemFluxSpline(
  const std::vector<double>& x_grid,
  const std::vector<double>& radius2,
  const std::vector<double>& source_function,
  const std::vector<double>& mean_intensity,
  const std::vector<double>& eddington_factor,
  const std::vector<double>& sphericality_factor,
  aux::TriDiagonalMatrix& m,
  std::vector<double>& rhs)
{
  static thread_local std::vector<double> a;
  a.assign(nb_grid_points, 0.);

  for (size_t i=0; i<nb_grid_points; ++i)
    a[i] = eddington_factor[i] * sphericality_factor[i] * radius2[i] * mean_intensity[i];
  
  const double hr = x_grid[1] - x_grid.front();
  
  m.c[0] = 2./hr;
  m.b[0] = 4./hr;
  rhs[0] = 6./hr/hr * (a[1] - a[0])
           - radius2[0]/sphericality_factor[0] * (mean_intensity[0] - source_function[0]); 

  for (size_t i=1; i<nb_grid_points-1; ++i)
  {
    const double hr = x_grid[i+1] - x_grid[i];
	  const double hl = x_grid[i] - x_grid[i-1];

    m.a[i] = 1./hl;
	  m.b[i] = 2. * (1./hl + 1./hr);
    m.c[i] = 1./hr;
    rhs[i] = 3./(hl*hl) * (a[i] - a[i-1]) + 3./(hr*hr) * (a[i+1] - a[i]);
  }

  const double hl = x_grid.back() - x_grid[nb_grid_points-2];

  m.a.back() = 2./hl;
  m.b.back() = 4./hl;
  rhs.back() = 6./(hl*hl) * (a.back() - a[nb_grid_points-2])
               + radius2.back()/sphericality_factor.back()
               * (mean_intensity.back() - source_function.back());
}



//Taylor discretisation
void RadiativeTransfer::assembleMomentSystemTaylor(
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
  std::vector<double>& rhs)
{
  const double hr = x_grid[1] - x_grid.front();

  static thread_local std::vector<double> a;
  a.assign(nb_grid_points, 0.);

  for (size_t i=0; i<nb_grid_points; ++i)
    a[i] = eddington_factor[i] * sphericality_factor[i] * radius2[i];

  auto b = [&](size_t i) {
    return radius2[i]/sphericality_factor[i] * (1 - scattering_coeff[i]/extinction_coeff[i]);};

  auto c = [&](size_t i) {
    return radius2[i]/sphericality_factor[i] / extinction_coeff[i];};

  m.b[0] = -1/hr * a[0] - hr/2 * b(0);
  m.c[0] = 1/hr * a[1];

  rhs[0] = radius2[0]/extinction_coeff[0] * boundary_flux_correction * boundary_planck_derivative
           - hr/2. * c(0) * emission_coeff[0];
  
  for (size_t i=1 ; i<nb_grid_points-1; ++i)
  {
    const double hr = x_grid[i+1] - x_grid[i];
	  const double hl = x_grid[i] - x_grid[i-1];

    m.a[i] = 2./(hl*(hl+hr)) * a[i-1];
	  m.b[i] = - 2./(hl*hr) * a[i] - b(i);
	  m.c[i] = 2./(hr*(hr+hl)) * a[i+1];

	  rhs[i] = - c(i) * emission_coeff[i];
  }

  const double hl = x_grid.back() - x_grid[nb_grid_points-2];

  m.a.back() = 1./hl * a[nb_grid_points-2];
  m.b.back() = -1./hl * a.back() - hl/2 * b(nb_grid_points-1) + radius2.back() * boundary_eddington_h;

  rhs.back() = -hl/2. * c(nb_grid_points-1) * emission_coeff.back();
}


void RadiativeTransfer::assembleMomentSystemFluxTaylor(
  const std::vector<double>& x_grid,
  const std::vector<double>& radius2,
  const std::vector<double>& source_function,
  const std::vector<double>& mean_intensity,
  const std::vector<double>& eddington_factor,
  const std::vector<double>& sphericality_factor,
  aux::TriDiagonalMatrix& m,
  std::vector<double>& rhs)
{
  static thread_local std::vector<double> a;
  a.assign(nb_grid_points, 0.);

  for (size_t i=0; i<nb_grid_points; ++i)
    a[i] = eddington_factor[i] * sphericality_factor[i] * radius2[i];
  
  const double hr = x_grid[1] - x_grid.front(); 
  
  m.c[0] = 0.;
  m.b[0] = 1.;
  
  rhs[0] = 1./hr * (a[1]*mean_intensity[1] - a[0]*mean_intensity[0])
           - hr/2. * radius2[0]/sphericality_factor[0] * (mean_intensity[0] - source_function[0]); 
  
  //rhs[0] = 6./hr/hr * (a[1] - a[0])
  //         - radius2[0]/sphericality_factor[0] * (mean_intensity[0] - source_function[0]); 

  for (size_t i=1; i<nb_grid_points-1; ++i)
  {
    const double hr = x_grid[i+1] - x_grid[i];
	  const double hl = x_grid[i] - x_grid[i-1];

    m.a[i] = 0;
	  m.b[i] = 1;
    m.c[i] = 0;
	  
    rhs[i] = - hr/(hl*(hl+hr)) * a[i-1] * mean_intensity[i-1]
	           + (hr-hl)/(hr*hl) * a[i] * mean_intensity[i]
			       + hl/(hr*(hr+hl)) * a[i+1] * mean_intensity[i+1];
  }

  const double hl = x_grid.back() - x_grid[nb_grid_points-2];

  m.a.back() = 0;
  m.b.back() = 1;
  rhs.back() = 1/hl*(a.back() * mean_intensity.back() - a[nb_grid_points-2] * mean_intensity[nb_grid_points-2])
	           + hl/2 * radius2.back()/sphericality_factor.back()
             * (mean_intensity.back() - source_function.back());
}


//RT half of the full linearisation (thesis eq. 3.51): assemble, for one frequency, the
//frozen Taylor moment operator V = M_nu, the diagonal thermal source dS_{s,i} =
//d rhs_i/d T_{s,i}, and the residual K_nu = rhs_nu - M_nu J_nu. The temperature-correction
//module consumes these (constraints, Rybicki elimination, dense Newton solve). The
//inner-boundary flux-injection term's T_gas[0] dependence is left frozen (not linearised).
void RadiativeTransfer::buildLinearisedMomentSystem(
  const size_t nu,
  const std::vector<double>& radius,
  const std::vector<double>& radius2,
  const double boundary_flux_correction,
  aux::TriDiagonalMatrix& m,
  std::vector<double>& source_gas,
  std::vector<double>& source_dust,
  std::vector<double>& residual_K)
{
  const size_t D = nb_grid_points;
  const double wl = spectral_grid->wavelength_list[nu];

  const double boundary_planck_derivative =
    aux::planckFunctionDerivWavelength(atmosphere->temperature_gas[0], wl);

  const std::vector<double>& x_grid = generateXGrid(
    extinction_coeff[nu], radius, sphericality_factor[nu]);

  //resize allocates the a/b/c diagonals AND the solveInto() scratch buffers the caller
  //needs for the Rybicki solves
  m.resize(D);
  static thread_local std::vector<double> rhs_nu;
  rhs_nu.assign(D, 0.);

  assembleMomentSystemTaylor(
    x_grid, radius, radius2,
    planck_emission[nu], extinction_coeff[nu], scattering_coeff[nu],
    eddington_factor_f[nu], boundary_eddington_factor_h[nu], sphericality_factor[nu],
    boundary_planck_derivative, boundary_flux_correction,
    m, rhs_nu);

  //diagonal source d rhs_i/d T_{s,i} = -c(i) kappa_abs,s dB/dT, c(i)=radius2/(q chi);
  //the boundary rows carry the extra hr/2, hl/2 factor of the Taylor assembly
  source_gas.assign(D, 0.);
  source_dust.assign(D, 0.);
  for (size_t i=0; i<D; ++i)
  {
    const double dB_gas  = aux::planckFunctionDerivWavelength(atmosphere->temperature_gas[i],  wl);
    const double dB_dust = aux::planckFunctionDerivWavelength(atmosphere->temperature_dust[i], wl);

    const double c_i = radius2[i] / sphericality_factor[nu][i] / extinction_coeff[nu][i];

    double face = 1.0;
    if (i == 0)       face = 0.5*(x_grid[1] - x_grid[0]);
    else if (i+1==D)  face = 0.5*(x_grid[D-1] - x_grid[D-2]);

    source_gas[i]  = - face * c_i * atmosphere->absorption_coeff_gas[i][nu]  * dB_gas;
    source_dust[i] = - face * c_i * atmosphere->absorption_coeff_dust[i][nu] * dB_dust;
  }

  //residual K = rhs_nu - M_nu J_nu (~0 when the RT is converged)
  residual_K.assign(D, 0.);
  residual_K[0]   = rhs_nu[0]   - (m.b[0]*radiation_field[0].mean_intensity[nu]
                                 + m.c[0]*radiation_field[1].mean_intensity[nu]);
  for (size_t i=1; i+1<D; ++i)
    residual_K[i] = rhs_nu[i]   - (m.a[i]*radiation_field[i-1].mean_intensity[nu]
                                 + m.b[i]*radiation_field[i].mean_intensity[nu]
                                 + m.c[i]*radiation_field[i+1].mean_intensity[nu]);
  residual_K[D-1] = rhs_nu[D-1] - (m.a[D-1]*radiation_field[D-2].mean_intensity[nu]
                                 + m.b[D-1]*radiation_field[D-1].mean_intensity[nu]);
}


//Diagnostic (env FLUX_CONSIST): contrast the node-centred conservative flux integral with a
//face-centred flux r^2 H_{i+1/2} = d(a J)/dX (a = f q r^2). Prints how far each r^2 H_int
//stays from the target L/(16 pi^2).
void RadiativeTransfer::fluxConsistencyDiagnostic()
{
  const size_t D = nb_grid_points;
  const size_t N = nb_spectral_points;
  const std::vector<double>& lambda = spectral_grid->wavelength_list;
  const double target = config->stellar_luminosity / (16. * constants::pi * constants::pi);

  std::vector<double> radius(D), radius2(D);
  for (size_t i=0; i<D; ++i) { radius[i] = atmosphere->radius[i]; radius2[i] = radius[i]*radius[i]; }

  std::vector<double> w_int(N, 0.);
  if (N >= 2)
  {
    w_int[0] = -0.5*(lambda[1]-lambda[0]);  w_int[N-1] = -0.5*(lambda[N-1]-lambda[N-2]);
    for (size_t n=1; n+1<N; ++n) w_int[n] = -0.5*(lambda[n+1]-lambda[n-1]);
  }

  std::vector<double> r2h_face(D-1, 0.);
  for (size_t n=0; n<N; ++n)
  {
    const std::vector<double>& x = generateXGrid(extinction_coeff[n], radius, sphericality_factor[n]);
    for (size_t i=0; i+1<D; ++i)
    {
      const double ai  = eddington_factor_f[n][i]   * sphericality_factor[n][i]   * radius2[i];
      const double ai1 = eddington_factor_f[n][i+1] * sphericality_factor[n][i+1] * radius2[i+1];
      const double Ji  = radiation_field[i].mean_intensity[n];
      const double Ji1 = radiation_field[i+1].mean_intensity[n];
      r2h_face[i] += w_int[n] * (ai1*Ji1 - ai*Ji) / (x[i+1] - x[i]);
    }
  }

  double node_lo=1e300, node_hi=-1e300, face_lo=1e300, face_hi=-1e300;
  for (size_t i=1; i+1<D; ++i)
  {
    const double node = radiation_field[i].eddington_flux_int_conservative * radius2[i];
    node_lo = std::min(node_lo, node); node_hi = std::max(node_hi, node);
  }
  for (size_t i=1; i+2<D; ++i)
  {
    face_lo = std::min(face_lo, r2h_face[i]); face_hi = std::max(face_hi, r2h_face[i]);
  }

  std::cout << "[FLUXCON] target r2H = " << target << "\n"
            << "[FLUXCON] conservative (eq2.59) r2H_int in ["
            << node_lo << ", " << node_hi << "]  spread = "
            << (node_hi-node_lo)/target*100. << " %\n"
            << "[FLUXCON] face-centred (eq2.59) r2H_int in ["
            << face_lo << ", " << face_hi << "]  spread = "
            << (face_hi-face_lo)/target*100. << " %\n";
}


}