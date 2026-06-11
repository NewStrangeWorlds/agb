#include "hydrodynamics.h"

#include "../additional/aux_functions.h"
#include "../spectral_grid/spectral_grid.h"
#include "../config/config.h"
#include "../additional/exceptions.h"
#include "../atmosphere/atmosphere.h"
#include "../additional/quadrature.h"
#include "../additional/physical_const.h"
#include "../radiative_transfer/radiative_transfer.h"
#include "structure_solver.h"

#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>


namespace agb {


Hydrodynamics::Hydrodynamics(
  ModelConfig* config_,
  SpectralGrid* spectral_grid_,
  Atmosphere* atmosphere_,
  std::vector<RadiationField>& radiation_field_)
  : config(config_)
  , spectral_grid(spectral_grid_)
  , atmosphere(atmosphere_)
  , radiation_field(radiation_field_)
  , nb_grid_points(atmosphere->nb_grid_points)
{
  isothermal_sound_speed.assign(nb_grid_points, 0);
  alpha.assign(nb_grid_points, 0);
  sound_speed_derivative.assign(nb_grid_points, 0);
  phi.assign(nb_grid_points, 0);
}


void Hydrodynamics::saveOutput(const std::string file_path)
{
  std::fstream file;
  file.open(file_path.c_str(), std::ios::out);

  std::cout << "Saving hydrodynamic structure to " << file_path << "\n\n";

  file << std::setw(16) << std::left << "#r/R*" << "\t"
       << std::setw(16) << std::left << "v(cm/s)" << "\t"
       << std::setw(16) << std::left << "c_T(cm/s)" << "\t"
       << std::setw(16) << std::left << "alpha" << "\t"
       << std::setw(16) << std::left << "phi" << "\n";
      
  for (size_t i=0; i<nb_grid_points; ++i)
  {
    file << std::setprecision(10) << std::scientific
         << atmosphere->radius_grid[i] << "\t"
         << atmosphere->velocity[i] << "\t" 
         << isothermal_sound_speed[i] << "\t" 
         << alpha[i] << "\t" 
         << phi[i] << "\n";
  }

  file.close();
}



void Hydrodynamics::calcWindVelocity()
{
  isothermal_sound_speed.assign(nb_grid_points, 0);
  alpha.assign(nb_grid_points, 0);
  sound_speed_derivative.assign(nb_grid_points, 0);
  phi.assign(nb_grid_points, 0);

  last_solve_rejected = false;

  speedOfSound();

  std::vector<double> flux_weighted_extinction(nb_grid_points, 0);

  for (size_t i=0; i<flux_weighted_extinction.size(); ++i)
    flux_weighted_extinction[i] = radiation_field[i].fluxWeightedExtinction(atmosphere->extinction_coeff[i]);
  
  calcAlpha(flux_weighted_extinction);

  sound_speed_derivative = soundSpeedDerivative();

  int critical_point = findCriticalPoint();

  //dumpCriticalPointDebug(flux_weighted_extinction, critical_point);

  //a valid transonic wind solution requires the critical point to lie in the
  //interior. If the search hit a boundary (e.g. a transiently super-Eddington
  //alpha during the iteration), the wind equation has no transonic solution and
  //calcPhi/windVelocity would produce a (super-luminal) garbage velocity. Keep
  //the previous structure for this step instead of poisoning the iteration.
  if (critical_point <= 0 || critical_point >= static_cast<int>(nb_grid_points) - 1)
  {
    std::cout << "Hydrodynamics: no interior critical point (index " << critical_point
              << ") - skipping velocity/density update this iteration\n";
    last_solve_rejected = true;
    return;
  }

  phi = calcPhi(critical_point);

  std::vector<double> wind_velocity = windVelocity(phi, critical_point);

  //Eigenvalue mass-loss rate implied by the critical-point regularity condition.
  //In Winters (1994) this would be the quantity used to fix the stellar radius.
  const double mass_loss_rate_eigenvalue = calcMassLossRate(
    critical_point,
    sound_speed_derivative[critical_point],
    flux_weighted_extinction[critical_point]);

  //Henyey path (default; config->use_henyey_solver): predict the mass-loss rate
  //self-consistently as a transonic eigenvalue (Setup A). Seed Mdot from the
  //shooting estimate above and warm-start Phi from the previous solve. On a
  //rejected solve we do NOT fall back to shooting: solveStructureHenyey keeps the
  //last good structure and the next outer iteration retries from the same
  //warm-start. (Overwriting with the shooting result on reject was observed to
  //knock the solution onto a different, non-settling Mdot branch -- shooting is
  //exactly the unstable method we are replacing.)
  if (config->use_henyey_solver)
  {
    last_solve_rejected = !solveStructureHenyey(
      flux_weighted_extinction, mass_loss_rate_eigenvalue);
    return;
  }

  //--- validation hook (set HENYEY_CHECK=1): on the FIRST call (converged RT,
  //physical chi_F), cross-check the CppAD Henyey eigenvalue solver against
  //calcPhi/calcMassLossRate on the real (chi_F, c_T) snapshot, WITHOUT modifying
  //the structure. The wind velocity is seeded from a transonic guess based on the
  //real sound speed (independent of the shooting solution); dust is decoupled so
  //this isolates the wind + Mdot eigenvalue.
  static bool henyey_check_done = false;

  bool chi_f_finite = true;
  for (const double x : flux_weighted_extinction)
    if (!std::isfinite(x)) chi_f_finite = false;

  if (std::getenv("HENYEY_CHECK") && !henyey_check_done
      && mass_loss_rate_eigenvalue > 0. && chi_f_finite)
  {
    henyey_check_done = true;

    StructureSolver solver(nb_grid_points);

    for (size_t i=0; i<nb_grid_points; ++i)
    {
      solver.radius[i] = atmosphere->radius[i];
      solver.sound_speed[i] = isothermal_sound_speed[i];
      solver.flux_mean_extinction[i] = flux_weighted_extinction[i];
      solver.nucleation_per_h[i] = 0.;          //decouple dust for the wind check
      solver.growth_timescale[i] = 1.;
    }
    solver.sound_speed2_deriv = sound_speed_derivative;  //d(c_T^2)/dr
    solver.stellar_luminosity = config->stellar_luminosity;
    solver.gravitational_mass = constants::gravitation_const * config->stellar_mass;
    solver.inner_velocity = 0.3 * isothermal_sound_speed[0];

    const size_t M = nb_grid_points;
    const size_t N = (1 + StructureSolver::nb_moments)*M + 1;

    //transonic seed based on the real sound speed (sub-sonic base -> super-sonic)
    std::vector<double> x(N, 0.);
    for (size_t i=0; i<M; ++i)
    {
      const double frac = 0.3 + 1.5*double(i)/double(M);
      const double v = frac * isothermal_sound_speed[i];
      x[i] = 0.5*(v + isothermal_sound_speed[i]*isothermal_sound_speed[i]/v);
    }
    x[(1 + StructureSolver::nb_moments)*M] = std::log(mass_loss_rate_eigenvalue);

    double final_res = 0.; int iters = 0;
    const std::vector<double> sol = solver.solveEigen(x, 300, 1e-8, &final_res, &iters);
    const double mdot_eigen = std::exp(sol[(1 + StructureSolver::nb_moments)*M]);

    std::cout << "[HENYEY_CHECK] solve rel-residual " << final_res << " in " << iters
              << " iters; critical point " << solver.critical_point << "/" << M << "\n"
              << "[HENYEY_CHECK] Mdot eigen " << mdot_eigen
              << " g/s  vs calcMassLossRate " << mass_loss_rate_eigenvalue
              << " g/s  (ratio " << mdot_eigen/mass_loss_rate_eigenvalue << ")\n";
  }

  //EXPERIMENT (Winters 1994): prescribe Mdot as a fixed model parameter and feed it
  //into the continuity equation, instead of deriving it from the critical point.
  //With the radius grid (hence R*) also fixed this formally over-constrains the
  //system; the eigenvalue/fixed ratio printed below measures the inconsistency
  //(how far the critical point is from being regular for this prescribed Mdot).
  //const double mass_loss_rate_cgs = config->stellar_mass_loss_rate;
  const double mass_loss_rate_cgs = mass_loss_rate_eigenvalue;

  //store the mass loss rate in solar masses per year for output/diagnostics
  mass_loss_rate = mass_loss_rate_cgs / constants::mass_sun * constants::year;


  //Under-relaxation of the velocity update. Via mass conservation the radiative
  //acceleration alpha ~ chi_F/rho ~ chi_F * v / Mdot, so an unrelaxed velocity is
  //an undamped feedback channel: when dust forms, alpha can jump from ~0 to >1 in
  //a single step and run away. We relax in log space (v spans many orders of
  //magnitude).
  //
  //The relaxation factor is adapted from the trend of the (log) residual between
  //the proposed and current velocity: when the residual shrinks the iteration is
  //converging and we accelerate (larger factor); when it grows we are in the
  //unstable transient and damp harder. This keeps the early dust-ignition step
  //stable while letting the late, near-converged steps move quickly.
  const double relaxation_min = 0.1;
  const double relaxation_max = 0.9;

  double velocity_residual = 0.;

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    const double v_new = wind_velocity[i];

    if (std::isfinite(v_new) && v_new > 0.)
      velocity_residual = std::max(
        velocity_residual,
        std::abs(std::log(v_new / atmosphere->velocity[i])));
  }

  if (prev_velocity_residual > 0.)
  {
    if (velocity_residual < prev_velocity_residual)
      velocity_relaxation = std::min(relaxation_max, velocity_relaxation * 1.3);
    else
      velocity_relaxation = std::max(relaxation_min, velocity_relaxation * 0.5);
  }

  prev_velocity_residual = velocity_residual;

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    const double v_new = wind_velocity[i];

    //keep the previous velocity for any point with a degenerate/invalid update
    if (std::isfinite(v_new) && v_new > 0.)
      atmosphere->velocity[i] = std::pow(atmosphere->velocity[i], 1. - velocity_relaxation)
                              * std::pow(v_new, velocity_relaxation);
  }

  std::cout << "Velocity relaxation factor: " << velocity_relaxation
            << "  (residual " << velocity_residual << ")\n";

  //the density follows from mass conservation (using the relaxed velocity) and
  //needs the CGS (g/s) mass loss rate
  std::vector<double> mass_density = massDensity(mass_loss_rate_cgs);

  std::cout << "Critical point: " << critical_point << "\n";
  std::cout << "Mass loss rate (fixed): " << mass_loss_rate
            << "  eigenvalue: " << mass_loss_rate_eigenvalue / constants::mass_sun * constants::year
            << "  eigenvalue/fixed: " << mass_loss_rate_eigenvalue / mass_loss_rate_cgs << "\n";
}



void Hydrodynamics::setDustState(
  const std::vector<double>& nucleation_rate,
  const std::vector<double>& growth_timescale,
  const std::vector<double>& second_moment)
{
  dust_nucleation_rate = nucleation_rate;
  dust_growth_timescale = growth_timescale;
  dust_second_moment = second_moment;
}



void Hydrodynamics::speedOfSound()
{
  //for (size_t i=0; i<nb_grid_points; ++i)
    //isothermal_sound_speed[i] = sqrt(atmosphere->pressure[i]/atmosphere->mass_density[i]);

  for (size_t i=0; i<nb_grid_points; ++i)
    isothermal_sound_speed[i] = std::sqrt(1.0 / (atmosphere->mean_molecuar_weight[i] * constants::mass_proton) 
      * constants::boltzmann_k * atmosphere->temperature_gas[i]);
}


void Hydrodynamics::calcAlpha(
  const std::vector<double>& flux_weighted_extinction)
{
  for (size_t i=0; i<nb_grid_points; ++i)
  {
    alpha[i] = flux_weighted_extinction[i]/atmosphere->mass_density[i] 
             * config->stellar_luminosity
             /(4*constants::pi * constants::light_c * constants::gravitation_const * config->stellar_mass);
  }
}



std::vector<double> Hydrodynamics::windVelocity(
  const std::vector<double>& phi,
  const int critical_point)
{
  std::vector<double> wind_velocity(nb_grid_points, 0);

  for (size_t i=0; i<critical_point; ++i)
  {
     wind_velocity[i] = phi[i] - std::sqrt(phi[i]*phi[i] - isothermal_sound_speed[i] * isothermal_sound_speed[i]);
    
    if (i > 0 && wind_velocity[i] < wind_velocity[i-1])
      wind_velocity[i] = wind_velocity[i-1];
  }
  
  
  for (size_t i=critical_point+1; i<nb_grid_points; ++i)
    wind_velocity[i] = phi[i] + std::sqrt(phi[i]*phi[i] - isothermal_sound_speed[i] * isothermal_sound_speed[i]);

  
  //to make the transition between the upward and downward branch smooth
  //we interpolate the wind velocity at the critical point between the two end points
  const double x_left = atmosphere->radius[critical_point-1];
  const double x_right = atmosphere->radius[critical_point+1];
  
  double y_left = wind_velocity[critical_point-1];
  double y_right = wind_velocity[critical_point+1];
  wind_velocity[critical_point] = y_left + (atmosphere->radius[critical_point] - x_left) * (y_right - y_left) / (x_right - x_left);

  return wind_velocity;
}


//integration of phi from the inner boundary to the critical point
//uses an explicit Euler scheme
void Hydrodynamics::integratePhiOutward(
  const int critical_point,
  const double start_velocity,
  std::vector<double>& phi)
{ 
  std::vector<double> vel (nb_grid_points, 0);

  double v = start_velocity;
  phi[0] = 0.5 * (v + isothermal_sound_speed[0]*isothermal_sound_speed[0]/v);

  for (size_t i=1; i<critical_point+1; ++i)
  {
    const double r = atmosphere->radius[i-1];

    double phi_deriv = - 0.5/v 
             * constants::gravitation_const * config->stellar_mass 
             /r/r * (1 - alpha[i-1])
             + isothermal_sound_speed[i-1] * isothermal_sound_speed[i-1] 
             / v/r;

    phi[i] = phi[i-1] + (atmosphere->radius[i] - r) * phi_deriv;
    
    v = phi[i] - std::sqrt(phi[i]*phi[i] - isothermal_sound_speed[i]*isothermal_sound_speed[i]);
  }
 
}


//integration of phi from the outer boundary inwards to the critical point
//uses an explicit Euler scheme
void Hydrodynamics::integratePhiInward(
  const int critical_point,
  const double end_velocity,
  std::vector<double>& phi)
{
  double v = end_velocity;
  
  phi.back() = 0.5 * (v + isothermal_sound_speed.back()*isothermal_sound_speed.back() / v);

  for (int i=phi.size()-2; i>critical_point-1; --i)
  {
    const double r = atmosphere->radius[i+1];

    double phi_deriv = - 0.5/v 
             * constants::gravitation_const * config->stellar_mass 
             /r/r * (1 - alpha[i+1])
             + isothermal_sound_speed[i+1] * isothermal_sound_speed[i+1] 
             / v/r;

    phi[i] = phi[i+1] - (r - atmosphere->radius[i]) * phi_deriv;

    v = phi[i] + std::sqrt(phi[i]*phi[i] - isothermal_sound_speed[i]*isothermal_sound_speed[i]);
  }

}


std::vector<double> Hydrodynamics::calcPhi(const int critical_point)
{
  std::vector<double> phi(nb_grid_points, 0.0);

  //first we solve phi from the inner boundary to the critical point
  //we use a shooting method to find the correct velocity at inner boundary
  double v_0 = 1e-2;
  double v_1 = isothermal_sound_speed[critical_point];

  for (int it=0; it<1000; ++it)
  {
    double v_start = (v_1 - v_0) * 0.5 + v_0;
    
    integratePhiOutward(critical_point, v_start, phi);

    if (std::isnan(phi[critical_point]))
      v_1 = v_start;
    else
    {
      if (phi[critical_point] > isothermal_sound_speed[critical_point])
        v_0 = v_start;
      else 
        v_1 = v_start;
    }

    if (std::abs(v_1 - v_0)/v_1 < 1e-5 && !std::isnan(phi[critical_point]))
      break;
  }

  //now we solve phi from the outer boundary to the critical point
  //we use a shooting method to find the correct velocity at outer boundary
  v_0 = isothermal_sound_speed[critical_point];
  v_1 = isothermal_sound_speed[critical_point] * 100000;

  for (int it=0; it<1000; ++it)
  {
    double v_start = (v_1 - v_0) * 0.5 + v_0;
    
    integratePhiInward(critical_point, v_start, phi);

    if (std::isnan(phi[critical_point]))
      v_0 = v_start;
    else
    {
       if (phi[critical_point] < isothermal_sound_speed[critical_point])
         v_0 = v_start;
       else
        v_1 = v_start;
    }
    
    if (std::abs(v_1 - v_0)/v_1 < 1e-5  && !std::isnan(phi[critical_point]))
      break;
  }
  
  //set phi at the critical point to its actual value
  phi[critical_point] = isothermal_sound_speed[critical_point];

  return phi;
}



int Hydrodynamics::findCriticalPoint()
{
  std::vector<double> phi_deriv(nb_grid_points, 0);

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    phi_deriv[i] = constants::gravitation_const * config->stellar_mass 
                                / (atmosphere->radius[i]*atmosphere->radius[i])
                                * (1. - alpha[i])
                                - 2 * isothermal_sound_speed[i] * isothermal_sound_speed[i]/atmosphere->radius[i]
                                + sound_speed_derivative[i];
  }

  int i = 0;
  // 1) skip the inner region until radiative driving is significant  
  while (i < nb_grid_points-1 && alpha[i] < 0.3)
     ++i;
  // 2) from there, find the first regularity crossing
  while (i < nb_grid_points-1 && phi_deriv[i] > 0)
     ++i;

  return i;
}



void Hydrodynamics::dumpCriticalPointDebug(
  const std::vector<double>& flux_weighted_extinction,
  const int critical_point)
{
  ++wind_call_count;

  //full per-gridpoint dump when the critical point is anomalously inner or jumped
  //sharply from the previous call; otherwise just a one-line summary
  const bool jumped = prev_critical_point >= 0
                    && std::abs(critical_point - prev_critical_point) > 10;
  const bool inner  = critical_point >= 0 && critical_point < 30;

  const std::string path = config->model_folder + "hydro_debug.dat";
  std::ofstream file(path.c_str(), std::ios::app);

  if (file.is_open())
  {
    file << std::setprecision(8) << std::scientific;
    file << "# call " << wind_call_count
         << "  critical_point " << critical_point
         << "  (prev " << prev_critical_point << ")"
         << (jumped ? "  <<< JUMP" : "")
         << (inner  ? "  <<< INNER" : "") << "\n";

    if (jumped || inner)
    {
      file << "# i\tr/R*\tT_gas\trho\tcs2\talpha\tchi_F\tgrav\tsound1\tsound2\tphi_deriv\n";

      const double gm = constants::gravitation_const * config->stellar_mass;

      for (size_t k=0; k<nb_grid_points; ++k)
      {
        const double r   = atmosphere->radius[k];
        const double cs2 = isothermal_sound_speed[k]*isothermal_sound_speed[k];
        const double grav   = gm/(r*r) * (1. - alpha[k]);
        const double sound1 = -2.*cs2/r;
        const double sound2 = sound_speed_derivative[k];

        file << k << "\t" << atmosphere->radius_grid[k] << "\t"
             << atmosphere->temperature_gas[k] << "\t"
             << atmosphere->mass_density[k] << "\t"
             << cs2 << "\t" << alpha[k] << "\t" << flux_weighted_extinction[k] << "\t"
             << grav << "\t" << sound1 << "\t" << sound2 << "\t"
             << (grav + sound1 + sound2) << "\n";
      }
      file << "\n";
    }
  }

  prev_critical_point = critical_point;
}



std::vector<double> Hydrodynamics::soundSpeedDerivative()
{
  std::vector<double> sound_speed_derivative(nb_grid_points, 0);

  for (size_t i=0; i<nb_grid_points-1; ++i)
  {
    sound_speed_derivative[i] = (isothermal_sound_speed[i+1]*isothermal_sound_speed[i+1] 
                              - isothermal_sound_speed[i]*isothermal_sound_speed[i])
                              / (atmosphere->radius[i+1] - atmosphere->radius[i]);
  }

  sound_speed_derivative[nb_grid_points-1] = sound_speed_derivative[nb_grid_points-2];

  return sound_speed_derivative;
}


double Hydrodynamics::calcMassLossRate(
  const int critical_point,
  const double sound_speed_derivative,
  const double flux_weighted_extinction)
{
  const double radius = atmosphere->radius[critical_point];
  const double sound_speed = isothermal_sound_speed[critical_point];
 
  //the denominator is the critical-point regularity condition (= G M / r_c^2 * alpha),
  //and must use the same sign for the sound-speed derivative term as findCriticalPoint:
  //  G M / r^2 (1 - alpha) - 2 c_T^2 / r + d(c_T^2)/dr = 0
  double mass_loss_rate = config->stellar_luminosity * flux_weighted_extinction * sound_speed / constants::light_c
                          * 1./(constants::gravitation_const * config->stellar_mass/radius/radius
                                - 2*sound_speed*sound_speed/radius + sound_speed_derivative);

  return mass_loss_rate;
}


std::vector<double> Hydrodynamics::massDensity(
  const double mass_loss_rate)
{
  std::vector<double> mass_density(nb_grid_points, 0);

  //const double mass_loss_rate = config->stellar_mass_loss_rate;
  
  for (size_t i=0; i<nb_grid_points; ++i)
    mass_density[i] = mass_loss_rate 
      / (4*constants::pi*atmosphere->radius[i]*atmosphere->radius[i] * atmosphere->velocity[i]);

  return mass_density;
}

}