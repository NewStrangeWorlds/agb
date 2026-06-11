 
#include "atmosphere.h"


#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <omp.h>
#include <stdlib.h>
#include <cmath>

#include "../config/config.h"
#include "../additional/physical_const.h"
#include "../additional/exceptions.h"
#include "../additional/interpolation.h"

namespace agb {

Atmosphere::Atmosphere(
  ModelConfig* config_,
  const size_t nb_spectral_points)
 : config(config_)
{
  if (config->starting_model_path == "grey")
  {
    //generate the radial grid and a static hydrostatic grey (Lucy) starting structure
    //(Dominik et al. 1990 / Winters 1994): no file needed
    buildGreyStart();
  }
  else
  {
    std::string file_path = config->model_folder + config->starting_model_path;
    readStructure(file_path);
  }

  std::cout << "\nInitial atmosphere structure:\n";
  for (size_t i=0; i<nb_grid_points; ++i)
  {
    std::cout << std::setprecision(5) << std::scientific  
              << radius_grid[i]
              << "\t" << radius[i] 
              << "\t" << mass_density[i] 
              << "\t" << pressure[i] 
              << "\t" << temperature_gas[i] 
              << "\t" << temperature_dust[i] << "\t" 
              << velocity[i] << "\n";
  }

  std::cout << "\n";

  nb_grid_points = radius_grid.size();

  number_densities.assign(nb_grid_points, std::vector<double>(nb_chemistry_species, 0.));
  mean_molecuar_weight.assign(nb_grid_points, 0.);
  total_element_density.assign(nb_grid_points, 0.);
  total_h_density.assign(nb_grid_points, 0.);

  absorption_coeff.assign(
    nb_grid_points, 
    std::vector<double>(nb_spectral_points, 0.));

  scattering_coeff.assign(
    nb_grid_points, 
    std::vector<double>(nb_spectral_points, 0.));

  extinction_coeff.assign(
    nb_grid_points, 
    std::vector<double>(nb_spectral_points, 0.));

  absorption_coeff_gas.assign(
    nb_grid_points, 
    std::vector<double>(nb_spectral_points, 0.));

  scattering_coeff_gas.assign(
    nb_grid_points, 
    std::vector<double>(nb_spectral_points, 0.));

  extinction_coeff_gas.assign(
    nb_grid_points, 
    std::vector<double>(nb_spectral_points, 0.));

  absorption_coeff_dust.assign(
    nb_grid_points, 
    std::vector<double>(nb_spectral_points, 0.));

  scattering_coeff_dust.assign(
    nb_grid_points, 
    std::vector<double>(nb_spectral_points, 0.));

  extinction_coeff_dust.assign(
    nb_grid_points, 
    std::vector<double>(nb_spectral_points, 0.));
}



//Static hydrostatic grey (Lucy) starting model (Dominik et al. 1990; Winters 1994, App. A),
//used when starting_model = "grey": no structure file needed. A logarithmic radial grid
//spans grid_inner_radius..grid_outer_radius (in R_*); on it we solve the coupled grey system
//  - grey extinction          chi(r) = kappa_gas * rho(r)
//  - geometrically diluted optical depth (Winters Eq. A.2/A.3)
//                             dtau_L/dr = (R_*/r)^2 chi,  tau_L(r_out) = 0
//  - Lucy spherical grey T    (Winters Eq. A.1, with J/P ~ 1, and 2 W = 1 - sqrt(1-(R*/r)^2))
//                             T^4 = T_eff^4 (W(r) + 3/4 tau_L)
//  - radial optical depth     dtau/dr = chi,  tau(r_out) = 0          (for hydrostatic balance)
//  - hydrostatic balance      dP/dtau = g_eff/kappa,  g_eff = (GM/r^2)(1-alpha),
//                             alpha = kappa_gas L/(4 pi c G M), integrated inward from P=0
//  - ideal gas                rho = P mu m_H/(k T)
//The map rho->tau->P->rho is homogeneous in rho, so the scale is fixed physically by the
//photospheric anchor (Winters Eq. A.4): the diluted optical depth reaches 2/3 at r_in = R_*,
//i.e. tau_L(R_*) = 2/3. Rescaling rho each step pins it. A density floor keeps the (unphys-
//ically thin) outer hydrostatic atmosphere numerically benign for chemistry/RT; the full
//chemistry/dust/hydro/temperature cycle turns this seed into the dust-driven wind.
void Atmosphere::buildGreyStart()
{
  const size_t M     = config->grid_nb_points;
  const double Rstar = config->stellar_radius;                 //R_* [cm]
  const double r_in  = config->grid_inner_radius * Rstar;
  const double r_out = config->grid_outer_radius * Rstar;
  const double L     = config->stellar_luminosity;             //[erg/s]
  const double GM    = constants::gravitation_const * config->stellar_mass;
  const double kappa = config->grey_gas_opacity;               //[cm^2/g]
  const double mu    = 1.3;                                    //initial mean molecular weight
  const double rho_match = 3.0e-14;                            //hydrostatic -> wind-tail handover

  radius.assign(M, 0.);
  radius_grid.assign(M, 0.);
  const double dlog = std::log(r_out / r_in) / static_cast<double>(M - 1);
  for (size_t i=0; i<M; ++i)
  {
    radius[i]      = r_in * std::exp(dlog * static_cast<double>(i));
    radius_grid[i] = radius[i] / r_in;
  }
  nb_grid_points = M;

  const double t_eff = std::pow(
    L / (4.*constants::pi * Rstar*Rstar * constants::stefan_boltzmann), 0.25);

  //alpha is constant for a grey gas opacity (kappa_gas const): radiative acceleration in
  //units of gravity = kappa_gas L / (4 pi c G M). With dust absent here it is tiny.
  const double alpha = kappa * L / (4.*constants::pi * constants::light_c * GM);

  //dilution factor and an initial barometric density guess (one scale height at T_eff)
  std::vector<double> W(M, 0.), T(M, t_eff), rho(M, 0.);
  const double g_in = GM / (r_in*r_in) * (1. - alpha);
  const double H    = constants::boltzmann_k * t_eff / (mu * constants::mass_proton * g_in);
  const double rho0 = 1.0e-9;   //rough photospheric guess; the tau_L=2/3 anchor sets the scale
  for (size_t i=0; i<M; ++i)
  {
    W[i]   = 0.5*(1. - std::sqrt(std::max(1. - (Rstar/radius[i])*(Rstar/radius[i]), 0.)));
    rho[i] = rho0 * std::exp(-(radius[i] - r_in)/H);
  }

  const double tau_phot = 2./3.;
  std::vector<double> tau(M, 0.), tau_L(M, 0.), pres(M, 0.);
  for (int iter=0; iter<500; ++iter)
  {
    //radial and geometrically diluted optical depths, inward from the outer edge (=0 at i=M-1)
    tau[M-1]   = 0.;
    tau_L[M-1] = 0.;
    for (int i=static_cast<int>(M)-2; i>=0; --i)
    {
      const double chi_i  = kappa*rho[i];
      const double chi_ip = kappa*rho[i+1];
      const double dr     = radius[i+1] - radius[i];
      tau[i]   = tau[i+1]   + 0.5*(chi_i + chi_ip)*dr;
      const double w_i  = (Rstar/radius[i])  *(Rstar/radius[i]);     //(R_*/r)^2 dilution
      const double w_ip = (Rstar/radius[i+1])*(Rstar/radius[i+1]);
      tau_L[i] = tau_L[i+1] + 0.5*(w_i*chi_i + w_ip*chi_ip)*dr;
    }

    //photospheric anchor (Winters Eq. A.4): rescale rho so tau_L(r_in) = 2/3. tau, tau_L are
    //both linear in rho, so the same scale factor pins the whole structure.
    if (tau_L[0] > 0.)
    {
      const double scale = tau_phot / tau_L[0];
      for (size_t i=0; i<M; ++i) { rho[i] *= scale; tau[i] *= scale; tau_L[i] *= scale; }
    }

    //Lucy spherical grey temperature (Winters Eq. A.1) from the diluted optical depth
    for (size_t i=0; i<M; ++i)
      T[i] = t_eff * std::pow(W[i] + 0.75*tau_L[i], 0.25);

    //hydrostatic pressure: dP/dtau = g_eff/kappa, integrated inward from P(tau=0) = 0
    pres[M-1] = 0.;
    for (int i=static_cast<int>(M)-2; i>=0; --i)
    {
      const double ge_i  = GM/(radius[i]*radius[i])     * (1. - alpha);
      const double ge_ip = GM/(radius[i+1]*radius[i+1]) * (1. - alpha);
      pres[i] = pres[i+1] + 0.5*(ge_i + ge_ip)/kappa * (tau[i] - tau[i+1]);
    }

    //updated density from the ideal-gas law, under-relaxed for stability
    double max_change = 0.;
    for (size_t i=0; i<M; ++i)
    {
      const double rho_new = pres[i] * mu * constants::mass_proton
                           / (constants::boltzmann_k * T[i]);
      if (rho[i] > 0.)
        max_change = std::max(max_change, std::abs(rho_new - rho[i])/rho[i]);
      rho[i] = 0.5*rho[i] + 0.5*rho_new;
    }
    if (max_change < 1e-6) break;
  }

  //Beyond the photospheric layers the bare hydrostatic atmosphere becomes unphysically
  //thin (~120 orders of magnitude over the grid), which both is unphysical (the real outer
  //structure is the wind, rho ~ 1/r^2) and ill-conditions the RT moment system (constant or
  //vanishing extinction with large r^2 -> negative mean intensities). Where the hydrostatic
  //density first drops below rho_match, continue with a smooth rho ~ (r_t/r)^2 wind tail
  //anchored there: monotonic, RT-friendly and close to the real stratification. Pressure is
  //kept consistent with the tail density via the ideal-gas law. The wind solver replaces it.
  size_t i_t = M;
  for (size_t i=0; i<M; ++i) { if (rho[i] < rho_match) { i_t = i; break; } }
  if (i_t > 0 && i_t < M)
  {
    const double r_t   = radius[i_t-1];
    const double rho_t = rho[i_t-1];
    for (size_t i=i_t; i<M; ++i)
    {
      rho[i]  = rho_t * (r_t/radius[i]) * (r_t/radius[i]);
      pres[i] = rho[i] * constants::boltzmann_k * T[i] / (mu * constants::mass_proton);
    }
  }

  temperature_gas  = T;
  temperature_dust = T;
  mass_density     = rho;
  pressure         = pres;
  pressure_bar     = cgsToBar(pressure);

  //Transonic velocity seed from mass continuity (Winters Sect. 4.3 / Dominik et al. 1990:
  //the stationary initial model is a grey wind with the mass-loss rate as the eigenvalue).
  //With the prescribed Mdot, v(r) = Mdot/(4 pi r^2 rho) is fixed by the density: ~tens of
  //cm/s at the (dense) photosphere, rising through the sonic point to a ~km/s terminal value
  //in the rho ~ 1/r^2 tail. A non-zero velocity is essential because the Gail-Sedlmayr dust
  //moments are advected outward (dK/dr ~ 1/v); with the static v=0 seed no dust forms, the
  //radiative acceleration alpha never exceeds the driving threshold, and the stationary wind
  //solver finds no interior critical point. The hydro solver then refines v/rho/Mdot.
  velocity.assign(M, 0.);
  const double mdot = config->stellar_mass_loss_rate;   //cgs (g/s)
  if (mdot > 0.)
    for (size_t i=0; i<M; ++i)
      velocity[i] = mdot / (4.*constants::pi * radius[i]*radius[i] * rho[i]);

  std::cout << "Grey wind start: T_eff = " << t_eff
            << " K, T(R_*) = " << T[0] << " K, photospheric rho ~ " << rho[0]
            << " g/cm^3 (tau_L,inner = " << tau_L[0] << "), v(R_*) ~ " << velocity[0]
            << " cm/s, v_out ~ " << velocity[M-1] << " cm/s\n";
}


std::vector<double> Atmosphere::cgsToBar(const std::vector<double> pressure_data)
{
  std::vector<double> pressure_data_bar = pressure_data;

  for (auto & p : pressure_data_bar)
    p *= 1e-6;

  return pressure_data_bar;
}


void Atmosphere::remapToGrid(const std::vector<double>& new_radius)
{
  //interpolate the persistent state from the OLD radius (still in `radius`) onto
  //the new grid; positive quantities (rho, p) in log space
  temperature_gas  = aux::interpolate(radius, temperature_gas,  new_radius);
  temperature_dust = aux::interpolate(radius, temperature_dust, new_radius);
  velocity         = aux::interpolate(radius, velocity,         new_radius);
  mass_density     = aux::interpolate(radius, mass_density,     new_radius, true);
  pressure         = aux::interpolate(radius, pressure,         new_radius, true);

  //adopt the new grid (radius_grid is r normalised to the inner radius)
  radius = new_radius;
  for (size_t i=0; i<nb_grid_points; ++i)
    radius_grid[i] = radius[i] / radius[0];

  pressure_bar = cgsToBar(pressure);
}


void Atmosphere::equationOfState()
{
  for (size_t i=0; i<nb_grid_points; ++i)
    pressure[i] = mass_density[i] / (mean_molecuar_weight[i] * constants::mass_proton) * constants::boltzmann_k * temperature_gas[i];

  pressure_bar = cgsToBar(pressure);
}


}