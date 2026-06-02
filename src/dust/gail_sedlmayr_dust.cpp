
#include "gail_sedlmayr_dust.h" 

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
#include <array>

#include "../additional/aux_functions.h"
#include "../spectral_grid/spectral_grid.h"
#include "../config/config.h"
#include "../additional/exceptions.h"
#include "../atmosphere/atmosphere.h"
#include "../additional/quadrature.h"
#include "../additional/physical_const.h"
#include "../chemistry/chem_species.h"



namespace agb{


GailSedlmayrDust::GailSedlmayrDust(
  ModelConfig* config_,
  SpectralGrid* spectral_grid_,
  Atmosphere* atmosphere_)
  : DustSpecies(config_, spectral_grid_, atmosphere_)
{

}



double GailSedlmayrDust::nucleationRate(
  const double temperature,
  const double number_density_c,
  const double number_density_c2,
  const double number_density_c2h,
  const double number_density_c2h2)
{
  
  const double ln_saturation_ratio = saturationRatio(
    temperature,
    number_density_c);

  if (ln_saturation_ratio <= 0)
    return 1e-100;


  const double n_star = criticalClusterSize(
    temperature,
    ln_saturation_ratio);

  const double delta_f = freeEnergyOfFormation(
    temperature,
    ln_saturation_ratio,
    n_star);

  const double beta = monomerGrowthRate(
    temperature,
    number_density_c,
    number_density_c2,
    number_density_c2h,
    number_density_c2h2);

  const double z = zeldovichFactor(
    temperature,
    n_star,
    ln_saturation_ratio);

  const double c0 = equilibriumClusterDistribution(
    temperature,
    delta_f,
    number_density_c,
    ln_saturation_ratio);

  const double a_star = monomer_surface_area * std::pow(n_star, 2./3.);
  double nucleation_rate = beta * a_star * z * c0;

  if (nucleation_rate < 1e-100)
    nucleation_rate = 1e-100;

  return nucleation_rate;
}



double GailSedlmayrDust::growthRate(
  const double temperature,
  const double number_density_c,
  const double number_density_c2,
  const double number_density_c2h,
  const double number_density_c2h2)
{
  double tau = monomer_surface_area * std::sqrt(constants::boltzmann_k * temperature/(2* constants::pi))
              * ( sticking_coeff[0]/std::sqrt(mass_c) * number_density_c
              + 2 * sticking_coeff[1]/std::sqrt(mass_c2) * number_density_c2
              + 2 * sticking_coeff[1]/std::sqrt(mass_c2h) * number_density_c2h
              + 2 * sticking_coeff[1]/std::sqrt(mass_c2h2) * number_density_c2h2);

  return 1./tau;
}


void GailSedlmayrDust::calcDistribution(const double condensable_carbon_abundance)
{
  const size_t M = nb_grid_points;

  //Gather the (full, un-depleted) gas-phase inputs on the radius grid.
  std::vector<double> nc(M), nc2(M), nc2h(M), nc2h2(M);
  std::vector<double> temp(M), vel(M), n_h(M), tau_full(M), ccond_per_h(M);

  for (size_t i=0; i<M; ++i)
  {
    nc[i]    = atmosphere->number_densities[i][_C];
    nc2[i]   = atmosphere->number_densities[i][_C2];
    nc2h[i]  = atmosphere->number_densities[i][_C2H];
    nc2h2[i] = atmosphere->number_densities[i][_C2H2];
    temp[i]  = atmosphere->temperature_gas[i];
    vel[i]   = atmosphere->velocity[i];
    n_h[i]   = atmosphere->total_h_density[i];

    //growthRate returns the growth timescale tau for the un-depleted gas
    tau_full[i] = growthRate(temp[i], nc[i], nc2[i], nc2h[i], nc2h2[i]);

    //condensable carbon per hydrogen nucleus, (eps_C - eps_O) * n_tot / n_H, with all
    //oxygen assumed locked in CO. The degree of condensation is fc = K3 / ccond_per_h
    //(Winters Eq. 5.4); the growth species available for further condensation is
    //throttled by (1 - fc).
    ccond_per_h[i] = (condensable_carbon_abundance > 0. && n_h[i] > 0.)
      ? condensable_carbon_abundance * atmosphere->total_element_density[i] / n_h[i]
      : 0.;
  }

  auto lin = [](const double a, const double b, const double f)
    { return a + (b - a)*f; };
  auto logterp = [](const double a, const double b, const double f)
    { return (a > 0. && b > 0.)
        ? std::exp(std::log(a)*(1.-f) + std::log(b)*f) : a + (b - a)*f; };

  //Right-hand side of the coupled moment ODE system dK_i/dr at one point, with the
  //growth species density throttled by the RUNNING degree of condensation fc=K3/ccond.
  //Both the growth rate (1/tau, linear in the carbon density) and the nucleation rate
  //J* (strongly nonlinear in the supersaturation) are evaluated from the depleted
  //densities, so condensation self-limits as the carbon is consumed. Because the sweep
  //marches outward and K3 only grows, fc rises monotonically and the throttle decreases
  //smoothly -- the resource limit is differential and causal (Gail & Sedlmayr book
  //Sec. 14.3), with no outer fixed-point iteration and hence no oscillation.
  auto rhs = [&](const double c, const double c2, const double c2h, const double c2h2,
                 const double t, const double v, const double nh,
                 const double tau, const double ccond,
                 const std::array<double,6>& K, std::array<double,6>& dK)
  {
    double depl = (ccond > 0.) ? 1. - K[3]/ccond : 0.;
    if (depl < 0.) depl = 0.;
    if (depl > 1.) depl = 1.;

    const double inv_tau_eff = (tau > 0.) ? depl/tau : 0.;
    const double j_star = nucleationRate(t, c*depl, c2*depl, c2h*depl, c2h2*depl);
    const double j_per_h = (nh > 0.) ? j_star/nh : 0.;
    const double v_inv = (v > 0.) ? 1./v : 0.;

    dK[0] = j_per_h * v_inv;
    for (int o=1; o<6; ++o)
      dK[o] = v_inv * ( (o/3.) * inv_tau_eff * K[o-1]
        + std::pow(static_cast<double>(minimum_monomer_number), o/3.) * j_per_h );
  };

  //Integrate the six moments (normalised to n<H>) outward with a coupled RK4 sweep.
  dust_moments.assign(nb_moments, std::vector<double>(M, 0.));

  std::array<double,6> K;
  K.fill(0.);  //no dust at the inner boundary

  for (size_t i=1; i<M; ++i)
  {
    const double h = atmosphere->radius[i] - atmosphere->radius[i-1];

    std::array<double,6> k1, k2, k3, k4, kt;

    rhs(nc[i-1], nc2[i-1], nc2h[i-1], nc2h2[i-1], temp[i-1], vel[i-1], n_h[i-1],
        tau_full[i-1], ccond_per_h[i-1], K, k1);

    const double nc_m   = logterp(nc[i-1],   nc[i],   0.5);
    const double nc2_m  = logterp(nc2[i-1],  nc2[i],  0.5);
    const double nc2h_m = logterp(nc2h[i-1], nc2h[i], 0.5);
    const double nc2h2_m= logterp(nc2h2[i-1],nc2h2[i],0.5);
    const double t_m    = lin(temp[i-1], temp[i], 0.5);
    const double v_m    = lin(vel[i-1],  vel[i],  0.5);
    const double nh_m   = logterp(n_h[i-1], n_h[i], 0.5);
    const double tau_m  = logterp(tau_full[i-1], tau_full[i], 0.5);
    const double cc_m   = lin(ccond_per_h[i-1], ccond_per_h[i], 0.5);

    for (int j=0; j<6; ++j) kt[j] = K[j] + 0.5*h*k1[j];
    rhs(nc_m, nc2_m, nc2h_m, nc2h2_m, t_m, v_m, nh_m, tau_m, cc_m, kt, k2);

    for (int j=0; j<6; ++j) kt[j] = K[j] + 0.5*h*k2[j];
    rhs(nc_m, nc2_m, nc2h_m, nc2h2_m, t_m, v_m, nh_m, tau_m, cc_m, kt, k3);

    for (int j=0; j<6; ++j) kt[j] = K[j] + h*k3[j];
    rhs(nc[i], nc2[i], nc2h[i], nc2h2[i], temp[i], vel[i], n_h[i],
        tau_full[i], ccond_per_h[i], kt, k4);

    for (int j=0; j<6; ++j) K[j] += h/6.*(k1[j] + 2.*k2[j] + 2.*k3[j] + k4[j]);
    for (int j=0; j<6; ++j) dust_moments[j][i] = K[j];
  }

  //scale the n<H>-normalised moments back to absolute units
  for (size_t i=0; i<M; ++i)
    for (auto& k : dust_moments)
    {
      k[i] *= atmosphere->total_h_density[i];
      if (k[i] < 1e-45) k[i] = 1e-45;
    }

  //store the (depleted) nucleation rate and growth timescale actually felt at each
  //grid point, for output and as the frozen kernels for the hydrodynamics Newton solver
  nucleation_rate.assign(M, 0.);
  growth_rate.assign(M, 0.);
  for (size_t i=0; i<M; ++i)
  {
    const double ccond_dim = condensable_carbon_abundance * atmosphere->total_element_density[i];
    double depl = (ccond_dim > 0.) ? 1. - dust_moments[3][i]/ccond_dim : 0.;
    if (depl < 0.) depl = 0.;
    if (depl > 1.) depl = 1.;

    nucleation_rate[i] = nucleationRate(temp[i], nc[i]*depl, nc2[i]*depl, nc2h[i]*depl, nc2h2[i]*depl);
    growth_rate[i]     = growthRate(temp[i], nc[i]*depl, nc2[i]*depl, nc2h[i]*depl, nc2h2[i]*depl);
  }

  number_density.assign(M, 0.);
  size_distribution.assign(M, std::vector<double>(1, 1.0));
  particle_radius.assign(M, std::vector<double>(1, 0.)); //in cm

  std::cout << "Gail&Sedlmayr dust calculation\n";

  for (size_t i=0; i<M; ++i)
  {
    const double ln_s = saturationRatio(temp[i], nc[i]);

    if (ln_s < 0)
    {
      number_density[i] = 0.;
      particle_radius[i][0] = 0.;
    }
    else
    {
      number_density[i] = dust_moments[0][i];

      //particle radius from the mean condensed volume per grain
      particle_radius[i][0] = monomer_radius * std::pow(dust_moments[3][i] / dust_moments[0][i], 1./3.);
    }

    //floor the number density and radius (otherwise RT/temperature calculations crash)
    if (number_density[i] < 1e-30)
    {
      number_density[i] = 1e-30;
      particle_radius[i][0] = 1e-8;
    }
  }

  std::cout << "\n";
}



std::vector<double> GailSedlmayrDust::degreeOfCondensation(
  const double condensable_carbon_abundance)
{
  std::vector<double> degree_of_condensation(atmosphere->nb_grid_points, 0);

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    //fc = (carbon condensed into grains) / (condensable carbon). The condensable
    //carbon abundance (eps_C - eps_O, all O locked in CO) is given with respect to
    //n<tot>, so the carbon number density is total_element_density * condensable.
    degree_of_condensation[i] =
      dust_moments[3][i] /(atmosphere->total_element_density[i] * condensable_carbon_abundance);

    if (degree_of_condensation[i] >=1) degree_of_condensation[i] = 0.9999999;
  }

  return degree_of_condensation;
}


void GailSedlmayrDust::saveOutput(const std::string file_path)
{
  std::fstream file;
  file.open(file_path.c_str(), std::ios::out);

  if (file.fail()) 
  {
    std::cout << "Couldn't open dust output file " << file_path << "\n";
    return;
  }

  std::cout << "Saving dust output to " << file_path << "\n\n";

  file << std::setprecision(10) << std::scientific << "#r/R*\tn_<H>(cm-3)\tnumber_density(cm-3)\tradius(micron)\tnucleation_rate\tgrowth_rate\tK0\tK1\tK2\tK3\n";
  
  for (size_t i=0; i<nb_grid_points; ++i)
  {  
     const double number_density_c = atmosphere->number_densities[i][_C];
     const double number_density_c2 = atmosphere->number_densities[i][_C2];
     const double number_density_c2h = atmosphere->number_densities[i][_C2H];
     const double number_density_c2h2 = atmosphere->number_densities[i][_C2H2];
     const double temperature = atmosphere->temperature_gas[i];

     file << atmosphere->radius_grid[i] << "\t"
          << atmosphere->total_h_density[i] << "\t"
          << temperature << "\t"
          << number_density_c << "\t"
          << number_density_c2 << "\t"
          << number_density_c2h << "\t"
          << number_density_c2h2 << "\t"
          << number_density[i] << "\t"
          << particle_radius[i][0]*1e4 << "\t"
          << nucleation_rate[i] << "\t"
          << growth_rate[i] << "\t"
          << dust_moments[0][i] << "\t"
          << dust_moments[1][i] << "\t"
          << dust_moments[2][i] << "\t"
          << dust_moments[3][i] << "\t"
          << dust_moments[3][i]/dust_moments[0][i] * 4./3. * constants::pi * std::pow(monomer_radius, 3) << "\t"
          << "\n";
  }

  file.close();
}



}