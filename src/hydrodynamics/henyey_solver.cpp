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

void Hydrodynamics::resetWarmStart()
{
  henyey_x_prev.clear();
  henyey_chi_prev.clear();
}

bool Hydrodynamics::solveStructureHenyey(
  const std::vector<double>& flux_weighted_extinction,
  const double seed_mass_loss_rate)
{
  const size_t M = nb_grid_points;
  const size_t s_idx = (1 + StructureSolver::nb_moments)*M;
  const size_t N = s_idx + 1;

  //reject if the (RT-derived) flux-mean extinction is not usable
  for (const double x : flux_weighted_extinction)
    if (!std::isfinite(x) || x < 0.)
    {
      std::cout << "Henyey solve skipped (non-finite flux-mean extinction)"
                   " - keeping previous structure\n";
      return false;
    }

  //Damp the change of the flux-mean extinction (hence the radiative acceleration
  //alpha) between outer iterations, in log space. alpha sets where the wind-equation
  //numerator N=0, i.e. the sonic radius r_c; an unrelaxed alpha makes r_c jump
  //dramatically between iterations. Relaxing alpha lets r_c migrate smoothly while
  //the critical point stays FREE (relocated each Newton step to where Phi is
  //minimal, consistent with dPhi/dr=0). This is Winters' "damping factors applied
  //to the radiative-acceleration corrections".
  const double w_chi = 0.1;
  std::vector<double> chi_used(M);
  for (size_t i=0; i<M; ++i)
  {
    const double chi_new = flux_weighted_extinction[i];
    chi_used[i] = (henyey_chi_prev.size()==M && henyey_chi_prev[i] > 0. && chi_new > 0.)
      ? std::pow(henyey_chi_prev[i], 1.-w_chi) * std::pow(chi_new, w_chi)
      : chi_new;
  }
  henyey_chi_prev = chi_used;

  //couple the dust moments in the Newton with the (frozen) Gail-Sedlmayr kernels.
  //GATED (HENYEY_DUST): the K0..K5 block spans ~20 orders of magnitude, so the
  //coupled Jacobian is ill-conditioned and the un-scaled Newton fails -- this needs
  //per-block residual/variable scaling first. Off by default (dust decoupled, the
  //opacity enters via the frozen chi_F) so the working Setup A is preserved.
  const bool have_dust_kernels = std::getenv("HENYEY_DUST")
    && dust_nucleation_rate.size() == M && dust_growth_timescale.size() == M;

  StructureSolver solver(M);
  for (size_t i=0; i<M; ++i)
  {
    solver.radius[i] = atmosphere->radius[i];
    solver.sound_speed[i] = isothermal_sound_speed[i];
    solver.flux_mean_extinction[i] = chi_used[i];

    if (have_dust_kernels && atmosphere->total_h_density[i] > 0.)
    {
      solver.nucleation_per_h[i] = dust_nucleation_rate[i] / atmosphere->total_h_density[i];
      solver.growth_timescale[i] = std::max(dust_growth_timescale[i], 1e-30);
    }
    else
    {
      solver.nucleation_per_h[i] = 0.;
      solver.growth_timescale[i] = 1.;
    }
  }
  solver.sound_speed2_deriv = sound_speed_derivative;
  solver.stellar_luminosity = config->stellar_luminosity;
  solver.gravitational_mass = constants::gravitation_const * config->stellar_mass;

  //Stage 2 dust radiative-acceleration feedback (with the dust coupled): split the
  //(damped) flux-mean extinction into gas + dust parts by their current flux-mean
  //fractions, and let the dust part scale with the in-Newton 2nd moment relative to
  //its reference value. This puts the local dust -> alpha -> wind loop in the Newton.
  //GATED (HENYEY_ALPHA_FB), and OFF by default: for a marginal dust-driven wind
  //(alpha ~ 1 at the sonic point) this feedback is near-critical -- any K2 excursion
  //pushes alpha super-Eddington and the Newton blows up (cp->1, Mdot->NaN). The
  //robust approach is to keep alpha frozen in the Newton and damp it in the outer
  //loop (Setup A), as Winters/Dominik do. Enabling this needs feedback-strength
  //continuation (ramp) or a trust-region step -- future work.
  if (have_dust_kernels && std::getenv("HENYEY_ALPHA_FB"))
  {
    solver.dust_alpha_feedback = true;
    solver.chi_gas.assign(M, 0.);
    solver.chi_dust_ref.assign(M, 0.);
    solver.k2_ref.assign(M, 1e-300);

    for (size_t i=0; i<M; ++i)
    {
      const double chi_total = flux_weighted_extinction[i];
      const double chi_gas_raw =
        radiation_field[i].fluxWeightedExtinction(atmosphere->extinction_coeff_gas[i]);
      const double gas_frac = (chi_total > 0.)
        ? std::min(std::max(chi_gas_raw/chi_total, 0.), 1.) : 1.;

      solver.chi_gas[i] = chi_used[i] * gas_frac;
      solver.chi_dust_ref[i] = chi_used[i] * (1. - gas_frac);

      if (dust_second_moment.size() == M && atmosphere->total_h_density[i] > 0.)
        solver.k2_ref[i] =
          std::max(dust_second_moment[i] / atmosphere->total_h_density[i], 1e-300);
    }
  }

  //inner-velocity boundary condition fixed at the PHYSICAL base velocity (from the
  //current structure), not an arbitrary fraction of c_T. A too-high inner velocity
  //gives a too-low base density (rho = Mdot/4pi r^2 v) and starves the dust driving.
  solver.inner_velocity = std::max(atmosphere->velocity[0], 1.0);

  //Setup B (Winters, HENYEY_FIX_MDOT=1): prescribe Mdot at the config value rather
  //than solving for it as the (collapse-prone) transonic eigenvalue. The inner BC
  //and Mdot are both fixed; the regularity is not imposed (the structure is the
  //forward solution of the wind equation for the prescribed Mdot).
  double seed_mdot = seed_mass_loss_rate;
  if (std::getenv("HENYEY_FIX_MDOT"))
  {
    solver.fix_mass_loss_rate = true;
    solver.log_mass_loss_rate_target = std::log(config->stellar_mass_loss_rate);
    seed_mdot = config->stellar_mass_loss_rate;
  }

  //seed: warm-start from the previous solve, else a transonic guess on c_T
  std::vector<double> x;
  if (henyey_x_prev.size() == N)
    x = henyey_x_prev;
  else
  {
    x.assign(N, 0.);
    for (size_t i=0; i<M; ++i)
    {
      const double frac = 0.3 + 1.5*double(i)/double(M);
      const double v = frac * isothermal_sound_speed[i];
      x[i] = 0.5*(v + isothermal_sound_speed[i]*isothermal_sound_speed[i]/v);
    }
    x[s_idx] = std::log(seed_mdot);
  }

  //Setup A: Mdot as the transonic eigenvalue, critical point FREE (relocated each
  //Newton step). The chi_F damping above keeps r_c migrating smoothly.
  //Setup B (fix_mass_loss_rate): Mdot pinned, regularity not imposed.
  const int max_newton = 300;
  double res = 0.; int iters = 0;
  x = solver.solveEigen(x, max_newton, 1e-8, &res, &iters);

  const double mdot_cgs = std::exp(x[s_idx]);
  const int cp = solver.critical_point;

  //reject a non-converged or degenerate solve (iteration cap, boundary critical
  //point, or non-finite Mdot): keep the last good structure and do not warm-start
  //from this x, so a diverged solve cannot poison the next iteration.
  if (iters >= max_newton || !std::isfinite(mdot_cgs) || mdot_cgs <= 0.
      || cp <= 0 || cp >= static_cast<int>(M)-1)
  {
    std::cout << "Henyey solve rejected (iters=" << iters << ", cp=" << cp
              << ", Mdot=" << mdot_cgs << ") - keeping previous structure\n";
    return false;
  }

  const std::vector<double> v = solver.velocityProfile(x);

  //Under-relax the structure write-back (log space). The mass-loss eigenvalue is
  //a steep function of the structure (dust-driving feedback), so an unrelaxed
  //Mdot collapses/explodes across outer iterations -- this is exactly why
  //Winters/Dominik fix Mdot. Damping the coupling keeps Setup A (predicted Mdot)
  //stable. The per-solve Newton is unaffected (it solves at frozen chi_F).
  const double w = 0.3;

  const double mdot_old = mass_loss_rate * constants::mass_sun / constants::year;  //CGS, prev step
  const double mdot_rel = (mdot_old > 0.)
    ? std::pow(mdot_old, 1.-w) * std::pow(mdot_cgs, w) : mdot_cgs;

  for (size_t i=0; i<M; ++i)
  {
    const double v_old = atmosphere->velocity[i];
    const double v_rel = (v_old > 0. && std::isfinite(v[i]) && v[i] > 0.)
      ? std::pow(v_old, 1.-w) * std::pow(v[i], w) : v[i];

    atmosphere->velocity[i] = v_rel;
    atmosphere->mass_density[i] = mdot_rel
      / (4.*constants::pi * atmosphere->radius[i]*atmosphere->radius[i] * v_rel);
  }

  mass_loss_rate = mdot_rel / constants::mass_sun * constants::year;
  henyey_x_prev = x;

  //raw eigenvalue Mdot (what the solve wants) vs the relaxed value actually used;
  //the raw/relaxed ratio approaching 1 signals a converged outer fixed point
  const double mdot_raw_msunyr = mdot_cgs / constants::mass_sun * constants::year;
  std::cout << "Henyey structure solve: rel-residual " << res << " in " << iters
            << " iters; critical point " << cp << "/" << M
            << "; Mdot raw " << mdot_raw_msunyr << " relaxed " << mass_loss_rate
            << " Msun/yr (raw/relaxed " << mdot_raw_msunyr/mass_loss_rate << ")\n";

  return true;
}


}