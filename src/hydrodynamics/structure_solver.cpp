
#include "structure_solver.h"

#include <cppad/cppad.hpp>
#include <Eigen/Dense>

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>


namespace agb {

//physical constants [cgs] (match src/additional/physical_const.h)
namespace {
  constexpr double light_c = 2.99792458e10;
  constexpr double sigma_sb = 5.670374e-5;     //Stefan-Boltzmann
  constexpr double boltzmann_k = 1.380649e-16;
  constexpr double mass_proton = 1.6726219e-24;
  constexpr double pi = 3.14159265358979;
}

//floor a value at a small positive number, in a way that is safe to tape with
//CppAD (CondExp keeps the operation differentiable through the branch)
namespace {
  inline double safePositive(const double x, const double floor)
  {
    return x > floor ? x : floor;
  }
  inline CppAD::AD<double> safePositive(const CppAD::AD<double>& x, const double floor)
  {
    return CppAD::CondExpGt(x, CppAD::AD<double>(floor), x, CppAD::AD<double>(floor));
  }

  //clamp to [lo, hi], taped-safe for both scalar types
  inline double clampRange(const double x, const double lo, const double hi)
  {
    return x < lo ? lo : (x > hi ? hi : x);
  }
  inline CppAD::AD<double> clampRange(const CppAD::AD<double>& x, const double lo, const double hi)
  {
    const CppAD::AD<double> a = CppAD::CondExpLt(x, CppAD::AD<double>(lo), CppAD::AD<double>(lo), x);
    return CppAD::CondExpGt(a, CppAD::AD<double>(hi), CppAD::AD<double>(hi), a);
  }

}


StructureSolver::StructureSolver(const size_t nb_points_)
  : nb_points(nb_points_)
{
  radius.assign(nb_points, 0.);
  sound_speed.assign(nb_points, 0.);
  flux_mean_extinction.assign(nb_points, 0.);
  sound_speed2_deriv.assign(nb_points, 0.);
  nucleation_per_h.assign(nb_points, 0.);
  growth_timescale.assign(nb_points, 1.);
}


template<class Scalar>
Scalar StructureSolver::windVelocity(
  const Scalar& phi, const double c_t, const bool supersonic) const
{
  using std::sqrt;  //CppAD::sqrt found by ADL for AD<double>

  //floor the discriminant just above zero so the velocity stays real and the
  //tape is NaN-free even if a Newton step momentarily pushes phi toward c_T
  const Scalar root = sqrt(safePositive(phi*phi - c_t*c_t, 1.0));

  return supersonic ? phi + root : phi - root;
}


template<class Scalar>
std::vector<Scalar> StructureSolver::residualPhi(
  const std::vector<Scalar>& phi) const
{
  std::vector<Scalar> residual(nb_points, Scalar(0.));

  const Scalar mdot(mass_loss_rate);

  //RHS of dPhi/dr = -G M /(2 v r^2) (1 - alpha) + c_T^2 /(v r),
  //with alpha = chi_F L r^2 v /(c G M Mdot)
  auto rhs = [&](const size_t i) -> Scalar {
    const bool supersonic = (static_cast<int>(i) > critical_point);
    const Scalar v = windVelocity(phi[i], sound_speed[i], supersonic);
    const double r = radius[i];
    const double c2 = sound_speed[i]*sound_speed[i];

    const Scalar alpha = flux_mean_extinction[i] * stellar_luminosity * r*r * v
                       / (light_c * gravitational_mass * mdot);

    return -gravitational_mass/(2.*v*r*r) * (1. - alpha) + c2/(v*r);
  };

  //inner boundary condition: Phi fixed by the prescribed inner velocity
  const double v_in = inner_velocity;
  const double phi_inner = 0.5*(v_in + sound_speed[0]*sound_speed[0]/v_in);
  residual[0] = phi[0] - phi_inner;

  //interior: trapezoidal discretization of the first-order ODE
  for (size_t i=1; i<nb_points; ++i)
  {
    const double h = radius[i] - radius[i-1];
    residual[i] = (phi[i] - phi[i-1]) - 0.5*h*(rhs(i) + rhs(i-1));
  }

  return residual;
}

template std::vector<double> StructureSolver::residualPhi<double>(
  const std::vector<double>&) const;


//fill the wind (Phi) + dust-moment residual blocks of a block-ordered unknown
//vector x = [Phi(M), K0(M), ..., K5(M), ...], for a given mass-loss rate mdot
//(passed as a Scalar so it may itself be an unknown of the eigenvalue problem).
template<class Scalar>
void StructureSolver::fillStructureResidual(
  const std::vector<Scalar>& x, const Scalar& mdot,
  std::vector<Scalar>& residual) const
{
  const size_t M = nb_points;

  auto phi = [&](const size_t i) -> const Scalar& { return x[i]; };
  auto K   = [&](const int j, const size_t i) -> const Scalar& {
    return x[(1 + j)*M + i]; };

  auto vel = [&](const size_t i) -> Scalar {
    const bool supersonic = (static_cast<int>(i) > critical_point);
    return windVelocity(phi(i), sound_speed[i], supersonic); };

  //flux-mean extinction at point i. With the dust feedback enabled it depends on
  //the in-Newton 2nd moment K2 = K(2,i): chi_F = chi_gas + chi_dust_ref*(K2/k2_ref).
  //dust-opacity ratio bounded to [0, 5]: a transient K2 overshoot in a Newton step
  //must not blow alpha up to super-Eddington/NaN; the bound is inactive at
  //convergence (where the ratio ~ 1 since k2_ref tracks the converged K2).
  auto chi_mean = [&](const size_t i) -> Scalar {
    if (!dust_alpha_feedback) return Scalar(flux_mean_extinction[i]);
    const Scalar ratio = clampRange(K(2, i) / std::max(k2_ref[i], 1e-300), 0., 5.);
    return chi_gas[i] + chi_dust_ref[i] * ratio; };

  //--- wind (Melia-Phi) residual, block 0 ---
  auto phi_rhs = [&](const size_t i) -> Scalar {
    const Scalar v = vel(i);
    const double r = radius[i];
    const double c2 = sound_speed[i]*sound_speed[i];
    const Scalar alpha = chi_mean(i) * stellar_luminosity * r*r * v
                       / (light_c * gravitational_mass * mdot);
    return -gravitational_mass/(2.*v*r*r) * (1. - alpha) + c2/(v*r); };

  const double v_in = inner_velocity;
  const double phi_inner = 0.5*(v_in + sound_speed[0]*sound_speed[0]/v_in);
  residual[0] = phi(0) - phi_inner;

  for (size_t i=1; i<M; ++i)
  {
    const double h = radius[i] - radius[i-1];
    residual[i] = (phi(i) - phi(i-1)) - 0.5*h*(phi_rhs(i) + phi_rhs(i-1));
  }

  //--- dust-moment residuals, blocks 1..nb_moments ---
  //dK0/dr = J*/v ;  dKj/dr = (1/v)[ (j/3)/tau * K_{j-1} + N_l^{j/3} * J* ]
  auto moment_rhs = [&](const int j, const size_t i) -> Scalar {
    const Scalar v = vel(i);
    const double jstar = nucleation_per_h[i];

    if (j == 0)
      return jstar / v;

    const double growth_coeff = (j/3.) / growth_timescale[i];
    const double nucleation_coeff = std::pow(min_monomer_number, j/3.);
    return (growth_coeff * K(j-1, i) + nucleation_coeff * jstar) / v; };

  for (int j=0; j<nb_moments; ++j)
  {
    const size_t base = (1 + j)*M;
    residual[base] = K(j, 0);  //inner boundary: no dust at the base

    for (size_t i=1; i<M; ++i)
    {
      const double h = radius[i] - radius[i-1];
      residual[base + i] = (K(j, i) - K(j, i-1))
                         - 0.5*h*(moment_rhs(j, i) + moment_rhs(j, i-1));
    }
  }
}


template<class Scalar>
std::vector<Scalar> StructureSolver::residualCoupled(
  const std::vector<Scalar>& x) const
{
  std::vector<Scalar> residual(x.size(), Scalar(0.));
  fillStructureResidual(x, Scalar(mass_loss_rate), residual);
  return residual;
}

template std::vector<double> StructureSolver::residualCoupled<double>(
  const std::vector<double>&) const;


//coupled system with Mdot promoted to a transonic eigenvalue. The unknown vector
//is x = [ Phi(M), K0..K5(M), s ] with s = ln(Mdot) (s last). The extra residual
//is the critical-point regularity: the wind-equation numerator vanishes at the
//sonic point r_c, which fixes Mdot.
template<class Scalar>
std::vector<Scalar> StructureSolver::residualEigen(
  const std::vector<Scalar>& x) const
{
  using std::exp;

  const size_t M = nb_points;
  const size_t base = (1 + nb_moments) * M;  //index of s = ln(Mdot)

  const Scalar mdot = exp(x[base]);

  std::vector<Scalar> residual(x.size(), Scalar(0.));
  fillStructureResidual(x, mdot, residual);  //sets residual[0] = inner-velocity BC

  //--- critical-point regularity ---
  //N = -G M / r_cp^2 (1 - alpha_cp) + 2 c_T^2 / r_cp - d(c_T^2)/dr |_cp = 0
  const size_t cp = static_cast<size_t>(critical_point);
  const bool supersonic = false;  //at r_c the branch is degenerate (v = c_T)
  const Scalar v_cp = windVelocity(x[cp], sound_speed[cp], supersonic);
  const double r = radius[cp];
  const double c2 = sound_speed[cp]*sound_speed[cp];

  //flux-mean extinction at the critical point, with the dust feedback (K2-dependent)
  Scalar chi_cp = Scalar(flux_mean_extinction[cp]);
  if (dust_alpha_feedback)
  {
    const Scalar ratio = clampRange(x[(1+2)*M + cp] / std::max(k2_ref[cp], 1e-300), 0., 5.);
    chi_cp = chi_gas[cp] + chi_dust_ref[cp] * ratio;
  }

  const Scalar alpha_cp = chi_cp * stellar_luminosity * r*r * v_cp
                        / (light_c * gravitational_mass * mdot);

  const Scalar regularity = -gravitational_mass/(r*r) * (1. - alpha_cp)
                          + 2.*c2/r - sound_speed2_deriv[cp];

  if (fix_mass_loss_rate)
  {
    //Setup B (Winters): Mdot prescribed. Keep the inner-velocity BC (already set in
    //residual[0] by fillStructureResidual) and pin Mdot; the regularity is NOT
    //imposed -- the structure is the forward solution of the wind equation for the
    //prescribed Mdot and inner boundary condition.
    (void) regularity;
    residual[base] = x[base] - log_mass_loss_rate_target;
  }
  else
  {
    //Setup A: Mdot is the eigenvalue, closed by the regularity; inner BC kept.
    //(A sub-grid sonic-radius evaluation was tried but proved fragile in the cold
    //phase -- r_c slammed to clamp edges and flipped discontinuously, diverging the
    //Newton -- so the node-based regularity is retained; see project memory.)
    residual[base] = regularity;
  }

  return residual;
}

template std::vector<double> StructureSolver::residualEigen<double>(
  const std::vector<double>&) const;


std::vector<double> StructureSolver::solvePhi(
  std::vector<double> phi,
  const int max_iterations,
  const double tolerance,
  double* final_residual,
  int* iterations)
{
  using CppAD::AD;

  //tape the residual once on the structure of the problem; re-zero/re-tape each
  //Newton step since the branch (critical_point) and coefficients are fixed here
  const size_t n = nb_points;

  double res_norm = 0.;
  int it = 0;

  for (; it < max_iterations; ++it)
  {
    //--- residual (double) ---
    std::vector<double> res = residualPhi(phi);

    res_norm = 0.;
    for (const double r : res) res_norm += r*r;
    res_norm = std::sqrt(res_norm);

    if (res_norm < tolerance) break;

    //--- CppAD Jacobian ---
    std::vector<AD<double>> aphi(n);
    for (size_t i=0; i<n; ++i) aphi[i] = phi[i];

    CppAD::Independent(aphi);
    std::vector<AD<double>> ares = residualPhi(aphi);
    CppAD::ADFun<double> f(aphi, ares);

    const std::vector<double> jac = f.Jacobian(phi);  //dense, row-major n*n

    //--- linear solve  J * delta = -res  (Eigen dense) ---
    Eigen::MatrixXd J(n, n);
    Eigen::VectorXd b(n);
    for (size_t i=0; i<n; ++i)
    {
      b(i) = -res[i];
      for (size_t k=0; k<n; ++k) J(i, k) = jac[i*n + k];
    }

    const Eigen::VectorXd delta = J.partialPivLu().solve(b);

    //--- damped update ---
    //limit the global step so (a) no component moves more than 50% and (b) no
    //Phi is driven below its local sound speed (which would make v complex)
    double damping = 1.0;
    for (size_t i=0; i<n; ++i)
    {
      const double max_rel = 0.5 * std::abs(phi[i]);
      if (std::abs(delta(i)) > max_rel && max_rel > 0.)
        damping = std::min(damping, max_rel/std::abs(delta(i)));

      //keep phi[i] + damping*delta >= 1.01*c_T[i]
      const double floor = 1.01 * sound_speed[i];
      if (delta(i) < 0.)
      {
        const double allowed = (floor - phi[i]) / delta(i);  //>0 since delta<0
        if (allowed > 0. && allowed < damping) damping = allowed;
      }
    }

    for (size_t i=0; i<n; ++i)
      phi[i] += damping * delta(i);
  }

  if (final_residual) *final_residual = res_norm;
  if (iterations) *iterations = it;

  return phi;
}


namespace {
  //relocate the critical point to the interior node where Phi is closest to c_T
  int locateCriticalPoint(
    const std::vector<double>& x, const std::vector<double>& c_t, const size_t M)
  {
    int best = 1;
    double best_gap = 1e300;
    for (size_t i=1; i<M-1; ++i)
    {
      const double gap = x[i] - c_t[i];
      if (gap < best_gap) { best_gap = gap; best = static_cast<int>(i); }
    }
    return best;
  }
}


double StructureSolver::residualEigenNorm(const std::vector<double>& x)
{
  critical_point = locateCriticalPoint(x, sound_speed, nb_points);

  const std::vector<double> res = residualEigen(x);

  double norm = 0.;
  for (const double r : res) norm += r*r;
  return std::sqrt(norm);
}


int StructureSolver::sonicNode(const std::vector<double>& x) const
{
  return locateCriticalPoint(x, sound_speed, nb_points);
}


std::vector<double> StructureSolver::velocityProfile(const std::vector<double>& x) const
{
  std::vector<double> v(nb_points, 0.);
  for (size_t i=0; i<nb_points; ++i)
  {
    const bool supersonic = (static_cast<int>(i) > critical_point);
    v[i] = windVelocity(x[i], sound_speed[i], supersonic);
  }
  return v;
}


std::vector<double> StructureSolver::solveEigen(
  std::vector<double> x,
  const int max_iterations,
  const double tolerance,
  double* final_residual,
  int* iterations)
{
  using CppAD::AD;

  const size_t M = nb_points;
  const size_t N = (1 + nb_moments)*M + 1;
  const size_t s_idx = (1 + nb_moments)*M;  //index of ln(Mdot)

  double res_norm = 0.;
  int it = 0;

  for (; it < max_iterations; ++it)
  {
    if (relocate_critical_point)
      critical_point = locateCriticalPoint(x, sound_speed, M);

    std::vector<double> res = residualEigen(x);

    //--- per-block scaling ---
    //Scale each variable by its magnitude (Phi by its block peak ~ inner Phi, each
    //K_j by its block peak, ln(Mdot) absolute) and each equation by its primary
    //variable's scale (the regularity row by ~G M / r_c^2). Without this the
    //K0..K5 block spans ~20 orders of magnitude and the dense solve is hopelessly
    //ill-conditioned. Solve (R^-1 J S) y = -R^-1 res, then delta = S y.
    double s_phi = 1.0;
    for (size_t i=0; i<M; ++i) s_phi = std::max(s_phi, std::abs(x[i]));
    std::vector<double> s_k(nb_moments, 1e-300);
    for (int j=0; j<nb_moments; ++j)
      for (size_t i=0; i<M; ++i)
        s_k[j] = std::max(s_k[j], std::abs(x[(1+j)*M + i]));
    const double s_mdot = std::max(std::abs(x[s_idx]), 1.0);

    auto var_scale = [&](const size_t idx) -> double {
      if (idx == s_idx) return s_mdot;
      const size_t blk = idx / M;            //0 = Phi, 1..nb_moments = K0..K5
      return (blk == 0) ? s_phi : s_k[blk-1]; };

    auto eq_scale = [&](const size_t k) -> double {
      if (k == s_idx)
        return fix_mass_loss_rate ? 1.0
             : gravitational_mass/(radius[critical_point]*radius[critical_point]);
      return var_scale(k); };

    //dimensionless (scaled) residual RMS as the convergence measure
    res_norm = 0.;
    for (size_t k=0; k<N; ++k) { const double e = res[k]/eq_scale(k); res_norm += e*e; }
    res_norm = std::sqrt(res_norm / N);

    if (res_norm < tolerance) break;

    //--- CppAD dense Jacobian ---
    std::vector<AD<double>> ax(N);
    for (size_t i=0; i<N; ++i) ax[i] = x[i];

    CppAD::Independent(ax);
    std::vector<AD<double>> ares = residualEigen(ax);
    CppAD::ADFun<double> f(ax, ares);

    const std::vector<double> jac = f.Jacobian(x);

    //--- scaled linear solve ---
    Eigen::MatrixXd J(N, N);
    Eigen::VectorXd b(N);
    for (size_t k=0; k<N; ++k)
    {
      const double rk = eq_scale(k);
      b(k) = -res[k]/rk;
      for (size_t j=0; j<N; ++j) J(k, j) = jac[k*N + j] * var_scale(j) / rk;
    }

    const Eigen::VectorXd y = J.partialPivLu().solve(b);
    Eigen::VectorXd delta(N);
    for (size_t j=0; j<N; ++j) delta(j) = var_scale(j) * y(j);

    //--- damped update ---
    //The residual mixes scales over ~20 orders of magnitude, so a residual-norm
    //line search is a poor merit function (it stalls) and an uncapped Newton step
    //overshoots. A fixed fractional cap + the physical safety constraints (Phi
    //above c_T, bounded ln(Mdot) step) converges reliably (empirically ~100 iters).
    double damping = 1.0;
    for (size_t i=0; i<M; ++i)
    {
      const double max_rel = 0.5 * std::abs(x[i]);
      if (std::abs(delta(i)) > max_rel && max_rel > 0.)
        damping = std::min(damping, max_rel/std::abs(delta(i)));

      if (delta(i) < 0.)
      {
        const double allowed = (sound_speed[i] - x[i]) / delta(i);  //>0
        if (allowed > 0. && allowed < damping) damping = allowed;
      }
    }
    if (std::abs(delta(s_idx)) > 0.5)
      damping = std::min(damping, 0.5/std::abs(delta(s_idx)));

    //apply the step and measure its (scale-relative) size. A residual-relative
    //stop alone fails for warm starts (the initial residual is already small, so
    //an 8-order further reduction is unreachable); converging on a negligible
    //Newton step is scale-robust and handles both cold and warm seeds.
    double max_rel_step = 0.;
    for (size_t i=0; i<N; ++i)
    {
      x[i] += damping * delta(i);
      if (std::abs(x[i]) > 0.)
        max_rel_step = std::max(max_rel_step, std::abs(damping*delta(i))/std::abs(x[i]));
    }
    if (max_rel_step < 1e-12) { ++it; break; }
  }

  if (final_residual) *final_residual = res_norm;
  if (iterations) *iterations = it;

  return x;
}


std::vector<double> StructureSolver::greyBootstrap(
  const GreyParameters& params,
  std::vector<double>& T_out,
  const int max_grey_iterations,
  const double tolerance,
  int* grey_iterations)
{
  const size_t M = nb_points;
  const size_t s_idx = (1 + nb_moments)*M;
  const size_t N = s_idx + 1;

  const double T_eff = std::pow(
    stellar_luminosity / (4.*pi * params.stellar_radius*params.stellar_radius * sigma_sb),
    0.25);

  T_out.assign(M, T_eff);

  //--- initial guess: a transonic velocity profile + a starting Mdot, no dust ---
  std::vector<double> x(N, 0.);
  {
    //provisional sound speed from T_eff to shape the initial Phi
    for (size_t i=0; i<M; ++i)
    {
      const double c0 = std::sqrt(boltzmann_k * T_eff / (params.mean_molecular_weight*mass_proton));
      const double frac = 0.3 + 1.5*double(i)/double(M);
      const double v = frac * c0;
      x[i] = 0.5*(v + c0*c0/v);
    }
    x[s_idx] = std::log(1.0e21);
  }

  double mdot_prev = std::exp(x[s_idx]);
  int it = 0;

  for (; it < max_grey_iterations; ++it)
  {
    const double mdot = std::exp(x[s_idx]);

    //reconstruct velocity (branch from current critical point) and density
    critical_point = locateCriticalPoint(x, sound_speed.empty() ? std::vector<double>(M,1.) : sound_speed, M);

    std::vector<double> rho(M, 0.), chi_grey(M, 0.);
    for (size_t i=0; i<M; ++i)
    {
      const bool supersonic = (static_cast<int>(i) > critical_point);
      //use the previous sound speed if available, else a nominal value
      const double c = sound_speed.empty() ? 1.0e5 : sound_speed[i];
      const double disc = std::max(x[i]*x[i] - c*c, 1.0);
      const double v = supersonic ? x[i] + std::sqrt(disc) : x[i] - std::sqrt(disc);

      rho[i] = mdot / (4.*pi * radius[i]*radius[i] * std::max(v, 1.0));

      const double n_h = rho[i] / (params.mean_molecular_weight * mass_proton);
      const double k2_norm = x[(1+2)*M + i];  //normalized 2nd dust moment
      const double chi_dust = pi * params.monomer_radius*params.monomer_radius
                            * params.dust_q_ext * std::max(k2_norm, 0.) * n_h;
      chi_grey[i] = params.gas_opacity * rho[i] + chi_dust;
    }

    //optical depth integrated inward from the outer edge
    std::vector<double> tau(M, 0.);
    for (int i=M-2; i>=0; --i)
      tau[i] = tau[i+1] + 0.5*(chi_grey[i]+chi_grey[i+1])*(radius[i+1]-radius[i]);

    //Lucy spherical-grey temperature, and the resulting sound speed
    sound_speed.assign(M, 0.);
    for (size_t i=0; i<M; ++i)
    {
      const double ratio = params.stellar_radius/radius[i];
      const double w = 0.5*(1.0 - std::sqrt(std::max(1.0 - ratio*ratio, 0.0)));
      const double t4 = T_eff*T_eff*T_eff*T_eff * (w + 0.75*tau[i]);
      T_out[i] = std::pow(t4, 0.25);
      sound_speed[i] = std::sqrt(boltzmann_k * T_out[i] / (params.mean_molecular_weight*mass_proton));
    }

    //d(c_T^2)/dr
    for (size_t i=0; i<M; ++i)
    {
      const size_t lo = (i==0)?0:i-1, hi = (i==M-1)?M-1:i+1;
      sound_speed2_deriv[i] = (sound_speed[hi]*sound_speed[hi] - sound_speed[lo]*sound_speed[lo])
                            / (radius[hi]-radius[lo]);
    }

    //grey dust source terms: nucleation switches on smoothly below T_cond
    for (size_t i=0; i<M; ++i)
    {
      const double s = 1.0/(1.0 + std::exp((T_out[i]-params.condensation_temperature)/50.0));
      nucleation_per_h[i] = params.nucleation_per_h_scale * s + 1.0e-30;
      growth_timescale[i] = params.dust_growth_timescale;
    }

    flux_mean_extinction = chi_grey;
    inner_velocity = 0.3 * sound_speed[0];

    //solve the wind+dust structure for this grey state
    double res = 0.; int newton_it = 0;
    x = solveEigen(x, 300, 1e-8, &res, &newton_it);

    const double mdot_new = std::exp(x[s_idx]);
    const double rel_change = std::abs(mdot_new - mdot_prev)/mdot_prev;
    mdot_prev = mdot_new;

    if (rel_change < tolerance) { ++it; break; }
  }

  if (grey_iterations) *grey_iterations = it;
  return x;
}


double StructureSolver::selfTestGreyBootstrap(int* cp_out, bool* monotone_T)
{
  const size_t M = nb_points;

  //IRC10216-scale setup
  const double r_inner = 9.0e13;
  stellar_luminosity = 3.0e4 * 3.828e33;
  gravitational_mass = 6.674e-8 * 0.7 * 1.98847e33;

  for (size_t i=0; i<M; ++i)
    radius[i] = r_inner * (1.0 + 0.06*i);

  GreyParameters params;
  params.stellar_radius = r_inner;
  params.gas_opacity = 1.0e-4;
  params.nucleation_per_h_scale = 1.0e-18;
  params.dust_q_ext = 2.0;
  params.dust_growth_timescale = 1.0e5;
  params.condensation_temperature = 1500.;

  std::vector<double> T;
  int grey_it = 0;
  const std::vector<double> x = greyBootstrap(params, T, 60, 1e-4, &grey_it);

  const size_t s_idx = (1 + nb_moments)*M;
  const double mdot = std::exp(x[s_idx]);

  bool monotone = true;
  for (size_t i=1; i<M; ++i)
    if (T[i] > T[i-1] + 1.0) monotone = false;  //decreasing outward (1 K slack)

  if (cp_out) *cp_out = critical_point;
  if (monotone_T) *monotone_T = monotone;

  return mdot;
}


void StructureSolver::setSyntheticState()
{
  const double r_inner = 9.0e13;
  stellar_luminosity = 3.0e4 * 3.828e33;            //3e4 Lsun
  gravitational_mass = 6.674e-8 * 0.7 * 1.98847e33; //G * 0.7 Msun
  mass_loss_rate = 5.0e21;                          //~1e-4 Msun/yr in g/s
  inner_velocity = 1.2e5;                           //sub-sonic (c_T base ~2e5)
  //all-sub-sonic for the infrastructure self-test (no branch switch): the
  //transonic critical point is introduced with the eigenvalue closure
  critical_point = static_cast<int>(nb_points) + 1;

  nucleation_per_h.assign(nb_points, 0.);
  growth_timescale.assign(nb_points, 0.);

  for (size_t i=0; i<nb_points; ++i)
  {
    //gentle radial range + near radiative/gravity force balance (alpha ~ 0.8 at
    //the base) so a smooth sub-sonic solution exists for this infrastructure test
    radius[i] = r_inner * (1.0 + 0.002*i);
    sound_speed[i] = 2.0e5 * std::pow(radius[0]/radius[i], 0.25);
    flux_mean_extinction[i] = 1.0e-13;

    //synthetic dust source coefficients (a smooth nucleation bump + a timescale)
    const double d = (double(i) - 0.5*nb_points) / (0.2*nb_points);
    nucleation_per_h[i] = 1.0e-20 * std::exp(-d*d);
    growth_timescale[i] = 1.0e7;
  }

  //d(c_T^2)/dr by central differences (frozen)
  sound_speed2_deriv.assign(nb_points, 0.);
  for (size_t i=0; i<nb_points; ++i)
  {
    const size_t lo = (i==0) ? 0 : i-1;
    const size_t hi = (i==nb_points-1) ? nb_points-1 : i+1;
    const double c2_hi = sound_speed[hi]*sound_speed[hi];
    const double c2_lo = sound_speed[lo]*sound_speed[lo];
    sound_speed2_deriv[i] = (c2_hi - c2_lo) / (radius[hi] - radius[lo]);
  }
}


double StructureSolver::selfTestPhiJacobian()
{
  using CppAD::AD;

  setSyntheticState();

  std::vector<double> phi(nb_points, 0.);
  for (size_t i=0; i<nb_points; ++i)
    phi[i] = 1.5*sound_speed[i] + 1.0e4*i;

  std::vector<AD<double>> aphi(nb_points);
  for (size_t i=0; i<nb_points; ++i) aphi[i] = phi[i];

  CppAD::Independent(aphi);
  std::vector<AD<double>> ares = residualPhi(aphi);
  CppAD::ADFun<double> f(aphi, ares);

  const std::vector<double> jac_ad = f.Jacobian(phi);

  double max_rel_error = 0.;
  for (size_t k=0; k<nb_points; ++k)
  {
    const double eps = 1e-6 * std::max(1.0, std::abs(phi[k]));

    std::vector<double> phi_p = phi, phi_m = phi;
    phi_p[k] += eps;
    phi_m[k] -= eps;

    const std::vector<double> res_p = residualPhi(phi_p);
    const std::vector<double> res_m = residualPhi(phi_m);

    for (size_t i=0; i<nb_points; ++i)
    {
      const double fd = (res_p[i] - res_m[i]) / (2.*eps);
      const double ad = jac_ad[i*nb_points + k];
      const double scale = std::max(1e-30, std::abs(ad) + std::abs(fd));
      max_rel_error = std::max(max_rel_error, std::abs(ad - fd)/scale);
    }
  }

  return max_rel_error;
}


double StructureSolver::selfTestCoupledJacobian()
{
  using CppAD::AD;

  setSyntheticState();

  const size_t M = nb_points;
  const size_t N = (1 + nb_moments) * M;

  //build a positive, smooth state vector [Phi, K0..K5]
  const double v_in = inner_velocity;
  const double phi_inner = 0.5*(v_in + sound_speed[0]*sound_speed[0]/v_in);

  std::vector<double> x(N, 0.);
  for (size_t i=0; i<M; ++i)
  {
    x[i] = phi_inner;
    for (int j=0; j<nb_moments; ++j)
      x[(1+j)*M + i] = 1.0e-13 * (1.0 + 0.01*i);
  }

  std::vector<AD<double>> ax(N);
  for (size_t i=0; i<N; ++i) ax[i] = x[i];

  CppAD::Independent(ax);
  std::vector<AD<double>> ares = residualCoupled(ax);
  CppAD::ADFun<double> f(ax, ares);

  const std::vector<double> jac_ad = f.Jacobian(x);  //dense N*N

  double max_rel_error = 0.;
  for (size_t k=0; k<N; ++k)
  {
    const double eps = 1e-6 * std::max(1e-20, std::abs(x[k]));

    std::vector<double> xp = x, xm = x;
    xp[k] += eps;
    xm[k] -= eps;

    const std::vector<double> rp = residualCoupled(xp);
    const std::vector<double> rm = residualCoupled(xm);

    for (size_t i=0; i<N; ++i)
    {
      const double fd = (rp[i] - rm[i]) / (2.*eps);
      const double ad = jac_ad[i*N + k];
      const double scale = std::max(1e-30, std::abs(ad) + std::abs(fd));
      max_rel_error = std::max(max_rel_error, std::abs(ad - fd)/scale);
    }
  }

  return max_rel_error;
}


double StructureSolver::selfTestEigenSolve(double* mdot_out, int* cp_out, int* iters_out)
{
  const size_t M = nb_points;
  const size_t s_idx = (1 + nb_moments)*M;
  const size_t N = s_idx + 1;

  //--- synthetic transonic state: dust-driven breeze crossing the sonic point ---
  const double r_inner = 9.0e13;
  stellar_luminosity = 3.0e4 * 3.828e33;
  gravitational_mass = 6.674e-8 * 0.7 * 1.98847e33;
  inner_velocity = 0.;  //unused by residualEigen except via the inner BC below
  min_monomer_number = 1000.;

  nucleation_per_h.assign(M, 0.);   //decouple dust for this wind+eigenvalue test
  growth_timescale.assign(M, 1.);

  for (size_t i=0; i<M; ++i)
  {
    radius[i] = r_inner * (1.0 + 0.06*i);
    sound_speed[i] = 2.0e5 * std::pow(radius[0]/radius[i], 0.25);
    flux_mean_extinction[i] = 1.0e-13;  //constant; alpha grows outward via r^2 v
  }

  sound_speed2_deriv.assign(M, 0.);
  for (size_t i=0; i<M; ++i)
  {
    const size_t lo = (i==0) ? 0 : i-1;
    const size_t hi = (i==M-1) ? M-1 : i+1;
    sound_speed2_deriv[i] = (sound_speed[hi]*sound_speed[hi] - sound_speed[lo]*sound_speed[lo])
                          / (radius[hi] - radius[lo]);
  }

  //transonic velocity guess: sub-sonic at base, crossing c_T near mid-grid
  std::vector<double> x(N, 0.);
  for (size_t i=0; i<M; ++i)
  {
    const double frac = 0.3 + 1.5*double(i)/double(M);   //0.3 -> ~1.8 of c_T
    const double v = frac * sound_speed[i];
    x[i] = 0.5*(v + sound_speed[i]*sound_speed[i]/v);
  }
  //inner velocity BC: anchor Phi[0] at the guessed inner velocity
  inner_velocity = 0.3 * sound_speed[0];
  x[s_idx] = std::log(5.0e21);  //Mdot guess [g/s]

  double final_res = 0.; int iters = 0;
  const std::vector<double> sol = solveEigen(x, 200, 1e-8, &final_res, &iters);

  const double seed_res = residualEigenNorm(x);
  const double mdot = std::exp(sol[s_idx]);

  if (mdot_out) *mdot_out = mdot;
  if (cp_out) *cp_out = critical_point;
  if (iters_out) *iters_out = iters;

  //relative residual reduction achieved
  return (seed_res > 0.) ? final_res/seed_res : final_res;
}


double StructureSolver::selfTestEigenJacobian()
{
  using CppAD::AD;

  setSyntheticState();
  critical_point = static_cast<int>(nb_points)/2;  //interior cp for the regularity row

  const size_t M = nb_points;
  const size_t N = (1 + nb_moments)*M + 1;  //+1 for s = ln(Mdot)

  const double v_in = inner_velocity;
  const double phi_inner = 0.5*(v_in + sound_speed[0]*sound_speed[0]/v_in);

  std::vector<double> x(N, 0.);
  for (size_t i=0; i<M; ++i)
  {
    x[i] = phi_inner;
    for (int j=0; j<nb_moments; ++j)
      x[(1+j)*M + i] = 1.0e-13 * (1.0 + 0.01*i);
  }
  x[(1+nb_moments)*M] = std::log(mass_loss_rate);  //s = ln(Mdot)

  std::vector<AD<double>> ax(N);
  for (size_t i=0; i<N; ++i) ax[i] = x[i];

  CppAD::Independent(ax);
  std::vector<AD<double>> ares = residualEigen(ax);
  CppAD::ADFun<double> f(ax, ares);

  const std::vector<double> jac_ad = f.Jacobian(x);

  double max_rel_error = 0.;
  for (size_t k=0; k<N; ++k)
  {
    const double eps = 1e-6 * std::max(1e-20, std::abs(x[k]));

    std::vector<double> xp = x, xm = x;
    xp[k] += eps;
    xm[k] -= eps;

    const std::vector<double> rp = residualEigen(xp);
    const std::vector<double> rm = residualEigen(xm);

    for (size_t i=0; i<N; ++i)
    {
      const double fd = (rp[i] - rm[i]) / (2.*eps);
      const double ad = jac_ad[i*N + k];
      const double scale = std::max(1e-30, std::abs(ad) + std::abs(fd));
      max_rel_error = std::max(max_rel_error, std::abs(ad - fd)/scale);
    }
  }

  return max_rel_error;
}


double StructureSolver::selfTestPhiSolve(int* iterations)
{
  setSyntheticState();

  //initial guess: the inner-boundary Phi held constant across the grid
  const double v_in = inner_velocity;
  const double phi_inner = 0.5*(v_in + sound_speed[0]*sound_speed[0]/v_in);

  std::vector<double> phi(nb_points, phi_inner);

  double final_residual = 0.;
  phi = solvePhi(phi, 100, 1e-10, &final_residual, iterations);

  return final_residual;
}


}
