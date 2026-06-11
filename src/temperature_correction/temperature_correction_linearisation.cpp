
#include "temperature_correction.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <omp.h>

#include "../radiative_transfer/radiative_transfer.h"
#include "../spectral_grid/spectral_grid.h"
#include "../config/config.h"
#include "../additional/physical_const.h"
#include "../additional/tri_diagonal_matrix.h"
#include "../additional/aux_functions.h"


namespace agb{


namespace {

//Dense Gaussian elimination with partial pivoting. Solves A x = b in place; A is
//row-major (n*n). Returns false if the matrix is singular.
bool solveDense(std::vector<double>& A, std::vector<double>& b, const size_t n)
{
  for (size_t k=0; k<n; ++k)
  {
    size_t piv = k;
    double pmax = std::abs(A[k*n+k]);
    for (size_t i=k+1; i<n; ++i)
    {
      const double v = std::abs(A[i*n+k]);
      if (v > pmax) { pmax = v; piv = i; }
    }
    if (pmax == 0.) return false;

    if (piv != k)
    {
      for (size_t j=k; j<n; ++j) std::swap(A[k*n+j], A[piv*n+j]);
      std::swap(b[k], b[piv]);
    }

    const double akk = A[k*n+k];
    for (size_t i=k+1; i<n; ++i)
    {
      const double f = A[i*n+k]/akk;
      if (f == 0.) continue;
      for (size_t j=k; j<n; ++j) A[i*n+j] -= f*A[k*n+j];
      b[i] -= f*b[k];
    }
  }

  for (int i=static_cast<int>(n)-1; i>=0; --i)
  {
    double s = b[i];
    for (size_t j=i+1; j<n; ++j) s -= A[i*n+j]*b[j];
    b[i] = s/A[i*n+i];
  }

  return true;
}

}



//Full-linearisation temperature correction (thesis 3.2.3 / App. B.1): one Newton step on
//(T_gas, T_dust) with the converged RT operator frozen. RadiativeTransfer supplies, per
//frequency, the moment operator V_nu, the diagonal thermal source dS_{s,i} = d rhs_i/d
//T_{s,i} and the residual K_nu (buildLinearisedMomentSystem). Here we own the radiative-
//equilibrium side: the per-frequency mean intensities are eliminated Rybicki-style (dJ_nu =
//V_nu^{-1}[dS_gas.dT_gas + dS_dust.dT_dust + K_nu]) by contracting against the local-RE
//residuals in ratio form (eq. 3.40/3.41, depth-stable) plus the composite flux-constancy
//term (eq. 3.64) on the gas equation, giving a dense (2D x 2D) Newton system in
//(dT_gas, dT_dust). This replaces the approximate Unsoeld-Lucy dJ/dT with the exact one.
void TemperatureCorrection::linearisedCorrection(
  RadiativeTransfer& rt,
  const std::vector<double>& temperature_gas,
  const std::vector<double>& temperature_dust,
  const std::vector<double>& radius,
  const std::vector<std::vector<double>>& extinction_coeff,
  const std::vector<std::vector<double>>& absorption_coeff_gas,
  const std::vector<std::vector<double>>& absorption_coeff_dust,
  std::vector<double>& delta_temperature_gas,
  std::vector<double>& delta_temperature_dust)
{
  std::cout << "Linearisation temperature correction.\n\n";

  if (config->use_spline_discretisation)
    throw std::logic_error(
      "linearisedCorrection: only the Taylor moment discretisation is supported "
      "(set use_spline_discretisation = false).\n");

  const size_t D = temperature_gas.size();
  const size_t N = spectral_grid->nbSpectralPoints();

  delta_temperature_gas.assign(D, 0.);
  delta_temperature_dust.assign(D, 0.);

  std::vector<double> radius2(D);
  for (size_t i=0; i<D; ++i) radius2[i] = radius[i]*radius[i];

  //linear weights of the wavelength integration: wavelengthIntegration(f) = sum_n w_n f_n
  const std::vector<double>& lambda = spectral_grid->wavelength_list;
  std::vector<double> w_int(N, 0.);
  if (N >= 2)
  {
    w_int[0]   = -0.5*(lambda[1] - lambda[0]);
    w_int[N-1] = -0.5*(lambda[N-1] - lambda[N-2]);
    for (size_t n=1; n+1<N; ++n)
      w_int[n] = -0.5*(lambda[n+1] - lambda[n-1]);
  }

  //inner-boundary flux injection (frequency-invariant), passed to the operator builder
  const double boundary_flux_correction = rt.boundaryFluxCorrection();

  const size_t M2 = 2*D;
  std::vector<double> A(M2*M2, 0.);
  std::vector<double> rhs_reduced(M2, 0.);

  std::vector<double> C_gas(D, 0.), C_dust(D, 0.);
  //ratio-form local-RE constraint (thesis 3.40/3.41): g_s,i = (sum w kappa_s J)/(sum w kappa_s B) - 1
  std::vector<double> num_gas(D, 0.), num_dust(D, 0.);   //sum_n w_n kappa_s J
  std::vector<double> den_gas(D, 0.), den_dust(D, 0.);   //sum_n w_n kappa_s B(T_s)

  //flux-constancy term (eq. 2.59 integral), gas equation only
  const bool use_flux = config->linearisation_flux_constraint;
  const double xi      = config->linearisation_xi;
  const double flux_target = config->stellar_luminosity / (16. * constants::pi * constants::pi);

  std::vector<double> A_resp(use_flux ? D*M2 : 0, 0.);   //kappa_abs-weighted (V^{-1}) response
  std::vector<double> rhs_flux(use_flux ? D : 0, 0.);

  //zeta(r): grey radial optical depth from the outer boundary (flux-mean extinction),
  //ramped to zeta = tau/(tau + tau_scale)  (~1 deep, ->0 at the thin outer edge)
  std::vector<double> zeta(D, 0.);
  if (use_flux)
  {
    std::vector<double> chi_h(D, 0.);
    for (size_t i=0; i<D; ++i)
      chi_h[i] = radiation_field[i].fluxWeightedExtinction(extinction_coeff[i]);

    std::vector<double> tau(D, 0.);
    for (int i=static_cast<int>(D)-2; i>=0; --i)
      tau[i] = tau[i+1] + 0.5*(chi_h[i] + chi_h[i+1])*(radius[i+1] - radius[i]);

    const double tau_scale = config->linearisation_zeta_tau_scale;
    for (size_t i=0; i<D; ++i)
      zeta[i] = tau[i] / (tau[i] + tau_scale);
  }

  //per-thread accumulators kept in thread-indexed buffers and reduced in fixed order for
  //run-to-run reproducibility (non-associative FP summation)
  const int n_threads = omp_get_max_threads();
  std::vector<std::vector<double>> A_th(n_threads), rhs_th(n_threads);
  std::vector<std::vector<double>> Cg_th(n_threads), Cd_th(n_threads);
  std::vector<std::vector<double>> Ng_th(n_threads), Nd_th(n_threads);
  std::vector<std::vector<double>> Dg_th(n_threads), Dd_th(n_threads);
  std::vector<std::vector<double>> Ar_th(n_threads), rf_th(n_threads);

  #pragma omp parallel
  {
    std::vector<double> A_loc(M2*M2, 0.);
    std::vector<double> rhs_loc(M2, 0.);
    std::vector<double> Cg_loc(D, 0.), Cd_loc(D, 0.);
    std::vector<double> Ng_loc(D, 0.), Nd_loc(D, 0.);
    std::vector<double> Dg_loc(D, 0.), Dd_loc(D, 0.);
    std::vector<double> Ar_loc(use_flux ? D*M2 : 0, 0.);
    std::vector<double> rf_loc(use_flux ? D : 0, 0.);

    aux::TriDiagonalMatrix m(D);
    std::vector<double> source_gas(D, 0.), source_dust(D, 0.);
    std::vector<double> Kvec(D, 0.), MinvK(D, 0.), unit(D, 0.), cinv(D, 0.);

    #pragma omp for schedule(static)
    for (size_t n=0; n<N; ++n)
    {
      const double wl = lambda[n];
      const double wn = w_int[n];

      //RT operator for this frequency: V = m, the source dS/dT, and K = rhs - M J
      rt.buildLinearisedMomentSystem(
        n, radius, radius2, boundary_flux_correction, m, source_gas, source_dust, Kvec);

      m.solveInto(Kvec, MinvK);   //MinvK = V^{-1} K

      //ratio-form local-RE numerators/denominators, dB/dT diagonal, and the K-residual term
      for (size_t i=0; i<D; ++i)
      {
        const double J  = radiation_field[i].mean_intensity[n];
        const double Bg = aux::planckFunctionWavelength(temperature_gas[i],  wl);
        const double Bd = aux::planckFunctionWavelength(temperature_dust[i], wl);
        const double dBg = aux::planckFunctionDerivWavelength(temperature_gas[i],  wl);
        const double dBd = aux::planckFunctionDerivWavelength(temperature_dust[i], wl);

        const double kg = absorption_coeff_gas[i][n];
        const double kd = absorption_coeff_dust[i][n];

        Ng_loc[i] += wn * kg * J;   Dg_loc[i] += wn * kg * Bg;
        Nd_loc[i] += wn * kd * J;   Dd_loc[i] += wn * kd * Bd;

        Cg_loc[i] += wn * kg * dBg;
        Cd_loc[i] += wn * kd * dBd;

        rhs_loc[i]   -= wn * kg * MinvK[i];
        rhs_loc[D+i] -= wn * kd * MinvK[i];

        if (use_flux) rf_loc[i] -= wn * (kg + kd) * MinvK[i];
      }

      //Rybicki contraction: column j of V^{-1} couples a unit dT at j to every node i
      for (size_t j=0; j<D; ++j)
      {
        unit[j] = 1.;
        m.solveInto(unit, cinv);   //cinv[i] = (V^{-1})_{ij}
        unit[j] = 0.;

        const double sg_j = source_gas[j];
        const double sd_j = source_dust[j];

        for (size_t i=0; i<D; ++i)
        {
          const double kg_i = absorption_coeff_gas[i][n];
          const double kd_i = absorption_coeff_dust[i][n];
          const double base = wn * cinv[i];

          const double bg = base * kg_i;
          const double bd = base * kd_i;

          A_loc[(i)*M2     + (j)]     += bg * sg_j;
          A_loc[(i)*M2     + (D+j)]   += bg * sd_j;
          A_loc[(D+i)*M2   + (j)]     += bd * sg_j;
          A_loc[(D+i)*M2   + (D+j)]   += bd * sd_j;
        }

        if (use_flux)
          for (size_t i=0; i<D; ++i)
          {
            const double kb = wn * (absorption_coeff_gas[i][n]
                                  + absorption_coeff_dust[i][n]) * cinv[i];
            Ar_loc[i*M2 + (j)]   += kb * sg_j;
            Ar_loc[i*M2 + (D+j)] += kb * sd_j;
          }
      }
    }

    const int tid = omp_get_thread_num();
    A_th[tid]   = std::move(A_loc);    rhs_th[tid] = std::move(rhs_loc);
    Cg_th[tid]  = std::move(Cg_loc);   Cd_th[tid]  = std::move(Cd_loc);
    Ng_th[tid]  = std::move(Ng_loc);   Nd_th[tid]  = std::move(Nd_loc);
    Dg_th[tid]  = std::move(Dg_loc);   Dd_th[tid]  = std::move(Dd_loc);
    Ar_th[tid]  = std::move(Ar_loc);   rf_th[tid]  = std::move(rf_loc);
  } //omp parallel

  //deterministic reduction: sum the thread partials in ascending thread order
  for (int t=0; t<n_threads; ++t)
  {
    if (A_th[t].empty()) continue;
    for (size_t k=0; k<M2*M2; ++k) A[k] += A_th[t][k];
    for (size_t k=0; k<M2; ++k)    rhs_reduced[k] += rhs_th[t][k];
    for (size_t i=0; i<D; ++i)
    {
      C_gas[i]   += Cg_th[t][i];  C_dust[i]  += Cd_th[t][i];
      num_gas[i] += Ng_th[t][i];  num_dust[i] += Nd_th[t][i];
      den_gas[i] += Dg_th[t][i];  den_dust[i] += Dd_th[t][i];
    }
    if (use_flux)
    {
      for (size_t k=0; k<D*M2; ++k) A_resp[k] += Ar_th[t][k];
      for (size_t i=0; i<D; ++i)    rhs_flux[i] += rf_th[t][i];
    }
  }

  //Composite reduced rows. Flux-constancy (eq. 3.64) is added to the GAS equation only (not
  //duplicated into both species: with the exact V^{-1} the local-RE diagonal cancels to ~0
  //in the diffusion limit, so duplicating the total-flux constraint would make the gas/dust
  //rows collinear there). Deep region where local RE degenerates carries no dust, so flux
  //constancy is the gas-temperature condition there; dust uses its own local RE.
  const double tiny = 1e-300;
  const double inv_target = 1.0 / flux_target;

  std::vector<double> Sbal(D, 0.);
  if (use_flux)
    for (size_t j=0; j<D; ++j)
    {
      Sbal[j] = (num_gas[j] + num_dust[j]) - (den_gas[j] + den_dust[j]);  //int kappa_abs(B-J)
      A_resp[j*M2 + (j)]   -= C_gas[j];
      A_resp[j*M2 + (D+j)] -= C_dust[j];
    }

  std::vector<double> cumJ(use_flux ? M2 : 0, 0.);
  double cumS = 0., cumRK = 0.;

  for (size_t i=0; i<D; ++i)
  {
    const double sg = 1.0 / (den_gas[i]  + std::copysign(tiny, den_gas[i]));
    const double sd = 1.0 / (den_dust[i] + std::copysign(tiny, den_dust[i]));

    if (use_flux && i > 0)
    {
      const double dr = radius[i] - radius[i-1];
      const double ri2 = radius2[i], rm2 = radius2[i-1];
      cumS  += 0.5*(ri2*Sbal[i]     + rm2*Sbal[i-1])     * dr;
      cumRK += 0.5*(ri2*rhs_flux[i] + rm2*rhs_flux[i-1]) * dr;
      for (size_t c=0; c<M2; ++c)
        cumJ[c] += 0.5*(ri2*A_resp[i*M2+c] + rm2*A_resp[(i-1)*M2+c]) * dr;
    }
    const double fz = use_flux ? zeta[i] * inv_target : 0.;

    const double rk_loc_g = rhs_reduced[i];
    const double rk_loc_d = rhs_reduced[D+i];

    for (size_t c=0; c<M2; ++c)
    {
      const double fl = use_flux ? fz * cumJ[c] : 0.;
      A[(i)*M2   + c] = xi * sg * A[(i)*M2   + c] + fl;
      A[(D+i)*M2 + c] = xi * sd * A[(D+i)*M2 + c];
    }

    A[(i)*M2     + (i)]   -= xi * (num_gas[i]  * sg * sg) * C_gas[i];
    A[(D+i)*M2   + (D+i)] -= xi * (num_dust[i] * sd * sd) * C_dust[i];

    rhs_reduced[i]   = xi*(1.0 - num_gas[i]  * sg) + xi*sg*rk_loc_g
                     + fz*(cumRK - cumS);
    rhs_reduced[D+i] = xi*(1.0 - num_dust[i] * sd) + xi*sd*rk_loc_d;
  }

  //Tikhonov regularisation scaled to the mean |diagonal|, to keep the dense solve
  //non-singular where a layer has negligible absorption (empty row)
  double diag_mean = 0.;
  for (size_t k=0; k<M2; ++k) diag_mean += std::abs(A[k*M2+k]);
  diag_mean /= static_cast<double>(M2);
  const double reg = 1e-6 * (diag_mean + 1e-300);
  for (size_t k=0; k<M2; ++k)
    A[k*M2+k] += (A[k*M2+k] >= 0. ? reg : -reg);

  std::vector<double> solution = rhs_reduced;
  const bool ok = solveDense(A, solution, M2);

  if (!ok)
  {
    std::cout << "[linearisation] reduced Newton system singular - skipping step\n";
    return;
  }

  const double omega = config->linearisation_relaxation;
  for (size_t i=0; i<D; ++i)
  {
    delta_temperature_gas[i]  = omega * solution[i];
    delta_temperature_dust[i] = omega * solution[D+i];
  }

  if (std::getenv("LIN_DEBUG"))
  {
    double mg = 0., md = 0.; size_t ig = 0, id = 0;
    double rg_max = 0., rd_max = 0.;
    for (size_t i=0; i<D; ++i)
    {
      const double rg = std::abs(delta_temperature_gas[i]) / temperature_gas[i];
      const double rd = std::abs(delta_temperature_dust[i]) / temperature_dust[i];
      if (rg > mg) { mg = rg; ig = i; }
      if (rd > md) { md = rd; id = i; }
      if (den_gas[i]  != 0.) rg_max = std::max(rg_max, std::abs(num_gas[i] /den_gas[i]  - 1.));
      if (den_dust[i] != 0.) rd_max = std::max(rd_max, std::abs(num_dust[i]/den_dust[i] - 1.));
    }
    double rf_max = 0.;
    if (use_flux)
    {
      double cs = 0.;
      for (size_t i=1; i<D; ++i)
      {
        cs += 0.5*(radius2[i]*Sbal[i] + radius2[i-1]*Sbal[i-1])*(radius[i]-radius[i-1]);
        rf_max = std::max(rf_max, std::abs(cs * inv_target));
      }
    }
    std::cout << "[lin] max |dT_gas|/T = " << mg << " @ " << ig
              << " ,  max |dT_dust|/T = " << md << " @ " << id
              << "  |  max residual: RE_gas = " << rg_max
              << " , RE_dust = " << rd_max
              << " , flux = " << rf_max << "\n";
  }
}


}
