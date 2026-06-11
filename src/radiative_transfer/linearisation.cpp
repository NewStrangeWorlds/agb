
#include "radiative_transfer.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <omp.h>

#include "../spectral_grid/spectral_grid.h"
#include "../config/config.h"
#include "../atmosphere/atmosphere.h"
#include "../additional/physical_const.h"
#include "../additional/tri_diagonal_matrix.h"
#include "../additional/aux_functions.h"


namespace agb{


//Diagnostic: contrast the node-centred flux from calcFlux (eq. 2.58, central first
//derivative of a J at the node) with a FACE-centred flux r^2 H_{i+1/2} = (a J)_{i+1} -
//(a J)_i over Delta X (a = f q r^2). The face-flux divergence equals the solveMomentSystem
//second-order stencil exactly (= eq. 2.59, the local balance), so the face flux is
//conservative by construction. Prints how far each r^2 H_int stays from L/(16 pi^2).
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

  std::vector<double> r2h_face(D-1, 0.);   //face-centred r^2 H_int at face i+1/2
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
  for (size_t i=1; i+2<D; ++i)   //interior faces only
  {
    face_lo = std::min(face_lo, r2h_face[i]); face_hi = std::max(face_hi, r2h_face[i]);
  }

  std::cout << "[FLUXCON] target r2H = " << target << "\n"
            << "[FLUXCON] node-centred (calcFlux, eq2.58): r2H_int in ["
            << node_lo << ", " << node_hi << "]  spread = "
            << (node_hi-node_lo)/target*100. << " %\n"
            << "[FLUXCON] face-centred (conservative, eq2.59): r2H_int in ["
            << face_lo << ", " << face_hi << "]  spread = "
            << (face_hi-face_lo)/target*100. << " %\n";
}


namespace {

//Dense Gaussian elimination with partial pivoting. Solves A x = b in place; A is
//row-major (n*n). Returns false if the matrix is singular.
bool solveDense(std::vector<double>& A, std::vector<double>& b, const size_t n)
{
  for (size_t k=0; k<n; ++k)
  {
    //partial pivot
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



//One Newton step of the full-linearisation temperature correction. The converged
//RT operator (Eddington factor f, sphericality factor q, opacities, geometry) is
//held fixed; only the thermal source term's explicit temperature dependence is
//linearised. For every frequency nu the moment system is
//   M_nu J_nu = rhs_nu(T),   rhs_nu = -c(i) [kappa_gas B(T_gas) + kappa_dust B(T_dust)]
//(+ a frozen inner-boundary flux-injection term), so the response of the mean
//intensity to a temperature change is the EXACT operator response
//   dJ_nu = M_nu^{-1} [ dS_gas .* dT_gas + dS_dust .* dT_dust + K_nu ],
//   dS_{s,i} = d rhs_i / d T_{s,i}  (diagonal),   K_nu = rhs_nu - M_nu J_nu (residual).
//The per-frequency intensities are eliminated (Rybicki) by contracting against the
//two local radiative-equilibrium residuals
//   E_{s,i} = \int kappa_{abs,s} (J_nu - B_nu(T_{s,i})) dnu = 0,   s in {gas, dust},
//giving a dense (2D x 2D) Newton system in (dT_gas, dT_dust), D = nb_grid_points.
//This replaces the approximate (tau-weighted) Unsoeld-Lucy estimate of dJ/dT with
//the exact one, which is what removes the Unsoeld-Lucy accuracy floor.
void RadiativeTransfer::linearisedTemperatureCorrection(
  std::vector<double>& delta_temperature_gas,
  std::vector<double>& delta_temperature_dust)
{
  std::cout << "Linearisation temperature correction.\n\n";

  if (config->use_spline_discretisation)
    throw std::logic_error(
      "linearisedTemperatureCorrection: only the Taylor moment discretisation is "
      "supported (set use_spline_discretisation = false).\n");

  const size_t D = nb_grid_points;
  const size_t N = nb_spectral_points;

  delta_temperature_gas.assign(D, 0.);
  delta_temperature_dust.assign(D, 0.);

  std::vector<double> radius(D), radius2(D);
  for (size_t i=0; i<D; ++i)
  {
    radius[i]  = atmosphere->radius[i];
    radius2[i] = radius[i]*radius[i];
  }

  //linear weights of the wavelength integration: wavelengthIntegration(f) = sum_n w_n f_n
  //(trapezoid on an ascending wavelength axis, with the overall -1 sign the module uses)
  const std::vector<double>& lambda = spectral_grid->wavelength_list;
  std::vector<double> w_int(N, 0.);
  if (N >= 2)
  {
    w_int[0]   = -0.5*(lambda[1] - lambda[0]);
    w_int[N-1] = -0.5*(lambda[N-1] - lambda[N-2]);
    for (size_t n=1; n+1<N; ++n)
      w_int[n] = -0.5*(lambda[n+1] - lambda[n-1]);
  }

  //inner-boundary flux injection is held fixed (its T_gas[0] dependence is not
  //linearised); we still need the values to reconstruct rhs_nu for the residual K_nu
  const double boundary_flux_correction = boundaryFluxCorrection();

  //reduced (2D x 2D) Newton matrix (row-major) and right-hand side, accumulated over
  //frequencies; layout: gas dofs [0,D), dust dofs [D,2D)
  const size_t M2 = 2*D;
  std::vector<double> A(M2*M2, 0.);
  std::vector<double> rhs_reduced(M2, 0.);

  //diagonal d/dT couplings C_{s,i} = \int kappa_{abs,s} dB/dT dnu and the energy-balance
  //residuals E_{s,i}; accumulated over frequencies
  std::vector<double> C_gas(D, 0.), C_dust(D, 0.);
  //numerator/denominator of the RATIO-form local-RE constraint (thesis eq. 3.40/3.41):
  //  g_s,i = (sum_n w_n kappa_s J)/(sum_n w_n kappa_s B(T_s)) - 1.
  //The ratio is O(1) at every depth, whereas the difference sum_n w_n kappa_s (J - B)
  //spans many orders of magnitude with depth and makes the dense system ill-scaled.
  std::vector<double> num_gas(D, 0.), num_dust(D, 0.);   //sum_n w_n kappa_s J
  std::vector<double> den_gas(D, 0.), den_dust(D, 0.);   //sum_n w_n kappa_s B(T_s)

  //Flux-constancy term in the INTEGRATED (eq. 2.59) form, J-consistent by construction:
  //  d(r^2 H_int)/dr = r^2 \int kappa_abs (B - J) dnu   (eq. 2.59 = the local balance),
  //so  r^2 H_int,i = r*^2 H* + \int_{r1}^{ri} r'^2 [\int kappa_abs (B-J) dnu] dr'  and
  //  g_flux,i = (r^2 H_int,i)/(r*^2 H*) - 1 = (1/target) * cumtrap_i( r^2 S ),
  //  S_k = sum_n w_n [kappa_g(J-B_g) + kappa_d(J-B_d)]  (= num_total - den_total).
  //This is the radial integral of the ABSOLUTE local balance, normalised to O(1) by the
  //target - well scaled, unlike S_k itself - and its divergence equals the moment-system
  //stencil exactly, so it does not fight the local RE. A_resp[k][col] holds the kappa_abs-
  //weighted (V^{-1}) intensity response at output node k (the dS_k/dT J-part).
  const bool use_flux = config->linearisation_flux_constraint;
  const double xi      = config->linearisation_xi;
  const double flux_target = config->stellar_luminosity / (16. * constants::pi * constants::pi);

  std::vector<double> A_resp(use_flux ? D*M2 : 0, 0.);  //sum_n w_n kappa_abs_k (V^{-1})_{kj} source_t
  std::vector<double> rhs_flux(use_flux ? D : 0, 0.);   //-sum_n w_n kappa_abs_k (M^{-1} K)_k

  //zeta(r): grey radial optical depth from the outer boundary (flux-mean extinction),
  //ramped to zeta = tau/(tau + tau_scale)  (~1 deep, ->0 at the thin outer edge)
  std::vector<double> zeta(D, 0.);
  if (use_flux)
  {
    std::vector<double> chi_h(D, 0.);
    for (size_t i=0; i<D; ++i)
      chi_h[i] = radiation_field[i].fluxWeightedExtinction(atmosphere->extinction_coeff[i]);

    std::vector<double> tau(D, 0.);
    for (int i=static_cast<int>(D)-2; i>=0; --i)
      tau[i] = tau[i+1] + 0.5*(chi_h[i] + chi_h[i+1])*(radius[i+1] - radius[i]);

    const double tau_scale = config->linearisation_zeta_tau_scale;
    for (size_t i=0; i<D; ++i)
      zeta[i] = tau[i] / (tau[i] + tau_scale);
  }

  //Per-thread accumulators kept in thread-indexed buffers and reduced in a FIXED order
  //below, rather than via an arrival-order omp-critical. Together with schedule(static)
  //this makes the frequency sum into A/rhs bit-reproducible run-to-run (floating-point
  //addition is non-associative, so the summation order must be fixed for determinism).
  const int n_threads = omp_get_max_threads();
  std::vector<std::vector<double>> A_th(n_threads), rhs_th(n_threads);
  std::vector<std::vector<double>> Cg_th(n_threads), Cd_th(n_threads);
  std::vector<std::vector<double>> Ng_th(n_threads), Nd_th(n_threads);
  std::vector<std::vector<double>> Dg_th(n_threads), Dd_th(n_threads);
  std::vector<std::vector<double>> Ar_th(n_threads), rf_th(n_threads);

  #pragma omp parallel
  {
    //per-thread accumulators (reduced into the shared ones at the end)
    std::vector<double> A_loc(M2*M2, 0.);
    std::vector<double> rhs_loc(M2, 0.);
    std::vector<double> Cg_loc(D, 0.), Cd_loc(D, 0.);
    std::vector<double> Ng_loc(D, 0.), Nd_loc(D, 0.);
    std::vector<double> Dg_loc(D, 0.), Dd_loc(D, 0.);
    std::vector<double> Ar_loc(use_flux ? D*M2 : 0, 0.);   //kappa_abs-weighted response
    std::vector<double> rf_loc(use_flux ? D : 0, 0.);

    aux::TriDiagonalMatrix m(D);
    m.resize(D);   //the size ctor does not allocate the solveInto() scratch buffers
    std::vector<double> rhs_nu(D, 0.), unit(D, 0.), cinv(D, 0.);
    std::vector<double> source_gas(D, 0.), source_dust(D, 0.);
    std::vector<double> dB_gas(D, 0.), dB_dust(D, 0.);
    std::vector<double> Jvec(D, 0.), Kvec(D, 0.), MinvK(D, 0.);

    #pragma omp for schedule(static)
    for (size_t n=0; n<N; ++n)
    {
      const double wl = lambda[n];
      const double wn = w_int[n];

      const double boundary_planck_derivative =
        aux::planckFunctionDerivWavelength(atmosphere->temperature_gas[0], wl);

      const std::vector<double>& x_grid = generateXGrid(
        extinction_coeff[n], radius, sphericality_factor[n]);

      //assemble M_nu and rhs_nu exactly as the forward solve does
      assembleMomentSystemTaylor(
        x_grid, radius, radius2,
        planck_emission[n], extinction_coeff[n], scattering_coeff[n],
        eddington_factor_f[n], boundary_eddington_factor_h[n], sphericality_factor[n],
        boundary_planck_derivative, boundary_flux_correction,
        m, rhs_nu);

      //per-node Planck derivatives and the diagonal source d rhs_i / d T_{s,i}.
      //c(i) = radius2[i] / (q[i] chi[i]); the boundary rows carry an extra hr/2, hl/2
      //factor (matching the Taylor assembly), the inner-boundary flux term is frozen.
      for (size_t i=0; i<D; ++i)
      {
        dB_gas[i]  = aux::planckFunctionDerivWavelength(atmosphere->temperature_gas[i],  wl);
        dB_dust[i] = aux::planckFunctionDerivWavelength(atmosphere->temperature_dust[i], wl);

        const double c_i = radius2[i] / sphericality_factor[n][i] / extinction_coeff[n][i];

        double face = 1.0;
        if (i == 0)     face = 0.5*(x_grid[1] - x_grid[0]);
        else if (i+1==D) face = 0.5*(x_grid[D-1] - x_grid[D-2]);

        source_gas[i]  = - face * c_i * atmosphere->absorption_coeff_gas[i][n]  * dB_gas[i];
        source_dust[i] = - face * c_i * atmosphere->absorption_coeff_dust[i][n] * dB_dust[i];
      }

      //current mean intensity and the moment-equation residual K = rhs - M J
      for (size_t i=0; i<D; ++i) Jvec[i] = radiation_field[i].mean_intensity[n];

      Kvec[0]   = rhs_nu[0]   - (m.b[0]*Jvec[0] + m.c[0]*Jvec[1]);
      for (size_t i=1; i+1<D; ++i)
        Kvec[i] = rhs_nu[i]   - (m.a[i]*Jvec[i-1] + m.b[i]*Jvec[i] + m.c[i]*Jvec[i+1]);
      Kvec[D-1] = rhs_nu[D-1] - (m.a[D-1]*Jvec[D-2] + m.b[D-1]*Jvec[D-1]);

      m.solveInto(Kvec, MinvK);

      //energy-balance residuals and the dB/dT diagonal couplings
      for (size_t i=0; i<D; ++i)
      {
        const double Bg = aux::planckFunctionWavelength(atmosphere->temperature_gas[i],  wl);
        const double Bd = aux::planckFunctionWavelength(atmosphere->temperature_dust[i], wl);

        const double kg = atmosphere->absorption_coeff_gas[i][n];
        const double kd = atmosphere->absorption_coeff_dust[i][n];

        Ng_loc[i] += wn * kg * Jvec[i];   Dg_loc[i] += wn * kg * Bg;
        Nd_loc[i] += wn * kd * Jvec[i];   Dd_loc[i] += wn * kd * Bd;

        Cg_loc[i] += wn * kg * dB_gas[i];
        Cd_loc[i] += wn * kd * dB_dust[i];

        //residual contribution R_{s,i} = sum_n w_n kappa_s (M^{-1} K)_i  -> moves to RHS
        //(still unscaled here; the 1/den_s row scaling is applied in the post-loop)
        rhs_loc[i]   -= wn * kg * MinvK[i];
        rhs_loc[D+i] -= wn * kd * MinvK[i];

        //kappa_abs-weighted moment residual for the flux term (kappa_abs = kappa_g + kappa_d)
        if (use_flux) rf_loc[i] -= wn * (kg + kd) * MinvK[i];
      }

      //Rybicki contraction: for each source node j, column j of M_nu^{-1} couples the
      //unit temperature perturbation there to the residual at every node i
      for (size_t j=0; j<D; ++j)
      {
        unit[j] = 1.;
        m.solveInto(unit, cinv);   //cinv[i] = (M_nu^{-1})_{ij}
        unit[j] = 0.;

        const double sg_j = source_gas[j];
        const double sd_j = source_dust[j];

        for (size_t i=0; i<D; ++i)
        {
          const double kg_i = atmosphere->absorption_coeff_gas[i][n];
          const double kd_i = atmosphere->absorption_coeff_dust[i][n];
          const double base = wn * cinv[i];

          const double bg = base * kg_i;   //gas residual row weight
          const double bd = base * kd_i;   //dust residual row weight

          A_loc[(i)*M2     + (j)]     += bg * sg_j;   //d E_gas,i / d T_gas,j
          A_loc[(i)*M2     + (D+j)]   += bg * sd_j;   //d E_gas,i / d T_dust,j
          A_loc[(D+i)*M2   + (j)]     += bd * sg_j;   //d E_dust,i / d T_gas,j
          A_loc[(D+i)*M2   + (D+j)]   += bd * sd_j;   //d E_dust,i / d T_dust,j
        }

        //kappa_abs-weighted intensity response at each output node i (the J-part of
        //dS_i/dT_t,j); the radial integral and the direct -C term are applied post-loop
        if (use_flux)
          for (size_t i=0; i<D; ++i)
          {
            const double kb = wn * (atmosphere->absorption_coeff_gas[i][n]
                                  + atmosphere->absorption_coeff_dust[i][n]) * cinv[i];
            Ar_loc[i*M2 + (j)]   += kb * sg_j;
            Ar_loc[i*M2 + (D+j)] += kb * sd_j;
          }
      }
    }

    //hand this thread's partials to its own slot (no race: unique tid); the sum across
    //threads is done in fixed order after the region for run-to-run reproducibility
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
    if (A_th[t].empty()) continue;   //a thread that ran no iterations
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

  //Assemble the composite reduced rows. The flux-constancy term (thesis eq. 3.64) is
  //added to the GAS equation ONLY, not duplicated into both species. Reason: with the
  //exact V^{-1} response the local-RE diagonal cancels to ~0 in the diffusion limit, so
  //if the same total-flux constraint were placed in both the gas and dust rows they would
  //become collinear there (det -> 0, divergent spurious dT split). Physically the deep
  //region where local RE degenerates carries no dust (too hot to condense), so flux
  //constancy is the gas-temperature condition there; dust is set by its own local RE
  //wherever it exists (the cooler outer region, where local RE is well conditioned).
  //  gas row : A = xi/den_g * (local J-response)  + zeta_i/(r*^2 H*) * (flux-op response)
  //  dust row: A = xi/den_d * (local J-response)   [pure ratio-form local RE, B.9-B.12]
  const double tiny = 1e-300;
  const double inv_target = 1.0 / flux_target;

  //fold the direct dB/dT term into the kappa_abs-weighted response so A_resp becomes the
  //full dS_k/dT_t,j (J-response minus the local Planck term at the diagonal)
  std::vector<double> Sbal(D, 0.);
  if (use_flux)
    for (size_t j=0; j<D; ++j)
    {
      Sbal[j] = (num_gas[j] + num_dust[j]) - (den_gas[j] + den_dust[j]);  //S_j = int kappa_abs(B-J)
      A_resp[j*M2 + (j)]   -= C_gas[j];
      A_resp[j*M2 + (D+j)] -= C_dust[j];
    }

  //running radial cumulative trapezoids (eq. 2.59 integral from the inner boundary):
  //  cumJ[c]  = int r^2 dS/dT_c dr,   cumS = int r^2 S dr,   cumRK = int r^2 (-RK_flux) dr
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

    const double rk_loc_g = rhs_reduced[i];     //= -RK_loc_gas,i
    const double rk_loc_d = rhs_reduced[D+i];   //= -RK_loc_dust,i

    for (size_t c=0; c<M2; ++c)
    {
      const double fl = use_flux ? fz * cumJ[c] : 0.;
      A[(i)*M2   + c] = xi * sg * A[(i)*M2   + c] + fl;   //gas row (local RE + flux)
      A[(D+i)*M2 + c] = xi * sd * A[(D+i)*M2 + c];        //dust row (local RE only)
    }

    A[(i)*M2     + (i)]   -= xi * (num_gas[i]  * sg * sg) * C_gas[i];
    A[(D+i)*M2   + (D+i)] -= xi * (num_dust[i] * sd * sd) * C_dust[i];

    //gas RHS: -R - (dR/dJ)V^{-1}K, with R = xi*g_loc + zeta*g_flux, g_flux = cumS/target
    rhs_reduced[i]   = xi*(1.0 - num_gas[i]  * sg) + xi*sg*rk_loc_g
                     + fz*(cumRK - cumS);                 //gas: local RE + flux
    rhs_reduced[D+i] = xi*(1.0 - num_dust[i] * sd) + xi*sd*rk_loc_d;   //dust: local RE only
  }

  //Tikhonov regularisation, scaled to the mean |diagonal|: with the ratio form the rows
  //are O(1)-scaled, but a layer with negligible absorption (kappa_abs -> 0) still leaves
  //an essentially empty row; the regularisation keeps the dense solve non-singular there
  //(that layer's temperature is then simply left almost unchanged) without affecting the
  //well-determined rows.
  double diag_mean = 0.;
  for (size_t k=0; k<M2; ++k) diag_mean += std::abs(A[k*M2+k]);
  diag_mean /= static_cast<double>(M2);
  const double reg = 1e-6 * (diag_mean + 1e-300);
  for (size_t k=0; k<M2; ++k)
    A[k*M2+k] += (A[k*M2+k] >= 0. ? reg : -reg);   //grow |diagonal| away from zero

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
    double rg_max = 0., rd_max = 0.;   //max |constraint residual| g_s = num/den - 1
    for (size_t i=0; i<D; ++i)
    {
      const double rg = std::abs(delta_temperature_gas[i]) / atmosphere->temperature_gas[i];
      const double rd = std::abs(delta_temperature_dust[i]) / atmosphere->temperature_dust[i];
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
        rf_max = std::max(rf_max, std::abs(cs * inv_target));   //|r^2 H_int/target - 1|
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
