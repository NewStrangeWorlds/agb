
#ifndef _structure_solver_h
#define _structure_solver_h

#include <vector>
#include <cmath>

namespace agb {


//Parameters for the grey (Lucy) bootstrap starting model. These collapse the
//frequency-resolved RT + Mie to a single grey opacity + analytic temperature so
//a smooth structure can be produced to seed the full Newton iteration.
struct GreyParameters{
  double stellar_radius = 0.;            //R_* [cm], for T_eff and dilution W(r)
  double gas_opacity = 1.0e-4;           //grey gas opacity kappa_gas [cm^2/g]
  double dust_q_ext = 1.5;               //grey dust extinction efficiency
  double monomer_radius = 1.28e-8;       //a_0 [cm]
  double mean_molecular_weight = 1.3;    //mu, for n<H> = rho/(mu m_H)
  double condensation_temperature = 1500.;  //dust switches on below this T [K]
  double nucleation_per_h_scale = 1.0e-12;  //J*/n<H> scale below T_cond [1/s]
  double dust_growth_timescale = 1.0e6;  //tau [s]
};





//Henyey-type relaxation solver for the stationary wind structure block.
//
//Built incrementally toward the full coupled system (Melia-Phi wind equation +
//Gail-Sedlmayr dust moments, with Mdot/R0 as transonic eigenvalues). The current
//stage solves the Melia-Phi wind equation by global Newton-Raphson relaxation at
//frozen flux-mean extinction chi_F(r) and sound speed c_T(r), with the mass-loss
//rate prescribed. The radiative acceleration is computed inside the residual from
//  alpha(r) = chi_F * L_star * r^2 * v / (c * G M_star * Mdot)
//(so that Mdot can later become a genuine eigenvalue), rather than frozen.
//
//The residual is a function template so the same code is evaluated in double
//precision (Newton residual) and in CppAD::AD<double> (to tape the Jacobian).
class StructureSolver{
  public:
    StructureSolver(const size_t nb_points);

    //frozen-coefficient inputs for one structure solve
    std::vector<double> radius;              //r [cm]
    std::vector<double> sound_speed;         //isothermal sound speed c_T(r) [cm/s]
    std::vector<double> flux_mean_extinction;//chi_F(r) [1/cm], frozen from the RT
    std::vector<double> sound_speed2_deriv;  //d(c_T^2)/dr [cm/s^2], frozen

    //frozen dust source coefficients (n<H>-normalized), from the Gail-Sedlmayr
    //nucleation/growth kernels evaluated at the current (frozen) T and chemistry
    std::vector<double> nucleation_per_h;    //J* / n<H>  [1/s]
    std::vector<double> growth_timescale;    //tau [s] (== GailSedlmayrDust growth_rate[])

    double stellar_luminosity = 0.;          //L_star [erg/s]
    double gravitational_mass = 0.;          //G * M_star [cgs]
    double mass_loss_rate = 0.;              //Mdot [g/s] (prescribed at this stage)
    double inner_velocity = 0.;              //v at the inner boundary [cm/s]

    int critical_point = 0;                  //grid index of the sonic point (branch split)

    static constexpr int nb_moments = 6;     //dust moments K0..K5
    double min_monomer_number = 1000.;       //N_l (lower integration bound)

    //if true, residualEigen pins Mdot to exp(log_mass_loss_rate_target) (replacing
    //the critical-point regularity row) -> prescribed-Mdot mode (Setup B), instead
    //of treating Mdot as the transonic eigenvalue (Setup A)
    bool fix_mass_loss_rate = false;
    double log_mass_loss_rate_target = 0.;

    //if false, solveEigen keeps critical_point fixed (does not relocate it each
    //Newton step). Used by the prescribed-Mdot solve, where the sonic node is held
    //fixed for stability and re-placed by an outer critical-point iteration.
    bool relocate_critical_point = true;

    //residual of the discretized Melia-Phi wind equation alone (M unknowns)
    template<class Scalar>
    std::vector<Scalar> residualPhi(const std::vector<Scalar>& phi) const;

    //residual of the coupled wind + dust-moment system. The unknown vector is
    //laid out in blocks: [ Phi(M), K0(M), K1(M), ..., K5(M) ]  -> 7*M entries.
    template<class Scalar>
    std::vector<Scalar> residualCoupled(const std::vector<Scalar>& x) const;

    //coupled system with Mdot promoted to a transonic eigenvalue. Unknown vector
    //x = [ Phi(M), K0..K5(M), s ] with s = ln(Mdot); the extra residual is the
    //critical-point regularity that closes Mdot.  -> 7*M + 1 entries.
    template<class Scalar>
    std::vector<Scalar> residualEigen(const std::vector<Scalar>& x) const;

    //damped Newton-Raphson relaxation of residualPhi. Returns the converged Phi;
    //writes the achieved residual norm and iteration count if pointers are given.
    std::vector<double> solvePhi(
      std::vector<double> phi,
      const int max_iterations = 100,
      const double tolerance = 1e-10,
      double* final_residual = nullptr,
      int* iterations = nullptr);

    //damped Newton-Raphson relaxation of the eigenvalue system residualEigen.
    //The unknown vector is [Phi(M), K0..K5(M), ln(Mdot)]; the critical point is
    //relocated each step to the interior node where Phi is closest to c_T.
    //Returns the converged unknown vector.
    std::vector<double> solveEigen(
      std::vector<double> x,
      const int max_iterations = 200,
      const double tolerance = 1e-8,
      double* final_residual = nullptr,
      int* iterations = nullptr);

    //evaluate the L2 norm of residualEigen at a given state (after relocating the
    //critical point), for consistency checks against an externally-provided seed
    double residualEigenNorm(const std::vector<double>& x);

    //reconstruct the wind velocity v(r) from a solution vector's Phi block, using
    //the current critical_point for the sub-/super-sonic branch split
    std::vector<double> velocityProfile(const std::vector<double>& x) const;

    //interior grid node where Phi is closest to c_T (the discrete sonic point)
    int sonicNode(const std::vector<double>& x) const;

    //Grey (Lucy) bootstrap: self-consistently iterate the grey opacity -> optical
    //depth -> Lucy temperature -> radiative acceleration -> wind+dust structure
    //(via solveEigen) until the structure converges. Produces a smooth starting
    //model. radius[] must be set beforehand. On return, sound_speed/flux_mean_
    //extinction/etc. hold the grey state and the converged unknowns are returned
    //(layout [Phi(M), K0..K5(M), ln Mdot]); temperature(r) is written to T_out.
    std::vector<double> greyBootstrap(
      const GreyParameters& params,
      std::vector<double>& T_out,
      const int max_grey_iterations = 50,
      const double tolerance = 1e-4,
      int* grey_iterations = nullptr);

    //run greyBootstrap on a synthetic IRC10216-scale setup; returns the converged
    //Mdot [g/s] and (via pointers) the interior critical-point index and a flag
    //for whether the temperature profile is monotonically decreasing outward
    double selfTestGreyBootstrap(int* cp_out = nullptr, bool* monotone_T = nullptr);

    //--- verification helpers ---
    //max relative discrepancy between the CppAD and finite-difference Jacobians
    double selfTestPhiJacobian();
    //solve a synthetic Phi system and return the final residual L2 norm
    double selfTestPhiSolve(int* iterations = nullptr);
    //max relative AD-vs-FD Jacobian discrepancy for the coupled (Phi+dust) system
    double selfTestCoupledJacobian();
    //max relative AD-vs-FD Jacobian discrepancy for the eigenvalue (Mdot) system
    double selfTestEigenJacobian();
    //solve a synthetic transonic eigenvalue problem; returns the relative residual
    //reduction and (via pointers) the converged Mdot, critical-point index, iters
    double selfTestEigenSolve(double* mdot_out = nullptr,
                              int* cp_out = nullptr, int* iters_out = nullptr);

  private:
    size_t nb_points = 0;

    //wind velocity from the Melia variable: v = phi -/+ sqrt(phi^2 - c_T^2),
    //sub-sonic (minus) inside the critical point, super-sonic (plus) outside
    template<class Scalar>
    Scalar windVelocity(const Scalar& phi, const double c_t, const bool supersonic) const;

    //fill the wind + dust-moment residual blocks for a given (possibly unknown)
    //mass-loss rate; shared by residualCoupled and residualEigen
    template<class Scalar>
    void fillStructureResidual(
      const std::vector<Scalar>& x, const Scalar& mdot,
      std::vector<Scalar>& residual) const;

    //fill a synthetic, physically-scaled frozen state (used by the self-tests)
    void setSyntheticState();
};


}

#endif
