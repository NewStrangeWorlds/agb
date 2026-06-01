

#include <iostream>
#include <omp.h>
#include <csignal>
#include <cstdlib>

#include "../agb_model/agb_model.h"
#include "../hydrodynamics/structure_solver.h"


int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    std::cout << "wrong number of command line parameters!\n";

    return 1;
  }

  std::string model_folder = argv[1];

  //temporary: validate the Henyey structure-solver autodiff + Newton wiring
  if (model_folder == "--selftest")
  {
    agb::StructureSolver solver(50);

    const double jac_err = solver.selfTestPhiJacobian();
    std::cout << "StructureSolver Phi-Jacobian AD-vs-FD max rel error: " << jac_err << "\n";

    int iters = 0;
    const double res = solver.selfTestPhiSolve(&iters);
    std::cout << "StructureSolver Phi Newton solve: residual " << res
              << " after " << iters << " iterations\n";

    //coupled system mixes scales over ~20 orders of magnitude (Phi ~1e5 vs
    //K_j ~1e-13), so the central finite-difference reference is only good to
    //~(machine_eps)^(1/3) ~ 1e-5 on the tiny cross-derivatives; AD is exact.
    const double coupled_jac_err = solver.selfTestCoupledJacobian();
    std::cout << "StructureSolver coupled (Phi+dust) Jacobian AD-vs-FD max rel error: "
              << coupled_jac_err << "\n";

    const double eigen_jac_err = solver.selfTestEigenJacobian();
    std::cout << "StructureSolver eigenvalue (Mdot) Jacobian AD-vs-FD max rel error: "
              << eigen_jac_err << "\n";

    agb::StructureSolver tsolver(80);
    double mdot = 0.; int cp = 0, eig_iters = 0;
    const double eig_red = tsolver.selfTestEigenSolve(&mdot, &cp, &eig_iters);
    std::cout << "StructureSolver transonic eigen solve: rel residual reduction "
              << eig_red << " in " << eig_iters << " iters; Mdot=" << mdot
              << " g/s, critical point index=" << cp << "/80\n";

    agb::StructureSolver gsolver(80);
    int gcp = 0; bool gmono = false;
    const double gmdot = gsolver.selfTestGreyBootstrap(&gcp, &gmono);
    std::cout << "StructureSolver grey bootstrap: Mdot=" << gmdot
              << " g/s, critical point index=" << gcp << "/80, T monotone="
              << (gmono ? "yes" : "no") << "\n";

    return (jac_err < 1e-6 && res < 1e-9
            && coupled_jac_err < 1e-4 && eigen_jac_err < 1e-4) ? 0 : 1;
  }


  agb::AGBStarModel agb_wind(model_folder);

  agb_wind.calcModel();

  return 0;
}