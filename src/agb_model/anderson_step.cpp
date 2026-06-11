 
#include "agb_model.h"

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "../additional/aux_functions.h"
#include "../additional/physical_const.h"
#include "../additional/solve_linear_system.h"


namespace agb{


std::vector<double> AGBStarModel::andersonStep(
  std::vector<std::vector<double>>& x_history,
  std::vector<std::vector<double>>& f_history,
  const std::vector<double>& x_k,
  const std::vector<double>& f_k)
{
  const size_t n = x_k.size();

  //record the current iterate/residual and trim to the window (m differences need
  //m+1 stored pairs)
  x_history.push_back(x_k);
  f_history.push_back(f_k);

  while (x_history.size() > config.anderson_window + 1)
  {
    x_history.erase(x_history.begin());
    f_history.erase(f_history.begin());
  }

  const size_t stored = x_history.size();

  //plain (damped) step G(x_k) = x_k + f_k until we have at least one difference
  std::vector<double> result(n);
  for (size_t i=0; i<n; ++i)
    result[i] = x_k[i] + f_k[i];

  if (stored < 2)
    return result;

  const size_t m = stored - 1; //number of difference columns

  //difference matrices: dF[j] = f_{j+1} - f_j, dX[j] = x_{j+1} - x_j
  std::vector<std::vector<double>> dF(m, std::vector<double>(n));
  std::vector<std::vector<double>> dX(m, std::vector<double>(n));

  for (size_t j=0; j<m; ++j)
    for (size_t i=0; i<n; ++i)
    {
      dF[j][i] = f_history[j+1][i] - f_history[j][i];
      dX[j][i] = x_history[j+1][i] - x_history[j][i];
    }

  //normal equations (dF^T dF) gamma = dF^T f_k  with light Tikhonov regularisation
  std::vector<std::vector<double>> a(m, std::vector<double>(m, 0.));
  std::vector<double> b(m, 0.);

  for (size_t j=0; j<m; ++j)
  {
    for (size_t k=0; k<m; ++k)
    {
      double sum = 0.;
      for (size_t i=0; i<n; ++i) sum += dF[j][i] * dF[k][i];
      a[j][k] = sum;
    }
    double sum = 0.;
    for (size_t i=0; i<n; ++i) sum += dF[j][i] * f_k[i];
    b[j] = sum;
  }

  double trace = 0.;
  for (size_t j=0; j<m; ++j) trace += a[j][j];
  const double reg = 1e-10 * (trace/m + 1e-300);
  for (size_t j=0; j<m; ++j) a[j][j] += reg;

  std::vector<double> gamma = aux::solveLinearSystem(a, b);

  //if the solve failed (singular), fall back to the plain damped step
  if (gamma.empty())
    return result;

  //x_{k+1} = (x_k + f_k) - sum_j gamma_j (dX_j + dF_j)
  for (size_t j=0; j<m; ++j)
    for (size_t i=0; i<n; ++i)
      result[i] -= gamma[j] * (dX[j][i] + dF[j][i]);

  return result;
}


}