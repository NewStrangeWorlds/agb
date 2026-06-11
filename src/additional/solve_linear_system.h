
#ifndef _interpolation_h
#define _interpolation_h

#include <vector>
#include <cmath>
#include <algorithm>

namespace agb { namespace aux {

//Dense Gaussian elimination with partial pivoting for the small (m x m) Anderson
//least-squares system. Returns an empty vector if the matrix is singular.
std::vector<double> solveLinearSystem(
  std::vector<std::vector<double>> a,
  std::vector<double> b)
{
  const size_t m = b.size();

  for (size_t col=0; col<m; ++col)
  {
    size_t pivot = col;
    for (size_t row=col+1; row<m; ++row)
      if (std::abs(a[row][col]) > std::abs(a[pivot][col])) pivot = row;

    if (std::abs(a[pivot][col]) < 1e-300)
      return std::vector<double>{};

    std::swap(a[col], a[pivot]);
    std::swap(b[col], b[pivot]);

    for (size_t row=col+1; row<m; ++row)
    {
      const double factor = a[row][col] / a[col][col];
      for (size_t k=col; k<m; ++k) a[row][k] -= factor * a[col][k];
      b[row] -= factor * b[col];
    }
  }

  std::vector<double> x(m, 0.);
  for (size_t row=m; row-- > 0; )
  {
    double sum = b[row];
    for (size_t k=row+1; k<m; ++k) sum -= a[row][k] * x[k];
    x[row] = sum / a[row][row];
  }

  return x;
}



}}

#endif
