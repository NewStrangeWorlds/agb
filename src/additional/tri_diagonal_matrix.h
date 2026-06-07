 
#ifndef _tridiagonal_matrix_h
#define _tridiagonal_matrix_h

#include <vector>
#include <stdexcept>

namespace agb{ namespace aux{

//a tri-diagonal matrix
//b is main diagonal,
//a the one below, c the one above
struct TriDiagonalMatrix{
  TriDiagonalMatrix(const size_t size)
  {
    a.assign(size, 0.);
    c.assign(size, 0.);
    b.assign(size, 0.);
  }

  std::vector<double> a;
  std::vector<double> b;
  std::vector<double> c;

  //reusable scratch for the allocation-free solveInto()
  std::vector<double> cp_scratch;
  std::vector<double> dp_scratch;

  //resize all diagonals (and scratch) to n; reused buffers keep their capacity
  //so repeated calls with n <= capacity do not reallocate
  void resize(const size_t n)
  {
    a.resize(n);
    b.resize(n);
    c.resize(n);
    cp_scratch.resize(n);
    dp_scratch.resize(n);
  }

  //solves M*x = rhs into the caller-provided buffer x, reusing cp_scratch /
  //dp_scratch instead of allocating. Same Thomas algorithm / arithmetic order as
  //solve(), so results are bit-identical.
  void solveInto(const std::vector<double>& rhs, std::vector<double>& x)
  {
    const size_t n = b.size();

    if (rhs.size() != n)
      throw std::logic_error("TriDiagonal matrix solveInto: matrix and rhs vector need to have the same size!\n");

    x.resize(n);

    double* cp = cp_scratch.data();
    double* dp = dp_scratch.data();

    cp[0] = c[0]/b[0];
    dp[0] = rhs[0]/b[0];

    for (size_t i=1; i<n; ++i)
    {
      const double m = b[i] - a[i]*cp[i-1];
      cp[i] = c[i]/m;
      dp[i] = (rhs[i] - a[i]*dp[i-1])/m;
    }

    x[n-1] = dp[n-1];

    for (int i=n-2; i>-1; --i)
      x[i] = dp[i] - cp[i]*x[i+1];
  }

  //solves M*x = rhs
  std::vector<double> solve(const std::vector<double>& rhs)
  {
    if (rhs.size() != a.size())
      throw std::logic_error("TriDiagonal matrix solve: matrix and rhs vector need to have the same size!\n");

    /*std::vector<double> d(a.size(), 0.);
    std::vector<double> v(a.size(), 0.);
    
    v[0] = b[0];
    d[0] = rhs[0];

    for (size_t i=1; i<a.size(); ++i)
    {
      v[i] = b[i] - a[i]/v[i-1] * c[i-1];
	    d[i] = rhs[i] - a[i]/v[i-1] * d[i-1];
    }
    
    std::vector<double> x(a.size(), 0.);

    x.back() = d.back()/v.back();

    for (int i=a.size()-2; i>-1; --i)
      x[i] = (d[i] - c[i]*x[i+1])/v[i];*/

    std::vector<double> cp(b.size(), 0);
    std::vector<double> dp(b.size(), 0);

    cp[0] = c[0]/b[0];
    dp[0] = rhs[0]/b[0];

    for (size_t i=1; i<b.size(); ++i)
    {
      cp[i] = c[i]/(b[i] - a[i]*cp[i-1]);
      dp[i] = (rhs[i] - a[i]*dp[i-1])/(b[i] - a[i]*cp[i-1]);
    }


    std::vector<double> x(a.size(), 0.);

    x.back() = dp.back();

    for (int i=a.size()-2; i>-1; --i)
      x[i] = dp[i] - cp[i]*x[i+1];

    return x;
  }
};


}
}

#endif