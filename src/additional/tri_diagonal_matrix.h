 
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