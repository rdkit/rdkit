
// solving A * X = B
// in two steps -- factor (getrf()) and solve (getrs())

//#define BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 
//#define BOOST_NO_FUNCTION_TEMPLATE_ORDERING

#include <cstddef>
#include <iostream>
#include <boost/numeric/bindings/atlas/cblas.hpp>
#include <boost/numeric/bindings/atlas/clapack.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/traits/std_vector.hpp>
#include "utils.h"

namespace ublas = boost::numeric::ublas;
namespace atlas = boost::numeric::bindings::atlas;

using std::size_t; 
using std::cout;
using std::endl; 

#ifndef F_ROW_MAJOR
typedef ublas::matrix<double, ublas::column_major> m_t;
#else
typedef ublas::matrix<double, ublas::row_major> m_t;
#endif

int main() {

  cout << endl; 

  size_t n = 5; 
  m_t a (n, n);   // system matrix 

  size_t nrhs = 2; 
  m_t x (n, nrhs); 
  // b -- right-hand side matrix:
  // .. see leading comments for `gesv()' in clapack.hpp
#ifndef F_ROW_MAJOR
  m_t b (n, nrhs);
#else
  m_t b (nrhs, n);
#endif

  init_symm (a); 
  //     [n   n-1 n-2  ... 1]
  //     [n-1 n   n-1  ... 2]
  // a = [n-2 n-1 n    ... 3]
  //     [        ...       ]
  //     [1   2   ...  n-1 n]
  ublas::matrix_row<m_t> ar1 (a, 0), ar3 (a, 3);
  swap (ar1, ar3);   // swap rows to force pivoting 

  ublas::matrix_column<m_t> xc0 (x, 0), xc1 (x, 1); 
  atlas::set (1., xc0);  // x[.,0] = 1
  atlas::set (2., xc1);  // x[.,1] = 2
#ifndef F_ROW_MAJOR
  atlas::gemm (a, x, b);  // b = a x, so we know the result ;o) 
#else
  // see leading comments for `gesv()' in clapack.hpp
  ublas::matrix_row<m_t> br0 (b, 0), br1 (b, 1); 
  atlas::gemv (a, xc0, br0);  // b[0,.] = a x[.,0]
  atlas::gemv (a, xc1, br1);  // b[1,.] = a x[.,1]  =>  b^T = a x
#endif 

  print_m (a, "A"); 
  cout << endl; 
  print_m (b, "B"); 
  cout << endl; 

  std::vector<int> ipiv (n);  // pivot vector

  atlas::getrf (a, ipiv);      // factor a
  atlas::getrs (a, ipiv, b);   // solve from factorization 
  print_m (b, "X"); 
  cout << endl; 

  print_v (ipiv, "pivots"); 

  cout << endl; 

}

