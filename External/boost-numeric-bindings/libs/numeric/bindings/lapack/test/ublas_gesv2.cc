
// solving A * X = B
// using driver function gesv()

#include <cstddef>
#include <iostream>
#include <vector>
#include <boost/numeric/bindings/lapack/gesv.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/std_vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "utils.h"

namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;

using std::size_t; 
using std::cout;
using std::endl; 

typedef ublas::matrix<double, ublas::column_major> m_t;

int main() {

  cout << endl; 

  size_t n = 5;   
  m_t a (n, n);   // system matrix 

  size_t nrhs = 2; 
  m_t x (n, nrhs), b (n, nrhs);  // b -- right-hand side matrix

  init_symm (a); 
  //     [n   n-1 n-2  ... 1]
  //     [n-1 n   n-1  ... 2]
  // a = [n-2 n-1 n    ... 3]
  //     [        ...       ]
  //     [1   2   ...  n-1 n]

  m_t aa (a); // copy of a, because a is `lost' after gesv()

  ublas::matrix_column<m_t> xc0 (x, 0), xc1 (x, 1); 
  for (int i = 0; i < xc0.size(); ++i) {
    xc0 (i) = 1.;
    xc1 (i) = 2.; 
  }
  b = prod (a, x); 

  print_m (a, "A"); 
  cout << endl; 
  print_m (b, "B"); 
  cout << endl; 

  lapack::gesv (a, b);  // solving the system, b contains x 

  print_m (b, "X");
  cout << endl; 

  x = prod (aa, b); 
  print_m (x, "B = A X"); 

  cout << endl; 

}

