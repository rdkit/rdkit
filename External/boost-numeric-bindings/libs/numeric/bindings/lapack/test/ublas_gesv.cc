
// solving A * X = B
// using driver function gesv()

#include <cstddef>
#include <iostream>
#include <complex>
#include <boost/numeric/bindings/lapack/gesv.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/std_vector.hpp>
#include "utils.h"

namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;

using std::size_t; 
using std::cout;
using std::endl; 

typedef std::complex<double> cmpx_t; 

typedef ublas::matrix<double, ublas::column_major> m_t;
typedef ublas::matrix<cmpx_t, ublas::column_major> cm_t;

int main() {

  cout << endl; 
  cout << "real system:" << endl << endl; 

  size_t n = 5;   
  m_t a (n, n);   // system matrix 

  size_t nrhs = 2; 
  m_t x (n, nrhs), b (n, nrhs);  // b -- right-hand side matrix

  std::vector<int> ipiv (n);  // pivot vector

  init_symm (a); 
  //     [n   n-1 n-2  ... 1]
  //     [n-1 n   n-1  ... 2]
  // a = [n-2 n-1 n    ... 3]
  //     [        ...       ]
  //     [1   2   ...  n-1 n]

  m_t aa (a);  // copy of a, because a is `lost' after gesv()

  for (int i = 0; i < x.size1(); ++i) {
    x (i, 0) = 1.;
    x (i, 1) = 2.; 
  }
  b = prod (a, x); 

  print_m (a, "A"); 
  cout << endl; 
  print_m (b, "B"); 
  cout << endl; 

  lapack::gesv (a, ipiv, b);   // solving the system, b contains x 
  print_m (b, "X"); 
  cout << endl; 

  x = prod (aa, b); 
  print_m (x, "B = A X"); 
  cout << endl; 

  ////////////////////////////////////////////////////////

  cout << endl; 
  cout << "complex system:" << endl << endl; 
  cm_t ca (3, 3), cb (3, 1), cx (3, 1);
  std::vector<int> ipiv2 (3); 

  ca (0, 0) = cmpx_t (3, 0);
  ca (0, 1) = cmpx_t (4, 2);
  ca (0, 2) = cmpx_t (-7, 5);
  ca (1, 0) = cmpx_t (4, -2);
  ca (1, 1) = cmpx_t (-5, 0);
  ca (1, 2) = cmpx_t (0, -3);
  ca (2, 0) = cmpx_t (-7, -5);
  ca (2, 1) = cmpx_t (0, 3);
  ca (2, 2) = cmpx_t (2, 0);
  print_m (ca, "CA"); 
  cout << endl; 

  for (int i = 0; i < cx.size1(); ++i) 
    cx (i, 0) = cmpx_t (1, -1); 
  cb = prod (ca, cx); 
  print_m (cb, "CB"); 
  cout << endl; 
  
  int ierr = lapack::gesv (ca, ipiv2, cb); 
  if (ierr == 0) 
    print_m (cb, "CX");
  else
    cout << "matrix is singular" << endl; 

  cout << endl; 

}

