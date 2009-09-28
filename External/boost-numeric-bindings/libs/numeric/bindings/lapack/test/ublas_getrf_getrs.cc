
// solving A * X = B
// in two steps -- factor (getrf()) and solve (getrs())

#include <cstddef>
#include <iostream>
#include <complex>
#include <boost/numeric/bindings/lapack/gesv.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/std_vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "utils.h"

namespace ublas = boost::numeric::ublas;
namespace traits = boost::numeric::bindings::traits;
namespace lapack = boost::numeric::bindings::lapack;

using std::size_t; 
using std::cin;
using std::cout;
using std::endl; 

typedef std::complex<double> cmpx; 
typedef ublas::matrix<double, ublas::column_major> m_t;
typedef ublas::matrix<cmpx, ublas::column_major> cm_t;


int main (int argc, char **argv) {
  size_t n = 0;
  if (argc > 1) {
    n = atoi(argv [1]);
  }

  cout << endl; 

  if (n <= 0) {
  cout << "n -> ";
  cin >> n; 
  }
  if (n < 5) {
    n = 5;
    cout << "min n = 5" << endl; 
  }
  cout << endl; 
  m_t a (n, n);   // system matrix 

  size_t nrhs = 2; 
  m_t x (n, nrhs), b (n, nrhs);  // b -- right-hand side matrix

  init_symm (a); 
  //     [n   n-1 n-2  ... 1]
  //     [n-1 n   n-1  ... 2]
  // a = [n-2 n-1 n    ... 3]
  //     [        ...       ]
  //     [1   2   ...  n-1 n]

  for (int i = 0; i < x.size1(); ++i) {
    x (i, 0) = 1.;
    x (i, 1) = 2.; 
  }
  b = prod (a, x); 
  m_t a2 (a);  // for part 2
  m_t b2 (b);

  print_m (a, "A"); 
  cout << endl; 
  print_m (b, "B"); 
  cout << endl; 

  ublas::matrix_row<m_t> ar1 (a, 1), ar3 (a, 3);
  ublas::matrix_row<m_t> br1 (b, 1), br3 (b, 3);
  swap (ar1, ar3);   // swap rows to force pivoting 
  swap (br1, br3);
  print_m (a, "A");  // print `new' system  
  cout << endl; 
  print_m (b, "B"); 
  cout << endl; 

  std::vector<int> ipiv (n);  // pivot vector

  lapack::getrf (a, ipiv);      // factor a
  m_t ia (a);
  lapack::getrs (a, ipiv, b);   // solve from factorization 
  print_m (b, "X"); 
  cout << endl; 
  lapack::getri (ia, ipiv);     // invert a
  print_m (ia, "InvA"); 
  cout << endl; 

  print_v (ipiv, "pivots"); 

  cout << endl; 

  ublas::matrix_column<m_t> a2c1 (a2, 1), a2c4 (a2, 4);
  ublas::matrix_row<m_t> b2r1 (b2, 1), b2r4 (b2, 4);
  swap (a2c1, a2c4);   // swap columns
  swap (b2r1, b2r4);
  print_m (a2, "A");  // print `new' system  
  cout << endl; 
  print_m (b2, "B"); 
  cout << endl; 
  
  lapack::getrf (a2, ipiv); // factor a
  m_t ia2 (a2);
  lapack::getrs ('T', a2, ipiv, b2); // solve 
  print_m (b2, "X"); 
  cout << endl; 
  lapack::getri (ia2, ipiv); // invert a2
  print_m (ia2, "InvA2"); 
  cout << endl; 

  print_v (ipiv, "pivots"); 

  cout << endl; 
}

