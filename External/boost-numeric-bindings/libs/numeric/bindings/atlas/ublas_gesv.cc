
// solving A * X = B
// using driver function gesv()

//#define BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 
//#define BOOST_NO_FUNCTION_TEMPLATE_ORDERING

#include <cstddef>
#include <iostream>
#include <complex>
#include <boost/numeric/bindings/atlas/cblas.hpp>
#include <boost/numeric/bindings/atlas/clapack.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/std_vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "utils.h"

namespace ublas = boost::numeric::ublas;
namespace atlas = boost::numeric::bindings::atlas;

using std::size_t; 
using std::cout;
using std::endl; 

typedef std::complex<double> cmpx; 

#ifndef F_ROW_MAJOR
typedef ublas::matrix<double, ublas::column_major> m_t;
typedef ublas::matrix<cmpx, ublas::column_major> cm_t;
#else
typedef ublas::matrix<double, ublas::row_major> m_t;
typedef ublas::matrix<cmpx, ublas::row_major> cm_t;
#endif

int main() {

  cout << endl; 
  cout << "real system:" << endl << endl; 

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

  std::vector<int> ipiv (n);  // pivot vector

  init_symm (a); 
  //     [n   n-1 n-2  ... 1]
  //     [n-1 n   n-1  ... 2]
  // a = [n-2 n-1 n    ... 3]
  //     [        ...       ]
  //     [1   2   ...  n-1 n]

  m_t aa (a);  // copy of a, because a is `lost' after gesv()

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

  atlas::gesv (a, ipiv, b);   // solving the system, b contains x 
  print_m (b, "X"); 
  cout << endl; 

#ifndef F_ROW_MAJOR
  atlas::gemm (aa, b, x);     // check the solution 
#else
  atlas::gemv (aa, br0, xc0);
  atlas::gemv (aa, br1, xc1);
#endif 
  print_m (x, "B = A X"); 
  cout << endl; 

  ////////////////////////////////////////////////////////

  cout << endl; 
  cout << "complex system:" << endl << endl; 
  cm_t ca (3, 3);
  cm_t cx (3, 1);
#ifndef F_ROW_MAJOR
  cm_t cb (3, 1);
#else
  cm_t cb (1, 3); 
#endif  

  ca (0, 0) = cmpx (3, 0);
  ca (0, 1) = cmpx (4, 2);
  ca (0, 2) = cmpx (-7, 5);
  ca (1, 0) = cmpx (4, -2);
  ca (1, 1) = cmpx (-5, 0);
  ca (1, 2) = cmpx (0, -3);
  ca (2, 0) = cmpx (-7, -5);
  ca (2, 1) = cmpx (0, 3);
  ca (2, 2) = cmpx (2, 0);
  print_m (ca, "CA"); 
  cout << endl; 
  cm_t caa (ca); 
  
  ublas::matrix_column<cm_t> cx0 (cx, 0);
  atlas::set (cmpx (1, -1), cx0);
#ifndef F_ROW_MAJOR
  ublas::matrix_column<cm_t> cb0 (cb, 0); 
#else
  ublas::matrix_row<cm_t> cb0 (cb, 0); 
#endif
  atlas::gemv (ca, cx0, cb0); 
  print_m (cb, "CB"); 
  cout << endl; 
  
  int ierr = atlas::gesv (ca, cb); // with `internal' pivot vector
  if (ierr == 0) {
    print_m (cb, "CX");
    cout << endl; 
    atlas::gemv (caa, cb0, cx0);
    print_m (cx, "CB");
  }
  else
    cout << "matrix is singular" << endl; 

  cout << endl; 
}

