
// inverting A
// using getrf() & getri() 

//#define BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 
//#define BOOST_NO_FUNCTION_TEMPLATE_ORDERING

#include <cstddef>
#include <iostream>
#include <complex>
#include <boost/numeric/bindings/atlas/cblas3.hpp>
#include <boost/numeric/bindings/atlas/clapack.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/std_vector.hpp>
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
  cout << "real matrix:" << endl << endl; 

  size_t n = 5; 
  m_t a (n, n);
  init_symm (a); 
  //     [n   n-1 n-2  ... 1]
  //     [n-1 n   n-1  ... 2]
  // a = [n-2 n-1 n    ... 3]
  //     [        ...       ]
  //     [1   2   ...  n-1 n]
  print_m (a, "A"); 
  cout << endl; 

  m_t aa (a);  // copy of a, for later use

  std::vector<int> ipiv (n);   // pivot vector 
  atlas::lu_factor (a, ipiv);  // alias for getrf()
  atlas::lu_invert (a, ipiv);  // alias for getri() 

  m_t i1 (n, n), i2 (n, n);  
  atlas::gemm (a, aa, i1);   // i1 should be (almost) identity matrix
  atlas::gemm (aa, a, i2);   // i2 should be (almost) identity matrix

  print_m (i1, "I = A^(-1) * A");  
  cout << endl; 
  print_m (i2, "I = A * A^(-1)"); 
  cout << endl; 
  
  cout << endl; 

  //////////////////////////////////////////////////////

  cout << "complex matrix:" << endl << endl; 
  cm_t ca (3, 3); 

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
  
  std::vector<int> ipiv2 (3); 
  
  int ierr = atlas::getrf (ca, ipiv2);
  if (ierr == 0) {
    atlas::getri (ca, ipiv2); 
    cm_t ii (3, 3); 
    atlas::gemm (ca, caa, ii);
    print_m (ii, "I = CA^(-1) * CA"); 
    cout << endl; 
    atlas::gemm (caa, ca, ii);
    print_m (ii, "I = CA * CA^(-1)"); 
    cout << endl; 
  }
  else
    cout << "matrix is singular" << endl; 
  
  cout << endl; 

}

