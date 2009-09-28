
// solving A * X = B
// using driver function gesv()

//#define BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 
//#define BOOST_NO_FUNCTION_TEMPLATE_ORDERING

#include <cstddef>
#include <iostream>
#include <boost/numeric/bindings/atlas/cblas.hpp>
#include <boost/numeric/bindings/atlas/clapack.hpp>
#include <boost/numeric/bindings/traits/tnt.hpp>
#ifndef F_FORTRAN 
#  include <tnt/tnt_array2d_utils.h>
#else
#  include <tnt/tnt_fortran_array2d_utils.h>
#endif 

namespace atlas = boost::numeric::bindings::atlas;

using std::size_t; 
using std::cout;
using std::endl; 

#ifndef F_FORTRAN 
typedef TNT::Array2D<double> m_t;
#else
typedef TNT::Fortran_Array2D<double> m_t;
#endif

int main() {

  cout << endl; 
  size_t n = 3, nrhs = 1; 

  m_t a (n, n);   // system matrix 
#ifndef F_FORTRAN 
  a[0][0] = 1.; a[0][1] = 1.; a[0][2] = 1.;
  a[1][0] = 2.; a[1][1] = 3.; a[1][2] = 1.;
  a[2][0] = 1.; a[2][1] = -1.; a[2][2] = -1.;
#else
  a(1,1) = 1.; a(1,2) = 1.; a(1,3) = 1.;
  a(2,1) = 2.; a(2,2) = 3.; a(2,3) = 1.;
  a(3,1) = 1.; a(3,2) = -1.; a(3,3) = -1.;
#endif 

// see leading comments for `gesv()' in clapack.hpp
#ifndef F_FORTRAN 
  m_t b (nrhs, n);  // right-hand side matrix
  b[0][0] = 4.; b[0][1] = 9.; b[0][2] = -2.; 
#else
  m_t b (n, nrhs);  
  b(1,1) = 4.; b(2,1) = 9.; b(3,1) = -2.; 
#endif 

  cout << "A: " << a << endl; 
  cout << "B: " << b << endl; 

  atlas::lu_solve (a, b);  
  cout << "X: " << b << endl; 

  cout << endl; 
}

