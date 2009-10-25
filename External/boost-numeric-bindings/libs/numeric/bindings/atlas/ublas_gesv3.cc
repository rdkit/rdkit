
// solving A * X = B
// using driver function gesv()

//#define BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 
//#define BOOST_NO_FUNCTION_TEMPLATE_ORDERING

#include <cstddef>
#include <iostream>
#include <boost/numeric/bindings/atlas/cblas.hpp>
#include <boost/numeric/bindings/atlas/clapack.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/ublas/io.hpp> 

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
  size_t n = 3, nrhs = 1; 

  m_t a (n, n);   // system matrix 
  a(0,0) = 1.; a(0,1) = 1.; a(0,2) = 1.;
  a(1,0) = 2.; a(1,1) = 3.; a(1,2) = 1.;
  a(2,0) = 1.; a(2,1) = -1.; a(2,2) = -1.;

// see leading comments for `gesv()' in clapack.hpp
#ifndef F_ROW_MAJOR
  m_t b (n, nrhs);  // right-hand side matrix
  b(0,0) = 4.; b(1,0) = 9.; b(2,0) = -2.; 
#else
  m_t b (nrhs, n);  
  b(0,0) = 4.; b(0,1) = 9.; b(0,2) = -2.; 
#endif 

  cout << "A: " << a << endl; 
  cout << "B: " << b << endl; 

  atlas::lu_solve (a, b);  
  cout << "X: " << b << endl; 

  cout << endl; 
}

