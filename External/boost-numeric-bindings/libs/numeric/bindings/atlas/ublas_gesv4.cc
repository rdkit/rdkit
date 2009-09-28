
// solving A * X = B
// using driver function gesv()
// with ublas::vector<> as RHS

//#define BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 
//#define BOOST_NO_FUNCTION_TEMPLATE_ORDERING

//#define BOOST_NUMERIC_BINDINGS_NO_SANITY_CHECK
//#define BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK

#include <cstddef>
#include <iostream>
#include <boost/numeric/bindings/atlas/cblas.hpp>
#include <boost/numeric/bindings/atlas/clapack.hpp>
#include <boost/numeric/bindings/traits/ublas_vector2.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/ublas/io.hpp> 

namespace ublas = boost::numeric::ublas;
namespace atlas = boost::numeric::bindings::atlas;

using std::size_t; 
using std::cout;
using std::endl; 

typedef ublas::matrix<double, ublas::column_major> m_t;
typedef ublas::vector<double> v_t; 

int main() {

  cout << endl; 
  size_t n = 3;

  m_t a (n, n);   // system matrix 
  a(0,0) = 1.; a(0,1) = 1.; a(0,2) = 1.;
  a(1,0) = 2.; a(1,1) = 3.; a(1,2) = 1.;
  a(2,0) = 1.; a(2,1) = -1.; a(2,2) = -1.;

  v_t b (n);  // right-hand side vector
  b(0) = 4.; b(1) = 9.; b(2) = -2.; 

  cout << "A: " << a << endl; 
  cout << "B: " << b << endl; 

  atlas::lu_solve (a, b);  
  cout << "X: " << b << endl; 
  cout << endl; 

}

