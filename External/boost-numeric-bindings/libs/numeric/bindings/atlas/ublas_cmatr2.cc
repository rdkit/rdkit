
// BLAS level 2 -- complex numbers

//#define BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 
//#define BOOST_NO_FUNCTION_TEMPLATE_ORDERING

#include <iostream>
#include <complex>
#include <boost/numeric/bindings/atlas/cblas1.hpp>
#include <boost/numeric/bindings/atlas/cblas2.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#ifdef F_USE_STD_VECTOR
#include <vector>
#include <boost/numeric/bindings/traits/std_vector.hpp> 
#endif 
#include "utils.h" 

namespace ublas = boost::numeric::ublas;
namespace atlas = boost::numeric::bindings::atlas;

using std::cout;
using std::endl; 

typedef double real_t;
typedef std::complex<real_t> cmplx_t; 

#ifndef F_USE_STD_VECTOR
typedef ublas::vector<cmplx_t> vct_t;
typedef ublas::matrix<cmplx_t, ublas::row_major> m_t;
#else
typedef ublas::vector<cmplx_t, std::vector<cmplx_t> > vct_t;
typedef ublas::matrix<cmplx_t, ublas::column_major, std::vector<cmplx_t> > m_t;
#endif 

int main() {

  cout << endl; 

  vct_t vx (2);
  atlas::set (cmplx_t (1., 0.), vx);
  print_v (vx, "vx"); 
  vct_t vy (4); // vector size can be larger 
                // than corresponding matrix size 
  atlas::set (cmplx_t (0., 0.), vy); 
  print_v (vy, "vy"); 
  cout << endl; 

  m_t m (3, 2);
  init_m (m, kpp (1)); 
  print_m (m, "m"); 
  cout << endl; 

  // vy = m vx
  atlas::gemv (CblasNoTrans, 1.0, m, vx, 0.0, vy);
  print_v (vy, "m vx"); 
  atlas::gemv (m, vx, vy);
  print_v (vy, "m vx"); 
  cout << endl; 

  m (0, 0) = cmplx_t (0., 1.);
  m (0, 1) = cmplx_t (0., 2.);
  m (1, 0) = cmplx_t (0., 3.);
  m (1, 1) = cmplx_t (0., 4.);
  m (2, 0) = cmplx_t (0., 5.);
  m (2, 1) = cmplx_t (0., 6.);
  print_m (m, "m"); 
  cout << endl; 

  // vy = m vx
  atlas::gemv (CblasNoTrans, 1.0, m, vx, 0.0, vy);
  print_v (vy, "m vx"); 
  atlas::gemv (m, vx, vy);
  print_v (vy, "m vx"); 
  cout << endl; 

  m (0, 0) = cmplx_t (-1., 1.);
  m (0, 1) = cmplx_t (-2., 2.);
  m (1, 0) = cmplx_t (-3., 3.);
  m (1, 1) = cmplx_t (-4., 4.);
  m (2, 0) = cmplx_t (-5., 5.);
  m (2, 1) = cmplx_t (-6., 6.);
  print_m (m, "m"); 
  cout << endl; 

  // vy = m vx
  atlas::gemv (CblasNoTrans, 1.0, m, vx, 0.0, vy);
  print_v (vy, "m vx"); 
  atlas::gemv (m, vx, vy);
  print_v (vy, "m vx"); 
  cout << endl; 

  atlas::set (cmplx_t (0., 1.), vx);
  print_v (vx, "vx"); 

  // vy = m vx
  atlas::gemv (CblasNoTrans, 1.0, m, vx, 0.0, vy);
  print_v (vy, "m vx"); 
  atlas::gemv (m, vx, vy);
  print_v (vy, "m vx"); 
  cout << endl; 

  atlas::set (cmplx_t (1., 1.), vx);
  print_v (vx, "vx"); 

  // vy = m vx
  atlas::gemv (CblasNoTrans, 1.0, m, vx, 0.0, vy);
  print_v (vy, "m vx"); 
  atlas::gemv (m, vx, vy);
  print_v (vy, "m vx"); 
  cout << endl; 

  // vx = m^H vy
  atlas::set (cmplx_t (-1., -1.), vy); 
  print_v (vy, "vy"); 
  atlas::gemv (CblasConjTrans, 1.0, m, vy, 0.0, vx);
  print_v (vx, "m^H vy"); 
  cout << endl; 

  m_t mx (2, 2); 
  m_t my (3, 2); 

  ublas::matrix_column<m_t> mxc0 (mx, 0), mxc1 (mx, 1); 
  ublas::matrix_column<m_t> myc0 (my, 0), myc1 (my, 1); 

  atlas::set (cmplx_t (1., 0.), mxc0);
  atlas::set (cmplx_t (0., 0.), mxc1);
  atlas::set (cmplx_t (0., 0.), myc0);
  atlas::set (cmplx_t (0., 0.), myc1);
  print_m (mx, "mx");
  cout << endl; 
  print_m (my, "my");
  cout << endl; 

  // my[.,0] = m mx[.,0] 
  atlas::gemv (m, mxc0, myc0); 
  print_m (my, "m mx[.,0]");

  cout << endl;

}
