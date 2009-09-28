
// BLAS level 3 -- complex numbers

//#define BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 
//#define BOOST_NO_FUNCTION_TEMPLATE_ORDERING

#include <iostream>
#include <complex>
#include <boost/numeric/bindings/atlas/cblas3.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
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
typedef ublas::matrix<cmplx_t, ublas::row_major> m_t;
#else
typedef ublas::matrix<cmplx_t, ublas::column_major, std::vector<cmplx_t> > m_t;
#endif 

int main() {

  cout << endl; 

  m_t a (2, 2);
  a (0, 0) = cmplx_t (1., 0.);
  a (0, 1) = cmplx_t (2., 0.);
  a (1, 0) = cmplx_t (3., 0.);
  a (1, 1) = cmplx_t (4., 0.);
  print_m (a, "A"); 
  cout << endl; 

  m_t b (2, 3);
  b (0, 0) = cmplx_t (1., 0.);
  b (0, 1) = cmplx_t (2., 0.);
  b (0, 2) = cmplx_t (3., 0.);
  b (1, 0) = cmplx_t (1., 0.);
  b (1, 1) = cmplx_t (2., 0.);
  b (1, 2) = cmplx_t (3., 0.);
  print_m (b, "B"); 
  cout << endl; 
  
  m_t c (2, 3);

  // c = a b
  atlas::gemm (a, b, c); 
  print_m (c, "A B"); 
  cout << endl; 

  a (0, 0) = cmplx_t (0., 1.);
  a (0, 1) = cmplx_t (0., 2.);
  a (1, 0) = cmplx_t (0., 3.);
  a (1, 1) = cmplx_t (0., 4.);
  print_m (a, "A"); 
  cout << endl; 
  
  // c = a b
  atlas::gemm (CblasNoTrans, CblasNoTrans, 1.0, a, b, 0.0, c); 
  print_m (c, "A B"); 
  cout << endl; 

  // c = a^T b
  atlas::gemm (CblasTrans, CblasNoTrans, 1.0, a, b, 0.0, c); 
  print_m (c, "A^T B"); 
  cout << endl; 

  // c = a^H b
  atlas::gemm (CblasConjTrans, CblasNoTrans, 1.0, a, b, 0.0, c); 
  print_m (c, "A^H B"); 

  cout << endl; 

}
