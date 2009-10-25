
//#define BOOST_NUMERIC_BINDINGS_LAPACK_2

#include <cstddef>
#include <iostream>
#include <algorithm> 
#include <complex>
#include <boost/numeric/bindings/lapack/gesvd.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/std_vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "utils.h"

namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;

using std::size_t; 
using std::cout;
using std::endl; 

typedef double real_t; 
typedef std::complex<real_t> cmplx_t; 
typedef ublas::matrix<cmplx_t, ublas::column_major> m_t;
typedef ublas::vector<cmplx_t> v_t;
typedef ublas::matrix<real_t, ublas::column_major> rm_t;
typedef ublas::vector<real_t> rv_t;

int main() {

  cout << endl; 

  size_t m = 1, n = 3;   
  size_t minmn = m < n ? m : n; 
  m_t a (m, n);  
  a(0,0) = cmplx_t (2, -1);
  a(0,1) = cmplx_t (1, 1);
  a(0,2) = cmplx_t (-2, 0);

  m_t a2 (a); // for part 2
  m_t a3 (a); // for part 3

  print_m (a, "A"); 
  cout << endl; 

  rv_t s (minmn); // singular values are real 
  m_t u (m, m);
  m_t vt (n, n);

  size_t lw; 

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_2
  lw = lapack::gesvd_work ('O', 'N', 'N', a); 
  cout << "opt NN lw: " << lw << endl; 
  lw = lapack::gesvd_work ('O', 'A', 'A', a); 
  cout << "opt AA lw: " << lw << endl; 
  lw = lapack::gesvd_work ('O', 'S', 'S', a); 
  cout << "opt SS lw: " << lw << endl; 
  lw = lapack::gesvd_work ('O', 'O', 'N', a); 
  cout << "opt ON lw: " << lw << endl; 
  lw = lapack::gesvd_work ('O', 'N', 'O', a); 
  cout << "opt NO lw: " << lw << endl; 
#endif 
  lw = lapack::gesvd_work ('M', 'A', 'A', a); 
  cout << "min lw: " << lw << endl << endl; 

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_2
  lw = lapack::gesvd_work ('O', 'A', 'A', a); 
#endif 

  std::vector<cmplx_t> w (lw); 

  size_t lrw = lapack::gesvd_rwork (a);
  cout << "lrw: " << lrw << endl << endl; 
  std::vector<real_t> rw (lrw);

  lapack::gesvd ('A', 'A', a, s, u, vt, w, rw);  

  print_v (s, "s"); 
  cout << endl; 
  print_m (u, "U"); 
  cout << endl; 
  print_m (vt, "V^T"); 
  cout << endl; 

  rm_t sm (m, n); 
  for (size_t i = 0; i < s.size(); ++i) 
    sm (i,i) = s (i); 
  print_m (sm, "S"); 
  cout << endl;

  a = ublas::prod (u, m_t (ublas::prod (sm, vt))); 
  print_m (a, "A == U S V^T"); 
  cout << endl; 

  // part 2 

  cout << endl << "part 2" << endl << endl; 
 
#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_2
  lapack::gesvd ('A', 'A', a2, s, u, vt);  
#else
  lapack::gesvd ('M', 'A', 'A', a2, s, u, vt);  
#endif

  print_v (s, "s"); 
  cout << endl; 
  print_m (u, "U"); 
  cout << endl; 
  print_m (vt, "V^T"); 
  cout << endl; 

  for (size_t i = 0; i < s.size(); ++i) 
    sm (i,i) = s (i); 
  print_m (sm, "S"); 
  cout << endl;

  a2 = ublas::prod (u, m_t (ublas::prod (sm, vt))); 
  print_m (a2, "A == U S V^T"); 
  cout << endl; 
 
  // part 3

  cout << endl << "part 3" << endl << endl; 

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_2
  cout << "opt lw: " << lapack::gesvd_work ('O', 'N', 'N', a3) 
       << endl << endl; 
  lapack::gesvd (a3, s);
#else 
  cout << "min lw: " << lapack::gesvd_work ('M', 'N', 'N', a3) 
       << endl << endl; 
  lapack::gesvd ('M', 'N', 'N', a3, s, u, vt);
#endif 

  print_v (s, "singular values only"); 
  cout << endl; 

  cout << endl; 
}

