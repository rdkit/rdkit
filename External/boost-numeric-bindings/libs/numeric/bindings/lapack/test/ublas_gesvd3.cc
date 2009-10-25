
//#define BOOST_NUMERIC_BINDINGS_LAPACK_2

#include <cstddef>
#include <iostream>
#include <algorithm> 
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
typedef ublas::matrix<real_t, ublas::column_major> m_t;
typedef ublas::vector<real_t> v_t;

int main() {

  cout << endl; 

  size_t m = 3, n = 4;   
  size_t minmn = m < n ? m : n; 
  m_t a (m, n);  
  a(0,0) = 2.;  a(0,1) = 2.;  a(0,2) = 2.;   a(0,3) = 2.;
  a(1,0) = 1.7; a(1,1) = 0.1; a(1,2) = -1.7; a(1,3) = -0.1;
  a(2,0) = 0.6; a(2,1) = 1.8; a(2,2) = -0.6; a(2,3) = -1.8;

  m_t a2 (a); // for part 2
  m_t a3 (a); // for parts 3 & 4

  print_m (a, "A"); 
  cout << endl; 

  m_t sm (m, n); 
  ublas::matrix_vector_range<m_t> 
    smr (sm, ublas::range (0, minmn), ublas::range (0, minmn)); 
  v_t s (minmn); 
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

  std::vector<real_t> w (lw); 

  lapack::gesvd ('A', 'A', a, s, u, vt, w);  

  print_v (s, "s"); 
  cout << endl; 
  smr = s; 
  print_m (sm, "S"); 
  cout << endl;
  print_m (u, "U"); 
  cout << endl; 
  print_m (vt, "V^T"); 
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
  smr = s; 
  print_m (sm, "S"); 
  cout << endl;
  print_m (u, "U"); 
  cout << endl; 
  print_m (vt, "V^T"); 
  cout << endl; 

  a2 = ublas::prod (u, m_t (ublas::prod (sm, vt))); 
  print_m (a2, "A == U S V^T"); 
  cout << endl; 

  // part 3 -- another check
  sm = ublas::prod (ublas::trans (u), 
                    m_t (ublas::prod (a3, ublas::trans (vt)))); 
  print_m (sm, "S == U^T A V"); 
  cout << endl << endl;

  // part 4 
  // singular values only

  cout << endl << "part 4" << endl << endl; 

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_2
  cout << "opt lw: " << lapack::gesvd_work ('O', 'N', 'N', a3) 
       << endl << endl; 
  lapack::gesvd ('N', 'N', a3, s, u, vt);
#else 
  cout << "min lw: " << lapack::gesvd_work ('M', 'N', 'N', a3) 
       << endl << endl; 
  lapack::gesvd ('M', 'N', 'N', a3, s, u, vt);
#endif 

  print_v (s, "singular values only"); 
  cout << endl; 

  cout << endl; 
}

