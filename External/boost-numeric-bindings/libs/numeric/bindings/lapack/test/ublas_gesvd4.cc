
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

  m_t a2 (a); // for parts 2 & 3

  print_m (a, "A"); 
  cout << endl; 

  v_t s (minmn); 
  m_t u (m, minmn);
  m_t vt (minmn, n);

  lapack::gesvd (a, s, u, vt);  

  print_v (s, "s"); 
  cout << endl; 
  print_m (u, "U"); 
  cout << endl; 
  print_m (vt, "V^T"); 
  cout << endl << endl;

  // part 2

  // singular values and singular vectors satisfy  A v_i == s_i u_i
  for (int i = 0; i < minmn; ++i) {
    cout << "A v_" << i << " == s_" << i << " u_" << i << endl; 
    v_t avi = ublas::prod (a2, ublas::row (vt, i));  
    print_v (avi);
    v_t siui = s[i] * ublas::column (u, i);
    print_v (siui);
    cout << endl; 
  }
  cout << endl; 
  
  // singular values and singular vectors satisfy  A^T u_i == s_i v_i
  for (int i = 0; i < minmn; ++i) {
    cout << "A^T u_" << i << " == s_" << i << " v_" << i << endl; 
    v_t atui = ublas::prod (trans (a2), ublas::column (u, i));  
    print_v (atui);
    v_t sivi = s[i] * ublas::row (vt, i);
    print_v (sivi);
    cout << endl; 
  }
  cout << endl; 
  
  // part 3 
  
  lapack::gesvd (a2, s);  
  
  print_v (s, "singular values only"); 
  cout << endl; 
}

