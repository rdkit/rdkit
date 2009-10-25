
// nothing to do with ATLAS bindings ;o) 
// various initializations from `utils.h' 

#include <iostream>
#include <complex> 
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include "utils.h"

#ifndef F_FLOAT
typedef double real_t;
#else
typedef float real_t;
#endif 

#ifndef F_COMPLEX
typedef real_t elem_t;
#else
typedef std::complex<real_t> elem_t;
#endif 

typedef boost::numeric::ublas::vector<elem_t> vct_t; 
#ifdef F_ROW_MAJOR
typedef boost::numeric::ublas::matrix<elem_t> matr_t; 
#else
typedef boost::numeric::ublas::matrix<
  elem_t,
  boost::numeric::ublas::column_major
> matr_t; 
#endif 

int main() {
  std::cout << std::endl; 

  vct_t v (10); 

  // v[i] = i || (i, 0)
  init_v<ident> (v); 
  print_v (v); 
  std::cout << std::endl; 

  // v[i] = i+1 || (i+1, 0)
  init_v (v, iplus1()); 
  print_v (v); 
  std::cout << std::endl; 

  // v[i] = k || (k, 0) ; i = 0, k = 5 ;  ++i, ++k
  init_v (v, kpp (5)); 
  print_v (v); 
  std::cout << std::endl; 

  // v[i] = 2 || (2, 0)
  init_v (v, const_val<elem_t> (2)); 
  print_v (v); 
  std::cout << std::endl; 

  // v[i] = 10 i || (10 i, 0)
  init_v (v, times_plus<elem_t> (10)); 
  print_v (v); 
  std::cout << std::endl; 

  // v[i] = i+0.1 || (i+0.1, 0)
  init_v (v, times_plus<elem_t> (1, 0.1)); 
  print_v (v); 
  std::cout << std::endl; 
  std::cout << std::endl; 

  matr_t m (5, 7); 

  // m[i, j] = k ; i = 0, j = 0, k = 0 ; ++j, ++k
  init_m<kpp> (m); 
  print_m (m); 
  std::cout << std::endl; 

#ifndef F_COMPLEX
  elem_t s (10); 
#else
  elem_t s (10, -1); 
#endif 

  // m[i, j] = s (i+1) + (j+0.1)
  init_m (m, times_plus<elem_t> (s, 1, 0.1)); 
  print_m (m); 
  std::cout << std::endl; 

#ifndef F_COMPLEX
  elem_t c (-2); 
#else
  elem_t c (2, -2); 
#endif 

  // m[i, j] = c
  init_m (m, const_val<elem_t> (c)); 
  print_m (m); 
  std::cout << std::endl; 

  // m[i, j] = i
  init_m<rws> (m); 
  print_m (m); 
  std::cout << std::endl; 

  // m[i, j] = j+1
  init_m (m, cls1()); 
  print_m (m); 
  std::cout << std::endl; 

  matr_t sm (4, 4); 

  // [4 3 2 1]
  // [0 4 3 2]
  // [0 0 4 3]
  // [0 0 0 4]
  init_symm (sm, 'u'); 
  print_m (sm); 
  std::cout << std::endl; 
  print_m_data (sm);
  std::cout << std::endl; 

  init_m (sm, const_val<elem_t> (0)); 
  // [4 3 2 1]
  // [3 4 3 2]
  // [2 3 4 3]
  // [1 2 3 4]
  init_symm (sm); 
  print_m (sm); 
  std::cout << std::endl; 
  print_m_data (sm);
  std::cout << std::endl; 

  init_m (sm, const_val<elem_t> (0)); 
  // [4 0 0 0]
  // [3 4 0 0]
  // [2 3 4 0]
  // [1 2 3 4]
  init_symm (sm, 'L'); 
  print_m (sm); 
  std::cout << std::endl; 
  print_m_data (sm);
  std::cout << std::endl; 

} 
