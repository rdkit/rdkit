
// BLAS level 3
// TNT arrays

#include <iostream>
#include <boost/numeric/bindings/atlas/cblas3.hpp>
#include <boost/numeric/bindings/traits/tnt.hpp>
#include "utils.h"
#include "tnt_utils.h"

namespace atlas = boost::numeric::bindings::atlas;

using std::cout;
using std::endl; 

#ifndef F_FORTRAN
typedef TNT::Array2D<double> m_t;
#else
typedef TNT::Fortran_Array2D<double> m_t;
#endif 

int main() {

  cout << endl; 

  m_t a (2, 2);
  init_m (a, kpp (1)); 
  print_m (a, "A"); 
  cout << endl; 

  m_t b (2, 3);
  init_m (b, cls1()); 
  print_m (b, "B"); 
  cout << endl; 
  
  m_t c (2, 3);
  init_m<const_val<double> > (c); 
  print_m (c, "C"); 
  cout << endl; 

  atlas::gemm (a, b, c); 
  print_m (c, "A B"); 
  cout << endl; 

  atlas::gemm (CblasNoTrans, CblasNoTrans, 1.0, a, b, 0.0, c); 
  print_m (c, "A B"); 
  cout << endl; 

  atlas::gemm (1.0, a, b, 0.0, c); 
  print_m (c, "A B"); 
  cout << endl; 

  init_m (c, const_val<double> (1.)); 
  print_m (c, "C"); 
  cout << endl; 

  atlas::gemm (1.0, a, b, 0.05, c);
  print_m (c, "0.05 C + A B"); 
  cout << endl; 

  m_t d (3, 2);

  atlas::gemm (CblasTrans, CblasTrans, 1.0, b, a, 0.0, d);
  print_m (d, "B^T A^T"); 
  cout << endl; 

  atlas::gemm (CblasTrans, CblasNoTrans, 1.0, a, b, 0.0, c); 
  print_m (c, "A^T B"); 
  cout << endl; 

  atlas::gemm (CblasTrans, CblasNoTrans, 1.0, b, a, 0.0, d);
  print_m (d, "B^T A"); 

  cout << endl; 

}
