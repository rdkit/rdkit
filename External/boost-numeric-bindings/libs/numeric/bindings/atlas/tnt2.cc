
// BLAS level 2
// TNT arrays

#include <iostream>
#include <boost/numeric/bindings/atlas/cblas1.hpp>
#include <boost/numeric/bindings/atlas/cblas2.hpp>
#include <boost/numeric/bindings/traits/tnt.hpp>
#include "utils.h"
#include "tnt_utils.h"

namespace atlas = boost::numeric::bindings::atlas;

using std::cout;
using std::endl; 

#ifndef F_FORTRAN
typedef TNT::Array1D<double> vct_t;
typedef TNT::Array2D<double> matr_t;
#else
typedef TNT::Fortran_Array1D<double> vct_t;
typedef TNT::Fortran_Array2D<double> matr_t;
#endif 

int main() {

  cout << endl; 

  vct_t vx (2);
  atlas::set (1., vx);
  print_v (vx, "vx"); 
  vct_t vy (3); 
  atlas::set (0., vy); 
  print_v (vy, "vy"); 
  cout << endl; 

  matr_t m (3, 2);
  init_m (m, kpp (1)); 
  print_m (m, "m"); 
  cout << endl; 

  atlas::gemv (CblasNoTrans, 1.0, m, vx, 0.0, vy);
  print_v (vy, "m vx"); 

  atlas::gemv (m, vx, vy);
  print_v (vy, "m vx"); 
  cout << endl; 

  atlas::set (0, vx); 
  atlas::set (1, vy); 
  atlas::gemv (CblasTrans, 1.0, m, vy, 0.0, vx);
  print_v (vx, "m^T vy"); 
  cout << endl; 

  atlas::set (1, vy); 
  atlas::gemv (CblasNoTrans, 1.0, m, vx, 1.0, vy);
  print_v (vy, "vy + m vx"); 
  cout << endl; 

  atlas::set (1, vy); 
  atlas::gemv (CblasNoTrans, 2.0, m, vx, 0.5, vy);
  print_v (vy, "0.5 vy + 2.0 m vx"); 
  cout << endl; 

}
