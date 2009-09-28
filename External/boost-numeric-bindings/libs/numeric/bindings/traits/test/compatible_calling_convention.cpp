/* This test serves to test if the calling conventions
 * between C and Fortran are compatible.
 * This serves to correctly make bindings for fortran libraries
 * like BLAS and LAPACK.
 */

#include <iostream>
#include <complex>

/*
 * signatures for functions defined in 'fortran_functions.f'
 */
extern "C" {
  double dfunction() ;
  void dzfunction(double*,double*,double*) ;
}

int main()
{
  {
    double d = dfunction() ;
    if ( ( d - 9.87 ) > 1e-5 ) {
      std::cerr << "return value was " << d << " instead of 9.87" << std::endl ;
      return 1 ;
    }
  }

  {
    std::complex< double > dc ; 
    double* p = (double*)&dc ; 
    double d = 9.87, i = 3.21 ;
    dzfunction(p,&d,&i) ;
    if ( ( dc.real() - 9.87 ) > 1e-5 && ( dc.imag() - 3.21 ) > 1e-5 ) {
      std::cerr << "return value was " << dc << " instead of (9.87,3.21)" << std::endl ;
      return 2 ;
    }
  }

  return 0 ;
} 
