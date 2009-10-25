#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/traits/ublas_vector2.hpp>
#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/bindings/mumps/mumps_driver.hpp>
#include <iostream>
#include <fstream>
#include <complex>

template <typename T>
int test() {
  namespace ublas = ::boost::numeric::ublas ;
  namespace mumps = ::boost::numeric::bindings::mumps ;


  int const n = 10 ;

  typedef ublas::coordinate_matrix<T, ublas::column_major, 1, ublas::unbounded_array<int> > coo_type ;

  coo_type coo( n, n, n + 6 ) ;

  for (int i=0; i<n; ++i) coo(i,i) = i+1.0 ;
  coo(2,3) = T(1.0) ;
  coo(2,4) = T(1.0) ;
  coo(5,6) = T(-1.0) ;
  coo(2,6) = T(1.0) ;
  coo(9,0) = T(1.0) ;
  coo(2,7) = T(-1.0) ;

  coo.sort() ;
  std::cout << "matrix " << coo << std::endl ;

  ublas::vector<T> v( 10 ) ;
  ublas::vector<T> w( 10 ) ;

  std::fill( w.begin(), w.end(), 1.0 ) ;

  for (int i=1; i<n; ++i) {
    w(i) += w(i-1) ;
  }

  for (int i=0; i<n; ++i) {
    v[i] = T(coo(i,i)) * w[i] ;
  }
  v[2] += T(coo(2,3)) * w[3] ;
  v[2] += T(coo(2,4)) * w[4] ;
  v[5] += T(coo(5,6)) * w[6] ;
  v[2] += T(coo(2,6)) * w[6] ;
  v[9] += T(coo(9,0)) * w[0] ;
  v[2] += T(coo(2,7)) * w[7] ;
  std::cout << "rhs : " << v << std::endl ;

  mumps::mumps< coo_type > mumps_coo ;

  mumps_coo.icntl[2]=mumps_coo.icntl[3] = 0 ;

  // Analysis
  mumps_coo.job = 1 ;
  matrix_integer_data( mumps_coo, coo ) ;
  driver( mumps_coo ) ;

  // Factorization
  mumps_coo.job = 2 ;
  matrix_value_data( mumps_coo, coo ) ;
  driver( mumps_coo ) ;

  // Solve
  mumps_coo.job = 3 ;
  rhs_sol_value_data( mumps_coo, v ) ;
  driver( mumps_coo ) ;

  std::cout << "w : " << w << std::endl ;
  std::cout << "v : " << v << std::endl ;

  if ( norm_2( v - w ) > 1.e-10 * norm_2( v ) ) return 1 ;

  return 0 ;
}

int main() {
  if ( test<float>() ) return 1 ;
  if ( test<double>() ) return 2 ;
  if ( test< std::complex<float> >() ) return 3 ;
  if ( test< std::complex<double> >() ) return 4 ;
  return 0 ;
}
