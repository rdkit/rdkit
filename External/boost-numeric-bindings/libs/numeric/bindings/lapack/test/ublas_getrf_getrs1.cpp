#include "../../blas/test/blas.hpp"
#include <boost/numeric/bindings/lapack/gesv.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>

template < typename matrix_type >
void test_getrf_getrs(matrix_type& lu, matrix_type& x)
{
  typedef typename matrix_type::value_type value_type ;

  numerics::matrix< value_type > a( lu ) ; // tmp to verify result
  numerics::matrix< value_type > b( x ) ;  // tmp to verify result
  std::vector< int > ipiv( x.size1() ) ;

  boost::numeric::bindings::lapack::getrf( lu, ipiv ) ;
  matrix_type ia( lu );
  boost::numeric::bindings::lapack::getrs( 'N', lu, ipiv, x ) ;
  boost::numeric::bindings::lapack::getri( ia, ipiv ) ;

  std::cout << prod(a,x) - b << std::endl ; 
  std::cout << prod(a,ia) << std::endl ; 
}

template < typename value_type, typename orientation, int size >
void test_getrf_getrs_matrix()
{
  numerics::matrix< value_type, orientation > a(size,size) ;
  random_initialise_matrix( a ) ;

  numerics::matrix< value_type, orientation > b(size,1) ;
  random_initialise_matrix( b ) ;

  test_getrf_getrs( a, b ) ;
}

int main()
{
  const int size = 5 ;

  test_getrf_getrs_matrix< double, numerics::column_major, size >() ;
  test_getrf_getrs_matrix< std::complex< double >, numerics::column_major, size >() ;

/*
  test_getrf_getrs_matrix< double, numerics::row_major, size >() ;
  test_getrf_getrs_matrix< std::complex< double >, numerics::row_major, size >() ;
*/

  return 0 ;
}
