//
//  Copyright Karl Meerbergen, 2008
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "../../blas/test/random.hpp"

#include <boost/numeric/bindings/lapack/ptsv.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <iostream>


namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;


template <typename B, typename X>
bool check_residual( B const& b, X const& x ) {
  typedef typename B::value_type value_type ;

  ublas::matrix<value_type, ublas::column_major> res( b ) ;
  row(res,0).minus_assign( value_type(2.0) * row(x,0) - row(x,1) ) ;
  for (int i=1; i<res.size1()-1; ++i) {
    row(res,i).minus_assign( value_type(2.0) * row(x,i) - row(x,i+1) - row(x,i-1) ) ;
  }
  row(res,res.size1()-1).minus_assign( value_type(2.0) * row(x,res.size1()-1) - row(x,res.size1()-2) ) ;

  return norm_frobenius(res)<norm_frobenius(b)*1.e-5 ;
} // check_residual()

template <typename T>
int do_value_type() {
   const int n = 8 ;
   typedef typename boost::numeric::bindings::traits::type_traits<T>::real_type real_type ;

   typedef ublas::matrix<T, ublas::column_major>     matrix_type ;

   // Set matrix
   int const nrhs = 1 ;
   matrix_type b( n, nrhs );
   ublas::vector< real_type > d( n );
   ublas::vector<T>           e( n-1 );

   std::fill( d.begin(), d.end(), 2.0 ) ;
   std::fill( e.begin(), e.end(), -1.0 ) ;

   for (int i=0; i<b.size1(); ++i) b(i,0) = random_value<T>() ;

   // Factorize and solve
   matrix_type x( b );
   if( lapack::ptsv( d, e, x ) ) return -1 ;
   if (!check_residual(b,x)) return 1 ;

   // Restart computations
   std::fill( d.begin(), d.end(), 2.0 ) ;
   std::fill( e.begin(), e.end(), -1.0 ) ;

   // Compute factorization.
   if( lapack::pttrf( d, e ) ) return -1 ;

   // Compute solve
   x.assign( b ) ;
   if( lapack::pttrs( 'U', d, e, x ) ) return -2 ;

   if (!check_residual(b,x)) return 1 ;

   x.assign( b ) ;
   if( lapack::pttrs( 'L', d, e, x ) ) return -3 ;

   if (!check_residual(b,x)) return 2 ;

   return 0 ;
} // do_value_type()


int main() {
   // Run tests for different value_types
   std::cout << "double\n" ;
   if (do_value_type< double >()) return 255;

   std::cout << "float\n" ;
   if (do_value_type< float >()) return 255;

   std::cout << "complex<double>\n" ;
   if (do_value_type< std::complex<double> >()) return 255;

   std::cout << "complex<float>\n" ;
   if (do_value_type< std::complex<float> >()) return 255;

   std::cout << "Regression test succeeded\n" ;
   return 0;
}

