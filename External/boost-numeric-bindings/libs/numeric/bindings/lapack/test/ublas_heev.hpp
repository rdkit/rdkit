//
//  Copyright Toon Knapen, Karl Meerbergen
//  Copyright Thomas Klimpel 2008
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_TEST_UBLAS_HEEV_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_TEST_UBLAS_HEEV_HPP

#include "../../blas/test/random.hpp"

#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <limits>


inline float conj(float v) { return v; }
inline double conj(double v) { return v; }
inline float real(float v) { return v; }
inline double real(double v) { return v; }

// Fill a banded matrix
template <typename M>
void fill_banded(M& m) {
   typedef typename M::size_type  size_type ;
   typedef typename M::value_type value_type ;

   int size = m.size2() ;
   int band = m.upper() ;

   for (int i=0; i<size; ++i) {
      for (int j=std::max(0,i-band); j<i; ++j) m(j,i) = random_value<value_type>();
      m(i,i) = real( random_value<value_type>() );
   }
} // randomize()

// Fill a matrix
template <typename M>
void fill(M& m) {
   typedef typename M::size_type  size_type ;
   typedef typename M::value_type value_type ;

   int size = m.size2() ;

   for (int i=0; i<size; ++i) {
      for (int j=0; j<i; ++j) {
         m(j,i) = random_value<value_type>();
         m(i,j) = conj( m(j,i) ) ;
      }
      m(i,i) = real( random_value<value_type>() );
   }
} // randomize()

template <typename H, typename E, typename Z>
int check_residual(H const& h, E const& e, Z const& z) {
   typedef typename H::value_type value_type ;
   typedef typename E::value_type real_type ;
   real_type safety_factor (1.5);

   // Check eigen decomposition
   int n = h.size1();
   boost::numeric::ublas::matrix<value_type> error( n, n ); error.clear();

   // Copy band matrix in error
   error.assign( h );
   assert( norm_frobenius( error - herm( error ) ) == 0.0 ) ;

   for (int i=0; i<n; ++i) {
      error .minus_assign( outer_prod( column(z, i), e(i) * conj( column(z, i) ) ) ) ;
   }
   return (norm_frobenius( error )
           >= safety_factor*n* norm_2( e ) * std::numeric_limits< real_type >::epsilon() ) ;
} // check_residual()

#endif
