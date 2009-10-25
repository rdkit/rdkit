//
//  Copyright Toon Knapen, Karl Meerbergen
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "../../blas/test/random.hpp"

#include <boost/numeric/bindings/lapack/gees.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <iostream>
#include <limits>


namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;


// Randomize a matrix
template <typename M>
void randomize(M& m) {
   typedef typename M::size_type  size_type ;
   typedef typename M::value_type value_type ;

   size_type size1 = m.size1() ;
   size_type size2 = m.size2() ;

   for (size_type i=0; i<size2; ++i) {
      for (size_type j=0; j<size1; ++j) {
         m(j,i) = random_value< value_type >() ;
      }
   }
} // randomize()




template <typename T, typename W>
int do_memory_type(int n, W workspace) {
   typedef typename boost::numeric::bindings::traits::type_traits<T>::real_type real_type ;
   typedef std::complex< real_type >                                            complex_type ;

   typedef ublas::matrix<T, ublas::column_major> matrix_type ;
   typedef ublas::vector<complex_type>           vector_type ;

   // Set matrix
   matrix_type a( n, n );
   matrix_type z( n, n );
   vector_type e1( n );
   vector_type e2( n );

   randomize( a );
   matrix_type a2( a );

   // Compute Schur decomposition.
   lapack::gees( a, e1, z, workspace ) ;

   // Check Schur factorization
   if (norm_frobenius( prod( a2, z ) - prod( z, a ) )
           >= 10.0* norm_frobenius( a2 ) * std::numeric_limits< real_type >::epsilon() ) return 255 ;

   lapack::gees( a2, e2, workspace ) ;
   if (norm_2( e1 - e2 ) > norm_2( e1 ) * std::numeric_limits< real_type >::epsilon()) return 255 ;

   if (norm_frobenius( a2 - a )
           >= 10.0* norm_frobenius( a2 ) * std::numeric_limits< real_type >::epsilon() ) return 255 ;


   return 0 ;
} // do_value_type()


template <typename T>
struct Workspace {
   typedef ublas::vector<T>                         array_type ;
   typedef lapack::detail::workspace1< array_type > type ;

   Workspace(size_t n)
   : work_( 3*n )
   {}

   type operator() () {
      return type( work_ );
   }

   array_type work_ ;
};


template <typename T>
struct Workspace< std::complex<T> > {
   typedef ublas::vector<T>                                                 real_array_type ;
   typedef ublas::vector< std::complex<T> >                                 complex_array_type ;
   typedef lapack::detail::workspace2< complex_array_type,real_array_type > type ;

   Workspace(size_t n)
   : work_( 2*n )
   , rwork_( n )
   {}

   type operator() () {
      return type( work_, rwork_ );
   }

   complex_array_type work_ ;
   real_array_type    rwork_ ;
};


template <typename T>
int do_value_type() {
   const int n = 8 ;
   
   if (do_memory_type<T,lapack::optimal_workspace>( n, lapack::optimal_workspace() ) ) return 255 ;
   if (do_memory_type<T,lapack::minimal_workspace>( n, lapack::minimal_workspace() ) ) return 255 ;

   Workspace<T> work( n );
   do_memory_type<T,typename Workspace<T>::type >( n, work() );
   return 0;
} // do_value_type()


int main() {
   // Run tests for different value_types
   if (do_value_type<float>()) return 255;
   if (do_value_type<double>()) return 255;
   if (do_value_type< std::complex<float> >()) return 255;
   if (do_value_type< std::complex<double> >()) return 255;

   std::cout << "Regression test succeeded\n" ;
   return 0;
}

