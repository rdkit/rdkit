//
//  Copyright Toon Knapen, Karl Meerbergen
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "ublas_heev.hpp"

#include <boost/numeric/bindings/lapack/hegv.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <iostream>


namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;


template <typename T, typename W, char UPLO>
int do_memory_uplo(int n, W& workspace ) {
   typedef typename boost::numeric::bindings::traits::type_traits<T>::real_type real_type ;

   typedef ublas::matrix<T, ublas::column_major>     matrix_type ;
   typedef ublas::vector<real_type>                  vector_type ;

   // Set matrix
   matrix_type a( n, n ); a.clear();
   vector_type e1( n );
   vector_type e2( n );

   fill( a );
   matrix_type a2( a );

   matrix_type b( n, n ); b.clear();
   for (int i = 0; i < n; ++i) b(i,i) = 1;

   // Compute eigen decomposition.
   lapack::hegv( 1, 'V', UPLO, a, b, e1, workspace ) ;

   if (check_residual( a2, e1, a )) return 255 ;

   lapack::hegv( 1, 'N', UPLO, a2, b, e2, workspace ) ;
   if (norm_2( e1 - e2 ) > n * norm_2( e1 ) * std::numeric_limits< real_type >::epsilon()) return 255 ;

   // Test for a matrix range
   fill( a ); a2.assign( a );

   typedef ublas::matrix_range< matrix_type > matrix_range ;

   ublas::range r(1,n-1) ;
   matrix_range a_r( a, r, r );
   ublas::vector_range< vector_type> e_r( e1, r );
   matrix_range b_r( b, r, r );

   lapack::hegv(1, 'V', UPLO,  a_r, b_r, e_r, workspace );

   matrix_range a2_r( a2, r, r );
   if (check_residual( a2_r, e_r, a_r )) return 255 ;

   return 0 ;
} // do_memory_uplo()


template <typename T, typename W>
int do_memory_type(int n, W workspace) {
   std::cout << "  upper\n" ;
   if (do_memory_uplo<T,W,'U'>(n, workspace)) return 255 ;
   std::cout << "  lower\n" ;
   if (do_memory_uplo<T,W,'L'>(n, workspace)) return 255 ;
   return 0 ;
}


template <typename T>
struct Workspace {
   typedef ublas::vector<T>                         array_type ;
   typedef lapack::detail::workspace1< array_type > type ;

   Workspace(size_t n)
   : work_( 3*n-1 )
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
   : work_( 2*n-1 )
   , rwork_( 3*n-2 )
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

   std::cout << " optimal workspace\n";
   if (do_memory_type<T,lapack::optimal_workspace>( n, lapack::optimal_workspace() ) ) return 255 ;

   std::cout << " minimal workspace\n";
   if (do_memory_type<T,lapack::minimal_workspace>( n, lapack::minimal_workspace() ) ) return 255 ;

   std::cout << " workspace array\n";
   Workspace<T> work( n );
   do_memory_type<T,typename Workspace<T>::type >( n, work() );
   return 0;
} // do_value_type()


int main() {
   // Run tests for different value_types
   std::cout << "float\n" ;
   if (do_value_type<float>()) return 255;

   std::cout << "double\n" ;
   if (do_value_type<double>()) return 255;

   std::cout << "complex<float>\n" ;
   if (do_value_type< std::complex<float> >()) return 255;

   std::cout << "complex<double>\n" ;
   if (do_value_type< std::complex<double> >()) return 255;

   std::cout << "Regression test succeeded\n" ;
   return 0;
}

