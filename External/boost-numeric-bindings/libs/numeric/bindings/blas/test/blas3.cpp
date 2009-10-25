//
//  Copyright Toon Knapen, Karl Meerbergen
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "random.hpp"

#include <boost/numeric/bindings/traits/transpose.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_symmetric.hpp>
#include <boost/numeric/bindings/traits/ublas_hermitian.hpp>
#include <boost/numeric/bindings/blas/blas3.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/hermitian.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <iostream>
#include <complex>


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


template <typename T>
struct RealDiagonal {
   template <typename M>
   void operator() (M& m) const {}
};


template <typename T>
struct RealDiagonal< std::complex<T> > {
   template <typename M>
   void operator() (M& m) const {
      for (typename M::size_type i=0; i<m.size1(); ++i) {
         m(i,i) = m(i,i).real();
      }
   }
};


template <typename M>
boost::numeric::ublas::triangular_adaptor<const M, boost::numeric::ublas::unit_upper> unit_upper_part(const M& m) {
   return boost::numeric::ublas::triangular_adaptor<const M, boost::numeric::ublas::unit_upper> (m);
}

template <typename M>
boost::numeric::ublas::triangular_adaptor<const M, boost::numeric::ublas::upper> upper_part(const M& m) {
   return boost::numeric::ublas::triangular_adaptor<const M, boost::numeric::ublas::upper> (m);
}

template <typename M>
boost::numeric::ublas::triangular_adaptor<const M, boost::numeric::ublas::lower> lower_part(const M& m) {
   return boost::numeric::ublas::triangular_adaptor<const M, boost::numeric::ublas::lower> (m);
}


// Run over all matrix types and aply a functor to it.
template <typename F>
int do_matrix_types(const F& f) {
   using namespace boost::numeric::ublas ;

   typedef typename F::value_type value_type ;

   {
      matrix<value_type,column_major> m( f.size1(), f.size2() );
      if (f(m)) return 255;
   }
   {
      typedef matrix<value_type,column_major> matrix_type ;
      matrix_type m( 2*f.size1(), 2*f.size2() );
      matrix_range<matrix_type> m_r( m, range(1,1+f.size1()), range(2,2+f.size2()) );
      if (f(m_r)) return 255 ;
   }

   return 0;
} // do_matrix_types()


// Functor to select syrk/herk calls given the two matrix_types
template <typename M1>
struct Syrk2 {
   typedef typename M1::value_type                                                         value_type ;
   typedef boost::numeric::ublas::matrix<value_type, boost::numeric::ublas::column_major>  ref_matrix_type ;
   typedef typename boost::numeric::bindings::traits::type_traits< value_type >::real_type real_type ;
   typedef typename ref_matrix_type::size_type                                             size_type ;

   Syrk2(const M1& a,
        const ref_matrix_type& a_ref)
   : a_( a ), a_ref_( a_ref )
   , c_ref_( a_ref_.size1(), a_ref_.size1() )
   {
      randomize( c_ref_ );
      RealDiagonal<value_type>() (c_ref_); // Must have real diagonal for herk.
   }

   template <typename M>
   int operator() (M& c) const {
      using namespace boost::numeric::ublas ;

      real_type alpha = 2.0 ;
      real_type beta = -3.0 ;

      c.assign( c_ref_ );
      boost::numeric::bindings::blas::syrk( 'U', 'N', value_type(alpha), a_, value_type(beta), c ) ;
      if ( norm_frobenius( upper_part( c - (beta*c_ref_ + alpha * prod( a_ref_, trans( a_ref_ ) ) ) ) )
          > std::numeric_limits< real_type >::epsilon() * norm_frobenius( upper_part(c) ) ) return 255;

      c.assign( c_ref_ );
      symmetric_adaptor<M, upper> c_s( c );
      boost::numeric::bindings::blas::syrk( 'U', 'N', value_type(alpha), a_, value_type(beta), c_s ) ;
      if ( norm_frobenius( upper_part( c_s - (beta*c_ref_ + alpha * prod( a_ref_, trans( a_ref_ ) ) ) ) )
          > std::numeric_limits< real_type >::epsilon() * norm_frobenius( upper_part(c_s) ) ) return 255;

      c.assign( c_ref_ );
      boost::numeric::bindings::blas::herk( 'U', 'N', alpha, a_, beta, c ) ;
      if ( norm_frobenius( upper_part( c - (beta*c_ref_ + alpha * prod( a_ref_, herm( a_ref_ ) ) ) ) )
          > std::numeric_limits< real_type >::epsilon() * norm_frobenius( upper_part(c) ) ) return 255;

      c.assign( c_ref_ );
      hermitian_adaptor<M, upper> c_h( c );
      boost::numeric::bindings::blas::herk( 'U', 'N', alpha, a_, beta, c_h ) ;
      if ( norm_frobenius( upper_part( c_h - (beta*c_ref_ + alpha * prod( a_ref_, herm( a_ref_ ) ) ) ) )
          > std::numeric_limits< real_type >::epsilon() * norm_frobenius( upper_part(c_h) ) ) return 255;

      return 0;
   }

   size_type size1() const {return c_ref_.size1() ;}
   size_type size2() const {return c_ref_.size2() ;}

   const M1& a_ ;
   const ref_matrix_type& a_ref_;
   ref_matrix_type  c_ref_;
}; // Syrk2


template <typename T>
struct Syrk1 {
   typedef T                                                                               value_type ;
   typedef boost::numeric::ublas::matrix<value_type, boost::numeric::ublas::column_major>  ref_matrix_type ;
   typedef typename boost::numeric::bindings::traits::type_traits< value_type >::real_type real_type ;
   typedef typename ref_matrix_type::size_type                                             size_type ;

   Syrk1()
   : a_ref_( 4, 7 )
   {
      randomize( a_ref_ );
   }

   template <typename M>
   int operator() (M& a) const {
      typedef Syrk2<M> functor_type ;
      a.assign( a_ref_ );
      return do_matrix_types( functor_type(a, a_ref_) );
   }

   size_type size1() const {return a_ref_.size1();}
   size_type size2() const {return a_ref_.size2();}

   ref_matrix_type a_ref_;
}; // Syrk1


// Functor to select gemm calls given the three matrix_types
template <typename M1, typename M2>
struct Gemm3 {
   typedef typename M1::value_type                                                         value_type ;
   typedef boost::numeric::ublas::matrix<value_type, boost::numeric::ublas::column_major>  ref_matrix_type ;
   typedef typename boost::numeric::bindings::traits::type_traits< value_type >::real_type real_type ;
   typedef typename ref_matrix_type::size_type                                             size_type ;

   Gemm3(const M1& a, const M2& b,
        const ref_matrix_type& a_ref, const ref_matrix_type& b_ref)
   : a_( a ), b_( b )
   , a_ref_( a_ref ), b_ref_( b_ref ), c_ref_( a_ref_.size1(), b_ref_.size2() )
   {
      randomize( c_ref_ );
   }

   template <typename M>
   int operator() (M& c) const {
      c.assign( c_ref_ );

      value_type alpha = 2.0 ;
      value_type beta = -3.0 ;

      boost::numeric::bindings::blas::gemm( 'N', 'N', alpha, a_, b_, beta, c ) ;
      if ( norm_frobenius( c - (beta*c_ref_ + alpha * prod( a_ref_, b_ref_ ) ) )
          > std::numeric_limits< real_type >::epsilon() * norm_frobenius( c ) ) return 255;

      return 0;
   }

   size_type size1() const {return c_ref_.size1() ;}
   size_type size2() const {return c_ref_.size2() ;}

   const M1& a_ ;
   const M2& b_ ;
   const ref_matrix_type& a_ref_, b_ref_;
   ref_matrix_type  c_ref_;
}; // Gemm3


template <typename M1>
struct Gemm2 {
   typedef typename M1::value_type                                                         value_type ;
   typedef boost::numeric::ublas::matrix<value_type, boost::numeric::ublas::column_major>  ref_matrix_type ;
   typedef typename boost::numeric::bindings::traits::type_traits< value_type >::real_type real_type ;
   typedef typename ref_matrix_type::size_type                                             size_type ;

   Gemm2(const M1& a, const ref_matrix_type& a_ref)
   : a_( a )
   , a_ref_( a_ref )
   , b_ref_( a_ref.size2(), 7 )
   {
      randomize( b_ref_ );
   }

   template <typename M>
   int operator() (M& b) const {
      typedef Gemm3<M1,M> functor_type ;
      b.assign( b_ref_ );
      return do_matrix_types( functor_type(a_, b, a_ref_, b_ref_) );
   }

   size_type size1() const {return b_ref_.size1();}
   size_type size2() const {return b_ref_.size2();}

   const M1&              a_ ;
   const ref_matrix_type& a_ref_;
   ref_matrix_type        b_ref_;
}; // Gemm2


template <typename T>
struct Gemm1 {
   typedef T                                                                               value_type ;
   typedef boost::numeric::ublas::matrix<value_type, boost::numeric::ublas::column_major>  ref_matrix_type ;
   typedef typename boost::numeric::bindings::traits::type_traits< value_type >::real_type real_type ;
   typedef typename ref_matrix_type::size_type                                             size_type ;

   Gemm1()
   : a_ref_( 4, 7 )
   {
      randomize( a_ref_ );
   }

   template <typename M>
   int operator() (M& a) const {
      typedef Gemm2<M> functor_type ;
      a.assign( a_ref_ );
      return do_matrix_types( functor_type(a, a_ref_) );
   }

   size_type size1() const {return a_ref_.size1();}
   size_type size2() const {return a_ref_.size2();}

   ref_matrix_type a_ref_;
}; // Gemm1


// Functor to select syrk/herk calls given the two matrix_types
template <typename M1>
struct Trsm2 {
   typedef typename M1::value_type                                                         value_type ;
   typedef boost::numeric::ublas::matrix<value_type, boost::numeric::ublas::column_major>  ref_matrix_type ;
   typedef typename boost::numeric::bindings::traits::type_traits< value_type >::real_type real_type ;
   typedef typename ref_matrix_type::size_type                                             size_type ;

   Trsm2(M1& c,
        const ref_matrix_type& b_ref, const ref_matrix_type& c_ref )
   : c_( c ), c_ref_( c_ref ), b_ref_( b_ref )
   , a_ref_( b_ref_.size1(), c_ref_.size2() )
   {
      randomize( a_ref_ ) ;
   }

   template <typename M>
   int operator() (M& a) const {
      using namespace boost::numeric::ublas ;

      real_type alpha = 2.0 ;

      a.assign( a_ref_ );
      boost::numeric::bindings::blas::trsm( 'L', 'U', 'N', 'N', value_type(alpha), b_ref_, a ) ;
      if ( norm_frobenius( alpha * a_ref_ - prod( upper_part(b_ref_), a ) )
          > std::numeric_limits< real_type >::epsilon() * norm_frobenius( upper_part(b_ref_) ) * norm_frobenius( a) ) return 255;

      a.assign( a_ref_ );
      boost::numeric::bindings::blas::trsm( 'R', 'U', 'N', 'N', value_type(alpha), c_, a ) ;
      if ( norm_frobenius( alpha * a_ref_ - prod( a, upper_part(c_) ) )
          > std::numeric_limits< real_type >::epsilon() * norm_frobenius( upper_part(c_) ) * norm_frobenius( a ) ) return 255;

      a.assign( a_ref_ );
      boost::numeric::bindings::blas::trsm( 'R', 'L', 'N', 'N', value_type(alpha), c_, a ) ;
      if ( norm_frobenius( alpha * a_ref_ - prod( a, lower_part(c_) ) )
          > std::numeric_limits< real_type >::epsilon() * norm_frobenius( lower_part(c_) ) * norm_frobenius( a ) ) return 255;

      a.assign( a_ref_ );
      boost::numeric::bindings::blas::trsm( 'L', 'U', 'T', 'N', value_type(alpha), b_ref_, a ) ;
      if ( norm_frobenius( alpha * a_ref_ - prod( trans(upper_part(b_ref_)), a ) )
          > std::numeric_limits< real_type >::epsilon() * norm_frobenius( upper_part(b_ref_) ) * norm_frobenius( a) ) return 255;

      a.assign( a_ref_ );
      boost::numeric::bindings::blas::trsm( 'L', 'U', 'N', 'U', value_type(alpha), b_ref_, a ) ;
      if ( norm_frobenius( alpha * a_ref_ - prod( unit_upper_part(b_ref_), a ) )
          > std::numeric_limits< real_type >::epsilon() * norm_frobenius( unit_upper_part(b_ref_) ) * norm_frobenius( a) ) return 255;

      return 0;
   }

   size_type size1() const {return a_ref_.size1() ;}
   size_type size2() const {return a_ref_.size2() ;}

   const M1& c_ ;
   const ref_matrix_type& c_ref_;
   const ref_matrix_type& b_ref_;
   ref_matrix_type        a_ref_;
} ; // Trsm2


template <typename T>
struct Trsm1 {
   typedef T                                                                               value_type ;
   typedef boost::numeric::ublas::matrix<value_type, boost::numeric::ublas::column_major>  ref_matrix_type ;
   typedef typename boost::numeric::bindings::traits::type_traits< value_type >::real_type real_type ;
   typedef typename ref_matrix_type::size_type                                             size_type ;

   Trsm1()
   : b_( 4, 4 )
   , c_( 7, 7 )
   {
      randomize( b_ );
      randomize( c_ );
      for ( std::size_t i=0; i<b_.size1(); ++i ) if ( b_(i,i)==value_type(0.0) ) b_(i,i) = value_type(1.0) ;
      for ( std::size_t i=0; i<c_.size1(); ++i ) if ( c_(i,i)==value_type(0.0) ) c_(i,i) = value_type(1.0) ;
   }

   template <typename M>
   int operator() (M& c) const {
      typedef Trsm2<M> functor_type ;
      return do_matrix_types( functor_type(c, b_, c_) );
   }

   size_type size1() const {return c_.size1();}
   size_type size2() const {return c_.size2();}

   ref_matrix_type b_;
   ref_matrix_type c_;
} ; // Trsm1



template <typename T>
int do_value_type() {
   // Gemm test
   if (do_matrix_types( Gemm1<T>() )) return 255 ;

   // Syrk and herk test
   if (do_matrix_types( Syrk1<T>() )) return 255 ;

   // Trsm test
   if (do_matrix_types( Trsm1<T>() )) return 255 ;

   return 0 ;
} // do_value_type()


int main() {
   // Test for various value_types
   if (do_value_type<float>()) return 255 ;
   if (do_value_type<double>()) return 255 ;
   if (do_value_type<std::complex<float> >()) return 255 ;
   if (do_value_type<std::complex<double> >()) return 255 ;

   std::cout << "Regression test succeeded\n" ;

   return 0 ;
}
