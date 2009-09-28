//
//  Copyright Toon Knapen, Karl Meerbergen
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "random.hpp"

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/std_vector.hpp>
#include <boost/numeric/bindings/blas/blas1.hpp>

#include <vector>
#include <complex>
#include <iostream>
#include <limits>
#include <cmath>



// Randomize a vector (using functions from random.hpp)
template <typename V>
void randomize(V& v) {
   for (typename V::size_type i=0; i<v.size(); ++i)
      v[i] = random_value< typename V::value_type >() ;
} // randomize()


float abs_sum_value( float const& f ) {
  using namespace std ;
  return abs(f) ;
}

double abs_sum_value( double const& f ) {
  using namespace std ;
  return abs(f) ;
}

float abs_sum_value( std::complex< float > const& f ) {
  using namespace std ;
  return abs(f.real()) + abs(f.imag()) ;
}

double abs_sum_value( std::complex< double > const& f ) {
  using namespace std ;
  return abs(f.real()) + abs(f.imag()) ;
}

template <typename V>
typename boost::numeric::bindings::traits::type_traits<typename V::value_type>::real_type abs_sum( V const& v) {
  typedef typename boost::numeric::bindings::traits::type_traits<typename V::value_type>::real_type real_type ;

  real_type sum( 0.0 ) ;
  for ( typename V::size_type i=0; i<v.size(); ++i ) {
    sum += abs_sum_value( v[i] ) ;
  }
  return sum ;
}


// Blas operations using one vector.
template <typename T>
struct OneVector {
  boost::numeric::ublas::vector<T> v_ref_ ;

  // Initialize : set reference vector (ublas)
  OneVector()
  : v_ref_( 10 )
  {
     randomize(v_ref_);
  }

  template <typename V>
  int operator()(V& v) const {
     using namespace boost::numeric::bindings::blas ;

     typedef typename V::value_type                                                        value_type ;
     typedef typename boost::numeric::bindings::traits::type_traits<value_type>::real_type real_type ;

     // Copy vector from reference
     for (typename V::size_type i=0; i<v_ref_.size(); ++i)
        v[i] = v_ref_(i);

     // Test blas routines and compare with reference
     real_type nrm = nrm2( v );
     if ( std::abs(nrm - norm_2(v_ref_)) > std::numeric_limits< real_type >::epsilon() * norm_2(v_ref_)) {
       std::cout << "nrm2 : " << std::abs(nrm - norm_2(v_ref_)) << " > " << std::numeric_limits< real_type >::epsilon() * norm_2(v_ref_) << std::endl ;
       return 255 ;
     }

     nrm = asum( v );
     if ( std::abs(nrm - abs_sum(v_ref_)) > std::numeric_limits< real_type >::epsilon() * abs_sum(v_ref_)) {
       std::cout << "asum : " << std::abs(nrm - abs_sum(v_ref_)) << " > " << std::numeric_limits< real_type >::epsilon() * abs_sum(v_ref_) << std::endl ;
       return 255 ;
     }

     scal( value_type(2.0), v );
     for (typename V::size_type i=0; i<v_ref_.size(); ++i)
        if (std::abs( v[i] - real_type(2.0)*v_ref_(i) ) > real_type(2.0)*std::abs(v_ref_(i))) return 255 ;

     return 0;
  }

  // Return the size of a vector.
  size_t size() const {return v_ref_.size();}
};


// Operations with two vectors.
template <typename T, typename V>
struct BaseTwoVectorOperations {
  typedef T                                                                             value_type ;
  typedef typename boost::numeric::bindings::traits::type_traits<value_type>::real_type real_type ;
  typedef boost::numeric::ublas::vector<T>                                              ref_vector_type ;

  // Initialize: select the first vector and set the reference vectors (ublas)
  BaseTwoVectorOperations(V& v, const ref_vector_type& v1_ref, const ref_vector_type& v2_ref)
  : v_( v )
  , v1_ref_( v1_ref )
  , v2_ref_( v2_ref )
  {}

  // Copy the 2nd reference vector into w.
  template <typename W>
  void copy_vector(W& w) const {
     for (size_t i=0; i<size(); ++i) {
        w[i] = v2_ref_(i);
     }
  } // copy_vector()

  // Get the size of a vector.
  size_t size() const {return v_.size();}

  // Data members.
  V&                     v_ ;
  const ref_vector_type& v1_ref_, v2_ref_ ;
};


template <typename T, typename V>
struct TwoVectorOperations { } ;


template <typename V>
struct TwoVectorOperations< float, V>
: BaseTwoVectorOperations<float,V> {
  typedef typename V::value_type                                                        value_type ;
  typedef typename boost::numeric::bindings::traits::type_traits<value_type>::real_type real_type ;
  typedef typename BaseTwoVectorOperations<float,V>::ref_vector_type                    ref_vector_type ;

  TwoVectorOperations(V& v, const ref_vector_type& v1_ref, const ref_vector_type& v2_ref)
  : BaseTwoVectorOperations<float,V>( v, v1_ref, v2_ref )
  {}

  // Perform the tests of blas functions and compare with reference
  template <typename W>
  int operator()(W& w) const {
     using namespace boost::numeric::bindings::blas ;

     copy_vector(w);

     // Test blas routines
     value_type prod = dot( this->v_, w );
     if ( std::abs(prod - inner_prod( this->v1_ref_, this->v2_ref_ ))
          > std::numeric_limits< real_type >::epsilon() * std::abs(prod)) return 255 ;

     axpy( value_type(2.0), this->v_, w );
     for (size_t i=0; i<this->size(); ++i)
        if ( std::abs(w[i] - (this->v2_ref_(i) + value_type(2.0)*this->v1_ref_(i)))
          > std::numeric_limits< real_type >::epsilon() * std::abs(w[i])) return 255 ;

     scal( value_type(0.0), w ) ;
     copy( this->v_, w ) ;
     for (size_t i=0; i<this->size(); ++i) {
        if ( std::abs( w[i] - this->v_[i] ) != 0.0 ) return 255 ;
     }

     return 0;
  }
};


template <typename V>
struct TwoVectorOperations< double, V>
: BaseTwoVectorOperations<double,V> {
  typedef typename V::value_type                                                        value_type ;
  typedef typename boost::numeric::bindings::traits::type_traits<value_type>::real_type real_type ;
  typedef typename BaseTwoVectorOperations<double,V>::ref_vector_type                   ref_vector_type ;

  TwoVectorOperations(V& v, const ref_vector_type& v1_ref, const ref_vector_type& v2_ref)
  : BaseTwoVectorOperations<double,V>( v, v1_ref, v2_ref )
  {}

  // Perform the tests of blas functions and compare with reference
  template <typename W>
  int operator()(W& w) const {
     using namespace boost::numeric::bindings::blas ;

     copy_vector( w );

     // Test blas routines
     value_type prod = dot( this->v_, w );
     if ( std::abs(prod - inner_prod( this->v1_ref_, this->v2_ref_ ))
          > std::numeric_limits< real_type >::epsilon() * std::abs(prod)) return 255 ;

     axpy( value_type(2.0), this->v_, w );
     for (size_t i=0; i<this->size(); ++i)
        if ( std::abs(w[i] - (this->v2_ref_(i) + value_type(2.0)*this->v1_ref_(i)))
          > std::numeric_limits< real_type >::epsilon() * std::abs(w[i])) return 255 ;

     copy_vector( w ) ;
     scal( value_type(-1.0), w ) ;
     ::boost::numeric::bindings::blas::copy( this->v_, w ) ;
     for (size_t i=0; i<this->size(); ++i) {
        if ( w[i] != this->v_[i] ) return 255 ;
     }

     return 0;
  }
};


template <typename V>
struct TwoVectorOperations< std::complex<float>, V>
: BaseTwoVectorOperations< std::complex<float>, V>
{
  typedef typename V::value_type                                                        value_type ;
  typedef typename boost::numeric::bindings::traits::type_traits<value_type>::real_type real_type ;
  typedef typename BaseTwoVectorOperations<std::complex<float>,V>::ref_vector_type      ref_vector_type ;

  TwoVectorOperations(V& v, const ref_vector_type& v1_ref, const ref_vector_type& v2_ref)
  : BaseTwoVectorOperations< std::complex<float>, V>( v, v1_ref, v2_ref )
  {}

  // Perform the tests of blas functions and compare with reference
  template <typename W>
  int operator()(W& w) const {
     using namespace boost::numeric::bindings::blas ;

     copy_vector( w );

     // Test blas routines
     value_type prod = dotc( this->v_, w );
     if ( std::abs(prod - inner_prod( conj(this->v1_ref_), this->v2_ref_ ))
          > std::numeric_limits< real_type >::epsilon() * std::abs(prod)) return 255 ;

     prod = dotu( this->v_, w );
     if ( std::abs(prod - inner_prod( this->v1_ref_, this->v2_ref_ ))
          > std::numeric_limits< real_type >::epsilon() * std::abs(prod)) return 255 ;

     axpy( value_type(2.0), this->v_, w );
     for (size_t i=0; i<this->size(); ++i)
        if ( std::abs(w[i] - (this->v2_ref_(i) + value_type(2.0)*this->v1_ref_(i)))
          > std::numeric_limits< real_type >::epsilon() * std::abs(w[i])) return 255 ;

     scal( value_type(0.0), w ) ;
     copy( this->v_, w ) ;
     for (size_t i=0; i<this->size(); ++i) {
        if ( std::abs( w[i] - this->v_[i] ) != 0.0 ) return 255 ;
     }

     return 0;
  }
};


template <typename V>
struct TwoVectorOperations< std::complex<double>, V>
: BaseTwoVectorOperations< std::complex<double>, V>
{
  typedef typename V::value_type                                                        value_type ;
  typedef typename boost::numeric::bindings::traits::type_traits<value_type>::real_type real_type ;
  typedef typename BaseTwoVectorOperations<std::complex<double>,V>::ref_vector_type     ref_vector_type ;

  TwoVectorOperations(V& v, const ref_vector_type& v1_ref, const ref_vector_type& v2_ref)
  : BaseTwoVectorOperations< std::complex<double>, V>( v, v1_ref, v2_ref )
  {}

  // Perform the tests of blas functions and compare with reference
  template <typename W>
  int operator()(W& w) const {
     using namespace boost::numeric::bindings::blas ;

     copy_vector( w );

     // Test blas routines
     value_type prod = dotc( this->v_, w );
     if ( std::abs(prod - inner_prod( conj(this->v1_ref_), this->v2_ref_ ))
          > std::numeric_limits< real_type >::epsilon() * std::abs(prod)) return 255 ;

     prod = dotu( this->v_, w );
     if ( std::abs(prod - inner_prod( this->v1_ref_, this->v2_ref_ ))
          > std::numeric_limits< real_type >::epsilon() * std::abs(prod)) return 255 ;

     axpy( value_type(2.0), this->v_, w );
     for (size_t i=0; i<this->size(); ++i)
        if ( std::abs(w[i] - (this->v2_ref_(i) + value_type(2.0)*this->v1_ref_(i)))
          > std::numeric_limits< real_type >::epsilon() * std::abs(w[i])) return 255 ;

     scal( value_type(0.0), w ) ;
     copy( this->v_, w ) ;
     for (size_t i=0; i<this->size(); ++i) {
        if ( std::abs( w[i] - this->v_[i] ) != 0.0 ) return 255 ;
     }

     return 0;
  }
};


// Run the tests for different types of vectors.
template <typename T, typename F>
int different_vectors(const F& f) {
   // Do test for different types of vectors
   {
      std::cout << "  ublas::vector\n" ;
      boost::numeric::ublas::vector< T > v(f.size());
      if (f( v )) return 255 ;
   }
   { 
      std::cout << "  std::vector\n" ;
      std::vector<T> v_ref(f.size());
      if (f( v_ref )) return 255 ;
   }
   {
      std::cout << "  ublas::vector_range\n" ;
      typedef boost::numeric::ublas::vector< T > vector_type ;
      vector_type v(f.size()*2);
      boost::numeric::ublas::vector_range< vector_type > vr(v, boost::numeric::ublas::range(1,1+f.size()));
      if (f( vr )) return 255 ;
   }
   {
      typedef boost::numeric::ublas::matrix< T, boost::numeric::ublas::column_major >  matrix_type ;
      matrix_type  m(f.size(),f.size()) ;

      std::cout << "  ublas::matrix_column\n" ;
      boost::numeric::ublas::matrix_column< matrix_type > m_c( m, 2 );
      if (f( m_c )) return 255 ;

      std::cout << "  ublas::matrix_row\n" ;
      boost::numeric::ublas::matrix_row< matrix_type > m_r( m, 1 );
      if (f( m_r )) return 255 ;
   }
   return 0;
} // different_vectors()


// This is the functor that selects the first vector of the tests that use two vectors.
template <typename T>
struct TwoVector {
   TwoVector()
   : v1_ref_( 10 )
   , v2_ref_( 10 )
   {}

   template <typename V>
   int operator() (V& v) const {
      for (size_t i=0; i<size(); ++i) v[i] = v1_ref_(i) ;
      return different_vectors<T,TwoVectorOperations<T,V> >( TwoVectorOperations<T,V>(v, v1_ref_, v2_ref_) ) ;
   }

   size_t size() const {
      return v1_ref_.size() ;
   }

   boost::numeric::ublas::vector<T> v1_ref_ ;
   boost::numeric::ublas::vector<T> v2_ref_ ;
}; // TwoVector


// Run the test for a specific value_type T.
template <typename T>
int do_value_type() {
   // Tests for functions with one vector argument.
   std::cout << " one argument\n";
   if (different_vectors<T,OneVector<T> >(OneVector<T> ())) return 255 ;

   // Tests for functions with two vector arguments.
   std::cout << " two arguments\n";
   if (different_vectors<T,TwoVector<T> >(TwoVector<T>())) return 255;
   return 0;
} // do_value_type()


int main() {
  // Run regression for Real/Complex
  std::cout << "float\n"; if (do_value_type<float>() ) return 255 ;
  std::cout << "double\n"; if (do_value_type<double>() ) return 255 ;
  std::cout << "complex<float>\n"; if (do_value_type<std::complex<float> >() ) return 255 ;
  std::cout << "complex<double>\n"; if (do_value_type<std::complex<double> >() ) return 255 ;

  std::cout << "Regression test successful\n" ;

  return 0 ;
}


