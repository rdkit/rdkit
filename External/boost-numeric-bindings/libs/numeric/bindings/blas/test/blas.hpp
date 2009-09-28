#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <iomanip>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/timer.hpp>

namespace numerics = boost::numeric::ublas ;

template < typename T >
bool eq (const T& a, const T& b, double abstol)
{
   return std::abs(a-b)<abstol;
}

template < typename T >
struct is_equal
{
   is_equal (double tol) : abstol(tol) {}
   bool operator()(const T& a, const T&b) { return eq(a,b,abstol);}
   double abstol;
};

template < typename T >
void random_initialise(T& v)
{ v = 10.0 * rand()/(RAND_MAX+1.0) ; }

template < typename T >
void random_initialise_vector(numerics::vector< T >& v)
{
  size_t size = v.size();
  for(size_t i = 0 ; i < size ; ++i ) random_initialise( v[i] ) ;
}

template < typename T, typename Orientation >
void random_initialise_matrix(numerics::matrix< T, Orientation >& m)
{
  size_t size1 = m.size1();
  size_t size2 = m.size2();
  for(size_t i = 0 ; i < size1 ; ++i ) 
    for(size_t j = 0 ; j < size2 ; ++j )
      random_initialise( m(i,j) ) ;
}

template < typename T >
struct assign_multiplier
{
   T operator()() const
   { return 1.0000002 ; } // otherwise the accumulated result of the range will result in 0.0 or over/under-flow
};

template < typename T >
struct assign_multiplier< std::complex< T > >
{
  std::complex< T > operator()() const
  { return std::complex< T >( cos(0.5),sin(0.5) ) + ( assign_multiplier< T >().operator() - 1.0 ) ; }
};

template<class T>
double flops(int multiplies, int plus, int runs, double elapsed) 
{ return ( multiplies * boost::numeric::ublas::type_traits<T>::multiplies_complexity + plus * boost::numeric::ublas::type_traits<T>::plus_complexity ) * runs / (1024 * 1024 * elapsed) ; }

template<class T>
struct peak_c_plus {
  typedef T value_type;
  
  void operator () (int runs) const {
    try {
      static T s (0);
      boost::timer t;
      for (int i = 0; i < runs; ++ i) {
        s += T (0);
      }
      std::cerr << flops<value_type>(0, 1, runs, t.elapsed ()) << " Mflops\n";
    }
    catch (std::exception &e) {
      std::cerr << e.what () << std::endl;
    }
    catch (...) {
      std::cerr << "unknown exception" << std::endl;
    }
  }
};

template<class T>
struct peak_c_multiplies {
  typedef T value_type;

  void operator () (int runs) const {
    try {
      static T s (1);
      boost::timer t;
      for (int i = 0; i < runs; ++ i) {
        s *= T (1);
      }
      std::cerr << flops<value_type>(0, 1, runs, t.elapsed ()) << " Mflops\n";
    }
    catch (std::exception &e) {
      std::cerr << e.what () << std::endl;
    }
    catch (...) {
      std::cerr << "unknown exception" << std::endl;
    }
  }
};

template<class T>
struct peak
{
  void operator () (int runs) 
  {
    std::cerr << "plus       :";
    peak_c_plus<T> () (runs);
    std::cerr << "multiplies :";
    peak_c_multiplies<T> () (runs);
  }
};

template < typename value_type >
void check(value_type a, value_type b)
{
  if ( ! eq< value_type >( a, b, 1e-5 ) ) {
    std::cerr << "\n\nregression test failure : results are not identical" << std::endl;
    std::cerr << a << " != " << b << std::endl;
    exit( 1 );  
  }
}

template < typename Iterator0, typename Iterator1 >
void check(Iterator0 begin0, Iterator0 end0, Iterator1 begin1)
{
  Iterator0 it0 = begin0 ;
  Iterator1 end1 = begin1 + std::distance(begin0,end0);
  for(bool fail = false ; it0 != end0 || fail ; ++it0 ) ; // fail = ! ( *it0 < std::numeric_limits< typename Iterator0::value_type >::max() );

  Iterator1 it1 = begin1 ;
  for(bool fail = false ; it1 != end1 || fail ; ++it1 ) ; // fail = ! ( *it1 < std::numeric_limits< typename Iterator0::value_type >::max() ) ;

  if ( it0 != end0 || it1 != end1 ) {
    std::cerr << "\n\nregression test failure : results overflowed" << std::endl;
    std::copy( begin0, end0, std::ostream_iterator< typename Iterator0::value_type >( std::cerr, " " ) ); std::cerr << std::endl;
    std::copy( begin1, end1, std::ostream_iterator< typename Iterator1::value_type >( std::cerr, " " ) ); std::cerr << std::endl;
    exit(1);
  }

  if ( ! std::equal( begin0, end0, begin1, is_equal< typename Iterator0::value_type >( std::abs( *begin0 ) * 1e-5 ) ) ) {
    std::cerr << "\n\nregression test failure : results are not identical" << std::endl;
    std::cerr << std::setprecision( 20 ) ;
    std::copy( begin0, end0, std::ostream_iterator< typename Iterator0::value_type >( std::cerr, " " ) ); std::cerr << std::endl;
    std::copy( begin1, end1, std::ostream_iterator< typename Iterator1::value_type >( std::cerr, " " ) ); std::cerr << std::endl;
    exit( 1 );  
  }
}

template < typename value_type >
void check(numerics::matrix< value_type > &a, numerics::matrix< value_type > &b)
{
  bool ret = true ;
  size_t num_rows = a.size1() ;
  size_t num_cols = a.size2() ;

  for(size_t i = 0 ; i < num_rows && ret ; ++i ) 
    for(size_t j = 0 ; j < num_cols && ret ; ++j ) 
      ret = eq< value_type >( a(i,j), b(i,j), 1e-5 );

  if ( ! ret ) {
    std::cerr << "\n\nregression test failure : matrices not identical" << std::endl;
    exit( 1 );  
  }
}

template < typename T >
void report(std::ostream& os, int runs, int runs_i, int size_i, double time)
{
  double normed_time = runs * time / ( runs_i * size_i ) ;
  double mflops = flops<T>(size_i,0,runs_i,time) ;
  std::cerr << std::setw(12) << normed_time << std::setw(12) << mflops ;
  os        << std::setw(12) << normed_time << std::setw(12) << mflops ;
}

template < typename FunctorType >
void loop(std::ostream& os, int start, int step, int stop, int runs, FunctorType functor)
{
  for(int size_i = start ; size_i <= stop ; size_i = std::max( static_cast< int >( size_i * (1 + 1.0 / step ) ), size_i + 1 ) ) {
    int runs_i = 10 * runs / size_i ;
    
    std::cerr << size_i << "\t";
    os        << size_i << "\t";
    
    functor.operator()( os, stop, size_i, runs, runs_i ) ;
    
    std::cerr << std::endl;
    os        << std::endl;
  }
}

template < typename FunctorType, typename CallFunctor >
void loop(std::ostream& os, int start, int step, int stop, int runs, FunctorType functor, CallFunctor ublas_call)
{
  for(int size_i = start ; size_i <= stop ; size_i = std::max( static_cast< int >( size_i * (1 + 1.0 / step ) ), size_i + 1 ) ) {
    int runs_i = 10 * runs / size_i ;
    
    std::cerr << size_i << "\t";
    os        << size_i << "\t";
    
    functor.operator()( os, stop, size_i, runs, runs_i, ublas_call ) ;
    
    std::cerr << std::endl;
    os        << std::endl;
  }
}
