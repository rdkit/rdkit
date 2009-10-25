#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>

int main()
{
  using namespace boost::numeric::ublas ;
  using namespace boost::numeric::bindings::traits ;

  vector< double > v( 10 ) ;  
  range r(2, 5) ;
  slice s(1,3,3) ;
  vector_range< vector< double > >  vr( v, r ) ;
  vector_slice< vector< double > >  vs( v, s ) ;
  
  if ( vector_storage( v ) != & v[0] )                   return 1 ;
  if ( vector_size( v ) != 10 )                          return 2 ;
  if ( vector_stride( v ) != 1 )                         return 3 ;

  if ( vector_storage( v ) + 2 != vector_storage( vr ) ) return 4 ;
  if ( vector_size( vr ) != 3 )                          return 5 ;
  if ( vector_stride( vr ) != 1 )                        return 6 ;

  if ( vector_storage( v ) + 1 != vector_storage( vs ) ) return 7 ;
  if ( vector_size( vs ) != 3 )                          return 8 ;
  if ( vector_stride( vs ) != 3 )                        return 9 ;

  return 0 ;
}
