#include <boost/numeric/ublas/vector.hpp>

template < class value_type >
void test_scal(int size)
{
  boost::numeric::ublas::vector< value_type > vector( size ) ;

  // random generate content for the vector
}


int main()
{
  int max_size = 1000 ;
  typedef std::complex< double > dcomplex ;

  for(int size = 1 ; size < max_size ; ++size ) test_scal< double >( size ) ;
  for(int size = 1 ; size < max_size ; ++size ) test_scal< dcomplex >( size ) ;

  return 0 ;
}
