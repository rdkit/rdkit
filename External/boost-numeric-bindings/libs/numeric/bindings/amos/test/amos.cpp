#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/amos/amos.hpp>

template <class T>
int do_value_type()
{
  typedef T value_type;
  //boost::numeric::ublas::vector< value_type > z(10) ; 
  value_type z ;
  boost::numeric::ublas::vector< value_type > cy(10) ; 
  int nz;
  boost::numeric::bindings::amos::besi( z, 1.0, 1, cy, nz ) ;

  return 0 ;
}

int main()
{
  if (do_value_type< std::complex<float> >()) return 255;
//  if (do_value_type< double >()) return 255;
  return 0;
}

