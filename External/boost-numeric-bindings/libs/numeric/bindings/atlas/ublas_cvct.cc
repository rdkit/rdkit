
// BLAS level 1 (vector) -- complex numbers

#include <iostream>
#include <cmath>
#include <complex> 

//#define BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 
//#define BOOST_NO_FUNCTION_TEMPLATE_ORDERING

#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/atlas/cblas1.hpp>
#include "utils.h" 

namespace atlas = boost::numeric::bindings::atlas;
namespace ublas = boost::numeric::ublas;

using std::cout;
using std::endl; 

typedef double real_t;
typedef std::complex<real_t> cmplx_t;  
typedef ublas::vector<cmplx_t> vct_t; 

#ifdef BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 
using ublas::inner_prod; 
using ublas::norm_2; 
using ublas::conj; 
#endif 

int main() {

  int n = 6; 

  cout << endl; 
  vct_t v (n); 
  init_v (v, times_plus<cmplx_t> (cmplx_t (1, -1), cmplx_t (0, .1))); 
  print_v (v, "v"); 

  atlas::scal (2.0, v); 
  print_v (v, "2.0 v"); 

  atlas::scal (cmplx_t (-1, 0), v); 
  print_v (v, "(-1, 0) v"); 

  atlas::scal (cmplx_t (0, 1), v);
  print_v (v, "(0, 1) v"); 

  atlas::set (cmplx_t (1, -1), v); 
  print_v (v, "v"); 
  vct_t v1 (n); 
  atlas::set (cmplx_t (0, -1), v1); 
  print_v (v1, "v1"); 

  cout << endl; 
  cout << "v^T v1 = " << atlas::dot (v, v1) << " == "
    << atlas::dotu (v, v1) << " == "
    << inner_prod (v, v1) << endl; 
  cout << "v^T v = " << atlas::dot (v, v) << " == "
    << atlas::dotu (v, v) << " == "
    << inner_prod (v, v) << endl; 
  cout << "v1^T v1 = " << atlas::dot (v1, v1) << " == "
    << atlas::dotu (v1, v1) << " == "
    << inner_prod (v1, v1) << endl; 

  cout << endl; 
  cout << "v^H v1 = " << atlas::dotc (v, v1) << " == "
    << inner_prod (conj (v), v1) << " == "
    << inner_prod (v1, conj (v)) << " != "
    << inner_prod (v, conj (v1)) << endl; 
  cout << "v^H v = " << atlas::dotc (v, v) << " == "
    << inner_prod (conj (v), v) << " == "
    << inner_prod (v, conj (v)) << endl; 
  cout << "v1^H v1 = " << atlas::dotc (v1, v1) << " == "
    << inner_prod (conj (v1), v1) << " == "
    << inner_prod (v1, conj (v1)) << endl; 

  
  cout << endl;
  cout << "||v||_1 = " << atlas::asum (v) << endl; 
  cout << "||v||_2 = " << atlas::nrm2 (v) << " == "
    << norm_2 (v) << endl; 
  
  cout << endl;
}
