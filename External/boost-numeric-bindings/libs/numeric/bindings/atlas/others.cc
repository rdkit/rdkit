
// BLAS level 1
// std::vector<>, std::valarray<>, boost::array<>, C array

// element type: float or double
#ifdef F_FLOAT
typedef float real_t; 
#else
typedef double real_t; 
#endif

#include <iostream>
#include <iterator>
#include <algorithm>
#include <complex>

#include <boost/numeric/bindings/traits/std_vector.hpp>
#include <boost/numeric/bindings/traits/std_valarray.hpp>
#include <boost/numeric/bindings/traits/boost_array.hpp>
#include <boost/numeric/bindings/traits/c_array.hpp>
#include <boost/numeric/bindings/atlas/cblas1.hpp>

#include "utils.h"

namespace atlas = boost::numeric::bindings::atlas;

using std::cout;
using std::endl;
using std::size_t; 

int main() {

  int n = 10; 

  cout << endl; 
  cout << "std::vector" << endl; 
  std::vector<real_t> sv (n); 
  init_v (sv, kpp (1)); 
  print_v (sv, "sv"); 
  cout << "std::valarray" << endl; 
  std::valarray<real_t> va (n); 
  atlas::set (0.1, va); 
  print_v (va, "va"); 
  cout << endl; 

  cout << "dot(): sv^t va: ";
  real_t d = 0;
  for (int i = 0; i < n; ++i)
    d += sv[i] * va[i]; 

  cout << "is " << d << " == " << atlas::dot (sv, va) << " ?" << endl; 
  cout << endl; 

#ifdef F_FLOAT
  cout << "sdsdot(): 10 + sv^T va = " << atlas::sdsdot (10, sv, va) << endl; 
  cout << endl;
#endif  

  atlas::scal (real_t(2), sv);
  print_v (sv, "scal(): 2 sv"); 

  cout << endl; 

  std::random_shuffle (sv.begin(), sv.end());
  cout << "shuffled sv: "; 
  std::copy (sv.begin(), sv.end(), std::ostream_iterator<real_t> (cout, " ")); 
  cout << endl; 
  int i = atlas::iamax (sv); 
  cout << "iamax():\n  index of max el = " << i 
       << "; max el = " << sv[i] << endl; 
  cout << endl; 

  cout << "asum():\n  ||sv||_1 =  " << atlas::asum (sv) 
       << "; ||va||_1 = " << atlas::asum (va) << endl; 
  cout << "nrm2():\n  ||sv||_2 = " << atlas::nrm2 (sv) 
       << "; ||va||_2 = " << atlas::nrm2 (va) << endl; 
  cout << endl; 

  cout << "boost::array" << endl;
  boost::array<double, 10> ba;
  atlas::set (0.1, ba);
  print_v (ba, "ba");
  cout << "C array" << endl; 
  typedef double double_array[10]; 
  double_array ca; 
  atlas::set (1., ca); 
  print_v (ca, "ca");
  cout << endl; 
  
  atlas::axpy (0.1, ba, ca); 
  print_v (ca, "axpy(): 0.1 ba + ca"); 

  atlas::axpby (0.1, ba, 2., ca); 
  print_v (ca, "axpby(): 0.1 ba + 2.0 ca"); 

  cout << endl;
}
