
// BLAS level 1 (matrix rows and columns) 

//#define BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 
//#define BOOST_NO_FUNCTION_TEMPLATE_ORDERING

#include <iostream>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/atlas/cblas1.hpp>
#include "utils.h"

namespace atlas = boost::numeric::bindings::atlas;
namespace traits = boost::numeric::bindings::traits;
namespace ublas = boost::numeric::ublas;

using std::cout;
using std::endl; 
using std::size_t; 

typedef ublas::vector<double> vct_t;
typedef ublas::matrix<double> matr_t;
typedef ublas::matrix_row<matr_t> mr_t;
typedef ublas::matrix_column<matr_t> mc_t;
typedef ublas::matrix_range<matr_t> mrng_t;
typedef ublas::matrix_slice<matr_t> msl_t;
typedef ublas::matrix_slice<matr_t const> cmsl_t;

int main() {

  cout << endl; 

  int r = 7; 
  int c = 8; 
  matr_t m (r, c);
  init_m (m, times_plus<double> (10, 1, 1)); 
  print_m (m, "m"); 
  cout << endl; 

  vct_t v (10);
  atlas::set (0., v); 
  print_v (v, "v"); 
  cout << endl; 

  // m[2,.] <- 0.1 m[2,.]
  mr_t mr2 (m, 2); 
  atlas::scal (0.1, mr2);
  print_m (m, "0.1 m[2,.]"); 
  cout << endl; 
  
  // m[2,.] <-> m[4,.]
  mr_t mr4 (m, 4);
  atlas::swap (mr2, mr4);
  print_m (m, "m[2,.] <-> m[4,.]"); 
  cout << endl; 

  // m[4,.] <- m[2,.]
  atlas::copy (mr2, mr4);
  print_m (m, "m[4,.] <- m[2,.]"); 
  cout << endl; 

  // v[2..6] <- 10 m[5,.][1..5]
  mr_t mr5 (m, 5); 
  ublas::vector_range<vct_t> vr (v, ublas::range (2, 6)); 
  ublas::vector_range<mr_t> mr5r (mr5, ublas::range (1, 5)); 
  atlas::axpy (10.0, mr5r, vr);
  print_v (v, "v[2..6] <- 10 m[5,.][1..5]"); 
  cout << endl; 

  // ||m[.,3]||_1, ||m[.,3]||_2
  mc_t mc3 (m, 3); 
  cout << "||m[.,3]||_1 = " << atlas::asum (mc3) << endl; 
  cout << "||m[.,3]||_2 = " << atlas::nrm2 (mc3) << endl; 
  cout << endl; 
  
  // m[.,5] <- 0.01 m[.,3] + m[.,5]
  mc_t mc5 (m, 5); 
  atlas::axpy (0.01, mc3, mc5); 
  print_m (m, "m[.,5] <- 0.01 m[.,3] + m[.,5]"); 
  cout << endl; 

  // 0.1 m[.,5][1:2:3]
  ublas::vector_slice<mc_t> mc5s (mc5, ublas::slice (1, 2, 3)); 
  atlas::scal (0.1, mc5s); 
  print_m (m, "0.1 m[.,5][1:2:3]"); 
  cout << endl; 

  // 0.1 m[4,.][1:2:4][1..3]
  ublas::vector_slice<mr_t> mr4s (mr4, ublas::slice (1, 2, 4));
  ublas::vector_range<ublas::vector_slice<mr_t> >
    mr4sr (mr4s, ublas::range (1, 3)); 
  atlas::scal (0.1, mr4sr); 
  print_m (m, "0.1 m[4,.][1:2:4][1..3]"); 
  cout << endl; 

  // new initialization
  init_m (m, times_plus<double> (10, 1, 1)); 
#ifndef F_USE_DETAIL
  for (int i = 0; i < m.size1(); ++i) {
    mr_t mri (m, i); 
    atlas::scal (0.1, mri);
  }
#else
  // cblas level 1 function applied to matrix
  atlas::detail::scal (traits::matrix_storage_size (m),
                       0.1, traits::matrix_storage (m), 1); 
#endif 
  matr_t const cm (m); 
  print_m (cm, "new m, cm == const m"); 
  cout << endl; 

  // m[2..6][1..8]
  mrng_t mrng (m, ublas::range (2, 6), ublas::range (1, 8)); 
  print_m (mrng, "mrng = m[2..6][1..8]"); 
  cout << endl; 
  
  // mrng[1,.] <-> mrng[2,.]
  ublas::matrix_row<mrng_t> mrngr1 (mrng, 1); 
  ublas::matrix_row<mrng_t> mrngr2 (mrng, 2); 
  atlas::swap (mrngr1, mrngr2); 
  print_m (m, "mrng[1,.] <-> mrng[2,.]"); 
  cout << endl; 

  // mrng[2,.] <-> mrng[1,.]
  atlas::swap (mrngr2, mrngr1); 
  print_m (m, "mrng[2,.] <-> mrng[1,.]"); 
  cout << endl; 

  // mrng[.,3] <- 0.01 mrng[.,5] + mrng[.,3]
  ublas::matrix_column<mrng_t> mrngc3 (mrng, 3);
  ublas::matrix_column<mrng_t> const mrngc5 (mrng, 5);
  atlas::axpy (0.01, mrngc5, mrngc3); 
  print_m (m, "mrng[.,3] <- 0.01 mrng[.,5] + mrng[.,3]"); 
  cout << endl; 

  // cm[1:2:3][2:3:2] 
  cmsl_t msl (cm, ublas::slice (1, 2, 3), ublas::slice (2, 3, 2)); 
  print_m (msl, "cmsl = cm[1:2:3][2:3:2]"); 
  cout << endl; 

  // ||cmsl[.,0]||_1
  ublas::matrix_column<cmsl_t> mslc0 (msl, 0); 
  cout << "||cmsl[.,0]||_1 = " << atlas::asum (mslc0) << endl;
  cout << endl; 

  // mrng[.,3][1..4] <= 0.0001 cmsl[.,0] + mrng[.,3][1..4]
  ublas::vector_range<ublas::matrix_column<mrng_t> >
    vrmrngc3 (mrngc3, ublas::range (1, 4));
  atlas::axpy (0.0001, mslc0, vrmrngc3); 
  print_m (m, "mrng[.,3][1..4] <= 0.0001 cmsl[.,0] + mrng[.,3][1..4]"); 
  cout << endl; 

}
