
// BLAS level 2
// symmetric matrices

//#define BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 
//#define BOOST_NO_FUNCTION_TEMPLATE_ORDERING

#include <stddef.h>
#include <iostream>
#include <boost/numeric/bindings/atlas/cblas1.hpp>
#include <boost/numeric/bindings/atlas/cblas2.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_symmetric.hpp>
#include "utils.h"

namespace ublas = boost::numeric::ublas;
namespace atlas = boost::numeric::bindings::atlas;
namespace traits = boost::numeric::bindings::traits;

using std::cout;
using std::cin;
using std::endl; 

typedef double real_t; 
typedef ublas::vector<real_t> vct_t;

typedef ublas::symmetric_matrix<
  real_t, ublas::upper, ublas::column_major
> ucsymm_t; 
typedef ublas::symmetric_matrix<
  real_t, ublas::lower, ublas::column_major
> lcsymm_t; 
typedef ublas::symmetric_matrix<
  real_t, ublas::upper, ublas::row_major
> ursymm_t; 
typedef ublas::symmetric_matrix<
  real_t, ublas::lower, ublas::row_major
> lrsymm_t; 

typedef ublas::matrix<real_t, ublas::column_major> cm_t;
typedef ublas::matrix<real_t, ublas::row_major> rm_t;
typedef ublas::symmetric_adaptor<cm_t, ublas::upper> ucsa_t; 
typedef ublas::symmetric_adaptor<cm_t, ublas::lower> lcsa_t; 
typedef ublas::symmetric_adaptor<rm_t, ublas::upper> ursa_t; 
typedef ublas::symmetric_adaptor<rm_t, ublas::lower> lrsa_t; 

int main() {
  
  cout << endl; 
  size_t n; 
  cout << "n -> "; 
  cin >> n;
  cout << endl; 

  vct_t vx (n), vy (n); 
  atlas::set (1., vx); 
  print_v (vx, "vx"); 
  vy (1) = 1.; 
  print_v (vy, "vy"); 

  // symmetric matrix

  ucsymm_t ucs (n, n); 
  lcsymm_t lcs (n, n); 
  ursymm_t urs (n, n); 
  lrsymm_t lrs (n, n); 

  init_symm (ucs, 'u'); 
  init_symm (lcs, 'l'); 
  init_symm (urs, 'u'); 
  init_symm (lrs, 'l'); 

  print_m (ucs, "ucs");
  cout << endl; 
  print_m (lcs, "lcs");
  cout << endl; 
  print_m (urs, "urs");
  cout << endl; 
  print_m (lrs, "lrs");
  cout << endl; 

  // m += x x^T 
  atlas::spr (vx, ucs); 
  print_m (ucs, "ucs += x x^T"); 
  cout << endl; 
  atlas::spr (1, vx, lcs); 
  print_m (lcs, "lcs += x x^T"); 
  cout << endl; 
  atlas::spr (1.0f, vx, urs); 
  print_m (urs, "urs += x x^T"); 
  cout << endl; 
  atlas::spr (vx, lrs); 
  print_m (lrs, "lrs += x x^T"); 
  cout << endl; 

  // m += x y^T + y x^T
  atlas::spr2 (1, vx, vy, ucs); 
  print_m (ucs, "ucs += x y^T + y x^T"); 
  cout << endl; 
  atlas::spr2 (vx, vy, lcs); 
  print_m (lcs, "lcs += x y^T + y x^T"); 
  cout << endl; 
  atlas::spr2 (vx, vy, urs); 
  print_m (urs, "urs += x y^T + y x^T"); 
  cout << endl; 
  atlas::spr2 (1., vx, vy, lrs); 
  print_m (lrs, "lrs += x y^T + y x^T"); 
  cout << endl; 


  // symmetric adaptor

  cm_t cmu (n, n); 
  cm_t cml (n, n); 
  rm_t rmu (n, n); 
  rm_t rml (n, n); 

  ucsa_t ucsa (cmu); 
  lcsa_t lcsa (cml); 
  ursa_t ursa (rmu); 
  lrsa_t lrsa (rml); 

  init_symm (ucsa, 'u'); 
  init_symm (lcsa, 'l'); 
  init_symm (ursa, 'u'); 
  init_symm (lrsa, 'l'); 

  print_m (ucsa, "ucsa");
  cout << endl; 
  print_m (lcsa, "lcsa");
  cout << endl; 
  print_m (ursa, "ursa");
  cout << endl; 
  print_m (lrsa, "lrsa");
  cout << endl; 

  // m += x x^T 
  atlas::syr (vx, ucsa); 
  print_m (ucsa, "ucsa += x x^T"); 
  cout << endl; 
  atlas::syr (1, vx, lcsa); 
  print_m (lcsa, "lcsa += x x^T"); 
  cout << endl; 
  atlas::syr (1.0f, vx, ursa); 
  print_m (ursa, "ursa += x x^T"); 
  cout << endl; 
  atlas::syr (vx, lrsa); 
  print_m (lrsa, "lrsa += x x^T"); 
  cout << endl; 

  // m += x y^T + y x^T
  atlas::syr2 (1, vx, vy, ucsa); 
  print_m (ucsa, "ucsa += x y^T + y x^T"); 
  cout << endl; 
  atlas::syr2 (vx, vy, lcsa); 
  print_m (lcsa, "lcsa += x y^T + y x^T"); 
  cout << endl; 
  atlas::syr2 (vx, vy, ursa); 
  print_m (ursa, "ursa += x y^T + y x^T"); 
  cout << endl; 
  atlas::syr2 (1., vx, vy, lrsa); 
  print_m (lrsa, "lrsa += x y^T + y x^T"); 
  cout << endl; 

}
