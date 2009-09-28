
// BLAS level 3
// hermitian matrices, her2k 

#include <stddef.h>
#include <iostream>
#include <complex>
#include <boost/numeric/bindings/atlas/cblas3.hpp>
#include <boost/numeric/bindings/traits/ublas_hermitian.hpp>
#include "utils.h"

namespace ublas = boost::numeric::ublas;
namespace atlas = boost::numeric::bindings::atlas;
namespace traits = boost::numeric::bindings::traits;

using std::cout;
using std::cin;
using std::endl; 

typedef float real_t; 
typedef std::complex<real_t> cmplx_t; 

typedef ublas::matrix<cmplx_t, ublas::column_major> ccm_t;
typedef ublas::matrix<cmplx_t, ublas::row_major> crm_t;
typedef ublas::hermitian_adaptor<ccm_t, ublas::upper> cucha_t; 
typedef ublas::hermitian_adaptor<ccm_t, ublas::lower> clcha_t; 
typedef ublas::hermitian_adaptor<crm_t, ublas::upper> curha_t; 
typedef ublas::hermitian_adaptor<crm_t, ublas::lower> clrha_t; 

int main() {

  // complex 

  const int n1 = 3;
  const int k1 = 2; 

  ccm_t cac (n1, k1); 
  crm_t car (n1, k1); 
  cac(0,0) = car(0,0) = cmplx_t (1., 1.);
  cac(1,0) = car(1,0) = cmplx_t (2., 1.);
  cac(2,0) = car(2,0) = cmplx_t (3., 1.);
  cac(0,1) = car(0,1) = cmplx_t (1., 1.);
  cac(1,1) = car(1,1) = cmplx_t (2., 1.);
  cac(2,1) = car(2,1) = cmplx_t (3., 1.);
  print_m (cac, "cac"); 
  cout << endl; 
  print_m (car, "car"); 
  cout << endl << endl;

  ccm_t cbc (n1, k1); 
  crm_t cbr (n1, k1); 
  cbc(0,0) = cbr(0,0) = cmplx_t (0., -1.);
  cbc(1,0) = cbr(1,0) = cmplx_t (0., -1.);
  cbc(2,0) = cbr(2,0) = cmplx_t (0., -1.);
  cbc(0,1) = cbr(0,1) = cmplx_t (0., -1.);
  cbc(1,1) = cbr(1,1) = cmplx_t (0., -1.);
  cbc(2,1) = cbr(2,1) = cmplx_t (0., -1.);
  print_m (cbc, "cbc"); 
  cout << endl; 
  print_m (cbr, "cbr"); 
  cout << endl << endl;

  ccm_t ccmu (n1, n1); 
  ccm_t ccml (n1, n1); 
  crm_t crmu (n1, n1); 
  crm_t crml (n1, n1); 
  cucha_t cucha (ccmu); 
  clcha_t clcha (ccml); 
  curha_t curha (crmu); 
  clrha_t clrha (crml); 

  atlas::her2k (CblasNoTrans, cac, cbc, cucha); 
  atlas::her2k (CblasNoTrans, cmplx_t(1,0), cac, cbc, 0., clcha); 
  atlas::her2k (CblasNoTrans, cmplx_t(1,0), car, cbr, 0., curha); 
  atlas::her2k (CblasNoTrans, car, cbr, clrha); 

  print_m (cucha, "cucha");
  cout << endl; 
  print_m (clcha, "clcha");
  cout << endl; 
  print_m (curha, "curha");
  cout << endl; 
  print_m (clrha, "clrha");
  cout << endl << endl; 

  // part 2

  ccm_t cact (ublas::herm (cac)); 
  crm_t cart (ublas::herm (car)); 
  print_m (cact, "cact"); 
  cout << endl; 
  print_m (cart, "cart"); 
  cout << endl << endl;

  ccm_t cbct (ublas::herm (cbc)); 
  crm_t cbrt (ublas::herm (cbr)); 
  print_m (cbct, "cbct"); 
  cout << endl; 
  print_m (cbrt, "cbrt"); 
  cout << endl << endl;

  init_m (ccmu, const_val<cmplx_t> (cmplx_t (0, 0)));
  init_m (ccml, const_val<cmplx_t> (cmplx_t (0, 0)));
  init_m (crmu, const_val<cmplx_t> (cmplx_t (0, 0)));
  init_m (crml, const_val<cmplx_t> (cmplx_t (0, 0)));

  atlas::her2k (CblasUpper, CblasConjTrans, 1.0, cact, cbct, 0.0, ccmu); 
  atlas::her2k (CblasLower, CblasConjTrans, 1.0, cact, cbct, 0.0, ccml); 
  atlas::her2k (CblasUpper, CblasConjTrans, 1.0, cart, cbrt, 0.0, crmu); 
  atlas::her2k (CblasLower, CblasConjTrans, 1.0, cart, cbrt, 0.0, crml); 

  print_m (ccmu, "ccmu");
  cout << endl; 
  print_m (ccml, "ccml");
  cout << endl; 
  print_m (crmu, "crmu");
  cout << endl; 
  print_m (crml, "crml");
  cout << endl; 
}

