
// BLAS level 3
// symmetric matrices, syr2k 

#include <stddef.h>
#include <iostream>
#include <complex>
#include <boost/numeric/bindings/atlas/cblas3.hpp>
#include <boost/numeric/bindings/traits/ublas_symmetric.hpp>
#include "utils.h"

namespace ublas = boost::numeric::ublas;
namespace atlas = boost::numeric::bindings::atlas;
namespace traits = boost::numeric::bindings::traits;

using std::cout;
using std::cin;
using std::endl; 

typedef double real_t; 
typedef std::complex<real_t> cmplx_t; 

typedef ublas::matrix<real_t, ublas::column_major> cm_t;
typedef ublas::matrix<real_t, ublas::row_major> rm_t;
typedef ublas::symmetric_adaptor<cm_t, ublas::upper> ucsa_t; 
typedef ublas::symmetric_adaptor<cm_t, ublas::lower> lcsa_t; 
typedef ublas::symmetric_adaptor<rm_t, ublas::upper> ursa_t; 
typedef ublas::symmetric_adaptor<rm_t, ublas::lower> lrsa_t; 

typedef ublas::matrix<cmplx_t, ublas::column_major> ccm_t;
typedef ublas::matrix<cmplx_t, ublas::row_major> crm_t;
typedef ublas::symmetric_adaptor<ccm_t, ublas::upper> cucsa_t; 
typedef ublas::symmetric_adaptor<ccm_t, ublas::lower> clcsa_t; 
typedef ublas::symmetric_adaptor<crm_t, ublas::upper> cursa_t; 
typedef ublas::symmetric_adaptor<crm_t, ublas::lower> clrsa_t; 

int main() {

  int n, k;
  cout << "n -> ";
  cin >> n;
  cout << "k -> ";
  cin >> k; 

  cm_t ac (n, k); 
  rm_t ar (n, k); 
  init_m (ac, rws1()); 
  init_m (ar, rws1());
  print_m (ac, "ac"); 
  cout << endl; 
  print_m (ar, "ar"); 
  cout << endl << endl;

  cm_t bc (n, k); 
  rm_t br (n, k); 
  init_m (bc, const_val<real_t> (1)); 
  init_m (br, const_val<real_t> (1));
  print_m (bc, "bc"); 
  cout << endl; 
  print_m (br, "br"); 
  cout << endl << endl;

  cm_t cmu (n, n); 
  cm_t cml (n, n); 
  rm_t rmu (n, n); 
  rm_t rml (n, n); 
  ucsa_t ucsa (cmu); 
  lcsa_t lcsa (cml); 
  ursa_t ursa (rmu); 
  lrsa_t lrsa (rml); 

  atlas::syr2k (CblasNoTrans, ac, bc, ucsa); 
  atlas::syr2k (CblasNoTrans, 1.0, ac, bc, 0.0, lcsa); 
  atlas::syr2k (CblasNoTrans, 1.0, ar, br, 0.0, ursa); 
  atlas::syr2k (CblasNoTrans, ar, br, lrsa); 

  print_m (ucsa, "ucsa");
  cout << endl; 
  print_m (lcsa, "lcsa");
  cout << endl; 
  print_m (ursa, "ursa");
  cout << endl; 
  print_m (lrsa, "lrsa");
  cout << endl << endl; 

  // part 2

  cm_t act (k, n); 
  rm_t art (k, n); 
  init_m (act, cls1()); 
  init_m (art, cls1());
  print_m (act, "act"); 
  cout << endl; 
  print_m (art, "art"); 
  cout << endl << endl;

  cm_t bct (k, n); 
  rm_t brt (k, n); 
  init_m (bct, rws1()); 
  init_m (brt, rws1());
  print_m (bct, "bct"); 
  cout << endl; 
  print_m (brt, "brt"); 
  cout << endl << endl;

  init_m (cmu, const_val<real_t> (0));
  init_m (cml, const_val<real_t> (0));
  init_m (rmu, const_val<real_t> (0));
  init_m (rml, const_val<real_t> (0));

  atlas::syr2k (CblasUpper, CblasTrans, 1.0, act, bct, 0.0, cmu); 
  atlas::syr2k (CblasLower, CblasTrans, 1.0, act, bct, 0.0, cml); 
  atlas::syr2k (CblasUpper, CblasTrans, 1.0, art, brt, 0.0, rmu); 
  atlas::syr2k (CblasLower, CblasTrans, 1.0, art, brt, 0.0, rml); 

  print_m (cmu, "cmu");
  cout << endl; 
  print_m (cml, "cml");
  cout << endl; 
  print_m (rmu, "rmu");
  cout << endl; 
  print_m (rml, "rml");
  cout << endl; 

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
  cucsa_t cucsa (ccmu); 
  clcsa_t clcsa (ccml); 
  cursa_t cursa (crmu); 
  clrsa_t clrsa (crml); 

  atlas::syr2k (CblasNoTrans, cac, cbc, cucsa); 
  atlas::syr2k (CblasNoTrans, cmplx_t(1,0), cac, cbc, cmplx_t(0,0), clcsa); 
  atlas::syr2k (CblasNoTrans, cmplx_t(1,0), car, cbr, cmplx_t(0,0), cursa); 
  atlas::syr2k (CblasNoTrans, car, cbr, clrsa); 

  print_m (cucsa, "cucsa");
  cout << endl; 
  print_m (clcsa, "clcsa");
  cout << endl; 
  print_m (cursa, "cursa");
  cout << endl; 
  print_m (clrsa, "clrsa");
  cout << endl << endl; 

}

