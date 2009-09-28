
// BLAS level 3
// symmetric matrices, syrk 

#include <stddef.h>
#include <iostream>
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

typedef ublas::matrix<real_t, ublas::column_major> cm_t;
typedef ublas::matrix<real_t, ublas::row_major> rm_t;
typedef ublas::symmetric_adaptor<cm_t, ublas::upper> ucsa_t; 
typedef ublas::symmetric_adaptor<cm_t, ublas::lower> lcsa_t; 
typedef ublas::symmetric_adaptor<rm_t, ublas::upper> ursa_t; 
typedef ublas::symmetric_adaptor<rm_t, ublas::lower> lrsa_t; 

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

  cm_t cmu (n, n); 
  cm_t cml (n, n); 
  rm_t rmu (n, n); 
  rm_t rml (n, n); 
  ucsa_t ucsa (cmu); 
  lcsa_t lcsa (cml); 
  ursa_t ursa (rmu); 
  lrsa_t lrsa (rml); 

  atlas::syrk (CblasNoTrans, ac, ucsa); 
  atlas::syrk (CblasNoTrans, 1.0, ac, 0.0, lcsa); 
  atlas::syrk (CblasNoTrans, 1.0, ar, 0.0, ursa); 
  atlas::syrk (CblasNoTrans, ar, lrsa); 

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

  init_m (cmu, const_val<real_t> (0));
  init_m (cml, const_val<real_t> (0));
  init_m (rmu, const_val<real_t> (0));
  init_m (rml, const_val<real_t> (0));

  atlas::syrk (CblasUpper, CblasTrans, 1.0, act, 0.0, cmu); 
  atlas::syrk (CblasLower, CblasTrans, 1.0, act, 0.0, cml); 
  atlas::syrk (CblasUpper, CblasTrans, 1.0, art, 0.0, rmu); 
  atlas::syrk (CblasLower, CblasTrans, 1.0, art, 0.0, rml); 

  print_m (cmu, "cmu");
  cout << endl; 
  print_m (cml, "cml");
  cout << endl; 
  print_m (rmu, "rmu");
  cout << endl; 
  print_m (rml, "rml");
  cout << endl; 
}

