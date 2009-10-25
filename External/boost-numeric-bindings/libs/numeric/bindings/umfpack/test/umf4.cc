/* ========================================================================= */
/* === umf4 ================================================================ */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.     */
/* Davis.  All Rights Reserved.  See ../README for License.                  */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.           */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                      */
/* ------------------------------------------------------------------------- */


/***********************************************************************/
/*         UMFPACK Copyright, License and Availability                 */
/***********************************************************************/
/*
 *
 * UMFPACK Version 4.1 (Apr. 30, 2003),  Copyright (c) 2003 by Timothy A.
 * Davis.  All Rights Reserved.
 *
 * UMFPACK License:
 *
 *   Your use or distribution of UMFPACK or any modified version of
 *   UMFPACK implies that you agree to this License.
 *
 *   THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
 *   EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
 *
 *   Permission is hereby granted to use or copy this program, provided
 *   that the Copyright, this License, and the Availability of the original
 *   version is retained on all copies.  User documentation of any code that
 *   uses UMFPACK or any modified version of UMFPACK code must cite the
 *   Copyright, this License, the Availability note, and "Used by permission."
 *   Permission to modify the code and to distribute modified code is granted,
 *   provided the Copyright, this License, and the Availability note are
 *   retained, and a notice that the code was modified is included.  This
 *   software was developed with support from the National Science Foundation,
 *   and is provided to you free of charge.
 *
 * Availability:
 *
 *   http://www.cise.ufl.edu/research/sparse/umfpack
 *
 */

/* Used by permission. */ 


/* modified by Kresimir Fresl, 2003              */

/* UMFPACK bindings 
 * ublas::compressed_matrix<> or ublas::coordinate_matrix<>
 * ... to use coordinate_matrix<> define COORDINATE
 */


/* Demo program for UMFPACK v4.1.  
 *
 * Reads in a triplet-form matrix and right-hand side vector. 
 * Input file format:
 *   num_rows_of_A  num_cols_of_A
 *   num_nonzeros_in_A
 *   row col val
 *   row col val
 *   ...
 *   num_nonzeros_in_B
 *   idx val
 *   ...
 *
 * Then calls UMFPACK to analyze, factor, and solve the system.
 *
 * Syntax:
 *
 *	umf4		default "auto" strategy, 1-norm row scaling
 *	umf4 a		default "auto" strategy, 1-norm row scaling
 *	umf4 u		unsymmetric strategy, 1-norm row scaling
 *	umf4 s		symmetric strategy, 1-norm row scaling
 *	umf4 2		2-by-2 strategy, maxnorm row scaling
 *	umf4 A		default "auto" strategy, maxnorm row scaling
 *	umf4 U		unsymmetric strategy, maxnorm row scaling
 *	umf4 S		symmetric strategy, maxnorm row scaling
 *	umf4 T		2-by-2 strategy , maxnorm row scaling
 *      umf4 ?n         ? can be in [aus2AUST], n: no aggressive absorption 
 */


#include <iostream>
#include <fstream> 
#include <cstdlib>
#include <string> 
#include <algorithm> 
#include <math.h>
#ifdef __ICC
#  include <mathimf.h> 
#endif  
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/bindings/umfpack/umfpack.hpp>

using std::max;
using std::cout;
using std::cin;
using std::endl; 
using std::string; 
using std::ifstream;
using std::ofstream;
using std::exit; 

namespace ublas = boost::numeric::ublas; 
namespace umf = boost::numeric::bindings::umfpack; 
namespace traits = boost::numeric::bindings::traits; 

#ifndef COORDINATE
typedef ublas::compressed_matrix<double, ublas::column_major, 0,
  ublas::unbounded_array<int>, ublas::unbounded_array<double> > cm_t; 
#else
typedef ublas::coordinate_matrix<double, ublas::column_major, 0,
  ublas::unbounded_array<int>, ublas::unbounded_array<double> > cm_t; 
#endif 
typedef ublas::vector<double> v_t; 


// wait for 'y'  

void wait4y() {
  char yes = 'n'; 
  while (1) {
    cout << "Continue (y/n) -> "; 
    cin >> yes; 
    if (yes == 'y') break;
    if (yes == 'n') exit (0);
  }
} 

// resid: compute the relative residual, ||Ax-b||/||b|| 

template <typename M, typename V> 
double resid (M const& m, V const& x, V const& b, V& r) {

  int const* Ap = traits::spmatrix_index1_storage (m);
  int const* Ai = traits::spmatrix_index2_storage (m);
  double const* Ax = traits::spmatrix_value_storage (m);

  double rnorm, bnorm, absr, absb;

  int n = traits::vector_size (r); 

  r = prod (m, x); 

#if 0
  if (transpose) {
    int i; 
    for (int j = 0; j < n; j++) 
      for (int p = Ap [j]; p < Ap [j+1]; p++) {
        i = Ai [p];
        r [j] += Ax [p] * x [i];
      }
  } else {
    int i; 
    for (int j = 0; j < n; j++)
      for (int p = Ap [j]; p < Ap [j+1]; p++) {
        i = Ai [p];
        r [i] += Ax [p] * x [j];
      }
  }
#endif 

  

  for (int i = 0; i < n; i++)
    r [i] -= b [i];

  rnorm = 0.;
  bnorm = 0.;
  for (int i = 0; i < n; i++) {
#ifndef NO_NAN
    if (isnan (r [i])){
      rnorm = r [i];
      break;
    }
#endif 
    absr = fabs (r [i]);
    rnorm = max (rnorm, absr);
  }
  for (int i = 0; i < n; i++) {
#ifndef NO_NAN
    if (isnan (b [i])){
      bnorm = b [i];
      break;
    }
#endif 
    absb = fabs (b [i]);
    bnorm = max (bnorm, absb);
  }
  if (bnorm == 0)
    bnorm = 1;
  return (rnorm / bnorm);
}


////////////////////////////////////////////////////////////////
// main program                                                       

int main (int argc, char **argv) {

  cout << "\n===========================================================\n"
       << "=== UMFPACK v4.1 ==========================================\n"
       << "===========================================================\n"
       << endl; 

  // set controls                                                      

  umf::control_type<> Control; 
  Control [UMFPACK_PRL] = 3;
  Control [UMFPACK_BLOCK_SIZE] = 32;

  if (argc > 1) {
    char *s = argv [1];
    // get the strategy 
    if (s [0] == 'u') 
      Control [UMFPACK_STRATEGY] = UMFPACK_STRATEGY_UNSYMMETRIC;
    else if (s [0] == 'a')
      Control [UMFPACK_STRATEGY] = UMFPACK_STRATEGY_AUTO;
    else if (s [0] == 's')
      Control [UMFPACK_STRATEGY] = UMFPACK_STRATEGY_SYMMETRIC;
    else if (s [0] == '2')
      Control [UMFPACK_STRATEGY] = UMFPACK_STRATEGY_2BY2;
    else if (s [0] == 'U') {
      Control [UMFPACK_STRATEGY] = UMFPACK_STRATEGY_UNSYMMETRIC;
      Control [UMFPACK_SCALE] = UMFPACK_SCALE_MAX;
    }
    else if (s [0] == 'A') {
      Control [UMFPACK_STRATEGY] = UMFPACK_STRATEGY_AUTO;
      Control [UMFPACK_SCALE] = UMFPACK_SCALE_MAX;
    }
    else if (s [0] == 'S') {
      Control [UMFPACK_STRATEGY] = UMFPACK_STRATEGY_SYMMETRIC;
      Control [UMFPACK_SCALE] = UMFPACK_SCALE_MAX;
    }
    else if (s [0] == 'T') {
      Control [UMFPACK_STRATEGY] = UMFPACK_STRATEGY_2BY2;
      Control [UMFPACK_SCALE] = UMFPACK_SCALE_MAX;
    }
    else
      printf ("unrecognized strategy: %s\n", argv [1]);

    if (s [1] == 'n')
      // no aggressive absorption 
      Control [UMFPACK_AGGRESSIVE] = 0;
  }
  bool call_wait4y = true;
  string ifs_name; 
  if (argc > 2) {
    call_wait4y = false;
    ifs_name = argv [2];
  }

  umf::report_control (Control);

  // open the matrix file 
  if (ifs_name.empty()) {
  cout << "input file -> ";
  cin >> ifs_name; 
  }
  ifstream f (ifs_name.c_str()); 
  if (!f) {
    cout << "unable to open file" << endl; 
    exit (1);
  }

  // initialize matrix A 
  cout << "reading A" << endl; 
  int nrow, ncol, nz; 
  f >> nrow >> ncol >> nz; 
  cm_t A (nrow, ncol, nz); 
  int r, c;
  for (int i = 0; i < nz; ++i) {
    f >> r >> c;
    double val;
    f >> val;
    A (r, c) = val;
  }
  Control[UMFPACK_PRL] = (nz > 20) ? 4 : 5; 
  cout << "matrix A: "; 
  umf::report_matrix (A, Control);

  if (call_wait4y) wait4y(); 

  // symbolic factorization 
  umf::info_type<> Info; 
  umf::symbolic_type<> Symbolic; 

  int status = umf::symbolic (A, Symbolic, Control, Info);

  cout << "symbolic status:" << endl; 
  umf::report_status (Control, status);
  if (status != UMFPACK_OK) {
    umf::report_info (Control, Info);
    cout << "umf::symbolic failed" << endl; 
    exit (1);
  }
  umf::report_symbolic (Symbolic, Control);

  if (call_wait4y) wait4y(); 

  // numeric factorization 
  umf::numeric_type<> Numeric; 

  status = umf::numeric (A, Symbolic, Numeric, Control, Info);

  cout << "numeric status:" << endl; 
  umf::report_status (Control, status);
  if (status < UMFPACK_OK) {
    umf::report_info (Control, Info);
    cout << "umf::numeric failed" << endl; 
    exit (1);
  }
  umf::report_numeric (Numeric, Control);

  if (call_wait4y) wait4y(); 

  // solve Ax=b 
  if (nrow == ncol && status == UMFPACK_OK) {
    
    // create right-hand side vector B 
    v_t B (ncol), X (ncol), R (ncol); 
    cout << "reading B" << endl; 
    int m2, ii; 
    f >> m2; 
    for (int i = 0; i < m2; ++i) {
      f >> ii; 
      f >> B[ii];
    }
    cout << "vector B: "; 
    umf::report_vector (B, Control); 

    status = umf::solve (A, X, B, Numeric, Control, Info);

    cout << "solve status:" << endl; 
    umf::report_status (Control, status);
    umf::report_info (Control, Info);
    if (status < UMFPACK_OK) {
      cout << "umf::solve failed" << endl;
      exit (1);
    }
    cout << "solution vector X: "; 
    umf::report_vector (X, Control);

    cout << "relative maxnorm of residual, ||Ax-b||/||b||: "
	 << resid (A, X, B, R) << endl; 

  } else {

    cout << "system not solved" << endl; 
    umf::report_info (Control, Info);

  }

  cout << endl; 

}
