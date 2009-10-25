/* ====================================================================== */
/* === umfpack_zi_demo ================================================== */
/* ====================================================================== */


/* ---------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.  */
/* Davis.  All Rights Reserved.  See ../README for License.               */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.        */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                   */
/* ---------------------------------------------------------------------- */


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

/* UMFPACK bindings & ublas::compressed_matrix<> */


/*
  A demo of UMFPACK:   umfpack_zi_* version.

  First, factor and solve a 5-by-5 system, Ax=b, using default parameters.
  Then solve A'x=b using the factors of A.   Modify one entry (A (1,4) = 0,
  where the row and column indices range from 0 to 4.  The pattern of A
  has not changed (it has explicitly zero entry), so a reanalysis with
  umfpack_zi_symbolic does not need to be done.  Refactorize (with
  umfpack_zi_numeric), and solve Ax=b.  Note that the pivot ordering has
  changed.  Next, change all of the entries in A, but not the pattern.

  Finally, compute C = A', and do the symbolic and numeric factorization of C.
  Factorizing A' can sometimes be better than factorizing A itself (less work
  and memory usage).  Solve C'x=b; the solution is the same as the
  solution to Ax=b.
*/


#include <iostream>
#include <fstream> 
#include <cstdlib>
#include <algorithm> 
#include <complex> 
#include <math.h>
#include <boost/mpl/and.hpp>
#include <boost/numeric/bindings/traits/c_array.hpp>
#include <boost/numeric/bindings/traits/std_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/bindings/umfpack/umfpack.hpp>

using std::max;
using std::cout;
using std::cin;
using std::endl; 
using std::exit; 

namespace ublas = boost::numeric::ublas; 
namespace umf = boost::numeric::bindings::umfpack; 
namespace traits = boost::numeric::bindings::traits; 

typedef ublas::compressed_matrix<
  std::complex<double>, 
  ublas::column_major, 
  0,
  ublas::unbounded_array<int>, 
  ublas::unbounded_array<std::complex<double> > 
> cm_t; 

typedef std::vector<std::complex<double> > v_t; 


/* --------------------------------------------------------------- */
/* wait for 'y'                                                    */
/* --------------------------------------------------------------- */

void wait4y() {
  char yes = 'n'; 
  while (1) {
    cout << "Continue (y/n) -> "; 
    cin >> yes; 
    if (yes == 'y') break;
    if (yes == 'n') exit (0);
  }
} 


/* --------------------------------------------------------------- */
/* error: print a message and exit                                 */
/* --------------------------------------------------------------- */

void error (char *message) {
  cout << "\n\n====== error: " << message << " =====\n" << endl; 
  exit (1);
}


/* ------------------------------------------------------------------------- */
/* resid: compute the residual, r = Ax-b or r = A'x=b and return maxnorm (r) */
/* A' is the complex conjugate transpose, not the array transpose            */
/* note: works only with compressed column matrices                          */
/* ------------------------------------------------------------------------- */

template <typename M, typename V> 
double resid (M const& A, V const& x, V const& b, V& r, int transpose = 0) {

  double norm;
  int const* Ap = traits::spmatrix_index1_storage (A);
  int const* Ai = traits::spmatrix_index2_storage (A);
  std::complex<double> const* Ax = traits::spmatrix_value_storage (A);

  int n = traits::vector_size (r); 

  for (int i = 0; i < n; i++)
    r [i] = -b [i];

  if (transpose) {
    int i; 
    for (int j = 0; j < n; j++)
      for (int p = Ap [j]; p < Ap [j+1]; p++) {
        i = Ai [p];
        r [j] += std::conj (Ax [p]) * x [i];
      }
  } else {
    int i; 
    for (int j = 0; j < n; j++) 
      for (int p = Ap [j]; p < Ap [j+1]; p++) {
        i = Ai [p];
        r [i] += Ax [p] * x [j];
      }
  }

  norm = 0. ;
  for (int i = 0 ; i < n ; i++)
    norm = max (std::abs (r [i]), norm);
  return (norm) ;
}

/* ------------------------------------------------------------------- */
/* main program                                                        */
/* ------------------------------------------------------------------- */

int main (int argc, char **argv) {
  bool call_wait4y = true;
  if (argc > 1) {
    char *s = argv [1];
    if (s [0] == '-' && s [1] == 'y') {
      call_wait4y = false;
    }
  }

  /* ---------------------------------------------------------------- */
  /* initialisations                                                  */
  /* ---------------------------------------------------------------- */

  cout << "\n" << UMFPACK_VERSION << " demo: _di_ version\n" << endl; 

  // get the default control parameters
  umf::control_type<std::complex<double> > Control; 

  // change the default print level for this demo 
  // (otherwise, nothing will print) 
  Control [UMFPACK_PRL] = 6;
  // print the license agreement 
  umf::report_status (Control, UMFPACK_OK);
  Control [UMFPACK_PRL] = 5;
  // print the control parameters 
  umf::report_control (Control);

  umf::info_type<std::complex<double> > Info;

  // matrix A                                         

  int    n = 5, nz = 12 ;
  int    Arow [ ] = { 0,  4,  1,  1,   2,   2,  0,  1,  2,  3,  4,  4} ;
  int    Acol [ ] = { 0,  4,  0,  2,   1,   2,  1,  4,  3,  2,  1,  2} ;
  double Aval [ ] = {2., 1., 3., 4., -1., -3., 3., 6., 2., 1., 4., 2.} ;
  double Avalz[ ] = {1., .4, .1, .2, -1., -.2, 0., 6., 3., 0., .3, .3} ;

  cm_t A (n, n, nz); 
  for (int i = 0; i < nz; ++i) 
    A (Arow[i], Acol[i]) = std::complex<double> (Aval[i], Avalz[i]); 

  cout << "matrix A: ";
  umf::report_matrix (A, Control);
    
  // right-hand-side etc. 

  double bx[ ] = {8., 45., -3., 3., 19.};
  double bz[ ] = {1., -5., -2., 0., 2.2}; 

  v_t b (n), x (n), r (n); 
  for (int i = 0; i < n; ++i) 
    b[i] = std::complex<double> (bx[i], bz[i]); 

  cout << "vector B: "; 
  umf::report_vector (b, Control); 

  /* ---------------------------------------------------------------- */
  /* symbolic factorisation                                           */
  /* ---------------------------------------------------------------- */

  if (call_wait4y) wait4y(); 

  cout << "\nSymbolic factorisation:" << endl; 

  umf::symbolic_type<std::complex<double> > Symbolic; 

  int status = umf::symbolic (A, Symbolic, Control, Info);

  umf::report_status (Control, status);
  if (status < 0) {
    umf::report_info (Control, Info);
    error ("umf::symbolic failed");
  }
  umf::report_symbolic (Symbolic, Control);

  /* ---------------------------------------------------------------- */
  /* numeric factorisation                                            */
  /* ---------------------------------------------------------------- */

  if (call_wait4y) wait4y(); 

  cout << "\nNumeric factorisation:" << endl; 

  umf::numeric_type<std::complex<double> > Numeric; 

  status = umf::numeric (A, Symbolic, Numeric,Control, Info);

  umf::report_status (Control, status);
  if (status < 0) {
    umf::report_info (Control, Info);
    error ("umf::numeric failed");
  }
  umf::report_numeric (Numeric, Control);

  /* ---------------------------------------------------------------- */
  /* solve Ax=b                                                       */
  /* ---------------------------------------------------------------- */

  if (call_wait4y) wait4y(); 

  cout << "\nSolve Ax = b:" << endl; 

  status = umf::solve (UMFPACK_A, A, x, b, Numeric, Control, Info);

  umf::report_status (Control, status);
  umf::report_info (Control, Info);
  if (status < 0) 
    error ("umf::solve failed");

  cout << "\nx (solution of Ax=b): "; 
  umf::report_vector (x, Control);

  double rnorm = resid (A, x, b, r);
  cout << "\nmaxnorm of residual: " << rnorm << endl << endl; 

  /* ---------------------------------------------------------------- */
  /* solve Ax=b, broken down into steps                               */
  /* ---------------------------------------------------------------- */

  if (call_wait4y) wait4y(); 

  cout << "\nSolve Ax = b (broken down into steps):" << endl; 

  // Rb = R*b 
  std::vector<std::complex<double> > Rb (n), y (n); 
  status = umf::scale (Rb, b, Numeric);
  if (status < 0) 
    error ("umf::scale failed");
    
  // solve Ly = P*(Rb) 
  status = umf::solve (UMFPACK_Pt_L, A, y, Rb, Numeric, Control, Info);
  if (status < 0) 
    error ("first umf::solve failed");

  // solve UQ'x=y 
  status = umf::solve (UMFPACK_U_Qt, A, x, y, Numeric, Control, Info);
  if (status < 0) 
    error ("second umf::solve failed");

  cout << "\nx (solution of Ax=b, solve is split into 3 steps): ";
  umf::report_vector (x, Control);

  rnorm = resid (A, x, b, r);
  cout << "\nmaxnorm of residual: " << rnorm << endl << endl; 

  /* ------------------------------------------------------------------- */
  /* solve A'x=b                                                         */
  /* note that this is the complex conjugate transpose, A'               */
  /* ------------------------------------------------------------------- */

  if (call_wait4y) wait4y(); 

  cout << "\nSolve A'x = b:" << endl; 

  status = umf::solve (UMFPACK_At, A, x, b, Numeric, Control, Info);

  umf::report_status (Control, status);
  umf::report_info (Control, Info);
  if (status < 0) 
    error ("umf::solve failed");

  cout << "\nx (solution of A'x=b): "; 
  umf::report_vector (x, Control);

  rnorm = resid (A, x, b, r, 1); // '1' -- transpose 
  cout << "\nmaxnorm of residual: " << rnorm << endl << endl; 

  /* ----------------------------------------------------------------- */
  /* modify one numerical value in the column-form of A                */
  /* ----------------------------------------------------------------- */

  if (call_wait4y) wait4y(); 

  cout << "\nModify A (1,4) if present:" << endl; 

  int row = 1, col = 4;
  cm_t::iterator1 i1 = A.begin1(); 
  ++i1; // go to row 1
  cm_t::iterator2 i2 = i1.begin();
  for (; i2 != i1.end(); ++i2) 
    if (i2.index2() == 4) break; 
  if (i2 != i1.end()) {
    if (i2.index1() != 1 || i2.index2() != 4)
      error ("something\'s very wrong"); 
    cout << "\nchanging A (1,4) to zero" << endl; 
    A (1, 4) = 0.0; 
  }

  cout << "\nmodified A: ";
  umf::report_matrix (A, Control);

  /* ----------------------------------------------------------------- */
  /* redo the numeric factorization                                    */
  /* ----------------------------------------------------------------- */

  if (call_wait4y) wait4y(); 

  cout << "\nRedo the numeric factorization:" << endl; 

  // The pattern (Ap and Ai) hasn't changed, so the symbolic factorisation 
  // doesn't have to be redone, no matter how much we change Ax. 

  // We don't need the Numeric object any more, so free it. 
  umf::free (Numeric);  // or: Numeric.free(); 
  // Note that a memory leak would have occured if the old Numeric 
  // had not been free'd.

  status = umf::numeric (A, Symbolic, Numeric, Control, Info);

  umf::report_status (Control, status);
  if (status < 0) {
    umf::report_info (Control, Info);
    error ("umf::numeric failed");
  }

  cout << "\nNumeric factorization of modified A: "; 
  umf::report_numeric (Numeric, Control);

  /* ------------------------------------------------------------------ */
  /* solve Ax=b, with the modified A                                    */
  /* ------------------------------------------------------------------ */

  if (call_wait4y) wait4y(); 

  cout << "\nSolve Ax = b with the modified A:" << endl; 
  
  status = umf::solve (A, x, b, Numeric, Control, Info);

  umf::report_status (Control, status);
  umf::report_info (Control, Info);
  if (status < 0) 
    error ("umf::solve failed");

  cout << "\nx (with modified A): ";
  umf::report_vector (x, Control);

  rnorm = resid (A, x, b, r);
  cout << "\nmaxnorm of residual: " << rnorm << endl << endl; 


  /* ------------------------------------------------------------------ */
  /* modify all of the numerical values of A, but not the pattern       */
  /* ------------------------------------------------------------------ */

  if (call_wait4y) wait4y(); 

  cout << "\nModify numerical values of A:" << endl; 

  for (i2 = A.begin2(); i2 != A.end2(); ++i2)
    for (i1 = i2.begin(); i1 != i2.end(); ++i1) {
      std::complex<double> c = *i1; 
      cout << "changing A (" 
           << i1.index1() << "," << i1.index2() << ") from " << c;
      double a = std::real (c) + 10 * i1.index2() - i1.index1(); 
      double ai = std::imag (c); 
      c = std::complex<double> (a, ai); 
      cout << " to " << c << endl; 
      *i1 = c;
    }
  cout << endl; 

  cout << "\ncompletely modified A (same pattern): ";
  umf::report_matrix (A, Control);

  /* ------------------------------------------------------------------ */
  /* redo the numeric factorization                                     */
  /* ------------------------------------------------------------------ */

  if (call_wait4y) wait4y(); 

  cout << "\nRedo the numeric factorization:" << endl; 

  umf::free (Numeric);
  status = umf::numeric (A, Symbolic, Numeric, Control, Info);
  umf::report_status (Control, status);
  if (status < 0) {
    umf::report_info (Control, Info);
    error ("umf::numeric failed");
  }
   
  cout << "\nNumeric factorization of completely modified A: ";
  umf::report_numeric (Numeric, Control);

  /* ------------------------------------------------------------------ */
  /* solve Ax=b, with the modified A                                    */
  /* ------------------------------------------------------------------ */

  if (call_wait4y) wait4y(); 

  cout << "\nSolve Ax = b with the modified A:" << endl; 
  
  status = umf::solve (A, x, b, Numeric, Control, Info);

  umf::report_status (Control, status);
  umf::report_info (Control, Info);
  if (status < 0) 
    error ("umf::solve failed");

  cout << "\nx (with completely modified A): ";
  umf::report_vector (x, Control);

  rnorm = resid (A, x, b, r);
  cout << "\nmaxnorm of residual: " << rnorm << endl << endl; 

  /* ---------------------------------------------------------------- */
  /* free the symbolic and numeric factorization                      */
  /* ---------------------------------------------------------------- */

  Symbolic.free();
  Numeric.free();

  /* ---------------------------------------------------------------- */
  /* C = transpose of A                                               */
  /* ---------------------------------------------------------------- */

  if (call_wait4y) wait4y(); 

  cout << "\nC = transpose of A:" << endl; 

  cm_t C (herm (A)); 

  umf::report_matrix (C, Control);

  /* ----------------------------------------------------------------- */
  /* symbolic factorization of C                                       */
  /* ----------------------------------------------------------------- */

  if (call_wait4y) wait4y(); 

  cout << "\nSymbolic factorization of C:" << endl; 

  status = umf::symbolic (C, Symbolic,Control, Info);

  umf::report_status (Control, status);
  if (status < 0) {
    umf::report_info (Control, Info);
    error ("umf::symbolic failed");
  }
  umf::report_symbolic (Symbolic, Control);

  /* ----------------------------------------------------------------- */
  /* numeric factorization of C                                        */
  /* ----------------------------------------------------------------- */

  if (call_wait4y) wait4y(); 

  cout << "\nNumeric factorization of C:" << endl; 

  status = umf::numeric (C, Symbolic, Numeric,Control, Info);

  umf::report_status (Control, status);
  if (status < 0) {
    umf::report_info (Control, Info);
    error ("umf::numeric failed");
  }
  umf::report_numeric (Numeric, Control);

  /* ----------------------------------------------------------------- */
  /* solve C'x=b                                                       */
  /* ----------------------------------------------------------------- */

  if (call_wait4y) wait4y(); 

  cout << "\nSolve C'x = b:" << endl; 

  status = umf::solve (UMFPACK_At, C, x, b, Numeric, Control, Info);
  umf::report_status (Control, status);
  umf::report_info (Control, Info);
  if (status < 0)
    error ("umf::solve failed");

  cout << "\nx (solution of C'x=b): ";
  umf::report_vector (x, Control);

  rnorm = resid (C, x, b, r, 1);
  cout << "\nmaxnorm of residual: " << rnorm << endl << endl; 
}
