
// solving A * X = B
// A symmetric/hermitian positive definite
// driver function posv()

// #define BOOST_UBLAS_STRICT_HERMITIAN
// .. doesn't work (yet?)  

//#define BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 
//#define BOOST_NO_FUNCTION_TEMPLATE_ORDERING

#include <cstddef>
#include <iostream>
#include <complex>
#include <boost/numeric/bindings/atlas/cblas1.hpp>
#include <boost/numeric/bindings/atlas/cblas2.hpp>
#include <boost/numeric/bindings/atlas/cblas3.hpp>
#include <boost/numeric/bindings/atlas/clapack.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_symmetric.hpp>
#include <boost/numeric/bindings/traits/ublas_hermitian.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "utils.h"

namespace ublas = boost::numeric::ublas;
namespace atlas = boost::numeric::bindings::atlas;

using std::size_t; 
using std::cout;
using std::endl; 

typedef std::complex<float> cmplx_t; 

#ifndef F_ROW_MAJOR
typedef ublas::matrix<double, ublas::column_major> m_t;
typedef ublas::matrix<cmplx_t, ublas::column_major> cm_t;
#else
typedef ublas::matrix<double, ublas::row_major> m_t;
typedef ublas::matrix<cmplx_t, ublas::row_major> cm_t;
#endif

#ifndef F_UPPER
typedef ublas::symmetric_adaptor<m_t, ublas::lower> symm_t; 
typedef ublas::hermitian_adaptor<cm_t, ublas::lower> herm_t; 
#else
typedef ublas::symmetric_adaptor<m_t, ublas::upper> symm_t; 
typedef ublas::hermitian_adaptor<cm_t, ublas::upper> herm_t; 
#endif 

int main() {

  cout << endl; 

  // symmetric 
  cout << "real symmetric\n" << endl; 

  size_t n = 5; 
  m_t a (n, n);    // matrix (storage)
  symm_t sa (a);   // symmetric adaptor 
#ifdef F_UPPER 
  init_symm (sa, 'u');
#else
  init_symm (sa, 'l');
#endif 
  // ifdef F_UPPER 
  //        [5 4 3 2 1]
  //        [0 5 4 3 2]
  //    a = [0 0 5 4 3]
  //        [0 0 0 5 4]
  //        [0 0 0 0 n]
  // else 
  //        [5 0 0 0 0]
  //        [4 5 0 0 0]
  //    a = [3 4 5 0 0]
  //        [2 3 4 5 0]
  //        [1 2 3 4 5]
  print_m (sa, "sa"); 
  cout << endl; 
  print_m (a, "a"); 
  cout << endl; 
  print_m_data (sa, "sa"); 
  cout << endl; 

  size_t nrhs = 2; 
  m_t x (n, nrhs); 
  // b -- right-hand side matrix:
  // .. see leading comments for `gesv()' in clapack.hpp
#ifndef F_ROW_MAJOR
  m_t b (n, nrhs);
#else
  m_t b (nrhs, n);
#endif
  ublas::matrix_column<m_t> xc0 (x, 0), xc1 (x, 1); 
  atlas::set (1., xc0);  // x[.,0] = 1
  atlas::set (2., xc1);  // x[.,1] = 2
#ifndef F_ROW_MAJOR
#ifndef BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 
  atlas::symm (sa, x, b);  // b = a x, so we know the result ;o)
#else
  atlas::symm (CblasLeft, 1.0, sa, x, 0.0, b); 
#endif 
#else
  // see leading comments for `gesv()' in clapack.hpp
  ublas::matrix_row<m_t> br0 (b, 0), br1 (b, 1); 
  atlas::symv (sa, xc0, br0);  // b[0,.] = a x[.,0]
  atlas::symv (sa, xc1, br1);  // b[1,.] = a x[.,1]  =>  b^T = a x
#endif 
  print_m (b, "b"); 
  cout << endl; 

  atlas::cholesky_solve (sa, b);  // same as posv() 
  print_m (b, "x"); 
  cout << endl; 

  //////////////////////////////////////////////////////////
  // hermitian 
  cout << "\n===========================\n" << endl; 
  cout << "complex hermitian (well, not really)\n" << endl; 

  cm_t ca (3, 3);   // matrix (storage)
  herm_t ha (ca);   // hermitian adaptor 
  cm_t cx (3, 1);
#ifndef F_ROW_MAJOR
  cm_t cb (3, 1);
#else
  cm_t cb (1, 3); 
#endif  

#ifdef F_UPPER
  init_symm (ha, 'u'); 
#else
  init_symm (ha, 'l'); 
#endif 
  print_m (ha, "ha"); 
  cout << endl; 
  print_m (ca, "ca"); 
  cout << endl; 
  print_m_data (ha, "ha"); 
  cout << endl; 

  ublas::matrix_column<cm_t> cx0 (cx, 0);
  atlas::set (cmplx_t (1, -1), cx0);
  print_m (cx, "cx"); 
  cout << endl; 
#ifndef F_ROW_MAJOR
  ublas::matrix_column<cm_t> cb0 (cb, 0); 
#else
  ublas::matrix_row<cm_t> cb0 (cb, 0); 
#endif
  atlas::hemv (ha, cx0, cb0); 
  print_m (cb, "cb"); 
  cout << endl; 
  
  int ierr = atlas::posv (ha, cb); 
  if (ierr == 0)
    print_m (cb, "cx"); 
  else 
    cout << "matrix is not positive definite: ierr = " 
         << ierr << endl;
  cout << endl; 

  cout << "\n===========================\n" << endl; 
  cout << "complex hermitian\n" << endl; 

  // regular (see ublas_gesv.cc), but not positive definite: 
#ifndef F_UPPER
  ha (0, 0) = cmplx_t (3, 0);
  ha (1, 0) = cmplx_t (4, -2);
  ha (1, 1) = cmplx_t (5, 0);
  ha (2, 0) = cmplx_t (-7, -5);
  ha (2, 1) = cmplx_t (0, 3);
  ha (2, 2) = cmplx_t (2, 0);
#else
  ha (0, 0) = cmplx_t (3, 0);
  ha (0, 1) = cmplx_t (4, 2);
  ha (0, 2) = cmplx_t (-7, 5);
  ha (1, 1) = cmplx_t (5, 0);
  ha (1, 2) = cmplx_t (0, -3);
  ha (2, 2) = cmplx_t (2, 0);
#endif
  print_m (ha, "ha"); 
  cout << endl; 
  print_m (ca, "ca"); 
  cout << endl; 
  print_m_data (ha, "ha"); 
  cout << endl; 

  atlas::set (cmplx_t (1, 1), cx0);
  print_m (cx, "cx"); 
  cout << endl; 
  atlas::hemv (ha, cx0, cb0); 
  print_m (cb, "cb"); 
  cout << endl; 
  
  ierr = atlas::cholesky_solve (ha, cb); 
  if (ierr == 0)
    print_m (cb, "cx"); 
  else 
    cout << "matrix is not positive definite: ierr = " 
         << ierr << endl;
  cout << endl; 

  cout << "\n===========================\n" << endl; 
  cout << "complex hermitian\n" << endl; 

  // positive definite: 
#ifndef F_UPPER
  ha (0, 0) = cmplx_t (25, 0);
  ha (1, 0) = cmplx_t (-5, 5);
  ha (1, 1) = cmplx_t (51, 0);
  ha (2, 0) = cmplx_t (10, -5);
  ha (2, 1) = cmplx_t (4, 6);
  ha (2, 2) = cmplx_t (71, 0);
#else
  ha (0, 0) = cmplx_t (25, 0);
  ha (0, 1) = cmplx_t (-5, -5);
  ha (0, 2) = cmplx_t (10, 5);
  ha (1, 1) = cmplx_t (51, 0);
  ha (1, 2) = cmplx_t (4, -6);
  ha (2, 2) = cmplx_t (71, 0);
#endif
  print_m (ha, "ha"); 
  cout << endl; 

#ifndef F_ROW_MAJOR
  cm_t cb32 (3, 2); 
  cb32 (0, 0) = cmplx_t (60, -55);
  cb32 (1, 0) = cmplx_t (34, 58);
  cb32 (2, 0) = cmplx_t (13, -152);
  cb32 (0, 1) = cmplx_t (70, 10);
  cb32 (1, 1) = cmplx_t (-51, 110);
  cb32 (2, 1) = cmplx_t (75, 63);
#else
  cm_t cb32 (2, 3); 
  cb32 (0, 0) = cmplx_t (60, -55);
  cb32 (0, 1) = cmplx_t (34, 58);
  cb32 (0, 2) = cmplx_t (13, -152);
  cb32 (1, 0) = cmplx_t (70, 10);
  cb32 (1, 1) = cmplx_t (-51, 110);
  cb32 (1, 2) = cmplx_t (75, 63);
#endif 
  print_m (cb32, "cb"); 
  cout << endl; 
  
  ierr = atlas::posv (ha, cb32); 

  if (ierr == 0)
    print_m (cb32, "cx"); 
  else 
    cout << "matrix is not positive definite: ierr = " 
         << ierr << endl << endl; 
  cout << endl; 

}

