
// solving A * X = B
// A symmetric/hermitian positive definite in packed format 
// factor (potrf()) and solve (potrs())

#include <cstddef>
#include <iostream>
#include <complex>
#include <boost/numeric/bindings/lapack/ppsv.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_symmetric.hpp>
#include <boost/numeric/bindings/traits/ublas_hermitian.hpp>
#include "utils.h"

namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;

using std::size_t; 
using std::cout;
using std::endl; 

typedef float real_t; 
typedef std::complex<real_t> cmplx_t; 

typedef ublas::matrix<real_t, ublas::column_major> m_t;
typedef ublas::matrix<cmplx_t, ublas::column_major> cm_t;

typedef 
  ublas::symmetric_matrix<real_t, ublas::lower, ublas::column_major> symml_t; 
typedef 
  ublas::hermitian_matrix<cmplx_t, ublas::lower, ublas::column_major> herml_t; 

typedef 
  ublas::symmetric_matrix<real_t, ublas::upper, ublas::column_major> symmu_t; 
typedef 
  ublas::hermitian_matrix<cmplx_t, ublas::upper, ublas::column_major> hermu_t; 

int main() {

  // for more descriptive comments see ublas_ppsv.cc 
  cout << endl; 

  // symmetric 
  cout << "real symmetric\n" << endl; 

  size_t n = 5; 
  size_t nrhs = 2; 
  symml_t sal (5, 5);   // symmetric matrix
  symmu_t sau (5, 5);   // symmetric matrix
  m_t x (n, nrhs);
  m_t bl (n, nrhs), bu (n, nrhs);  // RHS matrices

  init_symm (sal, 'l'); 
  init_symm (sau, 'u'); 

  print_m (sal, "sal"); 
  cout << endl; 
  print_m (sau, "sau"); 
  cout << endl; 

  for (int i = 0; i < x.size1(); ++i) {
    x (i, 0) = 1.;
    x (i, 1) = 2.; 
  }
  bl = prod (sal, x); 
  bu = prod (sau, x); 

  print_m (bl, "bl"); 
  cout << endl; 
  print_m (bu, "bu"); 
  cout << endl; 

  int ierr = lapack::pptrf (sal);  
  if (!ierr) {
    symml_t isal (sal);
    lapack::pptrs (sal, bl); 
    print_m (bl, "xl"); 
    cout << endl; 
    lapack::pptri (isal);
    print_m (isal, "isal");
  }
  cout << endl; 

  ierr = lapack::pptrf (sau);  
  if (!ierr) {
    symmu_t isau (sau);
    lapack::pptrs (sau, bu); 
    print_m (bu, "xu"); 
    cout << endl; 
    lapack::pptri (isau);
    print_m (isau, "isau");
  }
  cout << endl; 

  //////////////////////////////////////////////////////////
  // hermitian 
  cout << "\n==========================================\n" << endl; 
  cout << "complex hermitian\n" << endl; 

  herml_t hal (3, 3);   // hermitian matrix 
  hermu_t hau (3, 3);   // hermitian matrix
  cm_t cx (3, 1); 
  cm_t cbl (3, 1), cbu (3, 1);  // RHS

  hal (0, 0) = cmplx_t (25, 0);
  hal (1, 0) = cmplx_t (-5, 5);
  hal (1, 1) = cmplx_t (51, 0);
  hal (2, 0) = cmplx_t (10, -5);
  hal (2, 1) = cmplx_t (4, 6);
  hal (2, 2) = cmplx_t (71, 0);

  hau (0, 0) = cmplx_t (25, 0);
  hau (0, 1) = cmplx_t (-5, -5);
  hau (0, 2) = cmplx_t (10, 5);
  hau (1, 1) = cmplx_t (51, 0);
  hau (1, 2) = cmplx_t (4, -6);
  hau (2, 2) = cmplx_t (71, 0);

  print_m (hal, "hal"); 
  cout << endl; 
  print_m (hau, "hau"); 
  cout << endl; 

  cm_t cbl2 (3, 2); 
  cbl2 (0, 0) = cmplx_t (60, -55);
  cbl2 (1, 0) = cmplx_t (34, 58);
  cbl2 (2, 0) = cmplx_t (13, -152);
  cbl2 (0, 1) = cmplx_t (70, 10);
  cbl2 (1, 1) = cmplx_t (-51, 110);
  cbl2 (2, 1) = cmplx_t (75, 63);
  cm_t cbu2 (cbl2); 
  print_m (cbl2, "cbl"); 
  cout << endl; 
  
  ierr = lapack::pptrf (hal); 
  if (ierr == 0) {
    herml_t ihal (hal);
    lapack::pptrs (hal, cbl2); 
    print_m (cbl2, "cxl"); 
    cout << endl; 
    lapack::pptri (ihal);
    print_m (ihal, "ihal");
  }
  else 
    cout << "matrix is not positive definite: ierr = " 
         << ierr << endl << endl; 
  cout << endl; 

  ierr = lapack::pptrf (hau); 
  if (ierr == 0) {
    hermu_t ihau (hau);
    ierr = lapack::pptrs (hau, cbu2); 
    print_m (cbu2, "cxu"); 
    cout << endl; 
    lapack::pptri (ihau);
    print_m (ihau, "ihau");
  }
  else 
    cout << "matrix is not positive definite: ierr = " 
         << ierr << endl << endl; 
  cout << endl; 

}

