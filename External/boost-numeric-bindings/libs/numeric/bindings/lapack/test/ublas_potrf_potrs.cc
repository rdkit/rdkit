
// solving A * X = B
// A symmetric/hermitian positive definite
// factor (potrf()) and solve (potrs())

#include <cstddef>
#include <iostream>
#include <complex>
#include <boost/numeric/bindings/lapack/posv.hpp>
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

typedef ublas::symmetric_adaptor<m_t, ublas::lower> symml_t; 
typedef ublas::hermitian_adaptor<cm_t, ublas::lower> herml_t; 

typedef ublas::symmetric_adaptor<m_t, ublas::upper> symmu_t; 
typedef ublas::hermitian_adaptor<cm_t, ublas::upper> hermu_t; 

int main() {

  // for more descriptive comments see ublas_posv.cc 
  cout << endl; 

  // symmetric 
  cout << "real symmetric\n" << endl; 

  size_t n = 5; 
  size_t nrhs = 2; 
  m_t al (n, n), au (n, n);  // matrices (storage)
  symml_t sal (al);   // symmetric adaptor
  symmu_t sau (au);   // symmetric adaptor
  m_t x (n, nrhs);
  m_t bl (n, nrhs), bu (n, nrhs);  // RHS matrices

  init_symm (sal, 'l'); 
  init_symm (sau, 'u'); 

  print_m (sal, "sal"); 
  cout << endl; 
  print_m (al, "al"); 
  cout << endl; 

  print_m (sau, "sau"); 
  cout << endl; 
  print_m (au, "au"); 
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

  int ierr = lapack::potrf (sal);  
  if (!ierr) {
    lapack::potrs (sal, bl); 
    print_m (bl, "xl"); 
  }
  cout << endl; 

  ierr = lapack::potrf ('U', au);  
  if (!ierr) {
    lapack::potrs ('U', au, bu); 
    print_m (bu, "xu"); 
  }
  cout << endl; 

  //////////////////////////////////////////////////////////
  // hermitian 
  cout << "\n==========================================\n" << endl; 
  cout << "complex hermitian\n" << endl; 

  cm_t cal (3, 3), cau (3, 3);   // matrices (storage)
  herml_t hal (cal);   // hermitian adaptor 
  hermu_t hau (cau);   // hermitian adaptor 
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

  print_m (cal, "cal"); 
  cout << endl; 
  print_m (cau, "cau"); 
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
  
  ierr = lapack::potrf ('L', cal); 
  if (ierr == 0) {
    lapack::potrs ('L', cal, cbl2); 
    print_m (cbl2, "cxl"); 
  }
  else 
    cout << "matrix is not positive definite: ierr = " 
         << ierr << endl << endl; 
  cout << endl; 

  ierr = lapack::potrf (hau); 
  if (ierr == 0) {
    ierr = lapack::potrs (hau, cbu2); 
    print_m (cbu2, "cxu"); 
  }
  else 
    cout << "matrix is not positive definite: ierr = " 
         << ierr << endl << endl; 
  cout << endl; 


}

