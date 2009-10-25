
// solving A * X = B
// A symmetric/hermitian positive definite
// driver function posv()

#include <cstddef>
#include <iostream>
#include <complex>
#include <boost/numeric/bindings/lapack/posv.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_symmetric.hpp>
#include <boost/numeric/bindings/traits/ublas_hermitian.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
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
  //        [5 0 0 0 0]
  //        [4 5 0 0 0]
  //   al = [3 4 5 0 0]
  //        [2 3 4 5 0]
  //        [1 2 3 4 5]

  init_symm (sau, 'u'); 
  //        [5 4 3 2 1]
  //        [0 5 4 3 2]
  //   au = [0 0 5 4 3]
  //        [0 0 0 5 4]
  //        [0 0 0 0 5]

  print_m (sal, "sal"); 
  cout << endl; 
  print_m (al, "al"); 
  cout << endl; 
  print_m_data (sal, "sal"); 
  cout << endl; 

  print_m (sau, "sau"); 
  cout << endl; 
  print_m (au, "au"); 
  cout << endl; 
  print_m_data (sau, "sau"); 
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

  m_t al1 (al), au1 (au);  // for part 2
  m_t bl1 (bl), bu1 (bu); 

  lapack::posv (sal, bl);  
  print_m (bl, "xl"); 
  cout << endl; 

  lapack::posv (sau, bu);  
  print_m (bu, "xu"); 
  cout << endl; 

  // part 2 

  lapack::posv ('L', al1, bl1);  
  print_m (al1, "al1 factored"); 
  cout << endl; 
  print_m (bl1, "xl1"); 
  cout << endl; 

  lapack::posv ('U', au1, bu1);  
  print_m (au1, "au1 factored"); 
  cout << endl; 
  print_m (bu1, "xu"); 
  cout << endl; 


  //////////////////////////////////////////////////////////
  // hermitian 
  cout << "\n==========================================\n" << endl; 
  cout << "complex hermitian (well, not really ;o)\n" << endl; 

  cm_t cal (3, 3), cau (3, 3);   // matrices (storage)
  herml_t hal (cal);   // hermitian adaptor 
  hermu_t hau (cau);   // hermitian adaptor 
  cm_t cx (3, 1); 
  cm_t cbl (3, 1), cbu (3, 1);  // RHS

  init_symm (hal, 'l'); 
  init_symm (hau, 'u'); 

  print_m (hal, "hal"); 
  cout << endl; 
  print_m (cal, "cal"); 
  cout << endl; 
  print_m (hau, "hau"); 
  cout << endl; 
  print_m (cau, "cau"); 
  cout << endl; 

  for (int i = 0; i < cx.size1(); ++i) 
    cx (i, 0) = cmplx_t (1, -1); 
  print_m (cx, "cx"); 
  cout << endl; 
  cbl = prod (hal, cx);
  cbu = prod (hau, cx);
  print_m (cbl, "cbl"); 
  cout << endl; 
  print_m (cbu, "cbu"); 
  cout << endl; 

  int ierr = lapack::posv (hal, cbl); 
  if (ierr == 0)
    print_m (cbl, "cxl"); 
  else 
    cout << "matrix is not positive definite: ierr = " 
         << ierr << endl;
  cout << endl; 

  ierr = lapack::posv (hau, cbu); 
  if (ierr == 0)
    print_m (cbu, "cxu"); 
  else 
    cout << "matrix is not positive definite: ierr = " 
         << ierr << endl;
  cout << endl; 

  cout << "\n===========================\n" << endl; 
  cout << "complex hermitian\n" << endl; 

  // regular (see ublas_gesv.cc), but not positive definite: 

  hal (0, 0) = cmplx_t (3, 0);
  hal (1, 0) = cmplx_t (4, -2);
  hal (1, 1) = cmplx_t (5, 0);
  hal (2, 0) = cmplx_t (-7, -5);
  hal (2, 1) = cmplx_t (0, 3);
  hal (2, 2) = cmplx_t (2, 0);

  hau (0, 0) = cmplx_t (3, 0);
  hau (0, 1) = cmplx_t (4, 2);
  hau (0, 2) = cmplx_t (-7, 5);
  hau (1, 1) = cmplx_t (5, 0);
  hau (1, 2) = cmplx_t (0, -3);
  hau (2, 2) = cmplx_t (2, 0);

  print_m (cal, "cal"); 
  cout << endl; 
  print_m (cau, "cau"); 
  cout << endl; 

  for (int i = 0; i < cx.size1(); ++i) 
    cx (i, 0) = cmplx_t (1, 1); 
  print_m (cx, "cx"); 
  cout << endl; 
  cbl = prod (hal, cx);
  cbu = prod (hau, cx);
  print_m (cbl, "cbl"); 
  cout << endl; 
  print_m (cbu, "cbu"); 
  cout << endl; 

  ierr = lapack::posv (hal, cbl); 
  if (ierr == 0)
    print_m (cbl, "cxl"); 
  else 
    cout << "matrix is not positive definite: ierr = " 
         << ierr << endl;
  cout << endl; 

  ierr = lapack::posv (hau, cbu); 
  if (ierr == 0)
    print_m (cbu, "cxu"); 
  else 
    cout << "matrix is not positive definite: ierr = " 
         << ierr << endl;
  cout << endl; 

  cout << "\n===========================\n" << endl; 
  cout << "complex hermitian\n" << endl; 

  // positive definite: 

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
  
  ierr = lapack::posv ('L', cal, cbl2); 
  if (ierr == 0)
    print_m (cbl2, "cxl"); 
  else 
    cout << "matrix is not positive definite: ierr = " 
         << ierr << endl << endl; 
  cout << endl; 

  ierr = lapack::posv ('U', cau, cbu2); 
  if (ierr == 0)
    print_m (cbu2, "cxu"); 
  else 
    cout << "matrix is not positive definite: ierr = " 
         << ierr << endl << endl; 
  cout << endl; 


}

