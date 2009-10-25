
// solving A * X = B
// A hermitian in packed storage 
// hetrf() & hetrs()

#include <cstddef>
#include <iostream>
#include <complex>
#include <boost/numeric/bindings/lapack/hpsv.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_hermitian.hpp>
#include <boost/numeric/bindings/traits/std_vector.hpp>
#include "utils.h"

namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;

using std::size_t; 
using std::cin;
using std::cout;
using std::endl; 

typedef double real_t; 
typedef std::complex<real_t> cmplx_t; 

typedef ublas::matrix<cmplx_t, ublas::column_major> cm_t;
typedef 
  ublas::hermitian_matrix<cmplx_t, ublas::lower, ublas::column_major> cherml_t;
typedef 
  ublas::hermitian_matrix<cmplx_t, ublas::upper, ublas::column_major> chermu_t;

int main() {

  cherml_t hcal (3, 3);   // hermitian matrix 
  chermu_t hcau (3, 3);   // hermitian matrix 
  cm_t cx (3, 1);
  cm_t cbl (3, 1), cbu (3, 1);  // RHS

  hcal (0, 0) = cmplx_t (3, 0);
  hcal (1, 0) = cmplx_t (4, -2);
  hcal (1, 1) = cmplx_t (5, 0);
  hcal (2, 0) = cmplx_t (-7, -5);
  hcal (2, 1) = cmplx_t (0, 3);
  hcal (2, 2) = cmplx_t (2, 0);

  hcau (0, 0) = cmplx_t (3, 0);
  hcau (0, 1) = cmplx_t (4, 2);
  hcau (0, 2) = cmplx_t (-7, 5);
  hcau (1, 1) = cmplx_t (5, 0);
  hcau (1, 2) = cmplx_t (0, -3);
  hcau (2, 2) = cmplx_t (2, 0);

  print_m (hcal, "hcal"); 
  cout << endl; 
  print_m (hcau, "hcau"); 
  cout << endl; 

  for (int i = 0; i < cx.size1(); ++i) 
    cx (i, 0) = cmplx_t (1, -1); 
  print_m (cx, "cx"); 
  cout << endl; 
  cbl = prod (hcal, cx);
  cbu = prod (hcau, cx);
  print_m (cbl, "cbl"); 
  cout << endl; 
  print_m (cbu, "cbu"); 
  cout << endl; 

  std::vector<int> ipiv (3); 

  int ierr = lapack::hptrf (hcal, ipiv); 
  if (ierr == 0) {
    lapack::hptrs (hcal, ipiv, cbl); 
    print_v (ipiv, "ipiv"); 
    cout << endl; 
    print_m (cbl, "cxl"); 
  }
  else 
    cout << "matrix is not regular: ierr = " 
         << ierr << endl;
  cout << endl; 

  ierr = lapack::hptrf (hcau, ipiv); 
  if (ierr == 0) {
    lapack::hptrs (hcau, ipiv, cbu); 
    print_v (ipiv, "ipiv"); 
    cout << endl; 
    print_m (cbu, "cxu"); 
  }
  else 
    cout << "matrix is not regular: ierr = " 
         << ierr << endl;
  cout << endl; 
}

