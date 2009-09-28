
// solving A * X = B
// A symmetric in packed format 
// sytrf() & sytrs() 

#include <cstddef>
#include <iostream>
#include <complex>
#include <boost/numeric/bindings/lapack/spsv.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_symmetric.hpp>
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

typedef ublas::matrix<real_t, ublas::column_major> m_t;
typedef ublas::matrix<cmplx_t, ublas::column_major> cm_t;

typedef 
  ublas::symmetric_matrix<real_t, ublas::lower, ublas::column_major> symml_t; 
typedef 
  ublas::symmetric_matrix<cmplx_t, ublas::lower, ublas::column_major> csymml_t;

typedef 
  ublas::symmetric_matrix<real_t, ublas::upper, ublas::column_major> symmu_t; 
typedef 
  ublas::symmetric_matrix<cmplx_t, ublas::upper, ublas::column_major> csymmu_t;


template <typename M>
void init_symm2 (M& m) {
  for (int i = 0; i < m.size1(); ++i) 
    for (int j = i; j < m.size1(); ++j)
      m (i, j) = m (j, i) = 1 + j - i; 
}

int main (int argc, char **argv) {
  size_t n = 0;
  if (argc > 1) {
    n = atoi(argv [1]);
  }

  cout << endl; 

  cout << "real symmetric\n" << endl; 

  if (n <= 0) {
  cout << "n -> ";
  cin >> n;
  }
  if (n < 5) n = 5; 
  cout << "min n = 5" << endl << endl; 
  size_t nrhs = 2; 
  symml_t sal (n, n);   // symmetric matrix
  symmu_t sau (n, n);   // symmetric matrix
  m_t x (n, nrhs);
  m_t bl (n, nrhs), bu (n, nrhs);  // RHS matrices

  init_symm2 (sal); 
  print_m (sal, "sal"); 
  cout << endl; 

  init_symm2 (sau); 
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

  std::vector<int> ipiv (n); 
  
  int err = lapack::sptrf (sal, ipiv);  
  if (err == 0) {
    symml_t isal (sal);
    lapack::sptrs (sal, ipiv, bl); 
    print_m (bl, "xl"); 
    lapack::sptri (isal, ipiv);
    print_m (isal, "isal"); 
  } 
  cout << endl; 

  err = lapack::sptrf (sau, ipiv);  
  if (err == 0) {
    symmu_t isau (sau);
    lapack::sptrs (sau, ipiv, bu); 
    print_m (bu, "xu"); 
    lapack::sptri (isau, ipiv);
    print_m (isau, "isau"); 
  } 
  else 
    cout << "?" << endl; 
  cout << endl; 

  //////////////////////////////////////////////////////////
  cout << "\n==========================================\n" << endl; 
  cout << "complex symmetric\n" << endl; 

  csymml_t scal (n, n);   // symmetric matrix 
  csymmu_t scau (n, n);   // symmetric matrix
  cm_t cx (n, 1); 
  cm_t cbl (n, 1), cbu (n, 1);  // RHS

  init_symm2 (scal); 
  init_symm2 (scau); 
  scal *= cmplx_t (0.1, 0.25); 
  scau *= cmplx_t (-1, -0.5); 

  print_m (scal, "scal"); 
  cout << endl; 
  print_m (scau, "scau"); 
  cout << endl; 

  for (int i = 0; i < cx.size1(); ++i) 
    cx (i, 0) = cmplx_t (1, -1); 
  print_m (cx, "cx"); 
  cout << endl; 
  cbl = prod (scal, cx);
  cbu = prod (scau, cx);
  print_m (cbl, "cbl"); 
  cout << endl; 
  print_m (cbu, "cbu"); 
  cout << endl; 

  int ierr = lapack::sptrf (scal, ipiv); 
  if (ierr == 0) {
    csymml_t iscal (scal);
    lapack::sptrs (scal, ipiv, cbl); 
    print_m (cbl, "cxl"); 
    lapack::sptri (iscal, ipiv);
    print_m (iscal, "iscal"); 
  }
  else 
    cout << "?" << endl;
  cout << endl; 

  ierr = lapack::sptrf (scau, ipiv); 
  if (ierr == 0) {
    csymmu_t iscau (scau);
    lapack::sptrs (scau, ipiv, cbu); 
    print_m (cbu, "cxu"); 
    lapack::sptri (iscau, ipiv);
    print_m (iscau, "iscau"); 
  }
  else 
    cout << "?" << endl;
  cout << endl; 

}

