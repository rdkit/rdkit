
// solving A * X = B
// A symmetric
// sytrf() & sytrs() 

#include <cstddef>
#include <iostream>
#include <complex>
#include <boost/numeric/bindings/lapack/sysv.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_symmetric.hpp>
#include <boost/numeric/bindings/traits/std_vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
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

typedef ublas::symmetric_adaptor<m_t, ublas::lower> symml_t; 
typedef ublas::symmetric_adaptor<m_t, ublas::upper> symmu_t; 

typedef ublas::symmetric_adaptor<cm_t, ublas::lower> csymml_t; 
typedef ublas::symmetric_adaptor<cm_t, ublas::upper> csymmu_t; 

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

  // symmetric 
  cout << "real symmetric\n" << endl; 

  if (n <= 0) {
  cout << "n -> ";
  cin >> n;
  }
  if (n < 5) n = 5; 
  cout << "min n = 5" << endl << endl; 
  size_t nrhs = 2; 
  m_t al (n, n), au (n, n);  // matrices (storage)
  symml_t sal (al);   // symmetric adaptor
  symmu_t sau (au);   // symmetric adaptor
  m_t x (n, nrhs);
  m_t bl (n, nrhs), bu (n, nrhs);  // RHS matrices

  init_symm2 (al); 
  swap (row (al, 1), row (al, 4)); 
  swap (column (al, 1), column (al, 4)); 

  print_m (al, "al"); 
  cout << endl; 

  init_symm2 (au); 
  swap (row (au, 2), row (au, 3)); 
  swap (column (au, 2), column (au, 3)); 

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

  m_t al1 (al), au1 (au);  // for part 2
  m_t bl1 (bl), bu1 (bu); 

  std::vector<int> ipiv (n); 
  
  int err = lapack::sytrf (sal, ipiv);  
  if (err == 0) {
    symml_t isal (sal);
    lapack::sytrs (sal, ipiv, bl); 
    print_m (bl, "xl"); 
    lapack::sytri (isal, ipiv);
    print_m (isal, "isal"); 
  } 
  cout << endl; 

  err = lapack::sytrf (sau, ipiv);  
  if (err == 0) {
    symmu_t isau (sau);
    lapack::sytrs (sau, ipiv, bu); 
    print_m (bu, "xu"); 
    lapack::sytri (isau, ipiv);
    print_m (isau, "isau"); 
  } 
  else 
    cout << "?" << endl; 
  cout << endl; 

  // part 2 

  cout << endl << "part 2" << endl << endl; 
  
  int lw = lapack::sytrf_block ('O', 'L', al1);
  cout << "nb = " << lw << endl;
  lw *= n; 
  cout << "lw = " << lw << " == " 
       << lapack::sytrf_work ('O', 'L', al1) << endl;
  cout << "mw = " << lapack::sytrf_work ('M', 'L', al1) << endl;
  std::vector<real_t> work (lw); 
  int mb = lapack::sytrf_block ('M', 'L', al1); 
  cout << "mb = " << mb << endl << endl;

  err = lapack::sytrf ('L', al1, ipiv, work);  
  if (err == 0) {
    lapack::sytrs ('L', al1, ipiv, bl1); 
    print_m (al1, "al1 factored"); 
    cout << endl; 
    print_v (ipiv, "ipiv"); 
    cout << endl; 
    print_m (bl1, "xl1"); 
  }
  else 
    cout << "?" << endl; 
  cout << endl; 

  lw = lapack::sytrf_block ('O', 'U', au1); 
  cout << "nb = " << lw << endl;
  lw *= n; 
  cout << "lw = " << lw << " == " 
       << lapack::sytrf_work ('O', 'U', au1) << endl;
  cout << "mw = " << lapack::sytrf_work ('M', 'U', au1) << endl;
  if (lw != work.size())
    work.resize (lw); 
  mb = lapack::sytrf_block ('M', 'U', au1);
  cout << "mb = " << mb << endl << endl;

  err = lapack::sytrf ('U', au1, ipiv, work);  
  if (err == 0) {
    lapack::sytrs ('U', au1, ipiv, bu1); 
    print_m (au1, "au1 factored"); 
    cout << endl; 
    print_v (ipiv, "ipiv"); 
    cout << endl; 
    print_m (bu1, "xu1"); 
  }
  else 
    cout << "?" << endl; 
  cout << endl; 
  cout << endl; 

  //////////////////////////////////////////////////////////
  cout << "\n==========================================\n" << endl; 
  cout << "complex symmetric\n" << endl; 

  cm_t cal (n, n), cau (n, n);   // matrices (storage)
  csymml_t scal (cal);   // symmetric adaptor 
  csymmu_t scau (cau);   // symmetric adaptor 
  cm_t cx (n, 1); 
  cm_t cbl (n, 1), cbu (n, 1);  // RHS

  init_symm2 (cal); 
  init_symm2 (cau); 
  cal *= cmplx_t (1, 1); 
  cau *= cmplx_t (1, -0.5); 

  print_m (cal, "cal"); 
  cout << endl; 
  print_m (cau, "cau"); 
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

  int ierr = lapack::sytrf (scal, ipiv); 
  if (ierr == 0) {
    csymml_t iscal (scal);
    lapack::sytrs (scal, ipiv, cbl); 
    print_m (cbl, "cxl"); 
    lapack::sytri (iscal, ipiv);
    print_m (iscal, "iscal"); 
  }
  else 
    cout << "?" << endl;
  cout << endl; 

  lw = lapack::sytrf_block ('O', scau); 
  cout << "nb = " << lw << endl;
  lw *= n; 
  cout << "lw = " << lw << " == " 
       << lapack::sytrf_work ('O', scau) << endl;
  cout << "mw = " << lapack::sytrf_work ('M', scau) << endl;
  std::vector<cmplx_t> cwork (lw); 
  mb = lapack::sytrf_block ('M', scau); 
  cout << "mb = " << mb << endl << endl;

  ierr = lapack::sytrf (scau, ipiv, cwork); 
  if (ierr == 0) {
    csymmu_t iscau (scau);
    lapack::sytrs (scau, ipiv, cbu); 
    print_v (ipiv, "ipiv"); 
    cout << endl; 
    print_m (cbu, "cxu"); 
    lapack::sytri (iscau, ipiv);
    print_m (iscau, "iscau"); 
  }
  else 
    cout << "?" << endl;
  cout << endl; 

}

