
// solving A * X = B
// A symmetric
// driver function sysv()

#define BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 
#define BOOST_NO_FUNCTION_TEMPLATE_ORDERING

#include <cstddef>
#include <iostream>
#include <complex>
#include <boost/numeric/bindings/lapack/sysv.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_symmetric.hpp>
#include <boost/numeric/bindings/traits/ublas_hermitian.hpp>
#include "utils2.h"

namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;

using std::size_t; 
using std::cin;
using std::cout;
using std::endl; 

typedef double real; 
typedef std::complex<real> cmplx_t; 

typedef ublas::matrix<real, ublas::column_major> m_t;
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

int main() {

  cout << endl; 

  // symmetric 
  cout << "real symmetric\n" << endl; 

  size_t n;
  cout << "n -> ";
  cin >> n;
  if (n < 5) n = 5; 
  cout << "min n = 5" << endl << endl; 
  size_t nrhs = 2; 
  m_t al (n, n), au (n, n);  // matrices (storage)
  symml_t sal (al);   // symmetric adaptor
  symmu_t sau (au);   // symmetric adaptor
  m_t x (n, nrhs);
  m_t bl (n, nrhs), bu (n, nrhs);  // RHS matrices

  init_symm2 (al); 
  ublas::swap (row (al, 1), row (al, 4)); 
  ublas::swap (column (al, 1), column (al, 4)); 

  print_m (al, "al"); 
  cout << endl; 

  init_symm2 (au); 
  ublas::swap (row (au, 2), row (au, 3)); 
  ublas::swap (column (au, 2), column (au, 3)); 

  print_m (au, "au"); 
  cout << endl; 

  for (int i = 0; i < x.size1(); ++i) {
    x (i, 0) = 1.;
    x (i, 1) = 2.; 
  }
  bl = ublas::prod (sal, x); 
  bu = ublas::prod (sau, x); 

  print_m (bl, "bl"); 
  cout << endl; 
  print_m (bu, "bu"); 
  cout << endl; 

  m_t al1 (al), au1 (au);  // for part 2
  m_t bl1 (bl), bu1 (bu); 

  lapack::sysv (sal, bl);  
  print_m (bl, "xl"); 
  cout << endl; 

  lapack::sysv (sau, bu);  
  print_m (bu, "xu"); 
  cout << endl; 

  // part 2 

  std::vector<int> ipiv (n); 
  std::vector<real> work (1); 

  int err = lapack::sysv ('L', al1, ipiv, bl1, work);  
  print_m (al1, "al1 factored"); 
  cout << endl; 
  print_v (ipiv, "ipiv"); 
  cout << endl; 
  print_m (bl1, "xl1"); 
  cout << endl; 

  err = lapack::sysv ('U', au1, ipiv, bu1, work);  
  print_m (au1, "au1 factored"); 
  cout << endl; 
  print_v (ipiv, "ipiv"); 
  cout << endl; 
  print_m (bu1, "xu1"); 
  cout << endl; 
  cout << endl; 

  //////////////////////////////////////////////////////////
  cout << "\n==========================================\n" << endl; 
  cout << "complex symmetric\n" << endl; 

  cm_t cal (n, n), cau (n, n);   // matrices (storage)
  csymml_t scal (cal);   // hermitian adaptor 
  csymmu_t scau (cau);   // hermitian adaptor 
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
  cbl = ublas::prod (scal, cx);
  cbu = ublas::prod (scau, cx);
  print_m (cbl, "cbl"); 
  cout << endl; 
  print_m (cbu, "cbu"); 
  cout << endl; 

  int ierr = lapack::sysv (scal, cbl); 
  if (ierr == 0)
    print_m (cbl, "cxl"); 
  else 
    cout << "matrix is not regular: ierr = " 
         << ierr << endl;
  cout << endl; 

  std::vector<cmplx_t> cwork (n); 

  ierr = lapack::sysv (scau, ipiv, cbu, cwork); 
  if (ierr == 0) {
    print_v (ipiv, "ipiv"); 
    cout << endl; 
    print_m (cbu, "cxu"); 
  }
  else 
    cout << "matrix is not regular: ierr = " 
         << ierr << endl;
  cout << endl; 

}

