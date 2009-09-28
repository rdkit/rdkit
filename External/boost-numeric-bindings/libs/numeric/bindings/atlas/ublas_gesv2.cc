
// solving A * X = B
// using driver function gesv()

//#define BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 
//#define BOOST_NO_FUNCTION_TEMPLATE_ORDERING

#include <cstddef>
#include <iostream>
#include <boost/numeric/bindings/atlas/cblas.hpp>
#include <boost/numeric/bindings/atlas/clapack.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "utils.h"

namespace ublas = boost::numeric::ublas;
namespace atlas = boost::numeric::bindings::atlas;

using std::size_t; 
using std::cout;
using std::endl; 

#ifndef F_ROW_MAJOR
typedef ublas::matrix<double, ublas::column_major> m_t;
#else
typedef ublas::matrix<double, ublas::row_major> m_t;
#endif

int main() {

  cout << endl; 

  size_t n = 5;   
  m_t a (n, n);   // system matrix 

  size_t nrhs = 2; 
  m_t x (n, nrhs), bb (n, nrhs);  
  // b -- right-hand side matrix, see below 

  init_symm (a); 
  //     [n   n-1 n-2  ... 1]
  //     [n-1 n   n-1  ... 2]
  // a = [n-2 n-1 n    ... 3]
  //     [        ...       ]
  //     [1   2   ...  n-1 n]

  m_t const aa (a); // copy of a, because a is `lost' after gesv()

  ublas::matrix_column<m_t> xc0 (x, 0), xc1 (x, 1); 
  atlas::set (1., xc0);
  atlas::set (2., xc1);
  atlas::gemm (a, x, bb);  // bb = a x, so we know the result ;o) 

  print_m (a, "A"); 
  cout << endl; 
  print_m (bb, "B"); 
  cout << endl; 

  // see leading comments for `gesv()' in clapack.hpp
#ifndef F_ROW_MAJOR
  m_t b (bb); 
#else 
  m_t b (ublas::trans (bb)); 
#endif 
  print_m (b, "B for gesv()"); 
  cout << endl; 

  atlas::gesv (a, b);  // solving the system, b contains x 

#ifndef F_ROW_MAJOR
  print_m (b, "X");
  cout << endl; 
  atlas::gemm (aa, b, x); 
#else
  print_m (b, "X^T"); 
  cout << endl; 
  atlas::gemm (CblasNoTrans, CblasTrans, 1.0, aa, b, 0.0, x); 
#endif 
  print_m (x, "B = A X"); 

  cout << endl; 

}

