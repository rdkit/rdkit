
// BLAS level 2

//#define BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 
//#define BOOST_NO_FUNCTION_TEMPLATE_ORDERING

#include <iostream>
#include <boost/numeric/bindings/atlas/cblas1.hpp>
#include <boost/numeric/bindings/atlas/cblas2.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include "utils.h" 

namespace ublas = boost::numeric::ublas;
namespace atlas = boost::numeric::bindings::atlas;

using std::cout;
using std::endl; 

typedef ublas::vector<double> vct_t;
typedef ublas::matrix<double, ublas::row_major> rm_t;
typedef ublas::matrix<double, ublas::column_major> cm_t;

int main() {

  cout << endl; 

  vct_t vx (2);
  vct_t vy (4); 

  // row major matrix
  rm_t rm (2, 4);
  init_m (rm, const_val<double> (0)); 
  print_m (rm, "row major matrix m"); 
  cout << endl; 

  vx(0) = 1.; 
  vy(1) = 1.; 
  print_v (vx, "vx"); 
  cout << endl; 
  print_v (vy, "vy"); 
  cout << endl; 

  // m += x y^T
  atlas::ger (vx, vy, rm); 
  print_m (rm, "m += x y^T"); 
  cout << endl << endl; 

  init_m (rm, const_val<double> (1)); 
  print_m (rm, "m"); 
  cout << endl; 

  atlas::set (1., vx);
  atlas::set (1., vy);
  print_v (vx, "vx"); 
  cout << endl; 
  print_v (vy, "vy"); 
  cout << endl; 

  // m += 2 x y^T
  atlas::ger (2., vx, vy, rm); 
  print_m (rm, "m += 2 x y^T"); 
  cout << endl << endl; 

  init_v (vx, iplus1());
  init_v (vy, iplus1());
  print_v (vx, "vx"); 
  cout << endl; 
  print_v (vy, "vy"); 
  cout << endl; 

  // m += x y^T
  atlas::ger (vx, vy, rm); 
  print_m (rm, "m += x y^T"); 
  cout << endl << endl; 

  // column major matrix
  cm_t cm (2, 4);
  init_m (cm, const_val<double> (0)); 
  print_m (cm, "column major matrix m"); 
  cout << endl; 

  vx(0) = 1.; 
  vy(1) = 1.; 
  print_v (vx, "vx"); 
  cout << endl; 
  print_v (vy, "vy"); 
  cout << endl; 

  // m += x y^T
  atlas::ger (vx, vy, cm); 
  print_m (cm, "m += x y^T"); 
  cout << endl << endl; 

  init_m (cm, const_val<double> (1)); 
  print_m (cm, "m"); 
  cout << endl; 

  atlas::set (1., vx);
  atlas::set (1., vy);
  print_v (vx, "vx"); 
  cout << endl; 
  print_v (vy, "vy"); 
  cout << endl; 

  // m += 2 x y^T
  atlas::ger (2., vx, vy, cm); 
  print_m (cm, "m += 2 x y^T"); 
  cout << endl << endl; 

  init_v (vx, iplus1());
  init_v (vy, iplus1());
  print_v (vx, "vx"); 
  cout << endl; 
  print_v (vy, "vy"); 
  cout << endl; 

  // m += x y^T
  atlas::ger (vx, vy, cm); 
  print_m (cm, "m += x y^T"); 
  cout << endl << endl; 

}
