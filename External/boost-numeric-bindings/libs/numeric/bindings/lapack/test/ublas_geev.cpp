#include<iostream>
#include<cstdlib>
#include<complex>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/lapack/geev.hpp>
#include "utils.h"

using std::cout;
using std::endl;
using std::vector;
using std::complex;

namespace ublas =  boost::numeric::ublas;
namespace lapack =  boost::numeric::bindings::lapack;

void geev(int);
template <typename T>
void Hessenberg(ublas::matrix<T, ublas::column_major>& );

int main(){
    cout << "I'm testing uBlas." << endl;

    int n = 5;
    geev(n);

}
void geev(int n){
    cout << "\nCalculating eigenvalues using LAPACK's geev." << endl;
    ublas::matrix<double, ublas::column_major> A(n,n);
    Hessenberg(A);
    print_m(A);

    ublas::vector<complex<double> > values(n);
    ublas::matrix<complex<double>, ublas::column_major>* Vectors_left = 0;
    ublas::matrix<complex<double>, ublas::column_major> Vectors_right(n,n);

    lapack::geev(A, values, Vectors_left, &Vectors_right, lapack::optimal_workspace());
    print_v(values, "values"); cout << endl;
    print_m(Vectors_right, "Vectors_right"); cout << endl;

    Hessenberg(A);
    cout << "A*x = l*x." << endl;
    for( int i = 0; i < Vectors_right.size2(); ++i ){
        ublas::vector<complex<double> > tmp(n);
        tmp = ublas::prod( A, column(Vectors_right, i) );
        cout << tmp - values(i)*column(Vectors_right,i) << endl;
    }

}
template <typename T>
void Hessenberg(ublas::matrix<T, ublas::column_major>& H){
    T k = 1;
    for( unsigned int i = 0; i < H.size1(); ++i ){
        for( unsigned int j = i; j <= H.size2(); ++j ){
            if( j > 0 ){
                H(i,j-1) = k;
                k += 1;
            }
        }
    }
}

