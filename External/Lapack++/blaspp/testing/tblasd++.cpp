//
//              LAPACK++ 1.1 Linear Algebra Package 1.1
//               University of Tennessee, Knoxvilee, TN.
//            Oak Ridge National Laboratory, Oak Ridge, TN.
//        Authors: J. J. Dongarra, E. Greaser, R. Pozo, D. Walker
//                 (C) 1992-1996 All Rights Reserved
//
//                             NOTICE
//
// Permission to use, copy, modify, and distribute this software and
// its documentation for any purpose and without fee is hereby granted
// provided that the above copyright notice appear in all copies and
// that both the copyright notice and this permission notice appear in
// supporting documentation.
//
// Neither the Institutions (University of Tennessee, and Oak Ridge National
// Laboratory) nor the Authors make any representations about the suitability 
// of this software for any purpose.  This software is provided ``as is'' 
// without express or implied warranty.
//
// LAPACK++ was funded in part by the U.S. Department of Energy, the
// National Science Foundation and the State of Tennessee.




#include <iostream>
#include "lafnames.h"
//#include LA_VECTOR_DOUBLE_H // changed of VC++
#include "lavd.h"

#include "blas++.h"

int main(int argc, char *argv[])
{
    int M = atoi(argv[1]);
    int N = atoi(argv[2]);

    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " M N \n";
        exit(1);
    }

//  test Blas_Sum()

    LaVectorDouble Sum(M);
    Sum = N;
    double ans = Blas_Norm1(Sum);
    fprintf(stdout,"\nBlas_Sum() A:%dx1, value:%d\n",M,N);
    std::cout << "output:\n" << ans << std::endl;

//  test Blas_Add_Mult()

    LaVectorDouble X(M);
    LaVectorDouble Y(M);
    X = N;
    Y = N;
    double scalar = M;
    Blas_Add_Mult(Y, scalar,X);
    fprintf(stdout,"\nBlas_Add_Mult() alpha:%d, X:%dx1, Y:%dx1\n",\
                    M,N,N); 
    std::cout << "output:\n" << Y << std::endl;

//  test Blas_Copy()

    X = M*N;
    Blas_Copy(Y, X);
    std::cout <<"\nBlas_Copy():\n" << "X:\n" << X << "\nY:\n" << Y << std::endl; 

//  test Blas_Dot_Prod()

    X = M;
    Y = N;
    double cans = Blas_Dot_Prod(X,Y);
    fprintf(stdout,"\nBlas_Dot_Prod() X = %d, Y = %d\n",M,N); 
    fprintf(stdout,"  X is %dx1, Y is %dx1\n",M,M); 
    std::cout << "\nAns:\n" << cans << std::endl;

//  test Blas_Norm2()

    ans = Blas_Norm2(X);
    fprintf(stdout,"\nBlas_Norm2() X = %d\n",M); 
    fprintf(stdout,"  X is %dx1\n",M); 
    std::cout << "\nAns:\n" << ans << std::endl;

// see note in blas1++.cc
//#if 0
//  test Blas_Scale()

    double scale = 5.0;
    X = 1.1;
    fprintf(stdout,"\nBlas_Scale() scale = 5.0, X = 1.1\n");
    Blas_Scale(scale,X);
    std::cout <<"X:\n"<< X << std::endl;
// #endif

//  test Blas_Swap()

    LaVectorDouble A(5);
    LaVectorDouble B(5);
    A = 1.1;
    B = 2.0;
    fprintf(stdout,"\nBlas_Swap() A = 1.1, B = 2.0\n");
    Blas_Swap(A,B);
    std::cout <<"A:\n"<< A << "\nB:\n" << B << std::endl;


//  test Blas_Index_Max()


    int index;
    X = 8.0;
    X (M/2) = 64.0;
    fprintf(stdout,"\nBlas_Index_Max() X = 8.0\n");
    index = Blas_Index_Max(X);
    std::cout <<"index:\n"<< index << std::endl;

    return 0;
}
