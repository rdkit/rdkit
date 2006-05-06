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

#if 0
#include "lapack++.h"
#endif

#if 1
#include "lafnames.h"       /* macros for LAPACK++ filenames */
#include "lapack.h"
//#include LA_GEN_MAT_DOUBLE_H // changed of VC++
#include "gmd.h"
//#include LA_VECTOR_DOUBLE_H // changed of VC++
#include "lavd.h"
#include "blas++.h"
//#include LA_SOLVE_DOUBLE_H // changed of VC++
#include "laslv.h"
//#include LA_GENERATE_MAT_DOUBLE_H // changed of VC++
#include "genmd.h"
//#include LA_EXCEPTION_H // changed of VC++
#include "laexcp.h"
//#include LA_UTIL_H // changed of VC++
#include "lautil.h"
#endif

double residual(LaGenMatDouble &A, LaVectorDouble &x, 
    const LaVectorDouble& b)
{
    int M = A.size(0);
    int N = A.size(1);


    std::cout << "\tNorm_Inf(A*x-b)" << Norm_Inf(A*x-b) << std::endl;
    std::cout << "\tNorm_Inf(A) " << Norm_Inf(A) << std::endl;
    std::cout << "\tNorm_Inf(x) " << Norm_Inf(x) << std::endl;
    std::cout << "\tMacheps :" << Mach_eps_double() << std::endl;

    if (M>N)
    {
        LaVectorDouble Axb = A*x-b;
        LaVectorDouble R(M);

        Blas_Mat_Trans_Vec_Mult(A, Axb, R);
        return Norm_Inf(R) / 
            (Norm_Inf(A)* Norm_Inf(x) * N * Mach_eps_double());

    }
    else
    {
        return Norm_Inf(A*x-b ) /
                ( Norm_Inf(A)* Norm_Inf(x) * N * Mach_eps_double());
    }
}


int TestGenLinearSolve(int M,int N)
{
    LaGenMatDouble A(M,N);
    LaVectorDouble x(N), b(M);

#ifndef HPPA
    const char fname[] = "TestGenLinearSolve(LaGenMat, x, b) ";
#else
    char *fname = NULL;
#endif

    //char e = 'e';
    double norm;
    double res;

    LaGenerateMatDouble(A);

    // save a snapshot of what A looked like before the solution
    LaGenMatDouble old_A = A;



    b = 1.1;

    std::cerr << fname << ": testing LaLinearSolve(Gen,...) M= "<< M
        << " N = " << N << std::endl;

    LaLinearSolve(A, x, b);

    if ( (norm = Norm_Inf( old_A - A)) >  0.0)  // this is an exact test, not
                                         // necessary to worry about
                                         // round-off issues.  We
                                         // are testing to see A was
                                         // overwritten.
    {
        std::cerr << fname << ": overwrote 1st arg.\n";
        std::cerr << "       error norm: " << norm << std::endl;
        exit(1);
    }

    res = residual(A,x,b);
    if (res > 1)
    {
        std::cerr << fname << "resdiual " << res << " is to too high.\n";
        exit(1);
    }

    std::cerr << fname << ": LaLinearSolve() success.\n\n";


    // now try the in-place solver


    std::cerr << fname << ": testing LaLinearSolveIP(Gen,...) \n";
    LaLinearSolveIP(A, x, b);


    res = residual(old_A, x, b);

    if (res > 1)
    {
        std::cerr << fname << "resdiual " << res << " is to too high.\n";
        exit(1);
    }

    std::cerr << fname << ": LaLinearSolveIP() success.\n\n";
    return 0;
}


int main(int argc, char **argv)
{

    std::cout.precision(4);
    std::cout.setf(std::ios::scientific, std::ios::floatfield);

    if (argc < 2)
    {
        std::cerr << "Usage " << argv[0] << " N " << std::endl;
        exit(1);
    }
    int N = atoi(argv[1]);
    int M = N;

    std::cout << "Testing " << N << " x " << N << " system." << std::endl;

    TestGenLinearSolve(M,N);

	return 0;
}

