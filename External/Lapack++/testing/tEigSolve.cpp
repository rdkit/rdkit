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


#include "lafnames.h"       /* macros for LAPACK++ filenames */
#include "lapack.h"
//#include LA_GEN_MAT_DOUBLE_H // changed of VC++
#include "gmd.h"
//#include LA_VECTOR_DOUBLE_H // changed of VC++
#include "lavd.h"
//#include LA_SYMM_MAT_DOUBLE_H // changed of VC++
#include "symd.h"
#include "blas++.h"
//#include LA_SOLVE_DOUBLE_H // changed of VC++
#include "laslv.h"
//#include LA_GENERATE_MAT_DOUBLE_H // changed of VC++
#include "genmd.h"
//#include LA_EXCEPTION_H // changed of VC++
#include "laexcp.h"
//#include LA_UTIL_H // changed of VC++
#include "lautil.h"


double eig_residual(const LaSymmMatDouble &A, double lambda, 
        const LaVectorDouble &x)
{
    int N = A.size(0);

    return Norm_Inf(A*x-lambda*x) / 
        (Norm_Inf(A) * Norm_Inf(x) * N * Mach_eps_double());


}

void TestGenEigSolve(int N)
{
    LaSymmMatDouble A(N,N);
    LaVectorDouble  v(N);

#ifndef HPPA
    const char fname[] = "TestGenEigSolve() ";
#else
    char *fname = NULL;
#endif

    //char e = 'e';

    LaGenerateMatDouble(A);

    LaGenMatDouble Eigenvectors(N,N);

    std::cerr << fname << 
            ": testing LaEigSolve(LaSymmMat, eig_value, eig_vectors) \n";

    LaEigSolve(A, v, Eigenvectors);


    
    for (int i=0; i<A.size(0); i++)
    {
        LaIndex I(0,N-1);
        double res = eig_residual(A, v(i), Eigenvectors(I,i));
        if (res > 1)
        {
            std::cerr << fname << " residual " <<  res << " is too high.\n";
            exit(1);
        }
    }

    // if we've made this far, all eigenvalue/vector pairs have
    // been tested.

    std::cerr << fname << ": LaEigSolve() success.\n\n";

}


main(int argc, char **argv)
{

    std::cout.precision(4);
    std::cout.setf(std::ios::scientific, std::ios::floatfield);

    if (argc < 2)
    {
        std::cerr << "Usage " << argv[0] << " N " << std::endl;
        exit(1);
    }
    int N = atoi(argv[1]);

    std::cout << "Testing " << N << " x " << N << " system." << std::endl;

    TestGenEigSolve(N);

}

