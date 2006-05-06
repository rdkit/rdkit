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



#include "lafnames.h"
#include "lapack.h"
//#include LA_GEN_MAT_DOUBLE_H // changed of VC++
#include "gmd.h"
//#include LA_VECTOR_DOUBLE_H // changed of VC++
#include "lavd.h"
//#include LA_SYMM_MAT_DOUBLE_H // changed of VC++
#include "symd.h"
//#include LA_SOLVE_DOUBLE_H // changed of VC++
#include "laslv.h"
//#include LA_GENERATE_MAT_DOUBLE_H // changed of VC++
#include "genmd.h"
//#include LA_EXCEPTION_H // changed of VC++
#include "laexcp.h"
#include "blas++.h"




double residual(LaSymmMatDouble &A, LaVectorDouble &x, 
    const LaVectorDouble& b)
{
    // Symm matrices are always square
    int N = A.size(0);

    cout << "\tNorm_Inf(A*x-b)" << Norm_Inf(A*x-b) << endl;
    cout << "\tNorm_Inf(A) " << Norm_Inf(A) << endl;
    cout << "\tNorm_Inf(x) " << Norm_Inf(x) << endl;
    cout << "\tMacheps :" << Mach_eps_double() << endl;

        return Norm_Inf(A*x-b ) /
                ( Norm_Inf(A)* Norm_Inf(x) * N * Mach_eps_double());
}
    
int TestSymmLinearSolve(int N)
{
    LaSymmMatDouble A(N,N);
    LaVectorDouble x(N), b(N);

    const char fname[] = "TestSymmLinearSolve(LaGenMat, x, b) ";


    LaGenerateMatDouble(A);

    // save a snapshot of what A looked like before the solution
    LaSymmMatDouble old_A(A);


    b = 1.1;

    cerr << fname << ": testing LaLinearSolve(Symm,...) \n";
    LaLinearSolve(A, x, b);
    double norm;


    if ( (norm = Norm_Inf( old_A - A)) >  0.0)  // this is a hard test, not
                                         // necessary to worry about
                                         // round-off issues.  We
                                         // are testing to see A was
                                         // overwritten.
    {
        cerr << fname << ": overwrote 1st arg.\n";
        cerr << "       error norm: " << norm << endl;
        exit(1);
    }

#if 0
    cout << "A\n" << A << endl;
    cout << "old_A\n" << old_A << endl;
    cout << "x \n" << x << endl;
#endif

    double res = residual(A,x,b);
    if (res > 1)
    {
        cerr << fname << "resdiual " << res << " is to too high.\n";
        exit(1);
    }
    cerr << fname << ": LaLinearSolve(Symm) success.\n\n";

    // now try the in-place solver

    cerr << fname << ": testing LaLinearSolveIP(Symm,...) \n";
    LaLinearSolveIP(A, x, b);
    res = residual(old_A, x, b);
    if (res > 1)
    {
        cerr << fname << "resdiual " << res << " is to too high.\n";
        exit(1);
    }
    cerr << fname << ": LaLinearSolveIP(Symm) success.\n\n";
    return 0;
}


main(int argc, char **argv)
{

    cout.precision(8);
    cout.setf(ios::scientific, ios::floatfield);

    if (argc < 2)
    {
        cerr << "Usage " << argv[0] << " M N " << endl;
        exit(1);
    }

    int M = atoi(argv[1]);
    int N = atoi(argv[2]);

    TestSymmLinearSolve(M);

}

