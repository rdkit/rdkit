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
#include "lapack.h"
//#include LA_GEN_MAT_DOUBLE_H // changed of VC++
#include "gmd.h"
//#include LA_VECTOR_DOUBLE_H // changed of VC++
#include "lavd.h"
//#include LA_VECTOR_INT_H // changed of VC++
#include "lavi.h"
//#include LA_SPD_MAT_DOUBLE_H // changed of VC++
#include "spdmd.h"
//#include LA_SYMM_MAT_DOUBLE_H // changed of VC++
#include "symd.h"
//#include LA_EXCEPTION_H // changed of VC++
#include "laexcp.h"
//#include LA_SOLVE_DOUBLE_H // changed of VC++
#include "laslv.h"
//#include LA_UTIL_H // changed of VC++
#include "lautil.h"

void LaEigSolve(const LaSymmMatDouble &S, LaVectorDouble &eigvals)
{   
//#ifndef HPPA
//     const char fname[] = "LaEigSolve(LaSymmMatDouble &A, &eigvals)";
//#else
//    char *fname = NULL;  // HP C++ does not support string initalization!
//#endif

    LaSymmMatDouble tmp(S);

    LaEigSolveIP(tmp, eigvals);


}

void LaEigSolve(const LaSymmMatDouble &S, LaVectorDouble &eigvals, 
    LaGenMatDouble &eigvec)
{   
//#ifndef HPPA
//   const char fname[] = "LaEigSolve(LaSymmMatDouble &A, &eigvals, &eigvecs)";
//#else
//    char *fname = NULL;  // HP C++ does not support string initalization!
//#endif

    // view symmetric matrix S as a LaGenMatDouble
    // assumes S is contiguous...
    LaSymmMatDouble tmp(S);
    LaEigSolveVecIP(tmp, eigvals);

    LaGenMatDouble tmp2(&tmp(0,0), S.size(0), S.size(1));
    eigvec.ref(tmp2);


}

void LaEigSolveIP(LaSymmMatDouble &S, LaVectorDouble &eigvals)
{   
#ifndef HPPA
     const char fname[] = "LaEigSolveIP(LaGenMatDouble &A, &v)";
#else
    char *fname = NULL;  // HP C++ does not support string initalization!
#endif
    long int N = S.size(0);
    char jobz = 'N';
    char uplo = 'L';
    long int info;
    long int lda = S.gdim(0);

    if (eigvals.size() < N)
    {
        throw(LaException(fname, "Not enough room to store eigenvalues"));
    }
        

    long int w = (LaEnvBlockSize("SSYTRD", S) +2) * N;
    LaVectorDouble Work(w);


    F77NAME(dsyev)(&jobz, &uplo, &N, S.addr(), &lda, &eigvals(0), &Work(0),
        &w, &info);

    if (info != 0)
        throw(LaException(fname, "Internal error in LAPACK: SSYEV()"));

}

void LaEigSolveVecIP(LaSymmMatDouble &S, LaVectorDouble &eigvals)
{   
#ifndef HPPA
     const char fname[] = "LaEigSolveVecIP(LaGenMatDouble &A, &eigvals)";
#else
    char *fname = NULL;  // HP C++ does not support string initalization!
#endif
    long int N = S.size(0);
    char jobz = 'V';
    char uplo = 'L';
    long int info;
    long int lda = S.gdim(0);

    if (eigvals.size() < N)
    {
        throw(LaException(fname, "Not enough room to store eigenvalues"));
    }
        

    long int w = (LaEnvBlockSize("SSYTRD", S) +2) * N;
    LaVectorDouble Work(w);


    F77NAME(dsyev)(&jobz, &uplo, &N, S.addr(), &lda, &eigvals(0), &Work(0),
        &w, &info);

    if (info != 0)
        throw(LaException(fname, "Internal error in LAPACK: SSYEV()"));

}


