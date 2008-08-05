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
//#include LA_GEN_MAT_DOUBLE_H // changed of VC++
#include "gmd.h"
//#include LA_VECTOR_LONG_INT_H // changed of VC++
#include "lavli.h"
//#include LA_SYMM_MAT_DOUBLE_H // changed of VC++
#include "symd.h"

//#include LA_UTIL_H // changed of VC++
#include "lautil.h"
#include "lapackd.h"
#include <cstring>     // strlen()

int LaEnvBlockSize(const char *fname, const LaGenMatDouble &A)
{
    char *opts = "U";

    int one = 1;
    int M = A.size(0);
    int N = A.size(1);
    int junk = -1;

    return ilaenv_(&one, fname, opts, &M, &N, &junk , &junk,
        strlen(fname), strlen(opts));
}

int LaEnvBlockSize(const char *fname, const LaSymmMatDouble &A)
{
    char *opts = "U";
    int one = 1;
    int M = A.size(0);
    int N = A.size(1);
    int junk = -1;

    return ilaenv_(&one, fname, opts, &M, &N, &junk , &junk,
        strlen(fname), strlen(opts));
}


double Mach_eps_double()
{
    char e= 'e';
    return F77NAME(dlamch)(&e);
}


void LaSwap(LaGenMatDouble &A, LaVectorLongInt &ipiv)
{
    long int lda = A.gdim(0),  n = A.size(1);
    long int k1 = ipiv.start(), k2 = ipiv.end(), incx = ipiv.inc();

    F77NAME(dlaswp)(&n, &A(0,0), &lda, &k1, &k2, &ipiv(0), &incx);
}


