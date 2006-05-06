//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.

#ifndef _LA_FACTORIZE_MAT_DOUBLE_H
#define _LA_FACTORIZE_MAT_DOUBLE_H

#include "lapack.h"

#include "lafnames.h"
//#include LA_GEN_MAT_DOUBLE_H // changed of VC++
#include "gmd.h"
//#include LA_GEN_FACT_DOUBLE_H // changed of VC++
#include "gfd.h"
//#include LA_BAND_MAT_DOUBLE_H // changed of VC++
#include "bmd.h"
//#include LA_BAND_FACT_DOUBLE_H // changed of VC++
#include "bfd.h"
//#include LA_TRIDIAG_MAT_DOUBLE_H // changed of VC++
#include "trmd.h"
//#include LA_TRIDIAG_FACT_DOUBLE_H // changed of VC++
#include "trfd.h"
//#include LA_SPD_MAT_DOUBLE_H // changed of VC++
#include "spdmd.h"
//#include LA_SYMM_MAT_DOUBLE_H // changed of VC++
#include "symd.h"
//#include LA_SYMM_FACT_DOUBLE_H // changed of VC++
#include "syfd.h"
//#include LA_SYMM_TRIDIAG_MAT_DOUBLE_H // changed of VC++
#include "sytrmd.h"
//#include LA_SYMM_BAND_MAT_DOUBLE_H // changed of VC++
#include "sybmd.h"



inline void LaGenMatFactorize(LaGenMatDouble &GM, LaGenFactDouble &GF)
{
    integer m = GM.size(0), n = GM.size(1), lda = GM.gdim(0);
    integer info=0;

    F77NAME(dgetrf)(&m, &n, &GM(0,0), &lda, &(GF.pivot()(0)), &info);
}


inline void LaGenMatSolve(LaGenMatDouble &A, LaVectorDouble &b,
                         LaVectorLongInt &piv)
{
    char trans = 'N';
    integer n = A.size(1), lda = A.gdim(0), nrhs = 1, 
            ldb = b.size(), info = 0;

    F77NAME(dgetrs)(&trans, &n, &nrhs, &A(0,0), &lda, &piv(0), 
                    &b(0), &ldb, &info);
}

#if 0
inline int LaBandMatFactorize(LaBandMatDouble &A, LaBandFactDouble &AF)
{
    integer n = A.size(1), m = n, LDA = A.gdim(0);
    integer KL = A.subdiags(), KU = A.superdiags(), info=0;

    F77NAME(dgbtrf)(&m, &n, &KL, &KU, &A(0,0), &LDA, &(AF.pivot())(0), &info);
    return 1;
}


inline void LaBandMatSolve(LaBandMatDouble &AB, LaVectorDouble &b,
                        LaVectorLongInt &piv)
{
    char trans = 'N';
    integer n = AB.size(1), ldab = AB.gdim(0), nrhs = 1, ldb = b.size();
    integer kl = AB.subdiags(), ku = AB.superdiags(), info=0;

    F77NAME(dgbtrs)(&trans, &n, &kl, &ku, &nrhs, &AB(0,0), &ldab, &piv(0), 
                    &b(0), &ldb, &info);
}


inline void LaTridiagMatFactorize(LaTridiagMatDouble &TD,
                                 LaTridiagFactDouble &TDF)
{
    integer N = TD.size(), info = 0;
    double *DL = & TD.diag(-1)(0), *D = &TD.diag(0)(0),
         *DU =  &TD.diag(1)(0), *DU2 = &TD.diag(2)(0);

    F77NAME(dgttrf)(&N, DL, D, DU, DU2, &(TDF.pivot())(0), &info);
}


inline void LaTridiagMatSolve(LaTridiagMatDouble &TD, LaVectorDouble &b,
                        LaVectorLongInt &piv)
{
    char trans = 'N';
    integer N = TD.size(), nrhs = 1, ldb = b.size(), info = 0;
    double *DL = &TD.diag(-1)(0), *D = &TD.diag(0)(0),
         *DU = &TD.diag(1)(0), *DU2 = &TD.diag(2)(0);

    F77NAME(dgttrs)(&trans, &N, &nrhs, DL, D, DU, DU2, &piv(0), 
                    &b(0), &ldb, &info);
}


inline void LaSpdMatFactorize(LaSpdMatDouble &SPD)
{
    char UPLO = 'L';
    integer N = SPD.size(0), LDSPD = SPD.gdim(0), info = 0;

    F77NAME(dpotrf)(&UPLO, &N, &SPD(0,0), &LDSPD, &info);
}


inline void LaSpdMatSolve(LaSpdMatDouble &SPD, LaVectorDouble &b)
{
    char UPLO = 'L';
    integer N = SPD.size(0), nrhs = 1, ldspd = SPD.gdim(0), ldb = b.size(),
            info = 0;

    F77NAME(dpotrs)(&UPLO, &N, &nrhs, &SPD(0,0), &ldspd, &b(0),
                    &ldb, &info);
}


inline void LaSymmMatFactorize(LaSymmMatDouble &A, LaSymmFactDouble &AF)
{
    char UPLO = 'L';
    integer N = A.size(0), LDA = A.gdim(0), info = 0;
    integer M = DSYTRF;
    int NB = F77NAME(get_nb)(&N,&M);

    integer LWORK = N*NB;
    double *WORK = new double[LWORK];

    F77NAME(dsytrf)(&UPLO, &N, &A(0,0), &LDA, &(AF.pivot()(0)), WORK, 
                    &LWORK, &info);

    delete [] WORK;
}


inline void LaSymmMatSolve(LaSymmMatDouble &A, LaVectorDouble &b,
                           LaVectorLongInt &piv)
{
    char uplo = 'L';
    integer N = A.size(1), nrhs = 1, lda = A.gdim(0),
            ldb = b.size(), info = 0;

    F77NAME(dsytrs)(&uplo, &N, &nrhs, &A(0,0), &lda,
            &piv(0), &b(0), &ldb, &info);

}

inline int LaSymmBandMatFactorize(LaSymmBandMatDouble &AB)
{
    char UPLO = 'L';
    integer N = AB.size(0), KD = AB.subdiags(), LDAB = AB.gdim(0),
            info = 0;

    F77NAME(dpbtrf)(&UPLO, &N, &KD, &AB(0,0), &LDAB, &info);
    return 1;
}


inline void LaSymmBandMatSolve(LaSymmBandMatDouble &AB, LaVectorDouble &b)
{
    char UPLO = 'L';
    integer N = AB.size(0), KD = AB.subdiags(), LDAB = AB.gdim(0),
            info = 0, nrhs = 1, ldb = b.size();

    F77NAME(dpbtrs)(&UPLO, &N, &KD, &nrhs, &AB(0,0), &LDAB, &b(0),
                    &ldb, &info);
}


inline int LaSymmTridiagMatFactorize(LaSymmTridiagMatDouble &STD)
{
    integer N = STD.size(), info = 0;
    double *D = (double*) STD.diag(0), *E = (double*) STD.diag(-1);

    F77NAME(dpttrf)(&N, D, E, &info);
    return 1;
}


inline void LaSymmTridiagMatSolve(LaSymmTridiagMatDouble &STD, 
                                    LaVectorDouble &b)
{
    integer N = STD.size(), nrhs = 1, ldb = b.size(), info = 0;
    double *D = (double*) STD.diag(0), *E = (double*) STD.diag(-1);

    F77NAME(dpttrs)(&N, &nrhs, D, E, &b(0), &ldb, &info);
}


inline void LaSwap(LaGenMatDouble &A, LaVectorLongInt &ipiv)
{
    integer lda = A.gdim(0),  n = A.size(1);
    integer k1 = ipiv.start(), k2 = ipiv.end(), incx = ipiv.inc();

    F77NAME(dlaswp)(&n, &A(0,0), &lda, &k1, &k2, ipiv, &incx);
}


inline void LaSwap_Trans(LaGenMatDouble &A, LaVectorLongInt &ipiv)
{
    integer lda = A.gdim(0),  n = A.size(1);
    integer k1 = ipiv.start(), k2 = ipiv.end(), incx = -(ipiv.inc());

    F77NAME(dlaswp)(&n, &A(0,0), &lda, &k1, &k2, ipiv, &incx);
}


inline int LaLUFactorDoubleUnblocked(LaGenMatDouble &A, LaGenFactDouble &F)
{
    integer m = A.size(0), n=A.size(1), lda = A.gdim(0);
    integer info=0;

    F77NAME(dgetf2)(&m, &n, &A(0,0), &lda, &(F.pivot()(0)), &info);
    return 1;
}


inline double LaDopla(char *subname, int m, int n, int kl, int ku, int nb)
{
    return F77NAME(dopla)(subname, &m, &n, &kl, &ku, &nb);
}

#endif
#endif 
