//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.

//      Complex (complex precision) Lapack routines

#ifndef _DLAPACK_H_
#define _DLAPACK_H_


#ifndef _ARCH_H_
#include "arch.h"
#endif

//#include <complex.h>
#include <complex>

extern "C"
{
// *************************** Utility Routines **********************

    float F77NAME(slamch)(char *t);

    COMPLEX F77NAME(zlamch)(char *t);


    void F77NAME(zswap)(int *n, COMPLEX *x, int *incx, COMPLEX *y, int *incy);

    void F77NAME(zgesv)(int *n, int *k, COMPLEX *A, int *lda, int *ipiv,
            COMPLEX *X, int *ldx, int *info);

    void F77NAME(zposv)(char *uplo, int *m, int *k , COMPLEX *A, int *lda,
        COMPLEX *X, int *ldx, int *info);

    void F77NAME(zgels)(char *trans, int *m, int *n, int *nrhs, COMPLEX *A,
        int *lda, COMPLEX *B, int *ldb, COMPLEX *work, int *lwork, int *info);

    void F77NAME(ztimmg)(int *iflag, int *m, int *n, COMPLEX *A, int *lda,
                int *kl, int *ku);

    void F77NAME(zlaswp)(int *n, COMPLEX *A, int *lda, int *k1, int *k2,
                int *ipiv, int *incx);

    COMPLEX F77NAME(zopla)(char *subname, int *m, int *n, int *kl, int *ku,
            int *nb);

// ******************* LU Factorization Routines **********************

    void F77NAME(zgetrf)(int *m, int *n, COMPLEX *A, int *lda, int *ipiv,
                int *info);

    void F77NAME(zgetf2)(int *m, int *n, COMPLEX *A, int *lda, int *ipiv,
                int *info);

    void F77NAME(zgbtrf)(int *m, int *n, int *KL, int *KU, COMPLEX *BM,
                int *LDBM, int *ipiv, int *info);

    void F77NAME(zgttrf)(int *N, COMPLEX *DL, COMPLEX *D, COMPLEX *DU,
                COMPLEX *DU2, int *ipiv, int *info);

    void F77NAME(zpotrf)(char *UPLO, int *N, COMPLEX *SM, int *LDSM,
                int *info);

    void F77NAME(zsytrf)(char *UPLO, int *N, COMPLEX *SM, int *LDSM,
                int *ipiv, COMPLEX *WORK, int *LWORK, int *info);

    void F77NAME(zpbtrf)(char *UPLO, int *N, int *KD, COMPLEX *SBM,
                int *LDSM, int *info);

    void F77NAME(zpttrf)(int *N, COMPLEX *D, COMPLEX *E, int *info);

// ********************* LU Solve Routines ***************************

    void F77NAME(zgetrs)(char *trans, int *N, int *nrhs, COMPLEX *A, int *lda, 
            int * ipiv, COMPLEX *b, int *ldb, int *info);

    void F77NAME(zgbtrs)(char *trans, int *N, int *kl, int *ku, int *nrhs,
            COMPLEX *AB, int *ldab, int *ipiv, COMPLEX *b, int *ldb, int *info);

    void F77NAME(zsytrs)(char *uplo, int *N, int *nrhs, COMPLEX *A, int *lda, 
            int *ipiv, COMPLEX *b, int *ldb, int *info);

    void F77NAME(zgttrs)(char *trans, int *N, int *nrhs, COMPLEX *DL, 
                COMPLEX *D, COMPLEX *DU, COMPLEX *DU2, int *ipiv, COMPLEX *b, 
                int *ldb, int *info);

    void F77NAME(zpotrs)(char *UPLO, int *N, int *nrhs, COMPLEX *A, int *LDA,
                COMPLEX *b, int *ldb, int *info);

    void F77NAME(zpttrs)(int *N, int *nrhs, COMPLEX *D, COMPLEX *E, 
                COMPLEX *b, int *ldb, int *info);

    void F77NAME(zpbtrs)(char *UPLO, int *N, int *KD, int *nrhs, COMPLEX *AB,
                int *LDAB, COMPLEX *b, int *ldb, int *info);

// *******************************
}

#endif 
