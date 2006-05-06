//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.

//      Double precision Lapack routines

#ifndef _DLAPACK_H_
#define _DLAPACK_H_


#ifndef _ARCH_H_
#include "arch.h"
#endif

#include "f2c.h"

extern "C"
{
// *************************** Utility Routines **********************


    double F77NAME(dlamch)(char *t);



//******************** Linear Equation Solvers *************************
    void F77NAME(dgesv)(integer *n, integer *k, doublereal *A, integer *lda, integer *ipiv,
            doublereal *X, integer *ldx, integer *info);

    void F77NAME(dposv)(char *uplo, integer *m, integer *k , doublereal *A, integer *lda,
        doublereal *X, integer *ldx, integer *info);

    void F77NAME(dsysv)(const char *uplo, integer *n, integer *nrhs, doublereal *A, 
            integer *lda, integer *ipv, integer *lidb, doublereal *work, integer *lwork, integer *info);

//******************** Lapack Utility Routines ************************

    void F77NAME(dgels)(char *trans, integer *m, integer *n, integer *nrhs, doublereal *A,
        integer *lda, doublereal *B, integer *ldb, doublereal *work, integer *lwork, integer *info);

    void F77NAME(dtimmg)(integer *iflag, integer *m, integer *n, doublereal *A, integer *lda,
                integer *kl, integer *ku);

    void F77NAME(dlaswp)(integer *n, doublereal *A, integer *lda, integer *k1, integer *k2,
                integer *ipiv, integer *incx);

    doublereal F77NAME(dopla)(char *subname, integer *m, integer *n, integer *kl, integer *ku,
            integer *nb);

// ******************* LU Factorization Routines **********************

    void F77NAME(dgetrf)(integer *m, integer *n, doublereal *A, integer *lda, integer *ipiv,
                integer *info);

    void F77NAME(dgetf2)(integer *m, integer *n, doublereal *A, integer *lda, integer *ipiv,
                integer *info);

    void F77NAME(dgbtrf)(integer *m, integer *n, integer *KL, integer *KU, doublereal *BM,
                integer *LDBM, integer *ipiv, integer *info);

    void F77NAME(dgttrf)(integer *N, doublereal *DL, doublereal *D, doublereal *DU,
                doublereal *DU2, integer *ipiv, integer *info);

    void F77NAME(dpotrf)(char *UPLO, integer *N, doublereal *SM, integer *LDSM,
                integer *info);

    void F77NAME(dsytrf)(char *UPLO, integer *N, doublereal *SM, integer *LDSM,
                integer *ipiv, doublereal *WORK, integer *LWORK, integer *info);

    void F77NAME(dpbtrf)(char *UPLO, integer *N, integer *KD, doublereal *SBM,
                integer *LDSM, integer *info);

    void F77NAME(dpttrf)(integer *N, doublereal *D, doublereal *E, integer *info);

// ********************* LU Solve Routines ***************************

    void F77NAME(dgetrs)(char *trans, integer *N, integer *nrhs, doublereal *A, integer *lda, 
            integer * ipiv, doublereal *b, integer *ldb, integer *info);

    void F77NAME(dgbtrs)(char *trans, integer *N, integer *kl, integer *ku, integer *nrhs,
            doublereal *AB, integer *ldab, integer *ipiv, doublereal *b, integer *ldb, integer *info);

    void F77NAME(dsytrs)(char *uplo, integer *N, integer *nrhs, doublereal *A, integer *lda, 
            integer *ipiv, doublereal *b, integer *ldb, integer *info);

    void F77NAME(dgttrs)(char *trans, integer *N, integer *nrhs, doublereal *DL, 
                doublereal *D, doublereal *DU, doublereal *DU2, integer *ipiv, doublereal *b, 
                integer *ldb, integer *info);

    void F77NAME(dpotrs)(char *UPLO, integer *N, integer *nrhs, doublereal *A, integer *LDA,
                doublereal *b, integer *ldb, integer *info);

    void F77NAME(dpttrs)(integer *N, integer *nrhs, doublereal *D, doublereal *E, 
                doublereal *b, integer *ldb, integer *info);

    void F77NAME(dpbtrs)(char *UPLO, integer *N, integer *KD, integer *nrhs, doublereal *AB,
                integer *LDAB, doublereal *b, integer *ldb, integer *info);

// ********************* Eigen Solve Routines ***************************

    void F77NAME(dsyev)(char *jobz, char *uplo, integer *N, doublereal *S,
    integer *lda, doublereal *eig, doublereal *work, integer *lwork, integer *info);

// *******************************
}

#endif 
