/*      LAPACK++ (V. 1.1)                                       */
/*      (C) 1992-1996 All Rights Reserved.                          */

#ifndef _BLAS2_H_
#define _BLAS2_H_

#include "arch.h"
#include "f2c.h"

extern "C"
{
     void F77NAME(dgemv)(char* trans, integer* M, integer* N, double* alpha, 
                    const double* A, integer* lda, const double* dx, 
                    integer* incx, double* beta, double* dy, integer* incy);

     void F77NAME(dgbmv)(char* trans, integer* M, integer* N, integer* kl, 
                    integer* ku, double* alpha, const double* A, integer* lda, 
                    const double* dx, integer* incx, double* beta, 
                    double* dy, integer* incy);

     void F77NAME(dsymv)(char* uplo, integer* N, double* alpha, const double* A, 
                    integer* lda, const double* dx, integer* incx, double* beta,
                    double* dy, integer* incy);

     void F77NAME(dsbmv)(char* uplo, integer* N, integer* k, double* alpha, 
                    const double* A, integer* lda, const double* dx, 
                    integer* incx, double* beta, double* dy, integer* incy);

     void F77NAME(dspmv)(char* uplo, integer* N, double* alpha, double* AP, 
                    double* dx, integer* incx, double* beta, double* dy, 
                    integer* incy);

     void F77NAME(dtrmv)(char* uplo, char* trans, char* diag, const integer* N, 
                    const double* A, integer* lda, const double* dx, 
                    integer* incx);

     // currently not implemented.
     //F77NAME(dtbmv) ( UPLO, TRANS, DIAG, N, K, A, LDA, dx, INCX )

     void F77NAME(dtrsv)(char* uplo, char* trans, char* diag, const integer* N, 
                    double* A, integer* lda, double* dx, integer* incx);

     // currently not implemented.
     //F77NAME(dtbsv) ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )

     // currently not implemented.
     //F77NAME(dtpsv) ( UPLO, TRANS, DIAG, N, AP, X, INCX )

     void F77NAME(dger)(integer* M, integer* N, double* alpha, 
                     double* dx, integer* incx, double* dy, integer* incy, 
                     double* A, integer* lda);

     void F77NAME(dsyr)(char* uplo, integer* N, double* alpha, double* dx, 
                   integer* incx, double* A, integer* lda);

     void F77NAME(dspr)(char* uplo, integer* N, double* alpha, double* dx, 
                   integer* incx, double* AP);

     void F77NAME(dsyr2)(char* uplo, integer* N, double* alpha, double* dx, 
                    integer* incx, double* dy, integer* incy, double* A, 
                    integer* lda);

     void F77NAME(dspr2)(char* uplo, integer* N, double* alpha, double* dx, 
                    integer* incx, double* dy, integer* incy, double* AP);

}

#endif
// _BLAS2_H_
