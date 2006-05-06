#include "lafnames.h"
//#include LA_VECTOR_DOUBLE_H // changed of VC++
#include "lavd.h"
//#include LA_SYMM_MAT_DOUBLE_H // changed of VC++
#include "symd.h"
//#include LA_UNIT_UPPER_TRIANG_MAT_DOUBLE_H // changed of VC++
#include "uutgmd.h"
//#include LA_UPPER_TRIANG_MAT_DOUBLE_H // changed of VC++
#include "utgmd.h"
//#include LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H // changed of VC++
#include "ultgmd.h"
//#include LA_LOWER_TRIANG_MAT_DOUBLE_H // changed of VC++
#include "ltgmd.h"
//#include LA_SPD_MAT_DOUBLE_H // changed of VC++
#include "spdmd.h"
//#include LA_SYMM_BAND_MAT_DOUBLE_H // changed of VC++
#include "sybmd.h"
//#include LA_TRIDIAG_MAT_DOUBLE_H // changed of VC++
#include "trmd.h"

#include "blas3++.h"


void Blas_Mat_Mat_Mult(const LaGenMatDouble &A, 
            const LaGenMatDouble &B, LaGenMatDouble &C, 
            double alpha , double beta )
{
        char t = 'N';
        integer m = A.size(0), k = A.size(1), n = B.size(1);
        integer lda = A.gdim(0), ldb = B.gdim(0), ldc = C.gdim(0);

  F77NAME(dgemm)(&t, &t, &m, &n, &k, &alpha, &A(0,0), &lda, &B(0,0), 
                &ldb, &beta, &C(0,0), &ldc);
    
}

void Blas_Mat_Trans_Mat_Mult(const LaGenMatDouble &A, 
            const LaGenMatDouble &B, LaGenMatDouble &C, 
            double alpha, double beta )
{
        char nt = 'N';
        char t = 'T';
        integer m = A.size(0), k = A.size(1), n = B.size(1);
        integer lda = A.gdim(0), ldb = B.gdim(0), ldc = C.gdim(0);

  F77NAME(dgemm)(&t, &nt, &m, &n, &k, &alpha, &A(0,0), &lda, &B(0,0), 
                &ldb, &beta, &C(0,0), &ldc);
    
}

void Blas_Mat_Mat_Trans_Mult(const LaGenMatDouble &A, 
            const LaGenMatDouble &B, LaGenMatDouble &C, 
            double alpha,  double beta )
{
        char nt = 'N';
        char t = 'T';
        integer m = A.size(0), k = A.size(1), n = B.size(1);
        integer lda = A.gdim(0), ldb = B.gdim(0), ldc = C.gdim(0);

  F77NAME(dgemm)(&nt, &t, &m, &n, &k, &alpha, &A(0,0), &lda, &B(0,0), 
                &ldb, &beta, &C(0,0), &ldc);
    
}









#ifdef _LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H_

void Blas_Mat_Mat_Solve(LaUnitLowerTriangMatDouble &A, 
            LaGenMatDouble &B, double alpha)
{
        char side = 'L', uplo = 'L', transa = 'N', diag = 'U';
        integer m = B.size(0), n = B.size(1), 
                lda = A.gdim(0), ldb = B.gdim(0);

  F77NAME(dtrsm)(&side, &uplo, &transa, &diag, &m, &n, &alpha, 
                &A(0,0), &lda, &B(0,0), &ldb);
}

#endif

#ifdef _LA_UNIT_UPPER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Mat_Mult(LaUnitUpperTriangMatDouble &A,
            LaGenMatDouble &B, double alpha )
{
        char side = 'L', uplo = 'U', transa = 'N', diag = 'U';
        integer m = B.size(0), n = B.size(1),
                lda = A.gdim(0), ldb = B.gdim(0);

  F77NAME(dtrmm)(&side, &uplo, &transa, &diag, &m, &n, &alpha,
                &A(0,0), &lda, &B(0,0), &ldb);
}

void Blas_Mat_Mat_Solve(LaUnitUpperTriangMatDouble &A, 
            LaGenMatDouble &B, double alpha )
{
        char side = 'L', uplo = 'U', transa = 'N', diag = 'U';
        integer m = B.size(0), n = B.size(1), 
                lda = A.gdim(0), ldb = B.gdim(0);

  F77NAME(dtrsm)(&side, &uplo, &transa, &diag, &m, &n, &alpha, 
                &A(0,0), &lda, &B(0,0), &ldb);
}

#endif

#ifdef _LA_LOWER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Mat_Mult(LaLowerTriangMatDouble &A,
            LaGenMatDouble &B, double alpha )
{
        char side = 'L', uplo = 'L', transa = 'N', diag = 'N';
        integer m = B.size(0), n = B.size(1),
                lda = A.gdim(0), ldb = B.gdim(0);

  F77NAME(dtrmm)(&side, &uplo, &transa, &diag, &m, &n, &alpha,
                &A(0,0), &lda, &B(0,0), &ldb);
}

void Blas_Mat_Mat_Solve(LaLowerTriangMatDouble &A, 
            LaGenMatDouble &B, double alpha )
{
        char side = 'L', uplo = 'L', transa = 'N', diag = 'N';
        integer m = B.size(0), n = B.size(1), 
                lda = A.gdim(0), ldb = B.gdim(0);

  F77NAME(dtrsm)(&side, &uplo, &transa, &diag, &m, &n, &alpha, 
                &A(0,0), &lda, &B(0,0), &ldb);
}
#endif


#ifdef _LA_UPPER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Mat_Mult(LaUpperTriangMatDouble &A,
            LaGenMatDouble &B, double alpha )
{
        char side = 'L', uplo = 'U', transa = 'N', diag = 'N';
        integer m = B.size(0), n = B.size(1),
                lda = A.gdim(0), ldb = B.gdim(0);

  F77NAME(dtrmm)(&side, &uplo, &transa, &diag, &m, &n, &alpha,
                &A(0,0), &lda, &B(0,0), &ldb);
}

void Blas_Mat_Mat_Solve(LaUpperTriangMatDouble &A, 
            LaGenMatDouble &B, double alpha )
{
        char side = 'L', uplo = 'U', transa = 'N', diag = 'N';
        integer m = B.size(0), n = B.size(1), 
                lda = A.gdim(0), ldb = B.gdim(0);

  F77NAME(dtrsm)(&side, &uplo, &transa, &diag, &m, &n, &alpha, 
                &A(0,0), &lda, &B(0,0), &ldb);
}
#endif


#ifdef _LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Trans_Mat_Solve(LaUnitLowerTriangMatDouble &A,
            LaGenMatDouble &B, double alpha )
{
        char side = 'L', uplo = 'L', transa = 'T', diag = 'U';
        integer m = B.size(0), n = B.size(1),
                lda = A.gdim(0), ldb = B.gdim(0);

  F77NAME(dtrsm)(&side, &uplo, &transa, &diag, &m, &n, &alpha,
                &A(0,0), &lda, &B(0,0), &ldb);
}
#endif

#ifdef _LA_UNIT_UPPER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Trans_Mat_Solve(LaUnitUpperTriangMatDouble &A,
            LaGenMatDouble &B, double alpha )
{
        char side = 'L', uplo = 'U', transa = 'T', diag = 'U';
        integer m = B.size(0), n = B.size(1),
                lda = A.gdim(0), ldb = B.gdim(0);

  F77NAME(dtrsm)(&side, &uplo, &transa, &diag, &m, &n, &alpha,
                &A(0,0), &lda, &B(0,0), &ldb);
}
#endif

#ifdef LA_LOWER_TRIANG_MAT_DOUBLE_H
void Blas_Mat_Mat_Mult(LaUnitLowerTriangMatDouble &A,
            LaGenMatDouble &B, double alpha )
{
        char side = 'L', uplo = 'L', transa = 'N', diag = 'U';
        integer m = B.size(0), n = B.size(1),
                lda = A.gdim(0), ldb = B.gdim(0);

  F77NAME(dtrmm)(&side, &uplo, &transa, &diag, &m, &n, &alpha,
                &A(0,0), &lda, &B(0,0), &ldb);
}

void Blas_Mat_Trans_Mat_Solve(LaLowerTriangMatDouble &A,
            LaGenMatDouble &B, double alpha )
{
        char side = 'L', uplo = 'L', transa = 'T', diag = 'N';
        integer m = B.size(0), n = B.size(1),
                lda = A.gdim(0), ldb = B.gdim(0);

  F77NAME(dtrsm)(&side, &uplo, &transa, &diag, &m, &n, &alpha,
                &A(0,0), &lda, &B(0,0), &ldb);
}
#endif


#ifdef _LA_UPPER_TRIANG_MAT_DOUBLE_H 
void Blas_Mat_Trans_Mat_Solve(LaUpperTriangMatDouble &A,
            LaGenMatDouble &B, double alpha )
{
        char side = 'L', uplo = 'U', transa = 'T', diag = 'N';
        integer m = B.size(0), n = B.size(1),
                lda = A.gdim(0), ldb = B.gdim(0);

  F77NAME(dtrsm)(&side, &uplo, &transa, &diag, &m, &n, &alpha,
                &A(0,0), &lda, &B(0,0), &ldb);
}

#endif

#ifdef _LA_SYMM_MAT_DOUBLE_H_
void Blas_Mat_Mat_Mult(LaSymmMatDouble &A, LaGenMatDouble &B, 
            LaGenMatDouble &C, double alpha , double beta )
{
        char side = 'L', uplo = 'L';
        integer m = C.size(0), n = C.size(1), lda = A.gdim(0), 
                ldb = B.gdim(0), ldc = C.gdim(0);

  F77NAME(dsymm)(&side, &uplo, &m, &n, &alpha, &A(0,0), &lda, 
                &B(0,0), &ldb, &beta, &C(0,0), &ldc);
}

void Blas_R1_Update(LaSymmMatDouble &C, LaGenMatDouble &A,
            double alpha , double beta )
{
        char uplo = 'L', transa = 'N';
        integer n = C.size(0), k = A.size(1), 
                lda = A.gdim(0), ldc = C.gdim(0);

  F77NAME(dsyrk)(&uplo, &transa, &n, &k, &alpha, &A(0,0), &lda, 
                &beta, &C(0,0), &ldc);
}


void Blas_R2_Update(LaSymmMatDouble &C, LaGenMatDouble &A,
            LaGenMatDouble &B, double alpha , double beta )
{
        char uplo = 'L', transa = 'N';
        integer n = C.size(0), k = A.size(1), lda = A.gdim(0), 
                ldb = B.gdim(0), ldc = C.gdim(0);

  F77NAME(dsyr2k)(&uplo, &transa, &n, &k, &alpha, &A(0,0), &lda, 
                &B(0,0), &ldb, &beta, &C(0,0), &ldc);
}

#endif
