//      LAPACK++ (V. 1.1)

#ifndef _BLAS_PP_H_
#define _BLAS_PP_H_

// requires
//

#include "laexcp.h"
#include "blas1++.h"
#include "blas2++.h"
#include "blas3++.h"
#include <math.h>

//double fabs(double);

//-------------------------------------
// Vector/Vector operators
//-------------------------------------

#ifdef _LA_VECTOR_DOUBLE_H_

inline LaVectorDouble operator*(const LaVectorDouble &x, double a)
{
    int N = x.size();
    LaVectorDouble t(N);

    for (int i=0; i<N; i++)
    {
        t(i) = a * x(i);
    }

    return t;
}
    
inline LaVectorDouble operator*(double a, const LaVectorDouble &x)
{
    return operator*(x,a);
}

inline double operator*(const LaVectorDouble &dx, 
                                const LaVectorDouble &dy)
{
    assert(dx.size()==dy.size());
    integer incx = dx.inc(), incy = dy.inc(), n = dx.size();

    return F77NAME(ddot)(&n, &dx(0), &incx, &dy(0), &incy);
}
                      
inline LaVectorDouble operator+(const LaVectorDouble &dx, 
                                const LaVectorDouble &dy)
{
    assert(dx.size()==dy.size());
    integer incx = dx.inc(), incy = dx.inc(), n = dx.size();
    double da = 1.0;

    LaVectorDouble tmp((int) n);
    tmp = dy;

    F77NAME(daxpy)(&n, &da, &dx(0), &incx, &tmp(0), &incy);
    return tmp;
}

inline LaVectorDouble operator-(const LaVectorDouble &dx, 
                                const LaVectorDouble &dy)
{
    assert(dx.size()==dy.size());
    integer incx = dx.inc(), incy = dy.inc(), n = dx.size();
    double da = -1.0;

    LaVectorDouble tmp(n);
    tmp = dx;

    F77NAME(daxpy)(&n, &da, &dy(0), &incx, &tmp(0), &incy);
    return tmp;
}

//-------------------------------------
// Matrix/Vector operators
//-------------------------------------

#ifdef _LA_GEN_MAT_DOUBLE_H_
inline LaVectorDouble operator*(const LaGenMatDouble &A, 
                                const LaVectorDouble &dx)
{
    char trans = 'N';
    double alpha = 1.0, beta = 0.0;
    integer M = A.size(0), N = A.size(1), lda = A.gdim(0);

    LaVectorDouble dy(M);
    integer incx = dx.inc();
    integer incy = dy.inc();

    dy = 0.0;

    F77NAME(dgemv)(&trans, &M, &N, &alpha, &A(0,0), &lda, &dx(0), &incx, 
        &beta, &dy(0), &incy); 
    return dy; 
        
}
#endif

#ifdef _LA_BAND_MAT_DOUBLE_H_
inline LaVectorDouble operator*(const LaBandMatDouble &A, 
                                const LaVectorDouble &dx)
{
    char trans = 'N';
    double alpha = 1.0, beta = 0.0;
    integer M = A.size(0), N = A.size(1), lda = A.gdim(0),
        kl = A.subdiags(), ku = A.superdiags(); 

    LaVectorDouble dy(N);
    integer incx = dx.inc(), incy = dy.inc();

    F77NAME(dgbmv)(&trans, &M, &N, &kl, &ku, &alpha, &A(0,0), &lda,
                   &dx(0), &incx, &beta, &dy(0), &incy);
    return dy;
}
#endif

#ifdef _LA_SYMM_MAT_DOUBLE_H_
inline LaVectorDouble operator*(const LaSymmMatDouble &A, 
                                const LaVectorDouble &dx)
{
    char uplo = 'L';
    double alpha = 1.0, beta = 0.0;
    integer N = A.size(1), lda = A.gdim(0);

    LaVectorDouble dy(N);
    integer incx = dx.inc(), incy = dy.inc();

    F77NAME(dsymv)(&uplo, &N, &alpha, &A(0,0), &lda, &dx(0), &incx,
                   &beta, &dy(0), &incy);
    return dy;
}
#endif

#ifdef _LA_SYMM_BAND_MAT_DOUBLE_H_ 
inline LaVectorDouble operator*(const LaSymmBandMatDouble &A, 
        const LaVectorDouble &dx) 
{
    char uplo = 'L';
    double alpha = 1.0, beta = 0.0;
    integer N = A.size(1), lda = A.gdim(0), k = A.subdiags();

    LaVectorDouble dy(N);
    integer incx = dx.inc(), incy = dy.inc();

    F77NAME(dsbmv)(&uplo, &N, &k, &alpha, &A(0,0), &lda, &dx(0), &incx,
                   &beta, &dy(0), &incy);
    return dy;
}

#endif


#ifdef _LA_SPD_MAT_DOUBLE_H_
inline LaVectorDouble operator*(const LaSpdMatDouble &AP, 
                                const LaVectorDouble &dx)
{
    char uplo = 'L';
    double alpha = 1.0, beta = 0.0;
    integer N = AP.size(1), incx = dx.inc(); 
    integer lda = AP.gdim(0);

    LaVectorDouble dy(N);
    integer incy = dy.inc();

    F77NAME(dsymv)(&uplo, &N, &alpha, &AP(0,0), &lda, &dx(0), &incx, &beta,
                    &dy(0), &incy);
    return dy;
}
#endif

#ifdef _LA_LOWER_TRIANG_MAT_DOUBLE_H_
inline LaVectorDouble operator*(const LaLowerTriangMatDouble &A, 
                                const LaVectorDouble &dx)
{
    char uplo = 'L', trans = 'N', diag = 'N';
    integer N = A.size(1), lda = A.gdim(0),
        incx = dx.inc();

    LaVectorDouble dy(dx);

    F77NAME(dtrmv)(&uplo, &trans, &diag, &N, &A(0,0), &lda, &dy(0), &incx);

    return dy;
}
#endif

#ifdef _LA_UPPER_TRIANG_MAT_DOUBLE_H_
inline LaVectorDouble operator*(const LaUpperTriangMatDouble &A, 
                                const LaVectorDouble &dx)
{
    char uplo = 'U', trans = 'N', diag = 'N';
    integer N = A.size(1), lda = A.gdim(0),
        incx = dx.inc();

    LaVectorDouble dy(dx);

    F77NAME(dtrmv)(&uplo, &trans, &diag, &N, &A(0,0), &lda, &dy(0), &incx);

    return dy;
}
#endif

#ifdef _LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H_
inline LaVectorDouble operator*(const LaUnitLowerTriangMatDouble &A,
                                const LaVectorDouble &dx)
{
    char uplo = 'L', trans = 'N', diag = 'U';
    integer N = A.size(1), lda = A.gdim(0),
        incx = dx.inc();

    LaVectorDouble dy(dx);

    F77NAME(dtrmv)(&uplo, &trans, &diag, &N, &A(0,0), &lda, &dy(0), &incx);

    return dy;
}

#endif

#ifdef _LA_UNIT_UPPER_TRIANG_MAT_DOUBLE_H_
inline LaVectorDouble operator*(const LaUnitUpperTriangMatDouble &A,
                                const LaVectorDouble &dx)
{
    char uplo = 'U', trans = 'N', diag = 'U';
    integer N = A.size(1), lda = A.gdim(0),
        incx = dx.inc();

    LaVectorDouble dy(dx);

    F77NAME(dtrmv)(&uplo, &trans, &diag, &N, &A(0,0), &lda, &dy(0), &incx);

    return dy;
}
#endif


//-------------------------------------
// Matrix/Matrix operators
//-------------------------------------

inline LaGenMatDouble operator*(const LaGenMatDouble &A, 
                                const LaGenMatDouble &B)
{
    char t = 'N';
    integer m = A.size(0), k = A.size(1), n = B.size(1);
    integer lda = A.gdim(0), ldb = B.gdim(0);
    double alpha = 1.0, beta = 1.0;

    LaGenMatDouble C(m,n);
    integer ldc = A.gdim(0);

    C = 0.0;

  F77NAME(dgemm)(&t, &t, &m, &n, &k, &alpha, &A(0,0), &lda, &B(0,0), &ldb,
                &beta, &C(0,0), &ldc);

    return C;
}

#ifdef _LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H_
inline LaGenMatDouble operator*(const LaUnitLowerTriangMatDouble &A,
                                const LaGenMatDouble &B)
{
        char side = 'L', uplo = 'L', transa = 'N', diag = 'U';
        double alpha = 1.0;
        integer m = B.size(0), n = B.size(1),
                lda = A.gdim(0), ldb = B.gdim(0);

        LaGenMatDouble C(B);

  F77NAME(dtrmm)(&side, &uplo, &transa, &diag, &m, &n, &alpha,
                &A(0,0), &lda, &C(0,0), &ldb);

        return C;
}
#endif

#ifdef _LA_UNIT_UPPER_TRIANG_MAT_DOUBLE_H_
inline LaGenMatDouble operator*(const LaUnitUpperTriangMatDouble &A,
                                const LaGenMatDouble &B)
{
        char side = 'L', uplo = 'U', transa = 'N', diag = 'U';
        double alpha = 1.0;
        integer m = B.size(0), n = B.size(1),
                lda = A.gdim(0), ldb = B.gdim(0);

        LaGenMatDouble C(B);

  F77NAME(dtrmm)(&side, &uplo, &transa, &diag, &m, &n, &alpha,
                &A(0,0), &lda, &C(0,0), &ldb);

        return C;
}
#endif

#ifdef _LA_LOWER_TRIANG_MAT_DOUBLE_H_
inline LaGenMatDouble operator*(const LaLowerTriangMatDouble &A,
                                const LaGenMatDouble &B)
{
        char side = 'L', uplo = 'L', transa = 'N', diag = 'N';
        double alpha = 1.0;
        integer m = B.size(0), n = B.size(1),
                lda = A.gdim(0), ldb = B.gdim(0);

        LaGenMatDouble C(B);

  F77NAME(dtrmm)(&side, &uplo, &transa, &diag, &m, &n, &alpha,
                &A(0,0), &lda, &C(0,0), &ldb);

        return C;
}
#endif

#ifdef _LA_UPPER_TRIANG_MAT_DOUBLE_H_
inline LaGenMatDouble operator*(const LaUpperTriangMatDouble &A,
                                const LaGenMatDouble &B)
{
        char side = 'L', uplo = 'U', transa = 'N', diag = 'N';
        double alpha = 1.0;
        integer m = B.size(0), n = B.size(1),
                lda = A.gdim(0), ldb = B.gdim(0);

        LaGenMatDouble C(B);

  F77NAME(dtrmm)(&side, &uplo, &transa, &diag, &m, &n, &alpha,
                &A(0,0), &lda, &C(0,0), &ldb);

        return C;
}
#endif

#ifdef _LA_SYMM_MAT_DOUBLE_H_
inline LaGenMatDouble operator*(const LaSymmMatDouble &A, 
                                const LaGenMatDouble &B)
{
        char side = 'L', uplo = 'L';
        double alpha = 1.0, beta = 1.0;
        LaGenMatDouble C(B.size(1),A.size(1));
        integer m = C.size(0), n = C.size(1), lda = A.gdim(0),
                ldb = B.gdim(0), ldc = C.gdim(0);

  F77NAME(dsymm)(&side, &uplo, &m, &n, &alpha, &A(0,0), &lda,
                &B(0,0), &ldb, &beta, &C(0,0), &ldc);

        return C;
}
#endif

#ifdef _LA_SYMM_TRIDIAG_MAT_DOUBLE_H_
inline LaVectorDouble operator*(const LaSymmTridiagMatDouble& A, 
                                const LaVectorDouble& X)
{
    integer M = A.size();
    integer N = X.size();
    LaVectorDouble R(M);

    R(0) = ((A.diag(0)(0) * X(0)) + (A.diag(-1)(0) * X(1)));

    for (integer i = 1; i < M-2; i++)
    {
        R(i) = ((A.diag(-1)(i-1) * X(i-1)) +
                (A.diag(0)(i) * X(i)) +
                (A.diag(-1)(i) * X(i+1)));
    }

    R(M-1) = ((A.diag(0)(M-1) * X(N-1)) + (A.diag(-1)(M-2) * X(N-2)));

    return R;
}
#endif

#ifdef  _LA_TRIDIAG_MAT_DOUBLE_H_
inline LaVectorDouble operator*(const LaTridiagMatDouble& A, 
                                const LaVectorDouble& X)
{
    integer M = A.size();
    integer N = X.size();
    LaVectorDouble R(M);

    R(0) = ((A.diag(0)(0) * X(0)) + (A.diag(-1)(0) * X(1)));

    for (integer i = 1; i < M-2; i++)
    {
        R(i) = ((A.diag(-1)(i-1) * X(i-1)) +
                (A.diag(0)(i) * X(i)) +
                (A.diag(1)(i) * X(i+1)));
    }

    R(M-1) = ((A.diag(0)(M-1) * X(N-1)) + (A.diag(1)(M-2) * X(N-2)));

    return R;
}
#endif

inline LaGenMatDouble operator-(const LaGenMatDouble &A, 
                                const LaGenMatDouble &B)
{
#ifndef HPPA
    const char fname[] = "operator+(A,B)";
#else
    char *fname = NULL;
#endif

    integer M = A.size(0);
    integer N = A.size(1);

    if (M != B.size(0) || N != B.size(1))
    {
        throw(LaException(fname, "matrices non-conformant."));
    }

    LaGenMatDouble C(M,N);

    // slow mode
    // we'll hook the BLAS in later

    for (integer i=0;  i<M; i++)
        for(integer j=0; j<N; j++)
            C(i,j) = A(i,j) - B(i,j);

    return C;
}

inline LaGenMatDouble operator+(const LaGenMatDouble &A, 
                                const LaGenMatDouble &B)
{
#ifndef HPPA
    const char fname[] = "operator+(A,B)";
#else
    char *fname = NULL;
#endif

    integer M = A.size(0);
    integer N = A.size(1);

    if (M != B.size(0) || N != B.size(1))
    {
        throw(LaException(fname, "matrices non-conformant."));
    }

    LaGenMatDouble C(M,N);

    // slow mode
    // we'll hook the BLAS in later

    for (integer i=0;  i<M; i++)
        for(integer j=0; j<N; j++)
            C(i,j) = A(i,j) + B(i,j);

    return C;
}


//-------------------------------------
// Matrix/Vector Norms
//-------------------------------------

inline double Norm_Inf(const LaVectorDouble &x)
{   

    integer index = Blas_Index_Max(x);
    return fabs(x(index));
}

inline double Norm_Inf(const LaGenMatDouble &A)
{
    integer M=A.size(0);
    integer index;

    // max row-sum

    LaVectorDouble R(M);

    for (integer i=0; i<M; i++)
        R(i) = Blas_Norm1( A( i, LaIndex() ));

    index = Blas_Index_Max(R);
    return R(index);
}

#ifdef _LA_BAND_MAT_DOUBLE_H_
inline double Norm_Inf(const LaBandMatDouble &A)
{
    integer kl = A.subdiags(), ku = A.superdiags(); 
    integer N=A.size(1);
    integer M=N;

    // slow version

    LaVectorDouble R(M);
    integer i;
    integer j;
    for (i=0; i<M; i++)
    {
        R(i) = 0.0;
        for (j=0; j<N; j++)
            R(i) += fabs(A(i,j));
    }

    double max = R(0);

    // report back largest row sum
    for (i=1; i<M; i++)
        if (R(i) > max) max=R(i);

    return max;
}
#endif

#ifdef _LA_SYMM_MAT_DOUBLE_H_
inline double Norm_Inf(const LaSymmMatDouble &S)
{
    integer N = S.size(0); // square matrix

    // slow version

    LaVectorDouble R(N);
    integer i; 
    integer j;      

    for (i=0; i<N; i++)
    {
        R(i) = 0.0;
        for (j=0; j<N; j++)
            R(i) += fabs(S(i,j));
    }
     
    double max = R(0);

    // report back largest row sum
    for (i=1; i<N; i++)
        if (R(i) > max) max=R(i);

    return max;
}
#endif

#ifdef _LA_SPD_MAT_DOUBLE_H_
inline double Norm_Inf(const LaSpdMatDouble &S)
{
    integer N = S.size(0); //SPD matrices are square

    // slow version

    LaVectorDouble R(N);
    integer i; 
    integer j;      

    for (i=0; i<N; i++)
    {
        R(i) = 0.0;
        for (j=0; j<N; j++)
            R(i) += fabs(S(i,j));
    }
     
    double max = R(0);

    // report back largest row sum
    for (i=1; i<N; i++)
        if (R(i) > max) max=R(i);

    return max;
}
#endif

#ifdef _LA_SYMM_TRIDIAG_MAT_DOUBLE_H_
inline double Norm_Inf(const LaSymmTridiagMatDouble &S)
{
    integer N = S.size();   // S is square
    LaVectorDouble R(N);

    R(0) = fabs(S(0,0)) + fabs(S(0,1));

    for (integer i=1; i<N-1; i++)
    {
        R(i) = fabs(S(i,i-1)) + fabs(S(i,i)) + fabs(S(i,i+1));
    }

    R(N-1) = fabs(S(N-1,N-2)) + fabs(S(N-1,N-1));

    return Norm_Inf(R);
}
#endif

#ifdef _LA_TRIDIAG_MAT_DOUBLE_H_
inline double Norm_Inf(const LaTridiagMatDouble &T)
{
    integer N = T.size();   // T is square
    LaVectorDouble R(N);

    R(0) = fabs(T(0,0)) + fabs(T(0,1));

    for (int i=1; i<N-1; i++)
    {
        R(i) = fabs(T(i,i-1)) + fabs(T(i,i)) + fabs(T(i,i+1));
    }

    R(N-1) = fabs(T(N-1,N-2)) + fabs(T(N-1,N-1));

    return Norm_Inf(R);
}
#endif


#endif
    // LA_VECTOR_DOUBLE_H
#endif 
    // _BLAS_PP_H_
