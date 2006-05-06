//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.


#ifndef _LA_SYMM_FACT_DOUBLE_H_
#define _LA_SYMM_FACT_DOUBLE_H_

#include "lafnames.h"
//#include LA_VECTOR_INT_H // changed of VC++
#include "lavi.h"
//#include LA_SYMM_MAT_DOUBLE_H // changed of VC++
#include "symd.h"


#include "lapack.h"

class LaSymmFactDouble
{
    LaSymmMatDouble             S_;
    LaVectorLongInt             pivot_;
    int                      info_;
    char                     uplo_;
    int                     size_;
    int                     gdim_;

public:

    // constructor

    inline LaSymmFactDouble();
    inline LaSymmFactDouble(int,int);
    inline LaSymmFactDouble(const LaSymmFactDouble &);
    inline ~LaSymmFactDouble();

    // extraction functions for components

    inline LaSymmMatDouble& S() { return S_; }
    inline LaVectorLongInt& pivot() { return pivot_; }
    inline int info() { return info_; }
    inline char uplo(){ return uplo_; }
    inline int size() { return size_; }
    inline int gdim() { return gdim_; }

    // operators

    inline LaSymmFactDouble ref(LaSymmFactDouble &);
    inline LaSymmFactDouble ref(LaSymmMatDouble &);
    inline LaSymmFactDouble& copy(const LaSymmFactDouble &);
    inline LaSymmFactDouble& copy(const LaSymmMatDouble &);

};



    // constructor/destructor functions

inline LaSymmFactDouble::LaSymmFactDouble():S_(),pivot_(),info_(0),uplo_('L')
{}


inline LaSymmFactDouble::LaSymmFactDouble(int n, int m):S_(n,m),pivot_(n*m),
                    info_(0),uplo_('L')
{}


inline LaSymmFactDouble::LaSymmFactDouble(const LaSymmFactDouble &F)
{
    S_.copy(F.S_);
    pivot_.copy(F.pivot_);
    info_ = F.info_;
    uplo_ = F.uplo_;
    size_ = F.size_;
    gdim_ = F.gdim_;
}

inline LaSymmFactDouble::~LaSymmFactDouble()
{}

    // operators


inline LaSymmFactDouble LaSymmFactDouble::ref(LaSymmFactDouble& F)
{
    S_.ref(F.S_);
    pivot_.ref(F.pivot_);
    info_ = F.info_;
    uplo_ = F.uplo_;
    size_ = F.size_;
    gdim_ = F.gdim_;
    
    return *this;
}

inline LaSymmFactDouble& LaSymmFactDouble::copy(const LaSymmFactDouble& F)
{
    S_.copy(F.S_);
    pivot_.copy(F.pivot_);
    info_ = F.info_;
    uplo_ = F.uplo_;
    size_ = F.size_;
    gdim_ = F.gdim_;
    
    return *this;
}

inline LaSymmFactDouble LaSymmFactDouble::ref(LaSymmMatDouble &G)
{
    S_.ref(G);
    info_ = 0;
    uplo_ = 'L';
    size_ = G.size(0);
    gdim_ = G.gdim(0);

    return *this;
}

inline LaSymmFactDouble& LaSymmFactDouble::copy(const LaSymmMatDouble &G)
{
    S_.copy(G);
    info_ = 0;
    uplo_ = 'L';
    size_ = G.size(0);
    gdim_ = G.gdim(0);

    return *this;
}

#if 0
inline void LaSymmMatFactorize(LaSymmMatDouble &A, LaSymmFactDouble &AF)
{
    char UPLO = 'L';
    integer N = A.size(0), LDA = A.gdim(0), info = 0;
//    integer M = DSYTRF;
//    integer NB = F77NAME(get_nb)(&N,&M);

    integer LWORK = N*NB;
    double *WORK = new double[LWORK];
    LaVectorLongInt piv(N);
    AF.pivot().copy(piv); // make copies of A and pivot information
    AF.copy(A);

    F77NAME(dsytrf)(&UPLO, &N, &(AF.S()(0,0)), &LDA, &(AF.pivot()(0)), WORK,
                    &LWORK, &info);

    delete [] WORK;
}
#endif
inline void LaLinearSolve(LaSymmFactDouble &AF, LaGenMatDouble &X,
                           LaGenMatDouble &B)
{
    char uplo = 'L';
    integer N = AF.size(), nrhs = X.size(1), lda = AF.gdim(),
            ldb = B.size(0), info = 0;

    X.inject(B);
    F77NAME(dsytrs)(&uplo, &N, &nrhs, &(AF.S()(0,0)), &lda,
            &(AF.pivot()(0)), &X(0,0), &ldb, &info);

}

#endif 
// _LA_SYMM_FACT_DOUBLE_H_
