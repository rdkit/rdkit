//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.

#ifndef _LA_SPD_FACT_DOUBLE_H_
#define _LA_SPD_FACT_DOUBLE_H_

//#include LA_SPD_MAT_DOUBLE_H // changed of VC++
#include "spdmd.h"

#include "lapack.h"

class LaSpdFactDouble
{
    int size_;
    int gdim_;
    LaSpdMatDouble S_;

    public:

    LaSpdFactDouble();
    LaSpdFactDouble(int,int);
    LaSpdFactDouble(const LaSpdFactDouble&);
    ~LaSpdFactDouble();

    LaSpdFactDouble& ref(LaSpdFactDouble&);
    LaSpdFactDouble& ref(LaSpdMatDouble&);
    LaSpdFactDouble& copy(const LaSpdFactDouble&);
    LaSpdFactDouble& copy(const LaSpdMatDouble&);

    LaSpdMatDouble& S() { return S_; }
    int size() { return size_; }
    int gdim() { return gdim_; }
};


inline LaSpdFactDouble::LaSpdFactDouble(): S_(), size_(0), gdim_(0)
{}

inline LaSpdFactDouble::LaSpdFactDouble(int m,int n):S_(m,n), 
                size_(n), gdim_(m)
{}

inline LaSpdFactDouble::LaSpdFactDouble(const LaSpdFactDouble &X)
{
    size_ = X.size_;
    gdim_ = X.gdim_;
    S_.copy(X.S_);
}

inline LaSpdFactDouble::~LaSpdFactDouble()
{}

inline LaSpdFactDouble& LaSpdFactDouble::ref(LaSpdFactDouble &X)
{
    size_ = X.size_;
    gdim_ = X.gdim_;
    S_.ref(X.S_);

    return *this;
}

inline LaSpdFactDouble& LaSpdFactDouble::ref(LaSpdMatDouble &X)
{
    size_ = X.size(1);
    gdim_ = X.gdim(0);
    S_.ref(X);

    return *this;
}
    
inline LaSpdFactDouble& LaSpdFactDouble::copy(const LaSpdFactDouble &X)
{
    size_ = X.size_;
    gdim_ = X.gdim_;
    S_.copy(X.S_);

    return *this;
}

inline LaSpdFactDouble& LaSpdFactDouble::copy(const LaSpdMatDouble &X)
{
    size_ = X.size(1);
    gdim_ = X.gdim(0);
    S_.copy(X);

    return *this;
}

inline void LaSpdMatFactorize(LaSpdMatDouble &A, LaSpdFactDouble &AF)
{
    char uplo = 'L';
    integer N = A.size(1), lda = A.gdim(0), info = 0;
    AF.copy(A);

    F77NAME(dpotrf)(&uplo, &N, &(AF.S()(0,0)), &lda, &info);
}


inline void LaLinearSolve(LaSpdFactDouble &AF, LaGenMatDouble &X,
                            LaGenMatDouble &B)
{
    char uplo = 'L';
    integer N = AF.size(), nrhs = X.size(1), lda = AF.gdim(), 
            ldb = B.size(0), info = 0;

    X.inject(B);
    F77NAME(dpotrs)(&uplo, &N, &nrhs, &(AF.S()(0,0)), &lda, &X(0,0),
                    &ldb, &info);
}

#endif 
// _LA_SPD_FACT_DOUBLE_H_
