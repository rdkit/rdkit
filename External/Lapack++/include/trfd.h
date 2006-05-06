//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.

#ifndef _LA_TRIDIAG_FACT_DOUBLE_H_
#define _LA_TRIDIAG_FACT_DOUBLE_H_

//#include LA_VECTOR_LONG_INT_H // changed of VC++
#include "lavli.h"
//#include LA_TRIDIAG_MAT_DOUBLE_H // changed of VC++
#include "trmd.h"

#include "lapack.h"

class LaTridiagFactDouble
{
    int size_;
    LaTridiagMatDouble T_;
    LaVectorLongInt pivot_;

public:

    // constructor

    LaTridiagFactDouble();
    LaTridiagFactDouble(int);
    LaTridiagFactDouble(LaTridiagFactDouble &);
    ~LaTridiagFactDouble();

    LaTridiagMatDouble& T() { return T_; }
    LaVectorLongInt& pivot() { return pivot_; }
    int size() { return size_; }
    LaVectorDouble diag(int);

    // operators

    LaTridiagFactDouble& ref(LaTridiagMatDouble &);
    LaTridiagFactDouble& ref(LaTridiagFactDouble &);
    LaTridiagFactDouble& copy(const LaTridiagMatDouble &);
    LaTridiagFactDouble& copy(const LaTridiagFactDouble &);

};



    // constructor/destructor functions

inline LaTridiagFactDouble::LaTridiagFactDouble():T_(),pivot_(),size_(0)
{}


inline LaTridiagFactDouble::LaTridiagFactDouble(int N):T_(N),pivot_(N),size_(N)
{}


inline LaTridiagFactDouble::LaTridiagFactDouble(LaTridiagFactDouble &F)
{
  T_.copy(F.T_);
  pivot_.copy(F.pivot_);
  size_ = F.size_;
}

inline LaTridiagFactDouble::~LaTridiagFactDouble()
{}

    // member functions

inline LaVectorDouble LaTridiagFactDouble::diag(int k)
{
    return T_.diag(k);
}
    
    // operators


inline LaTridiagFactDouble& LaTridiagFactDouble::ref(LaTridiagFactDouble& F)
{
    T_.ref(F.T_);
    pivot_.ref(F.pivot_);
    size_ = F.size_;
    
    return *this;
}

inline LaTridiagFactDouble& LaTridiagFactDouble::ref(LaTridiagMatDouble& A)
{
    T_.ref(A);

    return *this;
}

inline LaTridiagFactDouble& LaTridiagFactDouble::copy(const LaTridiagFactDouble& F)
{
    T_.copy(F.T_);
    pivot_.copy(F.pivot_);
    size_ = F.size_;
    
    return *this;
}

inline LaTridiagFactDouble& LaTridiagFactDouble::copy(const LaTridiagMatDouble& A)
{
    T_.copy(A);

    return *this;
}

inline void LaTridiagMatFactorize(LaTridiagMatDouble &A,
                                 LaTridiagFactDouble &AF)
{
    integer N = A.size(), info = 0;
    AF.copy(A);
    double *DL = &AF.diag(-1)(0), *D = &AF.diag(0)(0),
         *DU = &AF.diag(1)(0), *DU2 = &AF.diag(2)(0);

cerr << " \t*\n";

    F77NAME(dgttrf)(&N, DL, D, DU, DU2, &(AF.pivot()(0)), &info);

cerr << " \t\t**\n";
}


inline void LaLinearSolve(LaTridiagFactDouble &AF, LaGenMatDouble &X,
                        LaGenMatDouble &B)
{
    char trans = 'N';
    integer N = AF.size(), nrhs = X.size(1), ldb = B.size(0), info = 0;
    double *DL = &AF.diag(-1)(0), *D = &AF.diag(0)(0),
         *DU =  &AF.diag(1)(0), *DU2 = &AF.diag(2)(0);

    X.inject(B);
    F77NAME(dgttrs)(&trans, &N, &nrhs, DL, D, DU, DU2, &(AF.pivot()(0)),
                    &X(0,0), &ldb, &info);
}

#endif 
// _LA_TRIDIAG_FACT_DOUBLE_H_
