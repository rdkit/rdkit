//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.


#ifndef _LA_BAND_FACT_DOUBLE_H
#define _LA_BAND_FACT_DOUBLE_H

#include "lafnames.h"
//#include LA_VECTOR_LONG_INT_H // changed of VC++
#include "lavli.h"
//#include LA_BAND_MAT_DOUBLE_H // changed of VC++
#include "bmd.h"
#include "lapack.h"


class LaBandFactDouble
{
    LaBandMatDouble         B_;
    LaVectorLongInt         pivot_;
    int                  info_;

public:

    // constructor

    inline LaBandFactDouble();
    inline LaBandFactDouble(int,int,int);
    inline LaBandFactDouble(LaBandMatDouble &);
    inline LaBandFactDouble(LaBandFactDouble &);
    inline ~LaBandFactDouble();

    // extraction functions for components

    inline LaBandMatDouble& B();
    inline LaVectorLongInt& pivot();
    inline int& info();

    // operators

    inline LaBandFactDouble& ref(LaBandFactDouble &);
    inline LaBandFactDouble& ref(LaBandMatDouble &);

};



    // constructor/destructor functions

inline LaBandFactDouble::LaBandFactDouble():B_(),pivot_()
{
#ifdef BandFactDouble_DEBUG 
    cout << " called LaBandFactDouble::LaBandFactDouble() " << endl; 
#endif 

    info_ = 0;
}


inline LaBandFactDouble::LaBandFactDouble(int N, int kl, int ku)
    : B_(N,kl,ku),pivot_(kl+ku+1)
{
#ifdef BandFactDouble_DEBUG 
    cout << " called LaBandFactDouble::LaBandFactDouble(int,int,int) " << endl; 
#endif 

    info_ = 0;
}


inline LaBandFactDouble::LaBandFactDouble(LaBandMatDouble &G):pivot_()
{
#ifdef BandFactDouble_DEBUG 
    cout << " called LaBandFactDouble::LaBandFactDouble(LaBandMatDouble &)"
        <<endl; 
#endif 

  B_.ref(G);
  info_ = 0;
}


inline LaBandFactDouble::LaBandFactDouble(LaBandFactDouble &F)
{
#ifdef BandFactDouble_DEBUG 
    cout << " called LaBandFactDouble::LaBandFactDouble(LaBandFactDouble &) " << endl; 
#endif 

  B_.ref(F.B_);
  pivot_.ref(F.pivot_);
  info_ = F.info_;
}

inline LaBandFactDouble::~LaBandFactDouble()
{
}

    // member functions

inline LaBandMatDouble& LaBandFactDouble::B()
{

    return B_;
}

inline LaVectorLongInt& LaBandFactDouble::pivot()
{

    return pivot_;
}

inline int& LaBandFactDouble::info()
{

    return info_;
}

    
    // operators


inline LaBandFactDouble& LaBandFactDouble::ref(LaBandFactDouble& F)
{

    B_.ref(F.B_);
    pivot_.ref(F.pivot_);
    info_ = F.info_;
    
    return *this;
}


inline LaBandFactDouble& LaBandFactDouble::ref(LaBandMatDouble& F)
{
    B_.ref(F);

    return *this;
}


inline void LaBandMatFactorize(LaBandMatDouble &A, LaBandFactDouble &AF)
{
    integer n = A.size(1), m = n, LDA = A.gdim(0);
    integer KL = A.subdiags(), KU = A.superdiags(), info=0;

    F77NAME(dgbtrf)(&m, &n, &KL, &KU, &A(0,0), &LDA, &(AF.pivot()(0)), &info);
}

inline void LaLinearSolve(LaBandFactDouble &AF, LaGenMatDouble &X,
                        LaGenMatDouble &B)
{
    char trans = 'N';
    integer n = AF.B().size(1), lda = AF.B().gdim(0), nrhs = X.size(1), 
            ldb = B.size(0), kl = AF.B().subdiags(), 
            ku = AF.B().superdiags(), info=0;

    X.inject(B);
    F77NAME(dgbtrs)(
        &trans, 
        &n, 
        &kl, &ku, &nrhs, &(AF.B()(-kl,0)), 
            &lda, &(AF.pivot()(0)), &X(0,0), &ldb, &info);
}


#endif
