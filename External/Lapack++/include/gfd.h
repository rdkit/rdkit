//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.


#ifndef _LA_GEN_FACT_DOUBLE_H
#define _LA_GEN_FACT_DOUBLE_H

#include "lafnames.h"
//#include LA_VECTOR_LONG_INT_H // changed of VC++
#include "lavli.h"
//#include LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H // changed of VC++
#include "ultgmd.h"
//#include LA_UPPER_TRIANG_MAT_DOUBLE_H // changed of VC++
#include "utgmd.h"

#include "lapack.h"

class LaGenFactDouble
{
    LaUnitLowerTriangMatDouble  L_;
    LaUpperTriangMatDouble      U_;
    LaVectorLongInt             pivot_;
    int                      info_;
    int                 transpose_;

public:

    // constructor

    inline LaGenFactDouble();
    inline LaGenFactDouble(int,int);
    inline LaGenFactDouble(LaGenFactDouble &);
    inline ~LaGenFactDouble();

    // extraction functions for components

    inline LaUnitLowerTriangMatDouble& L();
    inline LaUpperTriangMatDouble& U();
    inline LaVectorLongInt& pivot();
    inline int& info();
    inline int& transpose();

    // operators

    inline LaGenFactDouble& ref(LaGenFactDouble &);
    inline LaGenFactDouble& ref(LaGenMatDouble &);

};



    // constructor/destructor functions

inline LaGenFactDouble::LaGenFactDouble():L_(),U_(),pivot_()
{

    info_ = 0;
    transpose_ = 0;
}


inline LaGenFactDouble::LaGenFactDouble(int n, int m):L_(n,m),U_(n,m),pivot_(n*m)
{

    info_ = 0;
    transpose_ = 0;
}


inline LaGenFactDouble::LaGenFactDouble(LaGenFactDouble &F)
{

  L_.ref(F.L_);
  U_.ref(F.U_);
  pivot_.ref(F.pivot_);
  info_ = F.info_;
  transpose_ = F.transpose_;
}

inline LaGenFactDouble::~LaGenFactDouble()
{
}

    // member functions

inline LaUnitLowerTriangMatDouble& LaGenFactDouble::L()
{

    return L_;
}

inline LaUpperTriangMatDouble& LaGenFactDouble::U()
{

    return U_;
}

inline LaVectorLongInt& LaGenFactDouble::pivot()
{

    return pivot_;
}

inline int& LaGenFactDouble::info()
{

    return info_;
}

inline int& LaGenFactDouble::transpose()
{

    return transpose_;
}

    
    // operators


inline LaGenFactDouble& LaGenFactDouble::ref(LaGenFactDouble& F)
{

    L_.ref(F.L_);
    U_.ref(F.U_);
    pivot_.ref(F.pivot_);
    info_ = F.info_;
    transpose_ = F.transpose_;
    
    return *this;
}

inline LaGenFactDouble& LaGenFactDouble::ref(LaGenMatDouble &G)
{

  L_.ref(G);
  U_.ref(G);
  info_ = 0;
  transpose_ = 0;

  return *this;
}

#if 0
inline void LaLinearSolve(LaGenFactDouble &AF, LaGenMatDouble& X,
    LaGenMatDouble& B )
{
    char trans = 'N';
    integer n = AF.L().size(1), lda = AF.L().gdim(0), nrhs = X.size(1),
            ldb = B.size(0), info = 0;

    X.inject(B);
    F77NAME(dgetrs)(&trans, &n, &nrhs, &(AF.U()(0,0)), &lda, &(AF.pivot()(0)),
         &X(0,0), &ldb, &info);
}

inline void LaGenMatFactorize(LaGenMatDouble &GM, LaGenFactDouble &GF)
{
    integer m = GM.size(0), n = GM.size(1), lda = GM.gdim(0);
    integer info=0;

    F77NAME(dgetrf)(&m, &n, &GM(0,0), &lda, &(GF.pivot()(0)), &info);
}

inline void LaGenMatFactorizeUnblocked(LaGenMatDouble &A, LaGenFactDouble &F)
{
    integer m = A.size(0), n=A.size(1), lda = A.gdim(0);
    integer info=0;

    F77NAME(dgetf2)(&m, &n, &A(0,0), &lda, &(F.pivot()(0)), &info);
}
#endif

void LaLUFactorDouble(LaGenMatDouble &A, LaGenFactDouble &F, integer nb);
void LaLUFactorDouble(LaGenMatDouble &A, LaGenFactDouble &F);

#endif
