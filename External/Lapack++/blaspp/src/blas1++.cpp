//
//              LAPACK++ 1.1a Linear Algebra Package 1.1a
//               University of Tennessee, Knoxvilee, TN.
//            Oak Ridge National Laboratory, Oak Ridge, TN.
//        Authors: J. J. Dongarra, E. Greaser, R. Pozo, D. Walker
//                 (C) 1992-1996 All Rights Reserved
//
//                             NOTICE
//
// Permission to use, copy, modify, and distribute this software and
// its documentation for any purpose and without fee is hereby granted
// provided that the above copyright notice appear in all copies and
// that both the copyright notice and this permission notice appear in
// supporting documentation.
//
// Neither the Institutions (University of Tennessee, and Oak Ridge National
// Laboratory) nor the Authors make any representations about the suitability 
// of this software for any purpose.  This software is provided ``as is'' 
// without express or implied warranty.
//
// LAPACK++ was funded in part by the U.S. Department of Energy, the
// National Science Foundation and the State of Tennessee.

#include "lafnames.h"
//#include LA_VECTOR_DOUBLE_H // changed of VC++
#include "lavd.h"
#include "blas1++.h"

     
double Blas_Norm1(const LaVectorDouble &dx)
{
    integer n = dx.size();
    integer incx = dx.inc();

    return F77NAME(dasum)(&n, &dx(0), &incx);
}


void Blas_Add_Mult(LaVectorDouble &dy, double da, const LaVectorDouble &dx) 
{
    assert(dx.size()==dy.size());
    integer n = dx.size();
    integer incx = dx.inc(), incy = dy.inc();

    F77NAME(daxpy)(&n, &da, &dx(0), &incx, &dy(0), &incy);
}

void Blas_Mult(LaVectorDouble &dy, double da, LaVectorDouble &dx)
{
    assert(dx.size()==dy.size());
    integer n = dx.size();
    integer incx = dx.inc(), incy = dy.inc();
    dy = 0.0;

    F77NAME(daxpy)(&n, &da, &dx(0), &incx, &dy(0), &incy);
}


void Blas_Copy(LaVectorDouble &dy, LaVectorDouble &dx)
{
    assert(dx.size()==dy.size());
    integer n = dx.size();
    integer incx = dx.inc(), incy = dy.inc();

    F77NAME(dcopy)(&n, &dx(0), &incx, &dy(0), &incy);
}


double Blas_Dot_Prod(const LaVectorDouble &dx, const LaVectorDouble &dy)
{
    assert(dx.size()==dy.size());
    integer n = dx.size();
    integer incx = dx.inc(), incy = dy.inc();

    return F77NAME(ddot)(&n, &dx(0), &incx, &dy(0), &incy);
}


double Blas_Norm2(const LaVectorDouble &dx)
{
    integer n = dx.size();
    integer incx = dx.inc();

    return F77NAME(dnrm2)(&n, &dx(0), &incx);
}



void Blas_Apply_Plane_Rot(LaVectorDouble &dx, LaVectorDouble &dy, 
                double &c, double &s)
{
    assert(dx.size()==dy.size());
    integer n = dx.size();
    integer incx = dx.inc(), incy = dy.inc();

    F77NAME(drot)(&n, &dx(0), &incx, &dy(0), &incy, &c, &s);
}



void Blas_Gen_Plane_Rot(double &da, double &db, double &c, double &s)
{
    F77NAME(drotg)(&da, &db, &c, &s);
}

// Note: AT&T C++ 4.1.3 causes this to crash.  Works fine with other
// compilers.
//
//#if 0
void Blas_Scale(double da, LaVectorDouble &dx)
{
    integer n = dx.size();
    integer incx = dx.inc();


    F77NAME(dscal)(&n, &da, &dx(0), &incx);
}
// #endif


void Blas_Swap(LaVectorDouble &dx, LaVectorDouble &dy)
{
    assert(dx.size()==dy.size());
    integer n = dx.size();
    integer incx = dx.inc(), incy = dy.inc();

    F77NAME(dswap)(&n, &dx(0), &incx, &dy(0), &incy);
}



int Blas_Index_Max(const LaVectorDouble &dx)
{
    integer n = dx.size();
    integer incx = dx.inc();

    // subtract one from index since f77 starts at 1, not 0.
    return F77NAME(idamax)(&n, &dx(0), &incx) - 1;
}

// Complex Routines
#ifdef LA_COMPLEX_SUPPORT
//#include LA_VECTOR_COMPLEX_H // changed of VC++
#include "lavc.h"


COMPLEX Blas_H_Dot_Prod(const LaVectorComplex &cx, const LaVectorComplex &cy)
{
    assert(cx.size()==cy.size());
    integer n = cx.size();
    integer incx = cx.inc(), incy = cy.inc();
    COMPLEX tmp;

    F77NAME(zdotc)(&tmp, &n, &cx(0), &incx, &cy(0), &incy);
    return tmp;
}

COMPLEX Blas_U_Dot_Prod(const LaVectorComplex &cx, const LaVectorComplex &cy)
{
    assert(cx.size()==cy.size());
    integer n = cx.size();
    integer incx = cx.inc(), incy = cy.inc();
    COMPLEX tmp;

    F77NAME(zdotu)(&tmp, &n, &cx(0), &incx, &cy(0), &incy);
    return tmp;
}


double Blas_Norm1(const LaVectorComplex &dx)
{
    integer n = dx.size();
    integer incx = dx.inc();

    return F77NAME(dzasum)(&n, &dx(0), &incx);
}


void Blas_Add_Mult(LaVectorComplex &dy, COMPLEX da, const LaVectorComplex &dx) 
{
    assert(dx.size()==dy.size());
    integer n = dx.size();
    integer incx = dx.inc(), incy = dy.inc();

    F77NAME(zaxpy)(&n, &da, &dx(0), &incx, &dy(0), &incy);
}

void Blas_Mult(LaVectorComplex &dy, COMPLEX da, LaVectorComplex &dx)
{
    assert(dx.size()==dy.size());
    integer n = dx.size();
    integer incx = dx.inc(), incy = dy.inc();
    dy = COMPLEX(0.0,0.0);

    F77NAME(zaxpy)(&n, &da, &dx(0), &incx, &dy(0), &incy);
}


void Blas_Copy(LaVectorComplex &dy, LaVectorComplex &dx)
{
    assert(dx.size()==dy.size());
    integer n = dx.size();
    integer incx = dx.inc(), incy = dy.inc();

    F77NAME(zcopy)(&n, &dx(0), &incx, &dy(0), &incy);
}




double Blas_Norm2(const LaVectorComplex &dx)
{
    integer n = dx.size();
    integer incx = dx.inc();

    return F77NAME(dznrm2)(&n, &dx(0), &incx);
}



void Blas_Scale(COMPLEX za, LaVectorComplex &dx)
{
    integer n = dx.size();
    integer incx = dx.inc();

    F77NAME(zscal)(&n, &za, &dx(0), &incx);
}


void Blas_Swap(LaVectorComplex &dx, LaVectorComplex &dy)
{
    assert(dx.size()==dy.size());
    integer n = dx.size();
    integer incx = dx.inc(), incy = dy.inc();

    F77NAME(zswap)(&n, &dx(0), &incx, &dy(0), &incy);
}



int Blas_Index_Max(const LaVectorComplex &dx)
{
    integer n = dx.size();
    integer incx = dx.inc();

    // subtract one from index since f77 starts at 1, not 0.
    return F77NAME(izamax)(&n, &dx(0), &incx) - 1;
}

#endif
// LA_COMPLEX_SUPPORT
