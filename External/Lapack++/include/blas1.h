/*      LAPACK++ (V. 1.1)                                       */
/*      (C) 1992-1996 All Rights Reserved.                          */

#ifndef _BLAS1_H_
#define _BLAS1_H_

#include "arch.h"
#include "f2c.h"

extern "C"
{


    double F77NAME(dasum)(const integer *n, const double *dx, const integer *incx);


    void F77NAME(daxpy)(const integer *n, const double *da, const double *dx, 
            const integer *incx, double *dy, const integer *incy);

    void F77NAME(dcopy)(const integer *n, double *dx, const integer *incx, double *dy, 
                        const integer *incy);


    double F77NAME(ddot)(const integer *n, const double *dx, const integer *incx, 
                        const double *dy, const integer *incy);

    double F77NAME(dnrm2)(const integer *n, const double *dx, const integer *incx); 

    void F77NAME(drot)(const integer *n, double *dx, const integer *incx, double *dy, 
                        const integer *incy, const double *c, const double *s);

    void F77NAME(drotg)(double *da, double *db, double *c, double *s);

    void F77NAME(dscal)(const integer *n, double *da, double *dx, const integer *incx);

    void F77NAME(dswap)(const integer *n, double *dx, const integer *incx, double *dy, 
                        const integer *incy);

    integer F77NAME(idamax)(const integer *n, const double *dx, const integer *incx);

#if defined( LA_COMPLEX_SUPPORT)
#include "f2c.h"

    double F77NAME(zdotc)(doublecomplex *c, const integer *n, 
            const doublecomplex *cx, 
            const integer *incx, const doublecomplex *cy, const integer *incy);

    double F77NAME(zdotu)(doublecomplex *c, const integer *n, 
        const doublecomplex *cx, const integer *incx, 
        const doublecomplex *cy, const integer *incy);

    void F77NAME(zaxpy)(const integer *n, const doublecomplex *da, 
            const doublecomplex *dx, 
            const integer *incx, doublecomplex *dy, 
            const integer *incy);

    void F77NAME(zcopy)(const integer *n, doublecomplex *dx, const integer *incx, 
                doublecomplex *dy, const integer *incy);

    double  F77NAME(dzasum)(const integer *n, const doublecomplex *dx, const integer *incx);

    double  F77NAME(dznrm2)(const integer *n, const doublecomplex *dx, const integer *incx); 

    void F77NAME(zdscal)(const integer *n, const double *da, doublecomplex *dx, 
            const integer *incx);

    void F77NAME(zscal)(const integer *n, const doublecomplex *da, doublecomplex *dx, 
            const integer *incx);

    integer F77NAME(izamax)(const integer *n, const doublecomplex *dx, const integer *incx);

    void F77NAME(zswap)(const integer *n, doublecomplex *dx, const integer *incx, 
                doublecomplex *dy, integer *incy);

#endif
}

#endif

