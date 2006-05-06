//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.


#ifndef _BLAS2_PP_H_
#define _BLAS2_PP_H_


#include "blas2.h"

#ifdef _LA_VECTOR_DOUBLE_H_

#ifdef _LA_GEN_MAT_DOUBLE_H_
void Blas_Mat_Trans_Vec_Mult(LaGenMatDouble &A, LaVectorDouble &dx, 
            LaVectorDouble &dy, double alpha = 1.0, double beta = 1.0);
void Blas_Mat_Vec_Mult(LaGenMatDouble &A, LaVectorDouble &dx, 
            LaVectorDouble &dy, double alpha = 1.0, double beta = 1.0);
void Blas_R1_Update(LaGenMatDouble &A, LaVectorDouble &dx, 
                LaVectorDouble &dy, double alpha = 1.0);
#endif


#ifdef _LA_SYMM_MAT_DOUBLE_H_
void Blas_Mat_Vec_Mult(LaSymmMatDouble &A, LaVectorDouble &dx, 
            LaVectorDouble &dy, double alpha = 1.0, double beta = 1.0);
void Blas_R1_Update(LaSymmMatDouble &A, LaVectorDouble &dx,
                double alpha = 1.0);
void Blas_R2_Update(LaSymmMatDouble &A, LaVectorDouble &dx, 
                LaVectorDouble &dy, double alpha = 1.0);
#endif

#ifdef _LA_SYMM_BAND_MAT_DOUBLE_H_
void Blas_Mat_Vec_Mult(LaSymmBandMatDouble &A, LaVectorDouble &dx, 
            LaVectorDouble &dy, double alpha = 1.0, double beta = 1.0);
#endif

#ifdef _LA_SPD_MAT_DOUBLE_H_ 
void Blas_Mat_Vec_Mult(LaSpdMatDouble &AP, LaVectorDouble &dx, 
        LaVectorDouble &dy, double alpha = 1.0, double beta = 1.0);
void Blas_R1_Update(LaSpdMatDouble &AP, LaVectorDouble &dx,
                double alpha = 1.0);
void Blas_R2_Update(LaSpdMatDouble &AP, LaVectorDouble &dx, 
                LaVectorDouble &dy, double alpha = 1.0);
#endif

#ifdef _LA_BAND_MAT_DOUBLE_H_
void Blas_Mat_Vec_Mult(LaBandMatDouble &A, LaVectorDouble &dx, 
            LaVectorDouble &dy, double alpha = 1.0, double beta = 1.0);
#endif

#ifdef _LA_LOWER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Vec_Mult(LaLowerTriangMatDouble &A, LaVectorDouble &dx);
#endif



#ifdef _LA_UPPER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Vec_Mult(LaUpperTriangMatDouble &A, LaVectorDouble &dx);
#endif



#ifdef _LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H_ 
void Blas_Mat_Vec_Mult(LaUnitLowerTriangMatDouble &A, 
                                LaVectorDouble &dx);
#endif

#ifdef _LA_LOWER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Vec_Solve(LaLowerTriangMatDouble &A, LaVectorDouble &dx);
#endif



#ifdef _LA_UNIT_UPPER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Vec_Mult(LaUnitUpperTriangMatDouble &A, 
                                LaVectorDouble &dx);
#endif

#ifdef _LA_UPPER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Vec_Solve(LaUpperTriangMatDouble &A, LaVectorDouble &dx);
#endif


#ifdef _LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Vec_Solve(LaUnitLowerTriangMatDouble &A, 
                                LaVectorDouble &dx);
#endif

#ifdef _LA_UNIT_UPPER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Vec_Solve(LaUnitUpperTriangMatDouble &A, LaVectorDouble &dx);
#endif

#endif 
    // _LA_VECTOR_DOUBLE_H_

#endif 
    //_BLAS2_PP_H_
