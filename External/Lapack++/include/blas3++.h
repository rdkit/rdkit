//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.


#ifndef _BLAS3_PP_H_
#define _BLAS3_PP_H_

#include "blas3.h"


#ifdef _LA_GEN_MAT_DOUBLE_H_
void Blas_Mat_Mat_Mult(const LaGenMatDouble &A, 
            const LaGenMatDouble &B, LaGenMatDouble &C, 
            double alpha = 1.0, double beta = 0.0);

void Blas_Mat_Trans_Mat_Mult(const LaGenMatDouble &A, 
            const LaGenMatDouble &B, LaGenMatDouble &C, 
            double alpha = 1.0, double beta = 0.0);

void Blas_Mat_Mat_Trans_Mult(const LaGenMatDouble &A, 
            const LaGenMatDouble &B, LaGenMatDouble &C, 
            double alpha = 1.0, double beta = 0.0);



#ifdef _LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H_

void Blas_Mat_Mat_Solve(LaUnitLowerTriangMatDouble &A, 
            LaGenMatDouble &B, double alpha = 1.0);

#endif

#ifdef _LA_UNIT_UPPER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Mat_Mult(LaUnitUpperTriangMatDouble &A,
            LaGenMatDouble &B, double alpha = 1.0);

void Blas_Mat_Mat_Solve(LaUnitUpperTriangMatDouble &A, 
            LaGenMatDouble &B, double alpha = 1.0);

#endif

#ifdef _LA_LOWER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Mat_Mult(LaLowerTriangMatDouble &A,
            LaGenMatDouble &B, double alpha = 1.0);

void Blas_Mat_Mat_Solve(LaLowerTriangMatDouble &A, 
            LaGenMatDouble &B, double alpha = 1.0);
#endif


#ifdef _LA_UPPER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Mat_Mult(LaUpperTriangMatDouble &A,
            LaGenMatDouble &B, double alpha = 1.0);

void Blas_Mat_Mat_Solve(LaUpperTriangMatDouble &A, 
            LaGenMatDouble &B, double alpha = 1.0);
#endif


#ifdef _LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Trans_Mat_Solve(LaUnitLowerTriangMatDouble &A,
            LaGenMatDouble &B, double alpha = 1.0);
#endif

#ifdef _LA_UNIT_UPPER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Trans_Mat_Solve(LaUnitUpperTriangMatDouble &A,
            LaGenMatDouble &B, double alpha = 1.0);
#endif

#ifdef _LA_UNIT_LOWER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Mat_Mult(LaUnitLowerTriangMatDouble &A,
            LaGenMatDouble &B, double alpha = 1.0);

#endif

#ifdef _LA_LOWER_TRIANG_MAT_DOUBLE_H_
void Blas_Mat_Trans_Mat_Solve(LaLowerTriangMatDouble &A,
            LaGenMatDouble &B, double alpha = 1.0);
#endif


#ifdef _LA_UPPER_TRIANG_MAT_DOUBLE_H_ 
void Blas_Mat_Trans_Mat_Solve(LaUpperTriangMatDouble &A,
            LaGenMatDouble &B, double alpha = 1.0);

#endif

#ifdef _LA_SYMM_MAT_DOUBLE_H_
void Blas_Mat_Mat_Mult(LaSymmMatDouble &A, LaGenMatDouble &B, 
            LaGenMatDouble &C, double alpha = 1.0, double beta = 1.0);

void Blas_R1_Update(LaSymmMatDouble &C, LaGenMatDouble &A,
            double alpha = 1.0, double beta = 1.0);
void Blas_R2_Update(LaSymmMatDouble &C, LaGenMatDouble &A,
            LaGenMatDouble &B, double alpha = 1.0, double beta = 1.0);
#endif



#endif
    /* _LA_GEN_MAT_DOUBLE_H_ */


#endif 
    // _BLAS3_PP_H_
            
