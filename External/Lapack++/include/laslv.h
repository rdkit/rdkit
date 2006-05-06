//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.



#ifdef _LA_GEN_MAT_DOUBLE_H_
void LaLinearSolve( const LaGenMatDouble& A, LaGenMatDouble& X, 
    LaGenMatDouble& B );

void LaLinearSolveIP( LaGenMatDouble& A, LaGenMatDouble& X, LaGenMatDouble& B );


void LaLULinearSolve( const LaGenMatDouble& A, LaGenMatDouble& X, 
    LaGenMatDouble& B );
void LaLULinearSolveIP( LaGenMatDouble& A, LaGenMatDouble& X, 
        LaGenMatDouble& B );

void LaQRLinearSolve( const LaGenMatDouble& A, LaGenMatDouble& X, 
        LaGenMatDouble& B );
void LaQRLinearSolveIP( LaGenMatDouble& A, LaGenMatDouble& X, 
        LaGenMatDouble& B );



#ifdef _LA_SPD_MAT_DOUBLE_H_

void LaLinearSolve( const LaSpdMatDouble& A, LaGenMatDouble& X, 
    LaGenMatDouble& B );
void LaLinearSolveIP( LaSpdMatDouble& A, LaGenMatDouble& X, LaGenMatDouble& B );

void LaCholLinearSolve( const LaSpdMatDouble& A, LaGenMatDouble& X, 
        LaGenMatDouble& B );
void LaCholLinearSolveIP( LaSpdMatDouble& A, LaGenMatDouble& X, 
        LaGenMatDouble& B );
#endif

#ifdef _LA_SYMM_MAT_DOUBLE_H_

void LaLinearSolve( const LaSymmMatDouble& A, LaGenMatDouble& X, 
    LaGenMatDouble& B );
void LaLinearSolveIP( LaSymmMatDouble& A, LaGenMatDouble& X, LaGenMatDouble& B );

void LaCholLinearSolve( const LaSymmMatDouble& A, LaGenMatDouble& X, 
        LaGenMatDouble& B );
void LaCholLinearSolveIP( LaSymmMatDouble& A, LaGenMatDouble& X, 
        LaGenMatDouble& B );

// Eigenvalue problems 

void LaEigSolve(const LaSymmMatDouble &S, LaVectorDouble &eigvals);
void LaEigSolve(const LaSymmMatDouble &S, LaVectorDouble &eigvals, 
    LaGenMatDouble &eigvec);
void LaEigSolveIP(LaSymmMatDouble &S, LaVectorDouble &eigvals);
void LaEigSolveVecIP(LaSymmMatDouble &S, LaVectorDouble &eigvals);
#endif



#endif
