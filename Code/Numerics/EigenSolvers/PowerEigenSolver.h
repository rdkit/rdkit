//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef _RD_POWER_EIGENSOLVER_H
#define _RD_POWER_EIGENSOLVER_H

#include <Numerics/Vector.h>
#include <Numerics/Matrix.h>
#include <Numerics/SymmMatrix.h>

namespace RDNumeric {
  namespace EigenSolvers {
    //! Compute the \c numEig largest eigenvalues and, optionally,  the corresponding
    //! eigenvectors. 
    /*!
      
    \param numEig       the number of eigenvalues we are interested in
    \param mat          symmetric input matrix of dimension N*N
    \param eigenValues  Vector used to return the eigenvalues (size = numEig)
    \param eigenVectors Optional matrix used to return the eigenvectors (size = N*numEig)
    \param seed         Optional values to seed the random value generator used to 
                        initialize the eigen vectors
    \return a boolean indicating whether or not the calculation converged.
    
    <b>Notes:</b>
    - The matrix, \c mat, is changed in this function
    
    <b>Algorithm:</b>
    
    We use the iterative power method, which works like this:

    \verbatim
     u = arbitrary unit vector
     tol = 0.001
     currEigVal = 0.0;
     prevEigVal = -1.0e100
     while (abs(currEigVal - prevEigVal) > tol) :
         v = Au
         prevEigVal = currEigVal
         currEigVal = v[i] // where i is the id os the largest absolute component
         u = c*v
    \endverbatim

      
    */
    bool powerEigenSolver(unsigned int numEig, DoubleSymmMatrix &mat,
                          DoubleVector &eigenValues,
                          DoubleMatrix *eigenVectors=0, 
                          int seed=-1);
    //! \overload
    static bool powerEigenSolver(unsigned int numEig, DoubleSymmMatrix &mat,
                          DoubleVector &eigenValues,
                          DoubleMatrix &eigenVectors, 
                          int seed=-1) {
      return powerEigenSolver(numEig,mat,eigenValues,&eigenVectors,seed);
    }
  };
};

#endif
                           


