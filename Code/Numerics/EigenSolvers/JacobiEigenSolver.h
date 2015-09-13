/*
 * JacobiEigenSolver.h
 *
 *  Created on: Aug 28, 2014
 *      Author: hahnda6
 */

#ifndef JACOBIEIGENSOLVER_H_
#define JACOBIEIGENSOLVER_H_

#include <Numerics/Vector.h>
#include <Numerics/SquareMatrix.h>
#include <Numerics/SymmMatrix.h>

namespace RDNumeric {
  namespace EigenSolvers {
    //! Obtain the eigen vectors and eigen values
       /*!
         \param quad        symmetric NxN matrix of interest
         \param eigenVals   storage for eigen values std::vector of size N
         \param eigenVecs   storage for eigen vectors NxN Matrix
         \param maxIter     max number of iterations

         <b>Reference:<\b>
         This is essentailly a copy of the jacobi routine taken from the program
         quatfit.c available from the Computational Chemistry Archives.
         http://www.ccl.net/cca/software/SOURCES/C/quaternion-mol-fit/index.shtml
         E-mail jkl@osc.edu for details.
         It was written by:

         David J. Heisterberg
         The Ohio Supercomputer Center
         1224 Kinnear Rd.
         Columbus, OH  43212-1163
         (614)292-6036
         djh@osc.edu    djh@ohstpy.bitnet    ohstpy::djh
         Copyright: Ohio Supercomputer Center, David J. Heisterberg, 1990.
         The program can be copied and distributed freely, provided that
         this copyright in not removed. You may acknowledge the use of the
         program in published material as:
         David J. Heisterberg, 1990, unpublished results.

         Also see page 463 in Numerical Recipes in C (second edition)
       */
    unsigned int jacobiEigenSolver(const DoubleSymmMatrix &mat,
                          DoubleVector &eigenValues,
                          DoubleSquareMatrix &eigenVectors,
                          unsigned int maxIter=50);

  };
};


#endif /* JACOBIEIGENSOLVER_H_ */
