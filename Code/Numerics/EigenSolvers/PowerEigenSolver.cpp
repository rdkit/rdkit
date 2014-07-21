// $Id$
//
//  Copyright (C) 2004-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "PowerEigenSolver.h"
#include <Numerics/Vector.h>
#include <Numerics/Matrix.h>
#include <Numerics/SymmMatrix.h>
#include <RDGeneral/Invariant.h>
#include <time.h>

#define MAX_ITERATIONS 1000
#define TOLERANCE 0.001
#define HUGE_EIGVAL 1.0e10
#define TINY_EIGVAL 1.0e-10

namespace RDNumeric {
  namespace EigenSolvers {
    bool powerEigenSolver(unsigned int numEig, DoubleSymmMatrix &mat,
                          DoubleVector &eigenValues, DoubleMatrix *eigenVectors, int seed) {
      // first check all the sizes
      unsigned int N = mat.numRows();
      CHECK_INVARIANT(eigenValues.size() >= numEig, "");
      CHECK_INVARIANT(numEig <= N, "");
      if(eigenVectors){
        unsigned int evRows, evCols;
        evRows = eigenVectors->numRows();
        evCols = eigenVectors->numCols();
        CHECK_INVARIANT(evCols >= N, "");
        CHECK_INVARIANT(evRows >= numEig, "");
      }
      
      unsigned int ei;
      double eigVal, prevVal;
      bool converged=false;
      unsigned int i, j, id, iter, evalId;
      
      DoubleVector v(N), z(N);
      if(seed<=0) seed = clock();
      for (ei = 0; ei < numEig; ei++) {
        eigVal = -HUGE_EIGVAL;
        seed += ei;
        v.setToRandom(seed);

        converged = false;
        for (iter = 0; iter < MAX_ITERATIONS; iter++) {
          // z = mat*v
          multiply(mat, v, z);
          prevVal = eigVal;
          evalId = z.largestAbsValId();
          eigVal = z.getVal(evalId);
          
          if (fabs(eigVal) < TINY_EIGVAL) {
            break;
          }
          
          // compute the next estimate for the eigen vector
          v.assign(z);
          v /= eigVal;
          if (fabs(eigVal - prevVal) < TOLERANCE) {
            converged = true;
            break;
          }
        }
        if (!converged) {
          break;
        }
        v.normalize();
        
        // save this is a eigen vector and value
        // directly access the data instead of setVal so that we save time
        double *vdata = v.getData();
        if(eigenVectors){
          id = ei*eigenVectors->numCols();
          double *eigVecData = eigenVectors->getData();
          for (i = 0; i < N; i++) {
            eigVecData[id + i] = vdata[i];
          }
        }        
        eigenValues[ei]=eigVal;

        // now remove this eigen vector space out of the matrix
        double *matData = mat.getData();
        for (i = 0; i < N; i++) {
          id = i*(i+1)/2;
          for (j = 0; j < i+1; j++) {
            matData[id+j] -= (eigVal*vdata[i]*vdata[j]);
          }
        }
      }
      return converged;
    }
  }
}
