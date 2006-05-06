// $Id$
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include "BoundsMatrix.h"
#include "DistGeomUtils.h"
#include "DistViolationContrib.h"
#include <Numerics/Matrix.h>
#include <Numerics/SymmMatrix.h>
#include <Numerics/Vector.h>
#include <RDGeneral/Invariant.h>
#include <Numerics/EigenSolvers/PowerEigenSolver.h>
#include <RDGeneral/utils.h>
#include <ForceField/ForceField.h>

namespace DistGeom {
  void pickRandomDistMat(const BoundsMatrix &mmat, 
                         RDNumeric::SymmMatrix<double> &distMat, int seed) {
    RDKit::rng_type &generator = RDKit::getRandomGenerator();
    if (seed > 0) {
      //std::cerr << "   SEED: " << seed << std::endl;
      generator.seed(seed);
    }
    //std::cerr << " sample: " << RDKit::getRandomVal() << std::endl;

    // make sure the sizes match up
    unsigned int npt = mmat.numRows();
    CHECK_INVARIANT(npt == distMat.numRows(), "Size mismatch");
    unsigned int i, j, id;
    double *ddata = distMat.getData();
    double ub, lb, d;
    for (i = 1; i < npt; i++) {
      id = i*(i+1)/2;
      for (j = 0; j < i; j++) {
        ub = mmat.getUpperBound(i,j);
        lb = mmat.getLowerBound(i,j);
        CHECK_INVARIANT(ub >= lb, "");
        double rval = RDKit::getRandomVal();
        //std::cout << rval << "\n";
        d = lb + (rval)*(ub - lb);
        ddata[id+j] = d;
      }
    }
  }
  /*
  void _eigenFromLapack(const RDNumeric::SymmMatrix<double> &T, 
                        RDNumeric::DoubleVector &eigVals, RDNumeric::DoubleMatrix &eigVecs) {
    unsigned int N = T.size();
    LaSymmMatDouble A(N, N);
    LaGenMatDouble eVecs(N,N);
    LaVectorDouble eVals(N);
    unsigned int i,j;
    // this copying is unfortunately necessary because our SymmMatrix representation is
    // diff from that of lapack
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        A(i,j) = A(j,i) = T.getVal(i,j);
      }
    }
    
    LaEigSolveIP(A, eVals); //, eVecs);
    // copy them out 
    for (i = 0; i < 3; i++) {
      eigVals.setVal(i, eVals(i));
      for (j = 0; j < N; j++) {
        eigVecs.setVal(i, j, eVecs(i,j));
      }
    }
    }*/

  bool computeInitialCoords(const RDNumeric::SymmMatrix<double> &distMat,  
                            PointPtrVect &positions, bool randNegEig, 
                            unsigned int numZeroFail) {
    unsigned int N = distMat.numRows();
    unsigned int nPt = positions.size();
    CHECK_INVARIANT(nPt == N, "Size mismatch");
    
    const double *data = distMat.getData();
    RDNumeric::SymmMatrix<double> sqMat(N), T(N, 0.0); 
    RDNumeric::DoubleMatrix eigVecs(3,N);
    RDNumeric::DoubleVector eigVals(3);
    
    unsigned int i, j;
    double *sqDat = sqMat.getData();
    
    unsigned int dSize = distMat.getDataSize();
    double sumSqD2  = 0.0;      
    for (i = 0; i < dSize; i++) {
      sqDat[i] = data[i]*data[i];
      sumSqD2 += sqDat[i];
    }
    sumSqD2 /= (N*N);

    RDNumeric::DoubleVector sqD0i(N, 0.0);
    double *sqD0iData = sqD0i.getData();
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        sqD0iData[i] += sqMat.getVal(i,j);
      }
      sqD0iData[i] /= N;
      sqD0iData[i] -= sumSqD2;
      
      if ((sqD0iData[i] < EIGVAL_TOL) && (N > 3)){
        return false;
      }
            
    }

    double val;
    for (i = 0; i < N; i++) {
      for (j = 0; j <= i; j++) {
        val = 0.5*(sqD0iData[i] + sqD0iData[j] - sqMat.getVal(i,j));
        T.setVal(i,j, val);
      }
    }
    int nEigs = (3 < N) ? 3 : N;
    RDNumeric::EigenSolvers::powerEigenSolver(nEigs, T, eigVecs, eigVals, (int)(sumSqD2*N));
    
    double *eigData = eigVals.getData();
    bool foundNeg = false;
    unsigned int zeroEigs = 0;
    for (i = 0; i < 3; i++) {
      if (eigData[i] > EIGVAL_TOL) {
        eigData[i] = sqrt(eigData[i]);
      } else if (fabs(eigData[i]) < EIGVAL_TOL) {
        eigData[i] = 0.0;
        zeroEigs++;
      } else {
        foundNeg = true;
      }
    }
    if ((foundNeg) && (!randNegEig) ) {
      return false;
    }

    if ((zeroEigs >= numZeroFail) && (N > 3)) {
      // this probably happened because we have degenerate eigen values 
      // and the powereigen solver doesn't do good job as a result
      // Use lapack in this case
      /*
      _eigenFromLapack(T, eigVals, eigVecs);
      zeroEigs = 0;
      for (i = 0; i < 3; i++) {
        if (eigData[i] > EIGVAL_TOL) {
          eigData[i] = sqrt(eigData[i]);
        } else if (fabs(eigData[i]) < EIGVAL_TOL) {
          eigData[i] = 0.0;
          zeroEigs++;
        } else {
          foundNeg = true;
        }
      }
      if (zeroEigs > 1) {
        return false;
        }*/
      return false;
    }

    for (i = 0; i < N; i++) {
      RDGeom::Point *pt = positions[i];
      if (eigData[0] >= 0.0) {
        pt->x = eigData[0]*eigVecs.getVal(0,i);
      } else {
        pt->x = 1.0 - 2.0*RDKit::getRandomVal();
      }
      if (eigData[1] >= 0.0) {
        pt->y = eigData[1]*eigVecs.getVal(1,i);
      } else {
        pt->y = 1.0 - 2.0*RDKit::getRandomVal();
      }
      if (eigData[2] >= 0.0) {
        pt->z = eigData[2]*eigVecs.getVal(2,i);
      } else {
        pt->z = 1.0 - 2.0*RDKit::getRandomVal();
      }
    }
    return true;
  }

  ForceFields::ForceField *constructForceField(const BoundsMatrix &mmat,
					       PointPtrVect &positions,
					       std::map< std::pair<int,int>,double> *extraWeights,
					       double basinSizeTol) {
    unsigned int N = mmat.numRows();
    CHECK_INVARIANT(N == positions.size(), "");
    ForceFields::ForceField *field=new ForceFields::ForceField();
    unsigned int i, j;
    for(i=0; i < N; i++){
      field->positions().push_back(positions[i]);
    }
    
    for (i = 1; i < N; i++) {
      for (j = 0; j < i; j++) {
	double w = 1.0;
        double l = mmat.getLowerBound(i,j);
        double u = mmat.getUpperBound(i,j);
	bool includeIt=false;
	if(extraWeights){
	  std::map< std::pair<int,int>,double>::const_iterator mapIt;
	  mapIt = extraWeights->find(std::make_pair<int,int>(i,j));
	  if(mapIt != extraWeights->end()){
	    w = mapIt->second;
	    includeIt=true;
	  }
	}
	if(u-l <= basinSizeTol) {
	  includeIt=true;
	}
	if(includeIt){
	  DistViolationContrib *contrib = new DistViolationContrib(field, i, j, u, l, w);
	  field->contribs().push_back(ForceFields::ContribPtr(contrib));
	}
      }
    }
    return field;
  }
    
}
