// $Id$
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include "BoundsMatrix.h"
#include "TriangleSmooth.h"

namespace DistGeom {
  bool triangleSmoothBounds(BoundsMatPtr boundsMat) {
    return triangleSmoothBounds(boundsMat.get());
  }
  bool triangleSmoothBounds(BoundsMatrix *boundsMat) {
    int npt = boundsMat->numRows();
    int i, j, k;
    double Uik, Lik, Ukj, sumUikUkj, diffLikUjk, diffLjkUik;
    
    for (k = 0; k < npt; k++) {
      for (i = 0; i < npt-1; i++) {
        if (i == k) {
          continue;
        }
	Uik = boundsMat->getUpperBound(i,k);
	Lik = boundsMat->getLowerBound(i,k);
        for (j = i+1; j < npt; j++) {
          if (j == k) {
            continue;
          }
	  Ukj = boundsMat->getUpperBound(k,j);
	  sumUikUkj = Uik + Ukj;
          if (boundsMat->getUpperBound(i,j) > sumUikUkj) {
            boundsMat->setUpperBound(i,j, sumUikUkj);
          } 
          
	  diffLikUjk = Lik - Ukj;
	  diffLjkUik = boundsMat->getLowerBound(j,k) - Uik;
	  if (boundsMat->getLowerBound(i,j) < diffLikUjk) {
            boundsMat->setLowerBound(i,j, diffLikUjk);
          } else if (boundsMat->getLowerBound(i,j) < diffLjkUik) {
            boundsMat->setLowerBound(i,j, diffLjkUik);
          }
           
          if (boundsMat->getLowerBound(i,j) > boundsMat->getUpperBound(i,j)) {
            //std::cout << boundsMat->getLowerBound(i,j) << " " << boundsMat->getUpperBound(i,j) << "\n";
            return false;
          }
        }
      }
    }
    return true;
  }
}
 
