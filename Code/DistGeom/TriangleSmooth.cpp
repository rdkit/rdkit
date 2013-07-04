// $Id$
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "BoundsMatrix.h"
#include "TriangleSmooth.h"

namespace DistGeom {
  bool triangleSmoothBounds(BoundsMatPtr boundsMat,double tol) {
    return triangleSmoothBounds(boundsMat.get(),tol);
  }
  bool triangleSmoothBounds(BoundsMatrix *boundsMat,double tol) {
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
          double lBound=boundsMat->getLowerBound(i,j);
          double uBound=boundsMat->getUpperBound(i,j);
          if( tol>0. &&
              (lBound-uBound)/lBound>0. &&
              (lBound-uBound)/lBound<tol ){
            boundsMat->setUpperBound(i,j,lBound);
            uBound=lBound;
          }
          if (lBound - uBound>0.) {
            // std::cerr<<std::endl;
            // for(unsigned int ii=0;ii<npt;++ii){
            //   for(unsigned int jj=0;jj<npt;++jj){
            //     std::cerr<<"  "<<std::setprecision(3)<<boundsMat->getVal(ii,jj);
            //   }
            //   std::cerr<<std::endl;
            // }
            // std::cerr<<std::endl;
            // std::cerr<<" Fail: "<<i<<"-"<<j<<": " << boundsMat->getLowerBound(i,j) << " " << boundsMat->getUpperBound(i,j) << "\n";
            return false;
          }
        }
      }
    }
    return true;
  }
}
 
