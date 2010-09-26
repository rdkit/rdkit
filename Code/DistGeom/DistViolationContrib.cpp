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

#include "DistViolationContrib.h"
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>

namespace DistGeom {
  
  DistViolationContrib::DistViolationContrib(ForceFields::ForceField *owner,unsigned int idx1,
                                             unsigned int idx2,
                                             double ub, double lb, double weight) {
    PRECONDITION(owner,"bad owner");
    RANGE_CHECK(0,idx1,owner->positions().size()-1);
    RANGE_CHECK(0,idx2,owner->positions().size()-1);

    dp_forceField = owner;
    d_end1Idx = idx1;
    d_end2Idx = idx2;
    d_ub = ub;
    d_lb = lb;
    d_weight = weight;
  }

  double DistViolationContrib::getEnergy(double *pos) const {
    PRECONDITION(dp_forceField,"no owner");
    PRECONDITION(pos,"bad vector");
    
    double d = this->dp_forceField->distance(this->d_end1Idx,this->d_end2Idx,pos);
    double val=0.0;
    if (d > d_ub) {
      val = ((d*d)/(d_ub*d_ub)) - 1.0;
    } else if (d < d_lb) {
      val = ((2*d_lb*d_lb)/(d_lb*d_lb + d*d)) - 1.0;
    }
    double res;
    if(val>0.0){
      res=d_weight*val*val;
    }else{
      res=0;
    }
    return res;
  }

  void DistViolationContrib::getGrad(double *pos, double *grad) const {
    PRECONDITION(dp_forceField,"no owner");
    PRECONDITION(pos,"bad vector");
    PRECONDITION(grad,"bad vector");
    unsigned int dim=this->dp_forceField->dimension();
    double d = this->dp_forceField->distance(this->d_end1Idx,this->d_end2Idx,pos);
    double preFactor = 0.0;
    if (d > d_ub) {
      double u2 = d_ub*d_ub;
      preFactor = 4.*(((d*d)/u2) - 1.0)*(d/u2);
    } else if (d < d_lb) {
      double d2 = d*d;
      double l2 = d_lb*d_lb;
      double l2d2 = d2+l2;
      preFactor = 8.*l2*d * (1.-2*l2/l2d2) / (l2d2*l2d2);
      //preFactor = -8.*((l2-d2)/pow(l2+d2,3))*d*l2;
    } else {
      return;
    }

    double *end1Coords = &(pos[dim*this->d_end1Idx]);
    double *end2Coords = &(pos[dim*this->d_end2Idx]);

    for(unsigned int i=0;i<dim;i++){
      double dGrad;
      if(d>0.0){
        dGrad= d_weight*preFactor * (end1Coords[i]-end2Coords[i])/d;
      } else {
        // FIX: this likely isn't right
        dGrad= d_weight*preFactor * (end1Coords[i]-end2Coords[i]);
      }
      grad[dim*this->d_end1Idx+i] += dGrad;
      grad[dim*this->d_end2Idx+i] -= dGrad;
    }
  }
}
