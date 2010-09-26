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
#include "DistanceConstraint.h"
#include <math.h>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>

namespace ForceFields {
  namespace UFF {
      DistanceConstraintContrib::DistanceConstraintContrib(ForceField *owner,
							   unsigned int idx1,unsigned int idx2,
							   double minLen,double maxLen,
							   double forceConst) {
      PRECONDITION(owner,"bad owner");
      RANGE_CHECK(0,idx1,owner->positions().size()-1);
      RANGE_CHECK(0,idx2,owner->positions().size()-1);
      PRECONDITION(maxLen>=minLen,"bad bounds");

      dp_forceField = owner;
      d_end1Idx = idx1;
      d_end2Idx = idx2;
      d_minLen = minLen;
      d_maxLen = maxLen;
      d_forceConstant = forceConst;

    }

    double DistanceConstraintContrib::getEnergy(double *pos) const {
      PRECONDITION(dp_forceField,"no owner");
      PRECONDITION(pos,"bad vector");
      //return 0.0;
      double dist=this->dp_forceField->distance(this->d_end1Idx,this->d_end2Idx,pos);
      double distTerm=0.0;
      if(dist<this->d_minLen){
	distTerm=this->d_minLen-dist;
      } else if(dist>this->d_maxLen) {
	distTerm=dist-this->d_maxLen;
      }
      double res = 0.5*this->d_forceConstant*distTerm*distTerm;
      //std::cerr << "DIST(" << this->d_end1Idx << ","<<this->d_end2Idx << "): " << this->d_maxLen << " " << dist << " E=" << res << std::endl;
      return res;
    }
    void DistanceConstraintContrib::getGrad(double *pos,double *grad) const {
      PRECONDITION(dp_forceField,"no owner");
      PRECONDITION(pos,"bad vector");
      PRECONDITION(grad,"bad vector");

      double dist=this->dp_forceField->distance(this->d_end1Idx,this->d_end2Idx,pos);

      double preFactor = 0.0;
      if(dist<this->d_minLen){
	preFactor=dist-this->d_minLen;
      } else if(dist>this->d_maxLen) {
	preFactor=dist-this->d_maxLen;
      } else {
	return;
      }
      preFactor *= this->d_forceConstant;
    
      double *end1Coords = &(pos[3*this->d_end1Idx]);
      double *end2Coords = &(pos[3*this->d_end2Idx]);
      for(int i=0;i<3;i++){
	double dGrad;
	if(dist>0.0){
	  dGrad=preFactor * (end1Coords[i]-end2Coords[i])/dist;
	} else {
	  // FIX: this likely isn't right
	  dGrad=preFactor * (end1Coords[i]-end2Coords[i]);
	}
	grad[3*this->d_end1Idx+i] += dGrad;
	grad[3*this->d_end2Idx+i] -= dGrad;
      }    
    }
  
  }
}  
