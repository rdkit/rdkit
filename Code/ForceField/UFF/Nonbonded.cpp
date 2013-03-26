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
#include "Nonbonded.h"
#include "Params.h"
#include <math.h>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>

namespace ForceFields {
  namespace UFF {
    namespace Utils {
      double calcNonbondedMinimum(const AtomicParams *at1Params,
                                  const AtomicParams *at2Params){
        return sqrt(at1Params->x1*at2Params->x1);
      }
      double calcNonbondedDepth(const AtomicParams *at1Params,
                                const AtomicParams *at2Params){
        return sqrt(at1Params->D1*at2Params->D1);
      }

    } // end of namespace utils
    
    vdWContrib::vdWContrib(ForceField *owner,
                           unsigned int idx1,unsigned int idx2,
                           const AtomicParams *at1Params,
                           const AtomicParams *at2Params,
                           double threshMultiplier){
      PRECONDITION(owner,"bad owner");
      PRECONDITION(at1Params,"bad params pointer");
      PRECONDITION(at2Params,"bad params pointer");
      RANGE_CHECK(0,idx1,owner->positions().size()-1);
      RANGE_CHECK(0,idx2,owner->positions().size()-1);

      dp_forceField = owner;
      d_at1Idx = idx1;
      d_at2Idx = idx2;

      // UFF uses the geometric mean of the vdW parameters:
      this->d_xij = Utils::calcNonbondedMinimum(at1Params,at2Params);
      this->d_wellDepth = Utils::calcNonbondedDepth(at1Params,at2Params);
      this->d_thresh = threshMultiplier*this->d_xij;

      //std::cerr << "  non-bonded: " << idx1 << "-" << idx2 << " " << this->d_xij << " " << this->d_wellDepth << " " << this->d_thresh << std::endl;
    }

    double vdWContrib::getEnergy(double *pos) const {
      PRECONDITION(dp_forceField,"no owner");
      PRECONDITION(pos,"bad vector");

      double dist=this->dp_forceField->distance(this->d_at1Idx,this->d_at2Idx,pos);
      if(dist>this->d_thresh || dist<=0.0) return 0.0;

      double r=this->d_xij/dist;
      double r6=int_pow<6>(r);
      double r12=r6*r6;
      double res = this->d_wellDepth*(r12 - 2.0*r6);
      //if(d_at1Idx==12 && d_at2Idx==21 ) std::cerr << "     >: " << d_at1Idx << "-" << d_at2Idx << " " << r << " = " << res << std::endl;
      return res;
    }
    void vdWContrib::getGrad(double *pos,double *grad) const {
      PRECONDITION(dp_forceField,"no owner");
      PRECONDITION(pos,"bad vector");
      PRECONDITION(grad,"bad vector");

      double dist=this->dp_forceField->distance(this->d_at1Idx,this->d_at2Idx,pos);
      if(dist>this->d_thresh) return;

      if(dist<=0){
        for(int i=0;i<3;i++){
          // move in an arbitrary direction
          double dGrad=100.0;
          grad[3*this->d_at1Idx+i] += dGrad;
          grad[3*this->d_at2Idx+i] -= dGrad;
        }    
        return;
      }
      
      double r = this->d_xij/dist;
      double r7 = int_pow<7>(r);
      double r13= int_pow<13>(r);
      double preFactor = 12.*this->d_wellDepth/this->d_xij * (r7-r13);
    
      double *at1Coords = &(pos[3*this->d_at1Idx]);
      double *at2Coords = &(pos[3*this->d_at2Idx]);
      for(int i=0;i<3;i++){
        double dGrad=preFactor * (at1Coords[i]-at2Coords[i])/dist;
        grad[3*this->d_at1Idx+i] += dGrad;
        grad[3*this->d_at2Idx+i] -= dGrad;
      }    
    }
  
  }
}
