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

#include "BondStretch.h"
#include "Params.h"
#include <math.h>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>

namespace ForceFields {
  namespace UFF {
    namespace Utils {
      double calcBondRestLength(double bondOrder,
				const AtomicParams *end1Params,
				const AtomicParams *end2Params){
	PRECONDITION(bondOrder>0,"bad bond order");
	PRECONDITION(end1Params,"bad params pointer");
	PRECONDITION(end2Params,"bad params pointer");

	double ri=end1Params->r1,rj=end2Params->r1;

	// this is the pauling correction:
	double rBO = -Params::lambda * (ri + rj) * log(bondOrder);

	// O'Keefe and Breese electronegativity correction:
	double Xi=end1Params->GMP_Xi,Xj=end2Params->GMP_Xi;
	double rEN = ri*rj*(sqrt(Xi)-sqrt(Xj))*(sqrt(Xi)-sqrt(Xj)) /
	  (Xi*ri + Xj*rj);
    
	double res = ri + rj + rBO - rEN;
	return res;
      }
  
      double calcBondForceConstant(double restLength,
				   const AtomicParams *end1Params,
				   const AtomicParams *end2Params){
	double res = 2.0 * Params::G * end1Params->Z1 * end2Params->Z1 / (restLength*restLength*restLength);
	return res;
      }
    } // end of namespace Utils
  
    BondStretchContrib::BondStretchContrib(ForceField *owner,
					   unsigned int idx1,unsigned int idx2,
					   double bondOrder,
					   const AtomicParams *end1Params,
					   const AtomicParams *end2Params){
      PRECONDITION(owner,"bad owner");
      PRECONDITION(end1Params,"bad params pointer");
      PRECONDITION(end2Params,"bad params pointer");
      RANGE_CHECK(0,idx1,owner->positions().size()-1);
      RANGE_CHECK(0,idx2,owner->positions().size()-1);

      dp_forceField = owner;
      d_end1Idx = idx1;
      d_end2Idx = idx2;

      d_restLen = Utils::calcBondRestLength(bondOrder,
						  end1Params,end2Params);
      d_forceConstant = Utils::calcBondForceConstant(d_restLen,
							   end1Params,end2Params);
    }

    double BondStretchContrib::getEnergy(double *pos) const {
      PRECONDITION(dp_forceField,"no owner");
      PRECONDITION(pos,"bad vector");

      double distTerm=dp_forceField->distance(d_end1Idx,d_end2Idx,pos) -
	d_restLen;
      double res = 0.5*d_forceConstant*distTerm*distTerm;
      return res;
    }
    void BondStretchContrib::getGrad(double *pos,double *grad) const {
      PRECONDITION(dp_forceField,"no owner");
      PRECONDITION(pos,"bad vector");
      PRECONDITION(grad,"bad vector");


      double dist=dp_forceField->distance(d_end1Idx,d_end2Idx,pos);
      double preFactor = d_forceConstant*(dist-d_restLen);

      //std::cout << "\tDist("<<d_end1Idx<<","<<d_end2Idx<<") " << dist << std::endl;
      double *end1Coords = &(pos[3*d_end1Idx]);
      double *end2Coords = &(pos[3*d_end2Idx]);
      for(int i=0;i<3;i++){
	double dGrad;
	if(dist>0.0){
	  dGrad=preFactor * (end1Coords[i]-end2Coords[i])/dist;
	} else {
	  // move a small amount in an arbitrary direction
	  dGrad=d_forceConstant*.01;
	}
	grad[3*d_end1Idx+i] += dGrad;
	grad[3*d_end2Idx+i] -= dGrad;
      }    
    }
  
  }
}  
