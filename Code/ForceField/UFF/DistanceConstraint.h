//
//  Copyright (C) 2013 Paolo Tosco
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_UFFDISTANCECONSTRAINT_H__
#define __RD_UFFDISTANCECONSTRAINT_H__
#include <iostream>
#include <ForceField/Contrib.h>

namespace ForceFields {
  namespace UFF {

    //! A distance range constraint modelled after a BondStretchContrib
    class DistanceConstraintContrib : public ForceFieldContrib {
    public:
      DistanceConstraintContrib() : d_end1Idx(-1), d_end2Idx(-1) {};
      //! Constructor
      /*!
	\param owner       pointer to the owning ForceField
	\param idx1        index of end1 in the ForceField's positions
	\param idx2        index of end2 in the ForceField's positions
	\param minLen      minimum distance
	\param maxLen      maximum distance
	\param forceConst  force Constant
	
      */
      DistanceConstraintContrib(ForceField *owner, unsigned int idx1, unsigned int idx2,
				double minLen, double maxLen, double forceConst);
      DistanceConstraintContrib(ForceField *owner, unsigned int idx1, unsigned int idx2,
				bool relative, double minLen, double maxLen, double forceConst);

      ~DistanceConstraintContrib() {
	//std::cerr << " ==== Destroy constraint " << d_end1Idx << " " << d_end2Idx << std::endl;
      }
      double getEnergy(double *pos) const;

      void getGrad(double *pos, double *grad) const;
    private:
      int d_end1Idx, d_end2Idx; //!< indices of end points
      double d_minLen, d_maxLen;        //!< rest length of the bond
      double d_forceConstant;  //!< force constant of the bond

    };
  }
}
#endif
