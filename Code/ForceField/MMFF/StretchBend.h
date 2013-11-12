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
#ifndef __RD_MMFFSTRETCHBEND_H__
#define __RD_MMFFSTRETCHBEND_H__

#include <utility>
#include <ForceField/Contrib.h>

namespace ForceFields {
  namespace MMFF {
    class MMFFBond;
    class MMFFAngle;
    class MMFFStbn;
    class MMFFProp;

    //! The angle-bend term for MMFF
    class StretchBendContrib : public ForceFieldContrib {
    public:
      StretchBendContrib() : d_at1Idx(-1), d_at2Idx(-1), d_at3Idx(-1) {};
      //! Constructor
      /*!
	The angle is between atom1 - atom2 - atom3
	
	\param owner       pointer to the owning ForceField
	\param idx1        index of atom1 in the ForceField's positions
	\param idx2        index of atom2 in the ForceField's positions
	\param idx3        index of atom3 in the ForceField's positions
	\param angleType   MMFF type of the angle (as an unsigned int)
	
      */
    StretchBendContrib(ForceField *owner,
      const unsigned int idx1, const unsigned int idx2, const unsigned int idx3,
      const std::pair<bool, const MMFFStbn *> mmffStbnParams,
      const MMFFAngle *mmffAngleParams, const MMFFBond *mmffBondParams1,
      const MMFFBond *mmffBondParams2);
      
      double getEnergy(double *pos) const;
      void getGrad(double *pos,double *grad) const;
    
    private:
      int d_at1Idx, d_at2Idx, d_at3Idx;
      double d_restLen1, d_restLen2, d_theta0;
      std::pair<double, double> d_forceConstants;
    };
    namespace Utils {
      //! returns the std::pair of stretch-bend force constants for an angle 
      std::pair<double, double> calcStbnForceConstants
        (const std::pair<bool, const MMFFStbn *> mmffStbnParams);
      //! calculates and returns the stretch-bending MMFF energy
      std::pair<double, double> calcStretchBendEnergy
        (const double deltaDist1, const double deltaDist2,
        const double deltaTheta, const std::pair<double, double> forceConstants);
    }  
  }
}
#endif
