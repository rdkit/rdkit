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
#ifndef __RD_MMFFANGLEBEND_H__
#define __RD_MMFFANGLEBEND_H__

#include <ForceField/ForceField.h>
#include <ForceField/Contrib.h>

namespace ForceFields {
  namespace MMFF {
    class MMFFBond;
    class MMFFAngle;
    class MMFFProp;

    //! The angle-bend term for MMFF
    class AngleBendContrib : public ForceFieldContrib {
    public:
      AngleBendContrib() : d_at1Idx(-1), d_at2Idx(-1), d_at3Idx(-1) {};
      //! Constructor
      /*!
	The angle is between atom1 - atom2 - atom3
	
	\param owner       pointer to the owning ForceField
	\param idx1        index of atom1 in the ForceField's positions
	\param idx2        index of atom2 in the ForceField's positions
	\param idx3        index of atom3 in the ForceField's positions
	\param angleType   MMFF type of the angle (as an unsigned int)
	
      */
      AngleBendContrib(ForceField *owner,
        unsigned int idx1, unsigned int idx2, unsigned int idx3,
        const MMFFAngle *mmffAngleParams, const MMFFProp *mmffPropParamsCentralAtom);
      double getEnergy(double *pos) const;
      void getGrad(double *pos,double *grad) const;
    
    private:
      bool d_isLinear;
      int d_at1Idx, d_at2Idx, d_at3Idx;
      double d_ka, d_theta0;
    };
    namespace Utils {
      //! returns the MMFF rest value for an angle
      double calcAngleRestValue(const MMFFAngle *mmffAngleParams);
      //! returns the MMFF force constant for an angle 
      double calcAngleForceConstant(const MMFFAngle *mmffAngleParams);
      //! calculates and returns the cosine of the angle between points p1, p2, p3
      double calcCosTheta(RDGeom::Point3D p1, RDGeom::Point3D p2,
        RDGeom::Point3D p3, double dist1, double dist2);
      //! calculates and returns the angle bending MMFF energy
      double calcAngleBendEnergy(const double theta0,
        const double ka, bool isLinear, const double cosTheta);
      void calcAngleBendGrad(RDGeom::Point3D *r, double *dist,
        double **g, double &dE_dTheta, double &cosTheta, double &sinTheta);
    }
  }
}
#endif
