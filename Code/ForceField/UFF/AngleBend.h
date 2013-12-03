//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_ANGLEBEND_H__
#define __RD_ANGLEBEND_H__

#include <ForceField/Contrib.h>
#include <Geometry/point.h>

namespace ForceFields {
  namespace UFF {
    class AtomicParams;


    //! The angle-bend term for the Universal Force Field
    class AngleBendContrib : public ForceFieldContrib {
    public:
      AngleBendContrib() : d_at1Idx(-1), d_at2Idx(-1), d_at3Idx(-1),d_order(0) {};
      //! Constructor
      /*!
	The angle is between atom1 - atom2 - atom3
	
	\param owner       pointer to the owning ForceField
	\param idx1        index of atom1 in the ForceField's positions
	\param idx2        index of atom2 in the ForceField's positions
	\param idx3        index of atom3 in the ForceField's positions
	\param bondOrder12 order of the bond between atoms 1 and 2 (as a double)
	\param bondOrder23 order of the bond between atoms 2 and 3 (as a double)
	\param at1Params   pointer to the parameters for atom 1
	\param at2Params   pointer to the parameters for atom 2
	\param at3Params   pointer to the parameters for atom 3
	\param order (optional) the order of the angle term, this is for
  	   special cases and should adopt values:
	     - 0: not a special case, use the \c theta0 value from \c at2Params
	     - 2: linear coordination
	     - 3: trigonal planar coordination
	     - 4: square planar or tetrahedral coordination
	
      */
      AngleBendContrib(ForceField *owner,
		       unsigned int idx1,unsigned int idx2,unsigned int idx3,
		       double bondOrder12,double bondOrder23,
		       const AtomicParams *at1Params,
		       const AtomicParams *at2Params,
		       const AtomicParams *at3Params,
		       unsigned int order=0);
      double getEnergy(double *pos) const;
      void getGrad(double *pos,double *grad) const;
    
    private:
      int d_at1Idx,d_at2Idx,d_at3Idx;
      unsigned int d_order;
      double d_forceConstant,d_C0,d_C1,d_C2;

      double getEnergyTerm(double cosTheta,double sinThetaSq) const;
      double getThetaDeriv(double cosTheta,double sinTheta) const;
    };
  
    namespace Utils {
      //! Calculate the force constant for an angle bend
      /*!
	The angle is between atom1 - atom2 - atom3
	
	\param theta0      the preferred value of the angle (in radians)
	\param bondOrder12 order of the bond between atoms 1 and 2 (as a double)
	\param bondOrder23 order of the bond between atoms 2 and 3 (as a double)
	\param at1Params   pointer to the parameters for atom 1
	\param at2Params   pointer to the parameters for atom 2
	\param at3Params   pointer to the parameters for atom 3
	
      */
      double calcAngleForceConstant(double theta0,
				    double bondOrder12,double bondOrder23,
				    const AtomicParams *at1Params,
				    const AtomicParams *at2Params,
				    const AtomicParams *at3Params);
      void calcAngleBendGrad(RDGeom::Point3D *r, double *dist,
            double **g, double &dE_dTheta, double &cosTheta,
            double &sinTheta);
    }  
  }
}
#endif
