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
#ifndef __RD_UFFINVERSION_H__
#define __RD_UFFINVERSION_H__
#include <ForceField/Contrib.h>
#include <boost/tuple/tuple.hpp>
#include <Geometry/point.h>

namespace ForceFields {
  namespace UFF {
    class AtomicParams;

    //! The inversion term for the Universal Force Field
    class InversionContrib : public ForceFieldContrib {
    public:
      InversionContrib() :  d_at1Idx(-1), d_at2Idx(-1), d_at3Idx(-1), d_at4Idx(-1) {};
      //! Constructor
      /*!
	\param owner          pointer to the owning ForceField
	\param idx1           index of atom1 in the ForceField's positions
	\param idx2           index of atom2 in the ForceField's positions
	\param idx3           index of atom3 in the ForceField's positions
	\param idx4           index of atom4 in the ForceField's positions
	\param at2AtomicNum   atomic number for atom 2
	\param isCBoundToO    boolean flag; true if atom 2 is sp2 carbon bound to sp2 oxygen
	
      */
      InversionContrib(ForceField *owner,
			  unsigned int idx1, unsigned int idx2,
			  unsigned int idx3, unsigned int idx4,
			  int at2AtomicNum, bool isCBoundToO);

      double getEnergy(double *pos) const;

      void getGrad(double *pos, double *grad) const;
    
    private:
      int d_at1Idx, d_at2Idx, d_at3Idx, d_at4Idx;
      double d_forceConstant, d_C0, d_C1, d_C2;

    };
  
    namespace Utils {
      //! calculates and returns the cosine of the Y angle in an improper torsion
      //! (see UFF paper, equation 19)
      double calculateCosY(const RDGeom::Point3D &iPoint,
        const RDGeom::Point3D &jPoint, const RDGeom::Point3D &kPoint,
        const RDGeom::Point3D &lPoint);
      
      //! calculates and returns the UFF force constant for an improper torsion
      /*!

	\param at2AtomicNum   atomic number for atom 2
	\param isCBoundToO    boolean flag; true if atom 2 is sp2 carbon bound to sp2 oxygen

	\return the force constant
	
      */
      boost::tuple<double, double, double, double>
        calcInversionCoefficientsAndForceConstant
        (int at2AtomicNum, bool isCBoundToO);

    }  
  }
}
#endif
