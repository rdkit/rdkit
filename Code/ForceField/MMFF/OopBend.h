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
#ifndef __RD_MMFFOopBend_H__
#define __RD_MMFFOopBend_H__

#include <ForceField/Contrib.h>
#include <Geometry/point.h>


namespace ForceFields {
  namespace MMFF {
    class MMFFOop;

    //! the out-of-plane term for MMFF
    class OopBendContrib : public ForceFieldContrib {
    public:
      OopBendContrib() : d_at1Idx(-1), d_at2Idx(-1), d_at3Idx(-1), d_at4Idx(-1) {};
      //! Constructor
      /*!
	The Wilson angle is between the vector formed by atom2-atom4
  and the angle formed by atom1-atom2-atom3
	
	\param owner       pointer to the owning ForceField
	\param idx1        index of atom1 in the ForceField's positions
	\param idx2        index of atom2 in the ForceField's positions
	\param idx3        index of atom3 in the ForceField's positions
	\param idx4        index of atom4 in the ForceField's positions
      */
      OopBendContrib(ForceField *owner, unsigned int idx1, unsigned int idx2,
        unsigned int idx3, unsigned int idx4, const MMFFOop *mmffOopParams);
      double getEnergy(double *pos) const;
      void getGrad(double *pos, double *grad) const;
    private:
      int d_at1Idx, d_at2Idx, d_at3Idx, d_at4Idx;
      double d_koop;
    };

    namespace Utils {
      //! calculates and returns the Wilson angle (in degrees)
      double calcOopChi(const RDGeom::Point3D &iPoint, const RDGeom::Point3D &jPoint,
        const RDGeom::Point3D &kPoint, const RDGeom::Point3D &lPoint);
      //! returns the out-of-plane force constant koop
      double calcOopBendForceConstant(const MMFFOop *mmffOopParams);
      //! calculates and returns the out-of-plane MMFF energy
      double calcOopBendEnergy(const double chi, const double koop);
    }
  }
}
#endif
