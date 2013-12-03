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
#ifndef __RD_MMFFTORSIONANGLE_H__
#define __RD_MMFFTORSIONANGLE_H__

#include <ForceField/Contrib.h>
#include <boost/tuple/tuple.hpp>

namespace RDGeom {
  class Point3D;
}

namespace ForceFields {
  namespace MMFF {
    class MMFFTor;

    //! the torsion term for MMFF
    class TorsionAngleContrib : public ForceFieldContrib {
    public:
      TorsionAngleContrib() : d_at1Idx(-1), d_at2Idx(-1), d_at3Idx(-1), d_at4Idx(-1) {};
      //! Constructor
      /*!
	The torsion is between atom1 - atom2 - atom3 - atom4
	(i.e the angle between bond atom1-atom2 and bond atom3-atom4
	while looking down bond atom2-atom3)
	
	\param owner       pointer to the owning ForceField
	\param idx1        index of atom1 in the ForceField's positions
	\param idx2        index of atom2 in the ForceField's positions
	\param idx3        index of atom3 in the ForceField's positions
	\param idx4        index of atom4 in the ForceField's positions
	\param torsionType MMFF type of the torsional bond between atoms 2 and 3
      */
      TorsionAngleContrib(ForceField *owner, unsigned int idx1, unsigned int idx2,
        unsigned int idx3, unsigned int idx4, const MMFFTor *mmffTorParams);
      double getEnergy(double *pos) const;
      void getGrad(double *pos, double *grad) const;
    private:
      int d_at1Idx, d_at2Idx, d_at3Idx, d_at4Idx;
      double d_V1, d_V2, d_V3;
    };

    namespace Utils {
      //! calculates and returns the cosine of a torsion angle
      double calcTorsionCosPhi(const RDGeom::Point3D &iPoint,
        const RDGeom::Point3D &jPoint, const RDGeom::Point3D &kPoint,
        const RDGeom::Point3D &lPoint);
      //! returns the 3-tuple of a torsion angle force constants
      boost::tuple<double, double, double>
        calcTorsionForceConstant(const MMFFTor *mmffTorParams);
      //! calculates and returns the torsional MMFF energy
      double calcTorsionEnergy(const double V1,
        const double V2, const double V3, const double cosPhi);
      void calcTorsionGrad(RDGeom::Point3D *r, RDGeom::Point3D *t,
        double *d, double **g, double &sinTerm, double &cosPhi);
    }
  }
}
#endif
