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
#include <RDGeneral/export.h>
#ifndef __RD_UFFINVERSION_H__
#define __RD_UFFINVERSION_H__
#include <ForceField/Contrib.h>
#include <boost/tuple/tuple.hpp>
#include <Geometry/point.h>

namespace ForceFields {
namespace UFF {
class AtomicParams;

//! The inversion term for the Universal Force Field
class RDKIT_FORCEFIELD_EXPORT InversionContrib : public ForceFieldContrib {
 public:
  InversionContrib(){};
  //! Constructor
  /*!
    \param owner          pointer to the owning ForceField
    \param idx1           index of atom1 in the ForceField's positions
    \param idx2           index of atom2 in the ForceField's positions
    \param idx3           index of atom3 in the ForceField's positions
    \param idx4           index of atom4 in the ForceField's positions
    \param at2AtomicNum   atomic number for atom 2
    \param isCBoundToO    boolean flag; true if atom 2 is sp2 carbon bound to
    sp2 oxygen

  */
  InversionContrib(ForceField *owner, unsigned int idx1, unsigned int idx2,
                   unsigned int idx3, unsigned int idx4, int at2AtomicNum,
                   bool isCBoundToO, double oobForceScalingFactor = 1.0);

  double getEnergy(double *pos) const;

  void getGrad(double *pos, double *grad) const;
  virtual InversionContrib *copy() const {
    return new InversionContrib(*this);
  };

 private:
  int d_at1Idx{-1};
  int d_at2Idx{-1};
  int d_at3Idx{-1};
  int d_at4Idx{-1};
  double d_forceConstant, d_C0, d_C1, d_C2;
};

namespace Utils {
//! calculates and returns the cosine of the Y angle in an improper torsion
//! (see UFF paper, equation 19)
RDKIT_FORCEFIELD_EXPORT double calculateCosY(const RDGeom::Point3D &iPoint,
                                             const RDGeom::Point3D &jPoint,
                                             const RDGeom::Point3D &kPoint,
                                             const RDGeom::Point3D &lPoint);

//! calculates and returns the UFF force constant for an improper torsion
/*!

  \param at2AtomicNum   atomic number for atom 2
  \param isCBoundToO    boolean flag; true if atom 2 is sp2 carbon bound to sp2
  oxygen

  \return the force constant

*/
RDKIT_FORCEFIELD_EXPORT boost::tuple<double, double, double, double>
calcInversionCoefficientsAndForceConstant(int at2AtomicNum, bool isCBoundToO);
}  // namespace Utils
}  // namespace UFF
}  // namespace ForceFields
#endif
