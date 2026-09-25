//  Copyright (C) 2013-2025 Paolo Tosco and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef __RD_MMFFSTRETCHBEND_H__
#define __RD_MMFFSTRETCHBEND_H__

#include <cstdint>
#include <utility>
#include <vector>

#include <ForceField/Contrib.h>

namespace ForceFields {
namespace MMFF {
class MMFFBond;
class MMFFAngle;
class MMFFStbn;
class MMFFProp;

//! The angle-bend term for MMFF
class RDKIT_FORCEFIELD_EXPORT StretchBendContrib : public ForceFieldContrib {
 public:
  StretchBendContrib() {}
  //! Constructor
  StretchBendContrib(ForceField *owner);
  /*!
    Adds a stretch-bend term to the force field contrib.

    The angle is between atom1 - atom2 - atom3

    \param idx1        index of atom1 in the ForceField's positions
    \param idx2        index of atom2 in the ForceField's positions
    \param idx3        index of atom3 in the ForceField's positions
    \param angleType   MMFF type of the angle (as an unsigned int)
  */
  void addTerm(const unsigned int idx1, const unsigned int idx2,
               const unsigned int idx3, const MMFFStbn *mmffStbnParams,
               const MMFFAngle *mmffAngleParams,
               const MMFFBond *mmffBondParams1,
               const MMFFBond *mmffBondParams2);
  double getEnergy(double *pos) const override;
  void getGrad(double *pos, double *grad) const override;
  StretchBendContrib *copy() const override {
    return new StretchBendContrib(*this);
  }

 private:
  std::vector<int16_t> d_at1Idxs;
  std::vector<int16_t> d_at2Idxs;
  std::vector<int16_t> d_at3Idxs;
  std::vector<double> d_restLen1s;
  std::vector<double> d_restLen2s;
  std::vector<double> d_theta0s;
  std::vector<double> d_forceConstants1;
  std::vector<double> d_forceConstants2;
};
namespace Utils {
//! returns the std::pair of stretch-bend force constants for an angle
RDKIT_FORCEFIELD_EXPORT std::pair<double, double> calcStbnForceConstants(
    const MMFFStbn *mmffStbnParams);
//! calculates and returns the stretch-bending MMFF energy
RDKIT_FORCEFIELD_EXPORT std::pair<double, double> calcStretchBendEnergy(
    const double deltaDist1, const double deltaDist2, const double deltaTheta,
    const std::pair<double, double> forceConstants);
}  // namespace Utils
}  // namespace MMFF
}  // namespace ForceFields
#endif
