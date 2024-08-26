//
//  Copyright (C) 2004-2024 Paolo Tosco and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_ANGLECONSTRAINT_H
#define RD_ANGLECONSTRAINT_H
#include "Contrib.h"

namespace ForceFields {

//! An angle range constraint modelled after a AngleBendContrib
class RDKIT_FORCEFIELD_EXPORT AngleConstraintContrib
    : public ForceFieldContrib {
 public:
  AngleConstraintContrib() {}
  //! Constructor
  /*!
  \param owner       pointer to the owning ForceField
  \param idx1        index of atom1 in the ForceField's positions
  \param idx2        index of atom2 in the ForceField's positions
  \param idx3        index of atom3 in the ForceField's positions
  \param minAngle    minimum angle, must be between 0.0 and 180.0
  \param maxAngle    maximum angle, must be between 0.0 and 180.0
  \param forceConst  force Constant

  */
  AngleConstraintContrib(ForceField *owner, unsigned int idx1,
                         unsigned int idx2, unsigned int idx3,
                         double minAngleDeg, double maxAngleDeg,
                         double forceConst);
  AngleConstraintContrib(ForceField *owner, unsigned int idx1,
                         unsigned int idx2, unsigned int idx3, bool relative,
                         double minAngleDeg, double maxAngleDeg,
                         double forceConst);

  ~AngleConstraintContrib() override = default;
  double getEnergy(double *pos) const override;

  void getGrad(double *pos, double *grad) const override;

  AngleConstraintContrib *copy() const override {
    return new AngleConstraintContrib(*this);
  }

 private:
  double computeAngleTerm(double angle) const;
  int d_at1Idx{-1}, d_at2Idx{-1},
      d_at3Idx{-1};                     //!< indices of atoms forming the angle
  double d_minAngleDeg, d_maxAngleDeg;  //!< rest amplitudes of the angle
  double d_forceConstant;  //!< force constant of the angle constraint
};
}  // namespace ForceFields
#endif
