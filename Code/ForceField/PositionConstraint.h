//
//  Copyright (C) 2013-2024 Paolo Tosco and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_POSITIONCONSTRAINT_H
#define RD_POSITIONCONSTRAINT_H
#include "Contrib.h"
#include <Geometry/point.h>

namespace ForceFields {

//! A position constraint of the type 0.5k * deltaX^2
class RDKIT_FORCEFIELD_EXPORT PositionConstraintContrib
    : public ForceFieldContrib {
 public:
  PositionConstraintContrib() {}
  //! Constructor
  /*!
    \param owner       pointer to the owning ForceField
    \param idx         index of the atom in the ForceField's positions
    \param minDispl    minimum displacement
    \param maxDispl    maximum displacement
    \param forceConst  force constant

  */
  PositionConstraintContrib(ForceField *owner, unsigned int idx,
                            double maxDispl, double forceConst);

  ~PositionConstraintContrib() override = default;
  double getEnergy(double *pos) const override;

  void getGrad(double *pos, double *grad) const override;
  PositionConstraintContrib *copy() const override {
    return new PositionConstraintContrib(*this);
  }

 private:
  int d_atIdx{-1};         //!< index of the restrained atom
  double d_maxDispl;       //!< maximum allowed displacement
  RDGeom::Point3D d_pos0;  //!< reference position
  double d_forceConstant;  //!< force constant of the bond
};
}  // namespace ForceFields
#endif
