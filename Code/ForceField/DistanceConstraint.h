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
#ifndef RD_DISTANCECONSTRAINT_H
#define RD_DISTANCECONSTRAINT_H
#include "Contrib.h"

namespace ForceFields {

//! A distance range constraint modelled after a BondStretchContrib
class RDKIT_FORCEFIELD_EXPORT DistanceConstraintContrib
    : public ForceFieldContrib {
 public:
  DistanceConstraintContrib() {}
  //! Constructor
  /*!
    \param owner       pointer to the owning ForceField
    \param idx1        index of end1 in the ForceField's positions
    \param idx2        index of end2 in the ForceField's positions
    \param minLen      minimum distance
    \param maxLen      maximum distance
    \param forceConst  force Constant

  */
  DistanceConstraintContrib(ForceField *owner, unsigned int idx1,
                            unsigned int idx2, double minLen, double maxLen,
                            double forceConst);
  DistanceConstraintContrib(ForceField *owner, unsigned int idx1,
                            unsigned int idx2, bool relative, double minLen,
                            double maxLen, double forceConst);

  ~DistanceConstraintContrib() override {}
  double getEnergy(double *pos) const override;

  void getGrad(double *pos, double *grad) const override;
  DistanceConstraintContrib *copy() const override {
    return new DistanceConstraintContrib(*this);
  }

 private:
  int d_end1Idx{-1};          //!< indices of end points
  int d_end2Idx{-1};          //!< indices of end points
  double d_minLen, d_maxLen;  //!< rest length of the bond
  double d_forceConstant;     //!< force constant of the bond
};
}  // namespace ForceFields
#endif
