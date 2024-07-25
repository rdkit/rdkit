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
#ifndef __RD_DISTVIOLATIONCONTRIB_H__
#define __RD_DISTVIOLATIONCONTRIB_H__

#include <ForceField/Contrib.h>

namespace DistGeom {
//! DEPRECATED: use DistViolationContribs instead
//! A term to capture the violation of the upper and lower bounds by
//! distance between two points
class RDKIT_DISTGEOMETRY_EXPORT DistViolationContrib
    : public ForceFields::ForceFieldContrib {
 public:
  DistViolationContrib() {}

  //! Constructor
  /*!
    \param owner       pointer to the owning ForceField
    \param idx1        index of end1 in the ForceField's positions
    \param idx2        index of end2 in the ForceField's positions
    \param ub          Upper bound on the distance
    \param lb          Lower bound on the distance
    \param weight      optional weight for this contribution
  */
  DistViolationContrib(ForceFields::ForceField *owner, unsigned int idx1,
                       unsigned int idx2, double ub, double lb,
                       double weight = 1.0);

  double getEnergy(double *pos) const override;

  void getGrad(double *pos, double *grad) const override;
  DistViolationContrib *copy() const override {
    return new DistViolationContrib(*this);
  }

 private:
  unsigned int d_end1Idx{0}, d_end2Idx{0};  //!< indices of end points
  double d_ub{
      1000.0};  //!< upper bound on the distance between d_end1Idx, d_end2Idx
  double d_lb{
      0.0};  //!< lower bound on the distance between d_end1Idx, d_end2Idx
  double d_weight{1.0};  //!< used to adjust relative contribution weights
};
}  // namespace DistGeom

#endif
