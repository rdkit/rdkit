//
//  Copyright (C) 2024 Niels Maeder and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "DistanceConstraints.h"
#include "ForceField.h"
#include <cmath>
#include <RDGeneral/Invariant.h>

namespace ForceFields {
DistanceConstraintContribs::DistanceConstraintContribs(ForceField *owner) {
  PRECONDITION(owner, "bad owner");
  dp_forceField = owner;
}

void DistanceConstraintContribs::addContrib(unsigned int idx1,
                                            unsigned int idx2, double minLen,
                                            double maxLen,
                                            double forceConstant) {
  URANGE_CHECK(idx1, dp_forceField->positions().size());
  URANGE_CHECK(idx2, dp_forceField->positions().size());
  PRECONDITION(maxLen >= minLen, "bad bounds");
  d_contribs.emplace_back(idx1, idx2, minLen, maxLen, forceConstant);
}

void DistanceConstraintContribs::addContrib(unsigned int idx1,
                                            unsigned int idx2, bool relative,
                                            double minLen, double maxLen,
                                            double forceConstant) {
  const RDGeom::PointPtrVect &pos = dp_forceField->positions();
  URANGE_CHECK(idx1, pos.size());
  URANGE_CHECK(idx2, pos.size());
  PRECONDITION(maxLen >= minLen, "bad bounds");
  if (relative) {
    const RDGeom::Point3D p1 = *((RDGeom::Point3D *)pos[idx1]);
    const RDGeom::Point3D p2 = *((RDGeom::Point3D *)pos[idx2]);
    const auto distance = (p1 - p2).length();
    minLen = std::max(minLen + distance, 0.0);
    maxLen = std::max(maxLen + distance, 0.0);
  }
  d_contribs.emplace_back(idx1, idx2, minLen, maxLen, forceConstant);
}

double DistanceConstraintContribs::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");

  double accum = 0.0;
  for (const auto &contrib : d_contribs) {
    const auto distance2 =
        dp_forceField->distance2(contrib.idx1, contrib.idx2, pos);
    double difference = 0.0;
    if (distance2 < contrib.minLen * contrib.minLen) {
      difference = contrib.minLen - std::sqrt(distance2);
    } else if (distance2 > contrib.maxLen * contrib.maxLen) {
      difference = std::sqrt(distance2) - contrib.maxLen;
    } else {
      continue;
    }
    accum += 0.5 * contrib.forceConstant * difference * difference;
  }
  return accum;
}

void DistanceConstraintContribs::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");

  for (const auto &contrib : d_contribs) {
    double preFactor = 0.0;
    double distance = 0.0;
    const auto distance2 =
        dp_forceField->distance2(contrib.idx1, contrib.idx2, pos);
    if (distance2 < contrib.minLen * contrib.minLen) {
      distance = std::sqrt(distance2);
      preFactor = distance - contrib.minLen;
    } else if (distance2 > contrib.maxLen * contrib.maxLen) {
      distance = std::sqrt(distance2);
      preFactor = distance - contrib.maxLen;
    } else {
      continue;
    }
    preFactor *= contrib.forceConstant;
    preFactor /= std::max(1.0e-8, distance);
    const double *atom1Coords = &(pos[3 * contrib.idx1]);
    const double *atom2Coords = &(pos[3 * contrib.idx2]);
    for (unsigned int i = 0; i < 3; i++) {
      const double dGrad = preFactor * (atom1Coords[i] - atom2Coords[i]);
      grad[3 * contrib.idx1 + i] += dGrad;
      grad[3 * contrib.idx2 + i] -= dGrad;
    }
  }
}

}  // namespace ForceFields