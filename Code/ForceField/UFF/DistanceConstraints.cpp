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
#include <cmath>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>

namespace ForceFields {
namespace UFF {
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
  double distance = 0.0;
  if (relative) {
    RDGeom::Point3D p1 = *((RDGeom::Point3D *)pos[idx1]);
    RDGeom::Point3D p2 = *((RDGeom::Point3D *)pos[idx2]);
    distance = (p1 - p2).length();
  }
  minLen = std::max(minLen + distance, 0.0);
  maxLen = std::max(maxLen + distance, 0.0);
  d_contribs.emplace_back(idx1, idx2, minLen, maxLen, forceConstant);
}

double DistanceConstraintContribs::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");

  double accum = 0.0;
  for (const auto &contrib : d_contribs) {
    double distance = dp_forceField->distance(contrib.idx1, contrib.idx2, pos);
    double difference = 0.0;
    if (distance < contrib.minLen) {
      difference = contrib.minLen - distance;
    } else if (distance > contrib.maxLen) {
      difference = distance - contrib.maxLen;
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
    double distance = dp_forceField->distance(contrib.idx1, contrib.idx2, pos);
    if (distance < contrib.minLen) {
      preFactor = distance - contrib.minLen;
    } else if (distance > contrib.maxLen) {
      preFactor = distance - contrib.maxLen;
    }
    preFactor *= contrib.forceConstant;

    double *atom1Coords = &(pos[3 * contrib.idx1]);
    double *atom2Coords = &(pos[3 * contrib.idx2]);
    for (unsigned int i = 0; i < 3; i++) {
      double dGrad = preFactor * (atom1Coords[i] - atom2Coords[i]) /
                     std::max(distance, 1.0e-8);
      grad[3 * contrib.idx1 + i] += dGrad;
      grad[3 * contrib.idx2 + i] -= dGrad;
    }
  }
}

}  // namespace UFF
}  // namespace ForceFields