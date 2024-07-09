//
//  Copyright (C) 2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "DistViolationContribs.h"
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>

namespace DistGeom {

DistViolationContribs::DistViolationContribs(ForceFields::ForceField *owner) {
  PRECONDITION(owner, "bad owner");
  dp_forceField = owner;
}

double DistViolationContribs::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");

  double res = 0.0;
  for (const auto &c : d_contribs) {
    double d2 = dp_forceField->distance2(c.idx1, c.idx2, pos);
    double val = 0.0;
    if (d2 > c.ub * c.ub) {
      val = (d2 / (c.ub * c.ub)) - 1.0;
    } else if (d2 < c.lb * c.lb) {
      val = ((2 * c.lb * c.lb) / (c.lb * c.lb + d2)) - 1.0;
    }
    if (val > 0.0) {
      res += c.weight * val * val;
    }
  }
  return res;
}

void DistViolationContribs::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");
  const unsigned int dim = this->dp_forceField->dimension();
  for (const auto &contrib : d_contribs) {
    double d = dp_forceField->distance(contrib.idx1, contrib.idx2, pos);
    double preFactor = 0.0;
    if (d > contrib.ub) {
      double u2 = contrib.ub * contrib.ub;
      preFactor = 4. * (((d * d) / u2) - 1.0) * (d / u2);
    } else if (d < contrib.lb) {
      double d2 = d * d;
      double l2 = contrib.lb * contrib.lb;
      double l2d2 = d2 + l2;
      preFactor = 8. * l2 * d * (1. - 2 * l2 / l2d2) / (l2d2 * l2d2);
    } else {
      continue;
    }

    double *end1Coords = &(pos[dim * contrib.idx1]);
    double *end2Coords = &(pos[dim * contrib.idx2]);

    for (unsigned int i = 0; i < dim; i++) {
      double dGrad;
      if (d > 0.0) {
        dGrad =
            contrib.weight * preFactor * (end1Coords[i] - end2Coords[i]) / d;
      } else {
        // FIX: this likely isn't right
        dGrad = contrib.weight * preFactor * (end1Coords[i] - end2Coords[i]);
      }
      grad[dim * contrib.idx1 + i] += dGrad;
      grad[dim * contrib.idx2 + i] -= dGrad;
    }
  }
}
}  // namespace DistGeom
