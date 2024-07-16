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

auto distance2 = [](const unsigned int idx1, const unsigned int idx2,
                    const double *pos, const unsigned int dim) {
  const auto *end1Coords = &(pos[dim * idx1]);
  const auto *end2Coords = &(pos[dim * idx2]);
  double d2 = 0.0;
  for (unsigned int i = 0; i < dim; i++) {
    double d = end1Coords[i] - end2Coords[i];
    d2 += d * d;
  }
  return d2;
};
auto distance = [](const unsigned int idx1, const unsigned int idx2,
                   const double *pos, const unsigned int dim) {
  return sqrt(distance2(idx1, idx2, pos, dim));
};
double DistViolationContribs::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");

  double accum = 0.0;
  auto contrib = [&](const auto &c) {
    double d2 = distance2(c.idx1, c.idx2, pos, dp_forceField->dimension());
    double val = 0.0;
    if (d2 > c.ub * c.ub) {
      val = (d2 / (c.ub * c.ub)) - 1.0;
    } else if (d2 < c.lb * c.lb) {
      val = ((2 * c.lb * c.lb) / (c.lb * c.lb + d2)) - 1.0;
    }
    if (val > 0.0) {
      accum += c.weight * val * val;
    }
  };
  for (const auto &c : d_contribs) {
    contrib(c);
  }
  return accum;
}

void DistViolationContribs::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");
  const unsigned int dim = this->dp_forceField->dimension();

  auto contrib = [&](const auto &c) {
    double d = distance(c.idx1, c.idx2, pos, dp_forceField->dimension());
    double preFactor = 0.0;
    if (d > c.ub) {
      double u2 = c.ub * c.ub;
      preFactor = 4. * (((d * d) / u2) - 1.0) * (d / u2);
    } else if (d < c.lb) {
      double d2 = d * d;
      double l2 = c.lb * c.lb;
      double l2d2 = d2 + l2;
      preFactor = 8. * l2 * d * (1. - 2 * l2 / l2d2) / (l2d2 * l2d2);
    } else {
      return;
    }

    for (unsigned int i = 0; i < dim; i++) {
      const auto p1 = dim * c.idx1 + i;
      const auto p2 = dim * c.idx2 + i;
      double dGrad;
      if (d > 0.0) {
        dGrad = c.weight * preFactor * (pos[p1] - pos[p2]) / d;
      } else {
        // FIX: this likely isn't right
        dGrad = c.weight * preFactor * (pos[p1] - pos[p2]);
      }
      grad[p1] += dGrad;
      grad[p2] -= dGrad;
    }
  };
  for (const auto &c : d_contribs) {
    contrib(c);
  }
}
}  // namespace DistGeom
