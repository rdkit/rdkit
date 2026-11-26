//
//  Copyright (C) 2024-2025 Greg Landrum and other RDKit contributors
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

inline double distance2(const unsigned int idx1, const unsigned int idx2,
                        const double *pos, const unsigned int dim) {
  const auto *end1Coords = &(pos[dim * idx1]);
  const auto *end2Coords = &(pos[dim * idx2]);
  double d2 = 0.0;
  for (unsigned int i = 0; i < dim; i++) {
    double d = end1Coords[i] - end2Coords[i];
    d2 += d * d;
  }
  return d2;
}

inline double distance(const unsigned int idx1, const unsigned int idx2,
                       const double *pos, const unsigned int dim) {
  return sqrt(distance2(idx1, idx2, pos, dim));
}

namespace {
// When FixedDim is nonzero, the dimension is a compile-time constant, which
// allows the distance calculations to be fully unrolled for the common 3D and
// 4D cases. FixedDim == 0 means the dimension is taken from runtimeDim.
template <unsigned int FixedDim, typename Params>
double calcEnergy(const std::vector<Params> &contribs, const double *pos,
                  unsigned int runtimeDim) {
  const unsigned int dim = FixedDim ? FixedDim : runtimeDim;
  double accum = 0.0;
  for (const auto &c : contribs) {
    double d2 = distance2(c.idx1, c.idx2, pos, dim);
    double val = 0.0;
    if (d2 > c.ub2) {
      val = (d2 / (c.ub2)) - 1.0;
    } else if (d2 < c.lb2) {
      val = ((2 * c.lb2) / (c.lb2 + d2)) - 1.0;
    }
    if (val > 0.0) {
      accum += c.weight * val * val;
    }
  }
  return accum;
}

template <unsigned int FixedDim, typename Params>
void calcGrad(const std::vector<Params> &contribs, const double *pos,
              double *grad, unsigned int runtimeDim) {
  const unsigned int dim = FixedDim ? FixedDim : runtimeDim;
  for (const auto &c : contribs) {
    double d2 = distance2(c.idx1, c.idx2, pos, dim);
    double d;
    double preFactor = 0.0;
    if (d2 > c.ub2) {
      d = sqrt(d2);
      preFactor = 4. * (((d * d) / c.ub2) - 1.0) * (d / c.ub2);
    } else if (d2 < c.lb2) {
      d = sqrt(d2);
      double l2d2 = d2 + c.lb2;
      preFactor = 8. * c.lb2 * d * (1. - 2 * c.lb2 / l2d2) / (l2d2 * l2d2);
    } else {
      continue;
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
  }
}
}  // namespace

double DistViolationContribs::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");

  const unsigned int dim = dp_forceField->dimension();
  switch (dim) {
    case 3:
      return calcEnergy<3>(d_contribs, pos, dim);
    case 4:
      return calcEnergy<4>(d_contribs, pos, dim);
    default:
      return calcEnergy<0>(d_contribs, pos, dim);
  }
}

void DistViolationContribs::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");

  const unsigned int dim = dp_forceField->dimension();
  switch (dim) {
    case 3:
      calcGrad<3>(d_contribs, pos, grad, dim);
      break;
    case 4:
      calcGrad<4>(d_contribs, pos, grad, dim);
      break;
    default:
      calcGrad<0>(d_contribs, pos, grad, dim);
  }
}
}  // namespace DistGeom
