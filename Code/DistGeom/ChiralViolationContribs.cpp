//
//  Copyright (C) 2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "ChiralViolationContribs.h"
#include "ChiralSet.h"
#include <ForceField/ForceField.h>

namespace DistGeom {
double calcChiralVolume(const unsigned int idx1, const unsigned int idx2,
                        const unsigned int idx3, const unsigned int idx4,
                        const double *pos, const unsigned int dim) {
  // even if we are minimizing in higher dimension the chiral volume is
  // calculated using only the first 3 dimensions
  RDGeom::Point3D v1(pos[idx1 * dim] - pos[idx4 * dim],
                     pos[idx1 * dim + 1] - pos[idx4 * dim + 1],
                     pos[idx1 * dim + 2] - pos[idx4 * dim + 2]);

  RDGeom::Point3D v2(pos[idx2 * dim] - pos[idx4 * dim],
                     pos[idx2 * dim + 1] - pos[idx4 * dim + 1],
                     pos[idx2 * dim + 2] - pos[idx4 * dim + 2]);

  RDGeom::Point3D v3(pos[idx3 * dim] - pos[idx4 * dim],
                     pos[idx3 * dim + 1] - pos[idx4 * dim + 1],
                     pos[idx3 * dim + 2] - pos[idx4 * dim + 2]);

  RDGeom::Point3D v2xv3 = v2.crossProduct(v3);

  double vol = v1.dotProduct(v2xv3);
  return vol;
}

ChiralViolationContribs::ChiralViolationContribs(
    ForceFields::ForceField *owner) {
  PRECONDITION(owner, "bad owner");
  dp_forceField = owner;
}
void ChiralViolationContribs::addContrib(const ChiralSet *cset, double weight) {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(cset, "bad chiral set");

  URANGE_CHECK(cset->d_idx1, dp_forceField->positions().size());
  URANGE_CHECK(cset->d_idx2, dp_forceField->positions().size());
  URANGE_CHECK(cset->d_idx3, dp_forceField->positions().size());
  URANGE_CHECK(cset->d_idx4, dp_forceField->positions().size());

  d_contribs.emplace_back(cset->d_idx1, cset->d_idx2, cset->d_idx3,
                          cset->d_idx4, cset->getUpperVolumeBound(),
                          cset->getLowerVolumeBound(), weight);
}

double ChiralViolationContribs::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");

  const unsigned int dim = dp_forceField->dimension();
  double res = 0.0;
  for (const auto &c : d_contribs) {
    double vol = calcChiralVolume(c.idx1, c.idx2, c.idx3, c.idx4, pos, dim);
    if (vol < c.volLower) {
      res += c.weight * (vol - c.volLower) * (vol - c.volLower);
    } else if (vol > c.volUpper) {
      res += c.weight * (vol - c.volUpper) * (vol - c.volUpper);
    }
  }

  return res;
}

void ChiralViolationContribs::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");

  const unsigned int dim = dp_forceField->dimension();

  for (const auto &c : d_contribs) {
    // even if we are minimizing in higher dimension the chiral volume is
    // calculated using only the first 3 dimensions
    RDGeom::Point3D v1(pos[c.idx1 * dim] - pos[c.idx4 * dim],
                       pos[c.idx1 * dim + 1] - pos[c.idx4 * dim + 1],
                       pos[c.idx1 * dim + 2] - pos[c.idx4 * dim + 2]);

    RDGeom::Point3D v2(pos[c.idx2 * dim] - pos[c.idx4 * dim],
                       pos[c.idx2 * dim + 1] - pos[c.idx4 * dim + 1],
                       pos[c.idx2 * dim + 2] - pos[c.idx4 * dim + 2]);

    RDGeom::Point3D v3(pos[c.idx3 * dim] - pos[c.idx4 * dim],
                       pos[c.idx3 * dim + 1] - pos[c.idx4 * dim + 1],
                       pos[c.idx3 * dim + 2] - pos[c.idx4 * dim + 2]);

    RDGeom::Point3D v2xv3 = v2.crossProduct(v3);

    double vol = v1.dotProduct(v2xv3);
    double preFactor;

    if (vol < c.volLower) {
      preFactor = c.weight * (vol - c.volLower);
    } else if (vol > c.volUpper) {
      preFactor = c.weight * (vol - c.volUpper);
    } else {
      continue;
    }

    // now comes the hard part - there are a total of 12 variables involved
    // 4 x 3 - four points and 3 dimensions
    //
    grad[dim * c.idx1] += preFactor * ((v2.y) * (v3.z) - (v3.y) * (v2.z));
    grad[dim * c.idx1 + 1] += preFactor * ((v3.x) * (v2.z) - (v2.x) * (v3.z));
    grad[dim * c.idx1 + 2] += preFactor * ((v2.x) * (v3.y) - (v3.x) * (v2.y));

    grad[dim * c.idx2] += preFactor * ((v3.y) * (v1.z) - (v3.z) * (v1.y));
    grad[dim * c.idx2 + 1] += preFactor * ((v3.z) * (v1.x) - (v3.x) * (v1.z));
    grad[dim * c.idx2 + 2] += preFactor * ((v3.x) * (v1.y) - (v3.y) * (v1.x));

    grad[dim * c.idx3] += preFactor * ((v2.z) * (v1.y) - (v2.y) * (v1.z));
    grad[dim * c.idx3 + 1] += preFactor * ((v2.x) * (v1.z) - (v2.z) * (v1.x));
    grad[dim * c.idx3 + 2] += preFactor * ((v2.y) * (v1.x) - (v2.x) * (v1.y));

    grad[dim * c.idx4] +=
        preFactor * (pos[c.idx1 * dim + 2] *
                         (pos[c.idx2 * dim + 1] - pos[c.idx3 * dim + 1]) +
                     pos[c.idx2 * dim + 2] *
                         (pos[c.idx3 * dim + 1] - pos[c.idx1 * dim + 1]) +
                     pos[c.idx3 * dim + 2] *
                         (pos[c.idx1 * dim + 1] - pos[c.idx2 * dim + 1]));

    grad[dim * c.idx4 + 1] +=
        preFactor *
        (pos[c.idx1 * dim] * (pos[c.idx2 * dim + 2] - pos[c.idx3 * dim + 2]) +
         pos[c.idx2 * dim] * (pos[c.idx3 * dim + 2] - pos[c.idx1 * dim + 2]) +
         pos[c.idx3 * dim] * (pos[c.idx1 * dim + 2] - pos[c.idx2 * dim + 2]));

    grad[dim * c.idx4 + 2] +=
        preFactor *
        (pos[c.idx1 * dim + 1] * (pos[c.idx2 * dim] - pos[c.idx3 * dim]) +
         pos[c.idx2 * dim + 1] * (pos[c.idx3 * dim] - pos[c.idx1 * dim]) +
         pos[c.idx3 * dim + 1] * (pos[c.idx1 * dim] - pos[c.idx2 * dim]));
  }
}
}  // namespace DistGeom
