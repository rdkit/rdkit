//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "FiniteDifference.h"
#include "ForceField.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace ForceFields {

double calcFiniteDifference(ForceField &ff, double stepSize) {
  const unsigned int dim = ff.dimension();
  const unsigned int nPoints = ff.numPoints();
  const unsigned int nCoords = dim * nPoints;

  std::vector<double> analyticGrad(nCoords, 0.0);
  ff.calcGrad(analyticGrad.data());

  std::vector<double> pos(nCoords);
  for (unsigned int i = 0; i < nPoints; ++i) {
    for (unsigned int d = 0; d < dim; ++d) {
      pos[i * dim + d] = (*ff.positions()[i])[d];
    }
  }

  double maxDelta = 0.0;
  for (unsigned int i = 0; i < nCoords; ++i) {
    double orig = pos[i];

    pos[i] = orig + stepSize;
    double ePlus = ff.calcEnergy(pos.data());

    pos[i] = orig - stepSize;
    double eMinus = ff.calcEnergy(pos.data());

    pos[i] = orig;

    double fdGrad = (ePlus - eMinus) / (2.0 * stepSize);
    maxDelta = std::max(maxDelta, std::abs(fdGrad - analyticGrad[i]));
  }

  return maxDelta;
}

}  // namespace ForceFields
