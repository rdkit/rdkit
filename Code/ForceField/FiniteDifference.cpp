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
#include "Contrib.h"
#include "ForceField.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace ForceFields {

namespace {
void scatterPositions(const ForceField &ff, std::vector<double> &pos) {
  const unsigned int dim = ff.dimension();
  const unsigned int nPoints = ff.numPoints();
  pos.resize(dim * nPoints);
  for (unsigned int i = 0; i < nPoints; ++i) {
    for (unsigned int d = 0; d < dim; ++d) {
      pos[i * dim + d] = (*ff.positions()[i])[d];
    }
  }
}
}  // namespace

double calcFiniteDifference(ForceField &ff, double stepSize) {
  const unsigned int dim = ff.dimension();
  const unsigned int nCoords = dim * ff.numPoints();

  std::vector<double> analyticGrad(nCoords, 0.0);
  ff.calcGrad(analyticGrad.data());

  std::vector<double> pos;
  scatterPositions(ff, pos);

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

double calcContribFiniteDifference(const ForceFieldContrib &contrib,
                                   ForceField &ff, double stepSize) {
  const unsigned int dim = ff.dimension();
  const unsigned int nCoords = dim * ff.numPoints();

  std::vector<double> pos;
  scatterPositions(ff, pos);

  std::vector<double> analyticGrad(nCoords, 0.0);
  contrib.getGrad(pos.data(), analyticGrad.data());

  double maxDelta = 0.0;
  for (unsigned int i = 0; i < nCoords; ++i) {
    double orig = pos[i];

    pos[i] = orig + stepSize;
    ff.calcEnergy(pos.data());
    double ePlus = contrib.getEnergy(pos.data());

    pos[i] = orig - stepSize;
    ff.calcEnergy(pos.data());
    double eMinus = contrib.getEnergy(pos.data());

    pos[i] = orig;

    double fdGrad = (ePlus - eMinus) / (2.0 * stepSize);
    maxDelta = std::max(maxDelta, std::abs(fdGrad - analyticGrad[i]));
  }

  return maxDelta;
}

}  // namespace ForceFields
