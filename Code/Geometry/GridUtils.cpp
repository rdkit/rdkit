// $Id$
//
//   Copyright (C) 2005-2007 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "GridUtils.h"
#include "Grid3D.h"
#include "UniformGrid3D.h"
#include "point.h"
#include <RDGeneral/Exceptions.h>
#include <DataStructs/DiscreteValueVect.h>
#include <cmath>

using namespace RDKit;
namespace RDGeom {

template <class GRIDTYPE>
double tverskyIndex(const GRIDTYPE &grid1, const GRIDTYPE &grid2, double alpha,
                    double beta) {
  if (!grid1.compareParams(grid2)) {
    throw ValueErrorException("Grid parameters do not match");
  }
  const DiscreteValueVect *v1 = grid1.getOccupancyVect();
  const DiscreteValueVect *v2 = grid2.getOccupancyVect();
  unsigned int dist = computeL1Norm(*v1, *v2);
  unsigned int totv1 = v1->getTotalVal();
  unsigned int totv2 = v2->getTotalVal();
  double inter = 0.5 * (totv1 + totv2 - dist);
  //  double alpha = 1.0;
  //  double beta = 1.0;
  double tversky_res = inter / (alpha * (1.0 * totv1 - inter) +
                                beta * (1.0 * totv2 - inter) + inter);
  //  double res = dist / (dist + inter);
  return tversky_res;
}

template RDKIT_RDGEOMETRYLIB_EXPORT double tverskyIndex(
    const UniformGrid3D &grid1, const UniformGrid3D &grid2, double alpha,
    double beta);

template <class GRIDTYPE>
double tanimotoDistance(const GRIDTYPE &grid1, const GRIDTYPE &grid2) {
  if (!grid1.compareParams(grid2)) {
    throw ValueErrorException("Grid parameters do not match");
  }
  const DiscreteValueVect *v1 = grid1.getOccupancyVect();
  const DiscreteValueVect *v2 = grid2.getOccupancyVect();
  unsigned int dist = computeL1Norm(*v1, *v2);
  unsigned int totv1 = v1->getTotalVal();
  unsigned int totv2 = v2->getTotalVal();
  double inter = 0.5 * (totv1 + totv2 - dist);
  double res = dist / (dist + inter);
  return res;
}

template RDKIT_RDGEOMETRYLIB_EXPORT double tanimotoDistance(
    const UniformGrid3D &grid1, const UniformGrid3D &grid2);

template <class GRIDTYPE>
double protrudeDistance(const GRIDTYPE &grid1, const GRIDTYPE &grid2) {
  if (!grid1.compareParams(grid2)) {
    throw ValueErrorException("Grid parameters do not match");
  }
  const DiscreteValueVect *v1 = grid1.getOccupancyVect();
  const DiscreteValueVect *v2 = grid2.getOccupancyVect();
  unsigned int totv1 = v1->getTotalVal();
  unsigned int totv2 = v2->getTotalVal();
  unsigned int totProtrude = computeL1Norm(*v1, *v2);
  unsigned int intersectVolume = (totv1 + totv2 - totProtrude) / 2;
  double res = (1.0 * totv1 - intersectVolume) / (1.0 * totv1);
  return res;
}

template RDKIT_RDGEOMETRYLIB_EXPORT double protrudeDistance(
    const UniformGrid3D &grid1, const UniformGrid3D &grid2);

std::map<int, std::vector<int>> gridIdxCache;
std::vector<int> computeGridIndices(const UniformGrid3D &grid,
                                    double windowRadius) {
  double gridSpacing = grid.getSpacing();
  int radInGrid = static_cast<int>(ceil(windowRadius / gridSpacing));
  // if(gridIdxCache.count(radInGrid)>0){
  //  return gridIdxCache[radInGrid];
  //}
  unsigned int dX = grid.getNumX();
  unsigned int dY = grid.getNumY();
  std::vector<int> res;
  for (int i = -radInGrid; i <= radInGrid; ++i) {
    for (int j = -radInGrid; j <= radInGrid; ++j) {
      for (int k = -radInGrid; k <= radInGrid; ++k) {
        double r2 = i * i + j * j + k * k;
        int d = static_cast<int>(sqrt(r2));
        if (d <= radInGrid) {
          res.push_back((i * dY + j) * dX + k);
        }
      }
    }
  }
  gridIdxCache[radInGrid] = res;
  return res;
}

Point3D computeGridCentroid(const UniformGrid3D &grid, const Point3D &pt,
                            double windowRadius, double &weightSum) {
  weightSum = 0.0;
  const DiscreteValueVect *v1 = grid.getOccupancyVect();
  Point3D centroid(0.0, 0.0, 0.0);

  unsigned int idxI = grid.getGridPointIndex(pt);
  std::vector<int> indicesInSphere = computeGridIndices(grid, windowRadius);
  for (const auto it : indicesInSphere) {
    int idx = idxI + it;
    if (idx >= 0 && static_cast<unsigned int>(idx) < v1->getLength()) {
      unsigned int wt = v1->getVal(idx);
      centroid += grid.getGridPointLoc(static_cast<unsigned int>(idx)) * wt;
      weightSum += wt;
    }
  }
  return centroid / weightSum;
}

std::vector<Point3D> findGridTerminalPoints(const UniformGrid3D &grid,
                                            double windowRadius,
                                            double inclusionFraction) {
  std::vector<Point3D> res;
  std::vector<int> indicesInSphere = computeGridIndices(grid, windowRadius);
  const DiscreteValueVect *storage = grid.getOccupancyVect();
  unsigned int maxGridVal = (0x1 << storage->getNumBitsPerVal()) - 1;
  for (unsigned int i = 0; i < storage->getLength(); ++i) {
    if (storage->getVal(i) < maxGridVal) {
      continue;
    }

    // -----------
    // compute the weighted volume of the shape inside the sphere:
    double volInSphere = 0.0;
    unsigned int nPtsHere = 0;
    for (const auto it : indicesInSphere) {
      int idx = i + it;
      if (idx >= 0 && static_cast<unsigned int>(idx) < storage->getLength()) {
        volInSphere += storage->getVal(static_cast<unsigned int>(idx));
        ++nPtsHere;
      }
    }
    // -----
    // the shape may be cut off by the edge of the grid, so
    // the actual max volume in the sphere may well be less
    // than the theoretical max:
    double maxPossValInSphere = nPtsHere * maxGridVal;
    if (volInSphere / maxPossValInSphere <= inclusionFraction) {
      Point3D ptI = grid.getGridPointLoc(i);
      double weightSum;
      Point3D centroid =
          computeGridCentroid(grid, ptI, windowRadius, weightSum);
      res.push_back(centroid);
    }
  }
  return res;
}
}  // namespace RDGeom
