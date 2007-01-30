// $Id$
// 
//   Copyright (C) 2005-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include "GridUtils.h"
#include "Grid3D.h"
#include "UniformGrid3D.h"
#include <RDBoost/Exceptions.h>
#include <DataStructs/DiscreteValueVect.h>

using namespace RDKit;
namespace RDGeom {
  template<class GRIDTYPE> double tanimotoDistance(const GRIDTYPE &grid1, 
                                                   const GRIDTYPE &grid2) {
    if (!grid1.compareParams(grid2)) {
      throw ValueErrorException("Grid parameters do not match");
    }
    const DiscreteValueVect *v1 = grid1.getOccupancyVect();
    const DiscreteValueVect *v2 = grid2.getOccupancyVect();
    unsigned int dist = computeL1Norm(*v1, *v2);
    unsigned int totv1 = v1->getTotalVal();
    unsigned int totv2 = v2->getTotalVal();
    double inter = 0.5*(totv1 + totv2 - dist);
    double res = dist/(dist + inter);
    return res;
  }

  template double tanimotoDistance(const UniformGrid3D &grid1, 
                                   const UniformGrid3D &grid2);

  template<class GRIDTYPE> double protrudeDistance(const GRIDTYPE &grid1, 
                                                   const GRIDTYPE &grid2) {
    if (!grid1.compareParams(grid2)) {
      throw ValueErrorException("Grid parameters do not match");
    }
    const DiscreteValueVect *v1 = grid1.getOccupancyVect();
    const DiscreteValueVect *v2 = grid2.getOccupancyVect();
    unsigned int dist = computeL1Norm(*v1, *v2);
    unsigned int totv1 = v1->getTotalVal();
    unsigned int totv2 = v2->getTotalVal();
    double res = 1.0-(1.0*totv1-dist)/(1.0*totv1);
    return res;
  }

  template double protrudeDistance(const UniformGrid3D &grid1, 
                                   const UniformGrid3D &grid2);



}
