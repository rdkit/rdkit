//
//  Copyright (c) 2014-2025, Novartis Institutes for BioMedical Research and
//  other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>
#include "UniformRealValueGrid3D.h"
#include <GraphMol/GraphMol.h>
#include <Geometry/point.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/types.h>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ios>
#include <stdlib.h>

using namespace RDGeom;
using namespace RDKit;

TEST_CASE("test1UniformRealValueGrid3D") {
  RDGeom::UniformRealValueGrid3D grd(6.0, 5.0, 4.0);
  REQUIRE(grd.getSize() == 960);
  REQUIRE_THAT(grd.getSpacing(), Catch::Matchers::WithinAbs(.5, 1e-4));
  REQUIRE(grd.getNumX() == 12);
  REQUIRE(grd.getNumY() == 10);
  REQUIRE(grd.getNumZ() == 8);

  RDGeom::UniformRealValueGrid3D grd2(grd);
  REQUIRE(grd2.getSize() == 960);
  REQUIRE_THAT(grd2.getSpacing(), Catch::Matchers::WithinAbs(.5, 1e-4));
  REQUIRE(grd2.getNumX() == 12);
  REQUIRE(grd2.getNumY() == 10);
  REQUIRE(grd2.getNumZ() == 8);
  REQUIRE(grd2.getOccupancyVect()->getTotalVal() == 0);
  REQUIRE(grd.compareVectors(grd2));
  REQUIRE(grd.compareParams(grd2));

  grd.setVal(1, 1.0);
  REQUIRE(!grd.compareVectors(grd2));
  REQUIRE(grd.compareParams(grd2));
}

TEST_CASE("test2UniformRealValueGrid3DPickling") {
  RDGeom::UniformRealValueGrid3D grd(5.0, 5.0, 5.0, 0.1);
  grd.setVal(50, 2.3);
  const std::string pkl = grd.toString();

  RDGeom::UniformRealValueGrid3D grd2(pkl);
  REQUIRE(grd.compareVectors(grd2));
  REQUIRE(grd.compareParams(grd2));
}

TEST_CASE("test3UniformRealValueGrid3DIndexing") {
  RDGeom::UniformRealValueGrid3D grd(5.0, 5.0, 5.0, 0.1);
  {
    unsigned int xi = 3, yi = 3, zi = 3;
    unsigned int idx = grd.getGridIndex(xi, yi, zi);
    unsigned int nxi, nyi, nzi;
    grd.getGridIndices(idx, nxi, nyi, nzi);
    REQUIRE(nxi == xi);
    REQUIRE(nyi == yi);
    REQUIRE(nzi == zi);
  }
  {
    unsigned int xi = 3, yi = 3, zi = 5;
    unsigned int idx = grd.getGridIndex(xi, yi, zi);
    unsigned int nxi, nyi, nzi;
    grd.getGridIndices(idx, nxi, nyi, nzi);
    REQUIRE(nxi == xi);
    REQUIRE(nyi == yi);
    REQUIRE(nzi == zi);
  }
  {
    unsigned int xi = 3, yi = 6, zi = 3;
    unsigned int idx = grd.getGridIndex(xi, yi, zi);
    unsigned int nxi, nyi, nzi;
    grd.getGridIndices(idx, nxi, nyi, nzi);
    REQUIRE(nxi == xi);
    REQUIRE(nyi == yi);
    REQUIRE(nzi == zi);
  }
  {
    unsigned int xi = 0, yi = 0, zi = 0;
    unsigned int idx = grd.getGridIndex(xi, yi, zi);
    unsigned int nxi, nyi, nzi;
    grd.getGridIndices(idx, nxi, nyi, nzi);
    REQUIRE(nxi == xi);
    REQUIRE(nyi == yi);
    REQUIRE(nzi == zi);
  }
  {
    unsigned int xi = 8, yi = 2, zi = 1;
    unsigned int idx = grd.getGridIndex(xi, yi, zi);
    unsigned int nxi, nyi, nzi;
    grd.getGridIndices(idx, nxi, nyi, nzi);
    REQUIRE(nxi == xi);
    REQUIRE(nyi == yi);
    REQUIRE(nzi == zi);
  }

  RDGeom::Point3D pt = grd.getOffset();
  grd.setVal(pt, 2.3);
  unsigned int id = grd.getGridPointIndex(pt);

  REQUIRE_THAT(grd.getVal(pt), Catch::Matchers::WithinAbs(2.3, 1e-4));
  REQUIRE(id == 0);
  REQUIRE_THAT(grd.getVal(id), Catch::Matchers::WithinAbs(2.3, 1e-4));

  RDGeom::Point3D pt2 = grd.getGridPointLoc(id);
  REQUIRE_THAT(pt.x, Catch::Matchers::WithinAbs(pt2.x, 1e-4));
  REQUIRE_THAT(pt.y, Catch::Matchers::WithinAbs(pt2.y, 1e-4));
  REQUIRE_THAT(pt.z, Catch::Matchers::WithinAbs(pt2.z, 1e-4));
}

TEST_CASE("test4UniformRealValueGrid3DOps") {
  RDGeom::UniformRealValueGrid3D grd1(5.0, 5.0, 5.0, 0.1);
  RDGeom::UniformRealValueGrid3D grd2(5.0, 5.0, 5.0, 0.1);

  grd1.setVal(50, 37.37);
  grd2.setVal(50, 1.03);
  RDGeom::UniformRealValueGrid3D grd3(grd2);

  grd1 |= grd2;
  REQUIRE_THAT(grd1.getVal(50), Catch::Matchers::WithinAbs(37.37, 1e-4));
  REQUIRE_THAT(grd2.getVal(50), Catch::Matchers::WithinAbs(1.03, 1e-4));

  grd2 = grd2 | grd1;
  REQUIRE_THAT(grd1.getVal(50), Catch::Matchers::WithinAbs(37.37, 1e-4));
  REQUIRE_THAT(grd2.getVal(50), Catch::Matchers::WithinAbs(37.37, 1e-4));

  grd2 = grd2 & grd3;
  REQUIRE_THAT(grd2.getVal(50), Catch::Matchers::WithinAbs(1.03, 1e-4));
  REQUIRE_THAT(grd3.getVal(50), Catch::Matchers::WithinAbs(1.03, 1e-4));

  grd3 &= grd1;
  REQUIRE_THAT(grd1.getVal(50), Catch::Matchers::WithinAbs(37.37, 1e-4));
  REQUIRE_THAT(grd3.getVal(50), Catch::Matchers::WithinAbs(1.03, 1e-4));

  grd1 += grd2;
  REQUIRE_THAT(grd1.getVal(50), Catch::Matchers::WithinAbs(38.40, 1e-4));
  REQUIRE_THAT(grd2.getVal(50), Catch::Matchers::WithinAbs(1.03, 1e-4));

  grd1 -= grd2;
  REQUIRE_THAT(grd1.getVal(50), Catch::Matchers::WithinAbs(37.37, 1e-4));
  REQUIRE_THAT(grd2.getVal(50), Catch::Matchers::WithinAbs(1.03, 1e-4));

  grd2 = grd2 - grd1;
  REQUIRE_THAT(grd1.getVal(50), Catch::Matchers::WithinAbs(37.37, 1e-4));
  REQUIRE_THAT(grd2.getVal(50), Catch::Matchers::WithinAbs(-36.34, 1e-4));
}
