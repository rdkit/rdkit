//
//   Copyright (C) 2005-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>
#include "UniformGrid3D.h"
#include <DataStructs/DiscreteValueVect.h>
#include "point.h"
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>
#include "GridUtils.h"
#include <sstream>
#include <fstream>
#include <ios>
#include <cstdlib>

using namespace RDGeom;
using namespace RDKit;

TEST_CASE("testUniformGrid1") {
  UniformGrid3D grd(6.0, 5.0, 4.0);
  REQUIRE(grd.getSize() == 960);
  REQUIRE_THAT(grd.getSpacing(), Catch::Matchers::WithinAbs(.5, 1e-6));
  REQUIRE(grd.getNumX() == 12);
  REQUIRE(grd.getNumY() == 10);
  REQUIRE(grd.getNumZ() == 8);

  grd.setSphereOccupancy(Point3D(0.0, 0.0, 0.0), 1.5, 0.25);
  REQUIRE(grd.getOccupancyVect()->getTotalVal() == 523);

  UniformGrid3D grd2(grd);
  REQUIRE(grd2.getSize() == 960);
  REQUIRE_THAT(grd2.getSpacing(), Catch::Matchers::WithinAbs(.5, 1e-6));
  REQUIRE(grd2.getNumX() == 12);
  REQUIRE(grd2.getNumY() == 10);
  REQUIRE(grd2.getNumZ() == 8);
  REQUIRE(grd2.getOccupancyVect()->getTotalVal() == 523);

  grd.setSphereOccupancy(Point3D(1.0, 1.0, 0.0), 1.5, 0.25);
  REQUIRE(grd.getOccupancyVect()->getTotalVal() > 523);
  REQUIRE(grd2.getOccupancyVect()->getTotalVal() == 523);
}

TEST_CASE("testUniformGrid2") {
  UniformGrid3D grd(10.0, 10.0, 10.0);
  grd.setSphereOccupancy(Point3D(-2.0, -2.0, 0.0), 1.5, 0.25);
  grd.setSphereOccupancy(Point3D(-2.0, 2.0, 0.0), 1.5, 0.25);
  grd.setSphereOccupancy(Point3D(2.0, -2.0, 0.0), 1.5, 0.25);
  grd.setSphereOccupancy(Point3D(2.0, 2.0, 0.0), 1.5, 0.25);
  double dist = tanimotoDistance(grd, grd);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.0, 1e-6));
  dist = tverskyIndex(grd, grd, 1.0, 1.0);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(1.0, 1e-6));
  dist = tverskyIndex(grd, grd, 1.0, 0.0);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(1.0, 1e-6));
  dist = tverskyIndex(grd, grd, 0.0, 1.0);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(1.0, 1e-6));
  dist = tverskyIndex(grd, grd, 0.25, 0.75);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(1.0, 1e-6));
  dist = tverskyIndex(grd, grd, 0.75, 0.25);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(1.0, 1e-6));

  UniformGrid3D grd2(10.0, 10.0, 10.0);
  grd2.setSphereOccupancy(Point3D(-2.0, -2.0, 0.0), 1.5, 0.25);
  grd2.setSphereOccupancy(Point3D(-2.0, 2.0, 0.0), 1.5, 0.25);
  grd2.setSphereOccupancy(Point3D(2.0, -2.0, 0.0), 1.5, 0.25);
  dist = tanimotoDistance(grd, grd2);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.25, 1e-6));
  dist = tverskyIndex(grd, grd2, 1.0, 1.0);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.75, 1e-6));
  dist = tverskyIndex(grd, grd2, 1.0, 0.0);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.75, 1e-6));
  dist = tverskyIndex(grd, grd2, 0.0, 1.0);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(1.0, 1e-6));
  dist = tverskyIndex(grd, grd2, 0.25, 0.75);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.923, 1e-3));
  dist = tverskyIndex(grd, grd2, 0.75, 0.25);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.8, 1e-6));
  dist = protrudeDistance(grd, grd2);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.25, 1e-6));
  dist = protrudeDistance(grd2, grd);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.0, 1e-6));

  UniformGrid3D grd3(10.0, 10.0, 10.0);
  grd3.setSphereOccupancy(Point3D(-2.0, -2.0, 0.0), 1.5, 0.25);
  grd3.setSphereOccupancy(Point3D(-2.0, 2.0, 0.0), 1.5, 0.25);
  dist = tanimotoDistance(grd, grd3);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.5, 1e-6));
  dist = tverskyIndex(grd, grd3, 1.0, 1.0);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.5, 1e-6));
  dist = tverskyIndex(grd, grd3, 1.0, 0.0);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.5, 1e-6));
  dist = tverskyIndex(grd, grd3, 0.0, 1.0);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(1.0, 1e-6));
  dist = tverskyIndex(grd, grd3, 0.25, 0.75);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.8, 1e-6));
  dist = tverskyIndex(grd, grd3, 0.75, 0.25);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.5714, 1e-4));
  dist = protrudeDistance(grd, grd3);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.5, 1e-6));
  dist = protrudeDistance(grd3, grd);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.0, 1e-6));

  UniformGrid3D grd4(10.0, 10.0, 10.0);
  grd4.setSphereOccupancy(Point3D(2.0, 2.0, 0.0), 1.5, 0.25);
  grd4.setSphereOccupancy(Point3D(2.0, -2.0, 0.0), 1.5, 0.25);
  dist = tanimotoDistance(grd3, grd4);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(1.0, 1e-6));
  dist = tverskyIndex(grd3, grd4, 1.0, 1.0);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.0, 1e-6));
  dist = tverskyIndex(grd3, grd4, 1.0, 0.0);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.0, 1e-6));
  dist = tverskyIndex(grd3, grd4, 0.0, 1.0);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.0, 1e-6));
  dist = tverskyIndex(grd3, grd4, 0.25, 0.75);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.0, 1e-6));
  dist = tverskyIndex(grd3, grd4, 0.75, 0.25);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.0, 1e-6));

  UniformGrid3D grd5(10.0, 10.0, 10.0);
  grd5.setSphereOccupancy(Point3D(-2.0, -2.0, 0.0), 1.5, 0.25);
  dist = tanimotoDistance(grd, grd5);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.75, 1e-6));
  dist = tverskyIndex(grd, grd5, 1.0, 1.0);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.25, 1e-6));
  dist = tverskyIndex(grd, grd5, 1.0, 0.0);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.25, 1e-6));
  dist = tverskyIndex(grd, grd5, 0.0, 1.0);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(1.0, 1e-6));
  dist = tverskyIndex(grd, grd5, 0.25, 0.75);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.5714, 1e-4));
  dist = tverskyIndex(grd, grd5, 0.75, 0.25);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.3077, 1e-4));
  dist = protrudeDistance(grd, grd5);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.75, 1e-6));
  dist = protrudeDistance(grd5, grd);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.00, 1e-6));
}

TEST_CASE("testUniformGridPickling") {
  {
    UniformGrid3D grd(10.0, 10.0, 10.0);
    grd.setSphereOccupancy(Point3D(-2.0, -2.0, 0.0), 1.5, 0.25);
    grd.setSphereOccupancy(Point3D(-2.0, 2.0, 0.0), 1.5, 0.25);
    grd.setSphereOccupancy(Point3D(2.0, -2.0, 0.0), 1.5, 0.25);
    grd.setSphereOccupancy(Point3D(2.0, 2.0, 0.0), 1.5, 0.25);
    UniformGrid3D grd2(grd.toString());
    double dist = tanimotoDistance(grd, grd2);
    REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.0, 1e-6));
  }

  {
    std::string dirName = getenv("RDBASE");
    dirName += "/Code/Geometry/testData/";
    std::string pklName = dirName + "grid1.bin";
    std::ifstream inS;
    inS.open(pklName.c_str(), std::ios_base::binary);
    unsigned int length;
    inS >> length;
    auto *buff = new char[length];
    unsigned int nRead = 0;
    while (nRead < length) {
      nRead += inS.readsome(buff + nRead, length - nRead);
    }
    inS.close();
    std::string pkl(buff, length);
    delete[] buff;
    UniformGrid3D grd(pkl);

    UniformGrid3D grd2(10.0, 10.0, 10.0);
    grd2.setSphereOccupancy(Point3D(-2.0, -2.0, 0.0), 1.5, 0.25);
    grd2.setSphereOccupancy(Point3D(-2.0, 2.0, 0.0), 1.5, 0.25);
    grd2.setSphereOccupancy(Point3D(2.0, -2.0, 0.0), 1.5, 0.25);
    grd2.setSphereOccupancy(Point3D(2.0, 2.0, 0.0), 1.5, 0.25);

    std::string pkl2 = grd2.toString();
    REQUIRE(pkl2.length() == pkl.length());
    REQUIRE(pkl2 == pkl);

    REQUIRE(grd.getSize() == grd2.getSize());
    REQUIRE(grd.getNumX() == grd2.getNumX());
    REQUIRE(grd.getNumY() == grd2.getNumY());
    REQUIRE(grd.getNumZ() == grd2.getNumZ());
    REQUIRE(grd.compareParams(grd2));
    double dist = tanimotoDistance(grd, grd2);
    REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.0, 1e-6));
  }
}

TEST_CASE("testUniformGridOps") {
  UniformGrid3D grd(10.0, 10.0, 10.0);
  grd.setSphereOccupancy(Point3D(-2.0, -2.0, 0.0), 1.0, 0.25);
  grd.setSphereOccupancy(Point3D(-2.0, 2.0, 0.0), 1.0, 0.25);

  UniformGrid3D grd2(10.0, 10.0, 10.0);
  grd2.setSphereOccupancy(Point3D(2.0, -2.0, 0.0), 1.0, 0.25);
  grd2.setSphereOccupancy(Point3D(2.0, 2.0, 0.0), 1.0, 0.25);

  double dist = tanimotoDistance(grd, grd2);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(1.0, 1e-6));

  UniformGrid3D grd3(grd);
  grd3 |= grd2;

  dist = tanimotoDistance(grd3, grd);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.5, 1e-6));
  dist = tanimotoDistance(grd3, grd2);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.5, 1e-6));

  UniformGrid3D grd4(10.0, 10.0, 10.0);
  grd4.setSphereOccupancy(Point3D(-2.0, -2.0, 0.0), 1.0, 0.25);
  grd4.setSphereOccupancy(Point3D(-2.0, 2.0, 0.0), 1.0, 0.25);
  grd4.setSphereOccupancy(Point3D(2.0, -2.0, 0.0), 1.0, 0.25);

  UniformGrid3D grd5(grd4);
  grd5 &= grd2;

  dist = tanimotoDistance(grd5, grd);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(1.0, 1e-6));
  dist = tanimotoDistance(grd5, grd2);
  REQUIRE_THAT(dist, Catch::Matchers::WithinAbs(0.5, 1e-6));
}

TEST_CASE("testUniformGridIndexing") {
  UniformGrid3D grd(5.0, 5.0, 5.0);

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
}
