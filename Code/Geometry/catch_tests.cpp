//
//   Copyright (C) 2019-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "catch.hpp"
#include <Geometry/point.h>
#include <Geometry/UniformGrid3D.h>

TEST_CASE("construct Point2D from Point3D", "[point]") {
  SECTION("basics") {
    RDGeom::Point3D p3(1., 2., 3.);
    RDGeom::Point2D p2(p3);
    CHECK(p2.x == p3.x);
    CHECK(p2.y == p3.y);
  }
}

TEST_CASE("UniformGrid getGridIndex") {
  RDGeom::UniformGrid3D grd(6.0, 5.0, 4.0);
  CHECK(grd.getSize() == 960);
  CHECK(grd.getNumX() == 12);
  CHECK(grd.getNumY() == 10);
  CHECK(grd.getNumZ() == 8);
  CHECK(grd.getGridIndex(12, 1, 1) == -1);
  CHECK(grd.getGridIndex(1, 10, 1) == -1);
  CHECK(grd.getGridIndex(1, 1, 8) == -1);
  CHECK(grd.getGridIndex(100, 1, 1) == -1);
  CHECK(grd.getGridIndex(1, 100, 1) == -1);
  CHECK(grd.getGridIndex(1, 1, 100) == -1);
  unsigned int x, y, z;
  CHECK_THROWS_AS(grd.getGridIndices(960, x, y, z), IndexErrorException);
  CHECK_THROWS_AS(grd.getGridPointLoc(960), IndexErrorException);
}

TEST_CASE("UniformGrid copying") {
  RDGeom::UniformGrid3D grd(6.0, 5.0, 4.0);
  grd.setSphereOccupancy(RDGeom::Point3D(0.0, 0.0, 0.0), 1.5, 0.25);
  CHECK(grd.getOccupancyVect()->getTotalVal() == 523);
  SECTION("operator=") {
    RDGeom::UniformGrid3D grd2(3, 3, 3);
    grd2 = grd;
    CHECK(grd2.getSize() == grd.getSize());
    CHECK(grd2.getOccupancyVect()->getTotalVal() ==
          grd.getOccupancyVect()->getTotalVal());
  }
  SECTION("char * ctor") {
    auto pkl = grd.toString();
    RDGeom::UniformGrid3D grd2(pkl.c_str(), pkl.size());
    CHECK(grd2.getSize() == grd.getSize());
    CHECK(grd2.getOccupancyVect()->getTotalVal() ==
          grd.getOccupancyVect()->getTotalVal());
  }
}

TEST_CASE("UniformGrid get/setVal") {
  RDGeom::UniformGrid3D grd(6.0, 5.0, 4.0);
  SECTION("getVal()") {
    {
      RDGeom::Point3D pt(1, 0, 0);
      CHECK(grd.getGridPointIndex(pt) >= 0);
      CHECK(grd.getVal(pt) == 0);
    }
    {
      RDGeom::Point3D pt(10, 0, 0);
      CHECK(grd.getGridPointIndex(pt) == -1);
      CHECK(grd.getVal(pt) == -1);
    }
    {
      RDGeom::Point3D pt(0, 10, 0);
      CHECK(grd.getGridPointIndex(pt) == -1);
      CHECK(grd.getVal(pt) == -1);
    }
    {
      RDGeom::Point3D pt(0, 0, 10);
      CHECK(grd.getGridPointIndex(pt) == -1);
      CHECK(grd.getVal(pt) == -1);
    }
  }
  SECTION("setVal") {
    CHECK(grd.getOccupancyVect()->getTotalVal() == 0);
    grd.setVal(RDGeom::Point3D(1, 0, 0), 2);
    CHECK(grd.getOccupancyVect()->getTotalVal() == 2);
    // not on the grid, has no impact
    grd.setVal(RDGeom::Point3D(10, 0, 0), 2);
    CHECK(grd.getOccupancyVect()->getTotalVal() == 2);
    grd.setVal(3, 1);
    CHECK(grd.getOccupancyVect()->getTotalVal() == 3);
    // not on the grid
    CHECK_THROWS_AS(grd.setVal(grd.getSize() + 1, 1), IndexErrorException);
    CHECK(grd.getOccupancyVect()->getTotalVal() == 3);
    // value too large
    CHECK_THROWS_AS(grd.setVal(5, 8), ValueErrorException);
    CHECK(grd.getOccupancyVect()->getTotalVal() == 3);
  }
  SECTION("setSphereOccupancy out of range") {
    grd.setSphereOccupancy(RDGeom::Point3D(10.0, 0.0, 0.0), 1.5, 0.25);
    CHECK(grd.getOccupancyVect()->getTotalVal() == 0);
    int maxLayers = -1;
    bool ignoreOutOfBound = false;
    CHECK_THROWS_AS(grd.setSphereOccupancy(RDGeom::Point3D(10.0, 0.0, 0.0), 1.5,
                                           0.25, maxLayers, ignoreOutOfBound),
                    RDGeom::GridException);
  }
}
TEST_CASE("compareParams") {
  RDGeom::UniformGrid3D grd(6.0, 5.0, 4.0);
  {
    RDGeom::UniformGrid3D grd2(6.0, 5.0, 4.0);
    CHECK(grd.compareParams(grd2));
    CHECK(grd2.compareParams(grd));
  }
  {
    RDGeom::UniformGrid3D grd2(7.0, 5.0, 4.0);
    CHECK(!grd.compareParams(grd2));
    CHECK(!grd2.compareParams(grd));
  }
  {
    RDGeom::UniformGrid3D grd2(6.0, 6.0, 4.0);
    CHECK(!grd.compareParams(grd2));
    CHECK(!grd2.compareParams(grd));
  }
  {
    RDGeom::UniformGrid3D grd2(6.0, 5.0, 5.0);
    CHECK(!grd.compareParams(grd2));
    CHECK(!grd2.compareParams(grd));
  }
  {
    RDGeom::UniformGrid3D grd2(6.6, 5.5, 4.4, grd.getSpacing() + .05);
    CHECK(!grd.compareParams(grd2));
    CHECK(!grd2.compareParams(grd));
  }
  {
    RDGeom::Point3D offset = grd.getOffset();
    offset *= 1.5;
    RDGeom::UniformGrid3D grd2(6.0, 5.0, 4.0, grd.getSpacing(),
                               RDKit::DiscreteValueVect::TWOBITVALUE, &offset);
    CHECK(!grd.compareParams(grd2));
    CHECK(!grd2.compareParams(grd));
  }
}