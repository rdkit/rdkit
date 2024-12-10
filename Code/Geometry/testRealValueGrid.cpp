//
//  Copyright (c) 2014-2024, Novartis Institutes for BioMedical Research and
//  other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
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

void test1UniformRealValueGrid3D() {
  RDGeom::UniformRealValueGrid3D grd(6.0, 5.0, 4.0);
  CHECK_INVARIANT(grd.getSize() == 960, "");
  CHECK_INVARIANT(RDKit::feq(grd.getSpacing(), .5), "");
  CHECK_INVARIANT(grd.getNumX() == 12, "");
  CHECK_INVARIANT(grd.getNumY() == 10, "");
  CHECK_INVARIANT(grd.getNumZ() == 8, "");

  RDGeom::UniformRealValueGrid3D grd2(grd);
  CHECK_INVARIANT(grd2.getSize() == 960, "");
  CHECK_INVARIANT(RDKit::feq(grd2.getSpacing(), .5), "");
  CHECK_INVARIANT(grd2.getNumX() == 12, "");
  CHECK_INVARIANT(grd2.getNumY() == 10, "");
  CHECK_INVARIANT(grd2.getNumZ() == 8, "");
  CHECK_INVARIANT(grd2.getOccupancyVect()->getTotalVal() == 0, "");
  TEST_ASSERT(grd.compareVectors(grd2));
  TEST_ASSERT(grd.compareParams(grd2));

  grd.setVal(1, 1.0);
  TEST_ASSERT(!grd.compareVectors(grd2));
  TEST_ASSERT(grd.compareParams(grd2))
  // make sure the data are actually decoupled:
  // grd.setSphereOccupancy(Point3D(1.0, 1.0, 0.0), 1.5, 0.25);
  // CHECK_INVARIANT(grd.getOccupancyVect()->getTotalVal()>523, "" );
  // CHECK_INVARIANT(grd2.getOccupancyVect()->getTotalVal() == 523, "" );
}

void test2UniformRealValueGrid3DPickling() {
  RDGeom::UniformRealValueGrid3D grd(5.0, 5.0, 5.0, 0.1);
  grd.setVal(50, 2.3);
  const std::string pkl = grd.toString();

  RDGeom::UniformRealValueGrid3D grd2(pkl);
  TEST_ASSERT(grd.compareVectors(grd2));
  TEST_ASSERT(grd.compareParams(grd2));
}

void test3UniformRealValueGrid3DIndexing() {
  RDGeom::UniformRealValueGrid3D grd(5.0, 5.0, 5.0, 0.1);
  {
    unsigned int xi = 3, yi = 3, zi = 3;
    unsigned int idx = grd.getGridIndex(xi, yi, zi);
    unsigned int nxi, nyi, nzi;
    grd.getGridIndices(idx, nxi, nyi, nzi);
    TEST_ASSERT(nxi == xi);
    TEST_ASSERT(nyi == yi);
    TEST_ASSERT(nzi == zi);
  }
  {
    unsigned int xi = 3, yi = 3, zi = 5;
    unsigned int idx = grd.getGridIndex(xi, yi, zi);
    unsigned int nxi, nyi, nzi;
    grd.getGridIndices(idx, nxi, nyi, nzi);
    TEST_ASSERT(nxi == xi);
    TEST_ASSERT(nyi == yi);
    TEST_ASSERT(nzi == zi);
  }
  {
    unsigned int xi = 3, yi = 6, zi = 3;
    unsigned int idx = grd.getGridIndex(xi, yi, zi);
    unsigned int nxi, nyi, nzi;
    grd.getGridIndices(idx, nxi, nyi, nzi);
    TEST_ASSERT(nxi == xi);
    TEST_ASSERT(nyi == yi);
    TEST_ASSERT(nzi == zi);
  }
  {
    unsigned int xi = 0, yi = 0, zi = 0;
    unsigned int idx = grd.getGridIndex(xi, yi, zi);
    unsigned int nxi, nyi, nzi;
    grd.getGridIndices(idx, nxi, nyi, nzi);
    TEST_ASSERT(nxi == xi);
    TEST_ASSERT(nyi == yi);
    TEST_ASSERT(nzi == zi);
  }
  {
    unsigned int xi = 8, yi = 2, zi = 1;
    unsigned int idx = grd.getGridIndex(xi, yi, zi);
    unsigned int nxi, nyi, nzi;
    grd.getGridIndices(idx, nxi, nyi, nzi);
    TEST_ASSERT(nxi == xi);
    TEST_ASSERT(nyi == yi);
    TEST_ASSERT(nzi == zi);
  }

  RDGeom::Point3D pt = grd.getOffset();
  grd.setVal(pt, 2.3);
  unsigned int id = grd.getGridPointIndex(pt);

  TEST_ASSERT(RDKit::feq(grd.getVal(pt), 2.3));
  TEST_ASSERT(id == 0);
  TEST_ASSERT(RDKit::feq(grd.getVal(id), 2.3));

  RDGeom::Point3D pt2 = grd.getGridPointLoc(id);
  TEST_ASSERT(RDKit::feq(pt.x, pt2.x));
  TEST_ASSERT(RDKit::feq(pt.y, pt2.y));
  TEST_ASSERT(RDKit::feq(pt.z, pt2.z));
}

void test4UniformRealValueGrid3DOps() {
  RDGeom::UniformRealValueGrid3D grd1(5.0, 5.0, 5.0, 0.1);
  RDGeom::UniformRealValueGrid3D grd2(5.0, 5.0, 5.0, 0.1);

  grd1.setVal(50, 37.37);
  grd2.setVal(50, 1.03);
  RDGeom::UniformRealValueGrid3D grd3(grd2);

  grd1 |= grd2;
  TEST_ASSERT(RDKit::feq(grd1.getVal(50), 37.37));
  TEST_ASSERT(RDKit::feq(grd2.getVal(50), 1.03));

  grd2 = grd2 | grd1;
  TEST_ASSERT(RDKit::feq(grd1.getVal(50), 37.37));
  TEST_ASSERT(RDKit::feq(grd2.getVal(50), 37.37));

  grd2 = grd2 & grd3;
  TEST_ASSERT(RDKit::feq(grd2.getVal(50), 1.03));
  TEST_ASSERT(RDKit::feq(grd3.getVal(50), 1.03));

  grd3 &= grd1;
  TEST_ASSERT(RDKit::feq(grd1.getVal(50), 37.37));
  TEST_ASSERT(RDKit::feq(grd3.getVal(50), 1.03));

  grd1 += grd2;
  TEST_ASSERT(RDKit::feq(grd1.getVal(50), 38.40));
  TEST_ASSERT(RDKit::feq(grd2.getVal(50), 1.03));

  grd1 -= grd2;
  TEST_ASSERT(RDKit::feq(grd1.getVal(50), 37.37));
  TEST_ASSERT(RDKit::feq(grd2.getVal(50), 1.03));

  grd2 = grd2 - grd1;
  TEST_ASSERT(RDKit::feq(grd1.getVal(50), 37.37));
  TEST_ASSERT(RDKit::feq(grd2.getVal(50), -36.34));
}

int main() {
  std::cout << "***********************************************************\n";
  std::cout << "Testing RealValueGrid3D\n";

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test1UniformRealValueGrid3D \n\n";
  test1UniformRealValueGrid3D();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test2UniformRealValueGrid3DPickling \n\n";
  test2UniformRealValueGrid3DPickling();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test3UniformRealValueGrid3DIndexing \n\n";
  test3UniformRealValueGrid3DIndexing();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test4UniformRealValueGrid3DOps \n\n";
  test4UniformRealValueGrid3DOps();

  return 0;
}
