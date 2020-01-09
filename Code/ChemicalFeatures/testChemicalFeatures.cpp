//
//  Copyright (C) 2005-2019 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//

#include <RDGeneral/test.h>
#include "FreeChemicalFeature.h"
#include <Geometry/point.h>
#include <iostream>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>

using namespace ChemicalFeatures;

void test1() {
  std::cout << "-----------------------------------------" << std::endl;
  std::cout << "Test1" << std::endl;

  FreeChemicalFeature f1("foo", "bar", RDGeom::Point3D(1, 1, 1));
  TEST_ASSERT(f1.getId() == -1);
  TEST_ASSERT(f1.getFamily() == "foo");
  TEST_ASSERT(f1.getType() == "bar");
  TEST_ASSERT(RDKit::feq(f1.getPos().x, 1.0));
  TEST_ASSERT(RDKit::feq(f1.getPos().y, 1.0));
  TEST_ASSERT(RDKit::feq(f1.getPos().z, 1.0));

  FreeChemicalFeature f2("foo", "bar", RDGeom::Point3D(1, 1, 1), 123);
  TEST_ASSERT(f2.getId() == 123);
  TEST_ASSERT(f2.getFamily() == "foo");
  TEST_ASSERT(f2.getType() == "bar");
  TEST_ASSERT(RDKit::feq(f2.getPos().x, 1.0));
  TEST_ASSERT(RDKit::feq(f2.getPos().y, 1.0));
  TEST_ASSERT(RDKit::feq(f2.getPos().z, 1.0));

  FreeChemicalFeature f3;
  f3.initFromString(f2.toString());
  TEST_ASSERT(f3.getId() == 123);
  TEST_ASSERT(f3.getFamily() == "foo");
  TEST_ASSERT(f3.getType() == "bar");
  TEST_ASSERT(RDKit::feq(f3.getPos().x, 1.0));
  TEST_ASSERT(RDKit::feq(f3.getPos().y, 1.0));
  TEST_ASSERT(RDKit::feq(f3.getPos().z, 1.0));

  FreeChemicalFeature f4(f2);
  TEST_ASSERT(f4.getId() == 123);
  TEST_ASSERT(f4.getFamily() == "foo");
  TEST_ASSERT(f4.getType() == "bar");
  TEST_ASSERT(RDKit::feq(f4.getPos().x, 1.0));
  TEST_ASSERT(RDKit::feq(f4.getPos().y, 1.0));
  TEST_ASSERT(RDKit::feq(f4.getPos().z, 1.0));

  std::cout << "Done" << std::endl;
}

int main() {
  test1();
  return 0;
}
