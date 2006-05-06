// $Id: testChemicalFeatures.cpp 4943 2006-02-17 01:21:27Z glandrum $
//
//  Copyright (C) 2005-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
//

#include "FreeChemicalFeature.h"
#include <Geometry/point.h>
#include <iostream>
#include <RDGeneral/Invariant.h>

using namespace ChemicalFeatures;


void test1() {
  std::cout << "-----------------------------------------" << std::endl;
  std::cout << "Test1" << std::endl;

  FreeChemicalFeature f1("foo","bar",RDGeom::Point3D(0,0,0));
  TEST_ASSERT(f1.getFamily()=="foo");
  TEST_ASSERT(f1.getType()=="bar");

  FreeChemicalFeature f2;
  f2.initFromString(f1.toString());
  TEST_ASSERT(f2.getFamily()=="foo");
  TEST_ASSERT(f2.getType()=="bar");

  FreeChemicalFeature f3(f1.toString());
  TEST_ASSERT(f3.getFamily()=="foo");
  TEST_ASSERT(f3.getType()=="bar");

  
  std::cout << "Done" << std::endl;

}

int main() {
  test1();
  return 0;
}

