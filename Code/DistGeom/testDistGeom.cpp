// $Id$
//
// Copyright (C)  2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//


#include "BoundsMatrix.h"
#include "TriangleSmooth.h"
#include <iostream>
#include <boost/smart_ptr.hpp>
#include <math.h>
#include <Numerics/SymmMatrix.h>
#include "DistGeomUtils.h"
#include <RDGeneral/utils.h>

using namespace DistGeom;
using namespace RDNumeric;

void test1() {
  // test triangle smoothing
  unsigned int npt = 5;
  double x= sqrt(3.0); 
  BoundsMatrix * mmat = new BoundsMatrix(npt);

  mmat->setUpperBound(0, 1, 1.0); mmat->setLowerBound(0, 1, 1.0);
  mmat->setUpperBound(0, 2, x); mmat->setLowerBound(0, 2, x);
  mmat->setUpperBound(0, 3, 10.0); mmat->setLowerBound(0, 3, 0.0);
  mmat->setUpperBound(0, 4, 10.0); mmat->setLowerBound(0, 4, 0.0);
  mmat->setUpperBound(1, 2, 1.0); mmat->setLowerBound(1, 2, 1.0);
  mmat->setUpperBound(1, 3, x); mmat->setLowerBound(1, 3, x);
  mmat->setUpperBound(1, 4, 10.0); mmat->setLowerBound(1, 4, 0.0);
  mmat->setUpperBound(2, 3, 1.0); mmat->setLowerBound(2, 3, 1.0);
  mmat->setUpperBound(2, 4, x); mmat->setLowerBound(2, 4, x);
  mmat->setUpperBound(3, 4, 1.0); mmat->setLowerBound(3, 4, 1.0);
  
  BoundsMatPtr mptr(mmat);

  triangleSmoothBounds(mptr);
  CHECK_INVARIANT(RDKit::feq(mmat->getUpperBound(0, 1), 1.0, 0.001), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getLowerBound(0, 1), 1.0, 0.001), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getUpperBound(0, 2), 1.732, 0.001), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getLowerBound(0, 2), 1.732, 0.001), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getUpperBound(0, 3), 2.732, 0.001), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getLowerBound(0, 3), 0.732, 0.001), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getUpperBound(0, 4), 3.464, 0.001), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getLowerBound(0, 4), 0.0, 0.001), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getUpperBound(1, 2), 1.0, 0.001), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getLowerBound(1, 2), 1.0, 0.001), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getUpperBound(1, 3), 1.732, 0.001), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getLowerBound(1, 3), 1.732, 0.001), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getUpperBound(1, 4), 2.732, 0.001), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getLowerBound(1, 4), 0.732, 0.001), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getUpperBound(2, 3), 1.0, 0.001), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getLowerBound(2, 3), 1.0, 0.001), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getUpperBound(2, 4), 1.732, 0.001), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getLowerBound(2, 4), 1.732, 0.001), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getUpperBound(3, 4), 1.0, 0.001), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getLowerBound(3, 4), 1.0, 0.001), "");

  DoubleSymmMatrix dmat(npt, 0.0);
  RDKit::rng_type generator(42u);
  generator.seed(100);
  RDKit::uniform_double distrib(0,1.0);
  RDKit::double_source_type rng(generator,distrib);
  pickRandomDistMat(*mmat, dmat, rng);

  double sumElem = 0.0;
  for (unsigned int i = 0; i < dmat.getDataSize(); i++) {
    sumElem += dmat.getData()[i];
  }
  CHECK_INVARIANT(RDKit::feq(sumElem, 14.3079, 0.001), "");
}

void testIssue216() {
  RDNumeric::DoubleSymmMatrix dmat(4);
  dmat.setVal(0,0, 0.0); dmat.setVal(0,1, 1.0); dmat.setVal(0,2, 1.0); dmat.setVal(0,3, 1.0);
  dmat.setVal(1,1, 0.0); dmat.setVal(1,2, 1.0); dmat.setVal(1,3, 1.0);
  dmat.setVal(2,2, 0.0); dmat.setVal(2,3, 1.0);
  dmat.setVal(3,3, 0.0);
  
  std::cout << dmat;
  RDGeom::PointPtrVect pos;
  for (int i = 0; i < 4; i++) {
    RDGeom::Point3D *pt = new RDGeom::Point3D();
    pos.push_back(pt);
  }
  
  bool gotCoords = DistGeom::computeInitialCoords(dmat, pos);
  CHECK_INVARIANT(gotCoords, "");

  for (int i = 1; i < 4; i++) {
    RDGeom::Point3D pti = *(RDGeom::Point3D *)pos[i];
    for (int j = 0; j < i; j++) {
      RDGeom::Point3D ptj = *(RDGeom::Point3D *)pos[j];
      ptj -= pti;
      CHECK_INVARIANT(RDKit::feq(ptj.length(), 1.0, 0.02), "");
    }
  }
}
int main() {
  std::cout << "***********************************************************\n";
  std::cout << "   test1 \n";
  test1();

  std::cout << "***********************************************************\n";
  std::cout << "   testIssue216 \n";
  testIssue216();
  std::cout << "***********************************************************\n\n";
  return 0;
}
