// $Id$
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "PowerEigenSolver.h"
#include <Numerics/Matrix.h>
#include <Numerics/SquareMatrix.h>
#include <Numerics/SymmMatrix.h>
#include <Numerics/Vector.h>
#include <RDGeneral/utils.h>

using namespace RDNumeric;
using namespace RDNumeric::EigenSolvers;

void testPowerSolver() {
  unsigned int N = 5;
  DoubleSymmMatrix mat(N, 0.0);
  double x = 1.732; 
  double y = 2.268;
  double z = 3.268;
  mat.setVal(1,0,1.0);
  mat.setVal(2,0,x); mat.setVal(2,1,1.0);
  mat.setVal(3,0,y); mat.setVal(3,1,x); mat.setVal(3,2,1.0);
  mat.setVal(4,0,z); mat.setVal(4,1,y); mat.setVal(4,2,x); mat.setVal(4,3,1.0);

  DoubleSymmMatrix nmat(mat);
  DoubleMatrix eigVecs(N, N);
  DoubleVector eigVals(N);
  bool converge = powerEigenSolver(N, mat, eigVals, eigVecs,23);
  TEST_ASSERT(converge);
  DoubleVector ev1(N), ev2(N);
  eigVecs.getRow(0, ev1);
  eigVecs.getRow(1, ev2);
  TEST_ASSERT(RDKit::feq(ev1.dotProduct(ev2), 0.0, 0.001));
  
  TEST_ASSERT(RDKit::feq(eigVals[0],6.981,0.001));
  TEST_ASSERT(RDKit::feq(eigVals[1],-3.982,0.001));
  TEST_ASSERT(RDKit::feq(eigVals[2],-1.395,0.001));
  TEST_ASSERT(RDKit::feq(eigVals[3],-1.016,0.001));
  TEST_ASSERT(RDKit::feq(eigVals[4],-0.586,0.001));
  TEST_ASSERT(RDKit::feq(eigVecs.getVal(0,0),0.523,0.001));
  TEST_ASSERT(RDKit::feq(eigVecs.getVal(2,1),0.201,0.001));
  TEST_ASSERT(RDKit::feq(eigVecs.getVal(4,0),-.230,0.001));

}

void test2PowerSolver() {
  unsigned int N = 5;
  DoubleSymmMatrix mat(N, 0.0);
  double x = 1.0;
  mat.setVal(1,0,x);
  mat.setVal(2,0,x); mat.setVal(2,1,x);
  mat.setVal(3,0,x); mat.setVal(3,1,x); mat.setVal(3,2,x);
  mat.setVal(4,0,x); mat.setVal(4,1,x); mat.setVal(4,2,x); mat.setVal(4,3,x);
  DoubleSymmMatrix nmat(mat);
  
  DoubleVector eigVals(N);
  DoubleSquareMatrix eigVecs(N);
  DoubleVector ev1(N), ev2(N);
  bool converge = powerEigenSolver(N, mat, eigVals, eigVecs, 100);
  CHECK_INVARIANT(converge, "");
  CHECK_INVARIANT(RDKit::feq(eigVals.getVal(0), 4.0, 0.001), "");


  TEST_ASSERT(RDKit::feq(eigVals[0],4.000,0.001));
  TEST_ASSERT(RDKit::feq(eigVals[1],-1.0,0.001));
  TEST_ASSERT(RDKit::feq(eigVals[2],-1.0,0.001));
  TEST_ASSERT(RDKit::feq(eigVals[3],-1.0,0.001));
  TEST_ASSERT(RDKit::feq(eigVals[4],-1.0,0.001));
  TEST_ASSERT(RDKit::feq(eigVecs.getVal(0,0),0.447,0.001));
  TEST_ASSERT(RDKit::feq(eigVecs.getVal(2,1),0.028,0.001));
  TEST_ASSERT(RDKit::feq(eigVecs.getVal(4,0),0.193,0.001));
}

int main() {
  std::cout << "-----------------------------------------\n";
  std::cout << "Testing EigenSolvers code\n";

  std::cout << "---------------------------------------\n";
  std::cout << "\t testPowerSolver\n";
  testPowerSolver();

  std::cout << "---------------------------------------\n";
  std::cout << "\t test2PowerSolver\n";
  test2PowerSolver();
  std::cout << "---------------------------------------\n";
  return (0);
}

