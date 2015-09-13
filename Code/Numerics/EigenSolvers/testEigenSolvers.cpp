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
#include "JacobiEigenSolver.h"
#include <Numerics/Matrix.h>
#include <Numerics/SquareMatrix.h>
#include <Numerics/SymmMatrix.h>
#include <Numerics/Vector.h>
#include <RDGeneral/utils.h>

#ifdef RDK_USE_EIGEN3
#include <Eigen/Dense>
#endif

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

void testJacobiSolver() {
  {
  unsigned int N = 2;
  DoubleSymmMatrix mat(N, 0.0);
  mat.setVal(0,0,2.0); mat.setVal(0,1,1.0);
  mat.setVal(1,0,1.0); mat.setVal(1,1,2.0);

  DoubleVector eigVals(N);
  DoubleSquareMatrix eigVecs(N);
  jacobiEigenSolver(mat, eigVals, eigVecs, 100);

  TEST_ASSERT(RDKit::feq(eigVals[0],1.0,0.001));
  TEST_ASSERT(RDKit::feq(eigVals[1],3.0,0.001));
  TEST_ASSERT(RDKit::feq(eigVecs.getVal(0,0),0.707,0.001));
  TEST_ASSERT(RDKit::feq(eigVecs.getVal(0,1),-0.707,0.001));
  TEST_ASSERT(RDKit::feq(eigVecs.getVal(1,0),0.707,0.001));
  TEST_ASSERT(RDKit::feq(eigVecs.getVal(1,1),0.707,0.001));
  }

  {
  unsigned int N = 3;
  DoubleSymmMatrix mat(N, 0.0);
  mat.setVal(0,0,2.0); mat.setVal(0,1,0.0); mat.setVal(0,2,0.0);
  mat.setVal(1,0,0.0); mat.setVal(1,1,3.0); mat.setVal(1,2,4.0);
  mat.setVal(2,0,0.0); mat.setVal(2,1,4.0); mat.setVal(2,2,9.0);

  DoubleVector eigVals(N);
  DoubleSquareMatrix eigVecs(N);
  jacobiEigenSolver(mat, eigVals, eigVecs, 100);

  TEST_ASSERT(RDKit::feq(eigVals[0],1.0,0.001));
  TEST_ASSERT(RDKit::feq(eigVals[1],2.0,0.001));
  TEST_ASSERT(RDKit::feq(eigVals[2],11.0,0.001));
  TEST_ASSERT(RDKit::feq(eigVecs.getVal(0,0),0.0,0.001));
  TEST_ASSERT(RDKit::feq(fabs(eigVecs.getVal(0,1)),0.894,0.001));
  TEST_ASSERT(RDKit::feq(fabs(eigVecs.getVal(0,2)),0.447,0.001));
  TEST_ASSERT(RDKit::feq(eigVecs.getVal(1,0),1.0,0.001));
  TEST_ASSERT(RDKit::feq(eigVecs.getVal(1,1),0.0,0.001));
  TEST_ASSERT(RDKit::feq(eigVecs.getVal(1,2),0.0,0.001));
  TEST_ASSERT(RDKit::feq(eigVecs.getVal(2,0),0.0,0.001));
  TEST_ASSERT(RDKit::feq(fabs(eigVecs.getVal(2,1)),0.447,0.001));
  TEST_ASSERT(RDKit::feq(fabs(eigVecs.getVal(2,2)),0.894,0.001));
  }

}

void testJacobiSolver2() {
  int seed = 31415;
  std::srand(seed);
  for ( unsigned int i = 2; i < 51; ++i) {
    int dataSize = i*(i+1)/2;
    int conv;
    double *data = new double[dataSize];
    DATA_SPTR dataPtr;
    dataPtr.reset(data);
    DoubleSymmMatrix mat(i, dataPtr);
    DoubleVector eigVals(i);
    DoubleSquareMatrix eigVecs(i);
#ifdef RDK_USE_EIGEN3
    Eigen::MatrixXd matE(i,i);
    Eigen::MatrixXd eigValsE(i, 1);
    Eigen::MatrixXd eigVecsE(i, i);
#endif

    for (unsigned int num = 0; num < 10; ++num) {
      for (unsigned int j = 0; j < i; j++) {
        for (unsigned int k = 0; k <= j; k++) {
          double res = static_cast<double>(rand())/RAND_MAX*1000.0;
          mat.setVal(j, k, res);
#ifdef RDK_USE_EIGEN3
          matE(j, k) = res;
          matE(k, j) = res;
#endif
        }
      }
      conv = jacobiEigenSolver(mat, eigVals, eigVecs, 1000);

#ifdef RDK_USE_EIGEN3
      Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > eigensolver(matE);
      eigVecsE = eigensolver.eigenvectors();
      eigValsE = eigensolver.eigenvalues();
      for (unsigned int j = 0; j < i; ++j) {
        for (unsigned int k = 0; k < i; ++k) {
          TEST_ASSERT(RDKit::feq(mat.getVal(j,k),matE(j,k), 0.001));
          if ((eigVecs.getVal(j, 0) * eigVecsE(0, j)) > 0 ) {
            TEST_ASSERT(RDKit::feq(eigVecs.getVal(j, k), eigVecsE(k, j), 0.01));
          }
          else {
            TEST_ASSERT(RDKit::feq(eigVecs.getVal(j, k), -eigVecsE(k, j), 0.01));
          }
        }TEST_ASSERT(RDKit::feq(eigVals.getVal(j), eigValsE(j), 0.001));
      }
#endif
    }
  }
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
  std::cout << "\t testJacobiSolver\n";
  testJacobiSolver();

  std::cout << "---------------------------------------\n";
  std::cout << "\t testJacobiSolver2\n";
  testJacobiSolver2();

  std::cout << "---------------------------------------\n";
  return (0);
}

