// $Id$
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include "PowerEigenSolver.h"
#include <Numerics/Matrix.h>
#include <Numerics/SquareMatrix.h>
#include <Numerics/SymmMatrix.h>
#include <Numerics/Vector.h>
#include <RDGeneral/utils.h>

//lapack ++ includes
#include <lafnames.h>
#include <lapack.h>
#include <symd.h>
#include <lavd.h>
#include <laslv.h>

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
  bool converge = powerEigenSolver(N, mat, eigVals, eigVecs );
  CHECK_INVARIANT(converge, "");
  DoubleVector ev1(N), ev2(N);
  eigVecs.getRow(0, ev1);
  eigVecs.getRow(1, ev2);
  CHECK_INVARIANT(RDKit::feq(ev1.dotProduct(ev2), 0.0, 0.001), "");
  
  // compare this solvers to the values we get from the lapack solver
  double data[25];
  double *leigs = new double[5];
  unsigned int i, j;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      data[i*N + j] = nmat.getVal(i,j);
    }
  }
  LaSymmMatDouble laMat(data, N, N);
  LaVectorDouble laEigs(leigs, N);
  LaGenMatDouble laEigVecs(N, N); //(leigvecs, N, N);
  
  LaEigSolve(laMat, laEigs, laEigVecs);
  
  DoubleVector nEigVals(N, Vector<double>::DATA_SPTR(leigs));
  double *ndata = new double[N*N];
  memcpy(static_cast<void *>(ndata), 
	 static_cast<const void *>(laEigVecs.addr()),
	 N*N*sizeof(double));
  Matrix<double>::DATA_SPTR sdata(ndata);
  DoubleSquareMatrix nEigVecs(N, sdata);
  
  // compare eigen values and eigen vectors
  for (i = 0; i < N; i++) {
    unsigned int ei = eigVals.largestValId();
    unsigned int nei = nEigVals.largestValId();
    CHECK_INVARIANT(RDKit::feq(eigVals.getVal(ei), nEigVals.getVal(nei), 0.01), "");
    eigVecs.getRow(ei, ev1);
    nEigVecs.getRow(nei, ev2);
    CHECK_INVARIANT(RDKit::feq(fabs(ev1.dotProduct(ev2)), 1.0, 0.01), "" );
    eigVals.setVal(ei, 0.0);
    nEigVals.setVal(nei, 0.0);
  }
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
  bool converge = powerEigenSolver(N, mat, eigVals, eigVecs, 100);
  CHECK_INVARIANT(converge, "");
  CHECK_INVARIANT(RDKit::feq(eigVals.getVal(0), 4.0, 0.001), "");

  double data[25];//, leigvecs[25];
  double *leigs = new double[5];
  unsigned int i, j;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      data[i*N + j] = nmat.getVal(i,j);
    }
  }
  LaSymmMatDouble laMat(data, N, N);
  LaVectorDouble laEigs(leigs, N);
  LaGenMatDouble laEigVecs(N, N); //(leigvecs, N, N);
  
  LaEigSolve(laMat, laEigs, laEigVecs);
  DoubleVector ev1(N), ev2(N);
  DoubleVector nEigVals(N, Vector<double>::DATA_SPTR(leigs));
  double *ndata = new double[N*N];
  memcpy(static_cast<void *>(ndata), 
	 static_cast<const void *>(laEigVecs.addr()),
	 N*N*sizeof(double));
  Matrix<double>::DATA_SPTR sdata(ndata);
  DoubleSquareMatrix nEigVecs(N, sdata);//laEigVecs.addr());
  for (i = 0; i < N; i++) {
    unsigned int ei = eigVals.largestValId();
    unsigned int nei = nEigVals.largestValId();
    CHECK_INVARIANT(RDKit::feq(eigVals.getVal(ei), nEigVals.getVal(nei), 0.001), "");
    
    if (i == 0) {
      eigVecs.getRow(ei, ev1);
      nEigVecs.getRow(nei, ev2);
      CHECK_INVARIANT(RDKit::feq(fabs(ev1.dotProduct(ev2)), 1.0, 0.001), "" );
    }
    eigVals.setVal(ei, 0.0);
    nEigVals.setVal(nei, 0.0);
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
  return (0);
}

