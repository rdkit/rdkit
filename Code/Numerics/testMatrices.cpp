//
//  Copyright (C) 2004-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>
#include "Matrix.h"
#include "SquareMatrix.h"
#include "SymmMatrix.h"
#include <iostream>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <cmath>
#include <RDGeneral/utils.h>

using namespace RDNumeric;

TEST_CASE("test1Vector") {
  Vector<double> v1(3, 1.0);
  v1.setVal(0, 2.0);
  v1.setVal(2, -4.0);
  REQUIRE_THAT(v1.normL1(), Catch::Matchers::WithinAbs(7.0, 1e-4));
  REQUIRE_THAT(v1.normLinfinity(), Catch::Matchers::WithinAbs(4.0, 1e-4));
  REQUIRE_THAT(v1.normL2(), Catch::Matchers::WithinAbs(sqrt(21.0), 1e-4));

  v1.setVal(1, 2.0);
  REQUIRE_THAT(v1.getVal(1), Catch::Matchers::WithinAbs(2.0, 1e-4));
  REQUIRE_THAT(v1.normL1(), Catch::Matchers::WithinAbs(8.0, 1e-4));

  auto *data = new double[3];
  data[0] = 1.0;
  data[1] = 2.0;
  data[2] = 3.0;
  Vector<double>::DATA_SPTR sdata(data);
  Vector<double> *v2 = new Vector<double>(3, sdata);
  REQUIRE_THAT(v2->normL1(), Catch::Matchers::WithinAbs(6.0, 1e-4));
  Vector<double> v3(v1);
  unsigned int i;
  for (i = 0; i < v1.size(); i++) {
    REQUIRE_THAT(v1.getVal(i), Catch::Matchers::WithinAbs(v3.getVal(i), 1e-4));
  }

  delete v2;

  Vector<double> vr1(5);
  Vector<double> vr2(5);
  vr1.setToRandom();
  REQUIRE_THAT(vr1.normL2(), Catch::Matchers::WithinAbs(1.0, 1e-4));
  vr2.setToRandom(120);
  REQUIRE_THAT(vr2.normL2(), Catch::Matchers::WithinAbs(1.0, 1e-4));
}

TEST_CASE("test2Matrix") {
  Matrix<double> A(2, 3);
  A.setVal(0, 0, 1.0);
  A.setVal(0, 1, 0.5);
  A.setVal(0, 2, 2.0);
  A.setVal(1, 0, 0.5);
  A.setVal(1, 1, 1.0);
  A.setVal(1, 2, 3.0);

  Vector<double> v1(3, 1.0);
  v1.setVal(1, 2.0);
  v1.setVal(2, 3.0);

  Vector<double> v2(2);
  multiply(A, v1, v2);
  REQUIRE_THAT(v2.getVal(0), Catch::Matchers::WithinAbs(8.0, 1e-4));
  REQUIRE_THAT(v2.getVal(1), Catch::Matchers::WithinAbs(11.5, 1e-4));

  double *data = A.getData();
  data[2] = 3.0;
  REQUIRE_THAT(A.getVal(0, 2), Catch::Matchers::WithinAbs(3.0, 1e-4));
  multiply(A, v1, v2);
  REQUIRE_THAT(v2.getVal(0), Catch::Matchers::WithinAbs(11.0, 1e-4));
  REQUIRE_THAT(v2.getVal(1), Catch::Matchers::WithinAbs(11.5, 1e-4));

  data = new double[6];
  Matrix<double> *B = new Matrix<double>(2, 3, Matrix<double>::DATA_SPTR(data));
  Matrix<double> E(A);
  multiply(E, v1, v2);
  REQUIRE_THAT(v2.getVal(0), Catch::Matchers::WithinAbs(11.0, 1e-4));
  REQUIRE_THAT(v2.getVal(1), Catch::Matchers::WithinAbs(11.5, 1e-4));
  delete B;

  Matrix<double> D(3, 2);
  A.transpose(D);
  Matrix<double> C(2, 2);
  multiply(A, D, C);

  REQUIRE_THAT(C.getVal(0, 0), Catch::Matchers::WithinAbs(10.25, 1e-4));
  REQUIRE_THAT(C.getVal(0, 1), Catch::Matchers::WithinAbs(10.0, 1e-4));
  REQUIRE_THAT(C.getVal(1, 0), Catch::Matchers::WithinAbs(10.0, 1e-4));
  REQUIRE_THAT(C.getVal(1, 1), Catch::Matchers::WithinAbs(10.25, 1e-4));

  auto Ccp(C);
  Ccp += C;
  REQUIRE_THAT(Ccp.getVal(0, 0), Catch::Matchers::WithinAbs(20.5, 1e-4));
  REQUIRE_THAT(Ccp.getVal(0, 1), Catch::Matchers::WithinAbs(20.0, 1e-4));
  REQUIRE_THAT(Ccp.getVal(1, 0), Catch::Matchers::WithinAbs(20.0, 1e-4));
  REQUIRE_THAT(Ccp.getVal(1, 1), Catch::Matchers::WithinAbs(20.5, 1e-4));

  Ccp -= C;
  REQUIRE_THAT(Ccp.getVal(0, 0), Catch::Matchers::WithinAbs(10.25, 1e-4));
  REQUIRE_THAT(Ccp.getVal(0, 1), Catch::Matchers::WithinAbs(10.0, 1e-4));
  REQUIRE_THAT(Ccp.getVal(1, 0), Catch::Matchers::WithinAbs(10.0, 1e-4));
  REQUIRE_THAT(Ccp.getVal(1, 1), Catch::Matchers::WithinAbs(10.25, 1e-4));

  C *= 2.;
  REQUIRE_THAT(C.getVal(0, 0), Catch::Matchers::WithinAbs(20.5, 1e-4));
  REQUIRE_THAT(C.getVal(0, 1), Catch::Matchers::WithinAbs(20.0, 1e-4));
  REQUIRE_THAT(C.getVal(1, 0), Catch::Matchers::WithinAbs(20.0, 1e-4));
  REQUIRE_THAT(C.getVal(1, 1), Catch::Matchers::WithinAbs(20.5, 1e-4));

  C /= 2.;
  REQUIRE_THAT(C.getVal(0, 0), Catch::Matchers::WithinAbs(10.25, 1e-4));
  REQUIRE_THAT(C.getVal(0, 1), Catch::Matchers::WithinAbs(10.0, 1e-4));
  REQUIRE_THAT(C.getVal(1, 0), Catch::Matchers::WithinAbs(10.0, 1e-4));
  REQUIRE_THAT(C.getVal(1, 1), Catch::Matchers::WithinAbs(10.25, 1e-4));

  Vector<double> tRow(A.numCols());
  A.getRow(1, tRow);
  for (unsigned int i = 0; i < A.numCols(); ++i) {
    REQUIRE_THAT(A.getVal(1, i),
                 Catch::Matchers::WithinAbs(tRow.getVal(i), 1e-4));
  }
  Vector<double> tCol(A.numRows());
  A.getCol(1, tCol);
  for (unsigned int i = 0; i < A.numRows(); ++i) {
    REQUIRE_THAT(A.getVal(i, 1),
                 Catch::Matchers::WithinAbs(tCol.getVal(i), 1e-4));
  }
}

TEST_CASE("test3SquareMatrix") {
  SquareMatrix<double> A(2);
  A.setVal(0, 0, 1.0);
  A.setVal(0, 1, 2.0);
  A.setVal(1, 0, 3.0);
  A.setVal(1, 1, 4.0);
  SquareMatrix<double> B(A), C(2);

  multiply(A, B, C);
  REQUIRE_THAT(C.getVal(0, 0), Catch::Matchers::WithinAbs(7.0, 1e-4));
  REQUIRE_THAT(C.getVal(0, 1), Catch::Matchers::WithinAbs(10.0, 1e-4));
  REQUIRE_THAT(C.getVal(1, 0), Catch::Matchers::WithinAbs(15.0, 1e-4));
  REQUIRE_THAT(C.getVal(1, 1), Catch::Matchers::WithinAbs(22.0, 1e-4));
  B *= A;
  REQUIRE_THAT(B.getVal(0, 0), Catch::Matchers::WithinAbs(7.0, 1e-4));
  REQUIRE_THAT(B.getVal(0, 1), Catch::Matchers::WithinAbs(10.0, 1e-4));
  REQUIRE_THAT(B.getVal(1, 0), Catch::Matchers::WithinAbs(15.0, 1e-4));
  REQUIRE_THAT(B.getVal(1, 1), Catch::Matchers::WithinAbs(22.0, 1e-4));

  auto *data = new double[4];
  data[0] = 1.0;
  data[1] = 2.0;
  data[2] = 3.0;
  data[3] = 4.0;
  SquareMatrix<double> *D =
      new SquareMatrix<double>(2, Matrix<double>::DATA_SPTR(data));
  SquareMatrix<double> E(*D);
  multiply((*D), E, A);

  unsigned int i, j;
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      REQUIRE_THAT(B.getVal(i, j),
                   Catch::Matchers::WithinAbs(A.getVal(i, j), 1e-4));
    }
  }
  D->transposeInplace();
  REQUIRE_THAT(D->getVal(0, 0), Catch::Matchers::WithinAbs(1.0, 1e-4));
  REQUIRE_THAT(D->getVal(0, 1), Catch::Matchers::WithinAbs(3.0, 1e-4));
  REQUIRE_THAT(D->getVal(1, 0), Catch::Matchers::WithinAbs(2.0, 1e-4));
  REQUIRE_THAT(D->getVal(1, 1), Catch::Matchers::WithinAbs(4.0, 1e-4));
  delete D;
}

TEST_CASE("test4SymmMatrix") {
  SymmMatrix<double> A(3);
  A.setVal(0, 0, 1.0);
  A.setVal(0, 1, 2.0);
  A.setVal(1, 0, 2.0);
  A.setVal(1, 1, 1.0);
  A.setVal(1, 2, 3.0);
  A.setVal(2, 1, 3.0);
  A.setVal(2, 0, 1.5);
  A.setVal(0, 2, 1.5);
  A.setVal(2, 2, 1.0);

  SymmMatrix<double> B(A);
  SymmMatrix<double> C(3);

  multiply(A, B, C);

  B *= A;
  unsigned int i, j;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      REQUIRE_THAT(B.getVal(i, j),
                   Catch::Matchers::WithinAbs(C.getVal(i, j), 1e-4));
    }
  }

  Vector<double> x(3), y(3);
  x.setVal(0, 1.0);
  x.setVal(1, 2.0);
  x.setVal(2, 3.0);

  multiply(A, x, y);

  REQUIRE_THAT(y.getVal(0), Catch::Matchers::WithinAbs(9.5, 1e-4));
  REQUIRE_THAT(y.getVal(1), Catch::Matchers::WithinAbs(13.0, 1e-4));
  REQUIRE_THAT(y.getVal(2), Catch::Matchers::WithinAbs(10.5, 1e-4));

  REQUIRE(A.getDataSize() == 6);
  A.setToIdentity();
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      if (i != j) {
        REQUIRE_THAT(A.getVal(i, j), Catch::Matchers::WithinAbs(0.0, 1e-6));
      } else {
        REQUIRE_THAT(A.getVal(i, j), Catch::Matchers::WithinAbs(1.0, 1e-6));
      }
    }
  }
}
