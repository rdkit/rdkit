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
#include "Matrix.h"
#include "SquareMatrix.h"
#include "SymmMatrix.h"
#include <iostream>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <math.h>
#include <RDGeneral/utils.h>

using namespace RDNumeric;


void test1Vector() {
  Vector<double> v1(3, 1.0);
  v1.setVal(0, 2.0); v1.setVal(2, -4.0);
  CHECK_INVARIANT(RDKit::feq(v1.normL1(), 7.0), "");
  CHECK_INVARIANT(RDKit::feq(v1.normLinfinity(), 4.0), "");
  CHECK_INVARIANT(RDKit::feq(v1.normL2(), sqrt(21.0)), "");

  
  v1.setVal(1, 2.0);
  CHECK_INVARIANT(RDKit::feq(v1.getVal(1), 2.0), "");
  CHECK_INVARIANT(RDKit::feq(v1.normL1(), 8.0), "");

  double *data = new double[3];
  data[0] = 1.0;
  data[1] = 2.0;
  data[2] = 3.0;
  Vector<double>::DATA_SPTR sdata(data);
  Vector<double> *v2 = new Vector<double>(3,sdata);
  CHECK_INVARIANT(RDKit::feq(v2->normL1(), 6.0), "");
  Vector<double> v3(v1);
  unsigned int i;
  for (i = 0; i < v1.size(); i++) {
    CHECK_INVARIANT(RDKit::feq(v1.getVal(i), v3.getVal(i)), "");
  }

  delete v2;
  //delete [] data;

  Vector<double> vr1(5);
  Vector<double> vr2(5);
  vr1.setToRandom();
  CHECK_INVARIANT(RDKit::feq(vr1.normL2(), 1.0), "");
  vr2.setToRandom(120);
  CHECK_INVARIANT(RDKit::feq(vr2.normL2(), 1.0), "");
}

void test2Matrix() {
  
  Matrix<double> A(2,3);
  A.setVal(0,0, 1.0); A.setVal(0,1, 0.5); A.setVal(0,2, 2.0);
  A.setVal(1,0, 0.5); A.setVal(1,1, 1.0); A.setVal(1,2, 3.0);

  Vector<double> v1(3, 1.0);
  v1.setVal(1, 2.0); v1.setVal(2, 3.0);
  
  Vector<double> v2(2);
  multiply(A, v1, v2);
  CHECK_INVARIANT(RDKit::feq(v2.getVal(0), 8.0), "");
  CHECK_INVARIANT(RDKit::feq(v2.getVal(1), 11.5), "");

  double *data = A.getData();
  data[2] = 3.0;
  CHECK_INVARIANT(RDKit::feq(A.getVal(0,2), 3.0), "");
  multiply(A, v1, v2);
  CHECK_INVARIANT(RDKit::feq(v2.getVal(0), 11.0), "");
  CHECK_INVARIANT(RDKit::feq(v2.getVal(1), 11.5), "");

  data = new double[6];
  Matrix<double> *B = new Matrix<double>(2,3, Matrix<double>::DATA_SPTR(data));
  Matrix<double>E(A); //(*B) = A;
  multiply(E, v1, v2);
  CHECK_INVARIANT(RDKit::feq(v2.getVal(0), 11.0), "");
  CHECK_INVARIANT(RDKit::feq(v2.getVal(1), 11.5), "");
  delete B;
  //delete [] data;

  Matrix<double> D(3,2);
  A.transpose(D);
  Matrix<double> C(2,2);
  multiply(A, D, C);
  
  CHECK_INVARIANT(RDKit::feq(C.getVal(0,0), 10.25), "");
  CHECK_INVARIANT(RDKit::feq(C.getVal(0,1), 10), "");
  CHECK_INVARIANT(RDKit::feq(C.getVal(1,0), 10), "");
  CHECK_INVARIANT(RDKit::feq(C.getVal(1,1), 10.25), "");
}

void test3SquareMatrix() {
  SquareMatrix<double> A(2);
  A.setVal(0,0, 1.0); A.setVal(0,1, 2.0);
  A.setVal(1,0, 3.0); A.setVal(1,1, 4.0);
  SquareMatrix<double> B(A), C(2);

  multiply(A, B,C);
  CHECK_INVARIANT(RDKit::feq(C.getVal(0,0), 7.0), "");
  CHECK_INVARIANT(RDKit::feq(C.getVal(0,1), 10.0), "");
  CHECK_INVARIANT(RDKit::feq(C.getVal(1,0), 15.0), "");
  CHECK_INVARIANT(RDKit::feq(C.getVal(1,1), 22.0), "");
  B *= A;
  CHECK_INVARIANT(RDKit::feq(B.getVal(0,0), 7.0), "");
  CHECK_INVARIANT(RDKit::feq(B.getVal(0,1), 10.0), "");
  CHECK_INVARIANT(RDKit::feq(B.getVal(1,0), 15.0), "");
  CHECK_INVARIANT(RDKit::feq(B.getVal(1,1), 22.0), "");

  double *data = new double[4];
  data[0] = 1.0; data[1] = 2.0;
  data[2] = 3.0; data[3] = 4.0;
  SquareMatrix<double>* D = new SquareMatrix<double>(2, Matrix<double>::DATA_SPTR(data));
  SquareMatrix<double> E(*D);
  multiply((*D), E, A);
  
  unsigned int i, j;
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      CHECK_INVARIANT(RDKit::feq(B.getVal(i,j), A.getVal(i,j)), "");
    }
  }
  D->transposeInplace();
  CHECK_INVARIANT(RDKit::feq(D->getVal(0,0), 1.0), "");
  CHECK_INVARIANT(RDKit::feq(D->getVal(0,1), 3.0), "");
  CHECK_INVARIANT(RDKit::feq(D->getVal(1,0), 2.0), "");
  CHECK_INVARIANT(RDKit::feq(D->getVal(1,1), 4.0), "");
  delete D;
  //delete [] data;
}

void test4SymmMatrix() {
  SymmMatrix<double> A(3);
  A.setVal(0,0, 1.0); A.setVal(0,1, 2.0); A.setVal(1,0, 2.0);
  A.setVal(1,1, 1.0); A.setVal(1,2, 3.0); A.setVal(2,1, 3.0);
  A.setVal(2,0, 1.5); A.setVal(0,2, 1.5); A.setVal(2,2, 1.0);
  
  SymmMatrix<double> B(A);
  SymmMatrix<double> C(3);
  
  multiply(A, B, C);

  B *= A;
  unsigned int i, j;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      CHECK_INVARIANT(RDKit::feq(B.getVal(i,j), C.getVal(i,j)), "");
    }
  }

  Vector<double> x(3), y(3);
  x.setVal(0, 1.0); x.setVal(1, 2.0); x.setVal(2, 3.0);
  
  multiply(A, x, y);
  
  CHECK_INVARIANT(RDKit::feq(y.getVal(0), 9.5), "");
  CHECK_INVARIANT(RDKit::feq(y.getVal(1), 13.0), "");
  CHECK_INVARIANT(RDKit::feq(y.getVal(2), 10.5), "");
  
  CHECK_INVARIANT(A.getDataSize() == 6, "");
  A.setToIdentity();
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      if (i != j) {
        CHECK_INVARIANT(RDKit::feq(A.getVal(i,j), 0.0, 0.000001), "");
      } else {
        CHECK_INVARIANT(RDKit::feq(A.getVal(i,j), 1.0, 0.000001), "");
      }
    }
  }
    
}

int main() {
  RDLog::InitLogs();

  BOOST_LOG(rdErrorLog) <<"-----------------------------------------\n";
  BOOST_LOG(rdErrorLog) <<"Testing RDNumerics: vectors and matrices Code\n";
  
  BOOST_LOG(rdErrorLog) <<"---------------------------------------\n";
  BOOST_LOG(rdErrorLog) <<"\t test1Vector\n";
  test1Vector();
  
  BOOST_LOG(rdErrorLog) <<"---------------------------------------\n";
  BOOST_LOG(rdErrorLog) <<"\t test2Matrix\n";
  test2Matrix();

  BOOST_LOG(rdErrorLog) <<"---------------------------------------\n";
  BOOST_LOG(rdErrorLog) <<"\t test3SquareMatrix\n";
  test3SquareMatrix();

  BOOST_LOG(rdErrorLog) <<"---------------------------------------\n";
  BOOST_LOG(rdErrorLog) <<"\t test4SymmMatrix\n";
  test4SymmMatrix();
  return 0;
}

