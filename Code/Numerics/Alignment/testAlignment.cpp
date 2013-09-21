// $Id$
//  Copyright (C) 2004-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "AlignPoints.h"
#include <Numerics/Vector.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <Geometry/Transform3D.h>
#include <Geometry/point.h>

#include <math.h>

using namespace RDNumeric;
using namespace RDNumeric::Alignments;

void testBasic() {
  RDGeom::Point3DConstPtrVect rpts;
  RDGeom::Point3D rpt1(0.0, 0.0, 0.0); rpts.push_back(&rpt1);
  RDGeom::Point3D rpt2(1.0, 0.0, 0.0); rpts.push_back(&rpt2);

  RDGeom::Point3DConstPtrVect qpts;
  RDGeom::Point3D qpt1(2.0, 2.0, 0.0); qpts.push_back(&qpt1);
  RDGeom::Point3D qpt2(2.0, 3.0, 0.0); qpts.push_back(&qpt2);
  
  RDGeom::Transform3D trans;
  double ssr = AlignPoints(rpts, qpts, trans);
  CHECK_INVARIANT(RDKit::feq(ssr, 0.0), "");

  // transform qpts and see if we match the rpts
  trans.TransformPoint(qpt1); trans.TransformPoint(qpt2);
  qpt1 -= rpt1;
  qpt2 -= rpt2;
  CHECK_INVARIANT(RDKit::feq(qpt1.length(), 0.0), "");
  CHECK_INVARIANT(RDKit::feq(qpt2.length(), 0.0), "");
}

void testTraingle() {
  // try 3 point two equilateral triangles of different edge lengths
  RDGeom::Point3DConstPtrVect rpts;
  RDGeom::Point3D rpt1(-cos(M_PI/6), -sin(M_PI/6), 0.0); rpts.push_back(&rpt1);
  RDGeom::Point3D rpt2(cos(M_PI/6), -sin(M_PI/6), 0.0); rpts.push_back(&rpt2);
  RDGeom::Point3D rpt3(0.0, 1.0, 0.0); rpts.push_back(&rpt3);

  RDGeom::Point3DConstPtrVect qpts;
  RDGeom::Point3D qpt1(-2*sin(M_PI/6) + 3.0, 2*cos(M_PI/6), 4.0); qpts.push_back(&qpt1);
  RDGeom::Point3D qpt2(-2*sin(M_PI/6) + 3.0, -2*cos(M_PI/6), 4.0); qpts.push_back(&qpt2);
  RDGeom::Point3D qpt3(5.0, 0.0, 4.0); qpts.push_back(&qpt3);
  
  RDGeom::Transform3D trans;
  double ssr = AlignPoints(rpts, qpts, trans);
  CHECK_INVARIANT(RDKit::feq(ssr, 3.0), "");
  RDGeom::Point3D nqpt1, nqpt2, nqpt3;
  nqpt1 = qpt1;  nqpt2 = qpt2; nqpt3 = qpt3;
  trans.TransformPoint(nqpt1); 
  trans.TransformPoint(nqpt2);
  trans.TransformPoint(nqpt3);
  nqpt1 -= rpt1; nqpt2 -= rpt2; nqpt3 -= rpt3;
  
  CHECK_INVARIANT(RDKit::feq(nqpt1.length(), 1.0), "");
  CHECK_INVARIANT(RDKit::feq(nqpt2.length(), 1.0), "");
  CHECK_INVARIANT(RDKit::feq(nqpt3.length(), 1.0), "");

  DoubleVector wts(3);
  wts.setVal(0, 1.0); wts.setVal(1, 1.0); wts.setVal(2, 2.0);
  ssr = AlignPoints(rpts, qpts, trans, &wts);
  CHECK_INVARIANT(RDKit::feq(ssr, 3.75), "");
  nqpt1 = qpt1;  nqpt2 = qpt2; nqpt3 = qpt3;
  trans.TransformPoint(nqpt1); 
  trans.TransformPoint(nqpt2);
  trans.TransformPoint(nqpt3);
  
  nqpt1 -= rpt1; nqpt2 -= rpt2; nqpt3 -= rpt3;
  
  CHECK_INVARIANT(RDKit::feq(nqpt1.length(), 1.14564), "");
  CHECK_INVARIANT(RDKit::feq(nqpt2.length(), 1.14564), "");
  CHECK_INVARIANT(RDKit::feq(nqpt3.length(), 0.75), "");

  wts.setVal(0, 1.0); wts.setVal(1, 2.0); wts.setVal(2, 2.0);
 
  ssr = AlignPoints(rpts, qpts, trans, &wts);
  CHECK_INVARIANT(RDKit::feq(ssr, 4.8), "");
  nqpt1 = qpt1;  nqpt2 = qpt2; nqpt3 = qpt3;
  trans.TransformPoint(nqpt1); 
  trans.TransformPoint(nqpt2);
  trans.TransformPoint(nqpt3);
  qpt1 -= rpt1; qpt2 -= rpt2; qpt3 -= rpt3;
  nqpt1 -= rpt1; nqpt2 -= rpt2; nqpt3 -= rpt3;
  CHECK_INVARIANT(RDKit::feq(nqpt1.length(), 1.2), "");
  CHECK_INVARIANT(RDKit::feq(nqpt2.length(), 0.9165), "");
  CHECK_INVARIANT(RDKit::feq(nqpt3.length(), 0.9165), "");
}

void testFourPts() {
  // lets test most point 4 points 
  RDGeom::Point3DConstPtrVect rpts;
  RDGeom::Point3D rpt1(0.0, 0.0, 0.0); rpts.push_back(&rpt1);
  RDGeom::Point3D rpt2(1.0, 0.0, 0.0); rpts.push_back(&rpt2);
  RDGeom::Point3D rpt3(0.0, 1.0, 0.0); rpts.push_back(&rpt3);
  RDGeom::Point3D rpt4(0.0, 0.0, 1.0); rpts.push_back(&rpt4);

  RDGeom::Point3DConstPtrVect qpts;
  RDGeom::Point3D qpt1(2.0, 2.0, 3.0); qpts.push_back(&qpt1);
  RDGeom::Point3D qpt2(3.0, 2.0, 3.0); qpts.push_back(&qpt2);
  RDGeom::Point3D qpt3(2.0, 3.0, 3.0); qpts.push_back(&qpt3);
  RDGeom::Point3D qpt4(2.0, 2.0, 4.0); qpts.push_back(&qpt4);
  
  RDGeom::Transform3D trans;
  double ssr = AlignPoints(rpts, qpts, trans);
  CHECK_INVARIANT(RDKit::feq(ssr, 0.0), "");
  trans.TransformPoint(qpt1); 
  trans.TransformPoint(qpt2);
  trans.TransformPoint(qpt3);
  trans.TransformPoint(qpt4);
  qpt1 -= rpt1; qpt2 -= rpt2; qpt3 -= rpt3; qpt4 -= rpt4;
  CHECK_INVARIANT(RDKit::feq(qpt1.length(), 0.0), "");
  CHECK_INVARIANT(RDKit::feq(qpt2.length(), 0.0), "");
  CHECK_INVARIANT(RDKit::feq(qpt3.length(), 0.0), "");
  CHECK_INVARIANT(RDKit::feq(qpt4.length(), 0.0), "");
}

void testReflection() {
  // lets test most point 4 points 
  RDGeom::Point3DConstPtrVect rpts;
  RDGeom::Point3D rpt1(0.0, 0.0, 0.0); rpts.push_back(&rpt1);
  RDGeom::Point3D rpt2(1.0, 0.0, 0.0); rpts.push_back(&rpt2);
  RDGeom::Point3D rpt3(0.0, 1.0, 0.0); rpts.push_back(&rpt3);
  RDGeom::Point3D rpt4(0.0, 0.0, 1.0); rpts.push_back(&rpt4);

  RDGeom::Point3DConstPtrVect qpts;
  RDGeom::Point3D qpt1(2.0, 2.0, 3.0); qpts.push_back(&qpt1);
  RDGeom::Point3D qpt2(3.0, 2.0, 3.0); qpts.push_back(&qpt2);
  RDGeom::Point3D qpt3(2.0, 2.0, 4.0); qpts.push_back(&qpt3);
  RDGeom::Point3D qpt4(2.0, 3.0, 3.0); qpts.push_back(&qpt4);

  RDGeom::Transform3D trans;
  double ssr = AlignPoints(rpts, qpts, trans);
  CHECK_INVARIANT(RDKit::feq(ssr, 1.0), "");
  
  ssr = AlignPoints(rpts, qpts, trans, 0, true);
  CHECK_INVARIANT(RDKit::feq(ssr, 0.0), "");
  
  trans.TransformPoint(qpt1); 
  trans.TransformPoint(qpt2);
  trans.TransformPoint(qpt3);
  trans.TransformPoint(qpt4);
  qpt1 -= rpt1; qpt2 -= rpt2; qpt3 -= rpt3; qpt4 -= rpt4;
  
  CHECK_INVARIANT(RDKit::feq(qpt1.length(), 0.0), "");
  CHECK_INVARIANT(RDKit::feq(qpt2.length(), 0.0), "");
  CHECK_INVARIANT(RDKit::feq(qpt3.length(), 0.0), "");
  CHECK_INVARIANT(RDKit::feq(qpt4.length(), 0.0), "");
  
}

int main() {
  std::cout << "-----------------------------------------\n";
  std::cout << "Testing Alignment Code\n";

  std::cout << "---------------------------------------\n";
  std::cout << "\t testBasic\n";
  testBasic();

  std::cout << "---------------------------------------\n";
  std::cout << "\t testTriangle\n";
  testTraingle();
  
  std::cout << "---------------------------------------\n";
  std::cout << "\t testFourPts\n";
  testFourPts();

  std::cout << "---------------------------------------\n";
  std::cout << "\t testReflection\n";
  testReflection();
  std::cout << "---------------------------------------\n";
  return (0);
}

