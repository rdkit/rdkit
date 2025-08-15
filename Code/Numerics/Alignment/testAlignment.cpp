//  Copyright (C) 2004-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>
#include "AlignPoints.h"
#include <Numerics/Vector.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <Geometry/Transform3D.h>
#include <Geometry/point.h>

#include <cmath>

using namespace RDNumeric;
using namespace RDNumeric::Alignments;

TEST_CASE("testBasic") {
  RDGeom::Point3DConstPtrVect rpts;
  RDGeom::Point3D rpt1(0.0, 0.0, 0.0);
  rpts.push_back(&rpt1);
  RDGeom::Point3D rpt2(1.0, 0.0, 0.0);
  rpts.push_back(&rpt2);

  RDGeom::Point3DConstPtrVect qpts;
  RDGeom::Point3D qpt1(2.0, 2.0, 0.0);
  qpts.push_back(&qpt1);
  RDGeom::Point3D qpt2(2.0, 3.0, 0.0);
  qpts.push_back(&qpt2);

  RDGeom::Transform3D trans;
  double ssr = AlignPoints(rpts, qpts, trans);
  REQUIRE_THAT(ssr, Catch::Matchers::WithinAbs(0.0, 1e-4));

  // transform qpts and see if we match the rpts
  trans.TransformPoint(qpt1);
  trans.TransformPoint(qpt2);
  qpt1 -= rpt1;
  qpt2 -= rpt2;
  REQUIRE_THAT(qpt1.length(), Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(qpt2.length(), Catch::Matchers::WithinAbs(0.0, 1e-4));
}

TEST_CASE("testTriangle") {
  // try 3 point two equilateral triangles of different edge lengths
  RDGeom::Point3DConstPtrVect rpts;
  RDGeom::Point3D rpt1(-cos(M_PI / 6), -sin(M_PI / 6), 0.0);
  rpts.push_back(&rpt1);
  RDGeom::Point3D rpt2(cos(M_PI / 6), -sin(M_PI / 6), 0.0);
  rpts.push_back(&rpt2);
  RDGeom::Point3D rpt3(0.0, 1.0, 0.0);
  rpts.push_back(&rpt3);

  RDGeom::Point3DConstPtrVect qpts;
  RDGeom::Point3D qpt1(-2 * sin(M_PI / 6) + 3.0, 2 * cos(M_PI / 6), 4.0);
  qpts.push_back(&qpt1);
  RDGeom::Point3D qpt2(-2 * sin(M_PI / 6) + 3.0, -2 * cos(M_PI / 6), 4.0);
  qpts.push_back(&qpt2);
  RDGeom::Point3D qpt3(5.0, 0.0, 4.0);
  qpts.push_back(&qpt3);

  RDGeom::Transform3D trans;
  double ssr = AlignPoints(rpts, qpts, trans);
  REQUIRE_THAT(ssr, Catch::Matchers::WithinAbs(3.0, 1e-4));
  RDGeom::Point3D nqpt1, nqpt2, nqpt3;
  nqpt1 = qpt1;
  nqpt2 = qpt2;
  nqpt3 = qpt3;
  trans.TransformPoint(nqpt1);
  trans.TransformPoint(nqpt2);
  trans.TransformPoint(nqpt3);
  nqpt1 -= rpt1;
  nqpt2 -= rpt2;
  nqpt3 -= rpt3;

  REQUIRE_THAT(nqpt1.length(), Catch::Matchers::WithinAbs(1.0, 1e-4));
  REQUIRE_THAT(nqpt2.length(), Catch::Matchers::WithinAbs(1.0, 1e-4));
  REQUIRE_THAT(nqpt3.length(), Catch::Matchers::WithinAbs(1.0, 1e-4));

  DoubleVector wts(3);
  wts.setVal(0, 1.0);
  wts.setVal(1, 1.0);
  wts.setVal(2, 2.0);
  ssr = AlignPoints(rpts, qpts, trans, &wts);
  REQUIRE_THAT(ssr, Catch::Matchers::WithinAbs(3.75, 1e-4));
  nqpt1 = qpt1;
  nqpt2 = qpt2;
  nqpt3 = qpt3;
  trans.TransformPoint(nqpt1);
  trans.TransformPoint(nqpt2);
  trans.TransformPoint(nqpt3);

  nqpt1 -= rpt1;
  nqpt2 -= rpt2;
  nqpt3 -= rpt3;

  REQUIRE_THAT(nqpt1.length(), Catch::Matchers::WithinAbs(1.14564, 1e-4));
  REQUIRE_THAT(nqpt2.length(), Catch::Matchers::WithinAbs(1.14564, 1e-4));
  REQUIRE_THAT(nqpt3.length(), Catch::Matchers::WithinAbs(0.75, 1e-4));

  wts.setVal(0, 1.0);
  wts.setVal(1, 2.0);
  wts.setVal(2, 2.0);

  ssr = AlignPoints(rpts, qpts, trans, &wts);
  REQUIRE_THAT(ssr, Catch::Matchers::WithinAbs(4.8, 1e-4));
  nqpt1 = qpt1;
  nqpt2 = qpt2;
  nqpt3 = qpt3;
  trans.TransformPoint(nqpt1);
  trans.TransformPoint(nqpt2);
  trans.TransformPoint(nqpt3);
  qpt1 -= rpt1;
  qpt2 -= rpt2;
  qpt3 -= rpt3;
  nqpt1 -= rpt1;
  nqpt2 -= rpt2;
  nqpt3 -= rpt3;
  REQUIRE_THAT(nqpt1.length(), Catch::Matchers::WithinAbs(1.2, 1e-4));
  REQUIRE_THAT(nqpt2.length(), Catch::Matchers::WithinAbs(0.9165, 1e-4));
  REQUIRE_THAT(nqpt3.length(), Catch::Matchers::WithinAbs(0.9165, 1e-4));
}

TEST_CASE("testFourPts") {
  // lets test most point 4 points
  RDGeom::Point3DConstPtrVect rpts;
  RDGeom::Point3D rpt1(0.0, 0.0, 0.0);
  rpts.push_back(&rpt1);
  RDGeom::Point3D rpt2(1.0, 0.0, 0.0);
  rpts.push_back(&rpt2);
  RDGeom::Point3D rpt3(0.0, 1.0, 0.0);
  rpts.push_back(&rpt3);
  RDGeom::Point3D rpt4(0.0, 0.0, 1.0);
  rpts.push_back(&rpt4);

  RDGeom::Point3DConstPtrVect qpts;
  RDGeom::Point3D qpt1(2.0, 2.0, 3.0);
  qpts.push_back(&qpt1);
  RDGeom::Point3D qpt2(3.0, 2.0, 3.0);
  qpts.push_back(&qpt2);
  RDGeom::Point3D qpt3(2.0, 3.0, 3.0);
  qpts.push_back(&qpt3);
  RDGeom::Point3D qpt4(2.0, 2.0, 4.0);
  qpts.push_back(&qpt4);

  RDGeom::Transform3D trans;
  double ssr = AlignPoints(rpts, qpts, trans);
  REQUIRE_THAT(ssr, Catch::Matchers::WithinAbs(0.0, 1e-4));
  trans.TransformPoint(qpt1);
  trans.TransformPoint(qpt2);
  trans.TransformPoint(qpt3);
  trans.TransformPoint(qpt4);
  qpt1 -= rpt1;
  qpt2 -= rpt2;
  qpt3 -= rpt3;
  qpt4 -= rpt4;
  REQUIRE_THAT(qpt1.length(), Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(qpt2.length(), Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(qpt3.length(), Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(qpt4.length(), Catch::Matchers::WithinAbs(0.0, 1e-4));
}

TEST_CASE("testReflection") {
  // lets test most point 4 points
  RDGeom::Point3DConstPtrVect rpts;
  RDGeom::Point3D rpt1(0.0, 0.0, 0.0);
  rpts.push_back(&rpt1);
  RDGeom::Point3D rpt2(1.0, 0.0, 0.0);
  rpts.push_back(&rpt2);
  RDGeom::Point3D rpt3(0.0, 1.0, 0.0);
  rpts.push_back(&rpt3);
  RDGeom::Point3D rpt4(0.0, 0.0, 1.0);
  rpts.push_back(&rpt4);

  RDGeom::Point3DConstPtrVect qpts;
  RDGeom::Point3D qpt1(2.0, 2.0, 3.0);
  qpts.push_back(&qpt1);
  RDGeom::Point3D qpt2(3.0, 2.0, 3.0);
  qpts.push_back(&qpt2);
  RDGeom::Point3D qpt3(2.0, 2.0, 4.0);
  qpts.push_back(&qpt3);
  RDGeom::Point3D qpt4(2.0, 3.0, 3.0);
  qpts.push_back(&qpt4);

  RDGeom::Transform3D trans;
  double ssr = AlignPoints(rpts, qpts, trans);
  REQUIRE_THAT(ssr, Catch::Matchers::WithinAbs(1.0, 1e-4));

  ssr = AlignPoints(rpts, qpts, trans, nullptr, true);
  REQUIRE_THAT(ssr, Catch::Matchers::WithinAbs(0.0, 1e-4));

  trans.TransformPoint(qpt1);
  trans.TransformPoint(qpt2);
  trans.TransformPoint(qpt3);
  trans.TransformPoint(qpt4);
  qpt1 -= rpt1;
  qpt2 -= rpt2;
  qpt3 -= rpt3;
  qpt4 -= rpt4;

  REQUIRE_THAT(qpt1.length(), Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(qpt2.length(), Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(qpt3.length(), Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(qpt4.length(), Catch::Matchers::WithinAbs(0.0, 1e-4));
}
