//  $Id$
//
//   Copyright (C) 2005-2013 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Trajectory.h"
#include "point.h"
#include <RDGeneral/Invariant.h>

using namespace RDGeom;
using namespace RDKit;

void testSnapshot() {
  {
    boost::shared_array<double> pos;
    Snapshot s(pos);
    bool e = false;
    try {
      s.getPoint2D(12);
    }
    catch (...) {
      e = true;
    }
    TEST_ASSERT(e);
  }
  {
    boost::shared_array<double> pos(new double[3]());
    Snapshot s(pos);
    bool e = false;
    try {
      s.getPoint2D(0);
    }
    catch (...) {
      e = true;
    }
    TEST_ASSERT(e);
  }
}

void testTrajectory2D() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "testTrajectory2D" << std::endl;
  const unsigned int dim = 2;
  const unsigned int np = 10;
  const unsigned int ns = 5;
  Trajectory traj(dim, np);
  CHECK_INVARIANT(traj.dimension() == dim, "");
  CHECK_INVARIANT(traj.numPoints() == np, "");
  boost::shared_array<double> c(new double [np * dim]);
  for (unsigned int i = 0; i < np * dim; ++i)
    c[i] = static_cast<double>(i);
  for (unsigned int i = 0; i < ns; ++i)
    traj.addSnapshot(new Snapshot(c, static_cast<double>(i)));
  TEST_ASSERT(traj.size() == ns);
  {
    bool e = false;
    try {
      traj.getSnapshot(ns);
    }
    catch (...) {
      e = true;
    }
    TEST_ASSERT(e);
  }
  {
    bool e = false;
    try {
      traj.getSnapshot(0)->getPoint2D(np);
    }
    catch (...) {
      e = true;
    }
    TEST_ASSERT(e);
  }
  for (unsigned int i = 0; i < np; ++i) {
    TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(0)->getPoint2D(i).x), static_cast<double>(i * dim)));
    TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(0)->getPoint2D(i).y), static_cast<double>(i * dim + 1)));
    bool e = false;
    try {
      TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(0)->getPoint3D(i).z), 0.0));
    }
    catch (...) {
      e = true;
    }
    TEST_ASSERT(!e);
  }
  for (unsigned int i = 0; i < ns; ++i)
    TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(i)->getEnergy()), static_cast<double>(i)));
  traj.removeSnapshot(0);
  TEST_ASSERT(traj.size() == ns - 1);
  for (unsigned int i = 0; i < ns - 1; ++i)
    TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(i)->getEnergy()), static_cast<double>(i + 1)));
  traj.insertSnapshot(0, new Snapshot(c, 999.0));
  TEST_ASSERT(traj.size() == ns);
  Snapshot *copySnapshot = new Snapshot(*(traj.getSnapshot(0)));
  traj.addSnapshot(copySnapshot);
  TEST_ASSERT(traj.size() == ns + 1);
  TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(0)->getEnergy()), 999.0));
  TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(1)->getEnergy()), 1.0));
  TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(traj.size() - 1)->getEnergy()), 999.0));
  Trajectory traj2(traj);
  BOOST_LOG(rdErrorLog) << "done" << std::endl;
}

void testTrajectory3D() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "testTrajectory3D" << std::endl;
  const unsigned int dim = 3;
  const unsigned int np = 10;
  const unsigned int ns = 5;
  Trajectory traj(dim, np);
  CHECK_INVARIANT(traj.dimension() == dim, "");
  CHECK_INVARIANT(traj.numPoints() == np, "");
  boost::shared_array<double> c(new double [np * dim]);
  for (unsigned int i = 0; i < np * dim; ++i)
    c[i] = static_cast<double>(i);
  for (unsigned int i = 0; i < ns; ++i)
    traj.addSnapshot(new Snapshot(c, static_cast<double>(i)));
  TEST_ASSERT(traj.size() == ns);
  {
    bool e = false;
    try {
      traj.getSnapshot(ns);
    }
    catch (...) {
      e = true;
    }
    TEST_ASSERT(e);
  }
  {
    bool e = false;
    try {
      traj.getSnapshot(0)->getPoint2D(np);
    }
    catch (...) {
      e = true;
    }
    TEST_ASSERT(e);
  }
  for (unsigned int i = 0; i < np; ++i) {
    TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(0)->getPoint3D(i).x), static_cast<double>(i * dim)));
    TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(0)->getPoint3D(i).y), static_cast<double>(i * dim + 1)));
    TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(0)->getPoint3D(i).z), static_cast<double>(i * dim + 2)));
    if (!i) {
      bool e = false;
      try {
        traj.getSnapshot(0)->getPoint2D(i);
      }
      catch (...) {
        e = true;
      }
      TEST_ASSERT(e);
    }
  }
  for (unsigned int i = 0; i < ns; ++i)
    TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(i)->getEnergy()), static_cast<double>(i)));
  traj.removeSnapshot(0);
  TEST_ASSERT(traj.size() == ns - 1);
  for (unsigned int i = 0; i < ns - 1; ++i)
    TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(i)->getEnergy()), static_cast<double>(i + 1)));
  traj.insertSnapshot(0, new Snapshot(c, 999.0));
  TEST_ASSERT(traj.size() == ns);
  Snapshot *copySnapshot = new Snapshot(*(traj.getSnapshot(0)));
  traj.addSnapshot(copySnapshot);
  TEST_ASSERT(traj.size() == ns + 1);
  TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(0)->getEnergy()), 999.0));
  TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(1)->getEnergy()), 1.0));
  TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(traj.size() - 1)->getEnergy()), 999.0));
  Trajectory traj2(traj);
  BOOST_LOG(rdErrorLog) << "done" << std::endl;
}

void testReadAmber() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "testReadAmber" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName = rdbase + "/Code/Geometry/testData/water_coords_bad.trx";
  {
    Trajectory traj(2, 0);
    bool ok = false;
    try {
      traj.readAmber(fName);
    }
    catch (...) {
      ok = true;
    }
    TEST_ASSERT(ok);
  }
  {
    Trajectory traj(3, 3);
    bool ok = false;
    try {
      traj.readAmber(fName);
    }
    catch (ValueErrorException &e) {
      BOOST_LOG(rdErrorLog) << e.message() << std::endl;
      ok = true;
    }
    TEST_ASSERT(ok);
  }
  fName = rdbase + "/Code/Geometry/testData/water_coords_bad2.trx";
  {
    bool ok = false;
    try {
      Trajectory traj(3, 3);
      traj.readAmber(fName);
    }
    catch (ValueErrorException &e) {
      BOOST_LOG(rdErrorLog) << e.message() << std::endl;
      ok = true;
    }
    TEST_ASSERT(ok);
  }
  fName = rdbase + "/Code/Geometry/testData/water_coords.trx";
  {
    Trajectory traj(3, 3);
    traj.readAmber(fName);
    TEST_ASSERT(traj.size() == 1);
  }
  fName = rdbase + "/Code/Geometry/testData/water_coords2.trx";
  {
    Trajectory traj(3, 3);
    traj.readAmber(fName);
    TEST_ASSERT(traj.size() == 2);
    Trajectory trajCopy(traj);
  }
  BOOST_LOG(rdErrorLog) << "done" << std::endl;
}

void testReadGromos() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "testReadGromos" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName = rdbase + "/Code/Geometry/testData/water_coords_bad.trc";
  {
    Trajectory traj(2, 0);
    bool ok = false;
    try {
      traj.readGromos(fName);
    }
    catch (...) {
      ok = true;
    }
    TEST_ASSERT(ok);
  }
  {
    Trajectory traj(3, 3);
    bool ok = false;
    try {
      traj.readGromos(fName);
    }
    catch (ValueErrorException &e) {
      BOOST_LOG(rdErrorLog) << e.message() << std::endl;
      ok = true;
    }
    TEST_ASSERT(ok);
  }
  fName = rdbase + "/Code/Geometry/testData/water_coords_bad2.trc";
  {
    bool ok = false;
    try {
      Trajectory traj(3, 3);
      traj.readGromos(fName);
    }
    catch (ValueErrorException &e) {
      BOOST_LOG(rdErrorLog) << e.message() << std::endl;
      ok = true;
    }
    TEST_ASSERT(ok);
  }
  fName = rdbase + "/Code/Geometry/testData/water_coords.trc";
  {
    Trajectory traj(3, 3);
    traj.readGromos(fName);
    TEST_ASSERT(traj.size() == 1);
  }
  fName = rdbase + "/Code/Geometry/testData/water_coords2.trc";
  {
    Trajectory traj(3, 3);
    traj.readGromos(fName);
    TEST_ASSERT(traj.size() == 2);
    Trajectory trajCopy(traj);
  }
  BOOST_LOG(rdErrorLog) << "done" << std::endl;
}

int main() {
  BOOST_LOG(rdErrorLog) << "***********************************************************\n";
  BOOST_LOG(rdErrorLog) << "Testing Trajectory\n";

  testSnapshot();
  testTrajectory2D();
  testTrajectory3D();
  testReadAmber();
  testReadGromos();

  return 0;
}
