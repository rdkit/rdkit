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

void testTrajectory2D() {
  const unsigned int dim = 2;
  const unsigned int np = 10;
  const unsigned int ns = 5;
  Trajectory traj(dim, np);
  CHECK_INVARIANT(traj.dimension() == dim, "");
  CHECK_INVARIANT(traj.numPoints() == np, "");
  double *c = new double [np * dim];
  for (unsigned int i = 0; i < np * dim; ++i)
    c[i] = static_cast<double>(i);
  for (unsigned int i = 0; i < ns; ++i)
    traj.addSnapshot(Snapshot(c, static_cast<double>(i)));
  TEST_ASSERT(traj.size() == ns);
  for (unsigned int i = 0; i < np; ++i) {
    TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(0).getPoint2D(i).x), static_cast<double>(i * dim)));
    TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(0).getPoint2D(i).y), static_cast<double>(i * dim + 1)));
    bool e = false;
    try {
      TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(0).getPoint3D(i).z), 0.0));
    }
    catch (...) {
      e = true;
    }
    TEST_ASSERT(!e);
  }
  for (unsigned int i = 0; i < ns; ++i)
    TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(i).getEnergy()), static_cast<double>(i)));
  traj.removeSnapshot(0);
  TEST_ASSERT(traj.size() == ns - 1);
  for (unsigned int i = 0; i < ns - 1; ++i)
    TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(i).getEnergy()), static_cast<double>(i + 1)));
  traj.insertSnapshot(0, Snapshot(c, 999.0));
  TEST_ASSERT(traj.size() == ns);
  TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(0).getEnergy()), 999.0));
  TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(1).getEnergy()), 1.0));
}

void testTrajectory3D() {
  const unsigned int dim = 3;
  const unsigned int np = 10;
  const unsigned int ns = 5;
  Trajectory traj(dim, np);
  CHECK_INVARIANT(traj.dimension() == dim, "");
  CHECK_INVARIANT(traj.numPoints() == np, "");
  double *c = new double [np * dim];
  for (unsigned int i = 0; i < np * dim; ++i)
    c[i] = static_cast<double>(i);
  for (unsigned int i = 0; i < ns; ++i)
    traj.addSnapshot(Snapshot(c, static_cast<double>(i)));
  TEST_ASSERT(traj.size() == ns);
  for (unsigned int i = 0; i < np; ++i) {
    TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(0).getPoint3D(i).x), static_cast<double>(i * dim)));
    TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(0).getPoint3D(i).y), static_cast<double>(i * dim + 1)));
    TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(0).getPoint3D(i).z), static_cast<double>(i * dim + 2)));
    if (!i) {
      bool e = false;
      try {
        traj.getSnapshot(0).getPoint2D(i);
      }
      catch (...) {
        e = true;
      }
      TEST_ASSERT(e);
    }
  }
  for (unsigned int i = 0; i < ns; ++i)
    TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(i).getEnergy()), static_cast<double>(i)));
  traj.removeSnapshot(0);
  TEST_ASSERT(traj.size() == ns - 1);
  for (unsigned int i = 0; i < ns - 1; ++i)
    TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(i).getEnergy()), static_cast<double>(i + 1)));
  traj.insertSnapshot(0, Snapshot(c, 999.0));
  TEST_ASSERT(traj.size() == ns);
  TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(0).getEnergy()), 999.0));
  TEST_ASSERT(RDKit::feq(RDKit::round(traj.getSnapshot(1).getEnergy()), 1.0));
}

int main() {
  std::cout << "***********************************************************\n";
  std::cout << "Testing Trajectory\n";

  std::cout << "\t---------------------------------\n";
  std::cout << "\t testTrajectory2D \n\n";
  testTrajectory2D();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t testTrajectory3D \n\n";
  testTrajectory3D();

  return 0;
}
