//  $Id$
//
//   Copyright (C) 2016 Sereina Riniker, Paolo Tosco
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/test.h>
#include "Trajectory.h"
#include <Geometry/point.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <ForceField/ForceField.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <RDGeneral/Invariant.h>

using namespace RDKit;

void testSnapshot() {
  {
    boost::shared_array<double> pos;
    Snapshot s(pos);
    bool e = false;
    try {
      s.getPoint2D(12);
    } catch (...) {
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
    } catch (...) {
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
  boost::shared_array<double> c(new double[np * dim]);
  for (unsigned int i = 0; i < np * dim; ++i) {
    c[i] = static_cast<double>(i);
  }
  for (unsigned int i = 0; i < ns; ++i) {
    traj.addSnapshot(Snapshot(c, static_cast<double>(i)));
  }
  TEST_ASSERT(traj.size() == ns);
  {
    bool e = false;
    try {
      traj.getSnapshot(ns);
    } catch (...) {
      e = true;
    }
    TEST_ASSERT(e);
  }
  {
    bool e = false;
    try {
      traj.getSnapshot(0).getPoint2D(np);
    } catch (...) {
      e = true;
    }
    TEST_ASSERT(e);
  }
  for (unsigned int i = 0; i < np; ++i) {
    TEST_ASSERT(RDKit::feq(std::round(traj.getSnapshot(0).getPoint2D(i).x),
                           static_cast<double>(i * dim)));
    TEST_ASSERT(RDKit::feq(std::round(traj.getSnapshot(0).getPoint2D(i).y),
                           static_cast<double>(i * dim + 1)));
    bool e = false;
    try {
      TEST_ASSERT(
          RDKit::feq(std::round(traj.getSnapshot(0).getPoint3D(i).z), 0.0));
    } catch (...) {
      e = true;
    }
    TEST_ASSERT(!e);
  }
  for (unsigned int i = 0; i < ns; ++i) {
    TEST_ASSERT(RDKit::feq(std::round(traj.getSnapshot(i).getEnergy()),
                           static_cast<double>(i)));
  }
  traj.removeSnapshot(0);
  TEST_ASSERT(traj.size() == ns - 1);
  for (unsigned int i = 0; i < ns - 1; ++i) {
    TEST_ASSERT(RDKit::feq(std::round(traj.getSnapshot(i).getEnergy()),
                           static_cast<double>(i + 1)));
  }
  traj.insertSnapshot(0, Snapshot(c, 999.0));
  TEST_ASSERT(traj.size() == ns);
  Snapshot copySnapshot(traj.getSnapshot(0));
  traj.addSnapshot(copySnapshot);
  TEST_ASSERT(traj.size() == ns + 1);
  TEST_ASSERT(RDKit::feq(std::round(traj.getSnapshot(0).getEnergy()), 999.0));
  TEST_ASSERT(RDKit::feq(std::round(traj.getSnapshot(1).getEnergy()), 1.0));
  TEST_ASSERT(RDKit::feq(
      std::round(traj.getSnapshot(traj.size() - 1).getEnergy()), 999.0));
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
  boost::shared_array<double> c(new double[np * dim]);
  for (unsigned int i = 0; i < np * dim; ++i) {
    c[i] = static_cast<double>(i);
  }
  for (unsigned int i = 0; i < ns; ++i) {
    traj.addSnapshot(Snapshot(c, static_cast<double>(i)));
  }
  TEST_ASSERT(traj.size() == ns);
  {
    bool e = false;
    try {
      traj.getSnapshot(ns);
    } catch (...) {
      e = true;
    }
    TEST_ASSERT(e);
  }
  {
    bool e = false;
    try {
      traj.getSnapshot(0).getPoint2D(np);
    } catch (...) {
      e = true;
    }
    TEST_ASSERT(e);
  }
  for (unsigned int i = 0; i < np; ++i) {
    TEST_ASSERT(RDKit::feq(std::round(traj.getSnapshot(0).getPoint3D(i).x),
                           static_cast<double>(i * dim)));
    TEST_ASSERT(RDKit::feq(std::round(traj.getSnapshot(0).getPoint3D(i).y),
                           static_cast<double>(i * dim + 1)));
    TEST_ASSERT(RDKit::feq(std::round(traj.getSnapshot(0).getPoint3D(i).z),
                           static_cast<double>(i * dim + 2)));
    if (!i) {
      bool e = false;
      try {
        traj.getSnapshot(0).getPoint2D(i);
      } catch (...) {
        e = true;
      }
      TEST_ASSERT(e);
    }
  }
  for (unsigned int i = 0; i < ns; ++i) {
    TEST_ASSERT(RDKit::feq(std::round(traj.getSnapshot(i).getEnergy()),
                           static_cast<double>(i)));
  }
  traj.removeSnapshot(0);
  TEST_ASSERT(traj.size() == ns - 1);
  for (unsigned int i = 0; i < ns - 1; ++i) {
    TEST_ASSERT(RDKit::feq(std::round(traj.getSnapshot(i).getEnergy()),
                           static_cast<double>(i + 1)));
  }
  traj.insertSnapshot(0, Snapshot(c, 999.0));
  TEST_ASSERT(traj.size() == ns);
  Snapshot copySnapshot(traj.getSnapshot(0));
  traj.addSnapshot(copySnapshot);
  TEST_ASSERT(traj.size() == ns + 1);
  TEST_ASSERT(RDKit::feq(std::round(traj.getSnapshot(0).getEnergy()), 999.0));
  TEST_ASSERT(RDKit::feq(std::round(traj.getSnapshot(1).getEnergy()), 1.0));
  TEST_ASSERT(RDKit::feq(
      std::round(traj.getSnapshot(traj.size() - 1).getEnergy()), 999.0));
  Trajectory traj2(traj);
  BOOST_LOG(rdErrorLog) << "done" << std::endl;
}

void testReadAmber() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "testReadAmber" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName = rdbase + "/Code/GraphMol/test_data/water_coords_bad.trx";
  {
    Trajectory traj(2, 0);
    bool ok = false;
    try {
      readAmberTrajectory(fName, traj);
    } catch (...) {
      ok = true;
    }
    TEST_ASSERT(ok);
  }
  {
    Trajectory traj(3, 3);
    bool ok = false;
    try {
      readAmberTrajectory(fName, traj);
    } catch (ValueErrorException &e) {
      BOOST_LOG(rdErrorLog) << e.what() << std::endl;
      ok = true;
    }
    TEST_ASSERT(ok);
  }
  fName = rdbase + "/Code/GraphMol/test_data/water_coords_bad2.trx";
  {
    bool ok = false;
    try {
      Trajectory traj(3, 3);
      readAmberTrajectory(fName, traj);
    } catch (ValueErrorException &e) {
      BOOST_LOG(rdErrorLog) << e.what() << std::endl;
      ok = true;
    }
    TEST_ASSERT(ok);
  }
  fName = rdbase + "/Code/GraphMol/test_data/water_coords.trx";
  {
    Trajectory traj(3, 3);
    readAmberTrajectory(fName, traj);
    TEST_ASSERT(traj.size() == 1);
  }
  fName = rdbase + "/Code/GraphMol/test_data/water_coords2.trx";
  {
    Trajectory traj(3, 3);
    readAmberTrajectory(fName, traj);
    TEST_ASSERT(traj.size() == 2);
    Trajectory trajCopy(traj);
  }
  BOOST_LOG(rdErrorLog) << "done" << std::endl;
}

void testReadGromos() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "testReadGromos" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName = rdbase + "/Code/GraphMol/test_data/water_coords_bad.trc";
  {
    Trajectory traj(2, 0);
    bool ok = false;
    try {
      readGromosTrajectory(fName, traj);
    } catch (...) {
      ok = true;
    }
    TEST_ASSERT(ok);
  }
  {
    Trajectory traj(3, 3);
    bool ok = false;
    try {
      readGromosTrajectory(fName, traj);
    } catch (ValueErrorException &e) {
      BOOST_LOG(rdErrorLog) << e.what() << std::endl;
      ok = true;
    }
    TEST_ASSERT(ok);
  }
  fName = rdbase + "/Code/GraphMol/test_data/water_coords_bad2.trc";
  {
    bool ok = false;
    try {
      Trajectory traj(3, 3);
      readGromosTrajectory(fName, traj);
    } catch (ValueErrorException &e) {
      BOOST_LOG(rdErrorLog) << e.what() << std::endl;
      ok = true;
    }
    TEST_ASSERT(ok);
  }
  fName = rdbase + "/Code/GraphMol/test_data/water_coords.trc";
  {
    Trajectory traj(3, 3);
    readGromosTrajectory(fName, traj);
    TEST_ASSERT(traj.size() == 1);
  }
  fName = rdbase + "/Code/GraphMol/test_data/water_coords2.trc";
  {
    Trajectory traj(3, 3);
    readGromosTrajectory(fName, traj);
    TEST_ASSERT(traj.size() == 2);
    Trajectory trajCopy(traj);
  }
  BOOST_LOG(rdErrorLog) << "done" << std::endl;
}

void testAddConformersFromTrajectory() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n";
  BOOST_LOG(rdInfoLog) << "Testing adding conformers from a trajectory"
                       << std::endl;
  std::string molBlock =
      "\n"
      "     RDKit          3D\n"
      "\n"
      " 71 74  0  0  0  0  0  0  0  0999 V2000\n"
      "    8.2543    3.1901   -0.3005 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    7.4558    1.9712    0.0938 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    7.3934    1.0441   -0.9483 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    6.6660   -0.0533   -0.4641 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    5.1928    0.2346   -0.4609 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    4.3713   -0.9410   -0.5770 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    3.1852   -1.0034   -1.2291 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    2.2914    0.1276   -1.6316 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    0.9308   -0.4468   -1.9908 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    0.1417   -0.7821   -0.7545 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -0.1848    0.3695    0.0456 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -1.5661    0.7686   -0.0745 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -2.4768   -0.0640    0.8206 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -3.8874    0.1143    0.3941 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -4.6333   -0.9984    0.0264 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -6.0127   -0.9516   -0.0400 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -6.7062    0.1599    0.3963 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -8.0408    0.4828   -0.1977 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -7.7914    1.1180   -1.5591 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -8.7622    1.4403    0.7265 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -8.8409   -0.7397   -0.4395 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -8.9121   -1.6637    0.4258 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -9.7414   -0.7636   -1.5059 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -5.9736    1.2357    0.8565 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -4.5843    1.2252    0.8530 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    0.6263    1.4884   -0.3942 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    2.0541    1.0258   -0.4230 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    2.9225   -2.3317   -1.2963 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    3.6061   -2.9745   -0.3180 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    3.3554   -4.1536    0.3735 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    3.7653   -4.2712    1.6948 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    4.8254   -3.4613    2.0796 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    5.1978   -2.3436    1.3419 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    4.5694   -2.0799    0.1305 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    9.3138    3.1372    0.0031 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    7.8117    4.0754    0.1798 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    8.2358    3.3535   -1.4074 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    6.4027    2.2146    0.3634 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    7.9270    1.5444    1.0040 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    7.0677   -0.2415    0.5615 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    6.9530   -0.9105   -1.1025 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    4.9578    0.7259    0.5137 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    4.9985    0.9430   -1.3033 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    2.7171    0.7264   -2.4494 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    0.3994    0.2339   -2.6810 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    1.1342   -1.4171   -2.5076 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -0.7632   -1.3370   -1.0391 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    0.7845   -1.4394   -0.1311 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    0.0125    0.1989    1.0673 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -1.6672    1.8215    0.2925 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -1.8705    0.7271   -1.1337 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -2.3045    0.3159    1.8590 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -2.1980   -1.1367    0.7635 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -4.1513   -1.9468   -0.2114 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -6.6138   -1.7460   -0.4718 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -7.0727    0.4399   -2.0858 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -7.3144    2.1076   -1.4482 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -8.7609    1.1720   -2.1135 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -8.3137    2.4504    0.5729 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -8.6170    1.0817    1.7580 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -9.8244    1.4444    0.4200 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -6.4629    2.0541    1.3719 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -4.0445    2.0563    1.3058 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    0.3329    1.8224   -1.3991 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    0.4920    2.3164    0.3160 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    2.2025    0.3766    0.4766 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    2.7945    1.8369   -0.3969 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    2.4404   -4.6964    0.1303 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    3.3157   -5.0055    2.3587 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    5.4272   -3.7654    2.9380 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    5.5668   -1.5069    1.9380 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "  1  2  1  0\n"
      "  2  3  1  0\n"
      "  3  4  1  0\n"
      "  4  5  1  0\n"
      "  5  6  1  0\n"
      "  6  7  1  0\n"
      "  7  8  1  0\n"
      "  8  9  1  0\n"
      "  9 10  1  0\n"
      " 10 11  1  0\n"
      " 11 12  1  0\n"
      " 12 13  1  0\n"
      " 13 14  1  0\n"
      " 14 15  2  0\n"
      " 15 16  1  0\n"
      " 16 17  2  0\n"
      " 17 18  1  0\n"
      " 18 19  1  0\n"
      " 18 20  1  0\n"
      " 18 21  1  0\n"
      " 21 22  2  0\n"
      " 21 23  1  0\n"
      " 17 24  1  0\n"
      " 24 25  2  0\n"
      " 11 26  1  0\n"
      " 26 27  1  0\n"
      "  7 28  2  0\n"
      " 28 29  1  0\n"
      " 29 30  2  0\n"
      " 30 31  1  0\n"
      " 31 32  2  0\n"
      " 32 33  1  0\n"
      " 33 34  2  0\n"
      " 34  6  1  0\n"
      " 27  8  1  0\n"
      " 34 29  1  0\n"
      " 25 14  1  0\n"
      "  1 35  1  0\n"
      "  1 36  1  0\n"
      "  1 37  1  0\n"
      "  2 38  1  0\n"
      "  2 39  1  0\n"
      "  4 40  1  0\n"
      "  4 41  1  0\n"
      "  5 42  1  0\n"
      "  5 43  1  0\n"
      "  8 44  1  0\n"
      "  9 45  1  0\n"
      "  9 46  1  0\n"
      " 10 47  1  0\n"
      " 10 48  1  0\n"
      " 11 49  1  0\n"
      " 12 50  1  0\n"
      " 12 51  1  0\n"
      " 13 52  1  0\n"
      " 13 53  1  0\n"
      " 15 54  1  0\n"
      " 16 55  1  0\n"
      " 19 56  1  0\n"
      " 19 57  1  0\n"
      " 19 58  1  0\n"
      " 20 59  1  0\n"
      " 20 60  1  0\n"
      " 20 61  1  0\n"
      " 24 62  1  0\n"
      " 25 63  1  0\n"
      " 26 64  1  0\n"
      " 26 65  1  0\n"
      " 27 66  1  0\n"
      " 27 67  1  0\n"
      " 30 68  1  0\n"
      " 31 69  1  0\n"
      " 32 70  1  0\n"
      " 33 71  1  0\n"
      "M  CHG  2  11   1  23  -1\n"
      "M  END\n";
  RWMol *mol = MolBlockToMol(molBlock, true, false);
  const unsigned int everySteps = 10;
  const unsigned int maxIts = 1000;
  double gradTol = 0.01;
  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase + "/Code/GraphMol/test_data/bilastine_trajectory.sdf";
  SDWriter w(fName);
  ForceFields::ForceField *field = MMFF::constructForceField(*mol);
  field->initialize();
  auto *sv = new SnapshotVect;
  int res = field->minimize(everySteps, sv, maxIts, gradTol);
  TEST_ASSERT(res == 0);
  Trajectory traj(3, mol->getNumAtoms(), sv);
  mol->removeConformer(0);
  traj.addConformersToMol(*mol);
  for (unsigned int nConf = 0; nConf < mol->getNumConformers(); ++nConf) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(4)
       << traj.getSnapshot(nConf).getEnergy();
    mol->setProp("ENERGY", ss.str(), false);
    w.write(*mol, nConf);
  }
  w.close();
  traj.clear();
  unsigned int n1 = mol->getNumConformers();
  traj.addConformersToMol(*mol);
  unsigned int n2 = mol->getNumConformers();
  TEST_ASSERT(n1 == n2);
  // getSnapshot should raise exception after Clear()
  bool ok = false;
  try {
    traj.getSnapshot(0);
  } catch (Invar::Invariant &e) {
    BOOST_LOG(rdErrorLog) << e.what() << std::endl;
    ok = true;
  }
  TEST_ASSERT(ok);
  delete field;
  delete mol;
}

void testAddConformersFromAmberTrajectory() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "testAddConformersFromAmberTrajectory" << std::endl;

  ROMol *mol = SmilesToMol("CCC");
  std::string rdbase = getenv("RDBASE");
  std::string fName = rdbase + "/Code/GraphMol/test_data/water_coords.trx";
  {
    Trajectory traj(3, mol->getNumAtoms());
    readAmberTrajectory(fName, traj);
    TEST_ASSERT(traj.size() == 1);
    for (unsigned int i = 0; i < 2; ++i) {
      traj.addConformersToMol(*mol);
      TEST_ASSERT(mol->getNumConformers() == i + 1);
      TEST_ASSERT(mol->getConformer(i).getNumAtoms() == 3);
      TEST_ASSERT(RDKit::feq(mol->getConformer(i).getAtomPos(0).x, 0.1941767));
      TEST_ASSERT(RDKit::feq(mol->getConformer(i).getAtomPos(2).z, -0.4088006));
    }
    mol->clearConformers();
    bool e = false;
    try {
      traj.addConformersToMol(*mol, 1);
    } catch (...) {
      e = true;
    }
    TEST_ASSERT(e);
    TEST_ASSERT(mol->getNumConformers() == 0);
  }
  fName = rdbase + "/Code/GraphMol/test_data/water_coords2.trx";
  {
    Trajectory traj(3, mol->getNumAtoms());
    readAmberTrajectory(fName, traj);
    TEST_ASSERT(traj.size() == 2);
    traj.addConformersToMol(*mol);
    TEST_ASSERT(mol->getNumConformers() == 2);
    mol->clearConformers();
    traj.addConformersToMol(*mol, 0, 0);
    TEST_ASSERT(mol->getNumConformers() == 1);
    traj.addConformersToMol(*mol, 1);
    TEST_ASSERT(mol->getNumConformers() == 2);
  }
  delete mol;
  BOOST_LOG(rdErrorLog) << "done" << std::endl;
}

void testAddConformersFromGromosTrajectory() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "testAddConformersFromGromosTrajectory" << std::endl;

  ROMol *mol = SmilesToMol("CCC");
  std::string rdbase = getenv("RDBASE");
  std::string fName = rdbase + "/Code/GraphMol/test_data/water_coords.trc";
  {
    Trajectory traj(3, mol->getNumAtoms());
    readGromosTrajectory(fName, traj);
    TEST_ASSERT(traj.size() == 1);
    for (unsigned int i = 0; i < 2; ++i) {
      traj.addConformersToMol(*mol);
      TEST_ASSERT(mol->getNumConformers() == i + 1);
      TEST_ASSERT(mol->getConformer(i).getNumAtoms() == 3);
      TEST_ASSERT(RDKit::feq(mol->getConformer(i).getAtomPos(0).x, 1.941767));
      TEST_ASSERT(RDKit::feq(mol->getConformer(i).getAtomPos(2).z, -4.088006));
    }
    mol->clearConformers();
    bool e = false;
    try {
      traj.addConformersToMol(*mol, 1);
    } catch (...) {
      e = true;
    }
    TEST_ASSERT(e);
    TEST_ASSERT(mol->getNumConformers() == 0);
  }
  fName = rdbase + "/Code/GraphMol/test_data/water_coords2.trc";
  {
    Trajectory traj(3, mol->getNumAtoms());
    readGromosTrajectory(fName, traj);
    TEST_ASSERT(traj.size() == 2);
    traj.addConformersToMol(*mol);
    TEST_ASSERT(mol->getNumConformers() == 2);
    mol->clearConformers();
    traj.addConformersToMol(*mol, 0, 0);
    TEST_ASSERT(mol->getNumConformers() == 1);
    mol->clearConformers();
    traj.addConformersToMol(*mol, 1);
    TEST_ASSERT(mol->getNumConformers() == 1);
  }
  delete mol;
  BOOST_LOG(rdErrorLog) << "done" << std::endl;
}

int main() {
  BOOST_LOG(rdErrorLog)
      << "***********************************************************\n";
  BOOST_LOG(rdErrorLog) << "Testing Trajectory\n";

  testSnapshot();
  testTrajectory2D();
  testTrajectory3D();
  testReadAmber();
  testReadGromos();
  testAddConformersFromTrajectory();
  testAddConformersFromAmberTrajectory();
  testAddConformersFromGromosTrajectory();

  return 0;
}
