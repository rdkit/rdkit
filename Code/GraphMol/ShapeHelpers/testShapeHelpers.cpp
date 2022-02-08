//
//   Copyright (C) 2005-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/test.h>
#include <Geometry/UniformGrid3D.h>
#include "ShapeEncoder.h"
#include "ShapeUtils.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <Geometry/GridUtils.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/Substruct/SubstructMatch.h>

using namespace RDKit;

void test1Encode() {
  RDGeom::UniformGrid3D grd(30.0, 16.0, 10.0);
  std::string rdbase = getenv("RDBASE");
  std::string fname1 =
      rdbase + "/Code/GraphMol/ShapeHelpers/test_data/1oir.mol";
  ROMol *m = MolFileToMol(fname1);
  MolTransforms::canonicalizeMol(*m);

  MolShapes::EncodeShape(*m, grd, 0);
  delete m;

  CHECK_INVARIANT(grd.getOccupancyVect()->getTotalVal() == 7405, "");
}

void test2Compare() {
  RDGeom::UniformGrid3D grd(30.0, 16.0, 10.0);
  std::string rdbase = getenv("RDBASE");
  std::string fname1 =
      rdbase + "/Code/GraphMol/ShapeHelpers/test_data/1oir.mol";
  ROMol *m = MolFileToMol(fname1);
  MolTransforms::canonicalizeMol(*m);

  ROMol *mdup = MolFileToMol(fname1);
  MolTransforms::canonicalizeMol(*mdup);

  double dist = MolShapes::tanimotoDistance(*m, *mdup);
  CHECK_INVARIANT(dist == 0.0, "");
  dist = MolShapes::tverskyIndex(*m, *mdup, 1.0, 1.0);
  CHECK_INVARIANT(dist == 1.0, "");

  delete m;
  delete mdup;

  m = MolFileToMol(fname1);
  std::string fname2 =
      rdbase + "/Code/GraphMol/ShapeHelpers/test_data/1oir_conf.mol";
  ROMol *m2 = MolFileToMol(fname2);
  double rmsd = MolAlign::alignMol(*m, *m2);
  CHECK_INVARIANT(rmsd >= 0.0, "");
  dist = MolShapes::tanimotoDistance(*m, *m2);
  CHECK_INVARIANT(RDKit::feq(dist, 0.31, 0.01), "");
  dist = MolShapes::tverskyIndex(*m, *m2, 1.0, 1.0);
  CHECK_INVARIANT(RDKit::feq(dist, 0.68, 0.01), "");
  delete m2;

  m2 = MolFileToMol(fname2);
  MatchVectType atomMap;
  atomMap.push_back(std::pair<int, int>(18, 27));
  atomMap.push_back(std::pair<int, int>(13, 23));
  atomMap.push_back(std::pair<int, int>(21, 14));
  atomMap.push_back(std::pair<int, int>(24, 7));
  atomMap.push_back(std::pair<int, int>(9, 19));
  atomMap.push_back(std::pair<int, int>(16, 30));
  rmsd = MolAlign::alignMol(*m2, *m, 0, 0, &atomMap);
  dist = MolShapes::tanimotoDistance(*m, *m2);

  CHECK_INVARIANT(RDKit::feq(dist, 0.3593), "");
  delete m;
  delete m2;
}

void testGithub4364() {
  RDGeom::UniformGrid3D grd(30.0, 16.0, 10.0);
  std::string rdbase = getenv("RDBASE");
  std::string fname1 =
      rdbase + "/Code/GraphMol/ShapeHelpers/test_data/1oir.mol";
  std::unique_ptr<RWMol> m(MolFileToMol(fname1));
  MolTransforms::canonicalizeMol(*m);

  std::string fname2 =
      rdbase + "/Code/GraphMol/ShapeHelpers/test_data/1oir_conf.mol";
  std::unique_ptr<RWMol> m2(MolFileToMol(fname2));
  MolTransforms::canonicalizeMol(*m2);
  int cid1 = -1, cid2 = -1;
  auto dist = MolShapes::tanimotoDistance(*m, *m2, cid1, cid2);
  TEST_ASSERT(RDKit::feq(dist, 0.31, 0.01));
  double gridSpacing = 1.0;
  auto dist2 = MolShapes::tanimotoDistance(*m, *m2, cid1, cid2, gridSpacing);
  TEST_ASSERT(fabs(dist - dist2) > 0.001);
}

int main() {
#if 1
  std::cout << "***********************************************************\n";
  std::cout << "Testing ShapeHelpers\n";

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test1Encode \n\n";
  test1Encode();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test2Compare \n\n";
  test2Compare();
#endif
  std::cout << "\t---------------------------------\n";
  std::cout << "\t testGithub4364 \n\n";
  testGithub4364();
  std::cout << "***********************************************************\n";
  return 0;
}
