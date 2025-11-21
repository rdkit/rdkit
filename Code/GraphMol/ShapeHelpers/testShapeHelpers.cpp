//
//   Copyright (C) 2005-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

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

TEST_CASE("Encode") {
  RDGeom::UniformGrid3D grd(30.0, 16.0, 10.0);
  std::string rdbase = getenv("RDBASE");
  std::string fname1 =
      rdbase + "/Code/GraphMol/ShapeHelpers/test_data/1oir.mol";
  auto m = v2::FileParsers::MolFromMolFile(fname1);
  REQUIRE(m);
  MolTransforms::canonicalizeMol(*m);

  MolShapes::EncodeShape(*m, grd, 0);

  CHECK(grd.getOccupancyVect()->getTotalVal() == 7405);
}

TEST_CASE("Compare") {
  std::string rdbase = getenv("RDBASE");
  std::string fname1 =
      rdbase + "/Code/GraphMol/ShapeHelpers/test_data/1oir.mol";
  auto m = v2::FileParsers::MolFromMolFile(fname1);
  REQUIRE(m);
  MolTransforms::canonicalizeMol(*m);

  auto mdup = v2::FileParsers::MolFromMolFile(fname1);
  REQUIRE(mdup);
  MolTransforms::canonicalizeMol(*mdup);

  double dist = MolShapes::tanimotoDistance(*m, *mdup);
  CHECK(dist == 0.0);
  dist = MolShapes::tverskyIndex(*m, *mdup, 1.0, 1.0);
  CHECK(dist == 1.0);

  m = v2::FileParsers::MolFromMolFile(fname1);
  REQUIRE(m);
  std::string fname2 =
      rdbase + "/Code/GraphMol/ShapeHelpers/test_data/1oir_conf.mol";
  auto m2 = v2::FileParsers::MolFromMolFile(fname2);
  REQUIRE(m2);

  double rmsd = MolAlign::alignMol(*m, *m2);
  CHECK(rmsd >= 0.0);
  dist = MolShapes::tanimotoDistance(*m, *m2);
  CHECK_THAT(dist, Catch::Matchers::WithinAbs(0.31, 0.01));
  dist = MolShapes::tverskyIndex(*m, *m2, 1.0, 1.0);
  CHECK_THAT(dist, Catch::Matchers::WithinAbs(0.68, 0.01));

  m2 = v2::FileParsers::MolFromMolFile(fname2);
  MatchVectType atomMap;
  atomMap.push_back(std::pair<int, int>(18, 27));
  atomMap.push_back(std::pair<int, int>(13, 23));
  atomMap.push_back(std::pair<int, int>(21, 14));
  atomMap.push_back(std::pair<int, int>(24, 7));
  atomMap.push_back(std::pair<int, int>(9, 19));
  atomMap.push_back(std::pair<int, int>(16, 30));
  rmsd = MolAlign::alignMol(*m2, *m, 0, 0, &atomMap);
  dist = MolShapes::tanimotoDistance(*m, *m2);

  CHECK_THAT(dist, Catch::Matchers::WithinAbs(0.3561, 0.01));
}

TEST_CASE("GitHub #4364") {
  std::string rdbase = getenv("RDBASE");
  std::string fname1 =
      rdbase + "/Code/GraphMol/ShapeHelpers/test_data/1oir.mol";
  auto m = v2::FileParsers::MolFromMolFile(fname1);
  MolTransforms::canonicalizeMol(*m);

  std::string fname2 =
      rdbase + "/Code/GraphMol/ShapeHelpers/test_data/1oir_conf.mol";
  auto m2 = v2::FileParsers::MolFromMolFile(fname2);
  MolTransforms::canonicalizeMol(*m2);
  int cid1 = -1, cid2 = -1;
  auto dist = MolShapes::tanimotoDistance(*m, *m2, cid1, cid2);
  CHECK_THAT(dist, Catch::Matchers::WithinAbs(0.31, 0.01));
  double gridSpacing = 1.0;
  auto dist2 = MolShapes::tanimotoDistance(*m, *m2, cid1, cid2, gridSpacing);
  CHECK(dist2 - dist > 0.001);
}
