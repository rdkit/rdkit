
//
//  Copyright (C) 2024-2025 Niels Maeder and other RDKit contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <cmath>
#include <RDGeneral/test.h>
#include <catch2/catch_all.hpp>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>
#include <ForceField/MMFF/Params.h>
#include <ForceField/MMFF/BondStretch.h>
#include <GraphMol/MolTransforms/MolTransforms.h>

#include "FFConvenience.h"

using namespace RDKit;

TEST_CASE("Test empty force field") {
  auto mol =
      "CCCO |(-1.28533,-0.0567758,0.434662;-0.175447,0.695786,-0.299881;0.918409,-0.342619,-0.555572;1.30936,-0.801512,0.71705)|"_smiles;
  REQUIRE(mol);
  SECTION("basics") {
    auto forceField = ForceFieldsHelper::createEmptyForceFieldForMol(*mol);
    REQUIRE(forceField);
    forceField->initialize();
    CHECK(forceField->minimize() == 0);
    CHECK(forceField->calcEnergy() == 0.0);
    REQUIRE(forceField->numPoints() == mol->getNumAtoms());
    REQUIRE(forceField->positions().size() == mol->getNumAtoms());
    auto dist = forceField->distance(0, 1);
    auto pos1 = mol->getConformer().getAtomPos(0);
    auto pos2 = mol->getConformer().getAtomPos(1);
    CHECK(dist == (pos2 - pos1).length());
  }
  SECTION("add contrib and minimize") {
    auto forceField = ForceFieldsHelper::createEmptyForceFieldForMol(*mol);
    forceField->initialize();
    auto params = ForceFields::MMFF::MMFFBond{6.0, 100.0};
    auto *contrib = new ForceFields::MMFF::BondStretchContrib(forceField.get());
    contrib->addTerm(0, 1, &params);
    forceField->contribs().push_back(ForceFields::ContribPtr(contrib));
    CHECK(forceField->minimize() == 0);
    CHECK(std::round(forceField->distance(0, 1)) == 100);
    auto pos1 = mol->getConformer().getAtomPos(0);
    auto pos2 = mol->getConformer().getAtomPos(1);
    auto dist = (pos2 - pos1).length();
    CHECK(std::round(dist) == 100);
  }
}

TEST_CASE("github #7901") {
#if 1
  auto mb = R"CTAB(Acetone, enolate form
     RDKit          3D

  9  8  0  0  0  0  0  0  0  0999 V2000
   -1.2143   -0.3494    0.0962 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1072    0.3398    0.0041 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1682    1.7287    0.0554 O   0  0  0  0  0  1  0  0  0  0  0  0
    1.2210   -0.3737   -0.1263 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2311   -1.0203    0.9811 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3951   -0.9468   -0.8225 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0281    0.3987    0.2006 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.1862    0.1115   -0.1943 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.8861    0.1115   -0.1943 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  2  4  2  0
  1  5  1  0
  1  6  1  0
  1  7  1  0
  4  8  1  0
  4  9  1  0
M  CHG  1   3  -1
M  END)CTAB";
#else
  auto mb = R"CTAB(Acetone, enolate form
     RDKit          3D

  9  8  0  0  0  0  0  0  0  0999 V2000
   -1.2143   -0.3494    0.0962 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1072    0.3398    0.0041 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1682    1.7287    0.0554 O   0  0  0  0  0  1  0  0  0  0  0  0
    1.2210   -0.3737   -0.1263 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2311   -1.0203    0.9811 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3951   -0.9468   -0.8225 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0281    0.3987    0.2006 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.1862    0.1115   -0.1943 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.8861   -1.3115   -0.1943 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  2  4  2  0
  1  5  1  0
  1  6  1  0
  1  7  1  0
  4  8  1  0
  4  9  1  0
M  CHG  1   3  -1
M  END)CTAB";
#endif
  v2::FileParsers::MolFileParserParams params;
  params.removeHs = false;
  auto mol = v2::FileParsers::MolFromMolBlock(mb, params);
  REQUIRE(mol);
  auto &conf = mol->getConformer();
  std::unique_ptr<ForceFields::ForceField> ff{UFF::constructForceField(*mol)};
  ff->initialize();
  auto needsMore = ff->minimize(200);
  CHECK(!needsMore);
  CHECK((conf.getAtomPos(7) - conf.getAtomPos(8)).length() > 1.80);
  CHECK(MolTransforms::getAngleDeg(conf, 1, 3, 8) > 110);
  CHECK(MolTransforms::getAngleDeg(conf, 1, 3, 7) > 110);
  CHECK(MolTransforms::getAngleDeg(conf, 7, 3, 8) > 110);
}