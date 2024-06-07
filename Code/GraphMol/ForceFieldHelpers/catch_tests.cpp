
//
//  Copyright (C) 2024 Niels Maeder and other RDKit contributors
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
#include <ForceField/MMFF/Params.h>
#include <ForceField/MMFF/BondStretch.h>

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
    auto *contrib = new ForceFields::MMFF::BondStretchContrib(forceField.get(),
                                                              0, 1, &params);
    forceField->contribs().push_back(ForceFields::ContribPtr(contrib));
    CHECK(forceField->minimize() == 0);
    CHECK(std::round(forceField->distance(0, 1)) == 100);
    auto pos1 = mol->getConformer().getAtomPos(0);
    auto pos2 = mol->getConformer().getAtomPos(1);
    auto dist = (pos2 - pos1).length();
    CHECK(std::round(dist) == 100);
  }
}
