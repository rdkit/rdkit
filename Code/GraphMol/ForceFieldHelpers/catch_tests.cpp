
//
//  Copyright (C) 2024 Greg Landrum and other RDKit contributors
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
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <ForceField/MMFF/Params.h>
#include <ForceField/MMFF/BondStretch.h>

#include "FFConvenience.h"

using namespace RDKit;

TEST_CASE("Test empty force field") {
  auto mol = "CCCC"_smiles;
  REQUIRE(mol);
  REQUIRE_FALSE(DGeomHelpers::EmbedMolecule(*mol));
  SECTION("basics") {
    auto forceField = ForceFieldsHelper::constructEmptyForceField(*mol);
    REQUIRE(forceField);
    forceField->initialize();
    TEST_ASSERT(forceField->minimize() == 0);
    TEST_ASSERT(forceField->calcEnergy() == 0.0);
    TEST_ASSERT(forceField->numPoints() == mol->getNumAtoms());
    TEST_ASSERT(forceField->positions().size() == mol->getNumAtoms());
    auto dist = forceField->distance(0, 1);
    TEST_ASSERT(1.2 < dist && dist < 1.6);
    delete forceField;
  }
  SECTION("add contrib and minimize") {
    auto forceField = ForceFieldsHelper::constructEmptyForceField(*mol);
    forceField->initialize();
    auto params = ForceFields::MMFF::MMFFBond{6.0, 100.0};
    auto *contrib =
        new ForceFields::MMFF::BondStretchContrib(forceField, 0, 1, &params);
    forceField->contribs().push_back(ForceFields::ContribPtr(contrib));
    TEST_ASSERT(forceField->minimize() == 0);
    TEST_ASSERT(std::round(forceField->distance(0, 1)) == 100);
    auto pos1 = mol->getConformer().getAtomPos(0);
    auto pos2 = mol->getConformer().getAtomPos(1);
    auto dist = sqrt(pow(pos1.x - pos2.x, 2) + pow(pos1.y - pos2.y, 2) +
                     pow(pos1.z - pos2.z, 2));
    TEST_ASSERT(99.0 < dist && dist < 101.0);
    delete forceField;
  }
}
