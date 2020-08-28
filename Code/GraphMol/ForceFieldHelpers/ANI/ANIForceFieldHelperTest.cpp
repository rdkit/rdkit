//
//  Copyright (C) 2020 Manan Goel
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "RDGeneral/test.h"
#include "catch.hpp"

#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/FileParsers.h>

#include <ForceField/ANI/AtomicContrib.h>
#include <ForceField/ForceField.h>
#include <GraphMol/Descriptors/AtomicEnvironmentVector.h>
#include <Eigen/Dense>
#include <GraphMol/ForceFieldHelpers/ANI/Builder.h>
#include "ANI.h"

using namespace RDKit;

TEST_CASE("Check ANI Force Field builder") {
  SECTION("ANI-1ccx") {
    std::string pathName = getenv("RDBASE");
    std::string filePath =
        pathName + "/Code/GraphMol/Descriptors/test_data/CH4.mol";

    auto mol = MolFileToMol(filePath, true, false);
    REQUIRE(mol);
    int confId = -1;
    std::unique_ptr<ForceFields::ForceField> field(RDKit::ANI::constructForceField(*mol, "ANI-1ccx", 8));
    field->initialize();
    auto numAtoms = mol->getNumAtoms();
    double *pos;
    auto conf = mol->getConformer(confId);
    pos = new double[numAtoms * 3];
    for (unsigned int i = 0; i < mol->getNumAtoms(); i++) {
      auto atom = conf.getAtomPos(i);
      pos[3 * i] = atom.x;
      pos[3 * i + 1] = atom.y;
      pos[3 * i + 2] = atom.z;
    }
    CHECK(std::fabs(field->calcEnergy(pos) - (-40.0553)) < 0.05);
    field->minimize();
    CHECK(field->calcEnergy() - (-40.0553) < 0);

    delete[] pos;
  }
  SECTION("ANI-1x") {
    std::string pathName = getenv("RDBASE");
    std::string filePath =
        pathName + "/Code/GraphMol/Descriptors/test_data/CH4.mol";

    auto mol = MolFileToMol(filePath, true, false);
    REQUIRE(mol);
    int confId = -1;
    std::unique_ptr<ForceFields::ForceField> field(RDKit::ANI::constructForceField(*mol, "ANI-1ccx", 8));
    field->initialize();
    auto numAtoms = mol->getNumAtoms();
    double *pos;
    auto conf = mol->getConformer(confId);
    pos = new double[numAtoms * 3];
    for (unsigned int i = 0; i < mol->getNumAtoms(); i++) {
      auto atom = conf.getAtomPos(i);
      pos[3 * i] = atom.x;
      pos[3 * i + 1] = atom.y;
      pos[3 * i + 2] = atom.z;
    }
    CHECK(std::fabs(field->calcEnergy(pos) - (-40.0517)) < 0.05);
    field->minimize();
    CHECK(field->calcEnergy() - (-40.0517) < 0);
    delete[] pos;
  }
  SECTION("ANI Convenience Test") {
    std::string pathName = getenv("RDBASE");
    std::string filePath =
        pathName + "/Code/GraphMol/Descriptors/test_data/CH4.mol";

    auto mol = MolFileToMol(filePath, true, false);
    REQUIRE(mol);
    auto res = RDKit::ANI::ANIOptimizeMolecule(*mol, "ANI-1x", 8);
    CHECK(res.first == 0);
    CHECK(res.second < -40.0517);
    std::vector<std::pair<int, double>> res1;
    RDKit::ANI::ANIOptimizeMoleculeConfs(*mol, res1, "ANI-1x", 8);

    CHECK(res1.size() == 1);
    CHECK(res1[0].second < (-40.0517));
  }
  SECTION("Unsupported Atoms") {
    std::string pathName = getenv("RDBASE");
    std::string filePath =
        pathName + "/Code/GraphMol/Descriptors/test_data/SO2.mol";

    auto mol = MolFileToMol(filePath, true, false);
    REQUIRE(mol);
    REQUIRE_THROWS_AS(RDKit::ANI::constructForceField(*mol, "ANI-1x", 8),
                      ValueErrorException);
  }
}