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

using namespace RDKit;

TEST_CASE("Check ANI Force Field builder") {
  SECTION("ANI-1ccx") {
    std::string pathName = getenv("RDBASE");
    std::string filePath = pathName + "/Code/GraphMol/Descriptors/test_data/CH4.mol";

    auto mol = MolFileToMol(filePath, true, false);
    TEST_ASSERT(mol);
    int confId = -1;
    auto field = RDKit::ANI::constructForceField(*mol, "ANI-1ccx", 8);
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
  }
  SECTION("ANI-1x") {
    std::string pathName = getenv("RDBASE");
    std::string filePath = pathName + "/Code/GraphMol/Descriptors/test_data/CH4.mol";

    auto mol = MolFileToMol(filePath, true, false);
    TEST_ASSERT(mol);
    int confId = -1;
    auto field = RDKit::ANI::constructForceField(*mol, "ANI-1x", 8);
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
  }
}