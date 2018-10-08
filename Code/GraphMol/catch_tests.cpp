#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do
                           // this in one cpp file
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>

using namespace RDKit;

TEST_CASE("SMILES Parsing works", "[molops]") {
  std::unique_ptr<RWMol> mol(SmilesToMol("C1CC1"));
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms() == 3);
}

TEST_CASE("Sanitization tests", "[molops]") {
  std::unique_ptr<RWMol> mol(SmilesToMol("C1=CC=CC=C1Cc2ccccc2", false, false));
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms() == 13);

  SECTION("properties") {
    mol->updatePropertyCache();
    CHECK(mol->getAtomWithIdx(0)->getTotalNumHs() == 1);
    CHECK(!mol->getAtomWithIdx(0)->getIsAromatic());
    CHECK(mol->getAtomWithIdx(7)->getIsAromatic());
    SECTION("aromaticity") {
      unsigned int opThatFailed;
      MolOps::sanitizeMol(*mol, opThatFailed, MolOps::SANITIZE_SETAROMATICITY);
      // mol->debugMol(std::cerr);
      CHECK(mol->getAtomWithIdx(7)->getIsAromatic());
      // blocked by #1730
      // CHECK(mol->getAtomWithIdx(0)->getIsAromatic());
    }
    SECTION("kekulize") {
      unsigned int opThatFailed;
      MolOps::sanitizeMol(*mol, opThatFailed, MolOps::SANITIZE_KEKULIZE);
      CHECK(!mol->getAtomWithIdx(0)->getIsAromatic());
      CHECK(!mol->getAtomWithIdx(7)->getIsAromatic());
    }
  }
}

TEST_CASE("Github #2062", "[bug, molops]") {
  SmilesParserParams ps;
  ps.removeHs = false;
  ps.sanitize = true;
  std::unique_ptr<RWMol> mol(SmilesToMol("[C:1][C:2]([H:3])([H])[O:4][H]", ps));
  REQUIRE(mol);
  CHECK(mol->getNumAtoms() == 6);
  mol->getAtomWithIdx(1)->setProp("intProp", 42);
  MolOps::mergeQueryHs(*mol);
  CHECK(mol->getNumAtoms() == 3);
  SECTION("basics") { CHECK(mol->getAtomWithIdx(1)->getAtomMapNum() == 2); }
  SECTION("other props") {
    REQUIRE(mol->getAtomWithIdx(1)->hasProp("intProp"));
    CHECK(mol->getAtomWithIdx(1)->getProp<int>("intProp") == 42);
  }
}

TEST_CASE("Github #2086", "[bug, molops]") {
  SECTION("reported version") {
    auto mol = "C1CCCC1"_smiles;
    REQUIRE(mol);
    MolOps::addHs(*mol);
    REQUIRE(mol->getNumAtoms() == 15);
    mol->removeBond(4, 13);
    MolOps::removeHs(*mol);
    REQUIRE(mol->getNumAtoms() == 6);
  }
}
