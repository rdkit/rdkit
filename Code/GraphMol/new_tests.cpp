#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do
                           // this in one cpp file
#include <catch/catch.hpp>

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
  std::unique_ptr<RWMol> mol(SmilesToMol("C1=CC=CC=C1-c2ccccc2", false, false));
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms() == 12);

  SECTION("properties") {
    mol->updatePropertyCache();
    REQUIRE(mol->getAtomWithIdx(0)->getTotalNumHs() == 1);
    REQUIRE(!mol->getAtomWithIdx(0)->getIsAromatic());
    REQUIRE(mol->getAtomWithIdx(6)->getIsAromatic());
    SECTION("kekulize") {
      unsigned int opThatFailed;
      MolOps::sanitizeMol(*mol, opThatFailed, MolOps::SANITIZE_KEKULIZE);
      REQUIRE(!mol->getAtomWithIdx(0)->getIsAromatic());
      REQUIRE(!mol->getAtomWithIdx(6)->getIsAromatic());
    }
    SECTION("aromaticity") {
      unsigned int opThatFailed;
      MolOps::sanitizeMol(*mol, opThatFailed, MolOps::SANITIZE_SETAROMATICITY);
      mol->debugMol(std::cerr);
      REQUIRE(mol->getAtomWithIdx(0)->getIsAromatic());
      REQUIRE(mol->getAtomWithIdx(6)->getIsAromatic());
    }
  }
}
