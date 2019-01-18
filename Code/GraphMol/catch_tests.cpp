#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do
                           // this in one cpp file
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/FileParsers/FileParsers.h>
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

TEST_CASE("github #299", "[bug, molops, SSSR]"){
  SECTION("simplified"){
    auto mol = "C13%13%14.C124%18.C25%13%15.C368%17.C4679.C75%10%17.C8%11%14%16.C9%11%12%18.C%10%12%15%16"_smiles;
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms()==9);
  }

  SECTION("old example from molopstest"){
    auto mol = "C123C45C11C44C55C22C33C14C523"_smiles;
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms()==9);
  }

  SECTION("carborane"){
    std::unique_ptr<RWMol> mol(SmilesToMol("[B]1234[B]567[B]118[B]229[B]33%10[B]454[B]656[B]711[B]822[C]933[B]%1045[C]6123",0,false));
    REQUIRE(mol);
    CHECK(mol->getNumAtoms()==12);
    mol->updatePropertyCache(false);
    MolOps::findSSSR(*mol);
    REQUIRE(mol->getRingInfo()->isInitialized());
  }
  SECTION("original report from ChEbI"){
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    std::unique_ptr<RWMol> mol(MolFileToMol(pathName + "ChEBI_50252.mol",false));
    REQUIRE(mol);
    CHECK(mol->getNumAtoms()==80);
    mol->updatePropertyCache(false);
    MolOps::findSSSR(*mol);
    REQUIRE(mol->getRingInfo()->isInitialized());

  }
}

TEST_CASE("github #2224", "[bug, molops, removeHs, query]"){
  SECTION("the original report"){
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    std::unique_ptr<RWMol> mol(MolFileToMol(pathName + "github2224_1.mol"));
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms()==7);
  }
  SECTION("basics") {
  SmilesParserParams ps;
    ps.removeHs = false;
    ps.sanitize = true;
    std::unique_ptr<ROMol> mol(SmilesToMol("C[H]", ps));
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms()==2);
    { // The H without a query is removed
      std::unique_ptr<ROMol> m2(MolOps::removeHs(*mol));
      CHECK(m2->getNumAtoms()==1);
    }
    { // but if we add a query feature it's not removed
      RWMol m2(*mol);
      m2.replaceAtom(1,new QueryAtom(1));
      m2.getAtomWithIdx(1)->setAtomicNum(1);
      MolOps::removeHs(m2);
      CHECK(m2.getNumAtoms()==2);
    }
  }
}

