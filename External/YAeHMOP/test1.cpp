//
// Copyright (C) 2018 Greg Landrum
//
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do
                           // this in one cpp file per test
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include "EHTTools.h"

#include <GraphMol/FileParsers/FileParsers.h>

using namespace RDKit;

#if 1
TEST_CASE("methanol", "[basics]") {
  std::string pathName = getenv("RDBASE");
  pathName += "/External/YAeHMOP/test_data/";
  bool sanitize = true;
  bool removeHs = false;
  std::unique_ptr<RWMol> mol(
      MolFileToMol(pathName + "methanol.mol", sanitize, removeHs));
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms() == 6);
  EHTTools::EHTResults res;
  REQUIRE(EHTTools::runMol(*mol, res));
}
#endif
#if 1
TEST_CASE("benzene", "[basics]") {
  std::string pathName = getenv("RDBASE");
  pathName += "/External/YAeHMOP/test_data/";
  bool sanitize = true;
  bool removeHs = false;
  std::unique_ptr<RWMol> mol(
      MolFileToMol(pathName + "benzene.mol", sanitize, removeHs));
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms() == 12);
  EHTTools::EHTResults res;
  REQUIRE(EHTTools::runMol(*mol, res));
  for (unsigned int i = 0; i < 6; ++i) {
    CHECK(res.atomicCharges[i] == Approx(-0.026).margin(0.001));
  }
  for (unsigned int i = 6; i < 12; ++i) {
    CHECK(res.atomicCharges[i] == Approx(0.026).margin(0.001));
  }
  CHECK(res.totalEnergy == Approx(-535.026).margin(0.001));
  CHECK(res.fermiEnergy == Approx(-12.804).margin(0.001));
}
#endif
#if 1
TEST_CASE("phenol", "[basics]") {
  std::string pathName = getenv("RDBASE");
  pathName += "/External/YAeHMOP/test_data/";
  bool sanitize = true;
  bool removeHs = false;
  std::unique_ptr<RWMol> mol(
      MolFileToMol(pathName + "phenol.mol", sanitize, removeHs));
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms() == 13);
  EHTTools::EHTResults res;
  REQUIRE(EHTTools::runMol(*mol, res));
}
#endif