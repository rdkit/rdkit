//
// Copyright (C) 2018 Greg Landrum
//
#include <catch2/catch_all.hpp>

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
  int confId = -1;
  bool preserveMatrices = true;
  REQUIRE(EHTTools::runMol(*mol, res, confId, preserveMatrices));
  CHECK(res.numElectrons == 30);
  CHECK(res.numOrbitals == 30);
  CHECK(res.numAtoms == 12);
  for (unsigned int i = 0; i < 6; ++i) {
    CHECK(res.atomicCharges[i] == Catch::Approx(-0.026).margin(0.001));
  }
  for (unsigned int i = 6; i < 12; ++i) {
    CHECK(res.atomicCharges[i] == Catch::Approx(0.026).margin(0.001));
  }
  for (unsigned int i = 0; i < 6; ++i) {
    CHECK(res.reducedChargeMatrix[i * res.numOrbitals] ==
          Catch::Approx(0.1615).margin(0.001));
    CHECK(res.reducedChargeMatrix[i * res.numOrbitals + 6] ==
          Catch::Approx(0.1066).margin(0.001));
    CHECK(res.reducedChargeMatrix[i * res.numOrbitals + 9] ==
          Catch::Approx(0.1667).margin(0.001));
  }
  for (unsigned int i = 6; i < 12; ++i) {
    CHECK(res.reducedChargeMatrix[i * res.numOrbitals] ==
          Catch::Approx(0.0052).margin(0.001));
    CHECK(res.reducedChargeMatrix[i * res.numOrbitals + 6] ==
          Catch::Approx(0.0600).margin(0.001));
    CHECK(res.reducedChargeMatrix[i * res.numOrbitals + 9] ==
          Catch::Approx(0.0000).margin(0.001));
  }
  CHECK(res.orbitalEnergies[0] == Catch::Approx(-29.6302).margin(0.001));
  CHECK(res.orbitalEnergies[14] == Catch::Approx(-12.804).margin(0.001));
  CHECK(res.orbitalEnergies[29] == Catch::Approx(67.0404).margin(0.001));

  CHECK(res.totalEnergy == Catch::Approx(-535.026).margin(0.001));
  CHECK(res.fermiEnergy == Catch::Approx(-12.804).margin(0.001));

  CHECK(res.hamiltonianMatrix[0 * res.numOrbitals + 0] ==
        Catch::Approx(-21.4000).margin(0.001));
  CHECK(res.hamiltonianMatrix[0 * res.numOrbitals + 4] ==
        Catch::Approx(-15.3224).margin(0.001));
  CHECK(res.hamiltonianMatrix[4 * res.numOrbitals + 0] ==
        Catch::Approx(0.0000).margin(0.001));
  CHECK(res.overlapMatrix[0 * res.numOrbitals + 0] ==
        Catch::Approx(1.0000).margin(0.001));
  CHECK(res.overlapMatrix[0 * res.numOrbitals + 4] ==
        Catch::Approx(0.4091).margin(0.001));
  CHECK(res.overlapMatrix[4 * res.numOrbitals + 0] ==
        Catch::Approx(0.0000).margin(0.001));
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