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

TEST_CASE("basics", "[basics]") {
  EHTTools::stub();
  REQUIRE(1);
}

TEST_CASE("methanol", "[basics]") {
  std::string pathName = getenv("RDBASE");
  pathName += "/External/YAeHMOP/test_data/";
  std::unique_ptr<RWMol> mol(MolFileToMol(pathName + "methanol.mol",false));
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms()==6);
  REQUIRE(EHTTools::runMol(*mol));
}
