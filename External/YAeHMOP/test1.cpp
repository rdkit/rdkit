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
  std::unique_ptr<RWMol> mol(MolFileToMol(pathName + "methanol.mol",sanitize,removeHs));
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms()==6);
  REQUIRE(EHTTools::runMol(*mol));
}
#endif
#if 1
TEST_CASE("benzene", "[basics]") {
  std::string pathName = getenv("RDBASE");
  pathName += "/External/YAeHMOP/test_data/";
  bool sanitize = true;
  bool removeHs = false;
  std::unique_ptr<RWMol> mol(MolFileToMol(pathName + "benzene.mol",sanitize,removeHs));
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms()==12);
  std::cerr<<"RUNNING"<<std::endl;
  REQUIRE(EHTTools::runMol(*mol));
  std::cerr<<"back"<<std::endl;
}
#endif
#if 1
TEST_CASE("phenol", "[basics]") {
  std::string pathName = getenv("RDBASE");
  pathName += "/External/YAeHMOP/test_data/";
  bool sanitize = true;
  bool removeHs = false;
  std::unique_ptr<RWMol> mol(MolFileToMol(pathName + "phenol.mol",sanitize,removeHs));
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms()==13);
  REQUIRE(EHTTools::runMol(*mol));
}
#endif