//
//  Copyright (C) 2021 Greg Landrum
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/test.h>
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include "MolInterchange.h"

using namespace RDKit;

TEST_CASE("queries to JSON", "[query]") {
  SECTION("SMARTS 1") {
    auto mol = "C[cH]:[#7+]"_smarts;
    REQUIRE(mol);
    auto json = MolInterchange::MolToJSONData(*mol);
    // no implicit Hs added with queries
    CHECK(json.find("{\"impHs\":3}") == std::string::npos);
    CHECK(json.find("\"descr\":\"AtomHCount\"") != std::string::npos);
    auto nmols = MolInterchange::JSONDataToMols(json);
    CHECK(nmols.size() == 1);
    CHECK(MolToSmarts(*nmols[0]) == MolToSmarts(*mol));
  }
  SECTION("Recursive SMARTS 1") {
    auto mol = "[C;R2]=,#!@[$(CO),r6]"_smarts;
    REQUIRE(mol);
    auto json = MolInterchange::MolToJSONData(*mol);
    CHECK(json.find("\"subquery\":") != std::string::npos);
    auto nmols = MolInterchange::JSONDataToMols(json);
    CHECK(nmols.size() == 1);
    CHECK(MolToSmarts(*nmols[0]) == MolToSmarts(*mol));
  }
  SECTION("mol blocks") {
    auto mol = R"CTAB(
  Mrv2102 04092105442D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Q -1.4583 0.4583 0 0
M  V30 2 C -0.1247 1.2283 0 0
M  V30 3 C 1.209 0.4583 0 0
M  V30 4 C 2.5427 1.2283 0 0
M  V30 5 C 1.209 -1.0817 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 3 4
M  V30 3 5 2 3
M  V30 4 1 3 5 TOPO=2
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    auto json = MolInterchange::MolToJSONData(*mol);
    CHECK(json.find("\"type\":\"Q\"") != std::string::npos);
    auto ps = MolInterchange::JSONParseParameters();
    ps.useHCounts = false;
    auto nmols = MolInterchange::JSONDataToMols(json, ps);
    CHECK(nmols.size() == 1);
    auto mb = MolToV3KMolBlock(*nmols[0]);
    CHECK(mb.find(" Q ") != std::string::npos);
    CHECK(mb.find("TOPO=2") != std::string::npos);
  }
}

TEST_CASE("StereoGroups") {
  auto ormol = "C[C@H](O)C[C@@H](C)F |o1:1,4|"_smiles;
  REQUIRE(ormol);
  auto andmol = "C[C@H](O)C[C@@H](C)F |&1:1,4|"_smiles;
  REQUIRE(andmol);
  auto absmol = "CC(O)C[C@@H](C)F |a:4|"_smiles;
  REQUIRE(absmol);
  MolInterchange::JSONWriteParameters commonChemPs;
  commonChemPs.useRDKitExtensions = false;
  SECTION("writing") {
    {
      auto json = MolInterchange::MolToJSONData(*ormol);
      CHECK(json.find("stereoGroups") != std::string::npos);
      CHECK(json.find("\"or\"") != std::string::npos);
      CHECK(json.find("[1,4]") != std::string::npos);
      json = MolInterchange::MolToJSONData(*ormol, commonChemPs);
      CHECK(json.find("stereogroups") == std::string::npos);
      CHECK(json.find("\"or\"") == std::string::npos);
      CHECK(json.find("[1,4]") == std::string::npos);
    }
    {
      auto json = MolInterchange::MolToJSONData(*andmol);
      CHECK(json.find("stereoGroups") != std::string::npos);
      CHECK(json.find("\"and\"") != std::string::npos);
      CHECK(json.find("[1,4]") != std::string::npos);
      json = MolInterchange::MolToJSONData(*ormol, commonChemPs);
      CHECK(json.find("stereogroups") == std::string::npos);
      CHECK(json.find("\"and\"") == std::string::npos);
      CHECK(json.find("[1,4]") == std::string::npos);
    }
    {
      auto json = MolInterchange::MolToJSONData(*absmol);
      CHECK(json.find("stereoGroups") != std::string::npos);
      CHECK(json.find("\"abs\"") != std::string::npos);
      CHECK(json.find("[4]") != std::string::npos);
      json = MolInterchange::MolToJSONData(*ormol, commonChemPs);
      CHECK(json.find("stereogroups") == std::string::npos);
      CHECK(json.find("\"abs\"") == std::string::npos);
      CHECK(json.find("[4]") == std::string::npos);
    }
  }
  SECTION("reading") {
    {
      auto json = MolInterchange::MolToJSONData(*ormol);
      auto mols = MolInterchange::JSONDataToMols(json);
      REQUIRE(mols.size() == 1);
      auto sgs = mols[0]->getStereoGroups();
      REQUIRE(sgs.size() == 1);
      CHECK(sgs[0].getGroupType() == StereoGroupType::STEREO_OR);
      auto ats = sgs[0].getAtoms();
      REQUIRE(ats.size() == 2);
      CHECK(ats[0]->getIdx() == 1);
      CHECK(ats[1]->getIdx() == 4);
    }
    {
      auto json = MolInterchange::MolToJSONData(*andmol);
      auto mols = MolInterchange::JSONDataToMols(json);
      REQUIRE(mols.size() == 1);
      auto sgs = mols[0]->getStereoGroups();
      REQUIRE(sgs.size() == 1);
      CHECK(sgs[0].getGroupType() == StereoGroupType::STEREO_AND);
      auto ats = sgs[0].getAtoms();
      REQUIRE(ats.size() == 2);
      CHECK(ats[0]->getIdx() == 1);
      CHECK(ats[1]->getIdx() == 4);
    }
    {
      auto json = MolInterchange::MolToJSONData(*absmol);
      auto mols = MolInterchange::JSONDataToMols(json);
      REQUIRE(mols.size() == 1);
      auto sgs = mols[0]->getStereoGroups();
      REQUIRE(sgs.size() == 1);
      CHECK(sgs[0].getGroupType() == StereoGroupType::STEREO_ABSOLUTE);
      auto ats = sgs[0].getAtoms();
      REQUIRE(ats.size() == 1);
      CHECK(ats[0]->getIdx() == 4);
    }
  }
  SECTION("multiple groups") {
    auto mol =
        "C[C@H]1OC([C@H](O)F)[C@@H](C)[C@@H](C)C1[C@@H](O)F |a:1,o1:4,12,&1:7,&2:9|"_smiles;
    REQUIRE(mol);
    auto json = MolInterchange::MolToJSONData(*mol);
    CHECK(json.find("stereoGroups") != std::string::npos);
    CHECK(json.find("\"abs\"") != std::string::npos);
    CHECK(json.find("\"or\"") != std::string::npos);
    CHECK(json.find("\"and\"") != std::string::npos);

    auto mols = MolInterchange::JSONDataToMols(json);
    REQUIRE(mols.size() == 1);
    auto sgs = mols[0]->getStereoGroups();
    REQUIRE(sgs.size() == 4);
    // rather than worry about group order here, just check the CXSMILES:
    auto smi = MolToCXSmiles(*mols[0]);
    CHECK(
        smi ==
        "C[C@@H]1C([C@H](O)F)O[C@H](C)C([C@@H](O)F)[C@@H]1C |a:7,o1:3,10,&1:1,&2:13|");
  }
}