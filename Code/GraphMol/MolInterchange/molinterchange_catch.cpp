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
  auto andmol = "C[C@H](O)CC[C@@H](C)F |&1:1,5|"_smiles;
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
      CHECK(json.find("stereoGroups") == std::string::npos);
      CHECK(json.find("\"or\"") == std::string::npos);
      CHECK(json.find("[1,4]") == std::string::npos);
    }
    {
      auto json = MolInterchange::MolToJSONData(*andmol);
      CHECK(json.find("stereoGroups") != std::string::npos);
      CHECK(json.find("\"and\"") != std::string::npos);
      CHECK(json.find("[1,5]") != std::string::npos);
      json = MolInterchange::MolToJSONData(*ormol, commonChemPs);
      CHECK(json.find("stereoGroups") == std::string::npos);
      CHECK(json.find("\"and\"") == std::string::npos);
      CHECK(json.find("[1,5]") == std::string::npos);
    }
    {
      auto json = MolInterchange::MolToJSONData(*absmol);
      CHECK(json.find("stereoGroups") != std::string::npos);
      CHECK(json.find("\"abs\"") != std::string::npos);
      CHECK(json.find("[4]") != std::string::npos);
      json = MolInterchange::MolToJSONData(*ormol, commonChemPs);
      CHECK(json.find("stereoGroups") == std::string::npos);
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
      CHECK(ats[1]->getIdx() == 5);
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
        "C[C@H]1OC([C@H](O)F)[C@H](C)[C@H](C)C1[C@@H](O)F |a:1,o1:4,12,&1:7,&2:9|");
  }
  SECTION("writing multiple mols") {
    std::vector<ROMol *> mols{ormol.get(), andmol.get()};
    auto json = MolInterchange::MolsToJSONData(mols);
    CHECK(json.find("stereoGroups") != std::string::npos);
    CHECK(json.find("\"or\"") != std::string::npos);
    CHECK(json.find("[1,4]") != std::string::npos);
    CHECK(json.find("\"and\"") != std::string::npos);
    CHECK(json.find("[1,5]") != std::string::npos);
    json = MolInterchange::MolsToJSONData(mols, commonChemPs);
    CHECK(json.find("stereoGroups") == std::string::npos);
  }
}

TEST_CASE("SubstanceGroups") {
  auto polymol = R"CTAB(
  Mrv2219 12292206542D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 1 0 0
M  V30 BEGIN ATOM
M  V30 1 * -6.7083 3.2083 0 0
M  V30 2 C -5.3747 3.9783 0 0
M  V30 3 O -4.041 3.2083 0 0
M  V30 4 * -2.7073 3.9783 0 0
M  V30 5 C -5.3747 5.5183 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 2 5
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 0 ATOMS=(3 3 2 5) XBONDS=(2 3 1) BRKXYZ=(9 -3.9538 4.3256 0 -
M  V30 -3.0298 2.7252 0 0 0 0) BRKXYZ=(9 -5.4618 2.8611 0 -6.3858 4.4615 0 0 -
M  V30 0 0) CONNECT=HT LABEL=n
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
  REQUIRE(polymol);

  auto supmol = R"CTAB(example
 -ISIS-  10171405052D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 14 15 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 6.4292 -1.1916 0 0 CFG=3
M  V30 2 C 7.0125 -0.6042 0 0
M  V30 3 N 6.4292 -0.0250999 0 0
M  V30 4 C 5.8416 -0.6042 0 0
M  V30 5 C 5.8416 -1.7708 0 0
M  V30 6 N 6.4292 -2.3584 0 0 CFG=3
M  V30 7 C 7.0125 -1.7708 0 0
M  V30 8 O 5.7166 -3.5875 0 0
M  V30 9 C 5.7166 -4.4125 0 0 CFG=3
M  V30 10 C 4.8875 -4.4125 0 0
M  V30 11 C 6.5376 -4.4166 0 0
M  V30 12 C 5.7166 -5.2376 0 0
M  V30 13 C 6.4292 -3.175 0 0
M  V30 14 O 7.1375 -3.5875 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 1
M  V30 5 1 1 5
M  V30 6 1 5 6
M  V30 7 1 6 7
M  V30 8 1 7 1
M  V30 9 1 6 13
M  V30 10 1 8 9
M  V30 11 1 9 10
M  V30 12 1 9 11
M  V30 13 1 9 12
M  V30 14 2 13 14
M  V30 15 1 8 13
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SUP 0 ATOMS=(7 8 9 10 11 12 13 14) XBONDS=(1 9) BRKXYZ=(9 6.24 -2.9 0 -
M  V30 6.24 -2.9 0 0 0 0) CSTATE=(4 9 0 0.82 0) LABEL=Boc SAP=(3 13 6 1)
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
  REQUIRE(supmol);

  MolInterchange::JSONWriteParameters commonChemPs;
  commonChemPs.useRDKitExtensions = false;
  SECTION("writing") {
    {
      auto json = MolInterchange::MolToJSONData(*polymol);
      CHECK(json.find("substanceGroups") != std::string::npos);
      CHECK(json.find("\"TYPE\":\"SRU\"") != std::string::npos);
      CHECK(json.find("\"atoms\":[2,1,4]") != std::string::npos);
      CHECK(json.find("\"bonds\":[2,0]") != std::string::npos);
      CHECK(json.find("\"brackets\"") != std::string::npos);
      CHECK(json.find("\"cstates\"") == std::string::npos);
      CHECK(json.find("\"attachPoints\"") == std::string::npos);
    }
    {
      auto json = MolInterchange::MolToJSONData(*polymol, commonChemPs);
      CHECK(json.find("substanceGroups") == std::string::npos);
    }
    {
      auto json = MolInterchange::MolToJSONData(*supmol);
      CHECK(json.find("substanceGroups") != std::string::npos);
      CHECK(json.find("\"TYPE\":\"SUP\"") != std::string::npos);
      CHECK(json.find("\"atoms\":[7,8,9,10,11,12,13]") != std::string::npos);
      CHECK(json.find("\"bonds\":[8]") != std::string::npos);
      CHECK(json.find("\"brackets\"") != std::string::npos);
      CHECK(json.find("\"cstates\"") != std::string::npos);
      CHECK(json.find("\"bond\":8") != std::string::npos);
      CHECK(json.find("\"vector\"") != std::string::npos);
      CHECK(json.find("\"attachPoints\"") != std::string::npos);
      CHECK(json.find("{\"aIdx\":12,\"lvIdx\":5,\"id\":\"1\"}") !=
            std::string::npos);
    }
  }
  SECTION("parsing") {
    auto json = MolInterchange::MolToJSONData(*supmol);
    auto mols = MolInterchange::JSONDataToMols(json);
    REQUIRE(mols.size() == 1);
    auto sgs = getSubstanceGroups(*mols[0]);
    REQUIRE(sgs.size() == 1);
    CHECK(sgs[0].getAtoms() == std::vector<unsigned>{7, 8, 9, 10, 11, 12, 13});
    CHECK(sgs[0].getBonds() == std::vector<unsigned>{8});
    REQUIRE(sgs[0].getBrackets().size() == 1);
    REQUIRE(sgs[0].getBrackets()[0].size() == 3);
    CHECK(sgs[0].getBrackets()[0][0].x == Approx(6.24).margin(0.01));
    REQUIRE(sgs[0].getCStates().size() == 1);
    CHECK(sgs[0].getCStates()[0].bondIdx == 8);
    CHECK(sgs[0].getCStates()[0].vector.y == Approx(0.82).margin(0.01));
    REQUIRE(sgs[0].getAttachPoints().size() == 1);
    CHECK(sgs[0].getAttachPoints()[0].aIdx == 12);
    CHECK(sgs[0].getAttachPoints()[0].lvIdx == 5);
    CHECK(sgs[0].getAttachPoints()[0].id == "1");
    std::string pval;
    CHECK(sgs[0].getPropIfPresent("LABEL", pval));
    CHECK(pval == "Boc");
  }
}

TEST_CASE("do not crash with null molecules") {
  SECTION("only null") {
    std::vector<ROMol *> mols{nullptr};
    CHECK_THROWS_AS(MolInterchange::MolsToJSONData(mols), ValueErrorException);
  }
  SECTION("not just null") {
    auto tmol = "CCC"_smiles;
    std::vector<ROMol *> mols{tmol.get(), nullptr};
    CHECK_THROWS_AS(MolInterchange::MolsToJSONData(mols), ValueErrorException);
  }
}