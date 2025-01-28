//
//  Copyright (C) 2023 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include <GraphMol/RWMol.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;

TEST_CASE("convert haptic bond to explicit dative bonds") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/MolStandardize/test_data/";
  std::vector<std::tuple<std::string, int, int, int>> test_data{
      {"ferrocene.mol", 10, 26, 10},
      {"MOL_00002.mol", 14, 26, 10},
      {"MOL_00010.mol", 5, 46, 5}};
  {
    // doing in place
    for (const auto &td : test_data) {
      std::unique_ptr<RWMol> mol(MolFileToMol(pathName + std::get<0>(td)));
      REQUIRE(mol);
      MolOps::hapticBondsToDative(*mol);
      int numDats = 0, numMetDats = 0;
      for (const auto &b : mol->bonds()) {
        if (b->getBondType() == Bond::DATIVE) {
          ++numDats;
          if (b->getEndAtom()->getAtomicNum() == get<2>(td)) {
            ++numMetDats;
          }
        }
      }
      CHECK(numDats == std::get<1>(td));
      CHECK(numMetDats == std::get<3>(td));
    }
  }
  {
    // doing via temporary
    for (const auto &td : test_data) {
      std::unique_ptr<ROMol> mol(MolFileToMol(pathName + std::get<0>(td)));
      REQUIRE(mol);
      std::unique_ptr<ROMol> newMol(MolOps::hapticBondsToDative(*mol));
      int numDats = 0, numMetDats = 0;
      for (const auto &b : newMol->bonds()) {
        if (b->getBondType() == Bond::DATIVE) {
          ++numDats;
          if (b->getEndAtom()->getAtomicNum() == get<2>(td)) {
            ++numMetDats;
          }
        }
      }
      CHECK(numDats == std::get<1>(td));
      CHECK(numMetDats == std::get<3>(td));
    }
  }
}

TEST_CASE("convert explicit dative bonds to haptic bond") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/MolStandardize/test_data/";
  std::vector<std::string> test_data{"ferrocene.mol", "MOL_00002.mol",
                                     "MOL_00010.mol"};
  // MOL_00010.mol only has 1 dummy, so duplicating the test to fit in.
  std::vector<std::tuple<int, int, float, float, float, float>> exp_res{
      {11, 12, -36.5, 13.5, -36.4, 8.6},
      {43, 44, -36.5, 13.5, -36.4, 8.6},
      {82, 82, -1.0, -0.7, -1.0, -0.7}};
  {
    // doing in place
    for (size_t i = 0; i < 3; ++i) {
      std::unique_ptr<RWMol> mol(MolFileToMol(pathName + test_data[i]));
      REQUIRE(mol);
      auto initSmi = MolToSmiles(*mol);
      MolOps::hapticBondsToDative(*mol);
      MolOps::dativeBondsToHaptic(*mol);
      CHECK(initSmi == MolToSmiles(*mol));
      // Check the dummy atoms are in the right place.
      auto [d1, d2, d1x, d1y, d2x, d2y] = exp_res[i];
      auto dummy1pos = mol->getConformer().getAtomPos(d1);
      REQUIRE_THAT(dummy1pos.x, Catch::Matchers::WithinAbs(d1x, 0.1));
      REQUIRE_THAT(dummy1pos.y, Catch::Matchers::WithinAbs(d1y, 0.1));
      auto dummy2pos = mol->getConformer().getAtomPos(d2);
      REQUIRE_THAT(dummy2pos.x, Catch::Matchers::WithinAbs(d2x, 0.1));
      REQUIRE_THAT(dummy2pos.y, Catch::Matchers::WithinAbs(d2y, 0.1));
    }
  }
  {
    // doing via temporary
    for (const auto &td : test_data) {
      std::unique_ptr<ROMol> mol(MolFileToMol(pathName + td));
      REQUIRE(mol);
      auto initSmi = MolToSmiles(*mol);
      std::unique_ptr<ROMol> mol1{MolOps::hapticBondsToDative(*mol)};
      std::unique_ptr<ROMol> mol2{MolOps::dativeBondsToHaptic(*mol1)};
      CHECK(initSmi == MolToSmiles(*mol2));
    }
  }
}

TEST_CASE("get haptic bond end points") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/MolStandardize/test_data/";
  bool sanitize = false;
  std::unique_ptr<RWMol> mol(
      MolFileToMol(pathName + "MOL_00002.mol", sanitize));
  REQUIRE(mol);
  auto bond = mol->getBondWithIdx(0);
  auto endPts = MolOps::details::hapticBondEndpoints(bond);
  CHECK(std::vector<int>{1, 2, 3, 4, 0} == endPts);
  bond = mol->getBondWithIdx(11);
  endPts = MolOps::details::hapticBondEndpoints(bond);
  CHECK(std::vector<int>{6, 7, 8, 9, 5} == endPts);
  // dative but not haptic
  bond = mol->getBondWithIdx(49);
  endPts = MolOps::details::hapticBondEndpoints(bond);
  CHECK(std::vector<int>{} == endPts);
  // not dative
  bond = mol->getBondWithIdx(1);
  endPts = MolOps::details::hapticBondEndpoints(bond);
  CHECK(std::vector<int>{} == endPts);
}

TEST_CASE("Rings and dative bonds") {
  SECTION("ferrocene") {
    auto m =
        "C12->[Fe+2]3456789(<-C1=C->3[CH-]->4C->5=2)<-C1=C->6[CH-]->7C->8=C->91"_smiles;
    REQUIRE(m);
    std::vector<std::vector<int>> rings;
    auto nrings = MolOps::symmetrizeSSSR(*m, rings);
    CHECK(nrings == 2);
    CHECK(rings[0].size() == 5);
    CHECK(rings[1].size() == 5);
  }
}

TEST_CASE("aromaticity and dative bonds") {
  SECTION("ferrocene") {
    auto m =
        "C12->[Fe+2]3456789(<-C1=C->3[CH-]->4C->5=2)<-C1=C->6[CH-]->7C->8=C->91"_smiles;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(0)->getIsAromatic());
    CHECK(!m->getAtomWithIdx(1)->getIsAromatic());
  }
}

TEST_CASE("Github 6252 - wrong endpoints after dativeBondsToHaptic") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/MolStandardize/test_data/";
  std::string ferroccene = pathName + "ferrocene.mol";
  {
    std::unique_ptr<RWMol> mol(MolFileToMol(ferroccene));
    REQUIRE(mol);
    auto initSmi = MolToSmiles(*mol);
    MolOps::hapticBondsToDative(*mol);
    MolOps::dativeBondsToHaptic(*mol);
    CHECK(initSmi == MolToSmiles(*mol));
    std::string endpts;
    std::string attach;
    auto bond10 = mol->getBondWithIdx(10);
    CHECK(bond10->getBondType() == Bond::DATIVE);
    CHECK(bond10->getPropIfPresent(common_properties::_MolFileBondEndPts,
                                   endpts));
    CHECK(endpts == "(5 2 3 1 4 5)");
    auto bond11 = mol->getBondWithIdx(11);
    CHECK(bond11->getBondType() == Bond::DATIVE);
    CHECK(bond11->getPropIfPresent(common_properties::_MolFileBondEndPts,
                                   endpts));
    CHECK(endpts == "(5 7 8 6 9 10)");
  }
}