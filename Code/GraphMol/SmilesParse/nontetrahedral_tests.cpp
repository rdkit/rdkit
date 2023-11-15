//
//  Copyright (C) 2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include <iostream>
#include <fstream>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <boost/algorithm/string.hpp>

using namespace RDKit;

void inorganicSanitize(RDKit::RWMol &mol) {
  unsigned int failed = 0;
  unsigned int flags = RDKit::MolOps::SANITIZE_ALL;
  flags &= ~RDKit::MolOps::SANITIZE_CLEANUP;
  flags &= ~RDKit::MolOps::SANITIZE_PROPERTIES;

  mol.updatePropertyCache(false);
  RDKit::MolOps::sanitizeMol(mol, failed, flags);
}

TEST_CASE("bulk parse test") {
  Chirality::setAllowNontetrahedralChirality(true);
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/SmilesParse/test_data/inorganic_stereo.smi";
  std::ifstream ins(pathName);
  REQUIRE(ins.good());
  while (!ins.eof()) {
    auto inl = getLine(ins);
    if (inl.empty()) {
      continue;
    }
    std::vector<std::string> tokens;
    boost::split(tokens, inl, boost::is_any_of(" \t"));
    if (tokens.size() < 2) {
      continue;
    }
    SmilesParserParams ps;
    ps.sanitize = false;
    std::unique_ptr<RWMol> mol{SmilesToMol(tokens[0], ps)};
    if (!mol) {
      std::cerr << tokens[1] << std::endl;
    }
    REQUIRE(mol);
    inorganicSanitize(*mol);
    bool hasNontet = false;
    for (const auto atom : mol->atoms()) {
      if (Chirality::hasNonTetrahedralStereo(atom)) {
        hasNontet = true;
        break;
      }
    }
    if (!hasNontet) {
      std::cerr << tokens[0] << std::endl;
    }
    CHECK(hasNontet);
  }
}

TEST_CASE("TH and @ are equivalent") {
  Chirality::setAllowNontetrahedralChirality(true);
  SECTION("@TH") {
    auto m = "F[C@THH](C)O"_smiles;
    REQUIRE(m);
    CHECK(MolToSmiles(*m) == "C[C@H](O)F");
  }
  SECTION("@TH1") {
    auto m = "F[C@TH1H](C)O"_smiles;
    REQUIRE(m);
    CHECK(MolToSmiles(*m) == "C[C@H](O)F");
  }
  SECTION("@TH2") {
    auto m = "F[C@TH2H](C)O"_smiles;
    REQUIRE(m);
    CHECK(MolToSmiles(*m) == "C[C@@H](O)F");
  }
}

TEST_CASE("non-canonical non-tetrahedral output") {
  Chirality::setAllowNontetrahedralChirality(true);
  SECTION("no reordering") {
    // clang-format off
    std::vector<std::string> data = {
        "C[Pt@SP1](F)(O)Cl",        "C[Pt@SP2](F)(O)Cl",
        "C[Pt@TB1](F)(O)(N)Cl",     "C[Pt@TB2](F)(O)(N)Cl",
        "C[Pt@OH1](F)(O)(N)(Br)Cl", "C[Pt@OH2](F)(O)(N)(Br)Cl",
    };
    // clang-format on
    for (const auto &smi : data) {
      std::unique_ptr<RWMol> m{SmilesToMol(smi)};
      REQUIRE(m);
      SmilesWriteParams writeps;
      writeps.canonical = false;
      // be sure to skip stereo assignment
      m->setProp(common_properties::_StereochemDone, true);
      CHECK(MolToSmiles(*m, writeps) == smi);
    }
  }
  SECTION("reordering") {
    // clang-format off
    std::vector<std::pair<std::string,std::string>> data = {
        {"F[Pt@SP1](C)(O)Cl","C[Pt@SP3](O)(F)Cl",},
        {"F[Pt@SP2](C)(O)Cl","C[Pt@SP1](O)(F)Cl"},
        {"S[As@TB1](F)(Cl)(Br)N","N[As@TB6](F)(S)(Cl)Br"},
        {"C[Pt@OH1](F)(O)(N)(Br)Cl","C[Pt@OH16](N)(O)(F)(Cl)Br"},
    };
    // clang-format on
    for (const auto &pr : data) {
      std::unique_ptr<RWMol> m{SmilesToMol(pr.first)};
      REQUIRE(m);
      CHECK(MolToSmiles(*m) == pr.second);
    }
  }
}

TEST_CASE("SP getChiralAcrossBond et al.") {
  Chirality::setAllowNontetrahedralChirality(true);
  SECTION("basics") {
    {
      auto m = "C[Pt@SP1](F)(O)Cl"_smiles;
      REQUIRE(m);
      std::vector<std::pair<unsigned int, unsigned int>> bpairs = {{0, 2},
                                                                   {1, 3}};
      for (auto pr : bpairs) {
        CHECK(Chirality::getChiralAcrossBond(m->getAtomWithIdx(1),
                                             m->getBondWithIdx(pr.first))
                  ->getIdx() == pr.second);
        CHECK(Chirality::getChiralAcrossBond(m->getAtomWithIdx(1),
                                             m->getBondWithIdx(pr.second))
                  ->getIdx() == pr.first);
      }

      // we just need to check the other forms for one of the sets of pairs
      // since we know it's the same code underneath
      std::vector<std::pair<unsigned int, unsigned int>> apairs = {{0, 3},
                                                                   {2, 4}};
      for (auto pr : apairs) {
        CHECK(Chirality::getChiralAcrossAtom(m->getAtomWithIdx(1),
                                             m->getAtomWithIdx(pr.first))
                  ->getIdx() == pr.second);
        CHECK(Chirality::getChiralAcrossAtom(m->getAtomWithIdx(1),
                                             m->getAtomWithIdx(pr.second))
                  ->getIdx() == pr.first);
      }
      std::vector<std::pair<unsigned int, unsigned int>> abpairs = {
          {0, 2}, {2, 3}, {3, 0}, {4, 1}};
      for (auto pr : abpairs) {
        CHECK(Chirality::getChiralAcrossBond(m->getAtomWithIdx(1),
                                             m->getAtomWithIdx(pr.first))
                  ->getIdx() == pr.second);
        CHECK(Chirality::getChiralAcrossAtom(m->getAtomWithIdx(1),
                                             m->getBondWithIdx(pr.second))
                  ->getIdx() == pr.first);
      }
    }
    {
      auto m = "C[Pt@SP2](F)(O)Cl"_smiles;
      REQUIRE(m);
      std::vector<std::pair<unsigned int, unsigned int>> pairs = {{0, 1},
                                                                  {2, 3}};
      for (auto pr : pairs) {
        CHECK(Chirality::getChiralAcrossBond(m->getAtomWithIdx(1),
                                             m->getBondWithIdx(pr.first))
                  ->getIdx() == pr.second);
        CHECK(Chirality::getChiralAcrossBond(m->getAtomWithIdx(1),
                                             m->getBondWithIdx(pr.second))
                  ->getIdx() == pr.first);
      }
    }
    {
      auto m = "C[Pt@SP3](F)(O)Cl"_smiles;
      REQUIRE(m);
      std::vector<std::pair<unsigned int, unsigned int>> pairs = {{0, 3},
                                                                  {1, 2}};
      for (auto pr : pairs) {
        CHECK(Chirality::getChiralAcrossBond(m->getAtomWithIdx(1),
                                             m->getBondWithIdx(pr.first))
                  ->getIdx() == pr.second);
        CHECK(Chirality::getChiralAcrossBond(m->getAtomWithIdx(1),
                                             m->getBondWithIdx(pr.second))
                  ->getIdx() == pr.first);
      }
    }
  }
  SECTION("3 real ligands") {
    {
      auto m = "C[Pt@SP1](F)O"_smiles;
      REQUIRE(m);
      CHECK(Chirality::getChiralAcrossBond(m->getAtomWithIdx(1),
                                           m->getBondWithIdx(0))
                ->getIdx() == 2);
      CHECK(Chirality::getChiralAcrossBond(m->getAtomWithIdx(1),
                                           m->getBondWithIdx(2))
                ->getIdx() == 0);
      CHECK(Chirality::getChiralAcrossBond(m->getAtomWithIdx(1),
                                           m->getBondWithIdx(1)) == nullptr);
    }
  }
}
TEST_CASE("getChiralAcross edges") {
  Chirality::setAllowNontetrahedralChirality(true);
  SECTION("central atom isn't chiral") {
    auto m = "C[Pt](F)(O)CC"_smiles;
    REQUIRE(m);
    CHECK(Chirality::getChiralAcrossBond(m->getAtomWithIdx(1),
                                         m->getBondWithIdx(1)) == nullptr);
  }
  SECTION("others") {
    auto m = "C[Pt@SP1](F)(O)CC"_smiles;
    REQUIRE(m);
    // not the central atom:
    CHECK(Chirality::getChiralAcrossBond(m->getAtomWithIdx(0),
                                         m->getBondWithIdx(1)) == nullptr);
    // not a bond to central atom
    CHECK(Chirality::getChiralAcrossBond(m->getAtomWithIdx(1),
                                         m->getBondWithIdx(4)) == nullptr);
    CHECK(Chirality::getChiralAcrossAtom(m->getAtomWithIdx(1),
                                         m->getBondWithIdx(4)) == nullptr);
    // atom not connected to central atom
    CHECK(Chirality::getChiralAcrossBond(m->getAtomWithIdx(1),
                                         m->getAtomWithIdx(5)) == nullptr);
    CHECK(Chirality::getChiralAcrossAtom(m->getAtomWithIdx(1),
                                         m->getAtomWithIdx(5)) == nullptr);
  }
}

TEST_CASE("hasNonTetrahedralStereo") {
  Chirality::setAllowNontetrahedralChirality(true);
  std::vector<std::pair<std::string, bool>> data = {
      {"C[Pt@SP1](F)(O)CC", true},     {"C[Pt@SP](F)(O)CC", true},
      {"C[Pt](F)(O)CC", false},        {"C[Pt@](F)(O)CC", false},
      {"C[Pt@TH1](F)(O)CC", false},    {"C[Pt@TB1](N)(F)(O)CC", true},
      {"C[Pt@TB](N)(F)(O)CC", true},   {"C[Pt@OH1](N)(N)(F)(O)CC", true},
      {"C[Pt@OH](N)(N)(F)(O)CC", true}};
  for (const auto &pr : data) {
    std::unique_ptr<RWMol> m{SmilesToMol(pr.first)};
    REQUIRE(m);
    CHECK(Chirality::hasNonTetrahedralStereo(m->getAtomWithIdx(1)) ==
          pr.second);
  }
}

TEST_CASE("zero permutation is in SMILES") {
  Chirality::setAllowNontetrahedralChirality(true);
  std::vector<std::string> smis = {"CC[Pt@SP](C)(O)F", "C[Pt@TB](N)(F)(Cl)Br",
                                   "C[Pt@OH](N)(F)(Cl)(Br)I"};
  for (auto smi : smis) {
    std::unique_ptr<RWMol> m{SmilesToMol(smi)};
    REQUIRE(m);
    CHECK(MolToSmiles(*m) == smi);
  }
}

TEST_CASE("do not write bogus permutation values") {
  Chirality::setAllowNontetrahedralChirality(true);
  auto m = "CC[Pt@SP](C)(O)F"_smiles;
  REQUIRE(m);
  m->getAtomWithIdx(2)->setProp(common_properties::_chiralPermutation, 10);
  CHECK_THROWS_AS(MolToSmiles(*m), ValueErrorException);
  m->getAtomWithIdx(2)->setProp(common_properties::_chiralPermutation, -1);
  CHECK_THROWS_AS(MolToSmiles(*m), ValueErrorException);
}

TEST_CASE("do not read from/write to SMILES when disabled") {
  Chirality::setAllowNontetrahedralChirality(false);
  SECTION("parsing") {
    auto m = "CC[Pt@SP](C)(O)F"_smiles;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(2)->getChiralTag() ==
          Atom::ChiralType::CHI_UNSPECIFIED);
    CHECK(
        !m->getAtomWithIdx(2)->hasProp(common_properties::_chiralPermutation));
  }
  SECTION("writing") {
    auto m = "CC[Pt](C)(O)F"_smiles;
    REQUIRE(m);
    m->getAtomWithIdx(2)->setChiralTag(Atom::ChiralType::CHI_SQUAREPLANAR);
    m->getAtomWithIdx(2)->setProp(common_properties::_chiralPermutation, 1);
    auto smi = MolToSmiles(*m);
    CHECK(smi == "CC[Pt](C)(O)F");
  }
}