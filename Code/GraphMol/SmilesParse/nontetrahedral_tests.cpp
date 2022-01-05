//
//  Copyright (C) 2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "catch.hpp"

#include <iostream>
#include <fstream>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
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
    CHECK(mol != nullptr);
    inorganicSanitize(*mol);
  }
}

TEST_CASE("TH and @ are equivalent") {
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
  SECTION("no reordering") {
    // clang-format off
    std::vector<std::string> data = {
        "C[Pt@SP1](F)(O)Cl",
        "C[Pt@SP2](F)(O)Cl",
        "C[Pt@TB1](F)(O)(N)Cl",
        "C[Pt@TB2](F)(O)(N)Cl",
        "C[Pt@OH1](F)(O)(N)(Br)Cl",
        "C[Pt@OH2](F)(O)(N)(Br)Cl",
    };
    // clang-format on
    for (const auto &smi : data) {
      SmilesParserParams parseps;
      // we need to skip stereo assignment
      parseps.sanitize = false;
      parseps.removeHs = false;
      std::unique_ptr<RWMol> m{SmilesToMol(smi, parseps)};
      REQUIRE(m);
      m->updatePropertyCache(false);
      SmilesWriteParams writeps;
      writeps.canonical = false;
      // be sure to skip stereo assignment
      m->setProp(common_properties::_StereochemDone, true);
      CHECK(MolToSmiles(*m, writeps) == smi);
    }
  }
}