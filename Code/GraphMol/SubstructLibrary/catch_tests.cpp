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

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolBundle.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SubstructLibrary/SubstructLibrary.h>
#include <GraphMol/SubstructLibrary/PatternFactory.h>

using namespace RDKit;

TEST_CASE("querying with a molbundle") {
  std::vector<std::string> libSmiles = {"CCCC", "CCOC", "CCNC"};
  SubstructLibrary ssslib;
  for (const auto smi : libSmiles) {
    std::unique_ptr<RWMol> mol(SmilesToMol(smi));
    REQUIRE(mol);
    ssslib.addMol(*mol);
  }
  bool recursionPossible = true, useChirality = true,
       useQueryQueryMatches = false;
  SECTION("basics") {
    std::vector<std::string> qSmiles = {"CCC", "COC", "CNC"};
    MolBundle bundle;
    for (const auto smi : qSmiles) {
      boost::shared_ptr<ROMol> mol(SmilesToMol(smi));
      REQUIRE(mol);
      bundle.addMol(mol);
    }
    for (auto numThreads : std::vector<int>{1, -1}) {
      CHECK(ssslib.countMatches(bundle, recursionPossible, useChirality,
                                useQueryQueryMatches, numThreads) == 3);
      CHECK(ssslib.hasMatch(bundle, recursionPossible, useChirality,
                            useQueryQueryMatches, numThreads));
      auto libMatches =
          ssslib.getMatches(bundle, recursionPossible, useChirality,
                            useQueryQueryMatches, numThreads);
      CHECK(libMatches.size() == 3);
      std::sort(libMatches.begin(), libMatches.end());
      CHECK(libMatches == std::vector<unsigned int>{0, 1, 2});
    }
  }
  SECTION("some bundle members don't match") {
    std::vector<std::string> qSmiles = {"CCC", "CSC", "CNC"};
    MolBundle bundle;
    for (const auto smi : qSmiles) {
      boost::shared_ptr<ROMol> mol(SmilesToMol(smi));
      REQUIRE(mol);
      bundle.addMol(mol);
    }
    for (auto numThreads : std::vector<int>{1, -1}) {
      CHECK(ssslib.countMatches(bundle, recursionPossible, useChirality,
                                useQueryQueryMatches, numThreads) == 2);
      CHECK(ssslib.hasMatch(bundle, recursionPossible, useChirality,
                            useQueryQueryMatches, numThreads));
      auto libMatches =
          ssslib.getMatches(bundle, recursionPossible, useChirality,
                            useQueryQueryMatches, numThreads);
      CHECK(libMatches.size() == 2);
      std::sort(libMatches.begin(), libMatches.end());
      CHECK(libMatches == std::vector<unsigned int>{0, 2});
    }
  }
  SECTION("no dupes please") {
    std::vector<std::string> qSmiles = {"CCC", "CNC", "CC"};
    MolBundle bundle;
    for (const auto smi : qSmiles) {
      boost::shared_ptr<ROMol> mol(SmilesToMol(smi));
      REQUIRE(mol);
      bundle.addMol(mol);
    }
    for (auto numThreads : std::vector<int>{1, -1}) {
      CHECK(ssslib.countMatches(bundle, recursionPossible, useChirality,
                                useQueryQueryMatches, numThreads) == 3);
      CHECK(ssslib.hasMatch(bundle, recursionPossible, useChirality,
                            useQueryQueryMatches, numThreads));
      auto libMatches =
          ssslib.getMatches(bundle, recursionPossible, useChirality,
                            useQueryQueryMatches, numThreads);
      CHECK(libMatches.size() == 3);
      std::sort(libMatches.begin(), libMatches.end());
      CHECK(libMatches == std::vector<unsigned int>{0, 1, 2});
    }
  }
}