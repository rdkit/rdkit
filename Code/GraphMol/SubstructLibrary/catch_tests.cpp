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
  for (const auto &smi : libSmiles) {
    std::unique_ptr<RWMol> mol(SmilesToMol(smi));
    REQUIRE(mol);
    ssslib.addMol(*mol);
  }
  SubstructMatchParameters params;
  SECTION("basics") {
    std::vector<std::string> qSmiles = {"CCC", "COC", "CNC"};
    MolBundle bundle;
    for (const auto &smi : qSmiles) {
      boost::shared_ptr<ROMol> mol(SmilesToMol(smi));
      REQUIRE(mol);
      bundle.addMol(mol);
    }
    for (auto numThreads : std::vector<int>{1, -1}) {
      CHECK(ssslib.countMatches(bundle, params, numThreads) == 3);
      CHECK(ssslib.hasMatch(bundle, params, numThreads));
      auto libMatches = ssslib.getMatches(bundle, params, numThreads);
      CHECK(libMatches.size() == 3);
      std::sort(libMatches.begin(), libMatches.end());
      CHECK(libMatches == std::vector<unsigned int>{0, 1, 2});
    }
  }
  SECTION("some bundle members don't match") {
    std::vector<std::string> qSmiles = {"CCC", "CSC", "CNC"};
    MolBundle bundle;
    for (const auto &smi : qSmiles) {
      boost::shared_ptr<ROMol> mol(SmilesToMol(smi));
      REQUIRE(mol);
      bundle.addMol(mol);
    }
    for (auto numThreads : std::vector<int>{1, -1}) {
      CHECK(ssslib.countMatches(bundle, params, numThreads) == 2);
      CHECK(ssslib.hasMatch(bundle, params, numThreads));
      auto libMatches = ssslib.getMatches(bundle, params, numThreads);
      CHECK(libMatches.size() == 2);
      std::sort(libMatches.begin(), libMatches.end());
      CHECK(libMatches == std::vector<unsigned int>{0, 2});
    }
  }
  SECTION("no dupes please") {
    std::vector<std::string> qSmiles = {"CCC", "CNC", "CC"};
    MolBundle bundle;
    for (const auto &smi : qSmiles) {
      boost::shared_ptr<ROMol> mol(SmilesToMol(smi));
      REQUIRE(mol);
      bundle.addMol(mol);
    }
    for (auto numThreads : std::vector<int>{1, -1}) {
      CHECK(ssslib.countMatches(bundle, params, numThreads) == 3);
      CHECK(ssslib.hasMatch(bundle, params, numThreads));
      auto libMatches = ssslib.getMatches(bundle, params, numThreads);
      CHECK(libMatches.size() == 3);
      std::sort(libMatches.begin(), libMatches.end());
      CHECK(libMatches == std::vector<unsigned int>{0, 1, 2});
    }
  }
}

TEST_CASE("using modified query parameters") {
  std::vector<std::string> libSmiles = {"C[C@H](F)Cl", "C[C@@H](F)Cl",
                                        "C[CH](F)Cl"};
  SubstructLibrary ssslib;
  for (const auto &smi : libSmiles) {
    std::unique_ptr<RWMol> mol(SmilesToMol(smi));
    REQUIRE(mol);
    ssslib.addMol(*mol);
  }
  SECTION("basics") {
    auto qm = "C[C@@H](F)Cl"_smiles;
    {  // by default stereo is not used
      SubstructMatchParameters params;
      CHECK(ssslib.countMatches(*qm, params) == 3);
      CHECK(ssslib.hasMatch(*qm, params));
      auto libMatches = ssslib.getMatches(*qm, params);
      CHECK(libMatches.size() == 3);
      std::sort(libMatches.begin(), libMatches.end());
      CHECK(libMatches == std::vector<unsigned int>{0, 1, 2});
    }
    {  // use stereo
      SubstructMatchParameters params;
      params.useChirality = true;
      CHECK(ssslib.countMatches(*qm, params) == 1);
      CHECK(ssslib.hasMatch(*qm, params));
      auto libMatches = ssslib.getMatches(*qm, params);
      CHECK(libMatches.size() == 1);
      std::sort(libMatches.begin(), libMatches.end());
      CHECK(libMatches == std::vector<unsigned int>{1});
    }
  }
}

TEST_CASE("searchOrder") {
  std::vector<std::string> libSmiles = {"CCCOC", "CCCCOCC", "CCOC", "COC",
                                        "CCCCCOC"};
  boost::shared_ptr<MolHolder> mholder(new MolHolder());
  boost::shared_ptr<PatternHolder> fpholder(new PatternHolder());

  SubstructLibrary ssslib(mholder, fpholder);

  for (const auto &smi : libSmiles) {
    std::unique_ptr<RWMol> mol(SmilesToMol(smi));
    REQUIRE(mol);
    ssslib.addMol(*mol);
  }
  SECTION("basics") {
    auto qm = "COC"_smiles;
    for (auto numThreads : std::vector<int>{1, -1}) {
      {  // default search order
        ssslib.setSearchOrder(std::vector<unsigned int>());
        unsigned int maxResults = 2;
        SubstructMatchParameters params;
        CHECK(ssslib.countMatches(*qm, params, numThreads) == 5);
        CHECK(ssslib.hasMatch(*qm, params));
        auto libMatches =
            ssslib.getMatches(*qm, params, numThreads, maxResults);
        CHECK(libMatches.size() == 2);
        CHECK(libMatches == std::vector<unsigned int>{0, 1});
      }
      {  // specified search order
        unsigned int maxResults = 2;
        std::vector<unsigned int> searchOrder = {3, 2, 0, 1, 4};
        ssslib.setSearchOrder(searchOrder);
        SubstructMatchParameters params;
        CHECK(ssslib.countMatches(*qm, params, numThreads) == 5);
        CHECK(ssslib.hasMatch(*qm, params));
        auto libMatches =
            ssslib.getMatches(*qm, params, numThreads, maxResults);
        CHECK(libMatches.size() == 2);
        CHECK(libMatches == std::vector<unsigned int>{3, 2});
      }
      {  // search order not the full length
        unsigned int maxResults = 2;
        std::vector<unsigned int> searchOrder = {3, 2, 0};
        ssslib.setSearchOrder(searchOrder);
        SubstructMatchParameters params;
        CHECK(ssslib.countMatches(*qm, params, numThreads) == 3);
        CHECK(ssslib.hasMatch(*qm, params));
        auto libMatches =
            ssslib.getMatches(*qm, params, numThreads, maxResults);
        CHECK(libMatches.size() == 2);
        CHECK(libMatches == std::vector<unsigned int>{3, 2});
      }
    }
  }

  SECTION("overly long searchOrder") {
    auto qm = "COC"_smiles;
    unsigned int maxResults = 2;
    int numThreads = -1;
    std::vector<unsigned int> searchOrder = {3, 2, 0, 1, 4, 2, 1, 4, 0};
    ssslib.setSearchOrder(searchOrder);
    SubstructMatchParameters params;
    CHECK(ssslib.countMatches(*qm, params, numThreads) == 5);
    CHECK(ssslib.hasMatch(*qm, params));
    auto libMatches = ssslib.getMatches(*qm, params, numThreads, maxResults);
    CHECK(libMatches.size() == 2);
    CHECK(libMatches == std::vector<unsigned int>{3, 2});
  }
  SECTION("bogus indices in searchOrder") {
    auto qm = "COC"_smiles;
    std::vector<unsigned int> searchOrder = {3, 12, 0};
    CHECK_THROWS_AS(ssslib.setSearchOrder(searchOrder), IndexErrorException);
  }
#ifdef RDK_USE_BOOST_SERIALIZATION
  SECTION("serialization") {
    std::vector<unsigned int> searchOrder = {3, 0, 1};
    ssslib.setSearchOrder(searchOrder);
    std::string pickle = ssslib.Serialize();
    SubstructLibrary serialized;
    serialized.initFromString(pickle);
    CHECK(serialized.size() == ssslib.size());
    CHECK(serialized.getSearchOrder() == ssslib.getSearchOrder());
  }
#endif
}

void setSearchSmallestFirst(SubstructLibrary &ssslib) {
  const auto holder = ssslib.getMolHolder();
  std::vector<unsigned int> searchOrder(holder->size());
  std::iota(searchOrder.begin(), searchOrder.end(), 0);
  std::stable_sort(searchOrder.begin(), searchOrder.end(),
                   [holder](unsigned int i1, unsigned int i2) {
                     return holder->getMol(i1)->getNumAtoms() <
                            holder->getMol(i2)->getNumAtoms();
                   });
  ssslib.setSearchOrder(searchOrder);
}

TEST_CASE("searchOrderFunctionDemo") {
  std::vector<std::string> libSmiles = {"CCCOC", "CCCCOCC", "CCOC", "COC",
                                        "CCCCCOC"};
  boost::shared_ptr<MolHolder> mholder(new MolHolder());
  boost::shared_ptr<PatternHolder> fpholder(new PatternHolder());

  SubstructLibrary ssslib(mholder, fpholder);

  for (const auto &smi : libSmiles) {
    std::unique_ptr<RWMol> mol(SmilesToMol(smi));
    REQUIRE(mol);
    ssslib.addMol(*mol);
  }
  SECTION("basics") {
    auto qm = "COC"_smiles;
    setSearchSmallestFirst(ssslib);
    auto libMatches = ssslib.getMatches(*qm);
    CHECK(libMatches.size() == 5);
    CHECK(libMatches == std::vector<unsigned int>{3, 2, 0, 1, 4});
  }
}