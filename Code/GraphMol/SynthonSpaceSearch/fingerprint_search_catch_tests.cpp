//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#include <algorithm>
#include <fstream>

#include <GraphMol/SubstructLibrary/SubstructLibrary.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <catch2/catch_all.hpp>

using namespace RDKit;
using namespace RDKit::SynthonSpaceSearch;

using namespace RDKit::SynthonSpaceSearch::details;

const char *rdbase = getenv("RDBASE");

std::pair<std::set<std::string>, std::set<std::string>> bruteForceSearch(
    const std::string &enumLib, const ROMol &queryMol, double simCutoff) {
  v2::FileParsers::SmilesMolSupplierParams fileparams;
  fileparams.titleLine = false;
  v2::FileParsers::SmilesMolSupplier suppl(enumLib, fileparams);
  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGen;
  fpGen.reset(MorganFingerprint::getMorganGenerator<std::uint64_t>(2));
  auto queryFP =
      std::make_unique<ExplicitBitVect>(*fpGen->getFingerprint(queryMol));
  std::set<std::string> fullSmi;
  std::set<std::string> names;
  while (!suppl.atEnd()) {
    auto mol = suppl.next();
    auto fp = fpGen->getFingerprint(*mol);
    auto sim = TanimotoSimilarity(*fp, *queryFP);
    if (sim >= simCutoff) {
      fullSmi.insert(MolToSmiles(*mol));
      names.insert(mol->getProp<std::string>("_Name"));
    }
  }
  return std::make_pair(fullSmi, names);
}

TEST_CASE("Small tests") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string fullRoot(fName + "/Code/GraphMol/SynthonSpaceSearch/data/");
  std::vector<std::string> libNames{
      fullRoot + "amide_space.txt",
      fullRoot + "triazole_space.txt",
      fullRoot + "urea_space.txt",
  };

  std::vector<std::string> enumLibNames{
      fullRoot + "amide_space_enum.smi",
      fullRoot + "triazole_space_enum.smi",
      fullRoot + "urea_space_enum.smi",
  };
  std::vector<std::string> querySmis{
      "c1ccccc1C(=O)N1CCCC1",
      "CC1CCN(c2nnc(CO)n2C2CCCC2)C1",
      "C[C@@H]1CC(NC(=O)NC2COC2)CN(C(=O)c2nccnc2F)C1",
  };

  std::vector<size_t> expNumHits{2, 4, 4};

  for (size_t i = 0; i < libNames.size(); i++) {
    SynthonSpace synthonspace;
    synthonspace.readTextFile(libNames[i]);
    SynthonSpaceSearchParams params;
    params.maxBondSplits = 2;
    auto queryMol = v2::SmilesParse::MolFromSmiles(querySmis[i]);
    auto results = synthonspace.fingerprintSearch(*queryMol, params);
    CHECK(results.getHitMolecules().size() == expNumHits[i]);
    std::set<std::string> resSmi;
    for (const auto &r : results.getHitMolecules()) {
      resSmi.insert(MolToSmiles(*r));
      std::cout << MolToSmiles(*r) << std::endl;
    }

    // Do the enumerated library, just to check
    auto [fullSmi, names] =
        bruteForceSearch(enumLibNames[i], *queryMol, params.similarityCutoff);
    CHECK(resSmi == fullSmi);
  }
}

TEST_CASE("Biggy") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.csv";
  std::string enumLib =
      "/Users/david/Projects/FreedomSpace/Syntons_5567_space_a_enum.smi";
  SynthonSpace synthonspace;
  synthonspace.readTextFile(libName);

  const std::vector<std::string> smis{
      "c1ccccc1C(=O)N1CCCC1", "c1ccccc1NC(=O)C1CCN1",
      "c12ccccc1c(N)nc(N)n2", "c12ccc(C)cc1[nH]nc2C(=O)NCc1cncs1",
      "c1n[nH]cn1",           "C(=O)NC(CC)C(=O)N(CC)C"};
  const std::vector<size_t> numRes{46, 2, 48892, 130, 0, 0};
  const std::vector<size_t> maxRes{281273, 194725, 48893, 130, 427, 32817};
  SynthonSpaceSearchParams params;
  params.maxHits = -1;
  for (size_t i = 0; i < smis.size(); ++i) {
    auto queryMol = v2::SmilesParse::MolFromSmiles(smis[i]);
    auto results = synthonspace.fingerprintSearch(*queryMol, params);
    CHECK(results.getHitMolecules().size() == numRes[i]);
    CHECK(results.getMaxNumResults() == maxRes[i]);
    auto [fullSmi, names] =
        bruteForceSearch(enumLib, *queryMol, params.similarityCutoff);
    std::cout << "hits " << results.getHitMolecules().size() << " vs "
              << fullSmi.size() << std::endl;
    CHECK(results.getHitMolecules().size() == fullSmi.size());
  }
}

TEST_CASE("Random Hits") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.csv";
  SynthonSpace synthonspace;
  synthonspace.readTextFile(libName);

  auto queryMol = "c12ccc(C)cc1[nH]nc2C(=O)NCc1cncs1"_smiles;
  SynthonSpaceSearchParams params;
  params.maxBondSplits = 4;
  params.maxHits = 100;
  params.randomSample = true;
  params.randomSeed = 1;
  auto results = synthonspace.fingerprintSearch(*queryMol, params);
  std::map<std::string, int> libCounts;
  for (const auto &m : results.getHitMolecules()) {
    std::string lib(m->getProp<std::string>("_Name").substr(0, 2));
    if (const auto &c = libCounts.find(lib); c == libCounts.end()) {
      libCounts.insert(std::make_pair(lib, 1));
    } else {
      c->second++;
    }
  }
  CHECK(results.getHitMolecules().size() == 100);
  std::map<std::string, int> expCounts{{"a1", 100}};
  CHECK(expCounts == libCounts);
}
