//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>

#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/Fingerprints/RDKitFPGenerator.h>
#include <GraphMol/SubstructLibrary/SubstructLibrary.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>
#include <GraphMol/SynthonSpaceSearch/SearchResults.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <catch2/catch_all.hpp>

using namespace RDKit;
using namespace RDKit::SynthonSpaceSearch;

using namespace RDKit::SynthonSpaceSearch::details;

const char *rdbase = getenv("RDBASE");

std::map<std::string, std::unique_ptr<ExplicitBitVect>> getFingerprints(
    const std::string &molFilename,
    std::map<std::string, std::unique_ptr<RWMol>> &mols,
    const std::unique_ptr<FingerprintGenerator<std::uint64_t>> &fpGen) {
  v2::FileParsers::SmilesMolSupplierParams fileparams;
  fileparams.titleLine = false;
  v2::FileParsers::SmilesMolSupplier suppl(molFilename, fileparams);
  std::map<std::string, std::unique_ptr<ExplicitBitVect>> fps;
  while (!suppl.atEnd()) {
    auto mol = suppl.next();
    auto fp = fpGen->getFingerprint(*mol);
    fps.insert(std::make_pair(
        mol->getProp<std::string>(common_properties::_Name), fp));
    auto molName = mol->getProp<std::string>(common_properties::_Name);
    mols.insert(std::make_pair(molName, mol.release()));
  }

  return fps;
}

std::set<std::string> bruteForceSearch(
    std::map<std::string, std::unique_ptr<ExplicitBitVect>> &fps,
    const ROMol &queryMol, const double simCutoff) {
  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGen;
  fpGen.reset(MorganFingerprint::getMorganGenerator<std::uint64_t>(2));
  const auto queryFP =
      std::unique_ptr<ExplicitBitVect>(fpGen->getFingerprint(queryMol));
  std::set<std::string> fullSmi;
  std::set<std::string> names;
  for (auto &it : fps) {
    if (auto sim = TanimotoSimilarity(*it.second, *queryFP); sim >= simCutoff) {
      names.insert(it.first);
    }
  }
  return names;
}

TEST_CASE("FP Small tests") {
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

  std::vector<size_t> expNumHits{2, 3, 4};

  for (size_t i = 0; i < libNames.size(); i++) {
    SynthonSpace synthonspace;
    bool cancelled = false;
    synthonspace.readTextFile(libNames[i], cancelled);
    SynthonSpaceSearchParams params;
    params.randomSeed = 1;
    params.approxSimilarityAdjuster = 0.2;
    auto queryMol = v2::SmilesParse::MolFromSmiles(querySmis[i]);
    std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGen(
        MorganFingerprint::getMorganGenerator<std::uint64_t>(2));

    auto results = synthonspace.fingerprintSearch(*queryMol, *fpGen, params);
    CHECK(results.getHitMolecules().size() == expNumHits[i]);
    std::set<std::string> resSmis;
    for (const auto &r : results.getHitMolecules()) {
      resSmis.insert(MolToSmiles(*r));
    }

    // test with callback version
    std::set<std::string> cbSmis;
    auto cb = [&cbSmis](const std::vector<std::unique_ptr<ROMol>> &results) {
      for (auto &r : results) cbSmis.insert(MolToSmiles(*r));
      return false;
    };
    synthonspace.fingerprintSearch(*queryMol, *fpGen, cb, params);
    CHECK(resSmis == cbSmis);

    // Do the enumerated library, just to check
    std::map<std::string, std::unique_ptr<RWMol>> mols;
    auto fps = getFingerprints(enumLibNames[i], mols, fpGen);
    auto names = bruteForceSearch(fps, *queryMol, params.similarityCutoff);
    std::set<std::string> fullSmis;
    for (const auto &r : names) {
      fullSmis.insert(MolToSmiles(*mols[r]));
    }
    if (i != 1) {
      CHECK(resSmis == fullSmis);
    } else {
      // In the triazole library, one of the hits found by the brute force
      // method (triazole-1_1-1_2-2_3-1) is missed by the SynthonSpaceSearch
      // because it requires that the fragment [1*]n([3*])C1CCCC1 is similar
      // to synthon c1ccccc1-n([3*])[1*] which it isn't.  Instead, make sure
      // all the ones that are found are in the brute force results.
      for (const auto &rs : resSmis) {
        CHECK(fullSmis.find(rs) != fullSmis.end());
      }
    }
  }
}

TEST_CASE("FP Binary File") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  SynthonSpace synthonspace;
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/idorsia_toy_space_a.spc";
  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGen(
      RDKitFP::getRDKitFPGenerator<std::uint64_t>());
  SearchResults results;
  auto queryMol = "O=C(Nc1c(CNC=O)cc[s]1)c1nccnc1"_smiles;
  SynthonSpaceSearchParams params;
  for (auto numThreads : std::vector<int>{1, 2, -1}) {
    synthonspace.readDBFile(libName, numThreads);
    params.numThreads = numThreads;
    CHECK_NOTHROW(
        results = synthonspace.fingerprintSearch(*queryMol, *fpGen, params));
    CHECK(results.getHitMolecules().size() == 4);
    CHECK(results.getMaxNumResults() == 420);
  }

  // Make sure it rejects the wrong sort of fingerprint.
  synthonspace.readDBFile(libName);
  fpGen.reset(MorganFingerprint::getMorganGenerator<std::uint64_t>(2));
  CHECK_THROWS(results = synthonspace.fingerprintSearch(*queryMol, *fpGen));
}

TEST_CASE("Missing exact match") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  SynthonSpace synthonspace;
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/missing_hit.txt";
  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGen(
      RDKitFP::getRDKitFPGenerator<std::uint64_t>());
  SearchResults results;
  auto queryMol =
      "N#Cc1ccc(C(=O)Cn2cc(-c3cn[nH]c3-c3nc4ccccc4[nH]3)nn2)s1"_smiles;
  bool cancelled = false;
  synthonspace.readTextFile(libName, cancelled);
  synthonspace.buildSynthonFingerprints(*fpGen);
  CHECK_NOTHROW(results = synthonspace.fingerprintSearch(*queryMol, *fpGen));
  CHECK(results.getHitMolecules().size() == 1);
  CHECK(results.getHitMolecules()[0]->getProp<double>("Similarity") == 1.0);
}

TEST_CASE("Hit Filters") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  SynthonSpace synthonspace;
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/idorsia_toy_space_a.spc";
  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGen(
      RDKitFP::getRDKitFPGenerator<std::uint64_t>());
  SearchResults results;
  auto queryMol = "CCNC(=O)Cc1cncc(CCOC2c3ccccc3CC2)c1"_smiles;
  SynthonSpaceSearchParams params;
  params.similarityCutoff = 0.45;
  synthonspace.readDBFile(libName);
  results = synthonspace.fingerprintSearch(*queryMol, *fpGen, params);
  CHECK(results.getHitMolecules().size() == 18);
  {
    SynthonSpaceSearchParams params;
    params.minHitHeavyAtoms = 28;
    params.similarityCutoff = 0.45;
    results = synthonspace.fingerprintSearch(*queryMol, *fpGen, params);
    CHECK(results.getHitMolecules().size() == 13);
    params.maxHitHeavyAtoms = 29;
    results = synthonspace.fingerprintSearch(*queryMol, *fpGen, params);
    CHECK(results.getHitMolecules().size() == 12);
    for (const auto &r : results.getHitMolecules()) {
      auto numHeavies = Descriptors::calcNumHeavyAtoms(*r);
      CHECK((numHeavies == 28 || numHeavies == 29));
    }
  }
  {
    SynthonSpaceSearchParams params;
    params.similarityCutoff = 0.45;
    params.minHitMolWt = 375.0;
    results = synthonspace.fingerprintSearch(*queryMol, *fpGen, params);
    CHECK(results.getHitMolecules().size() == 13);
    params.maxHitMolWt = 390.0;
    results = synthonspace.fingerprintSearch(*queryMol, *fpGen, params);
    CHECK(results.getHitMolecules().size() == 4);
    for (const auto &r : results.getHitMolecules()) {
      auto molWt = Descriptors::calcExactMW(*r);
      CHECK((molWt >= 375.0 || molWt <= 390.0));
    }
  }
  {
    SynthonSpaceSearchParams params;
    params.similarityCutoff = 0.45;
    auto chiralQuery = "Cc1nccn1CCc1ccsc1COO[C@@H]1CCC[C@H](N)C1"_smiles;
    results = synthonspace.fingerprintSearch(*chiralQuery, *fpGen, params);
    CHECK(results.getHitMolecules().size() == 17);
    params.minHitChiralAtoms = 1;
    results = synthonspace.fingerprintSearch(*chiralQuery, *fpGen, params);
    CHECK(results.getHitMolecules().size() == 11);
    params.maxHitChiralAtoms = 1;
    results = synthonspace.fingerprintSearch(*chiralQuery, *fpGen, params);
    CHECK(results.getHitMolecules().size() == 4);
    for (const auto &r : results.getHitMolecules()) {
      auto numChiralAtoms = details::countChiralAtoms(*r);
      CHECK((numChiralAtoms >= 1 && numChiralAtoms <= 3));
    }
  }
}
