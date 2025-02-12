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
#include <fstream>

#include <GraphMol/SubstructLibrary/SubstructLibrary.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/Fingerprints/RDKitFPGenerator.h>
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
    synthonspace.readTextFile(libNames[i]);
    SynthonSpaceSearchParams params;
    params.maxBondSplits = 3;
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

TEST_CASE("FP Biggy") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.csv";
  SynthonSpace synthonspace;
  synthonspace.readTextFile(libName);

  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGen(
      MorganFingerprint::getMorganGenerator<std::uint64_t>(2));

  std::map<std::string, std::unique_ptr<RWMol>> mols;

  const std::vector<std::string> smis{
      "c1ccccc1C(=O)N1CCCC1", "c1ccccc1NC(=O)C1CCN1",
      "c12ccccc1c(N)nc(N)n2", "c12ccc(C)cc1[nH]nc2C(=O)NCc1cncs1",
      "c1n[nH]cn1",           "C(=O)NC(CC)C(=O)N(CC)C"};
  const std::vector<size_t> numRes{46, 2, 0, 123, 0, 0};
  const std::vector<size_t> maxRes{2408, 197, 0, 833, 0, 4};
  SynthonSpaceSearchParams params;
  params.approxSimilarityAdjuster = 0.2;
  params.maxHits = -1;
  for (size_t i = 0; i < smis.size(); ++i) {
    auto queryMol = v2::SmilesParse::MolFromSmiles(smis[i]);
    auto results = synthonspace.fingerprintSearch(*queryMol, *fpGen, params);
    CHECK(results.getHitMolecules().size() == numRes[i]);
    CHECK(results.getMaxNumResults() == maxRes[i]);
  }
}

TEST_CASE("FP Random Hits") {
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
  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGen(
      MorganFingerprint::getMorganGenerator<std::uint64_t>(2));
  auto results = synthonspace.fingerprintSearch(*queryMol, *fpGen, params);
  std::map<std::string, int> libCounts;
  for (const auto &m : results.getHitMolecules()) {
    std::string lib(
        m->getProp<std::string>(common_properties::_Name).substr(0, 2));
    if (const auto &c = libCounts.find(lib); c == libCounts.end()) {
      libCounts.insert(std::make_pair(lib, 1));
    } else {
      c->second++;
    }
  }
  CHECK(results.getHitMolecules().size() == 100);
  std::map<std::string, int> expCounts{{"a1", 100}};
  CHECK(expCounts == libCounts);
  CHECK(results.getHitMolecules().front()->getProp<double>("Similarity") ==
        Catch::Approx(0.711538));
  CHECK(results.getHitMolecules().back()->getProp<double>("Similarity") ==
        Catch::Approx(0.5));
}

TEST_CASE("Other Fingerprints") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.csv";
  SynthonSpace synthonspace;
  synthonspace.readTextFile(libName);
  SynthonSpaceSearchParams params;
  params.maxBondSplits = 3;
  params.maxHits = 100;
  params.randomSample = true;
  params.randomSeed = 1;
  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGen(
      MorganFingerprint::getMorganGenerator<std::uint64_t>(3));

  // This is Osimertinib (Tagrisso), which doesn't have any hits in this little
  // database.
  auto queryMol =
      "C=CC(=O)Nc1cc(Nc2nccc(-c3cn(C)c4ccccc34)n2)c(OC)cc1N(C)CCN(C)C"_smiles;
  auto results = synthonspace.fingerprintSearch(*queryMol, *fpGen, params);
  CHECK(results.getHitMolecules().empty());
}

TEST_CASE("Timeout") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.csv";
  SynthonSpace synthonspace;
  synthonspace.readTextFile(libName);
  SynthonSpaceSearchParams params;
  params.maxBondSplits = 3;
  params.maxHits = -1;
  params.similarityCutoff = 0.3;
  params.fragSimilarityAdjuster = 0.3;
  params.timeOut = 2;
  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGen(
      MorganFingerprint::getMorganGenerator<std::uint64_t>(3));

  auto queryMol = "c12ccc(C)cc1[nH]nc2C(=O)NCc1cncs1"_smiles;
  auto results = synthonspace.fingerprintSearch(*queryMol, *fpGen, params);
  CHECK(results.getTimedOut());

  // Make sure no timeout also works, but only on a short search.
  params.maxHits = 100;
  params.similarityCutoff = 0.3;
  params.fragSimilarityAdjuster = 0.2;
  params.timeOut = 0;
  auto results1 = synthonspace.fingerprintSearch(*queryMol, *fpGen, params);
  CHECK(!results1.getTimedOut());
}

TEST_CASE("FP Approx Similarity") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.csv";
  SynthonSpace synthonspace;
  synthonspace.readTextFile(libName);
  SynthonSpaceSearchParams params;
  // The addFP and subtractFP are built from a random selection of
  // products so do occasionally vary, so use a fixed seed.
  params.randomSeed = 1;
  params.similarityCutoff = 0.5;
  params.timeOut = 0;
  params.maxHits = 1000;

  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGen(
      RDKitFP::getRDKitFPGenerator<std::uint64_t>(3));
  auto queryMol = "c12ccc(C)cc1[nH]nc2C(=O)NCc1cncs1"_smiles;

  // With RDKit fingerprints, 0.05 gives a reasonable compromise
  // between speed and hits missed.
  params.approxSimilarityAdjuster = 0.05;
  auto results = synthonspace.fingerprintSearch(*queryMol, *fpGen, params);
  CHECK(results.getHitMolecules().size() == 482);
  CHECK(results.getMaxNumResults() == 1466);

  // A tighter adjuster misses more hits.
  params.approxSimilarityAdjuster = 0.01;
  results = synthonspace.fingerprintSearch(*queryMol, *fpGen, params);
  CHECK(results.getHitMolecules().size() == 124);

  // This is the actual number of hits achievable.
  params.approxSimilarityAdjuster = 0.25;
  results = synthonspace.fingerprintSearch(*queryMol, *fpGen, params);
  CHECK(results.getHitMolecules().size() == 914);
}
