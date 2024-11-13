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
    std::vector<std::unique_ptr<ExplicitBitVect>> fps;
    std::vector<std::unique_ptr<ROMol>> mols;
    v2::FileParsers::SmilesMolSupplierParams fileparams;
    fileparams.titleLine = false;
    v2::FileParsers::SmilesMolSupplier suppl(enumLibNames[i], fileparams);
    std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGen;
    fpGen.reset(MorganFingerprint::getMorganGenerator<std::uint64_t>(2));
    while (!suppl.atEnd()) {
      mols.emplace_back(suppl.next());
      fps.emplace_back(fpGen->getFingerprint(*mols.back()));
    }
    auto queryFP =
        std::make_unique<ExplicitBitVect>(*fpGen->getFingerprint(*queryMol));
    std::set<std::string> fullSmi;
    for (size_t i = 0; i < fps.size(); i++) {
      auto sim = TanimotoSimilarity(*fps[i], *queryFP);
      if (sim >= params.similarityCutoff) {
        fullSmi.insert(MolToSmiles(*mols[i]));
      }
    }
    CHECK(resSmi == fullSmi);
  }
}

TEST_CASE("Biggy") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.csv";
  SynthonSpace synthonspace;
  synthonspace.readTextFile(libName);

  const std::vector<std::string> smis{"c1ccccc1C(=O)N1CCCC1",
                                      "c1ccccc1NC(=O)C1CCN1",
                                      "c12ccccc1c(N)nc(N)n2",
                                      "c12ccc(C)cc1[nH]nc2C(=O)NCc1cncs1",
                                      "c1nncn1",
                                      "C(=O)NC(CC)C(=O)N(CC)C"};
  const std::vector<size_t> numRes{6785, 4544, 48892, 130, 29147, 5651};
  const std::vector<size_t> maxRes{6785, 4544, 48893, 130, 29312, 5869};
  SynthonSpaceSearchParams params;
  params.maxHits = -1;
  for (size_t i = 0; i < smis.size(); ++i) {
    if (i != 3) {
      continue;
    }
    auto queryMol = v2::SmilesParse::MolFromSmiles(smis[i]);
    auto results = synthonspace.fingerprintSearch(*queryMol, params);
    CHECK(results.getHitMolecules().size() == numRes[i]);
    CHECK(results.getMaxNumResults() == maxRes[i]);
  }
}
