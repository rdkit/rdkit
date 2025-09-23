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
#include <GraphMol/RascalMCES/RascalMCES.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>
#include <GraphMol/SynthonSpaceSearch/SearchResults.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <catch2/catch_all.hpp>

using namespace RDKit;
using namespace RDKit::SynthonSpaceSearch;
using namespace RDKit::RascalMCES;

const char *rdbase = getenv("RDBASE");

void getMols(const std::string &molFilename,
             std::map<std::string, std::unique_ptr<RWMol>> &mols) {
  v2::FileParsers::SmilesMolSupplierParams fileparams;
  fileparams.titleLine = false;
  v2::FileParsers::SmilesMolSupplier suppl(molFilename, fileparams);
  std::map<std::string, std::unique_ptr<ExplicitBitVect>> fps;
  while (!suppl.atEnd()) {
    auto mol = suppl.next();
    auto molName = mol->getProp<std::string>(common_properties::_Name);
    mols.insert(std::make_pair(molName, mol.release()));
  }
}

std::set<std::string> bruteForceSearch(
    const ROMol &queryMol,
    const std::map<std::string, std::unique_ptr<RWMol>> &mols,
    const RascalOptions &rascalOptions) {
  std::set<std::string> fullSmi;
  std::set<std::string> names;
  for (auto &[name, mol] : mols) {
    auto res = rascalMCES(queryMol, *mol, rascalOptions);
    if (!res.empty() &&
        res.front().getSimilarity() > rascalOptions.similarityThreshold) {
      mol->setProp<double>("Similarity", res.front().getSimilarity());
      names.insert(name);
    }
  }
  return names;
}

TEST_CASE("RASCAL Small tests") {
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

  std::vector<size_t> expNumHits{6, 4, 1};

  RascalOptions rascalOptions;

  for (size_t i = 0; i < libNames.size(); i++) {
    // if (i != 0) {
    // continue;
    // }
    SynthonSpace synthonspace;
    bool cancelled = false;
    synthonspace.readTextFile(libNames[i], cancelled);
    SynthonSpaceSearchParams params;
    auto queryMol = v2::SmilesParse::MolFromSmiles(querySmis[i]);

    auto results = synthonspace.rascalSearch(*queryMol, rascalOptions, params);
    CHECK(results.getHitMolecules().size() == expNumHits[i]);
    std::set<std::string> resSmis;
    for (const auto &r : results.getHitMolecules()) {
      resSmis.insert(MolToSmiles(*r));
    }

    // test with callback version
    std::set<std::string> cbSmis;
    auto cb = [&cbSmis](std::vector<std::unique_ptr<ROMol>>& results) {
      for (auto& r : results) cbSmis.insert(MolToSmiles(*r));
      return false;
    };
    synthonspace.rascalSearch(*queryMol, rascalOptions, cb, params);
    CHECK(resSmis == cbSmis);

    // Do the enumerated library, just to check
    std::map<std::string, std::unique_ptr<RWMol>> mols;
    getMols(enumLibNames[i], mols);
    auto names = bruteForceSearch(*queryMol, mols, rascalOptions);
    std::set<std::string> fullSmis;
    for (const auto &r : names) {
      fullSmis.insert(MolToSmiles(*mols[r]));
    }
    // As with fingerprints, we don't get all the hits with synthon search
    // that we would with a full search.
    for (const auto &rs : resSmis) {
      CHECK(fullSmis.find(rs) != fullSmis.end());
    }
  }
}
