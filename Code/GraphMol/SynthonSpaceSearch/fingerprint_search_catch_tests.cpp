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

TEST_CASE("Amide 1") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/amide_space.txt";
  std::string enumLibName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/amide_space_enum.smi";

  auto queryMol = "c1ccccc1C(=O)N1CCCC1"_smiles;
  SynthonSpace synthonspace;
  synthonspace.readTextFile(libName);
  SynthonSpaceSearchParams params;
  params.maxBondSplits = 2;
  auto results = synthonspace.fingerprintSearch(*queryMol, params);
  CHECK(results.getHitMolecules().size() == 2);
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
  v2::FileParsers::SmilesMolSupplier suppl(enumLibName, fileparams);
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
