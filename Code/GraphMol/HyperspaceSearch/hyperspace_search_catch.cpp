//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#include <GraphMol/SubstructLibrary/SubstructLibrary.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <catch2/catch_all.hpp>

#include "Hyperspace.h"
#include "HyperspaceSubstructureSearch.h"

using namespace RDKit;
using namespace RDKit::HyperspaceSSSearch;

namespace RDKit::HyperspaceSSSearch::details {
std::vector<std::vector<std::unique_ptr<ROMol>>> splitMolecule(
    const ROMol &query, unsigned int maxBondSplits);
}  // namespace RDKit::HyperspaceSSSearch::details

using namespace RDKit::HyperspaceSSSearch::details;

std::string fName = getenv("RDBASE");
const std::string LIB_NAME =
    fName + "/Code/GraphMol/HyperspaceSearch/data/idorsia_toy_space_a.txt";
const std::string ENUM_LIB_NAME =
    fName + "/Code/GraphMol/HyperspaceSearch/data/idorsia_toy_space_enum.smi";

std::unique_ptr<SubstructLibrary> loadSubstructLibrary(
    const std::string &smiFile = ENUM_LIB_NAME) {
  std::unique_ptr<SubstructLibrary> subsLib(new SubstructLibrary());
  v2::FileParsers::SmilesMolSupplierParams params;
  params.titleLine = false;
  v2::FileParsers::SmilesMolSupplier suppl(smiFile, params);
  while (!suppl.atEnd()) {
    auto mol = suppl.next();
    subsLib->addMol(*mol.release());
  }
  return subsLib;
}

TEST_CASE("Test splits 1", "[Test splits 1]") {
  std::vector<std::string> smiles{"c1ccccc1CN1CCN(CC1)C(-O)c1ncc(F)cc1",
                                  "CC(C)OCc1nnc(N2CC(C)CC2)n1C1CCCC1"};
  std::vector<std::vector<size_t>> expCounts{{6, 60, 350}, {8, 58, 326}};
  for (size_t i = 0; i < smiles.size(); ++i) {
    auto mol = v2::SmilesParse::MolFromSmiles(smiles[i]);
    REQUIRE(mol);
    auto fragments = splitMolecule(*mol, 3);
    REQUIRE(fragments.size() == 3);
    CHECK(fragments[0].size() == expCounts[i][0]);
    CHECK(fragments[1].size() == expCounts[i][1]);
    CHECK(fragments[2].size() == expCounts[i][2]);
  }
}

TEST_CASE("Amide 1", "[Amide 1]") {
  std::string libName =
      fName + "/Code/GraphMol/HyperspaceSearch/data/amide_space.txt";
  std::string enumLibName =
      fName + "/Code/GraphMol/HyperspaceSearch/data/amide_space_enum.smi";

  auto queryMol = "c1ccccc1C(=O)N1CCCC1"_smiles;
  auto results = SSSearch(*queryMol, 2, libName);
  CHECK(results.size() == 2);
  std::set<std::string> resSmi;
  for (const auto &r : results) {
    resSmi.insert(MolToSmiles(*r));
  }

  auto subsLib = loadSubstructLibrary(enumLibName);
  auto query = "c1ccccc1C(=O)N1CCCC1"_smarts;
  auto enumRes = subsLib->getMatches(*query);
  std::cout << "Number of enum results : " << enumRes.size() << " from "
            << subsLib->size() << " mols" << std::endl;
  std::set<std::string> enumSmi;
  for (auto i : enumRes) {
    std::cout << i << " : " << MolToSmiles(*subsLib->getMol(i)) << " : "
              << subsLib->getMol(i)->getProp<std::string>("_Name") << std::endl;
    enumSmi.insert(MolToSmiles(*subsLib->getMol(i)));
  }
  CHECK(resSmi == enumSmi);
}

TEST_CASE("Urea 1", "[Urea 1]") {
  std::string libName =
      fName + "/Code/GraphMol/HyperspaceSearch/data/urea_space.txt";
  Hyperspace hyperspace(libName);
  SECTION("Single fragged molecule") {
    auto fraggedMol =
        "O=C(NC1COC1)[1*].O=C(Nc1c(CN[1*])cc[s]1)[2*].Fc1nccnc1[2*]"_smiles;
    auto results = hyperspace.searchFragSet(*fraggedMol);
    CHECK(results.size() == 1);
  }
#if 1
  SECTION("Single molecule with fragging") {
    auto queryMol = "O=C(Nc1c(CNC=O)cc[s]1)c1nccnc1"_smiles;
    auto results = SSSearch(*queryMol, 3, hyperspace);
    CHECK(results.size() == 2);
  }
#endif
}

TEST_CASE("Simple query 1", "[Simple query 1]") {
  std::string libName =
      fName + "/Code/GraphMol/HyperspaceSearch/data/urea_3.txt";
  Hyperspace hyperspace(LIB_NAME);
  SECTION("Single fragged molecule") {
    auto fraggedMol = "c1ccccc1[1*].C1CCCN1C(=O)[1*]"_smiles;
    auto results = hyperspace.searchFragSet(*fraggedMol);
    CHECK(results.size() == 220);
  }
#if 1
  SECTION("Single molecule with fragging") {
    // should give 220 hits for urea-3
    auto queryMol = "c1ccccc1C(=O)N1CCCC1"_smiles;
    auto results = SSSearch(*queryMol, 3, hyperspace);
    CHECK(results.size() == 220);
  }
#endif
#if 0
  auto subsLib = loadSubstructLibrary();
  auto query = "c1ccccc1C(=O)N1CCCC1"_smarts;
  auto enumRes = subsLib->getMatches(*query);
  std::cout << "Number of enum results : " << enumRes.size() << " from "
            << subsLib->size() << " mols" << std::endl;
  for (auto i : enumRes) {
    std::cout << i << " : " << MolToSmiles(*subsLib->getMol(i)) << " : "
              << subsLib->getMol(i)->getProp<std::string>("_Name") << std::endl;
  }
#endif
}
