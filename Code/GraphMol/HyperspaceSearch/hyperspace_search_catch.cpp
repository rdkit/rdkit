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
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/HyperspaceSearch/Hyperspace.h>
#include <GraphMol/HyperspaceSearch/HyperspaceSubstructureSearch.h>
#include <GraphMol/HyperspaceSearch/ReactionSet.h>
#include <GraphMol/HyperspaceSearch/Reagent.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <catch2/catch_all.hpp>

using namespace RDKit;
using namespace RDKit::HyperspaceSSSearch;

namespace RDKit::HyperspaceSSSearch::details {
std::vector<std::vector<std::shared_ptr<ROMol>>> splitMolecule(
    const ROMol &query, unsigned int maxBondSplits);
}  // namespace RDKit::HyperspaceSSSearch::details

using namespace RDKit::HyperspaceSSSearch::details;

std::string fName = getenv("RDBASE");
const std::string TXT_LIB_NAME =
    fName + "/Code/GraphMol/HyperspaceSearch/data/idorsia_toy_space_a.txt";
const std::string BIN_LIB_NAME =
    fName + "/Code/GraphMol/HyperspaceSearch/data/idorsia_toy_space.spc";
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
  std::vector<std::vector<size_t>> expCounts{{1, 51, 345, 20},
                                             {1, 38, 298, 56}};
  for (size_t i = 0; i < smiles.size(); ++i) {
    auto mol = v2::SmilesParse::MolFromSmiles(smiles[i]);
    REQUIRE(mol);
    auto fragments = splitMolecule(*mol, 3);
    CHECK(fragments.size() ==
          std::reduce(expCounts[i].begin(), expCounts[i].end(), size_t(0)));
    // The first fragment set should just be the molecule itself.
    for (size_t j = 0; j < 4; ++j) {
      auto numFragSets = std::reduce(
          fragments.begin(), fragments.end(), size_t(0),
          [&](size_t prevRes, std::vector<std::shared_ptr<ROMol>> &frags) {
            if (frags.size() == j + 1) {
              return prevRes + 1;
            } else {
              return prevRes;
            }
          });
      CHECK(numFragSets == expCounts[i][j]);
    }
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
#if 1
  SECTION("Single fragged molecule") {
    std::vector<std::shared_ptr<ROMol>> fragSet{
        std::shared_ptr<ROMol>("O=C(NC1COC1)[1*]"_smiles),
        std::shared_ptr<ROMol>("O=C(Nc1c(CN[1*])cc[s]1)[2*]"_smiles),
        std::shared_ptr<ROMol>("Fc1nccnc1[2*]"_smiles),
    };
    auto results = hyperspace.searchFragSet(fragSet);
    CHECK(results.size() == 1);
  }
#endif
#if 1
  SECTION("Single molecule with fragging") {
    auto queryMol = "O=C(Nc1c(CNC=O)cc[s]1)c1nccnc1"_smiles;
    auto results = SSSearch(*queryMol, 3, hyperspace);
    CHECK(results.size() == 2);
  }
#endif
}

TEST_CASE("Simple query 1", "[Simple query 1]") {
  //  std::string libName =
  //      fName + "/Code/GraphMol/HyperspaceSearch/data/urea_3.txt";
  Hyperspace hyperspace(TXT_LIB_NAME);
  SECTION("Single fragged molecule") {
    std::vector<std::shared_ptr<ROMol>> fragSet{
        std::shared_ptr<ROMol>("c1ccccc1[1*]"_smiles),
        std::shared_ptr<ROMol>("C1CCCN1C(=O)[1*]"_smiles),
    };
    auto results = hyperspace.searchFragSet(fragSet);
    CHECK(results.size() == 220);
  }
  SECTION("Binary File") {
    Hyperspace hyperspace;
    hyperspace.readFromDBStream(BIN_LIB_NAME);
    // should give 220 hits for urea-3
    const auto start{std::chrono::steady_clock::now()};
    auto queryMol = "c1ccccc1C(=O)N1CCCC1"_smiles;
    auto results = SSSearch(*queryMol, 3, hyperspace);
    const auto end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_seconds{end - start};
    std::cout << "1 Elapsed time : " << elapsed_seconds.count() << std::endl;
    CHECK(results.size() == 220);
  }
#if 1
  SECTION("Single molecule with fragging") {
    {
      // should give 220 hits for urea-3
      const auto start{std::chrono::steady_clock::now()};
      auto queryMol = "c1ccccc1C(=O)N1CCCC1"_smiles;
      auto results = SSSearch(*queryMol, 3, hyperspace);
      const auto end{std::chrono::steady_clock::now()};
      const std::chrono::duration<double> elapsed_seconds{end - start};
      std::cout << "1 Elapsed time : " << elapsed_seconds.count() << std::endl;
      CHECK(results.size() == 220);
    }
    {
      const auto start{std::chrono::steady_clock::now()};
      auto queryMol = "O=C(Nc1c(CNC=O)cc[s]1)c1nccnc1"_smiles;
      auto results = SSSearch(*queryMol, 3, hyperspace);
      const auto end{std::chrono::steady_clock::now()};
      const std::chrono::duration<double> elapsed_seconds{end - start};
      std::cout << "2 Elapsed time : " << elapsed_seconds.count() << std::endl;
      CHECK(results.size() == 20);
    }
  }
#endif
#if 0
  SECTION("Brute Force") {
    const auto start{std::chrono::steady_clock::now()};
    auto subsLib = loadSubstructLibrary();
    auto query = "c1ccccc1C(=O)N1CCCC1"_smarts;
    auto enumRes = subsLib->getMatches(*query);
    const auto end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_seconds{end - start};
    CHECK(enumRes.size() == 220);
    std::cout << "3 Elapsed time : " << elapsed_seconds.count() << std::endl;
    std::cout << "Number of enum results : " << enumRes.size() << " from "
              << subsLib->size() << " mols" << std::endl;
    for (auto i : enumRes) {
      std::cout << i << " : " << MolToSmiles(*subsLib->getMol(i)) << " : "
                << subsLib->getMol(i)->getProp<std::string>("_Name")
                << std::endl;
    }
  }
#endif
}

TEST_CASE("Triazole", "[Triazole]") {
  std::string libName =
      fName + "/Code/GraphMol/HyperspaceSearch/data/triazole_space.txt";
  std::string enumLibName =
      fName + "/Code/GraphMol/HyperspaceSearch/data/triazole_space_enum.smi";

  Hyperspace hyperspace(libName);

#if 1
  SECTION("Fragged Mol") {
    auto queryMol =
        "OCC([1*])=NN=[2*].C1CCCC1N([3*])[1*].CC1CCN(C1)C(=[2*])[3*]"_smiles;
    REQUIRE(queryMol);
    std::vector<std::unique_ptr<ROMol>> tmpFrags;
    MolOps::getMolFrags(*queryMol, tmpFrags, false);
    std::vector<std::shared_ptr<ROMol>> queryFrags;
    std::transform(tmpFrags.begin(), tmpFrags.end(),
                   std::back_inserter(queryFrags),
                   [&](std::unique_ptr<ROMol> &m) -> std::shared_ptr<ROMol> {
                     return std::shared_ptr<ROMol>(m.release());
                   });
    auto results = hyperspace.searchFragSet(queryFrags);
    CHECK(results.size() == 4);
  }
#endif
#if 1
  SECTION("Full Molecule") {
    auto queryMol = "OCc1ncnn1"_smarts;
    REQUIRE(queryMol);
    auto results = SSSearch(*queryMol, 3, hyperspace);
    CHECK(results.size() == 8);
    std::set<std::string> resSmi;
    for (const auto &r : results) {
      resSmi.insert(MolToSmiles(*r));
    }

    auto subsLib = loadSubstructLibrary(enumLibName);
    auto enumRes = subsLib->getMatches(*queryMol);
    std::cout << "Number of enum results : " << enumRes.size() << " from "
              << subsLib->size() << " mols" << std::endl;
    std::set<std::string> enumSmi;
    for (auto i : enumRes) {
      std::cout << i << " : " << MolToSmiles(*subsLib->getMol(i)) << " : "
                << subsLib->getMol(i)->getProp<std::string>("_Name")
                << std::endl;
      enumSmi.insert(MolToSmiles(*subsLib->getMol(i)));
    }
    CHECK(resSmi == enumSmi);
  }
#endif
}

TEST_CASE("Quinoline", "[Quinoline]") {
  std::string libName =
      fName + "/Code/GraphMol/HyperspaceSearch/data/doebner_miller_space.txt";
  std::string enumLibName =
      fName +
      "/Code/GraphMol/HyperspaceSearch/data/doebner_miller_space_enum.smi";

  Hyperspace hyperspace(libName);
  {
    auto queryMol = "c1ccccn1"_smiles;
    auto results = SSSearch(*queryMol, 3, hyperspace);
    CHECK(results.size() == 12);
    std::set<std::string> resSmi;
    for (const auto &r : results) {
      std::cout << "Result : "
                << r->getProp<std::string>(common_properties::_Name) << " : "
                << MolToSmiles(*r) << std::endl;
      resSmi.insert(MolToSmiles(*r));
    }
    auto subsLib = loadSubstructLibrary(enumLibName);
    auto enumRes = subsLib->getMatches(*queryMol);
    std::cout << "Number of enum results : " << enumRes.size() << " from "
              << subsLib->size() << " mols" << std::endl;
    std::set<std::string> enumSmi;
    for (auto i : enumRes) {
      std::cout << i << " : " << MolToSmiles(*subsLib->getMol(i)) << " : "
                << subsLib->getMol(i)->getProp<std::string>("_Name")
                << std::endl;
      enumSmi.insert(MolToSmiles(*subsLib->getMol(i)));
    }
    CHECK(resSmi == enumSmi);
  }
}

TEST_CASE("Substructure in 1 reagent", "[Substructure in 1 reagent]") {
  // Making sure it works when the query is a complete substructure of 1
  // of the reagents in the library, so the whole library is a hit.
  {
    auto queryMol = "N1CCCC1"_smiles;
    std::string libName =
        fName + "/Code/GraphMol/HyperspaceSearch/data/triazole_space.txt";
    auto results = SSSearch(*queryMol, 3, libName);
    CHECK(results.size() == 8);
  }
  {
    auto queryMol = "N1CCC(C(F)(F)F)C1"_smiles;
    std::string libName =
        fName + "/Code/GraphMol/HyperspaceSearch/data/triazole_space.txt";
    auto results = SSSearch(*queryMol, 3, libName);
    CHECK(results.size() == 4);
  }
  {
    auto queryMol = "C1CCCCC1"_smiles;
    std::string libName =
        fName + "/Code/GraphMol/HyperspaceSearch/data/triazole_space.txt";
    auto results = SSSearch(*queryMol, 3, libName);
    CHECK(results.size() == 0);
  }
}

TEST_CASE("Connector Regions", "[Connector Regions]") {
  SECTION("Single tests") {
    auto m1 = "[1*]CN(C[2*])Cc1ccccc1"_smiles;
    REQUIRE(m1);
    CHECK(MolToSmiles(*getConnRegion(*m1)) == "[1*]CN(C)C[1*]");

    auto m2 = "[1*]CN(C[2*])Cc1ccc(CN(C[3*])C[1*])cc1"_smiles;
    REQUIRE(m2);
    CHECK(MolToSmiles(*getConnRegion(*m2)) == "[1*]CN(C)C[1*].[1*]CN(C)C[1*]");

    auto m3 = "[2*]C"_smiles;
    REQUIRE(m3);
    CHECK(MolToSmiles(*getConnRegion(*m3)) == "[1*]C");

    auto m4 = "[1*]c1cnccc1"_smiles;
    REQUIRE(m4);
    CHECK(MolToSmiles(*getConnRegion(*m4)) == "[1*]c(cc)cn");
  }

  SECTION("Built from file") {
    std::string libName =
        fName + "/Code/GraphMol/HyperspaceSearch/data/urea_3.txt";
    Hyperspace hyperspace(libName);
    const auto &rs = hyperspace.reactions().begin();
    CHECK(rs->second->connectorRegions().size() == 30);
  }
}

TEST_CASE("DB Writer", "[DB Writer]") {
  std::string libName =
      fName + "/Code/GraphMol/HyperspaceSearch/data/doebner_miller_space.txt";
  Hyperspace hyperspace(libName);
  CHECK(hyperspace.numReactions() == 1);
  hyperspace.writeToDBStream("doebner_miller_space.spc");

  Hyperspace newHyperspace;
  newHyperspace.readFromDBStream("doebner_miller_space.spc");
  auto it = newHyperspace.reactions().find("doebner-miller-quinoline");
  CHECK(it != newHyperspace.reactions().end());
  const auto &irxn = it->second;
  const auto &orxn =
      hyperspace.reactions().find("doebner-miller-quinoline")->second;
  CHECK(irxn->id() == orxn->id());
  CHECK(irxn->connectorRegions().size() == orxn->connectorRegions().size());
  CHECK(*irxn->connRegFP() == *orxn->connRegFP());
  CHECK(irxn->connectors() == orxn->connectors());
  CHECK(irxn->reagents().size() == orxn->reagents().size());
  for (size_t i = 0; i < irxn->reagents().size(); ++i) {
    CHECK(irxn->reagents()[i].size() == orxn->reagents()[i].size());
    for (size_t j = 0; j < irxn->reagents().size(); ++j) {
      CHECK(irxn->reagents()[i][j]->id() == orxn->reagents()[i][j]->id());
    }
  }
}
