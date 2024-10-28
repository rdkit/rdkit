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
#include <GraphMol/HyperspaceSearch/Reagent.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <catch2/catch_all.hpp>

using namespace RDKit;
using namespace RDKit::HyperspaceSearch;

namespace RDKit::HyperspaceSearch::details {
std::vector<std::vector<std::unique_ptr<ROMol>>> splitMolecule(
    const ROMol &query, unsigned int maxBondSplits);
}  // namespace RDKit::HyperspaceSearch::details

using namespace RDKit::HyperspaceSearch::details;

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
  std::cout << "Reading " << smiFile << std::endl;
  int mol_num = 0;
  while (!suppl.atEnd()) {
    //    auto mol = suppl.next();
    subsLib->addMol(*suppl.next());
    ++mol_num;
    if (mol_num && !(mol_num % 50000)) {
      std::cout << "Read " << mol_num << " molecules." << std::endl;
    }
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
          std::accumulate(expCounts[i].begin(), expCounts[i].end(), size_t(0)));
    // The first fragment set should just be the molecule itself.
    for (size_t j = 0; j < 4; ++j) {
      auto numFragSets = std::accumulate(
          fragments.begin(), fragments.end(), size_t(0),
          [&](size_t prevRes, std::vector<std::unique_ptr<ROMol>> &frags) {
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
  Hyperspace hyperspace;
  hyperspace.readTextFile(libName);
  HyperspaceSearchParams params;
  params.maxBondSplits = 2;
  auto results = hyperspace.substructureSearch(*queryMol, params);
  CHECK(results.hitMolecules().size() == 2);
  std::set<std::string> resSmi;
  for (const auto &r : results.hitMolecules()) {
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
  Hyperspace hyperspace;
  hyperspace.readTextFile(libName);
#if 1
  SECTION("Single fragged molecule") {
    std::vector<std::unique_ptr<ROMol>> fragSet;
    fragSet.emplace_back("O=C(NC1COC1)[1*]"_smiles);
    fragSet.emplace_back("O=C(Nc1c(CN[1*])cc[s]1)[2*]"_smiles);
    fragSet.emplace_back("Fc1nccnc1[2*]"_smiles);
    auto results = hyperspace.searchFragSet(fragSet);
    CHECK(results.size() == 1);
  }
#endif
#if 1
  SECTION("Single molecule with fragging") {
    auto queryMol = "O=C(Nc1c(CNC=O)cc[s]1)c1nccnc1"_smiles;
    auto results = hyperspace.substructureSearch(*queryMol);
    CHECK(results.hitMolecules().size() == 2);
  }
#endif
}

TEST_CASE("Simple query 1", "[Simple query 1]") {
  std::string libName =
      fName + "/Code/GraphMol/HyperspaceSearch/data/urea_3.txt";
  Hyperspace hyperspace;
  hyperspace.readTextFile(libName);
//  hyperspace.readTextFile(TXT_LIB_NAME);
#if 1
  SECTION("Single fragged molecule") {
    std::vector<std::unique_ptr<ROMol>> fragSet;
    fragSet.emplace_back("c1ccccc1[1*]"_smiles);
    fragSet.emplace_back("C1CCCN1C(=O)[1*]"_smiles);
    auto results = hyperspace.searchFragSet(fragSet);
    CHECK(results.size() == 1);
    CHECK(results.front().numHits == 220);
  }
#endif
#if 1
  SECTION("Binary File") {
    Hyperspace hyperspace;
    hyperspace.readDBFile(BIN_LIB_NAME);
    // should give 220 hits for urea-3
    const auto start{std::chrono::steady_clock::now()};
    auto queryMol = "c1ccccc1C(=O)N1CCCC1"_smiles;
    auto results = hyperspace.substructureSearch(*queryMol);
    const auto end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_seconds{end - start};
    std::cout << "1 Elapsed time : " << elapsed_seconds.count() << std::endl;
    CHECK(results.hitMolecules().size() == 220);
    CHECK(results.maxNumResults() == 220);
  }
#endif
#if 1
  SECTION("Single molecule with fragging") {
    {
      // should give 220 hits for urea-3
      const auto start{std::chrono::steady_clock::now()};
      auto queryMol = "c1ccccc1C(=O)N1CCCC1"_smiles;
      auto results = hyperspace.substructureSearch(*queryMol);
      const auto end{std::chrono::steady_clock::now()};
      const std::chrono::duration<double> elapsed_seconds{end - start};
      std::cout << "1 Elapsed time : " << elapsed_seconds.count() << std::endl;
      CHECK(results.hitMolecules().size() == 220);
      CHECK(results.maxNumResults() == 220);
    }
    {
      const auto start{std::chrono::steady_clock::now()};
      auto queryMol = "O=C(Nc1c(CNC=O)cc[s]1)c1nccnc1"_smiles;
      auto results = hyperspace.substructureSearch(*queryMol);
      const auto end{std::chrono::steady_clock::now()};
      const std::chrono::duration<double> elapsed_seconds{end - start};
      std::cout << "2 Elapsed time : " << elapsed_seconds.count() << std::endl;
      CHECK(results.hitMolecules().size() == 20);
    }
  }
#endif
}

TEST_CASE("Triazole", "[Triazole]") {
  std::string libName =
      fName + "/Code/GraphMol/HyperspaceSearch/data/triazole_space.txt";
  std::string enumLibName =
      fName + "/Code/GraphMol/HyperspaceSearch/data/triazole_space_enum.smi";

  Hyperspace hyperspace;
  hyperspace.readTextFile(libName);

#if 1
  SECTION("Fragged Mol") {
    auto queryMol =
        "OCC([1*])=NN=[2*].C1CCCC1N([3*])[1*].CC1CCN(C1)C(=[2*])[3*]"_smiles;
    REQUIRE(queryMol);
    std::vector<std::unique_ptr<ROMol>> queryFrags;
    MolOps::getMolFrags(*queryMol, queryFrags, false);
    auto results = hyperspace.searchFragSet(queryFrags);
    CHECK(results.size() == 1);
    CHECK(results.front().numHits == 4);
  }
#endif
#if 1
  SECTION("Full Molecule") {
    auto queryMol = "OCc1ncnn1"_smarts;
    REQUIRE(queryMol);
    auto results = hyperspace.substructureSearch(*queryMol);
    CHECK(results.hitMolecules().size() == 8);
    std::set<std::string> resSmi;
    for (const auto &r : results.hitMolecules()) {
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

  Hyperspace hyperspace;
  hyperspace.readTextFile(libName);
  {
    auto queryMol = "c1ccccn1"_smiles;
    auto results = hyperspace.substructureSearch(*queryMol);
    CHECK(results.hitMolecules().size() == 12);
    std::set<std::string> resSmi;
    for (const auto &r : results.hitMolecules()) {
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
  std::string libName =
      fName + "/Code/GraphMol/HyperspaceSearch/data/triazole_space.txt";
  Hyperspace hyperspace;
  hyperspace.readTextFile(libName);
  {
    auto queryMol = "N1CCCC1"_smiles;
    auto results = hyperspace.substructureSearch(*queryMol);
    CHECK(results.hitMolecules().size() == 8);
  }
  {
    auto queryMol = "N1CCC(C(F)(F)F)C1"_smiles;
    auto results = hyperspace.substructureSearch(*queryMol);
    CHECK(results.hitMolecules().size() == 4);
  }
  {
    auto queryMol = "C1CCCCC1"_smiles;
    auto results = hyperspace.substructureSearch(*queryMol);
    CHECK(results.hitMolecules().size() == 0);
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
    Hyperspace hyperspace;
    hyperspace.readTextFile(libName);
    const auto &rs = hyperspace.reactions().begin();
    CHECK(rs->second->connectorRegions().size() == 30);
  }
}

TEST_CASE("DB Writer", "[DB Writer]") {
  std::string libName =
      fName + "/Code/GraphMol/HyperspaceSearch/data/doebner_miller_space.txt";
  Hyperspace hyperspace;
  hyperspace.readTextFile(libName);
  CHECK(hyperspace.numReactions() == 1);
  hyperspace.writeDBFile("doebner_miller_space.spc");

  Hyperspace newHyperspace;
  newHyperspace.readDBFile("doebner_miller_space.spc");
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

TEST_CASE("Biggy", "[Biggy]") {
  std::string libName = "/Users/david/Projects/FreedomSpace/Syntons_5567.spc";
  std::string enumLibName =
      "/Users/david/Projects/FreedomSpace/Syntons_5567_space_a_enum.smi";
  Hyperspace hyperspace;
  hyperspace.readDBFile(libName);

#if 0
  SECTION("Fragged Mol") {
    std::vector<std::unique_ptr<ROMol>> fragSet;
    fragSet.emplace_back("c1ccccc1N[1*]"_smiles);
    fragSet.emplace_back("N1CCC1C(=O)[1*]"_smiles);
    auto results = hyperspace.searchFragSet(fragSet);
    CHECK(results.size() == 1);
  }
#endif

#if 1
  SECTION("WholeMol") {
    const std::vector<std::string> smis{"c1ccccc1C(=O)N1CCCC1",
                                        "c1ccccc1NC(=O)C1CCN1",
                                        "c12ccccc1c(N)nc(N)n2",
                                        "c12ccc(C)cc1[nH]nc2C(=O)NCc1cncs1",
                                        "c1nncn1",
                                        "C(=O)NC(CC)C(=O)N(CC)C"};
    const std::vector<size_t> numRes{6785, 4544, 48892, 1, 29147, 5651};
    const std::vector<size_t> maxRes{6785, 4544, 48893, 1, 29312, 5869};
    HyperspaceSearchParams params;
    params.maxHits = -1;
    for (size_t i = 0; i < smis.size(); ++i) {
      //      if (i != 2) {
      //        continue;
      //      }
      std::cout << "TTTTTTTTTTTTTTTTTTTTTTTTTTT : " << smis[i] << std::endl;
      auto queryMol = v2::SmilesParse::MolFromSmarts(smis[i]);
      const auto start{std::chrono::steady_clock::now()};
      auto results = hyperspace.substructureSearch(*queryMol, params);
      const auto end{std::chrono::steady_clock::now()};
      const std::chrono::duration<double> elapsed_seconds{end - start};
      std::cout << "Elapsed time : " << elapsed_seconds.count() << std::endl;
      std::cout << "Number of hits : " << results.hitMolecules().size()
                << std::endl;
      std::string outFile =
          std::string("/Users/david/Projects/FreedomSpace/hyperspace_hits_") +
          std::to_string(i) + ".smi";
      std::ofstream of(outFile);
      for (const auto &r : results.hitMolecules()) {
        of << MolToSmiles(*r) << " " << r->getProp<std::string>("_Name")
           << std::endl;
      }
      CHECK(results.hitMolecules().size() == numRes[i]);
      CHECK(results.maxNumResults() == maxRes[i]);
    }
  }
#endif
#if 0
  SECTION("Brute Force") {
    auto queryMol = "c1ccccc1C(=O)N1CCCC1"_smiles;
    const auto start{std::chrono::steady_clock::now()};
    auto subsLib = loadSubstructLibrary(enumLibName);
    std::cout << "Loaded SubstructLibrary" << std::endl;
    const auto iend{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> ielapsed_seconds{iend - start};
    std::cout << "Elapsed time : " << ielapsed_seconds.count() << std::endl;
    auto enumRes = subsLib->getMatches(*queryMol);
    std::cout << "Number of enum results : " << enumRes.size() << " from "
              << subsLib->size() << " mols" << std::endl;
    const auto end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_seconds{end - start};
    std::cout << "Elapsed time : " << elapsed_seconds.count() << std::endl;
  }
#endif
}

TEST_CASE("FreedomSpace", "[FreedomSpace]") {
  std::string libName =
      "/Users/david/Projects/FreedomSpace/2023-05_Freedom_synthons.spc";
  Hyperspace hyperspace;
  const auto rstart{std::chrono::steady_clock::now()};
  hyperspace.readDBFile(libName);
  const auto rend{std::chrono::steady_clock::now()};
  const std::chrono::duration<double> elapsed_seconds{rend - rstart};
  std::cout << "Time to read hyperspace : " << elapsed_seconds.count()
            << std::endl;

#if 0
  SECTION("Fragged Mol") {
    std::vector<std::unique_ptr<ROMol>> fragSet;
    fragSet.emplace_back("c1ccccc1N[1*]"_smiles);
    fragSet.emplace_back("N1CCC1C(=O)[1*]"_smiles);
    auto results = hyperspace.searchFragSet(fragSet);
    CHECK(results.size() == 1);
  }
#endif

#if 1
  SECTION("WholeMol") {
    const std::vector<std::string> smis{"c1ccccc1C(=O)N1CCCC1",
                                        "c1ccccc1NC(=O)C1CCN1",
                                        "c12ccccc1c(N)nc(N)n2",
                                        "c12ccc(C)cc1[nH]nc2C(=O)NCc1cncs1",
                                        "c1nncn1",
                                        "C(=O)NC(CC)C(=O)N(CC)C"};
    const std::vector<size_t> numRes{1000, 1000, 1000, 108, 1000, 1000};
    for (size_t i = 0; i < smis.size(); ++i) {
      //      if (i != 5) {
      //        continue;
      //      }
      std::cout << "TTTTTTTTTTTTTTTTTTTTTTTTTTT : " << smis[i] << std::endl;
      auto queryMol = v2::SmilesParse::MolFromSmarts(smis[i]);
      const auto start{std::chrono::steady_clock::now()};
      auto results = hyperspace.substructureSearch(*queryMol);
      const auto end{std::chrono::steady_clock::now()};
      const std::chrono::duration<double> elapsed_seconds{end - start};
      std::cout << "Elapsed time : " << elapsed_seconds.count() << std::endl;
      std::cout << "Number of hits : " << results.hitMolecules().size()
                << std::endl;
      std::string outFile =
          std::string("/Users/david/Projects/FreedomSpace/hyperspace_hits_") +
          std::to_string(i) + ".smi";
      std::ofstream of(outFile);
      for (const auto &r : results.hitMolecules()) {
        of << MolToSmiles(*r) << " " << r->getProp<std::string>("_Name")
           << std::endl;
      }
      CHECK(results.hitMolecules().size() == numRes[i]);
      std::cout << results.maxNumResults() << std::endl;
    }
  }
#endif
}

TEST_CASE("Small query", "[Small query]") {
  // Making sure it works when the query has fewer bonds than maxBondSplits.
  std::string libName =
      fName + "/Code/GraphMol/HyperspaceSearch/data/triazole_space.txt";
  Hyperspace hyperspace;
  hyperspace.readTextFile(libName);
  auto queryMol = "C=CC"_smiles;
  auto results = hyperspace.substructureSearch(*queryMol);
  // The number of results is immaterial, it just matters that the search
  // finished.
  CHECK(results.hitMolecules().size() == 0);
}

TEST_CASE("Random Hits", "[Random Hits]") {
  std::string libName = "/Users/david/Projects/FreedomSpace/Syntons_5567.spc";
  Hyperspace hyperspace;
  hyperspace.readDBFile(libName);

  auto queryMol = "c1ccccc1C(=O)N1CCCC1"_smiles;
  HyperspaceSearchParams params;
  params.maxBondSplits = 2;
  params.maxHits = 100;
  params.randomSample = true;
  params.randomSeed = 1;
  auto results = hyperspace.substructureSearch(*queryMol, params);
  std::map<std::string, int> libCounts;
  for (const auto &m : results.hitMolecules()) {
    std::string lib(m->getProp<std::string>("_Name").substr(0, 2));
    if (const auto &c = libCounts.find(lib); c == libCounts.end()) {
      libCounts.insert(std::make_pair(lib, 1));
    } else {
      c->second++;
    }
  }
  CHECK(results.hitMolecules().size() == 100);
  std::map<std::string, int> expCounts{{"a1", 73}, {"a6", 6}, {"a7", 21}};
  CHECK(expCounts == libCounts);
}

TEST_CASE("Later hits", "[Later Hits]") {
  // Test use of params.hitStart
  std::string libName = "/Users/david/Projects/FreedomSpace/Syntons_5567.spc";
  Hyperspace hyperspace;
  hyperspace.readDBFile(libName);

  auto queryMol = "c1ccccc1C(=O)N1CCCC1"_smiles;
  HyperspaceSearchParams params;
  params.maxBondSplits = 2;
  params.maxHits = 200;
  auto results = hyperspace.substructureSearch(*queryMol, params);
  std::vector<std::string> hitNames1;
  for (const auto &m : results.hitMolecules()) {
    hitNames1.push_back(m->getProp<std::string>("_Name"));
  }

  params.maxHits = 100;
  params.hitStart = 101;
  results = hyperspace.substructureSearch(*queryMol, params);
  std::vector<std::string> hitNames2;
  for (const auto &m : results.hitMolecules()) {
    hitNames2.push_back(m->getProp<std::string>("_Name"));
  }
  CHECK(hitNames1.size() == 200);
  CHECK(hitNames2.size() == 100);
  for (int i = 0; i < 100; ++i) {
    CHECK(hitNames1[100 + i] == hitNames2[i]);
  }

  params.hitStart = 6780;
  results = hyperspace.substructureSearch(*queryMol, params);
  CHECK(results.hitMolecules().size() == 6);

  params.hitStart = 7000;
  results = hyperspace.substructureSearch(*queryMol, params);
  CHECK(results.hitMolecules().empty());
}

TEST_CASE("Complex query", "[Complex query]") {
  std::string libName =
      "/Users/david/Projects/FreedomSpace/2023-05_Freedom_synthons.spc";
  Hyperspace hyperspace;
  hyperspace.readDBFile(libName);

  auto queryMol = v2::SmilesParse::MolFromSmarts(
      "[$(c1ccncc1),$(c1cnccc1)]C(=O)N1[C&!$(CC(=O))]CCC1");
  REQUIRE(queryMol);
  HyperspaceSearchParams params;
  params.maxBondSplits = 2;
  params.maxHits = 1000;
  auto results = hyperspace.substructureSearch(*queryMol, params);
  std::cout << "Number of hits : " << results.hitMolecules().size()
            << std::endl;
  std::cout << "Poss num hits : " << results.maxNumResults() << std::endl;
  CHECK(results.hitMolecules().size() == 1000);
}