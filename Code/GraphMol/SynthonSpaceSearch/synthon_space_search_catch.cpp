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
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SynthonSpaceSearch/Synthon.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <catch2/catch_all.hpp>

using namespace RDKit;
using namespace RDKit::SynthonSpaceSearch;

using namespace RDKit::SynthonSpaceSearch::details;

const char *rdbase = getenv("RDBASE");

std::unique_ptr<SubstructLibrary> loadSubstructLibrary(
    const std::string &smiFile) {
  std::unique_ptr<SubstructLibrary> subsLib(new SubstructLibrary());
  v2::FileParsers::SmilesMolSupplierParams params;
  params.titleLine = false;
  v2::FileParsers::SmilesMolSupplier suppl(smiFile, params);
  while (!suppl.atEnd()) {
    subsLib->addMol(*suppl.next());
  }
  return subsLib;
}

#if 1
TEST_CASE("Test splits 1") {
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
#endif

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
  auto results = synthonspace.substructureSearch(*queryMol, params);
  CHECK(results.getHitMolecules().size() == 2);
  std::set<std::string> resSmi;
  for (const auto &r : results.getHitMolecules()) {
    resSmi.insert(MolToSmiles(*r));
  }

  auto subsLib = loadSubstructLibrary(enumLibName);
  auto query = "c1ccccc1C(=O)N1CCCC1"_smarts;
  auto enumRes = subsLib->getMatches(*query);
  std::set<std::string> enumSmi;
  for (auto i : enumRes) {
    enumSmi.insert(MolToSmiles(*subsLib->getMol(i)));
  }
  CHECK(resSmi == enumSmi);
}

TEST_CASE("Urea 1") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/urea_space.txt";
  SynthonSpace synthonspace;
  synthonspace.readTextFile(libName);
  SECTION("Single fragged molecule") {
    std::vector<std::unique_ptr<ROMol>> fragSet;
    fragSet.emplace_back("O=C(NC1COC1)[1*]"_smiles);
    fragSet.emplace_back("O=C(Nc1c(CN[1*])cc[s]1)[2*]"_smiles);
    fragSet.emplace_back("Fc1nccnc1[2*]"_smiles);
    auto results = synthonspace.searchFragSet(fragSet);
    CHECK(results.size() == 1);
  }
  SECTION("Single molecule with fragging") {
    auto queryMol = "O=C(Nc1c(CNC=O)cc[s]1)c1nccnc1"_smiles;
    auto results = synthonspace.substructureSearch(*queryMol);
    CHECK(results.getHitMolecules().size() == 2);
  }
}

TEST_CASE("Simple query 1") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/urea_3.txt";
  SynthonSpace synthonspace;
  synthonspace.readTextFile(libName);
  SECTION("Single fragged molecule") {
    std::vector<std::unique_ptr<ROMol>> fragSet;
    fragSet.emplace_back("c1ccccc1[1*]"_smiles);
    fragSet.emplace_back("C1CCCN1C(=O)[1*]"_smiles);
    auto results = synthonspace.searchFragSet(fragSet);
    CHECK(results.size() == 1);
    CHECK(results.front().numHits == 220);
  }
  SECTION("Binary File") {
    SynthonSpace synthonspace;
    std::string libName =
        fName + "/Code/GraphMol/SynthonSpaceSearch/data/idorsia_toy_space.spc";
    synthonspace.readDBFile(libName);
    // should give 220 hits for urea-3
    auto queryMol = "c1ccccc1C(=O)N1CCCC1"_smiles;
    auto results = synthonspace.substructureSearch(*queryMol);
    CHECK(results.getHitMolecules().size() == 220);
    CHECK(results.getMaxNumResults() == 220);
  }
  SECTION("Single molecule with fragging") {
    {
      // should give 220 hits for urea-3
      auto queryMol = "c1ccccc1C(=O)N1CCCC1"_smiles;
      auto results = synthonspace.substructureSearch(*queryMol);
      CHECK(results.getHitMolecules().size() == 220);
      CHECK(results.getMaxNumResults() == 220);
    }
    {
      auto queryMol = "O=C(Nc1c(CNC=O)cc[s]1)c1nccnc1"_smiles;
      auto results = synthonspace.substructureSearch(*queryMol);
      CHECK(results.getHitMolecules().size() == 20);
    }
  }
}

TEST_CASE("Triazole") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/triazole_space.txt";
  std::string enumLibName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/triazole_space_enum.smi";

  SynthonSpace synthonspace;
  synthonspace.readTextFile(libName);

  SECTION("Fragged Mol") {
    auto queryMol =
        "OCC([1*])=NN=[2*].C1CCCC1N([3*])[1*].CC1CCN(C1)C(=[2*])[3*]"_smiles;
    REQUIRE(queryMol);
    std::vector<std::unique_ptr<ROMol>> queryFrags;
    MolOps::getMolFrags(*queryMol, queryFrags, false);
    auto results = synthonspace.searchFragSet(queryFrags);
    CHECK(results.size() == 1);
    CHECK(results.front().numHits == 4);
  }
  SECTION("Full Molecule") {
    auto queryMol = "OCc1ncnn1"_smarts;
    REQUIRE(queryMol);
    auto results = synthonspace.substructureSearch(*queryMol);
    CHECK(results.getHitMolecules().size() == 8);
    std::set<std::string> resSmi;
    for (const auto &r : results.getHitMolecules()) {
      resSmi.insert(MolToSmiles(*r));
    }

    auto subsLib = loadSubstructLibrary(enumLibName);
    auto enumRes = subsLib->getMatches(*queryMol);
    std::set<std::string> enumSmi;
    for (auto i : enumRes) {
      enumSmi.insert(MolToSmiles(*subsLib->getMol(i)));
    }
    CHECK(resSmi == enumSmi);
  }
}

TEST_CASE("Quinoline") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/doebner_miller_space.txt";
  std::string enumLibName =
      fName +
      "/Code/GraphMol/SynthonSpaceSearch/data/doebner_miller_space_enum.smi";

  SynthonSpace synthonspace;
  synthonspace.readTextFile(libName);
  {
    auto queryMol = "c1ccccn1"_smiles;
    auto results = synthonspace.substructureSearch(*queryMol);
    CHECK(results.getHitMolecules().size() == 12);
    std::set<std::string> resSmi;
    for (const auto &r : results.getHitMolecules()) {
      resSmi.insert(MolToSmiles(*r));
    }
    auto subsLib = loadSubstructLibrary(enumLibName);
    auto enumRes = subsLib->getMatches(*queryMol);
    std::set<std::string> enumSmi;
    for (auto i : enumRes) {
      enumSmi.insert(MolToSmiles(*subsLib->getMol(i)));
    }
    CHECK(resSmi == enumSmi);
  }
}

TEST_CASE("Substructure in 1 reagent") {
  // Making sure it works when the query is a complete substructure of 1
  // of the synthons in the library, so the whole library is a hit.
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/triazole_space.txt";
  SynthonSpace synthonspace;
  synthonspace.readTextFile(libName);
  {
    auto queryMol = "N1CCCC1"_smiles;
    auto results = synthonspace.substructureSearch(*queryMol);
    CHECK(results.getHitMolecules().size() == 8);
  }
  {
    auto queryMol = "N1CCC(C(F)(F)F)C1"_smiles;
    auto results = synthonspace.substructureSearch(*queryMol);
    CHECK(results.getHitMolecules().size() == 4);
  }
  {
    auto queryMol = "C1CCCCC1"_smiles;
    auto results = synthonspace.substructureSearch(*queryMol);
    CHECK(results.getHitMolecules().size() == 0);
  }
}

TEST_CASE("Connector Regions") {
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
    REQUIRE(rdbase);
    std::string fName(rdbase);
    std::string libName =
        fName + "/Code/GraphMol/SynthonSpaceSearch/data/urea_3.txt";
    SynthonSpace synthonspace;
    synthonspace.readTextFile(libName);
    const auto &rs = synthonspace.getReactions().begin();
    CHECK(rs->second->getConnectorRegions().size() == 30);
  }
}

TEST_CASE("DB Writer") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/doebner_miller_space.txt";
  SynthonSpace synthonspace;
  synthonspace.readTextFile(libName);
  CHECK(synthonspace.getNumReactions() == 1);
  synthonspace.writeDBFile("doebner_miller_space.spc");

  SynthonSpace newsynthonspace;
  newsynthonspace.readDBFile("doebner_miller_space.spc");
  auto it = newsynthonspace.getReactions().find("doebner-miller-quinoline");
  CHECK(it != newsynthonspace.getReactions().end());
  const auto &irxn = it->second;
  const auto &orxn =
      synthonspace.getReactions().find("doebner-miller-quinoline")->second;
  CHECK(irxn->getId() == orxn->getId());
  CHECK(irxn->getConnectorRegions().size() ==
        orxn->getConnectorRegions().size());
  CHECK(*irxn->getConnRegFP() == *orxn->getConnRegFP());
  CHECK(irxn->getConnectors() == orxn->getConnectors());
  CHECK(irxn->getSynthons().size() == orxn->getSynthons().size());
  for (size_t i = 0; i < irxn->getSynthons().size(); ++i) {
    CHECK(irxn->getSynthons()[i].size() == orxn->getSynthons()[i].size());
    for (size_t j = 0; j < irxn->getSynthons().size(); ++j) {
      CHECK(irxn->getSynthons()[i][j]->getId() ==
            orxn->getSynthons()[i][j]->getId());
    }
  }
}

TEST_CASE("Biggy") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.csv";
  SynthonSpace synthonspace;
  synthonspace.readTextFile(libName);

  SECTION("Fragged Mol") {
    std::vector<std::unique_ptr<ROMol>> fragSet;
    fragSet.emplace_back("c1ccccc1N[1*]"_smiles);
    fragSet.emplace_back("N1CCC1C(=O)[1*]"_smiles);
    auto results = synthonspace.searchFragSet(fragSet);
    CHECK(results.size() == 2);
  }

  SECTION("WholeMol") {
    const std::vector<std::string> smis{"c1ccccc1C(=O)N1CCCC1",
                                        "c1ccccc1NC(=O)C1CCN1",
                                        "c12ccccc1c(N)nc(N)n2",
                                        "c12ccc(C)cc1[nH]nc2C(=O)NCc1cncs1",
                                        "c1nncn1",
                                        "C(=O)NC(CC)C(=O)N(CC)C"};
    const std::vector<size_t> numRes{6785, 4544, 48892, 1, 29147, 5651};
    const std::vector<size_t> maxRes{6785, 4544, 48893, 1, 29312, 5869};
    SynthonSpaceSearchParams params;
    params.maxHits = -1;
    for (size_t i = 0; i < smis.size(); ++i) {
      auto queryMol = v2::SmilesParse::MolFromSmarts(smis[i]);
      auto results = synthonspace.substructureSearch(*queryMol, params);
      CHECK(results.getHitMolecules().size() == numRes[i]);
      CHECK(results.getMaxNumResults() == maxRes[i]);
    }
  }
}

#if 0
// The whole FreedomSpace synthon library is maybe a bit big for the
// repo, even in text format.
TEST_CASE("FreedomSpace", "[FreedomSpace]") {
  std::string libName =
      "/Users/david/Projects/FreedomSpace/2023-05_Freedom_synthons.spc";
  SynthonSpace synthonspace;
  const auto rstart{std::chrono::steady_clock::now()};
  synthonspace.readDBFile(libName);
  const auto rend{std::chrono::steady_clock::now()};
  const std::chrono::duration<double> elapsed_seconds{rend - rstart};
  std::cout << "Time to read synthonspace : " << elapsed_seconds.count()
            << std::endl;

  SECTION("Fragged Mol") {
    std::vector<std::unique_ptr<ROMol>> fragSet;
    fragSet.emplace_back("c1ccccc1N[1*]"_smiles);
    fragSet.emplace_back("N1CCC1C(=O)[1*]"_smiles);
    auto results = synthonspace.searchFragSet(fragSet);
    CHECK(results.size() == 1);
  }
  SECTION("WholeMol") {
    const std::vector<std::string> smis{"c1ccccc1C(=O)N1CCCC1",
                                        "c1ccccc1NC(=O)C1CCN1",
                                        "c12ccccc1c(N)nc(N)n2",
                                        "c12ccc(C)cc1[nH]nc2C(=O)NCc1cncs1",
                                        "c1nncn1",
                                        "C(=O)NC(CC)C(=O)N(CC)C"};
    const std::vector<size_t> numRes{1000, 1000, 1000, 108, 1000, 1000};
    for (size_t i = 0; i < smis.size(); ++i) {
      auto queryMol = v2::SmilesParse::MolFromSmarts(smis[i]);
      const auto start{std::chrono::steady_clock::now()};
      auto results = synthonspace.substructureSearch(*queryMol);
      const auto end{std::chrono::steady_clock::now()};
      const std::chrono::duration<double> elapsed_seconds{end - start};
      std::cout << "Elapsed time : " << elapsed_seconds.count() << std::endl;
      std::cout << "Number of hits : " << results.getHitMolecules().size()
                << std::endl;
      std::string outFile =
          std::string("/Users/david/Projects/FreedomSpace/synthonspace_hits_") +
          std::to_string(i) + ".smi";
      std::ofstream of(outFile);
      for (const auto &r : results.getHitMolecules()) {
        of << MolToSmiles(*r) << " " << r->getProp<std::string>("_Name")
           << std::endl;
      }
      CHECK(results.getHitMolecules().size() == numRes[i]);
      std::cout << results.getMaxNumResults() << std::endl;
    }
  }
}
#endif

TEST_CASE("Small query") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  // Making sure it works when the query has fewer bonds than maxBondSplits.
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/triazole_space.txt";
  SynthonSpace synthonspace;
  synthonspace.readTextFile(libName);
  auto queryMol = "C=CC"_smiles;
  auto results = synthonspace.substructureSearch(*queryMol);
  // The number of results is immaterial, it just matters that the search
  // finished.
  CHECK(results.getHitMolecules().size() == 0);
}

TEST_CASE("Random Hits") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.csv";
  SynthonSpace synthonspace;
  synthonspace.readTextFile(libName);

  auto queryMol = "c1ccccc1C(=O)N1CCCC1"_smiles;
  SynthonSpaceSearchParams params;
  params.maxBondSplits = 2;
  params.maxHits = 100;
  params.randomSample = true;
  params.randomSeed = 1;
  auto results = synthonspace.substructureSearch(*queryMol, params);
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
  std::map<std::string, int> expCounts{{"a1", 73}, {"a6", 6}, {"a7", 21}};
  CHECK(expCounts == libCounts);
}

TEST_CASE("Later hits") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  // Test use of params.hitStart
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.csv";
  SynthonSpace synthonspace;
  synthonspace.readTextFile(libName);

  auto queryMol = "c1ccccc1C(=O)N1CCCC1"_smiles;
  SynthonSpaceSearchParams params;
  params.maxBondSplits = 2;
  params.maxHits = 200;
  auto results = synthonspace.substructureSearch(*queryMol, params);
  std::vector<std::string> hitNames1;
  for (const auto &m : results.getHitMolecules()) {
    hitNames1.push_back(m->getProp<std::string>("_Name"));
  }

  params.maxHits = 100;
  params.hitStart = 101;
  results = synthonspace.substructureSearch(*queryMol, params);
  std::vector<std::string> hitNames2;
  for (const auto &m : results.getHitMolecules()) {
    hitNames2.push_back(m->getProp<std::string>("_Name"));
  }
  CHECK(hitNames1.size() == 200);
  CHECK(hitNames2.size() == 100);
  for (int i = 0; i < 100; ++i) {
    CHECK(hitNames1[100 + i] == hitNames2[i]);
  }

  params.hitStart = 6780;
  results = synthonspace.substructureSearch(*queryMol, params);
  CHECK(results.getHitMolecules().size() == 6);

  params.hitStart = 7000;
  results = synthonspace.substructureSearch(*queryMol, params);
  CHECK(results.getHitMolecules().empty());
}

TEST_CASE("Complex query") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  // Just to demonstrate that a complex query works.
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.csv";
  SynthonSpace synthonspace;
  synthonspace.readTextFile(libName);

  auto queryMol = v2::SmilesParse::MolFromSmarts(
      "[$(c1ccccc1),$(c1ccncc1),$(c1cnccc1)]C(=O)N1[C&!$(CC(=O))]CCC1");
  REQUIRE(queryMol);
  SynthonSpaceSearchParams params;
  params.maxBondSplits = 2;
  params.maxHits = 1000;
  auto results = synthonspace.substructureSearch(*queryMol, params);
  CHECK(results.getHitMolecules().size() == 1000);
  CHECK(results.getMaxNumResults() == 3257);
}

TEST_CASE("Map numbers in connectors") {
  // Map numbers might occur in the connectors, e.g. [1*:1] as well
  // as [1*].  This checks that that is the case.
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/map_numbers.txt";
  SynthonSpace synthonspace;
  synthonspace.readTextFile(libName);

  SECTION("Fragged Mol") {
    std::vector<std::unique_ptr<ROMol>> fragSet;
    fragSet.emplace_back("c1ccccc1C(=O)[1*]"_smiles);
    fragSet.emplace_back("[1*]N1CCCC1"_smiles);
    auto results = synthonspace.searchFragSet(fragSet);
    CHECK(results.size() == 1);
  }

  SECTION("Whole Mol") {
    auto queryMol = "c1ccccc1C(=O)N1CCCC1"_smarts;
    REQUIRE(queryMol);
    auto results = synthonspace.substructureSearch(*queryMol);
    // These were missing before map numbers were accommodated.
    std::set<std::string> missNames{
        "a7_67468_30577_29389",  "a7_67468_249279_29389",
        "a7_67468_24773_29389",  "a7_67468_29593_29389",
        "a7_67468_308698_29389", "a7_67468_56491_29389",
        "a7_67468_265474_29389", "a7_67468_15535_29389",
        "a7_67468_44908_29389",  "a7_67468_59597_29389",
        "a7_67468_45686_29389"};
    std::set<std::string> hitNames;
    for (const auto &hm : results.getHitMolecules()) {
      hitNames.insert(hm->getProp<std::string>("_Name"));
    }
    CHECK(results.getHitMolecules().size() == 11);
    CHECK(hitNames == missNames);
  }
}
