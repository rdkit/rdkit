//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#include <algorithm>
#include <cstdio>
#include <fstream>

#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SubstructLibrary/SubstructLibrary.h>
#include <GraphMol/SynthonSpaceSearch/SearchResults.h>
#include <GraphMol/SynthonSpaceSearch/Synthon.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSet.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
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
  const std::vector<std::string> smiles{"c1ccccc1CN1CCN(CC1)C(-O)c1ncc(F)cc1",
                                        "CC(C)OCc1nnc(N2CC(C)CC2)n1C1CCCC1"};
  std::vector<std::vector<size_t>> expCounts{{1, 47, 345, 20},
                                             {1, 37, 262, 41}};
  for (size_t i = 0; i < smiles.size(); ++i) {
    auto mol = v2::SmilesParse::MolFromSmiles(smiles[i]);
    REQUIRE(mol);
    bool timedOut = false;
    auto fragments = splitMolecule(*mol, 3, 100000, nullptr, timedOut);
    CHECK(fragments.size() ==
          std::accumulate(expCounts[i].begin(), expCounts[i].end(), size_t(0)));
    // The first fragment set should just be the molecule itself.
    for (size_t j = 0; j < 4; ++j) {
      const auto numFragSets = std::accumulate(
          fragments.begin(), fragments.end(), static_cast<size_t>(0),
          [&](size_t prevRes,
              const std::vector<std::unique_ptr<ROMol>> &frags) {
            if (frags.size() == j + 1) {
              return prevRes + 1;
            }
            return prevRes;
          });
      CHECK(numFragSets == expCounts[i][j]);
    }
  }
}
#endif

TEST_CASE("Enumerate") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  // Making sure it works when the query has fewer bonds than maxBondSplits.
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/triazole_space.txt";
  SynthonSpace synthonspace;
  bool cancelled = false;
  synthonspace.readTextFile(libName, cancelled);
  auto testName = std::tmpnam(nullptr);
  BOOST_LOG(rdInfoLog) << "Enumerating to " << testName << std::endl;
  synthonspace.writeEnumeratedFile(testName);

  std::string enumLibName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/triazole_space_enum.smi";

  auto loadLibrary =
      [](const std::string inFilename) -> std::map<std::string, std::string> {
    v2::FileParsers::SmilesMolSupplierParams params;
    params.titleLine = false;
    v2::FileParsers::SmilesMolSupplier suppl(inFilename, params);
    std::map<std::string, std::string> smiles;
    while (!suppl.atEnd()) {
      auto mol = suppl.next();
      if (mol) {
        smiles.insert(
            std::make_pair(mol->getProp<std::string>(common_properties::_Name),
                           std::string(MolToSmiles(*mol))));
      }
    }
    return smiles;
  };
  auto newSmiles = loadLibrary(testName);
  auto oldSmiles = loadLibrary(enumLibName);
  REQUIRE(newSmiles.size() == oldSmiles.size());
  for (const auto &[name, smiles] : oldSmiles) {
    REQUIRE(oldSmiles.find(name) != oldSmiles.end());
    REQUIRE(newSmiles.at(name) == oldSmiles.at(name));
  }
  std::remove(testName);
}

TEST_CASE("S Amide 1") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/amide_space.txt";
  std::string enumLibName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/amide_space_enum.smi";

  auto queryMol = "c1ccccc1C(=O)N1CCCC1"_smiles;
  SynthonSpace synthonspace;
  bool cancelled = false;
  synthonspace.readTextFile(libName, cancelled);
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

TEST_CASE("S Urea 1") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/urea_space.txt";
  SynthonSpace synthonspace;
  bool cancelled = false;
  synthonspace.readTextFile(libName, cancelled);
  auto queryMol = "O=C(Nc1c(CNC=O)cc[s]1)c1nccnc1"_smiles;
  auto results = synthonspace.substructureSearch(*queryMol);
  CHECK(results.getHitMolecules().size() == 2);
}

TEST_CASE("S Simple query 1") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  SynthonSpace synthonspace;
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/idorsia_toy_space_a.spc";
  synthonspace.readDBFile(libName);
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

TEST_CASE("S Triazole") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/triazole_space.txt";
  std::string enumLibName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/triazole_space_enum.smi";

  SynthonSpace synthonspace;
  bool cancelled = false;
  synthonspace.readTextFile(libName, cancelled);

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

TEST_CASE("S Quinoline") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/doebner_miller_space.txt";
  std::string enumLibName =
      fName +
      "/Code/GraphMol/SynthonSpaceSearch/data/doebner_miller_space_enum.smi";

  SynthonSpace synthonspace;
  bool cancelled = false;
  synthonspace.readTextFile(libName, cancelled);
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

TEST_CASE("S Substructure in 1 reagent") {
  // Making sure it works when the query is a complete substructure of 1
  // of the synthons in the library, so the whole library is a hit.
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/triazole_space.txt";
  SynthonSpace synthonspace;
  bool cancelled = false;
  synthonspace.readTextFile(libName, cancelled);
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
    CHECK(results.getHitMolecules().empty());
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
    bool cancelled = false;
    synthonspace.readTextFile(libName, cancelled);
    const auto &rnames = synthonspace.getReactionNames();
    const auto rs = synthonspace.getReaction(rnames.front());
    CHECK(rs->getConnectorRegions().size() == 30);
  }
}

TEST_CASE("DB Writer") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/doebner_miller_space.txt";
  SynthonSpace synthonspace;
  bool cancelled = false;
  synthonspace.readTextFile(libName, cancelled);
  CHECK(synthonspace.getNumReactions() == 1);
  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGen(
      MorganFingerprint::getMorganGenerator<std::uint64_t>(2));
  synthonspace.buildSynthonFingerprints(*fpGen);

  auto spaceName = std::tmpnam(nullptr);

  synthonspace.writeDBFile(spaceName);

  SynthonSpace newsynthonspace;
  newsynthonspace.readDBFile(spaceName);
  std::shared_ptr<SynthonSet> irxn;
  CHECK_NOTHROW(irxn = newsynthonspace.getReaction("doebner-miller-quinoline"));
  const auto &orxn = synthonspace.getReaction("doebner-miller-quinoline");
  CHECK(irxn->getId() == orxn->getId());
  CHECK(irxn->getConnectorRegions().size() ==
        orxn->getConnectorRegions().size());
  CHECK(*irxn->getConnRegFP() == *orxn->getConnRegFP());
  CHECK(irxn->getConnectors() == orxn->getConnectors());
  CHECK(irxn->getSynthons().size() == orxn->getSynthons().size());
  CHECK(newsynthonspace.hasFingerprints());
  for (size_t i = 0; i < irxn->getSynthons().size(); ++i) {
    CHECK(irxn->getSynthons()[i].size() == orxn->getSynthons()[i].size());
    for (size_t j = 0; j < irxn->getSynthons().size(); ++j) {
      CHECK(irxn->getSynthons()[i][j].first == orxn->getSynthons()[i][j].first);
      CHECK(*irxn->getSynthons()[i][j].second->getFP() ==
            *orxn->getSynthons()[i][j].second->getFP());
    }
  }
  std::remove(spaceName);
}

TEST_CASE("DB Converter") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/doebner_miller_space.txt";

  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGen(
      MorganFingerprint::getMorganGenerator<std::uint64_t>(2));

  auto spaceName = std::tmpnam(nullptr);
  bool cancelled;
  convertTextToDBFile(libName, spaceName, cancelled, fpGen.get());

  SynthonSpace synthonspace;
  synthonspace.readTextFile(libName, cancelled);
  CHECK(synthonspace.getNumReactions() == 1);
  synthonspace.buildSynthonFingerprints(*fpGen);
  synthonspace.buildAddAndSubstractFingerprints(*fpGen);

  SynthonSpace newsynthonspace;
  newsynthonspace.readDBFile(spaceName);

  std::shared_ptr<SynthonSet> irxn;
  CHECK_NOTHROW(irxn = newsynthonspace.getReaction("doebner-miller-quinoline"));
  const auto &orxn = synthonspace.getReaction("doebner-miller-quinoline");
  CHECK(irxn->getId() == orxn->getId());
  CHECK(irxn->getConnectorRegions().size() ==
        orxn->getConnectorRegions().size());
  CHECK(*irxn->getConnRegFP() == *orxn->getConnRegFP());
  CHECK(irxn->getConnectors() == orxn->getConnectors());
  CHECK(irxn->getSynthons().size() == orxn->getSynthons().size());
  CHECK(newsynthonspace.hasFingerprints());
  for (size_t i = 0; i < irxn->getSynthons().size(); ++i) {
    CHECK(irxn->getSynthons()[i].size() == orxn->getSynthons()[i].size());
    for (size_t j = 0; j < irxn->getSynthons().size(); ++j) {
      CHECK(irxn->getSynthons()[i][j].first == orxn->getSynthons()[i][j].first);
      CHECK(*irxn->getSynthons()[i][j].second->getFP() ==
            *orxn->getSynthons()[i][j].second->getFP());
    }
  }
  std::remove(spaceName);
  // Check it behaves gracefully with a missing file
  CHECK_THROWS(synthonspace.readDBFile(spaceName));
}

TEST_CASE("S Biggy") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.csv";
  SynthonSpace synthonspace;
  bool cancelled = false;
  synthonspace.readTextFile(libName, cancelled);
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
    std::ofstream ofs("new_hits.txt");
    for (size_t i = 0; i < results.getHitMolecules().size(); ++i) {
      ofs << i << " : " << MolToSmiles(*results.getHitMolecules()[i]) << "  "
          << results.getHitMolecules()[i]->getProp<std::string>(
                 common_properties::_Name)
          << std::endl;
    }
  }
}

TEST_CASE("S Small query") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  // Making sure it works when the query has fewer bonds than maxBondSplits.
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/triazole_space.txt";
  SynthonSpace synthonspace;
  bool cancelled = false;
  synthonspace.readTextFile(libName, cancelled);
  auto queryMol = "C=CC"_smiles;
  auto results = synthonspace.substructureSearch(*queryMol);
  // The number of results is immaterial, it just matters that the search
  // finished.
  CHECK(results.getHitMolecules().empty());
}

TEST_CASE("S Random Hits") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.csv";
  SynthonSpace synthonspace;
  bool cancelled = false;
  synthonspace.readTextFile(libName, cancelled);

  auto queryMol = "c1ccccc1C(=O)N1CCCC1"_smiles;
  SynthonSpaceSearchParams params;
  params.maxBondSplits = 2;
  params.maxHits = 100;
  params.randomSample = true;
  params.randomSeed = 1;
  auto results = synthonspace.substructureSearch(*queryMol, params);
  std::map<std::string, int> libCounts;
  for (const auto &m : results.getHitMolecules()) {
    auto molName = m->getProp<std::string>(common_properties::_Name);
    std::string lib(molName.substr(molName.length() - 2));
    if (const auto &c = libCounts.find(lib); c == libCounts.end()) {
      libCounts.insert(std::make_pair(lib, 1));
    } else {
      c->second++;
    }
  }
  CHECK(results.getHitMolecules().size() == 100);
  std::map<std::string, int> expCounts{{"a1", 61}, {"a6", 10}, {"a7", 29}};
  CHECK(expCounts == libCounts);
}

TEST_CASE("S Later hits") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  // Test use of params.hitStart
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.csv";
  SynthonSpace synthonspace;
  bool cancelled = false;
  synthonspace.readTextFile(libName, cancelled);

  auto queryMol = "c1ccccc1C(=O)N1CCCC1"_smiles;
  SynthonSpaceSearchParams params;
  params.maxBondSplits = 2;
  params.maxHits = 200;
  auto results = synthonspace.substructureSearch(*queryMol, params);
  std::vector<std::string> hitNames1;
  for (const auto &m : results.getHitMolecules()) {
    hitNames1.push_back(m->getProp<std::string>(common_properties::_Name));
  }

  params.maxHits = 100;
  params.hitStart = 101;
  results = synthonspace.substructureSearch(*queryMol, params);
  std::vector<std::string> hitNames2;
  for (const auto &m : results.getHitMolecules()) {
    hitNames2.push_back(m->getProp<std::string>(common_properties::_Name));
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

TEST_CASE("S Complex query") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  // Just to demonstrate that a complex query works.
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.csv";
  SynthonSpace synthonspace;
  bool cancelled = false;
  synthonspace.readTextFile(libName, cancelled);

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

TEST_CASE("S Map numbers in connectors") {
  // Map numbers might occur in the connectors, e.g. [1*:1] as well
  // as [1*].  This checks that that is the case.
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/map_numbers.txt";
  SynthonSpace synthonspace;
  bool cancelled = false;
  synthonspace.readTextFile(libName, cancelled);

  auto queryMol = "c1ccccc1C(=O)N1CCCC1"_smarts;
  REQUIRE(queryMol);
  auto results = synthonspace.substructureSearch(*queryMol);
  // These were missing before map numbers were accommodated.
  std::set<std::string> missNames{
      "67468;30577;29389;a7",  "67468;249279;29389;a7", "67468;24773;29389;a7",
      "67468;29593;29389;a7",  "67468;308698;29389;a7", "67468;56491;29389;a7",
      "67468;265474;29389;a7", "67468;15535;29389;a7",  "67468;44908;29389;a7",
      "67468;59597;29389;a7",  "67468;45686;29389;a7"};
  std::set<std::string> hitNames;
  for (const auto &hm : results.getHitMolecules()) {
    hitNames.insert(hm->getProp<std::string>(common_properties::_Name));
  }
  CHECK(results.getHitMolecules().size() == 11);
  CHECK(hitNames == missNames);
}

TEST_CASE("Greg Space Failure") {
  // This failed at one point due to the aliphatic synthon, aromatic
  // product issue.
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/gregs_space_fail.txt";
  SynthonSpace synthonspace;
  bool cancelled = false;
  synthonspace.readTextFile(libName, cancelled);

  auto queryMol =
      "Cc1nn(C)c(C)c1-c1nc(Cn2cc(CNC(C)C(=O)NC3CCCC3)nn2)no1"_smarts;
  REQUIRE(queryMol);
  SynthonSpaceSearchParams params;
  params.maxBondSplits = 2;

  auto results = synthonspace.substructureSearch(*queryMol, params);
  CHECK(results.getHitMolecules().size() == 1);
}

TEST_CASE("DOS File") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/amide_space_dos.txt";
  SynthonSpace synthonspace;
  bool cancelled = false;
  synthonspace.readTextFile(libName, cancelled);
  CHECK(synthonspace.getNumProducts() == 12);
}

TEST_CASE("Synthon Error") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  {
    std::string libName =
        fName + "/Code/GraphMol/SynthonSpaceSearch/data/amide_space_error.txt";
    SynthonSpace synthonspace;
    bool cancelled = false;
    CHECK_THROWS(synthonspace.readTextFile(libName, cancelled));
  }
  {
    std::string libName =
        fName + "/Code/GraphMol/SynthonSpaceSearch/data/synthon_error.txt";
    SynthonSpace synthonspace;
    bool cancelled = false;
    CHECK_THROWS(synthonspace.readTextFile(libName, cancelled));
  }
}

TEST_CASE("S Freedom Space") {
  std::string libName =
      "/Users/david/Projects/SynthonSpaceTests/FreedomSpace/2024-09_Freedom_synthons_rdkit_new.spc";
  libName =
      "/Users/david/Projects/SynthonSpaceTests/REAL/2024-09_RID-4-Cozchemix/2024-09_REAL_synthons_rdkit_3000.spc";
  SynthonSpace synthonspace;
  synthonspace.readDBFile(libName);

  auto m = "c12ccc(C)cc1[nH]nc2C(=O)NCc1cncs1"_smarts;
  SynthonSpaceSearchParams params;
  params.maxHits = 1000;
  SearchResults results;
  results = synthonspace.substructureSearch(*m, params);
  std::cout << "Number of results : " << results.getHitMolecules().size()
            << std::endl;
  int i = 0;
  int numHits = results.getHitMolecules().size();
  for (const auto &mol : results.getHitMolecules()) {
    if (i < 10 || i > numHits - 10) {
      std::cout << i << " : " << MolToSmiles(*mol) << " : "
                << mol->getProp<std::string>(common_properties::_Name)
                << std::endl;
    }
    ++i;
  }
}