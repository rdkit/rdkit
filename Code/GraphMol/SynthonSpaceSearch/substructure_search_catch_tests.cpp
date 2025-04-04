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
#include <cstdio>
#include <filesystem>
#include <fstream>

#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/GeneralizedSubstruct/XQMol.h>
#include <GraphMol/GenericGroups/GenericGroups.h>
#include <GraphMol/SubstructLibrary/SubstructLibrary.h>
#include <GraphMol/SynthonSpaceSearch/SearchResults.h>
#include <GraphMol/SynthonSpaceSearch/Synthon.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSet.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
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

std::map<std::string, std::string> loadLibrary(const std::string inFilename) {
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

TEST_CASE("Test splits 1") {
  const std::vector<std::string> smiles{"c1ccccc1CN1CCN(CC1)C(-O)c1ncc(F)cc1",
                                        "CC(C)OCc1nnc(N2CC(C)CC2)n1C1CCCC1",
                                        "c1ccccc1Oc1cccc2[nH]ccc12"};
  std::vector<std::vector<size_t>> expCounts{
      {1, 47, 1020, 0}, {1, 37, 562, 0}, {1, 29, 397, 0}};
  for (size_t i = 0; i < smiles.size(); ++i) {
    auto mol = v2::SmilesParse::MolFromSmiles(smiles[i]);
    REQUIRE(mol);
    bool timedOut = false;
    auto fragments = splitMolecule(*mol, 3, 100000, nullptr, 1, timedOut);
    CHECK(fragments.size() ==
          std::accumulate(expCounts[i].begin(), expCounts[i].end(), size_t(0)));
    // The first fragment set should just be the molecule itself.  There
    // shouldn't be any 4 fragment sets, but check.
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
  SubstructMatchParameters matchParams;
  SynthonSpaceSearchParams params;
  auto results =
      synthonspace.substructureSearch(*queryMol, matchParams, params);
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
  {
    // Test for multiple threads.
    auto queryMol = "c1ccccc1C(=O)N1CCCC1"_smiles;
    SynthonSpaceSearchParams params;
    params.numThreads = -1;
    SubstructMatchParameters matchParams;
    auto results =
        synthonspace.substructureSearch(*queryMol, matchParams, params);
    CHECK(results.getHitMolecules().size() == 220);
    CHECK(results.getMaxNumResults() == 220);
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
    CHECK(MolToSmiles(*buildConnRegion(*m1)) == "[1*]CN(C)C[1*]");

    auto m2 = "[1*]CN(C[2*])Cc1ccc(CN(C[3*])C[1*])cc1"_smiles;
    REQUIRE(m2);
    CHECK(MolToSmiles(*buildConnRegion(*m2)) ==
          "[1*]CN(C)C[1*].[1*]CN(C)C[1*]");

    auto m3 = "[2*]C"_smiles;
    REQUIRE(m3);
    CHECK(MolToSmiles(*buildConnRegion(*m3)) == "[1*]C");

    auto m4 = "[1*]c1cnccc1"_smiles;
    REQUIRE(m4);
    CHECK(MolToSmiles(*buildConnRegion(*m4)) == "[1*]c(cc)cn");
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
    CHECK(rs->getConnectorRegions().size() == 32);
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
  CHECK(irxn->getConnRegFPs().size() == orxn->getConnRegFPs().size());
  for (size_t i = 0; i < irxn->getConnRegFPs().size(); ++i) {
    CHECK(*irxn->getConnRegFPs()[i] == *orxn->getConnRegFPs()[i]);
  }
  CHECK(irxn->getConnectors() == orxn->getConnectors());
  CHECK(irxn->getSynthons().size() == orxn->getSynthons().size());
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

TEST_CASE("S Small query") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  // Making sure it works when the query has fewer bonds than the maximum
  // number of synthons.
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

  SubstructMatchParameters matchParams;
  SynthonSpaceSearchParams params;
  auto results =
      synthonspace.substructureSearch(*queryMol, matchParams, params);
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

TEST_CASE("Amino Acid") {
  // The issue here was that the SMARTS pattern should match just one synthon
  // in the "library" but doesn't because the connector is on the nitrogen
  // of the amino acid which says !$(N-[!#6;!#1]) i.e. the nitrogen can
  // only be attached to a carbon or hydrogen, and in the synthon it's
  // attached to a dummy atom.
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/amino_acid.txt";
  SynthonSpace synthonspace;
  bool cancelled = false;
  synthonspace.readTextFile(libName, cancelled);

  auto queryMol =
      "[$(C-[C;!$(C=[!#6])]-[N;!H0;!$(N-[!#6;!#1]);!$(N-C=[O,N,S])])](=O)([O;H,-])"_smarts;
  REQUIRE(queryMol);
  SubstructMatchParameters matchParams;
  SynthonSpaceSearchParams params;
  auto results =
      synthonspace.substructureSearch(*queryMol, matchParams, params);
  CHECK(results.getHitMolecules().size() == 1);
}

TEST_CASE("Extended Query") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/extended_query.csv";
  SynthonSpace synthonspace;
  bool cancelled = false;
  synthonspace.readTextFile(libName, cancelled);

  {
    auto queryMol =
        v2::SmilesParse::MolFromSmarts("[#6]-*.c1nc2cccnc2n1 |m:1:3.10|");
    REQUIRE(queryMol);
    auto xrq = GeneralizedSubstruct::createExtendedQueryMol(*queryMol);
    auto results = synthonspace.substructureSearch(xrq);
    CHECK(results.getHitMolecules().size() == 12);

    MolOps::AdjustQueryParameters aqps;
    aqps.adjustHeavyDegree = true;
    aqps.adjustHeavyDegreeFlags =
        MolOps::AdjustQueryWhichFlags::ADJUST_IGNORECHAINS;
    auto xrq1 = GeneralizedSubstruct::createExtendedQueryMol(*queryMol, true,
                                                             true, true, aqps);
    auto results1 = synthonspace.substructureSearch(xrq1);
    CHECK(results1.getHitMolecules().size() == 5);
  }

  {
    auto queryMol = R"CTAB(
  Mrv2401 02062512582D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 12 13 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -2.4167 7.8733 0 0
M  V30 2 C -3.7503 7.1033 0 0
M  V30 3 C -3.7503 5.5632 0 0
M  V30 4 N -2.4167 4.7932 0 0
M  V30 5 C -1.083 5.5632 0 0
M  V30 6 C -1.083 7.1033 0 0
M  V30 7 C 0.3973 7.5278 0 0
M  V30 8 N 0.3104 5.0376 0 0
M  V30 9 C 1.2585 6.251 0 0
M  V30 10 C 2.7975 6.1973 0 0
M  V30 11 N 3.6136 7.5032 0 0
M  V30 12 * -2.4167 9.4133 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 2 1 6
M  V30 7 1 7 6
M  V30 8 1 5 8
M  V30 9 1 8 9
M  V30 10 2 7 9
M  V30 11 1 9 10
M  V30 12 1 10 11
M  V30 13 1 1 12
M  V30 END BOND
M  V30 LINKNODE 1 3 2 10 9 10 11
M  V30 BEGIN SGROUP
M  V30 1 SUP 0 ATOMS=(1 12) SAP=(3 12 1 1) XBONDS=(1 13) LABEL=ARY ESTATE=E
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(queryMol);
    GenericGroups::setGenericQueriesFromProperties(*queryMol);
    auto xrq = GeneralizedSubstruct::createExtendedQueryMol(*queryMol);
    auto results = synthonspace.substructureSearch(xrq);
    CHECK(results.getHitMolecules().size() == 2);
  }

  {
    auto queryMol = R"CTAB(qry
  Mrv2305 09052314502D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 13 13 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N -4.75 1.9567 0 0
M  V30 2 C -6.0837 1.1867 0 0
M  V30 3 C -6.0837 -0.3534 0 0
M  V30 4 C -4.75 -1.1234 0 0
M  V30 5 C -3.4163 -0.3534 0 0
M  V30 6 C -3.4163 1.1867 0 0
M  V30 7 N -1.9692 1.7134 0 0
M  V30 8 N -1.8822 -0.7768 0 0
M  V30 9 C -1.0211 0.4999 0 0
M  V30 10 C 0.5179 0.5536 0 0
M  V30 11 N 1.2409 1.9133 0 0
M  V30 12 * -5.6391 -0.0967 0 0
M  V30 13 C -5.6391 -2.4067 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 2 1 6
M  V30 7 1 8 9
M  V30 8 1 7 6
M  V30 9 1 5 8
M  V30 10 2 7 9
M  V30 11 1 9 10
M  V30 12 1 10 11
M  V30 13 1 12 13 ENDPTS=(3 4 3 2) ATTACH=ANY
M  V30 END BOND
M  V30 LINKNODE 1 2 2 10 9 10 11
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(queryMol);
    GenericGroups::setGenericQueriesFromProperties(*queryMol);
    auto xrq = GeneralizedSubstruct::createExtendedQueryMol(*queryMol);
    auto results = synthonspace.substructureSearch(xrq);
    CHECK(results.getHitMolecules().size() == 12);
  }

  {
    // Check maxHits is working correctly.
    auto queryMol = v2::SmilesParse::MolFromSmarts(
        "[#6]-1-[#6]-c2ccccc2-[#7]-1 |LN:1:1.2|");
    REQUIRE(queryMol);
    auto xrq = GeneralizedSubstruct::createExtendedQueryMol(*queryMol);
    auto results = synthonspace.substructureSearch(xrq);
    CHECK(results.getHitMolecules().size() == 8);

    SynthonSpaceSearch::SynthonSpaceSearchParams params;
    params.maxHits = 5;
    SubstructMatchParameters mparams;
    auto results1 = synthonspace.substructureSearch(xrq, mparams, params);
    CHECK(results1.getHitMolecules().size() == 5);
  }

  {
    // Generic query check.
    auto queryMol = R"CTAB(
  Mrv2401 02062512582D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 12 13 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -2.4167 7.8733 0 0
M  V30 2 C -3.7503 7.1033 0 0
M  V30 3 C -3.7503 5.5632 0 0
M  V30 4 N -2.4167 4.7932 0 0
M  V30 5 C -1.083 5.5632 0 0
M  V30 6 C -1.083 7.1033 0 0
M  V30 7 N 0.3973 7.5278 0 0
M  V30 8 N 0.3104 5.0376 0 0
M  V30 9 C 1.2585 6.251 0 0
M  V30 10 C 2.7975 6.1973 0 0
M  V30 11 N 3.6136 7.5032 0 0
M  V30 12 * -2.4167 9.4133 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 2 1 6
M  V30 7 1 7 6
M  V30 8 1 5 8
M  V30 9 1 8 9
M  V30 10 2 7 9
M  V30 11 1 9 10
M  V30 12 1 10 11
M  V30 13 1 1 12
M  V30 END BOND
M  V30 LINKNODE 1 3 2 10 9 10 11
M  V30 BEGIN SGROUP
M  V30 1 SUP 0 ATOMS=(1 12) SAP=(3 12 1 1) XBONDS=(1 13) LABEL=ARY ESTATE=E
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(queryMol);
    GenericGroups::setGenericQueriesFromProperties(*queryMol);
    auto xrq = GeneralizedSubstruct::createExtendedQueryMol(*queryMol);
    SubstructMatchParameters mparams;
    mparams.useGenericMatchers = true;
    auto results1 = synthonspace.substructureSearch(xrq, mparams);
    CHECK(results1.getHitMolecules().size() == 2);
  }
}