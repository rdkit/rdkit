//
//  Copyright (C) 2022 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/test.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Fingerprints/AtomPairGenerator.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/Fingerprints/RDKitFPGenerator.h>
#include <GraphMol/Fingerprints/TopologicalTorsionGenerator.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/FileParsers.h>

using namespace RDKit;

TEST_CASE("includeRedundantEnvironments") {
  auto mol = "CC(=O)O"_smiles;
  REQUIRE(mol);
  SECTION("basics") {
    std::unique_ptr<FingerprintGenerator<std::uint32_t>> fpgen{
        MorganFingerprint::getMorganGenerator<std::uint32_t>(2)};
    REQUIRE(fpgen);
    {
      std::unique_ptr<SparseIntVect<std::uint32_t>> fp{
          fpgen->getCountFingerprint(*mol)};
      REQUIRE(fp);
      CHECK(fp->getTotalVal() == 8);
    }
    // turn on inclusion of redundant bits
    dynamic_cast<MorganFingerprint::MorganArguments *>(fpgen->getOptions())
        ->df_includeRedundantEnvironments = true;
    {
      std::unique_ptr<SparseIntVect<std::uint32_t>> fp{
          fpgen->getCountFingerprint(*mol)};
      REQUIRE(fp);
      CHECK(fp->getTotalVal() == 12);
    }
  }
}

TEST_CASE(
    "github #5838: FP count simulation is not accounted for when additional output is requested") {
  SECTION("as reported") {
    auto m = "OCCCCCCN"_smiles;
    REQUIRE(m);

    std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpgen(
        TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>());
    REQUIRE(fpgen);
    CHECK(fpgen->getOptions()->df_countSimulation);
    {
      AdditionalOutput ao;
      ao.allocateBitPaths();
      ao.allocateAtomToBits();
      ao.allocateAtomCounts();
      FingerprintFuncArguments args;
      args.additionalOutput = &ao;
      std::unique_ptr<ExplicitBitVect> fp{fpgen->getFingerprint(*m, args)};
      REQUIRE(fp);
      CHECK(fp->getNumOnBits() == ao.bitPaths->size());
      std::vector<int> obl;
      fp->getOnBits(obl);
      for (const auto bid : obl) {
        INFO(bid);
        CHECK(ao.bitPaths->find(bid) != ao.bitPaths->end());
      }
      std::vector<unsigned int> atomCounts{1, 2, 3, 4, 4, 3, 2, 1};
      CHECK(*ao.atomCounts == atomCounts);
      std::vector<std::vector<std::uint64_t>> atomToBits = {{0},
                                                            {0, 284, 285},
                                                            {0, 284, 285},
                                                            {0, 284, 285},
                                                            {284, 285, 384},
                                                            {284, 285, 384},
                                                            {284, 285, 384},
                                                            {384}};
      CHECK(*ao.atomToBits == atomToBits);
    }
    {
      AdditionalOutput ao;
      ao.allocateBitPaths();
      ao.allocateAtomToBits();
      ao.allocateAtomCounts();
      FingerprintFuncArguments args;
      args.additionalOutput = &ao;
      std::unique_ptr<SparseBitVect> fp{fpgen->getSparseFingerprint(*m, args)};
      REQUIRE(fp);
      CHECK(fp->getNumOnBits() == ao.bitPaths->size());
      std::vector<int> obl;
      fp->getOnBits(obl);
      for (const auto bid : obl) {
        INFO(bid);
        CHECK(ao.bitPaths->find(bid) != ao.bitPaths->end());
      }
      std::vector<unsigned int> atomCounts{1, 2, 3, 4, 4, 3, 2, 1};
      CHECK(*ao.atomCounts == atomCounts);
      std::vector<std::vector<std::uint64_t>> atomToBits = {
          {1046732804},
          {1046732804, 1046733088, 1046733089},
          {1046732804, 1046733088, 1046733089},
          {1046732804, 1046733088, 1046733089},
          {1046733088, 1046733089, 1046733188},
          {1046733088, 1046733089, 1046733188},
          {1046733088, 1046733089, 1046733188},
          {1046733188}};
      CHECK(*ao.atomToBits == atomToBits);
    }
  }
  SECTION("morgan") {
    auto m = "OCCCCCCN"_smiles;
    REQUIRE(m);

    std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpgen(
        MorganFingerprint::getMorganGenerator<std::uint64_t>(2));
    REQUIRE(fpgen);
    fpgen->getOptions()->df_countSimulation = true;
    {
      AdditionalOutput ao;
      ao.allocateBitInfoMap();
      ao.allocateAtomToBits();
      ao.allocateAtomCounts();
      FingerprintFuncArguments args;
      args.additionalOutput = &ao;
      std::unique_ptr<ExplicitBitVect> fp{fpgen->getFingerprint(*m, args)};
      REQUIRE(fp);
      CHECK(fp->getNumOnBits() == ao.bitInfoMap->size());
      std::vector<int> obl;
      fp->getOnBits(obl);
      for (const auto bid : obl) {
        INFO(bid);
        CHECK(ao.bitInfoMap->find(bid) != ao.bitInfoMap->end());
      }
      std::vector<unsigned int> atomCounts{2, 3, 3, 3, 3, 3, 3, 2};
      CHECK(*ao.atomCounts == atomCounts);
      std::vector<std::vector<std::uint64_t>> atomToBits = {
          {888, 1180},
          {320, 321, 322, 424, 1892},
          {116, 320, 321, 322, 1500, 1501, 1502},
          {320, 321, 322, 476, 477, 1500, 1501, 1502},
          {320, 321, 322, 476, 477, 1500, 1501, 1502},
          {232, 320, 321, 322, 1500, 1501, 1502},
          {320, 321, 322, 1216, 1972},
          {588, 1876},
      };
      CHECK(*ao.atomToBits == atomToBits);
    }
    {
      AdditionalOutput ao;
      ao.allocateBitInfoMap();
      ao.allocateAtomToBits();
      ao.allocateAtomCounts();
      FingerprintFuncArguments args;
      args.additionalOutput = &ao;
      std::unique_ptr<SparseBitVect> fp{fpgen->getSparseFingerprint(*m, args)};
      REQUIRE(fp);
      CHECK(fp->getNumOnBits() == ao.bitInfoMap->size());
      std::vector<int> obl;
      fp->getOnBits(obl);
      // there's an unfortunate bit of information loss happening here
      // due to the fact that the SparseBitVect uses ints, so we have to
      // do this test backwards:
      for (const auto &pr : *ao.bitInfoMap) {
        INFO(pr.first);
        CHECK(std::find(obl.begin(), obl.end(), (int)(pr.first)) != obl.end());
      }
      // for (auto i = 0u; i < m->getNumAtoms(); ++i) {
      //   std::cerr << " {";
      //   std::copy(ao.atomToBits->at(i).begin(), ao.atomToBits->at(i).end(),
      //             std::ostream_iterator<std::uint64_t>(std::cerr, ", "));
      //   std::cerr << " }," << std::endl;
      // }
      std::vector<unsigned int> atomCounts{2, 3, 3, 3, 3, 3, 3, 2};
      CHECK(*ao.atomCounts == atomCounts);
      std::vector<std::vector<std::uint64_t>> atomToBits = {
          {1845699452, 3458649244},
          {391602504, 391602505, 391602506, 3018396076, 3209717616},
          {391602504, 391602505, 391602506, 1746877920, 1746877921, 1746877922,
           3245328504},
          {391602504, 391602505, 391602506, 647852508, 647852509, 1746877920,
           1746877921, 1746877922},
          {391602504, 391602505, 391602506, 647852508, 647852509, 1746877920,
           1746877921, 1746877922},
          {391602504, 391602505, 391602506, 1746877920, 1746877921, 1746877922,
           4225765616},
          {10672060, 391602504, 391602505, 391602506, 3149364416},
          {1781206876, 3391828556}};
      CHECK(*ao.atomToBits == atomToBits);
    }
  }
}

TEST_CASE("numBitsPerFeature") {
  auto mol = "CC(=O)O"_smiles;
  REQUIRE(mol);
  SECTION("basics") {
    std::unique_ptr<FingerprintGenerator<std::uint32_t>> fpgen{
        MorganFingerprint::getMorganGenerator<std::uint32_t>(2)};
    REQUIRE(fpgen);
    {
      std::unique_ptr<SparseIntVect<std::uint32_t>> fp{
          fpgen->getCountFingerprint(*mol)};
      REQUIRE(fp);
      CHECK(fp->getTotalVal() == 8);
    }
    // turn on multiple bits per feature:
    fpgen->getOptions()->d_numBitsPerFeature = 2;
    {
      std::unique_ptr<SparseIntVect<std::uint32_t>> fp{
          fpgen->getCountFingerprint(*mol)};
      REQUIRE(fp);
      CHECK(fp->getTotalVal() == 16);
    }
    CHECK(fpgen->infoString().find("bitsPerFeature=2") != std::string::npos);
  }
}

TEST_CASE("multithreaded fp generation") {
  std::vector<std::string> smis = {"CC1CCC1", "CCC1CCC1", "CCCC1CCC1",
                                   "CCCC1CC(O)C1", "CCCC1CC(CO)C1"};
  std::vector<std::unique_ptr<ROMol>> ov;
  std::vector<const ROMol *> mols;
  for (const auto &smi : smis) {
    ov.emplace_back(SmilesToMol(smi));
    REQUIRE(ov.back());
    mols.push_back(ov.back().get());
  }
  for (auto i = 0u; i < 6; ++i) {
    auto n = mols.size();
    for (auto j = 0u; j < n; ++j) {
      mols.push_back(mols[j]);
    }
  }
  std::unique_ptr<FingerprintGenerator<std::uint32_t>> fpgen{
      MorganFingerprint::getMorganGenerator<std::uint32_t>(2)};
  REQUIRE(fpgen);
  SECTION("getFingerprints") {
    std::vector<std::unique_ptr<ExplicitBitVect>> ovs;
    FingerprintFuncArguments args;
    for (const auto mp : mols) {
      ovs.emplace_back(fpgen->getFingerprint(*mp, args));
    }
    mols.push_back(nullptr);  // make sure we handle this properly
    auto mtvs1 = fpgen->getFingerprints(mols, 1);
    CHECK(mtvs1.size() == ovs.size() + 1);
    for (auto fpi = 0u; fpi < ovs.size(); ++fpi) {
      CHECK(*ovs[fpi] == *mtvs1[fpi]);
    }
    CHECK(mtvs1.back().get() == nullptr);
#ifdef RDK_BUILD_THREADSAFE_SSS
    auto mtvs4 = fpgen->getFingerprints(mols, 4);
    CHECK(mtvs4.size() == ovs.size() + 1);
    for (auto fpi = 0u; fpi < ovs.size(); ++fpi) {
      CHECK(*ovs[fpi] == *mtvs4[fpi]);
    }
    CHECK(mtvs4.back().get() == nullptr);
#endif
  }
  SECTION("getSparseFingerprints") {
    std::vector<std::unique_ptr<SparseBitVect>> ovs;
    FingerprintFuncArguments args;
    for (const auto mp : mols) {
      ovs.emplace_back(fpgen->getSparseFingerprint(*mp, args));
    }
    mols.push_back(nullptr);  // make sure we handle this properly
    auto mtvs1 = fpgen->getSparseFingerprints(mols, 1);
    CHECK(mtvs1.size() == ovs.size() + 1);
    for (auto fpi = 0u; fpi < ovs.size(); ++fpi) {
      CHECK(*ovs[fpi] == *mtvs1[fpi]);
    }
    CHECK(mtvs1.back().get() == nullptr);
#ifdef RDK_BUILD_THREADSAFE_SSS
    auto mtvs4 = fpgen->getSparseFingerprints(mols, 4);
    CHECK(mtvs4.size() == ovs.size() + 1);
    for (auto fpi = 0u; fpi < ovs.size(); ++fpi) {
      CHECK(*ovs[fpi] == *mtvs4[fpi]);
    }
    CHECK(mtvs4.back().get() == nullptr);
#endif
  }
  SECTION("getCountFingerprints") {
    std::vector<std::unique_ptr<SparseIntVect<std::uint32_t>>> ovs;
    FingerprintFuncArguments args;
    for (const auto mp : mols) {
      ovs.emplace_back(fpgen->getCountFingerprint(*mp, args));
    }
    mols.push_back(nullptr);  // make sure we handle this properly
    auto mtvs1 = fpgen->getCountFingerprints(mols, 1);
    CHECK(mtvs1.size() == ovs.size() + 1);
    for (auto fpi = 0u; fpi < ovs.size(); ++fpi) {
      CHECK(*ovs[fpi] == *mtvs1[fpi]);
    }
    CHECK(mtvs1.back().get() == nullptr);
#ifdef RDK_BUILD_THREADSAFE_SSS
    auto mtvs4 = fpgen->getCountFingerprints(mols, 4);
    CHECK(mtvs4.size() == ovs.size() + 1);
    for (auto fpi = 0u; fpi < ovs.size(); ++fpi) {
      CHECK(*ovs[fpi] == *mtvs4[fpi]);
    }
    CHECK(mtvs4.back().get() == nullptr);
#endif
  }
  SECTION("getSparseCountFingerprints") {
    std::vector<std::unique_ptr<SparseIntVect<std::uint32_t>>> ovs;
    FingerprintFuncArguments args;
    for (const auto mp : mols) {
      ovs.emplace_back(fpgen->getSparseCountFingerprint(*mp, args));
    }
    mols.push_back(nullptr);  // make sure we handle this properly
    auto mtvs1 = fpgen->getSparseCountFingerprints(mols, 1);
    CHECK(mtvs1.size() == ovs.size() + 1);
    for (auto fpi = 0u; fpi < ovs.size(); ++fpi) {
      CHECK(*ovs[fpi] == *mtvs1[fpi]);
    }
    CHECK(mtvs1.back().get() == nullptr);
#ifdef RDK_BUILD_THREADSAFE_SSS
    auto mtvs4 = fpgen->getSparseCountFingerprints(mols, 4);
    CHECK(mtvs4.size() == ovs.size() + 1);
    for (auto fpi = 0u; fpi < ovs.size(); ++fpi) {
      CHECK(*ovs[fpi] == *mtvs4[fpi]);
    }
    CHECK(mtvs4.back().get() == nullptr);
#endif
  }
}

TEST_CASE("countBounds edge cases") {
  auto mol = "CC"_smiles;
  REQUIRE(mol);
  SECTION("just zeros") {
    std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGenerator(
        MorganFingerprint::getMorganGenerator<std::uint64_t>(2));
    REQUIRE(fpGenerator);
    fpGenerator->getOptions()->df_countSimulation = true;
    fpGenerator->getOptions()->d_countBounds = {0, 0, 0, 0};
    std::unique_ptr<ExplicitBitVect> fp(fpGenerator->getFingerprint(*mol));
    REQUIRE(fp);
    CHECK(fp->getNumBits() == 2048);
  }
  SECTION("empty bounds") {
    std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGenerator(
        MorganFingerprint::getMorganGenerator<std::uint64_t>(2));
    REQUIRE(fpGenerator);
    fpGenerator->getOptions()->df_countSimulation = true;
    fpGenerator->getOptions()->d_countBounds.clear();
    REQUIRE_THROWS_AS(fpGenerator->getFingerprint(*mol), ValueErrorException);
  }
  SECTION("really big") {
    std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGenerator(
        MorganFingerprint::getMorganGenerator<std::uint64_t>(2));
    REQUIRE(fpGenerator);
    fpGenerator->getOptions()->df_countSimulation = true;
    fpGenerator->getOptions()->d_countBounds =
        std::vector<unsigned int>((1 << 11) + 1, 0);
    REQUIRE_THROWS_AS(fpGenerator->getFingerprint(*mol), ValueErrorException);
  }
  SECTION("edge case") {
    std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGenerator(
        MorganFingerprint::getMorganGenerator<std::uint64_t>(2));
    REQUIRE(fpGenerator);
    fpGenerator->getOptions()->df_countSimulation = true;
    fpGenerator->getOptions()->d_countBounds =
        std::vector<unsigned int>(2047, 0);
    std::unique_ptr<ExplicitBitVect> fp(fpGenerator->getFingerprint(*mol));
    REQUIRE(fp);
    CHECK(fp->getNumBits() == 2048);
  }
}