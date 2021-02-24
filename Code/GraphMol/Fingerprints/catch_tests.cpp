//
//  Copyright (C) 2019 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "catch.hpp"

#include <memory>
#include <RDGeneral/test.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolBundle.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <RDGeneral/Exceptions.h>
#include <GraphMol/Fingerprints/RDKitFPGenerator.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>
#include <string>

using namespace RDKit;

TEST_CASE("Github 2051", "[patternfp][bug]") {
  auto mol = "CCC1CC1"_smiles;
  std::unique_ptr<ExplicitBitVect> mfp(PatternFingerprintMol(*mol));

  REQUIRE(mfp);
  SECTION("basics1") {
    auto qmol = "**"_smarts;
    std::unique_ptr<ExplicitBitVect> qfp(PatternFingerprintMol(*qmol));
    REQUIRE(qfp);

    CHECK(AllProbeBitsMatch(*qfp, *mfp));
  }
  SECTION("basics2") {
    auto qmol = "*"_smarts;
    std::unique_ptr<ExplicitBitVect> qfp(PatternFingerprintMol(*qmol));
    REQUIRE(qfp);
    CHECK(AllProbeBitsMatch(*qfp, *mfp));
  }
}

TEST_CASE("Github 2614", "[patternfp][bug]") {
  SECTION("basics") {
    auto mol =
        "F[P-](F)(F)(F)(F)F.F[P-](F)(F)(F)(F)F.F[P-](F)(F)(F)(F)F.F[P-](F)(F)(F)(F)F.F[P-](F)(F)(F)(F)F.F[P-](F)(F)(F)(F)F.F[P-](F)(F)(F)(F)F.F[P-](F)(F)(F)(F)F.F[P-](F)(F)(F)(F)F.c1ccc2ccccc2c1"_smiles;
    REQUIRE(mol);
    std::unique_ptr<ExplicitBitVect> mfp(PatternFingerprintMol(*mol));
    REQUIRE(mfp);
    auto query = "c1ccc2ccccc2c1"_smiles;
    REQUIRE(query);
    std::unique_ptr<ExplicitBitVect> qfp(PatternFingerprintMol(*query));
    REQUIRE(qfp);
    CHECK(AllProbeBitsMatch(*qfp, *mfp));
  }
}

TEST_CASE("Github 1761", "[patternfp][bug]") {
  SECTION("throw ValueErrorException") {
    auto mol = "CCC1CC1"_smiles;
    try {
      RDKit::MorganFingerprints::getHashedFingerprint(*mol, 0, 0);
      FAIL("Expected ValueErrorException");
    } catch (const ValueErrorException &e) {
      REQUIRE(std::string(e.what()) == "nBits can not be zero");
    }
  }
}

TEST_CASE("RDKit bits per feature", "[fpgenerator][rdkit]") {
  auto m1 = "CCCO"_smiles;
  REQUIRE(m1);
  SECTION("defaults") {
    unsigned int minPath = 1;
    unsigned int maxPath = 2;
    std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGenerator(
        RDKitFP::getRDKitFPGenerator<std::uint64_t>(minPath, maxPath));
    REQUIRE(fpGenerator);
    std::unique_ptr<ExplicitBitVect> fp(fpGenerator->getFingerprint(*m1));
    REQUIRE(fp);
    CHECK(fp->getNumBits() == 2048);
    CHECK(fp->getNumOnBits() == 8);
    CHECK(fpGenerator->infoString().find("bitsPerFeature=2") !=
          std::string::npos);
  }
  SECTION("change numBitsPerFeature") {
    // I won't lie: having to do this makes my head hurt, but fixing it to
    // create a ctor that takes a Parameters object is more effort than I can
    // devote at the moment
    unsigned int minPath = 1;
    unsigned int maxPath = 2;
    bool useHs = true;
    bool branchedPaths = true;
    bool useBondOrder = true;
    AtomInvariantsGenerator *atomInvariantsGenerator = nullptr;
    bool countSimulation = false;
    const std::vector<std::uint32_t> countBounds = {1, 2, 4, 8};
    std::uint32_t fpSize = 2048;
    std::uint32_t numBitsPerFeature = 1;
    std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGenerator(
        RDKitFP::getRDKitFPGenerator<std::uint64_t>(
            minPath, maxPath, useHs, branchedPaths, useBondOrder,
            atomInvariantsGenerator, countSimulation, countBounds, fpSize,
            numBitsPerFeature));
    REQUIRE(fpGenerator);
    std::unique_ptr<ExplicitBitVect> fp(fpGenerator->getFingerprint(*m1));
    REQUIRE(fp);
    CHECK(fp->getNumBits() == 2048);
    CHECK(fp->getNumOnBits() == 4);
    CHECK(fpGenerator->infoString().find("bitsPerFeature=1") !=
          std::string::npos);
  }
}

TEST_CASE("pattern fingerprints for MolBundles", "[patternfp]") {
  SECTION("basics") {
    boost::shared_ptr<ROMol> q1{SmilesToMol("OCCO")};
    REQUIRE(q1);
    boost::shared_ptr<ROMol> q2{SmilesToMol("OCCCO")};
    REQUIRE(q2);
    std::unique_ptr<ExplicitBitVect> pfp1{PatternFingerprintMol(*q1)};
    REQUIRE(pfp1);
    std::unique_ptr<ExplicitBitVect> pfp2{PatternFingerprintMol(*q2)};
    REQUIRE(pfp2);

    MolBundle bundle;
    bundle.addMol(q1);
    bundle.addMol(q2);
    std::unique_ptr<ExplicitBitVect> pfp{PatternFingerprintMol(bundle)};
    REQUIRE(pfp);
    CHECK(((*pfp1) & (*pfp2)).getNumOnBits() > 0);
    CHECK(((*pfp1) & (*pfp2)).getNumOnBits() == pfp->getNumOnBits());
    CHECK(((*pfp1) & (*pfp2)) == *pfp);
  }
}

TEST_CASE("MorganGenerator bit info", "[fpgenerator][morgan]") {
  auto m1 = "CCC(CC)CO"_smiles;
  REQUIRE(m1);
  unsigned radius = 1;
  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGenerator(
      MorganFingerprint::getMorganGenerator<std::uint64_t>(radius));
  REQUIRE(fpGenerator);
  const std::vector<std::uint32_t> *fromAtoms = nullptr;
  const std::vector<std::uint32_t> *ignoreAtoms = nullptr;
  const int confId = -1;

  SECTION("folded bitInfoMap") {
    AdditionalOutput::bitInfoMapType expected = {
        {1, {{2, 0}}}, {80, {{1, 0}, {3, 0}, {5, 0}}}, {294, {{0, 1}, {4, 1}}}};

    {
      AdditionalOutput ao;
      ao.allocateBitInfoMap(m1->getNumAtoms());
      std::unique_ptr<ExplicitBitVect> fp(fpGenerator->getFingerprint(
          *m1, fromAtoms, ignoreAtoms, confId, &ao));
      CHECK(fp->getNumOnBits() == ao.bitInfoMap->size());

      for (const auto &elem : expected) {
        CHECK((*ao.bitInfoMap)[elem.first] == elem.second);
        CHECK(fp->getBit(elem.first));
      }
    }
    {
      AdditionalOutput ao;
      ao.allocateBitInfoMap(m1->getNumAtoms());
      std::unique_ptr<SparseIntVect<std::uint32_t>> fp(
          fpGenerator->getCountFingerprint(*m1, fromAtoms, ignoreAtoms, confId,
                                           &ao));
      CHECK(fp->getNonzeroElements().size() == ao.bitInfoMap->size());

      for (const auto &elem : expected) {
        CHECK((*ao.bitInfoMap)[elem.first] == elem.second);
        CHECK(fp->getVal(elem.first));
      }
    }
  }

  SECTION("folded atomToBits atomCounts") {
    AdditionalOutput::atomToBitsType expected1 = {
        {1057, 294}, {80, 1544}, {1, 1420}, {80, 1544},
        {1057, 294}, {80, 482},  {807, 222}};
    AdditionalOutput::atomCountsType expected2 = {2, 2, 2, 2, 2, 2, 2};

    {
      AdditionalOutput ao;
      ao.allocateAtomCounts(m1->getNumAtoms());
      ao.allocateAtomToBits(m1->getNumAtoms());
      std::unique_ptr<ExplicitBitVect> fp(fpGenerator->getFingerprint(
          *m1, fromAtoms, ignoreAtoms, confId, &ao));
      CHECK(*ao.atomToBits == expected1);
      CHECK(*ao.atomCounts == expected2);
    }
    {
      AdditionalOutput ao;
      ao.allocateAtomCounts(m1->getNumAtoms());
      ao.allocateAtomToBits(m1->getNumAtoms());
      std::unique_ptr<SparseIntVect<std::uint32_t>> fp(
          fpGenerator->getCountFingerprint(*m1, fromAtoms, ignoreAtoms, confId,
                                           &ao));
      CHECK(*ao.atomToBits == expected1);
      CHECK(*ao.atomCounts == expected2);
    }
  }

  SECTION("unfolded bitInfoMap") {
    AdditionalOutput::bitInfoMapType expected = {
        {2245273601, {{2, 0}}},
        {2245384272, {{1, 0}, {3, 0}, {5, 0}}},
        {3542456614, {{0, 1}, {4, 1}}}};
    {
      AdditionalOutput ao;
      ao.allocateBitInfoMap(m1->getNumAtoms());
      std::unique_ptr<SparseBitVect> fp(fpGenerator->getSparseFingerprint(
          *m1, fromAtoms, ignoreAtoms, confId, &ao));
      CHECK(fp->getNumOnBits() == ao.bitInfoMap->size());

      for (const auto &elem : expected) {
        CHECK((*ao.bitInfoMap)[elem.first] == elem.second);
        CHECK(fp->getBit(elem.first));
      }
    }
    {
      AdditionalOutput ao;
      ao.allocateBitInfoMap(m1->getNumAtoms());
      std::unique_ptr<SparseIntVect<std::uint64_t>> fp(
          fpGenerator->getSparseCountFingerprint(*m1, fromAtoms, ignoreAtoms,
                                                 confId, &ao));
      CHECK(fp->getNonzeroElements().size() == ao.bitInfoMap->size());

      for (const auto &elem : expected) {
        CHECK((*ao.bitInfoMap)[elem.first] == elem.second);
        CHECK(fp->getVal(elem.first));
      }
    }
  }
  SECTION("unfolded atomToBits atomCounts") {
    AdditionalOutput::atomToBitsType expected1 = {
        {2246728737, 3542456614}, {2245384272, 1506563592},
        {2245273601, 3098934668}, {2245384272, 1506563592},
        {2246728737, 3542456614}, {2245384272, 4022716898},
        {864662311, 1535166686}};
    AdditionalOutput::atomCountsType expected2 = {2, 2, 2, 2, 2, 2, 2};

    {
      AdditionalOutput ao;
      ao.allocateAtomCounts(m1->getNumAtoms());
      ao.allocateAtomToBits(m1->getNumAtoms());
      std::unique_ptr<SparseBitVect> fp(fpGenerator->getSparseFingerprint(
          *m1, fromAtoms, ignoreAtoms, confId, &ao));
      CHECK(*ao.atomToBits == expected1);
      CHECK(*ao.atomCounts == expected2);
    }
    {
      AdditionalOutput ao;
      ao.allocateAtomCounts(m1->getNumAtoms());
      ao.allocateAtomToBits(m1->getNumAtoms());
      std::unique_ptr<SparseIntVect<std::uint64_t>> fp(
          fpGenerator->getSparseCountFingerprint(*m1, fromAtoms, ignoreAtoms,
                                                 confId, &ao));
      CHECK(*ao.atomToBits == expected1);
      CHECK(*ao.atomCounts == expected2);
    }
  }

  //   for (const auto &pr : *ao.bitInfoMap) {
  //     std::cerr << pr.first << std::endl;
  //     for (const auto &env : pr.second) {
  //       std::cerr << "   " << env.first << " " << env.second << std::endl;
  //     }
  //   }
  // }
}

TEST_CASE("RDKitGenerator bit info", "[fpgenerator][RDKit]") {
  auto m1 = "CCCO"_smiles;
  REQUIRE(m1);
  unsigned int minPath = 1;
  unsigned int maxPath = 3;
  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGenerator(
      RDKitFP::getRDKitFPGenerator<std::uint64_t>(minPath, maxPath));
  REQUIRE(fpGenerator);
  const std::vector<std::uint32_t> *fromAtoms = nullptr;
  const std::vector<std::uint32_t> *ignoreAtoms = nullptr;
  const int confId = -1;

  SECTION("folded bitInfo") {
    AdditionalOutput::bitPathsType expected = {{1233,
                                                {{
                                                    0,
                                                    1,
                                                    2,
                                                }}},
                                               {1308, {{0}, {1}}},
                                               {1339, {{2}}},
                                               {1728, {{1, 2}}},
                                               {1813,
                                                {{
                                                    0,
                                                    1,
                                                }}}};

    {
      AdditionalOutput ao;
      ao.allocateBitPaths(m1->getNumAtoms());
      std::unique_ptr<ExplicitBitVect> fp(fpGenerator->getFingerprint(
          *m1, fromAtoms, ignoreAtoms, confId, &ao));
      CHECK(*ao.bitPaths == expected);
    }
    {
      AdditionalOutput ao;
      ao.allocateBitPaths(m1->getNumAtoms());
      std::unique_ptr<SparseIntVect<std::uint32_t>> fp(
          fpGenerator->getCountFingerprint(*m1, fromAtoms, ignoreAtoms, confId,
                                           &ao));
      CHECK(*ao.bitPaths == expected);
    }
  }

  SECTION("folded atomToBits atomCounts") {
    AdditionalOutput::atomToBitsType expected1 = {
        {1308, 1813, 1233},
        {1308, 1308, 1813, 1728, 1233},
        {1308, 1339, 1813, 1728, 1233},
        {1339, 1728, 1233}};
    AdditionalOutput::atomCountsType expected2 = {3, 5, 5, 3};

    {
      AdditionalOutput ao;
      ao.allocateAtomCounts(m1->getNumAtoms());
      ao.allocateAtomToBits(m1->getNumAtoms());
      std::unique_ptr<ExplicitBitVect> fp(fpGenerator->getFingerprint(
          *m1, fromAtoms, ignoreAtoms, confId, &ao));
      CHECK(*ao.atomToBits == expected1);
      CHECK(*ao.atomCounts == expected2);
    }
    {
      AdditionalOutput ao;
      ao.allocateAtomCounts(m1->getNumAtoms());
      ao.allocateAtomToBits(m1->getNumAtoms());
      std::unique_ptr<SparseIntVect<std::uint32_t>> fp(
          fpGenerator->getCountFingerprint(*m1, fromAtoms, ignoreAtoms, confId,
                                           &ao));
      CHECK(*ao.atomToBits == expected1);
      CHECK(*ao.atomCounts == expected2);
    }
  }

  SECTION("unfolded bitInfo") {
    AdditionalOutput::bitPathsType expected = {{3977409745,
                                                {{
                                                    0,
                                                    1,
                                                    2,
                                                }}},
                                               {4275705116, {{0}, {1}}},
                                               {4274652475, {{2}}},
                                               {1524090560, {{1, 2}}},
                                               {1940446997,
                                                {{
                                                    0,
                                                    1,
                                                }}}};

    {
      AdditionalOutput ao;
      ao.allocateBitPaths(m1->getNumAtoms());
      std::unique_ptr<SparseBitVect> fp(fpGenerator->getSparseFingerprint(
          *m1, fromAtoms, ignoreAtoms, confId, &ao));
      CHECK(*ao.bitPaths == expected);
    }
    {
      AdditionalOutput ao;
      ao.allocateBitPaths(m1->getNumAtoms());
      std::unique_ptr<SparseIntVect<std::uint64_t>> fp(
          fpGenerator->getSparseCountFingerprint(*m1, fromAtoms, ignoreAtoms,
                                                 confId, &ao));
      CHECK(*ao.bitPaths == expected);
    }
  }

  SECTION("unfolded atomToBits atomCounts") {
    AdditionalOutput::atomToBitsType expected1 = {
        {4275705116, 1940446997, 3977409745},
        {4275705116, 4275705116, 1940446997, 1524090560, 3977409745},
        {4275705116, 4274652475, 1940446997, 1524090560, 3977409745},
        {4274652475, 1524090560, 3977409745}};
    AdditionalOutput::atomCountsType expected2 = {3, 5, 5, 3};

    {
      AdditionalOutput ao;
      ao.allocateAtomCounts(m1->getNumAtoms());
      ao.allocateAtomToBits(m1->getNumAtoms());
      std::unique_ptr<SparseBitVect> fp(fpGenerator->getSparseFingerprint(
          *m1, fromAtoms, ignoreAtoms, confId, &ao));
      CHECK(*ao.atomToBits == expected1);
      CHECK(*ao.atomCounts == expected2);
    }
    {
      AdditionalOutput ao;
      ao.allocateAtomCounts(m1->getNumAtoms());
      ao.allocateAtomToBits(m1->getNumAtoms());
      std::unique_ptr<SparseIntVect<std::uint64_t>> fp(
          fpGenerator->getSparseCountFingerprint(*m1, fromAtoms, ignoreAtoms,
                                                 confId, &ao));
      CHECK(*ao.atomToBits == expected1);
      CHECK(*ao.atomCounts == expected2);
    }
  }
}