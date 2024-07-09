//
//  Copyright (C) 2019-2021 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

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
#include <GraphMol/Fingerprints/TopologicalTorsionGenerator.h>
#include <GraphMol/Fingerprints/AtomPairGenerator.h>
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
    unsigned int minPath = 1;
    unsigned int maxPath = 2;
    std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGenerator(
        RDKitFP::getRDKitFPGenerator<std::uint64_t>(minPath, maxPath));
    REQUIRE(fpGenerator);
    fpGenerator->getOptions()->d_numBitsPerFeature = 1;
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
      ao.allocateBitInfoMap();
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
      ao.allocateBitInfoMap();
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
      ao.allocateAtomCounts();
      ao.allocateAtomToBits();
      std::unique_ptr<ExplicitBitVect> fp(fpGenerator->getFingerprint(
          *m1, fromAtoms, ignoreAtoms, confId, &ao));
      CHECK(*ao.atomToBits == expected1);
      CHECK(*ao.atomCounts == expected2);
    }
    {
      AdditionalOutput ao;
      ao.allocateAtomCounts();
      ao.allocateAtomToBits();
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
      ao.allocateBitInfoMap();
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
      ao.allocateBitInfoMap();
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
      ao.allocateAtomCounts();
      ao.allocateAtomToBits();
      std::unique_ptr<SparseBitVect> fp(fpGenerator->getSparseFingerprint(
          *m1, fromAtoms, ignoreAtoms, confId, &ao));
      CHECK(*ao.atomToBits == expected1);
      CHECK(*ao.atomCounts == expected2);
    }
    {
      AdditionalOutput ao;
      ao.allocateAtomCounts();
      ao.allocateAtomToBits();
      std::unique_ptr<SparseIntVect<std::uint64_t>> fp(
          fpGenerator->getSparseCountFingerprint(*m1, fromAtoms, ignoreAtoms,
                                                 confId, &ao));
      CHECK(*ao.atomToBits == expected1);
      CHECK(*ao.atomCounts == expected2);
    }
  }
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
    // clang-format off
    AdditionalOutput::bitPathsType expected = {
          {562, {{2}}},
          {709, {{0, 1}}},
          {1118, {{0, 1, 2}}},
          {1183, {{1, 2}}},
          {1233, {{0, 1, 2}}},
          {1308, {{0}, {1}}},
          {1339, {{2}}},
          {1728, {{1, 2}}},
          {1772, {{0}, {1}}},
          {1813, {{0, 1}}},
        };
    // clang-format on

    {
      AdditionalOutput ao;
      ao.allocateBitPaths();
      std::unique_ptr<ExplicitBitVect> fp(fpGenerator->getFingerprint(
          *m1, fromAtoms, ignoreAtoms, confId, &ao));

      CHECK(*ao.bitPaths == expected);
    }
    {
      AdditionalOutput ao;
      ao.allocateBitPaths();
      std::unique_ptr<SparseIntVect<std::uint32_t>> fp(
          fpGenerator->getCountFingerprint(*m1, fromAtoms, ignoreAtoms, confId,
                                           &ao));
      CHECK(*ao.bitPaths == expected);
    }
  }

  SECTION("folded atomToBits atomCounts") {
    AdditionalOutput::atomToBitsType expected1 = {
        {1308, 1772, 1813, 709, 1233, 1118},
        {1308, 1772, 1813, 709, 1728, 1183, 1233, 1118},
        {1308, 1772, 1339, 562, 1813, 709, 1728, 1183, 1233, 1118},
        {1339, 562, 1728, 1183, 1233, 1118}};
    AdditionalOutput::atomCountsType expected2 = {6, 10, 10, 6};

    {
      AdditionalOutput ao;
      ao.allocateAtomCounts();
      ao.allocateAtomToBits();
      std::unique_ptr<ExplicitBitVect> fp(fpGenerator->getFingerprint(
          *m1, fromAtoms, ignoreAtoms, confId, &ao));
      CHECK(*ao.atomToBits == expected1);
      CHECK(*ao.atomCounts == expected2);
    }
    {
      AdditionalOutput ao;
      ao.allocateAtomCounts();
      ao.allocateAtomToBits();
      std::unique_ptr<SparseIntVect<std::uint32_t>> fp(
          fpGenerator->getCountFingerprint(*m1, fromAtoms, ignoreAtoms, confId,
                                           &ao));
      CHECK(*ao.atomToBits == expected1);
      CHECK(*ao.atomCounts == expected2);
    }
  }

  SECTION("unfolded bitInfo") {
    // clang-format off
    AdditionalOutput::bitPathsType expected = {
        {54413874, {{2}}},
        {257418334, {{0, 1, 2}}},
        {986785516, {{0}, {1}}},
        {1135572127, {{1, 2}}},
        {1433230021, {{0, 1}}},
        {1524090560, {{1, 2}}},
        {1940446997, {{0, 1}}},
        {3977409745, {{0, 1, 2}}},
        {4274652475, {{2}}},
        {4275705116, {{0}, {1}}},
      };
    // clang-format on

    {
      AdditionalOutput ao;
      ao.allocateBitPaths();
      std::unique_ptr<SparseBitVect> fp(fpGenerator->getSparseFingerprint(
          *m1, fromAtoms, ignoreAtoms, confId, &ao));
      CHECK(*ao.bitPaths == expected);
    }
    {
      AdditionalOutput ao;
      ao.allocateBitPaths();
      std::unique_ptr<SparseIntVect<std::uint64_t>> fp(
          fpGenerator->getSparseCountFingerprint(*m1, fromAtoms, ignoreAtoms,
                                                 confId, &ao));
      CHECK(*ao.bitPaths == expected);
    }
  }

  SECTION("unfolded atomToBits atomCounts") {
    AdditionalOutput::atomToBitsType expected1 = {
        {4275705116, 986785516, 1940446997, 1433230021, 3977409745, 257418334},
        {4275705116, 986785516, 1940446997, 1433230021, 1524090560, 1135572127,
         3977409745, 257418334},
        {4275705116, 986785516, 4274652475, 54413874, 1940446997, 1433230021,
         1524090560, 1135572127, 3977409745, 257418334},
        {4274652475, 54413874, 1524090560, 1135572127, 3977409745, 257418334}};
    AdditionalOutput::atomCountsType expected2 = {6, 10, 10, 6};

    {
      AdditionalOutput ao;
      ao.allocateAtomCounts();
      ao.allocateAtomToBits();
      std::unique_ptr<SparseBitVect> fp(fpGenerator->getSparseFingerprint(
          *m1, fromAtoms, ignoreAtoms, confId, &ao));
      CHECK(*ao.atomToBits == expected1);
      CHECK(*ao.atomCounts == expected2);
    }
    {
      AdditionalOutput ao;
      ao.allocateAtomCounts();
      ao.allocateAtomToBits();
      std::unique_ptr<SparseIntVect<std::uint64_t>> fp(
          fpGenerator->getSparseCountFingerprint(*m1, fromAtoms, ignoreAtoms,
                                                 confId, &ao));
      CHECK(*ao.atomToBits == expected1);
      CHECK(*ao.atomCounts == expected2);
    }
  }
}

TEST_CASE("TopologicalTorsionGenerator bit info", "[fpgenerator][TT]") {
  auto m1 = "CCCCS"_smiles;
  REQUIRE(m1);
  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGenerator(
      TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>());
  REQUIRE(fpGenerator);
  const std::vector<std::uint32_t> *fromAtoms = nullptr;
  const std::vector<std::uint32_t> *ignoreAtoms = nullptr;
  const int confId = -1;

  SECTION("folded bitInfo") {
    {
      // clang-format off
      AdditionalOutput::bitPathsType expected = {
          {0, {{0, 1, 2, 3}}},
          {1536, {{1, 2, 3, 4}}}
      };
      // clang-format on
      AdditionalOutput ao;
      ao.allocateBitPaths();
      std::unique_ptr<ExplicitBitVect> fp(fpGenerator->getFingerprint(
          *m1, fromAtoms, ignoreAtoms, confId, &ao));
      CHECK(*ao.bitPaths == expected);
    }
    {
      // clang-format off
      AdditionalOutput::bitPathsType expected = {
          {0, {{0, 1, 2, 3}}},
          {1920, {{1, 2, 3, 4}}}
      };
      // clang-format on

      AdditionalOutput ao;
      ao.allocateBitPaths();
      std::unique_ptr<SparseIntVect<std::uint32_t>> fp(
          fpGenerator->getCountFingerprint(*m1, fromAtoms, ignoreAtoms, confId,
                                           &ao));
      CHECK(*ao.bitPaths == expected);
    }
  }

  SECTION("folded atomToBits atomCounts") {
    AdditionalOutput::atomCountsType expected2 = {1, 2, 2, 2, 1};

    {
      AdditionalOutput::atomToBitsType expected1 = {
          {0}, {0, 1536}, {0, 1536}, {0, 1536}, {1536}};
      AdditionalOutput ao;
      ao.allocateAtomCounts();
      ao.allocateAtomToBits();
      std::unique_ptr<ExplicitBitVect> fp(fpGenerator->getFingerprint(
          *m1, fromAtoms, ignoreAtoms, confId, &ao));
      CHECK(*ao.atomToBits == expected1);
      CHECK(*ao.atomCounts == expected2);
    }
    {
      AdditionalOutput::atomToBitsType expected1 = {
          {0}, {0, 1920}, {0, 1920}, {0, 1920}, {1920}};
      AdditionalOutput ao;
      ao.allocateAtomCounts();
      ao.allocateAtomToBits();
      std::unique_ptr<SparseIntVect<std::uint32_t>> fp(
          fpGenerator->getCountFingerprint(*m1, fromAtoms, ignoreAtoms, confId,
                                           &ao));
      CHECK(*ao.atomToBits == expected1);
      CHECK(*ao.atomCounts == expected2);
    }
  }

  SECTION("unfolded bitInfo") {
    {
      // clang-format off
      AdditionalOutput::bitPathsType expected = {
          {1046740484, {{1, 2, 3, 4}}},
          {1048313860, {{0, 1, 2, 3}}},
      };
      // clang-format on
      AdditionalOutput ao;
      ao.allocateBitPaths();
      std::unique_ptr<SparseBitVect> fp(fpGenerator->getSparseFingerprint(
          *m1, fromAtoms, ignoreAtoms, confId, &ao));
      CHECK(*ao.bitPaths == expected);
    }
    {
      // clang-format off
      AdditionalOutput::bitPathsType expected = {
          {4437590048, {{0, 1, 2, 3}}},
          {30073176097, {{1, 2, 3, 4}}}
      };
      // clang-format on
      AdditionalOutput ao;
      ao.allocateBitPaths();
      std::unique_ptr<SparseIntVect<std::uint64_t>> fp(
          fpGenerator->getSparseCountFingerprint(*m1, fromAtoms, ignoreAtoms,
                                                 confId, &ao));

      CHECK(*ao.bitPaths == expected);
    }
  }

  SECTION("unfolded atomToBits atomCounts") {
    AdditionalOutput::atomCountsType expected2 = {1, 2, 2, 2, 1};

    {
      AdditionalOutput::atomToBitsType expected1 = {{1048313860},
                                                    {1046740484, 1048313860},
                                                    {1046740484, 1048313860},
                                                    {1046740484, 1048313860},
                                                    {1046740484}};
      AdditionalOutput ao;
      ao.allocateAtomCounts();
      ao.allocateAtomToBits();
      std::unique_ptr<SparseBitVect> fp(fpGenerator->getSparseFingerprint(
          *m1, fromAtoms, ignoreAtoms, confId, &ao));
      CHECK(*ao.atomToBits == expected1);
      CHECK(*ao.atomCounts == expected2);
    }
    {
      AdditionalOutput::atomToBitsType expected1 = {{4437590048},
                                                    {4437590048, 30073176097},
                                                    {4437590048, 30073176097},
                                                    {4437590048, 30073176097},
                                                    {30073176097}};
      AdditionalOutput ao;
      ao.allocateAtomCounts();
      ao.allocateAtomToBits();
      std::unique_ptr<SparseIntVect<std::uint64_t>> fp(
          fpGenerator->getSparseCountFingerprint(*m1, fromAtoms, ignoreAtoms,
                                                 confId, &ao));
      CHECK(*ao.atomToBits == expected1);
      CHECK(*ao.atomCounts == expected2);
    }
  }
}

TEST_CASE("AtomPairGenerator bit info", "[fpgenerator][AP]") {
  auto m1 = "CCO"_smiles;
  REQUIRE(m1);
  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGenerator(
      AtomPair::getAtomPairGenerator<std::uint64_t>());
  REQUIRE(fpGenerator);
  const std::vector<std::uint32_t> *fromAtoms = nullptr;
  const std::vector<std::uint32_t> *ignoreAtoms = nullptr;
  const int confId = -1;

  SECTION("folded bitInfo") {
    {
      AdditionalOutput::bitInfoMapType expected = {
          {1404, {{0, 1}}},
          {1916, {{0, 2}}},
          {1596, {{1, 2}}},
      };
      AdditionalOutput ao;
      ao.allocateBitInfoMap();
      std::unique_ptr<ExplicitBitVect> fp(fpGenerator->getFingerprint(
          *m1, fromAtoms, ignoreAtoms, confId, &ao));
      CHECK(*ao.bitInfoMap == expected);
    }
    {
      AdditionalOutput::bitInfoMapType expected = {
          {1375, {{0, 1}}},
          {1503, {{0, 2}}},
          {1423, {{1, 2}}},
      };
      AdditionalOutput ao;
      ao.allocateBitInfoMap();
      std::unique_ptr<SparseIntVect<std::uint32_t>> fp(
          fpGenerator->getCountFingerprint(*m1, fromAtoms, ignoreAtoms, confId,
                                           &ao));
      CHECK(*ao.bitInfoMap == expected);
    }
  }

  SECTION("folded atomToBits atomCounts") {
    AdditionalOutput::atomCountsType expected2 = {2, 2, 2};

    {
      AdditionalOutput::atomToBitsType expected1 = {
          {1404, 1916},
          {1404, 1596},
          {1596, 1916},
      };
      AdditionalOutput ao;
      ao.allocateAtomCounts();
      ao.allocateAtomToBits();
      std::unique_ptr<ExplicitBitVect> fp(fpGenerator->getFingerprint(
          *m1, fromAtoms, ignoreAtoms, confId, &ao));
      CHECK(*ao.atomToBits == expected1);
      CHECK(*ao.atomCounts == expected2);
    }
    {
      AdditionalOutput::atomToBitsType expected1 = {
          {1375, 1503},
          {1375, 1423},
          {1503, 1423},
      };
      AdditionalOutput ao;
      ao.allocateAtomCounts();
      ao.allocateAtomToBits();
      std::unique_ptr<SparseIntVect<std::uint32_t>> fp(
          fpGenerator->getCountFingerprint(*m1, fromAtoms, ignoreAtoms, confId,
                                           &ao));
      CHECK(*ao.atomToBits == expected1);
      CHECK(*ao.atomCounts == expected2);
    }
  }

  SECTION("unfolded bitInfo") {
    {
      AdditionalOutput::bitInfoMapType expected = {
          {7918972, {{0, 1}}},
          {7919484, {{0, 2}}},
          {8066620, {{1, 2}}},
      };
      AdditionalOutput ao;
      ao.allocateBitInfoMap();
      std::unique_ptr<SparseBitVect> fp(fpGenerator->getSparseFingerprint(
          *m1, fromAtoms, ignoreAtoms, confId, &ao));
      CHECK(*ao.bitInfoMap == expected);
    }
    {
      AdditionalOutput::bitInfoMapType expected = {
          {558113, {{0, 1}}},
          {1590306, {{0, 2}}},
          {1590337, {{1, 2}}},
      };
      AdditionalOutput ao;
      ao.allocateBitInfoMap();
      std::unique_ptr<SparseIntVect<std::uint64_t>> fp(
          fpGenerator->getSparseCountFingerprint(*m1, fromAtoms, ignoreAtoms,
                                                 confId, &ao));
      CHECK(*ao.bitInfoMap == expected);
    }
  }

  SECTION("unfolded atomToBits atomCounts") {
    AdditionalOutput::atomCountsType expected2 = {2, 2, 2};

    {
      AdditionalOutput::atomToBitsType expected1 = {
          {7918972, 7919484}, {7918972, 8066620}, {7919484, 8066620}};
      AdditionalOutput ao;
      ao.allocateAtomCounts();
      ao.allocateAtomToBits();
      std::unique_ptr<SparseBitVect> fp(fpGenerator->getSparseFingerprint(
          *m1, fromAtoms, ignoreAtoms, confId, &ao));
      CHECK(*ao.atomToBits == expected1);
      CHECK(*ao.atomCounts == expected2);
    }
    {
      AdditionalOutput::atomToBitsType expected1 = {
          {558113, 1590306}, {558113, 1590337}, {1590306, 1590337}};
      AdditionalOutput ao;
      ao.allocateAtomCounts();
      ao.allocateAtomToBits();
      std::unique_ptr<SparseIntVect<std::uint64_t>> fp(
          fpGenerator->getSparseCountFingerprint(*m1, fromAtoms, ignoreAtoms,
                                                 confId, &ao));
      CHECK(*ao.atomToBits == expected1);
      CHECK(*ao.atomCounts == expected2);
    }
  }
}

TEST_CASE("Generators, bit info, and multiple calls", "[fpgenerator]") {
  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGenerator(
      AtomPair::getAtomPairGenerator<std::uint64_t>());
  REQUIRE(fpGenerator);
  const std::vector<std::uint32_t> *fromAtoms = nullptr;
  const std::vector<std::uint32_t> *ignoreAtoms = nullptr;
  const int confId = -1;

  AdditionalOutput ao;
  ao.allocateBitInfoMap();
  ao.allocateAtomCounts();
  ao.allocateAtomToBits();

  auto m1 = "CCO"_smiles;
  REQUIRE(m1);
  std::unique_ptr<ExplicitBitVect> fp1(
      fpGenerator->getFingerprint(*m1, fromAtoms, ignoreAtoms, confId, &ao));
  CHECK(ao.bitInfoMap->size() == fp1->getNumOnBits());
  CHECK(ao.atomCounts->size() == m1->getNumAtoms());
  CHECK(ao.atomToBits->size() == m1->getNumAtoms());

  auto m2 = "CCON"_smiles;
  REQUIRE(m2);
  fp1.reset(
      fpGenerator->getFingerprint(*m2, fromAtoms, ignoreAtoms, confId, &ao));
  CHECK(ao.bitInfoMap->size() == fp1->getNumOnBits());
  CHECK(ao.atomCounts->size() == m2->getNumAtoms());
  CHECK(ao.atomToBits->size() == m2->getNumAtoms());

  fp1.reset(
      fpGenerator->getFingerprint(*m1, fromAtoms, ignoreAtoms, confId, &ao));
  CHECK(ao.bitInfoMap->size() == fp1->getNumOnBits());
  CHECK(ao.atomCounts->size() == m1->getNumAtoms());
  CHECK(ao.atomToBits->size() == m1->getNumAtoms());
}

TEST_CASE(
    "github #4212: UnfoldedRDKFingerprintCountBased returns a different "
    "fingerprint length for every molecule") {
  auto m1 = "c1ccccc1"_smiles;
  REQUIRE(m1);
  auto m2 = "CCCC"_smiles;
  REQUIRE(m2);
  std::unique_ptr<SparseIntVect<boost::uint64_t>> fp1{
      getUnfoldedRDKFingerprintMol(*m1)};
  REQUIRE(fp1);
  std::unique_ptr<SparseIntVect<boost::uint64_t>> fp2{
      getUnfoldedRDKFingerprintMol(*m2)};
  REQUIRE(fp2);
  CHECK(fp1->getLength() == fp2->getLength());
}

TEST_CASE("RDKit set countBounds", "[fpgenerator][rdkit]") {
  auto m1 = "COc1ccc(CCNC(=O)c2ccccc2C(=O)NCCc2ccc(OC)cc2)cc1"_smiles;
  REQUIRE(m1);
  SECTION("change countBounds") {
    unsigned int minPath = 1;
    unsigned int maxPath = 7;
    std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGenerator(
        RDKitFP::getRDKitFPGenerator<std::uint64_t>(minPath, maxPath));
    REQUIRE(fpGenerator);
    fpGenerator->getOptions()->df_countSimulation = true;
    std::unique_ptr<ExplicitBitVect> fp1(fpGenerator->getFingerprint(*m1));
    REQUIRE(fp1);
    CHECK(fp1->getNumBits() == 2048);

    fpGenerator->getOptions()->d_countBounds = {2, 8, 16, 32};
    std::unique_ptr<ExplicitBitVect> fp2(fpGenerator->getFingerprint(*m1));
    REQUIRE(fp2);
    CHECK(fp2->getNumBits() == 2048);
    CHECK(fp2->getNumOnBits() != fp1->getNumOnBits());
  }
}

TEST_CASE(
    "github #5036: Neighboring Hs not taken into account in connectivity "
    "invariants",
    "[morgan]") {
  SECTION("basics") {
    auto mol = "CC[2H]"_smiles;
    REQUIRE(mol);
    std::vector<std::uint32_t> invars(mol->getNumAtoms());
    MorganFingerprints::getConnectivityInvariants(*mol, invars);
    CHECK(invars[1] == invars[0]);
  }
}

TEST_CASE("topological torsions shorted paths") {
  SECTION("basics") {
    auto mol = "CC1CCC1"_smiles;
    REQUIRE(mol);
    std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGenerator(
        TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>());
    REQUIRE(fpGenerator);
    static_cast<TopologicalTorsion::TopologicalTorsionArguments *>(
        fpGenerator->getOptions())
        ->df_countSimulation = false;
    std::unique_ptr<SparseBitVect> fp(fpGenerator->getSparseFingerprint(*mol));
    CHECK(fp->getNumOnBits() == 3);
    static_cast<TopologicalTorsion::TopologicalTorsionArguments *>(
        fpGenerator->getOptions())
        ->df_onlyShortestPaths = true;
    fp.reset(fpGenerator->getSparseFingerprint(*mol));
    CHECK(fp->getNumOnBits() == 1);
  }
}

TEST_CASE(
    "GitHub #7318: Utils.AtomPairs.NumPiElectrons fails on atoms with dative bonds",
    "[bug]") {
  auto mol =
      "O=C1[O-]->[Cr+3]23(<-[O-]C(=O)C4=CC=CC=N->24)(<-[O-]C(=O)C2=CC=CC=N->32)<-N2=CC=CC=C12"_smiles;
  REQUIRE(mol);

  for (auto bond_idx : {2, 3, 12, 21, 28, 31}) {
    INFO("bond = " << bond_idx);
    auto bond = mol->getBondWithIdx(bond_idx);
    REQUIRE(bond->getBondType() == Bond::DATIVE);

    const auto atom = bond->getBeginAtom();
    CHECK(atom->getHybridization() == Atom::SP2);

    if (atom->getAtomicNum() == 8) {
      CHECK(numPiElectrons(*atom) == 0);
    } else {
      CHECK(numPiElectrons(*atom) == 1);
    }
  }
}

TEST_CASE("github #7533: IndexError while computing fingerprint") {
  auto mol = "CC1C(B(C)C)S1(C)(C)=O"_smiles;
  REQUIRE(mol);
  auto fp = MorganFingerprints::getFingerprint(*mol, 2);
  REQUIRE(fp);
  CHECK(fp->getLength() == std::numeric_limits<unsigned>::max());
}