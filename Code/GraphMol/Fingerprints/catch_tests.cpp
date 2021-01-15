//
//  Copyright (C) 2019 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do
                           // this in one cpp file
#include "catch.hpp"

#include <memory>
#include <RDGeneral/test.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <RDGeneral/Exceptions.h>
#include <GraphMol/Fingerprints/RDKitFPGenerator.h>
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
        try
        {
            RDKit::MorganFingerprints::getHashedFingerprint(*mol, 0, 0);
            FAIL("Expected ValueErrorException");
        }
        catch (const ValueErrorException &e)
        {
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
    // I won't lie: having to do this makes my head hurt, but fixing it to create a 
    // ctor that takes a Parameters object is more effort than I can devote at the moment
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