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
#include <GraphMol/Substruct/SubstructMatch.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>
#include <string>

using namespace RDKit;

TEST_CASE("Github 2051", "[patternfp,bug]") {
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

TEST_CASE("Github 2614", "[patternfp,bug]") {
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
