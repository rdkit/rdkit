//
//  Copyright (C) 2024 Greg Landrum and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include <RDGeneral/test.h>
#include "BitVects.h"
#include "BitOps.h"
#include "BitVectUtils.h"
#include "ExplicitBitVect.h"
#include "SparseIntVect.h"
#include <limits>

using namespace RDKit;

TEST_CASE("special cases for the limits of sparse vectors") {
  SECTION("SparseBitVect") {
    SparseBitVect sbv(std::numeric_limits<unsigned int>::max());
    CHECK(sbv.getNumBits() == std::numeric_limits<unsigned int>::max());
    CHECK(!sbv.setBit(std::numeric_limits<unsigned int>::max()));
    CHECK(sbv.getBit(std::numeric_limits<unsigned int>::max()) == 1);
  }
}

TEST_CASE("github #9033: tversky is 1 when no bits are set") {
  ExplicitBitVect bv1(8);
  ExplicitBitVect bv2(8);
  CHECK(TverskySimilarity(bv1, bv2, 0.5, 0.5) == 0.0);
  CHECK(TanimotoSimilarity(bv1, bv2) == 0.0);
  CHECK(CosineSimilarity(bv1, bv2) == 0.0);
  CHECK(KulczynskiSimilarity(bv1, bv2) == 0.0);
  CHECK(SokalSimilarity(bv1, bv2) == 0.0);
  CHECK(McConnaugheySimilarity(bv1, bv2) == 0.0);
  CHECK(BraunBlanquetSimilarity(bv1, bv2) == 0.0);
  CHECK(RusselSimilarity(bv1, bv2) == 0.0);
  CHECK(RogotGoldbergSimilarity(bv1, bv2) == 0.0);
}