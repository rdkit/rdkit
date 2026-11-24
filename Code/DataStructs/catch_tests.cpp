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
#include "base64.h"
#include <RDGeneral/Exceptions.h>
#include <limits>
#include <memory>
#include <string>
#include <vector>

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

TEST_CASE("Base64Decode and badly padded input") {
  SECTION("properly padded input decodes") {
    struct Case {
      std::string in;
      std::string out;
    };
    std::vector<Case> cases{{"", ""},
                            {"TW==", "M"},
                            {"TWE=", "Ma"},
                            {"TWFu", "Man"},
                            {"TWFuTWE=", "ManMa"},
                            // 0xff is junk and has to be skipped: transTable
                            // needs to be initialized all the way to 255
                            {std::string("TW") + char(0xff) + "E=", "Ma"}};
    for (const auto &c : cases) {
      INFO(c.in);
      unsigned int len = 0;
      std::unique_ptr<char[]> res(Base64Decode(c.in.c_str(), &len));
      CHECK(std::string(res.get(), len) == c.out);
    }
  }
  SECTION("input which isn't properly padded throws") {
    unsigned int len = 0;
    for (const auto &in : {"T", "TW", "TWE", "TWFuTWE", "TWFu=", "TW="}) {
      INFO(in);
      CHECK_THROWS_AS(Base64Decode(in, &len), ValueErrorException);
    }
  }
}
