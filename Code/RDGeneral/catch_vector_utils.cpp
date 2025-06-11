//
//  Copyright (C) 2025 NVIDIA Corporation & Affiliates and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>

#include "vector_utils.h"

using namespace RDKit;

TEST_CASE("Test eraseMultipleIndices") {
  std::vector<int> vec = {0, 1, 2, 3, 4};
  SECTION("Empty indices") {
    std::vector<size_t> indices;
    std::vector<int> wantVec = vec;

    eraseMultipleIndices(vec, indices);
    CHECK(vec == wantVec);

    std::vector<double> empty;
    eraseMultipleIndices(empty, indices);
    CHECK(empty.empty());
  }

  SECTION("Out of order indices") {
    std::vector<size_t> indices = {3, 1};
    CHECK_THROWS_AS(eraseMultipleIndices(vec, indices), Invar::Invariant);
  }

  SECTION("Out of range indices") {
    std::vector<size_t> indices = {6};
    CHECK_THROWS_AS(eraseMultipleIndices(vec, indices), Invar::Invariant);
  }

  SECTION("Normal case") {
    std::vector<size_t> indices = {1, 3};
    std::vector<int> wantVec = {0, 2, 4};

    eraseMultipleIndices(vec, indices);
    CHECK(vec == wantVec);
  }

  SECTION("Edge case: Erase last element") {
    std::vector<size_t> indices = {1, 4};
    std::vector<int> wantVec = {0, 2, 3};
    eraseMultipleIndices(vec, indices);
    CHECK(vec == wantVec);
  }

  SECTION("Edge case: Erase first element") {
    std::vector<size_t> indices = {0, 3};
    std::vector<int> wantVec = {1, 2, 4};
    eraseMultipleIndices(vec, indices);
    CHECK(vec == wantVec);
  }

  SECTION("Edge case: erase only one element") {
    std::vector<size_t> indices = {2};
    std::vector<int> wantVec = {0, 1, 3, 4};
    eraseMultipleIndices(vec, indices);
    CHECK(vec == wantVec);
  }

  SECTION("Edge case: Erase adjacent elements") {
    std::vector<size_t> indices = {1, 2};
    std::vector<int> wantVec = {0, 3, 4};
    eraseMultipleIndices(vec, indices);
    CHECK(vec == wantVec);
  }

  SECTION("Edge case: size 1 input") {
    std::vector<int> input = {0};
    std::vector<size_t> indices = {0};

    eraseMultipleIndices(input, std::vector<size_t>{});
    std::vector<int> wantFromEmptyIndices = input;
    REQUIRE(input == wantFromEmptyIndices);
    eraseMultipleIndices(input, indices);
    REQUIRE(input == std::vector<int>{});
  }

  SECTION("Edge case: all elements erased") {
    std::vector<size_t> indices = {0, 1, 2, 3, 4};
    std::vector<int> wantVec = {};

    eraseMultipleIndices(vec, indices);
    CHECK(vec == wantVec);
  }
}