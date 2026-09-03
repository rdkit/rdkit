//
//  Copyright (C) 2025 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include "MacrocycleGenerator.h"

using namespace RDDepict;

TEST_CASE("macrocycle turn constraints", "[macrocycle]") {
  SECTION("fixed turn patterns are preserved") {
    MacrocycleGenerator generator(14);
    generator.addConstraint(
        {0, ConstraintType::FIXED, {1, -1, 1}, "three-atom fusion"});

    REQUIRE(generator.solve());
    const auto &turns = generator.getTurns();
    CHECK(turns[0] == 1);
    CHECK(turns[1] == -1);
    CHECK(turns[2] == 1);
  }

  SECTION("same and opposite constraints are preserved") {
    MacrocycleGenerator generator(16);
    generator.addConstraint({2, ConstraintType::SAME, {}, "cis bond"});
    generator.addConstraint(
        {7, ConstraintType::OPPOSITE, {}, "trans bond"});

    REQUIRE(generator.solve());
    const auto &turns = generator.getTurns();
    CHECK(turns[2] == turns[3]);
    CHECK(turns[7] == -turns[8]);
  }

  SECTION("conflicting constraints fail") {
    MacrocycleGenerator generator(12);
    generator.addConstraint({3, ConstraintType::FIXED, {1}, "force right"});
    generator.addConstraint({3, ConstraintType::FIXED, {-1}, "force left"});
    CHECK_FALSE(generator.solve());
  }

  SECTION("an over-constrained turn sequence fails") {
    MacrocycleGenerator generator(12);
    for (size_t position = 0; position < 12; ++position) {
      generator.addConstraint(
          {position, ConstraintType::FIXED, {1}, "force all right"});
    }
    CHECK_FALSE(generator.solve());
  }
}

TEST_CASE("unsupported shared endpoint topology is ignored", "[macrocycle]") {
  // A spiro ring has the same first and last macrocycle endpoint, so endpoint
  // tracking can present the same ring twice. Its full atom overlap is neither
  // a one-atom nor a two-atom fusion and must not create a zero-angle
  // constraint.
  const EndpointInfo spiroEndpoint{5, {0, 1, 2, 3, 4}};
  const auto constraint =
      computeSharedEndpointConstraint({spiroEndpoint, spiroEndpoint}, 3);

  CHECK(constraint.position == SIZE_MAX);
}
