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

TEST_CASE("basic macrocycle coordinate generation", "[macrocycle]") {
  SECTION("even macrocycles close exactly") {
    for (const size_t ringSize : {10, 12, 14}) {
      CAPTURE(ringSize);
      MacrocycleGenerator generator(ringSize);
      REQUIRE(generator.solve());
      CHECK(generator.getClosureError() == Catch::Approx(0.0).margin(1e-6));

      const auto coords = generator.generateCoordinates();
      CHECK(coords.size() == ringSize);
      CHECK(generator.calculateMinDistance(generator.getTurns()) > 0);
    }
  }

  SECTION("odd macrocycles are closed during coordinate generation") {
    constexpr size_t ringSize = 13;
    MacrocycleGenerator generator(ringSize);
    REQUIRE(generator.solve());
    CHECK(generator.getClosureError() == Catch::Approx(1.5).margin(1e-6));

    const auto coords = generator.generateCoordinates();
    CHECK(coords.size() == ringSize);
  }

  SECTION("coordinates require a solved turn sequence") {
    MacrocycleGenerator generator(12);
    CHECK_THROWS_AS(generator.generateCoordinates(), std::runtime_error);
  }
}
