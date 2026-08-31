//
//  Copyright (C) 2025 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <set>
#include <catch2/catch_all.hpp>

#include "MacrocycleGenerator.h"

using namespace RDDepict;

TEST_CASE("macrocycle self-crossing detection", "[macrocycle]") {
  SECTION("solved macrocycles have unique hex-grid positions") {
    for (const size_t ringSize : {12, 14}) {
      CAPTURE(ringSize);
      MacrocycleGenerator generator(ringSize);
      REQUIRE(generator.solve());

      const auto coords = generator.calculateHexCoords(generator.getTurns());
      const std::set<HexCoord> uniqueCoords(coords.begin(), coords.end() - 1);
      CHECK(uniqueCoords.size() == ringSize);
      CHECK(generator.calculateMinDistance(generator.getTurns()) > 0);
    }
  }

  SECTION("a crossing sequence is detected") {
    MacrocycleGenerator generator(12);
    const std::vector<int> crossingPattern(12, 1);
    CHECK(generator.calculateMinDistance(crossingPattern) == 0);
  }
}
