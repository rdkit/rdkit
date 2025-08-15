//
//  Copyright (C) 2005-2025 Greg Landrum and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//

#include <catch2/catch_all.hpp>
#include <Geometry/point.h>
#include "FreeChemicalFeature.h"
#include <RDGeneral/utils.h>

using namespace ChemicalFeatures;

TEST_CASE("FreeChemicalFeature Tests") {
  FreeChemicalFeature f1("foo", "bar", RDGeom::Point3D(1, 1, 1));
  REQUIRE(f1.getId() == -1);
  REQUIRE(f1.getFamily() == "foo");
  REQUIRE(f1.getType() == "bar");
  REQUIRE_THAT(f1.getPos().x, Catch::Matchers::WithinAbs(1.0, 1e-6));
  REQUIRE_THAT(f1.getPos().y, Catch::Matchers::WithinAbs(1.0, 1e-6));
  REQUIRE_THAT(f1.getPos().z, Catch::Matchers::WithinAbs(1.0, 1e-6));

  FreeChemicalFeature f2("foo", "bar", RDGeom::Point3D(1, 1, 1), 123);
  REQUIRE(f2.getId() == 123);
  REQUIRE(f2.getFamily() == "foo");
  REQUIRE(f2.getType() == "bar");
  REQUIRE_THAT(f2.getPos().x, Catch::Matchers::WithinAbs(1.0, 1e-6));
  REQUIRE_THAT(f2.getPos().y, Catch::Matchers::WithinAbs(1.0, 1e-6));
  REQUIRE_THAT(f2.getPos().z, Catch::Matchers::WithinAbs(1.0, 1e-6));

  FreeChemicalFeature f3;
  f3.initFromString(f2.toString());
  REQUIRE(f3.getId() == 123);
  REQUIRE(f3.getFamily() == "foo");
  REQUIRE(f3.getType() == "bar");
  REQUIRE_THAT(f3.getPos().x, Catch::Matchers::WithinAbs(1.0, 1e-6));
  REQUIRE_THAT(f3.getPos().y, Catch::Matchers::WithinAbs(1.0, 1e-6));
  REQUIRE_THAT(f3.getPos().z, Catch::Matchers::WithinAbs(1.0, 1e-6));

  FreeChemicalFeature f4(f2);
  REQUIRE(f4.getId() == 123);
  REQUIRE(f4.getFamily() == "foo");
  REQUIRE(f4.getType() == "bar");
  REQUIRE_THAT(f4.getPos().x, Catch::Matchers::WithinAbs(1.0, 1e-6));
  REQUIRE_THAT(f4.getPos().y, Catch::Matchers::WithinAbs(1.0, 1e-6));
  REQUIRE_THAT(f4.getPos().z, Catch::Matchers::WithinAbs(1.0, 1e-6));
}
