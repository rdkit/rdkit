//
// Copyright (C)  2004-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>
#include <Geometry/point.h>

#include <Features/Feature.h>

using namespace RDKit;
using namespace RDGeom;
using namespace RDFeatures;

typedef enum {
  fooType,
  barType,
  bazType,
  grnType,
} TypeMarker;

TEST_CASE("Basics for ExplicitFeatures") {
  ExplicitFeature<TypeMarker> f1;
  f1.setFamily(bazType);
  REQUIRE(f1.getFamily() == bazType);
  f1.setType(grnType);
  REQUIRE(f1.getType() == grnType);
  REQUIRE_THAT(f1.getLoc().x, Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(f1.getLoc().y, Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(f1.getLoc().z, Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE(f1.getDirs().size() == 0);

  f1 = ExplicitFeature<TypeMarker>(barType, fooType);
  REQUIRE(f1.getFamily() == barType);
  REQUIRE(f1.getType() == fooType);
  REQUIRE_THAT(f1.getLoc().x, Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(f1.getLoc().y, Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(f1.getLoc().z, Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE(f1.getDirs().size() == 0);

  f1 = ExplicitFeature<TypeMarker>(barType, fooType, Point3D(1.0, 2.0, 3.0));
  REQUIRE(f1.getFamily() == barType);
  REQUIRE(f1.getType() == fooType);
  REQUIRE_THAT(f1.getLoc().x, Catch::Matchers::WithinAbs(1.0, 1e-4));
  REQUIRE_THAT(f1.getLoc().y, Catch::Matchers::WithinAbs(2.0, 1e-4));
  REQUIRE_THAT(f1.getLoc().z, Catch::Matchers::WithinAbs(3.0, 1e-4));
  REQUIRE(f1.getDirs().size() == 0);

  f1.setLoc(Point3D(-1.0, -2.0, -3.0));
  REQUIRE_THAT(f1.getLoc().x, Catch::Matchers::WithinAbs(-1.0, 1e-4));
  REQUIRE_THAT(f1.getLoc().y, Catch::Matchers::WithinAbs(-2.0, 1e-4));
  REQUIRE_THAT(f1.getLoc().z, Catch::Matchers::WithinAbs(-3.0, 1e-4));
  REQUIRE(f1.getDirs().size() == 0);
}

TEST_CASE("Basics for ImplicitFeatures") {
  ImplicitFeature<TypeMarker> f1;
  f1.setType(fooType);
  REQUIRE(f1.getType() == fooType);
  f1.setType(grnType);
  REQUIRE(f1.getType() == grnType);
  REQUIRE_THAT(f1.getLoc().x, Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(f1.getLoc().y, Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(f1.getLoc().z, Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE(f1.getDirs().size() == 0);

  f1 = ImplicitFeature<TypeMarker>(barType, fooType);
  REQUIRE(f1.getFamily() == barType);
  REQUIRE(f1.getType() == fooType);
  REQUIRE_THAT(f1.getLoc().x, Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(f1.getLoc().y, Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(f1.getLoc().z, Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE(f1.getDirs().size() == 0);

  Point3D p1(0, 0, 0), p2(1, 0, 0), p3(0, 1, 0);
  f1.addPoint(&p1);
  f1.addPoint(&p2);
  REQUIRE_THAT(f1.getLoc().x, Catch::Matchers::WithinAbs(0.50, 1e-4));
  REQUIRE_THAT(f1.getLoc().y, Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(f1.getLoc().z, Catch::Matchers::WithinAbs(0.0, 1e-4));

  f1.addPoint(&p3);
  REQUIRE_THAT(f1.getLoc().x, Catch::Matchers::WithinAbs(0.3333, 1e-4));
  REQUIRE_THAT(f1.getLoc().y, Catch::Matchers::WithinAbs(0.3333, 1e-4));
  REQUIRE_THAT(f1.getLoc().z, Catch::Matchers::WithinAbs(0.0, 1e-4));

  f1.reset();
  REQUIRE_THAT(f1.getLoc().x, Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(f1.getLoc().y, Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(f1.getLoc().z, Catch::Matchers::WithinAbs(0.0, 1e-4));
}

TEST_CASE("ExplicitFeatures 2D, string type.") {
  typedef ExplicitFeature<std::string, std::string, Point2D> LocalFeature;
  LocalFeature f1;
  f1.setType("foo");
  REQUIRE(f1.getType() == "foo");
  f1.setFamily("foob");
  REQUIRE(f1.getFamily() == "foob");
  REQUIRE_THAT(f1.getLoc().x, Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(f1.getLoc().y, Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE(f1.getDirs().size() == 0);

  f1 = LocalFeature("foo", "bar");
  REQUIRE(f1.getFamily() == "foo");
  REQUIRE(f1.getType() == "bar");
  REQUIRE_THAT(f1.getLoc().x, Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE_THAT(f1.getLoc().y, Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE(f1.getDirs().size() == 0);

  f1 = LocalFeature("grm", "grn", Point2D(1.0, 2.0));
  REQUIRE(f1.getFamily() == "grm");
  REQUIRE(f1.getType() == "grn");
  REQUIRE_THAT(f1.getLoc().x, Catch::Matchers::WithinAbs(1.0, 1e-4));
  REQUIRE_THAT(f1.getLoc().y, Catch::Matchers::WithinAbs(2.0, 1e-4));
  REQUIRE(f1.getDirs().size() == 0);

  f1.setLoc(Point2D(-1.0, -2.0));
  REQUIRE_THAT(f1.getLoc().x, Catch::Matchers::WithinAbs(-1.0, 1e-4));
  REQUIRE_THAT(f1.getLoc().y, Catch::Matchers::WithinAbs(-2.0, 1e-4));
  REQUIRE(f1.getDirs().size() == 0);
}
