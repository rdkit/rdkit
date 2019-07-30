#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do
                           // this in one cpp file
#include "catch.hpp"
#include <Geometry/point.h>

TEST_CASE("construct Point2D from Point3D", "[point]") {
  SECTION("basics") {
    RDGeom::Point3D p3(1., 2., 3.);
    RDGeom::Point2D p2(p3);
    CHECK(p2.x == p3.x);
    CHECK(p2.y == p3.y);
  }
}
