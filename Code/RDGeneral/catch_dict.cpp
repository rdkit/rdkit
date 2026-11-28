//
//  Copyright (C) 2021-2026 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <string>
#include <boost/numeric/conversion/cast.hpp>
#include <catch2/catch_all.hpp>
#include "Dict.h"
#include "RDProps.h"
#include "Exceptions.h"
#include "StreamOps.h"
#include <cstdint>
#include <limits>
#include <sstream>
#include <vector>
using namespace std::string_literals;

TEST_CASE("Dict move semantics") {
  RDKit::Dict d1;
  d1.setVal("foo"s, 1);
  d1.setVal("bar"s, "yep");
  SECTION("move constructor") {
    CHECK(d1.hasVal("foo"s));
    CHECK(d1.hasVal("bar"s));
    CHECK(d1.getVal<std::string>("bar"s) == "yep"s);
    auto d2 = std::move(d1);
    CHECK(d2.hasVal("foo"s));
    CHECK(d2.hasVal("bar"s));
    CHECK(d2.getVal<std::string>("bar"s) == "yep"s);
    CHECK(!d1.hasVal("foo"s));
    CHECK(!d1.hasVal("bar"s));
  }
  SECTION("move assignment") {
    CHECK(d1.hasVal("foo"s));
    CHECK(d1.hasVal("bar"s));
    CHECK(d1.getVal<std::string>("bar"s) == "yep"s);
    RDKit::Dict d2;
    d2 = std::move(d1);
    CHECK(d2.hasVal("foo"s));
    CHECK(d2.hasVal("bar"s));
    CHECK(d2.getVal<std::string>("bar"s) == "yep"s);
    CHECK(!d1.hasVal("foo"s));
    CHECK(!d1.hasVal("bar"s));
  }
}

TEST_CASE("RDProps move semantics") {
  RDKit::RDProps d1;
  d1.setProp("foo"s, 1);
  d1.setProp("bar"s, "yep");
  SECTION("move constructor") {
    CHECK(d1.hasProp("foo"s));
    CHECK(d1.hasProp("bar"s));
    CHECK(d1.getProp<std::string>("bar"s) == "yep"s);
    auto d2 = std::move(d1);
    CHECK(d2.hasProp("foo"s));
    CHECK(d2.hasProp("bar"s));
    CHECK(d2.getProp<std::string>("bar"s) == "yep"s);
    CHECK(!d1.hasProp("foo"s));
    CHECK(!d1.hasProp("bar"s));
  }
  SECTION("move assignment") {
    CHECK(d1.hasProp("foo"s));
    CHECK(d1.hasProp("bar"s));
    CHECK(d1.getProp<std::string>("bar"s) == "yep"s);
    RDKit::RDProps d2;
    d2 = std::move(d1);
    CHECK(d2.hasProp("foo"s));
    CHECK(d2.hasProp("bar"s));
    CHECK(d2.getProp<std::string>("bar"s) == "yep"s);
    CHECK(!d1.hasProp("foo"s));
    CHECK(!d1.hasProp("bar"s));
  }
}
TEST_CASE("github #9068: properties with empty names") {
  RDKit::RDProps props;
  SECTION("setProp with empty key") {
    CHECK_THROWS_AS(props.setProp("", 1), ValueErrorException);
  }
  SECTION("getProp with empty key") {
    CHECK_THROWS_AS(props.getProp<int>(""), KeyErrorException);
  }
  SECTION("hasProp with empty key") { CHECK(!props.hasProp("")); }
  SECTION("clearProp with empty key") { CHECK_NOTHROW(props.clearProp("")); }
}
TEST_CASE("github #9068: dicts with empty keys") {
  RDKit::Dict dict;
  SECTION("setVal with empty key") {
    CHECK_THROWS_AS(dict.setVal("", 1), ValueErrorException);
  }
  SECTION("getVal with empty key") {
    CHECK_THROWS_AS(dict.getVal<int>(""), KeyErrorException);
  }
  SECTION("hasVal with empty key") { CHECK(!dict.hasVal("")); }
  SECTION("clearVal with empty key") { CHECK_NOTHROW(dict.clearVal("")); }
}
TEST_CASE("streamReadStringVec rejects element counts that do not fit in size_t") {
  auto makeStream = [](std::uint64_t count) {
    std::stringstream ss;
    ss.write(reinterpret_cast<const char *>(&count), sizeof(count));
    // a few empty strings (zero length prefix) so a bad count is not caught
    // by an early end of stream
    for (int i = 0; i < 4; ++i) {
      std::uint32_t len = 0;
      ss.write(reinterpret_cast<const char *>(&len), sizeof(len));
    }
    return ss;
  };
  std::vector<std::string> v;
  SECTION("largest possible count") {
    auto ss = makeStream(std::numeric_limits<std::uint64_t>::max());
    CHECK_THROWS(RDKit::streamReadStringVec(ss, v, 0));
  }
  if constexpr (sizeof(size_t) == 4) {
    // this count used to be silently truncated to 1 on 32-bit platforms,
    // after which the read loop ran past the end of the vector
    SECTION("count that wraps around on 32-bit platforms") {
      auto ss = makeStream((std::uint64_t(1) << 32) + 1);
      CHECK_THROWS_AS(RDKit::streamReadStringVec(ss, v, 0),
                      boost::numeric::bad_numeric_cast);
    }
  }
}
