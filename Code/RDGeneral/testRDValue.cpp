#include <catch2/catch_all.hpp>
#include "RDValue.h"
#include "RDProps.h"
#include "Invariant.h"
#include "StreamOps.h"
#include <limits>
#include <vector>
#include <list>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <any>

using namespace RDKit;

template <class T>
void testLimits() {
  // check numeric limits
  {
    RDValue v(std::numeric_limits<T>::min());
    std::cerr << "min: " << std::numeric_limits<T>::min() << " "
              << rdvalue_cast<T>(v) << std::endl;
    REQUIRE(rdvalue_cast<T>(v) == std::numeric_limits<T>::min());
    REQUIRE(rdvalue_cast<T>(RDValue(v)) == std::numeric_limits<T>::min());
    v = std::numeric_limits<T>::max();
    REQUIRE(rdvalue_cast<T>(v) == std::numeric_limits<T>::max());
    REQUIRE(rdvalue_cast<T>(RDValue(v)) == std::numeric_limits<T>::max());
  }
  {
    RDValue v(std::numeric_limits<T>::max());
    REQUIRE(rdvalue_cast<T>(v) == std::numeric_limits<T>::max());
    RDValue vv(v);
    REQUIRE(rdvalue_cast<T>(vv) == std::numeric_limits<T>::max());

    v = std::numeric_limits<T>::min();
    RDValue vvv(v);
    REQUIRE(rdvalue_cast<T>(v) == std::numeric_limits<T>::min());
    REQUIRE(rdvalue_cast<T>(vvv) == std::numeric_limits<T>::min());
  }
}

TEST_CASE("testLimits") {
  {
    RDValue v(std::numeric_limits<int>::min());
    REQUIRE(rdvalue_cast<int>(v) == std::numeric_limits<int>::min());
    REQUIRE(rdvalue_cast<int>(RDValue(v)) == std::numeric_limits<int>::min());
    v = std::numeric_limits<int>::max();
    REQUIRE(rdvalue_cast<int>(v) == std::numeric_limits<int>::max());
    REQUIRE(rdvalue_cast<int>(RDValue(v)) == std::numeric_limits<int>::max());
  }
  {
    RDValue v(std::numeric_limits<int>::max());
    REQUIRE(rdvalue_cast<int>(v) == std::numeric_limits<int>::max());
    RDValue vv(v);
    REQUIRE(rdvalue_cast<int>(vv) == std::numeric_limits<int>::max());

    v = std::numeric_limits<int>::min();
    RDValue vvv(v);
    REQUIRE(rdvalue_cast<int>(v) == std::numeric_limits<int>::min());
    REQUIRE(rdvalue_cast<int>(vvv) == std::numeric_limits<int>::min());
  }
}

TEST_CASE("testPOD") {
  testLimits<int>();
  testLimits<unsigned int>();
  testLimits<double>();
  testLimits<float>();
  testLimits<bool>();
}

template <class T>
void testVector() {
  T minv = std::numeric_limits<T>::min();
  T maxv = std::numeric_limits<T>::max();
  std::vector<T> data;
  data.push_back(minv);
  data.push_back(maxv);
  data.push_back(T());

  RDValue v(data);
  REQUIRE(rdvalue_cast<std::vector<T>>(v) == data);
  RDValue vv;
  copy_rdvalue(vv, v);
  REQUIRE(rdvalue_cast<std::vector<T>>(vv) == data);
  RDValue::cleanup_rdvalue(v);  // desctructor...
  RDValue::cleanup_rdvalue(vv);
}

TEST_CASE("testPODVectors") {
  testVector<int>();
  testVector<unsigned int>();
  testVector<double>();
  testVector<float>();
  testVector<long double>();  // stored in anys
}

TEST_CASE("testStringVect") {
  std::vector<std::string> vecs;
  vecs.emplace_back("my");
  vecs.emplace_back("dog");
  vecs.emplace_back("has");
  vecs.emplace_back("fleas");
  RDValue v(vecs);
  REQUIRE(rdvalue_cast<std::vector<std::string>>(v) == vecs);
  RDValue vc;
  copy_rdvalue(vc, v);
  REQUIRE(rdvalue_cast<std::vector<std::string>>(vc) == vecs);
  RDValue vv = vecs;
  RDValue vvc;
  copy_rdvalue(vvc, vv);
  REQUIRE(rdvalue_cast<std::vector<std::string>>(vv) == vecs);
  REQUIRE(rdvalue_cast<std::vector<std::string>>(vvc) == vecs);

  RDValue::cleanup_rdvalue(v);
  RDValue::cleanup_rdvalue(vc);
  RDValue::cleanup_rdvalue(vv);
  RDValue::cleanup_rdvalue(vvc);
}

TEST_CASE("testMapsAndLists") {
  {
    typedef std::map<std::string, int> listtype;
    listtype m;
    m["foo"] = 1;
    m["bar"] = 2;
    RDValue v(m);
    REQUIRE(rdvalue_cast<listtype>(v) == m);
    RDValue::cleanup_rdvalue(v);
  }
  {
    std::list<std::string> m;
    m.emplace_back("foo");
    m.emplace_back("bar");
    RDValue v(m);
    REQUIRE(rdvalue_cast<std::list<std::string>>(v) == m);
    RDValue::cleanup_rdvalue(v);
  }
}

TEST_CASE("testNaN") {
  double nan = sqrt(-1.0);
  RDValue v(nan);
  REQUIRE(v.getTag() == RDTypeTag::DoubleTag);
  REQUIRE(std::isnan(rdvalue_cast<double>(v)));

  RDValue vv(2.0);
  REQUIRE(rdvalue_is<double>(vv));
  REQUIRE(vv.getTag() == RDTypeTag::DoubleTag);
}

template <class T>
std::vector<T> makeVec() {
  std::vector<T> vec;
  vec.push_back((T)0);
  vec.push_back((T)1);
  vec.push_back((T)2);
  vec.push_back((T)3);
  vec.push_back((T)4);
  vec.push_back((T)5);
  return vec;
}

template <class T>
void testProp(T val) {
  std::stringstream ss;

  {
    RDProps p;
    p.setProp<T>("foo", val);
    TEST_ASSERT(streamWriteProps(ss, p));
  }

  {
    RDProps p2;
    streamReadProps(ss, p2);
    TEST_ASSERT(p2.getProp<T>("foo") == val);
  }
};

TEST_CASE("testPropertyPickler") {
  testProp<int>(1234);
  testProp<double>(1234.);
  testProp<float>(1234.0f);
  testProp<unsigned int>(1234u);
  testProp<bool>(true);
  testProp<std::string>(
      std::string("the quick brown fox jumps over the lazy dog"));

  testProp(0);
  testProp(0.);
  testProp(0.0f);
  testProp(0u);
  testProp(false);
}

TEST_CASE("testPickleBinaryString") {
  char buf[10];
  for (int i = 0; i < 10; ++i) {
    buf[i] = (char)i;
  }
  std::string str(buf, 10);
  std::stringstream ss;

  {
    RDProps p;
    p.setProp<std::string>("foo", str);
    REQUIRE(streamWriteProps(ss, p));
  }

  {
    RDProps p2;
    streamReadProps(ss, p2);
    REQUIRE(p2.getProp<std::string>("foo") == str);
  }
}

TEST_CASE("testIntConversions") {
  RDProps p;
  p.setProp<int>("foo", 1);
  p.getProp<std::int64_t>("foo");
  p.getProp<std::int8_t>("foo");
  p.getProp<std::int16_t>("foo");
  p.getProp<std::uint16_t>("foo");

  p.setProp<int64_t>("foo", 1);
  p.getProp<int64_t>("foo");

  p.setProp<unsigned int>("foo", 1);
  p.getProp<std::uint64_t>("foo");
  p.getProp<std::uint8_t>("foo");
  p.getProp<std::uint16_t>("foo");

  p.getProp<std::int16_t>("foo");

  p.setProp<unsigned int>("foo", 0);
  p.getProp<std::uint8_t>("foo");
  p.getProp<std::uint16_t>("foo");

  p.setProp<unsigned int>("foo", 255);
  p.getProp<std::uint8_t>("foo");

  p.setProp<unsigned int>("foo", 65535);
  p.getProp<std::uint16_t>("foo");

  p.setProp<int>("foo", -128);
  p.getProp<std::int8_t>("foo");

  p.setProp<int>("foo", -32768);
  p.getProp<std::int16_t>("foo");

  p.setProp<int>("foo", 127);
  p.getProp<std::int8_t>("foo");

  p.setProp<int>("foo", 32767);
  p.getProp<std::int16_t>("foo");

  p.setProp<int>("foo", 32767 + 1);
  REQUIRE_THROWS_AS(p.getProp<std::int8_t>("foo"),
                    boost::numeric::positive_overflow);
  REQUIRE_THROWS_AS(p.getProp<std::uint8_t>("foo"),
                    boost::numeric::positive_overflow);
  REQUIRE_THROWS_AS(p.getProp<std::int16_t>("foo"),
                    boost::numeric::positive_overflow);
  p.setProp<int>("foo", 65535 + 1);
  REQUIRE_THROWS_AS(p.getProp<std::uint16_t>("foo"),
                    boost::numeric::positive_overflow);

  p.setProp<int>("foo", -1);
  REQUIRE_THROWS_AS(p.getProp<std::uint8_t>("foo"),
                    boost::numeric::negative_overflow);

  p.getProp<std::int16_t>("foo");
  REQUIRE_THROWS_AS(p.getProp<std::uint16_t>("foo"),
                    boost::numeric::negative_overflow);
}

TEST_CASE("testStringToDouble") {
  RDProps p;
  p.setProp<std::string>("foo", "123.0 ");
  p.setProp<std::string>("bar", " 123.0 ");
  REQUIRE(p.getProp<double>("foo") == 123.0);
  REQUIRE(p.getProp<float>("foo") == 123.0f);
  REQUIRE_THROWS_AS(p.getProp<double>("bar"),
		    std::bad_any_cast);
}
