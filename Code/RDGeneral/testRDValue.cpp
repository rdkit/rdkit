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
#include <boost/math/special_functions/fpclassify.hpp>

using namespace RDKit;

template<class T>
void testLimits() {
  // check numeric limits
  {
    RDValue v(std::numeric_limits<T>::min());
    std::cerr << "min: " << std::numeric_limits<T>::min() << " " << rdvalue_cast<T>(v) <<
        std::endl;
    CHECK_INVARIANT(rdvalue_cast<T>(v) == std::numeric_limits<T>::min(), "bad min");
    CHECK_INVARIANT(rdvalue_cast<T>(RDValue(v)) == std::numeric_limits<T>::min(), "bad min");
    v = std::numeric_limits<T>::max();
    CHECK_INVARIANT(rdvalue_cast<T>(v) == std::numeric_limits<T>::max(), "bad max");
    CHECK_INVARIANT(rdvalue_cast<T>(RDValue(v)) == std::numeric_limits<T>::max(), "bad max");
    
  }
  {
    RDValue v(std::numeric_limits<T>::max());
    CHECK_INVARIANT(rdvalue_cast<T>(v) == std::numeric_limits<T>::max(), "bad max");
    RDValue vv(v);
    CHECK_INVARIANT(rdvalue_cast<T>(vv) == std::numeric_limits<T>::max(), "bad max");

    v = std::numeric_limits<T>::min();
    RDValue vvv(v);
    CHECK_INVARIANT(rdvalue_cast<T>(v) == std::numeric_limits<T>::min(), "bad min");
    CHECK_INVARIANT(rdvalue_cast<T>(vvv) == std::numeric_limits<T>::min(), "bad min");
  }
}

void testPOD() {
  testLimits<int>();
  testLimits<unsigned int>();
  testLimits<double>();
  testLimits<float>();
  testLimits<bool>();
}



template<class T>
void testVector() {
  T minv = std::numeric_limits<T>::min();
  T maxv = std::numeric_limits<T>::max();
  std::vector<T> data;
  data.push_back(minv);
  data.push_back(maxv);
  data.push_back(T());

  RDValue v(data);
  CHECK_INVARIANT(rdvalue_cast<std::vector<T> >(v) == data, "bad vec");
  RDValue vv; copy_rdvalue(vv,v);
  CHECK_INVARIANT(rdvalue_cast<std::vector<T> >(vv) == data, "bad copy constructor");
  RDValue::cleanup_rdvalue(v); // desctructor...
  RDValue::cleanup_rdvalue(vv);
}

void testPODVectors() {
  testVector<int>();
  testVector<unsigned int>();
  testVector<double>();
  testVector<float>();
  testVector<long double>(); // stored in anys
}

void testStringVect() {
  std::vector<std::string> vecs;
  vecs.push_back("my");
  vecs.push_back("dog");
  vecs.push_back("has");
  vecs.push_back("fleas");
  RDValue v(vecs);
  CHECK_INVARIANT(rdvalue_cast<std::vector<std::string> >(v) == vecs, "bad vect");
  RDValue vc; copy_rdvalue(vc, v);
  CHECK_INVARIANT(rdvalue_cast<std::vector<std::string> >(vc) == vecs, "bad vect");
  RDValue vv = vecs;
  RDValue vvc; copy_rdvalue(vvc, vv);
  CHECK_INVARIANT(rdvalue_cast<std::vector<std::string> >(vv) == vecs, "bad vect");          
  CHECK_INVARIANT(rdvalue_cast<std::vector<std::string> >(vvc) == vecs, "bad vect");

  RDValue::cleanup_rdvalue(v); // desctructor...
  RDValue::cleanup_rdvalue(vc); // desctructor...
  RDValue::cleanup_rdvalue(vv); // desctructor...
  RDValue::cleanup_rdvalue(vvc); // desctructor...
}

void testMapsAndLists() {
  {
    typedef std::map<std::string, int> listtype;
    listtype m;
    m["foo"] = 1;
    m["bar"] = 2;
    RDValue v(m);
    CHECK_INVARIANT(rdvalue_cast<listtype>(v) == m, "bad map cast");
    RDValue::cleanup_rdvalue(v);
  }
  {
    std::list<std::string> m;
    m.push_back("foo");
    m.push_back("bar");
    RDValue v(m);
    CHECK_INVARIANT(rdvalue_cast<std::list<std::string> >(v) == m, "bad map cast");
    RDValue::cleanup_rdvalue(v);
  }
}

void testNaN() {
  // make a NaN
  double nan=sqrt(-1.0);
  RDValue v(nan);
  TEST_ASSERT(v.getTag() == RDTypeTag::DoubleTag);

  CHECK_INVARIANT(boost::math::isnan(rdvalue_cast<double>(v)), "Oops, can't store NaNs!");

  RDValue vv(2.0);
  TEST_ASSERT(rdvalue_is<double>(vv));
  TEST_ASSERT(vv.getTag() == RDTypeTag::DoubleTag);
}

template<class T>
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

template<class T>
void testProp(T val) {
  std::string pklName = getenv("RDBASE");
  pklName += "/Code/GraphMol/test_data/prop.pkl";
  
  {
    std::ofstream ss(pklName.c_str(), std::ios_base::binary);
    RDProps p;
    p.setProp<T>("foo", val);
    TEST_ASSERT(streamWriteProps(ss, p));
  }
  
  {
    std::ifstream ss(pklName.c_str(), std::ios_base::binary);
    RDProps p2;
    streamReadProps(ss, p2);
    TEST_ASSERT(p2.getProp<T>("foo") == val);
  }
};

void testPropertyPickler() {
  std::cerr << "== int" << std::endl;
  testProp<int>(1234);
  std::cerr << "== double" << std::endl;
  testProp<double>(1234.);
  std::cerr << "== float" << std::endl;
  testProp<float>(1234.0f);
  std::cerr << "== unsigned int" << std::endl;
  testProp<unsigned int>(1234u);
  std::cerr << "== bool" << std::endl;
  testProp<bool>(true);
  std::cerr << "== std::string" << std::endl;
  testProp<std::string>(std::string("the quick brown fox jumps over the lazy dog"));
  
  testProp(0);
  testProp(0.);
  testProp(0.0f);
  testProp(0u);
  testProp(false);

  /*
  testProp(makeVec<int>());
  testProp(makeVec<int>());
  testProp(makeVec<int>());
  testProp(makeVec<unsigned int>());

  {
    std::vector<std::string> v;
    v.push_back("a");
    v.push_back("b");
    v.push_back("c");
    v.push_back("d");
    v.push_back("e");
    testProp(v);
  }
  */
}


int main() {
  std::cerr << "-- running tests -- " << std::endl;
  testPOD();
  testPODVectors();
  testStringVect();
  testNaN();
  testPropertyPickler();
}
