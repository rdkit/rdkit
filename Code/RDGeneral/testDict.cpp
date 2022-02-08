//
// Copyright 2001-2008 Randal M. Henne, Greg Landrum and
//                     Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//

#include <RDGeneral/test.h>
#include "types.h"
#include "StreamOps.h"
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDAny.h>
#include <RDGeneral/Dict.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>
#include <boost/shared_ptr.hpp>
#include <vector>
#include <ctime>

using namespace RDKit;
using namespace std;

struct Foo {
  int bar{0};
  float baz{0.f};
  Foo() {}
  Foo(int bar, float baz) : bar(bar), baz(baz) {}
  ~Foo() { std::cerr << "deleted!" << std::endl; }
};

void testGithub940() {
  BOOST_LOG(rdErrorLog)
      << "Testing Github940: property dictionaries leaking memory" << std::endl;

  // a couple small tests to check for memory leaks. Only useful with valgrind
  {  // tests computed props
    STR_VECT computed;
    auto *d = new Dict();
    d->setVal(RDKit::detail::computedPropName, computed);
    computed.push_back("foo");
    d->setVal(RDKit::detail::computedPropName, computed);
    delete d;
  }
  {  // tests computed props
    STR_VECT computed;
    auto *d = new Dict();
    d->setVal(RDKit::detail::computedPropName, computed);
    computed.push_back("foo");
    d->setVal(RDKit::detail::computedPropName, computed);
    d->clearVal(RDKit::detail::computedPropName);
    delete d;
  }
  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void testRDAny() {
  std::cerr << "Testing RDValue" << std::endl;
  { RDAny v(-2147450880); }

  {
    int vi = 0;
    RDValue v(0);
    for (int i = 0; i < 100; ++i) {
      vi += i;
      v = rdvalue_cast<int>(v) + i;
      TEST_ASSERT(vi == rdvalue_cast<int>(v));
    }
  }
  std::cerr << "Testing RDAny" << std::endl;
  {
    RDAny a(1);
    RDAny b = a;
    TEST_ASSERT(rdany_cast<int>(a) == 1);
    TEST_ASSERT(rdany_cast<int>(b) == 1);
  }

  {
    RDAny a(1);
    RDAny b = a;
    TEST_ASSERT(rdany_cast<int>(a) == 1);
    TEST_ASSERT(rdany_cast<int>(b) == rdany_cast<int>(a));
    std::map<std::string, RDAny> foo;
    foo["foo"] = a;
    foo["bar"] = std::string("This is a test");
    TEST_ASSERT(rdany_cast<int>(foo["foo"]) == 1);
    TEST_ASSERT(rdany_cast<int>(foo["foo"]) == rdany_cast<int>(a));
    TEST_ASSERT(rdany_cast<std::string>(foo["bar"]) == "This is a test");
  }

  {
    bool a = true;
    RDValue v(a);
    TEST_ASSERT(rdvalue_cast<bool>(v) == true);
    v = (int)10;
    TEST_ASSERT(rdvalue_cast<int>(v) == 10);
  }

  {
    Dict d;
    bool a = true;
    d.setVal("foo", a);
    d.getVal<bool>("foo");
  }

  {  // tests computed props
    STR_VECT computed;
    Dict d;
    d.setVal(RDKit::detail::computedPropName, computed);
    computed.push_back("foo");
    d.setVal(RDKit::detail::computedPropName, computed);
    STR_VECT computed2 = d.getVal<STR_VECT>(RDKit::detail::computedPropName);
    TEST_ASSERT(computed2[0] == "foo");
    Dict d2(d);
    computed2 = d2.getVal<STR_VECT>(RDKit::detail::computedPropName);
    TEST_ASSERT(computed2[0] == "foo");
  }

  {
    Dict d;
    // int v=1;
    // d.setVal("foo", v);
    // TEST_ASSERT(d.getVal<int>("foo") == 1, "bad getval");

    std::vector<int> fooV;
    fooV.resize(3);
    fooV[0] = 1;
    fooV[1] = 2;
    fooV[2] = 3;
    if (0) {
      std::vector<int> fooV2;
      std::cerr << "send int vect" << std::endl;
      RDAny a(fooV);
      std::cerr << "retrieve int vect" << std::endl;
      fooV2 = rdany_cast<std::vector<int>>(a);
      TEST_ASSERT(fooV == fooV2);
    }

    {
      std::vector<int> fooV2;
      std::cerr << "dict set int vect" << std::endl;
      d.setVal("bar", fooV);
      std::cerr << "dict get int vect" << std::endl;
      d.getVal("bar", fooV2);
      TEST_ASSERT(fooV == fooV2);
    }
  }

  {
    std::vector<int> v;
    for (int i = 0; i < 4; ++i) {
      v.push_back(i);
    }

    RDAny foo(v);
    RDAny bar = foo;
    RDAny baz(foo);

    for (int i = 0; i < 4; ++i) {
      TEST_ASSERT(rdany_cast<std::vector<int>>(foo)[i] == i);
      TEST_ASSERT(rdany_cast<std::vector<int>>(bar)[i] == i);
      TEST_ASSERT(rdany_cast<std::vector<int>>(baz)[i] == i);
    }
  }

  {
    std::vector<double> v;
    for (double i = 0; i < 4; ++i) {
      v.push_back(i);
    }

    RDAny foo(v);

    for (int i = 0; i < 4; ++i) {
      TEST_ASSERT(rdany_cast<std::vector<double>>(foo)[i] == i);
    }

    RDAny b = foo;

    for (int i = 0; i < 4; ++i) {
      TEST_ASSERT(rdany_cast<std::vector<double>>(b)[i] == i);
    }
  }
  const int loops = 10000000;
  {
    std::clock_t clock1 = std::clock();
    boost::any v;
    for (int i = 0; i < loops; ++i) {
      v = i;
    }
    std::clock_t clock2 = std::clock();

    std::cout << "static boost any:"
              << (double)(clock2 - clock1) / CLOCKS_PER_SEC << " s"
              << std::endl;
  }
  {
    std::clock_t clock1 = std::clock();
    boost::any *v = nullptr, *vv;
    for (int i = 0; i < loops; ++i) {
      vv = new boost::any(v ? boost::any_cast<int>(*v) + i : i);
      delete v;
      v = vv;
    }
    delete vv;
    std::clock_t clock2 = std::clock();

    std::cout << "dynamic boost any:"
              << (double)(clock2 - clock1) / CLOCKS_PER_SEC << " s"
              << std::endl;
  }

  {
    std::clock_t clock1 = std::clock();
    RDAny v;
    for (int i = 0; i < loops; ++i) {
      v = i;
    }
    std::clock_t clock2 = std::clock();

    std::cout << "static RDAny:" << (double)(clock2 - clock1) / CLOCKS_PER_SEC
              << " s" << std::endl;
  }

  {
    std::clock_t clock1 = std::clock();
    RDAny *v = nullptr, *vv;
    for (int i = 0; i < loops; ++i) {
      vv = new RDAny(v ? rdany_cast<int>(*v) + i : i);
      delete v;
      v = vv;
    }
    delete vv;
    std::clock_t clock2 = std::clock();

    std::cout << "dynamic RDAny:" << (double)(clock2 - clock1) / CLOCKS_PER_SEC
              << " s" << std::endl;
  }

  {
    std::clock_t clock1 = std::clock();
    RDValue v;
    for (int i = 0; i < loops; ++i) {
      v = i;
    }
    std::clock_t clock2 = std::clock();

    std::cout << "static RDValue:" << (double)(clock2 - clock1) / CLOCKS_PER_SEC
              << " s" << std::endl;
  }

  {
    std::clock_t clock1 = std::clock();
    RDValue v(0);
    for (int i = 0; i < loops; ++i) {
      v = RDValue(rdvalue_cast<int>(v) + i);
    }

    std::clock_t clock2 = std::clock();

    std::cout << "dynamic RDValue:"
              << (double)(clock2 - clock1) / CLOCKS_PER_SEC << " s"
              << std::endl;
  }

  {  // checks replacement with vector
    RDAny vv(2.0);
    TEST_ASSERT(rdany_cast<double>(vv) == 2.0);

    std::vector<int> vect;
    vect.push_back(1);
    vv = vect;
    TEST_ASSERT(rdany_cast<std::vector<int>>(vv)[0] == 1);

    // tests copy
    RDAny vvv(vv);

    TEST_ASSERT(rdany_cast<std::vector<int>>(vvv)[0] == 1);
  }

  {
    // Checks fallback to Any
    std::vector<std::pair<int, int>> pvect;
    pvect.push_back(std::make_pair<int, int>(2, 2));
    boost::any any1(pvect);
    boost::any_cast<std::vector<std::pair<int, int>>>(any1);
    boost::any_cast<std::vector<std::pair<int, int>> &>(any1);
    boost::any_cast<const std::vector<std::pair<int, int>> &>(any1);

    RDAny vv(pvect);
    auto &any = rdany_cast<boost::any &>(vv);
    boost::any_cast<std::vector<std::pair<int, int>>>(any);
    boost::any_cast<std::vector<std::pair<int, int>> &>(any);
    boost::any_cast<const std::vector<std::pair<int, int>> &>(any);

    const std::vector<std::pair<int, int>> &pv =
        rdany_cast<std::vector<std::pair<int, int>>>(vv);
    TEST_ASSERT(pv[0].first == 2);
    RDAny vvv(vv);
    TEST_ASSERT(
        (rdany_cast<std::vector<std::pair<int, int>>>(vvv)[0].first == 2));
  }

  {
    // Check pointers -- RDAny doesn't delete these, must do them manually
    auto *p = new std::vector<int>();
    p->push_back(100);
    RDAny v(p);
    RDAny vv(v);
    try {
      rdany_cast<std::vector<int>>(v);
#ifndef UNSAFE_RDVALUE
      PRECONDITION(0, "Should throw bad cast");
#endif
    } catch (boost::bad_any_cast &) {
    }

    TEST_ASSERT((*rdany_cast<std::vector<int> *>(vv))[0] == 100);
    TEST_ASSERT((*rdany_cast<std::vector<int> *>((const RDAny &)vv))[0] == 100);
    delete p;

    auto *m = new std::map<int, int>();
    (*m)[0] = 1;
    RDAny mv(m);
    // leaks
    auto *anym = rdany_cast<std::map<int, int> *>(mv);
    TEST_ASSERT(anym->find(0) != anym->end());
    delete anym;
  }

  {
    // check shared ptrs -- boost::any deletes these :)
    typedef boost::shared_ptr<std::vector<int>> vptr;
    vptr p(new std::vector<int>());
    p->push_back(100);
    RDAny v(p);
    RDAny vv(v);
    TEST_ASSERT((*rdany_cast<vptr>(v))[0] == 100);
    TEST_ASSERT((*rdany_cast<vptr>(vv))[0] == 100);
    TEST_ASSERT((*rdany_cast<vptr>((const RDAny &)vv))[0] == 100);

    typedef boost::shared_ptr<std::map<int, int>> mptr;
    mptr m(new std::map<int, int>());
    (*m)[0] = 1;
    RDAny mv(m);
    // leaks
    mptr anym = rdany_cast<mptr>(mv);
    TEST_ASSERT(anym->find(0) != anym->end());

    RDAny any3(boost::shared_ptr<Foo>(new Foo(1, 2.f)));
    TEST_ASSERT(any3.m_value.getTag() == RDTypeTag::AnyTag);
  }
}

class DictCon {
 public:
  DictCon() { d.reset(); };
  DictCon(const DictCon &other) { d = other.d; };
  DictCon &operator=(const DictCon &other) {
    d = other.d;
    return *this;
  };
  Dict *getDict() { return &d; };

 private:
  Dict d;
};

void testStringVals() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing String Pickle Roundtrips." << std::endl;
  {
    Dict d;
    std::string sv;
    sv = "1";
    d.setVal("foo", sv);
    int iv;
    d.getVal("foo", iv);
    TEST_ASSERT(iv == 1);
  }
  {
    Dict d;
    d.setVal("foo", "1");
    int iv;
    d.getVal("foo", iv);
    TEST_ASSERT(iv == 1);
  }
  {
    Dict d;
    std::string sv;
    sv = "1.3";
    d.setVal("foo", sv);
    double dv;
    d.getVal("foo", dv);
    TEST_ASSERT(feq(dv, 1.3));
  }

  {
    Dict d;
    int iv = 1;
    d.setVal("foo", iv);
    std::string sv;
    d.getVal("foo", sv);
    TEST_ASSERT(sv == "1");
    sv = d.getVal<std::string>("foo");
    TEST_ASSERT(sv == "1");
  }

  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void testVectToString() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing conversion of vect to string." << std::endl;
  {
    Dict d;
    std::vector<int> v;
    v.push_back(1);
    v.push_back(0);
    d.setVal("foo", v);
    std::string sv;
    d.getVal("foo", sv);
    TEST_ASSERT(sv == "[1,0,]");
  }
  {
    Dict d;
    std::vector<unsigned int> v;
    v.push_back(1);
    v.push_back(0);
    d.setVal("foo", v);
    std::string sv;
    d.getVal("foo", sv);
    TEST_ASSERT(sv == "[1,0,]");
    sv = d.getVal<std::string>("foo");
    TEST_ASSERT(sv == "[1,0,]");
  }
  {
    Dict d;
    std::vector<double> v;
    v.push_back(1.2);
    v.push_back(0);
    d.setVal("foo", v);
    std::string sv;
    d.getVal("foo", sv);
    TEST_ASSERT(sv == "[1.2,0,]");
    sv = d.getVal<std::string>("foo");
    TEST_ASSERT(sv == "[1.2,0,]");
  }
  {
    Dict d;
    std::vector<float> v;
    v.push_back(10001.f);
    v.push_back(0);
    d.setVal("foo", v);
    std::string sv;
    d.getVal("foo", sv);
    TEST_ASSERT(sv == "[10001,0,]");
    sv = d.getVal<std::string>("foo");
    TEST_ASSERT(sv == "[10001,0,]");
  }

  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void testConstReturns() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing returning const references." << std::endl;
  {
    std::string v = "foo";
    RDAny anyv(v);

    std::string tgt = rdany_cast<std::string>(anyv);
    const std::string &ctgt = rdany_cast<std::string>(anyv);
    TEST_ASSERT(ctgt != "");
  }

  {
    Dict d;
    std::string v = "foo";
    d.setVal("foo", v);

    // const std::string nv=d.getVal<const std::string &>("foo");
    std::string nv = d.getVal<std::string>("foo");
    TEST_ASSERT(nv == "foo");
  }

#if 0
  {
    Dict d;
    std::string v="foo";
    d.setVal("foo",v);

    double ls=0;
    BOOST_LOG(rdErrorLog) << "copy" << std::endl;
    for(int i=0;i<100000000;++i){
      std::string nv=d.getVal<std::string>("foo");
      ls+= nv.size();
    }
    BOOST_LOG(rdErrorLog) << "done: "<<ls << std::endl;
    ls=0;
    BOOST_LOG(rdErrorLog) << "ref" << std::endl;
    for(int i=0;i<100000000;++i){
      const std::string &nv=d.getVal<std::string>("foo");
      ls+= nv.size();
    }
    BOOST_LOG(rdErrorLog) << "done: "<<ls << std::endl;
    //std::string nv=d.getVal<std::string>("foo");
  }
#else
  {
    // int nreps=100000000;
    int nreps = 100000;
    Dict d;
    std::string v = "foo";
    RDAny anyv(v);
    d.setVal("foo", v);

    std::clock_t start, end;

    double ls = 0;
    BOOST_LOG(rdErrorLog) << "any cast" << std::endl;
    start = std::clock();
    for (int i = 0; i < nreps; ++i) {
      const std::string &nv = rdany_cast<std::string>(anyv);
      ls += nv.size();
    }
    end = std::clock();
    BOOST_LOG(rdErrorLog) << "done: "
                          << (end - start) / (double)(CLOCKS_PER_SEC) << " "
                          << ls << std::endl;

    ls = 0;
    BOOST_LOG(rdErrorLog) << "copy" << std::endl;
    start = std::clock();
    for (int i = 0; i < nreps; ++i) {
      std::string nv = rdany_cast<std::string>(anyv);
      ls += nv.size();
    }
    end = std::clock();
    BOOST_LOG(rdErrorLog) << "done: "
                          << (end - start) / (double)(CLOCKS_PER_SEC) << " "
                          << ls << std::endl;

    ls = 0;
    BOOST_LOG(rdErrorLog) << "ref" << std::endl;
    start = std::clock();
    for (int i = 0; i < nreps; ++i) {
      const std::string &nv = rdany_cast<std::string>(anyv);
      ls += nv.size();
    }
    end = std::clock();
    BOOST_LOG(rdErrorLog) << "done: "
                          << (end - start) / (double)(CLOCKS_PER_SEC) << " "
                          << ls << std::endl;

    ls = 0;
    BOOST_LOG(rdErrorLog) << "dict" << std::endl;
    start = std::clock();
    for (int i = 0; i < nreps; ++i) {
      const std::string &nv = d.getVal<std::string>("foo");
      ls += nv.size();
    }
    end = std::clock();
    BOOST_LOG(rdErrorLog) << "done: "
                          << (end - start) / (double)(CLOCKS_PER_SEC) << " "
                          << ls << std::endl;

    ls = 0;
    BOOST_LOG(rdErrorLog) << "ref with hasVal" << std::endl;
    start = std::clock();
    std::string k = "foo";
    for (int i = 0; i < nreps; ++i) {
      if (d.hasVal(k)) {
        const std::string &nv = rdany_cast<std::string>(anyv);
        ls += nv.size();
      }
    }
    end = std::clock();
    BOOST_LOG(rdErrorLog) << "done: "
                          << (end - start) / (double)(CLOCKS_PER_SEC) << " "
                          << ls << std::endl;

    // std::string nv=d.getVal<std::string>("foo");
  }

#endif

  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void testUpdate() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing dict update." << std::endl;

  {
    Dict d;
    std::string sv;
    sv = "1.3";
    d.setVal("foo", sv);
    double dv = 3.0;
    d.setVal("foo2", dv);
    std::vector<int> f;
    f.push_back(1);
    f.push_back(2);
    d.setVal("foo3", f);

    Dict d2;
    d2.update(d);
    TEST_ASSERT(d.getVal<std::string>("foo") == d2.getVal<std::string>("foo"));
    TEST_ASSERT(d.getVal<double>("foo2") == d2.getVal<double>("foo2"));
    TEST_ASSERT(d.getVal<std::vector<int>>("foo3") ==
                d2.getVal<std::vector<int>>("foo3"));
  }

  {  // a few tests to make sure copying/updating with nonPOD data is ok
    Dict d;
    std::string sv;
    sv = "1.3";
    d.setVal("foo", sv);
    double dv = 3.0;
    d.setVal("foo2", dv);
    std::vector<int> f;
    f.push_back(1);
    f.push_back(2);
    d.setVal("foo3", f);
    {
      Dict d2;
      d2.setVal("foo", 1);
      d2.update(d, true);
      TEST_ASSERT(d2.getVal<std::string>("foo") == "1.3");
      TEST_ASSERT(d.getVal<double>("foo2") == d2.getVal<double>("foo2"));
      TEST_ASSERT(d.getVal<std::vector<int>>("foo3") ==
                  d2.getVal<std::vector<int>>("foo3"));
    }

    {
      Dict d2 = d;
      d2.setVal("foo", 1);
      TEST_ASSERT(1 == d2.getVal<int>("foo"));
      TEST_ASSERT(d.getVal<double>("foo2") == d2.getVal<double>("foo2"));
      TEST_ASSERT(d.getVal<std::vector<int>>("foo3") ==
                  d2.getVal<std::vector<int>>("foo3"));
    }

    {
      Dict d2(d);
      TEST_ASSERT(d.getVal<double>("foo2") == d2.getVal<double>("foo2"));
      TEST_ASSERT(d.getVal<std::vector<int>>("foo3") ==
                  d2.getVal<std::vector<int>>("foo3"));
    }
  }
  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

class FooHandler : public CustomPropHandler {
 public:
  const char *getPropName() const override { return "Foo"; }
  bool canSerialize(const RDValue &value) const override {
    return rdvalue_is<Foo>(value);
  }
  bool read(std::istream &ss, RDValue &value) const override {
    int version = 0;
    streamRead(ss, version);
    Foo f;
    streamRead(ss, f.bar);
    streamRead(ss, f.baz);
    value = f;
    return true;
  }

  bool write(std::ostream &ss, const RDValue &value) const override {
    try {
      const Foo &f = rdvalue_cast<const Foo &>(value);
      const int version = 0;
      streamWrite(ss, version);
      streamWrite(ss, f.bar);
      streamWrite(ss, f.baz);
    } catch (boost::bad_any_cast &) {
      return false;
    }
    return true;
  }

  CustomPropHandler *clone() const override { return new FooHandler; }
};

void testCustomProps() {
  Foo f(1, 2.f);
  Dict d;
  d.setVal<Foo>("foo", f);
  RDValue &value = d.getData()[0].val;
  FooHandler foo_handler;
  std::vector<CustomPropHandler *> handlers = {&foo_handler,
                                               foo_handler.clone()};
  for (auto handler : handlers) {
    TEST_ASSERT(handler->canSerialize(value));
    RDValue bad_value = 1;
    TEST_ASSERT(!handler->canSerialize(bad_value));
    std::stringstream ss;
    TEST_ASSERT(handler->write(ss, value));
    RDValue newValue;
    TEST_ASSERT(handler->read(ss, newValue));
    TEST_ASSERT(from_rdvalue<const Foo &>(newValue).bar == f.bar);
    TEST_ASSERT(from_rdvalue<const Foo &>(newValue).baz == f.baz);
    newValue.destroy();
  }
  delete handlers[1];
}

void testGithub2910() {
  Dict d;
  d.setVal("foo", 1.0);
  d.clearVal("foo");
  d.clearVal("bar");
  d.clearVal("foo");
}

int main() {
  RDLog::InitLogs();
  testGithub940();

#if 1
  testRDAny();
  Dict d;
  INT_VECT fooV;
  fooV.resize(3);
  BOOST_LOG(rdInfoLog) << "dict test" << std::endl;
  TEST_ASSERT(!d.hasVal("foo"));
  int x = 1;
  d.setVal("foo", x);
  TEST_ASSERT(d.hasVal("foo"));
  TEST_ASSERT(!d.hasVal("bar"));
  int v, v2;
  d.getVal("foo", v);
  TEST_ASSERT(v == 1);
  v2 = d.getVal<int>("foo");
  TEST_ASSERT(v2 == v);
  d.setVal("bar", fooV);
  d.getVal("foo", v);
  TEST_ASSERT(v == 1);
  v2 = d.getVal<int>("foo");
  TEST_ASSERT(v2 == v);
  INT_VECT fooV2, fooV3;
  d.getVal("bar", fooV2);
  fooV3 = d.getVal<INT_VECT>("bar");
  TEST_ASSERT(fooV == fooV2);
  TEST_ASSERT(fooV2 == fooV3);

  VECT_INT_VECT fooV4;
  fooV4.resize(3);
  TEST_ASSERT(!d.hasVal("baz"));
  d.setVal("baz", fooV4);
  TEST_ASSERT(d.hasVal("baz"));

  DictCon dc1;
  TEST_ASSERT(!dc1.getDict()->hasVal("foo"));
  int y = 1;
  dc1.getDict()->setVal("foo", y);
  TEST_ASSERT(dc1.getDict()->hasVal("foo"));
  TEST_ASSERT(!dc1.getDict()->hasVal("bar"));
  dc1.getDict()->setVal("bar", fooV);
  dc1.getDict()->getVal("foo", v);
  TEST_ASSERT(v == 1);
  dc1.getDict()->getVal("bar", fooV2);
  TEST_ASSERT(fooV == fooV2);
  fooV4.resize(3);
  TEST_ASSERT(!dc1.getDict()->hasVal("baz"));
  dc1.getDict()->setVal("baz", fooV4);
  TEST_ASSERT(dc1.getDict()->hasVal("baz"));

  dc1.getDict()->reset();

  DictCon dc2 = dc1;
  TEST_ASSERT(!dc2.getDict()->hasVal("foo"));
  int z = 1;
  dc2.getDict()->setVal("foo", z);
  TEST_ASSERT(dc2.getDict()->hasVal("foo"));
  TEST_ASSERT(!dc2.getDict()->hasVal("bar"));
  dc2.getDict()->setVal("bar", fooV);
  dc2.getDict()->getVal("foo", v);
  TEST_ASSERT(v == 1);
  dc2.getDict()->getVal("bar", fooV2);
  TEST_ASSERT(fooV == fooV2);
  fooV4.resize(3);
  TEST_ASSERT(!dc2.getDict()->hasVal("baz"));
  dc2.getDict()->setVal("baz", fooV4);
  TEST_ASSERT(dc2.getDict()->hasVal("baz"));

  DictCon dc3(dc2);
  TEST_ASSERT(dc3.getDict()->hasVal("foo"));
  dc3.getDict()->getVal("foo", v);
  TEST_ASSERT(v == 1);
  dc3.getDict()->getVal("bar", fooV2);
  TEST_ASSERT(fooV == fooV2);
  fooV4.resize(3);
  TEST_ASSERT(dc3.getDict()->hasVal("baz"));

  TEST_ASSERT(dc3.getDict()->hasVal("foo"));
  dc3.getDict()->getVal("foo", v);
  TEST_ASSERT(v == 1);
  dc3.getDict()->getVal("bar", fooV2);
  TEST_ASSERT(fooV == fooV2);
  fooV4.resize(3);
  TEST_ASSERT(dc3.getDict()->hasVal("baz"));

  testStringVals();
  testVectToString();
#endif
  testConstReturns();
  testUpdate();
  testCustomProps();
  testGithub2910();
  return 0;
}
