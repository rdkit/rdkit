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

#include "types.h"
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
  int bar;
  float baz;
  ~Foo() { std::cerr << "deleted!" << std::endl; }
};

void testRDAny() {
  std::cerr << "Testing RDValue" << std::endl;
  {
    RDAny v(-2147450880);
  }
  
  {
    int vi = 0;
    RDValue v(0);
    for (int i=0; i<100; ++i) {
        vi += i;
        v = rdvalue_cast<int>(v) + i;
        CHECK_INVARIANT(vi == rdvalue_cast<int>(v), "Opps, bad variant");
      }
               
  }
  std::cerr << "Testing RDAny" << std::endl;
  {
    RDAny a(1);
    RDAny b = a;
    CHECK_INVARIANT(rdany_cast<int>(a) == 1, "Should be 1");
    CHECK_INVARIANT(rdany_cast<int>(b) == 1, "Should be 1");
  }
  
  
  {
    RDAny a(1);
    RDAny b = a;
    CHECK_INVARIANT(rdany_cast<int>(a) == 1, "should be one");
    CHECK_INVARIANT(rdany_cast<int>(b) == rdany_cast<int>(a), "Bad Any");
    std::map<std::string, RDAny> foo;
    foo["foo"] = a;
    foo["bar"] = std::string("This is a test");
    CHECK_INVARIANT(rdany_cast<int>(foo["foo"]) == 1, "should be one");
    CHECK_INVARIANT(rdany_cast<int>(foo["foo"]) == rdany_cast<int>(a), "Bad Any");
    CHECK_INVARIANT(rdany_cast<std::string>(foo["bar"]) == "This is a test", "Bad Any");
  }

  {
    bool a = true;
    RDValue v(a);
    CHECK_INVARIANT(rdvalue_cast<bool>(v) == true, "bad value cast");
    v = (int) 10;
    CHECK_INVARIANT(rdvalue_cast<int>(v) == 10, "bad value cast");
  }
  
  {
    Dict d;
    bool a=true;
    d.setVal("foo", a);
    d.getVal<bool>("foo");
  }

  { // tests computed props
    STR_VECT computed;
    Dict d;
    d.setVal(detail::computedPropName, computed);
    computed.push_back("foo");
    d.setVal(detail::computedPropName, computed);
    STR_VECT computed2 = d.getVal<STR_VECT>(detail::computedPropName);
    CHECK_INVARIANT(computed2[0] == "foo", "bad STR_VECT");
    Dict d2(d);
    computed2 = d2.getVal<STR_VECT>(detail::computedPropName);
    CHECK_INVARIANT(computed2[0] == "foo", "bad STR_VECT");
  }
  
  {
    Dict d;
    //int v=1;
    //d.setVal("foo", v);
    //CHECK_INVARIANT(d.getVal<int>("foo") == 1, "bad getval");
    
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
      fooV2 = rdany_cast<std::vector<int> >(a);
      CHECK_INVARIANT(fooV == fooV2, "bad getVal");
    }
    
    {
      std::vector<int> fooV2;
      std::cerr << "dict set int vect" << std::endl;
      d.setVal("bar", fooV);
      std::cerr << "dict get int vect" << std::endl;
      d.getVal("bar", fooV2);
      CHECK_INVARIANT(fooV == fooV2, "bad getVal");
    }
  }
  
  {
    std::vector<int> v;
    for(int i=0;i<4;++i)
      v.push_back(i);

    RDAny foo(v);
    RDAny bar = foo;
    RDAny baz(foo);
    
    for(int i=0;i<4;++i) { 
      CHECK_INVARIANT(rdany_cast<std::vector<int> >(foo)[i] == i, "Failed check");
      CHECK_INVARIANT(rdany_cast<std::vector<int> >(bar)[i] == i, "Failed check");
      CHECK_INVARIANT(rdany_cast<std::vector<int> >(baz)[i] == i, "Failed check");
    }
    
  }
  
  {
    std::vector<double> v;
    for(double i=0;i<4;++i)
      v.push_back(i);

    RDAny foo(v);

    for(int i=0;i<4;++i) { 
      CHECK_INVARIANT(rdany_cast<std::vector<double> >(foo)[i] == i, "Failed check");
    }
    
    RDAny b = foo;

    for(int i=0;i<4;++i) { 
      CHECK_INVARIANT(rdany_cast<std::vector<double> >(b)[i] == i, "Failed check");
    }
  }
  const int loops = 10000000;
  {
    std::clock_t clock1 = std::clock();
    boost::any v;
    for(int i=0;i<loops;++i) {
      v = i;
    }
    std::clock_t clock2 = std::clock();

    std::cout << "static boost any:" << (double)(clock2-clock1)/CLOCKS_PER_SEC << " s" << std::endl;
  }
  {
    std::clock_t clock1 = std::clock();
    boost::any *v=0, *vv;
    for(int i=0;i<loops;++i) {
      vv = new boost::any(v?boost::any_cast<int>(*v) + i: i);
      delete v;
      v = vv;
    }
    delete vv;
    std::clock_t clock2 = std::clock();

    std::cout << "dynamic boost any:" << (double)(clock2-clock1)/CLOCKS_PER_SEC << " s" << std::endl;
  }

  {
    std::clock_t clock1 = std::clock();
    RDAny v;
    for(int i=0;i<loops;++i) {
      v = i;
    }
    std::clock_t clock2 = std::clock();

    std::cout << "static RDAny:" << (double)(clock2-clock1)/CLOCKS_PER_SEC << " s" << std::endl;
  }

  {
    std::clock_t clock1 = std::clock();
    RDAny *v=0, *vv;
    for(int i=0;i<loops;++i) {
      vv = new RDAny(v ?rdany_cast<int>(*v) + i : i);
      delete v;
      v = vv;
    }
    delete vv;
    std::clock_t clock2 = std::clock();

    std::cout << "dynamic RDAny:" << (double)(clock2-clock1)/CLOCKS_PER_SEC << " s" << std::endl;
  }

  {
    std::clock_t clock1 = std::clock();
    RDValue v;
    for(int i=0;i<loops;++i) {
      v = i;
    }
    std::clock_t clock2 = std::clock();

    std::cout << "static RDValue:" << (double)(clock2-clock1)/CLOCKS_PER_SEC << " s" << std::endl;
  }

  {
    std::clock_t clock1 = std::clock();
    RDValue v(0);
    for(int i=0;i<loops;++i) {
      v = RDValue(rdvalue_cast<int>(v) + i);
    }

    std::clock_t clock2 = std::clock();

    std::cout << "dynamic RDValue:" << (double)(clock2-clock1)/CLOCKS_PER_SEC << " s" << std::endl;
  }
  
  { // checks replacement with vector
    RDAny vv(2.0);
    CHECK_INVARIANT(rdany_cast<double>(vv) == 2.0, "Bad double");
    
    
    std::vector<int> vect;
    vect.push_back(1);
    vv = vect;
    CHECK_INVARIANT(rdany_cast<std::vector<int> >(vv)[0] == 1, "Bad cast");
    
    // tests copy
    RDAny vvv(vv);
    
    CHECK_INVARIANT(rdany_cast<std::vector<int> >(vvv)[0] == 1, "Bad cast");    
  }

  {
    // Checks fallback to Any
    std::vector<std::pair<int,int> > pvect;
    pvect.push_back(std::make_pair<int,int>(2,2));
    boost::any any1(pvect);
    boost::any_cast<std::vector<std::pair<int,int> > >(any1);
    boost::any_cast<std::vector<std::pair<int,int> > &>(any1);    
    boost::any_cast<const std::vector<std::pair<int,int> > &>(any1);
    
    RDAny vv(pvect);
    boost::any &any = rdany_cast<boost::any&>(vv);
    boost::any_cast<std::vector<std::pair<int,int> > >(any);
    boost::any_cast<std::vector<std::pair<int,int> > &>(any);    
    boost::any_cast<const std::vector<std::pair<int,int> > &>(any);
    
    const std::vector<std::pair<int,int> > &pv = rdany_cast<std::vector<std::pair<int, int> > >(vv);
    CHECK_INVARIANT(pv[0].first == 2,
                    "Bad cast");
    RDAny vvv(vv);
    CHECK_INVARIANT((rdany_cast<std::vector<std::pair<int, int> > >(vvv)[0].first == 2),
                    "Bad cast");
    
  }

  {
    // Check pointers -- RDAny doesn't delete these, must do them manually
    std::vector<int> *p = new std::vector<int>();
    p->push_back(100);
    RDAny v(p);
    RDAny vv(v);
    try {
      rdany_cast<std::vector<int> >(v);
#ifndef UNSAFE_RDVALUE
      PRECONDITION(0, "Should throw bad cast");
#endif
    } catch (boost::bad_any_cast &e) {
    }
    
    CHECK_INVARIANT((*rdany_cast<std::vector<int> *>(vv))[0] == 100, "Bad cast");
    CHECK_INVARIANT((*rdany_cast<std::vector<int> *>((const RDAny&)vv))[0] == 100, "Bad cast");
    delete p;

    std::map<int,int> *m = new std::map<int,int>();
    (*m)[0] = 1;
    RDAny mv(m);
    // leaks
    std::map<int,int> *anym = rdany_cast<std::map<int,int> *>(mv);
    CHECK_INVARIANT(anym->find(0) != anym->end(),
                    "Bad cast");    
    delete anym;
  }

  {
    // check shared ptrs -- boost::any deletes these :)
    typedef boost::shared_ptr<std::vector<int> > vptr;
    vptr p(new std::vector<int>());
    p->push_back(100);
    RDAny v(p);
    RDAny vv(v);
    CHECK_INVARIANT((*rdany_cast<vptr>(v))[0] == 100, "Bad cast");    
    CHECK_INVARIANT((*rdany_cast<vptr>(vv))[0] == 100, "Bad cast");
    CHECK_INVARIANT((*rdany_cast<vptr>((const RDAny&)vv))[0] == 100, "Bad cast");

    typedef boost::shared_ptr<std::map<int,int> > mptr;
    mptr m(new std::map<int,int>());
    (*m)[0] = 1;
    RDAny mv(m);
    // leaks
    mptr anym = rdany_cast<mptr>(mv);
    CHECK_INVARIANT(anym->find(0) != anym->end(),
                    "Bad cast");

    RDAny any3(boost::shared_ptr<Foo>( new Foo ));
    CHECK_INVARIANT(any3.m_value.getTag() == RDTypeTag::AnyTag, "Wrong type");
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
  }
  {
    Dict d;
    std::vector<double> v;
    v.push_back(1.2);
    v.push_back(0);
    d.setVal("foo", v);
    std::string sv;
    d.getVal("foo", sv);
    std::cerr << sv << std::endl;
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
    std::cerr << sv << std::endl;
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

int main() {
  RDLog::InitLogs();
    testRDAny();

#if 1
  Dict d;
  INT_VECT fooV;
  fooV.resize(3);
  BOOST_LOG(rdInfoLog) << "dict test" << std::endl;
  CHECK_INVARIANT(!d.hasVal("foo"), "bad init");
  int x = 1;
  d.setVal("foo", x);
  CHECK_INVARIANT(d.hasVal("foo"), "should be there");
  CHECK_INVARIANT(!d.hasVal("bar"), "bad other key");
  int v, v2;
  d.getVal("foo", v);
  CHECK_INVARIANT(v == 1, "bad val");
  v2 = d.getVal<int>("foo");
  CHECK_INVARIANT(v2 == v, "bad val");
  d.setVal("bar", fooV);
  d.getVal("foo", v);
  CHECK_INVARIANT(v == 1, "bad val");
  v2 = d.getVal<int>("foo");
  CHECK_INVARIANT(v2 == v, "bad val");
  INT_VECT fooV2, fooV3;
  d.getVal("bar", fooV2);
  fooV3 = d.getVal<INT_VECT>("bar");
  CHECK_INVARIANT(fooV == fooV2, "bad getVal");
  CHECK_INVARIANT(fooV2 == fooV3, "bad getVal");

  VECT_INT_VECT fooV4;
  fooV4.resize(3);
  CHECK_INVARIANT(!d.hasVal("baz"), "bad get");
  d.setVal("baz", fooV4);
  CHECK_INVARIANT(d.hasVal("baz"), "bad get");

  DictCon dc1;
  CHECK_INVARIANT(!dc1.getDict()->hasVal("foo"), "bad init");
  int y = 1;
  dc1.getDict()->setVal("foo", y);
  CHECK_INVARIANT(dc1.getDict()->hasVal("foo"), "should be there");
  CHECK_INVARIANT(!dc1.getDict()->hasVal("bar"), "bad other key");
  dc1.getDict()->setVal("bar", fooV);
  dc1.getDict()->getVal("foo", v);
  CHECK_INVARIANT(v == 1, "bad val");
  dc1.getDict()->getVal("bar", fooV2);
  CHECK_INVARIANT(fooV == fooV2, "bad getVal");
  fooV4.resize(3);
  CHECK_INVARIANT(!dc1.getDict()->hasVal("baz"), "bad get");
  dc1.getDict()->setVal("baz", fooV4);
  CHECK_INVARIANT(dc1.getDict()->hasVal("baz"), "bad get");

  dc1.getDict()->reset();

  DictCon dc2 = dc1;
  CHECK_INVARIANT(!dc2.getDict()->hasVal("foo"), "bad init");
  int z = 1;
  dc2.getDict()->setVal("foo", z);
  CHECK_INVARIANT(dc2.getDict()->hasVal("foo"), "should be there");
  CHECK_INVARIANT(!dc2.getDict()->hasVal("bar"), "bad other key");
  dc2.getDict()->setVal("bar", fooV);
  dc2.getDict()->getVal("foo", v);
  CHECK_INVARIANT(v == 1, "bad val");
  dc2.getDict()->getVal("bar", fooV2);
  CHECK_INVARIANT(fooV == fooV2, "bad getVal");
  fooV4.resize(3);
  CHECK_INVARIANT(!dc2.getDict()->hasVal("baz"), "bad get");
  dc2.getDict()->setVal("baz", fooV4);
  CHECK_INVARIANT(dc2.getDict()->hasVal("baz"), "bad get");

  DictCon dc3(dc2);
  CHECK_INVARIANT(dc3.getDict()->hasVal("foo"), "should be there");
  dc3.getDict()->getVal("foo", v);
  CHECK_INVARIANT(v == 1, "bad val");
  dc3.getDict()->getVal("bar", fooV2);
  CHECK_INVARIANT(fooV == fooV2, "bad getVal");
  fooV4.resize(3);
  CHECK_INVARIANT(dc3.getDict()->hasVal("baz"), "bad get");

  CHECK_INVARIANT(dc3.getDict()->hasVal("foo"), "should be there");
  dc3.getDict()->getVal("foo", v);
  CHECK_INVARIANT(v == 1, "bad val");
  dc3.getDict()->getVal("bar", fooV2);
  CHECK_INVARIANT(fooV == fooV2, "bad getVal");
  fooV4.resize(3);
  CHECK_INVARIANT(dc3.getDict()->hasVal("baz"), "bad get");

  testStringVals();
  testVectToString();
#endif
  testConstReturns();
  
  return 0;
}
