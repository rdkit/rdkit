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
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>
#include <ctime>

using namespace RDKit;
using namespace std;

class DictCon{
public:
  DictCon() { d.reset(); };
  DictCon(const DictCon &other) {d = other.d;};
  DictCon &operator=(const DictCon &other) { d= other.d; return *this;};
  Dict *getDict() { return &d;};
private:
  Dict d;
};

void testStringVals(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing String Pickle Roundtrips." << std::endl;
  {
    Dict d;
    std::string sv;
    sv="1";
    d.setVal("foo",sv);
    int iv;
    d.getVal("foo",iv);
    TEST_ASSERT(iv==1);
  }
  {
    Dict d;
    d.setVal("foo","1");
    int iv;
    d.getVal("foo",iv);
    TEST_ASSERT(iv==1);
  }
  {
    Dict d;
    std::string sv;
    sv="1.3";
    d.setVal("foo",sv);
    double dv;
    d.getVal("foo",dv);
    TEST_ASSERT(feq(dv,1.3));
  }
  
  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void testVectToString(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing conversion of vect to string." << std::endl;
  {
    Dict d;
    std::vector<int> v;
    v.push_back(1);
    v.push_back(0);
    d.setVal("foo",v);
    std::string sv;
    d.getVal("foo",sv);
    TEST_ASSERT(sv=="[1,0,]");
  }
  {
    Dict d;
    std::vector<unsigned int> v;
    v.push_back(1);
    v.push_back(0);
    d.setVal("foo",v);
    std::string sv;
    d.getVal("foo",sv);
    TEST_ASSERT(sv=="[1,0,]");
  }
  {
    Dict d;
    std::vector<double> v;
    v.push_back(1);
    v.push_back(0);
    d.setVal("foo",v);
    bool ok=false;
    try{
      std::string sv;
      d.getVal("foo",sv);
    } catch (const boost::bad_any_cast &) {
      ok=true;
    }
    TEST_ASSERT(ok);
  }
  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}


void testConstReturns(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing returning const references." << std::endl;
  {
    std::string v="foo";
    boost::any anyv(v);

    std::string tgt=boost::any_cast<std::string>(anyv);
    const std::string &ctgt=boost::any_cast<const std::string &>(anyv);
    TEST_ASSERT(ctgt!="");
  }

  {
    Dict d;
    std::string v="foo";
    d.setVal("foo",v);
    
    //const std::string nv=d.getVal<const std::string &>("foo");
    std::string nv=d.getVal<std::string>("foo");
    TEST_ASSERT(nv=="foo");
    
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
      const std::string &nv=d.getVal<const std::string &>("foo");
      ls+= nv.size();
    }
    BOOST_LOG(rdErrorLog) << "done: "<<ls << std::endl;    
    //std::string nv=d.getVal<std::string>("foo");
  }
#else
  {
    //int nreps=100000000;
    int nreps=100000;
    Dict d;
    std::string v="foo";
    boost::any anyv(v);
    d.setVal("foo",v);

    std::clock_t start,end;

    double ls=0;
    BOOST_LOG(rdErrorLog) << "any cast" << std::endl;    
    start = std::clock();
    for(int i=0;i<nreps;++i){
      const std::string &nv=boost::any_cast<const std::string &>(anyv);
      ls+= nv.size();
    }
    end = std::clock();
    BOOST_LOG(rdErrorLog) << "done: "<<(end-start)/(double)(CLOCKS_PER_SEC)<<" "<<ls << std::endl;    
    
    ls=0;
    BOOST_LOG(rdErrorLog) << "copy" << std::endl;    
    start = std::clock();
    for(int i=0;i<nreps;++i){
      std::string nv=d.fromany<std::string>(anyv);
      ls+= nv.size();
    }
    end = std::clock();
    BOOST_LOG(rdErrorLog) << "done: "<<(end-start)/(double)(CLOCKS_PER_SEC)<<" "<<ls << std::endl;    

    ls=0;
    BOOST_LOG(rdErrorLog) << "ref" << std::endl;    
    start = std::clock();
    for(int i=0;i<nreps;++i){
      const std::string &nv=d.fromany<const std::string &>(anyv);
      ls+= nv.size();
    }
    end = std::clock();
    BOOST_LOG(rdErrorLog) << "done: "<<(end-start)/(double)(CLOCKS_PER_SEC)<<" "<<ls << std::endl;    


    ls=0;
    BOOST_LOG(rdErrorLog) << "dict" << std::endl;    
    start = std::clock();
    for(int i=0;i<nreps;++i){
      const std::string &nv=d.getVal<const std::string &>("foo");
      ls+= nv.size();
    }
    end = std::clock();
    BOOST_LOG(rdErrorLog) << "done: "<<(end-start)/(double)(CLOCKS_PER_SEC)<<" "<<ls << std::endl;    

    ls=0;
    BOOST_LOG(rdErrorLog) << "ref with hasVal" << std::endl;    
    start = std::clock();
    std::string k="foo";
    for(int i=0;i<nreps;++i){
      if(d.hasVal(k)){
        const std::string &nv=d.fromany<const std::string &>(anyv);
        ls+= nv.size();
      }
    }
    end = std::clock();
    BOOST_LOG(rdErrorLog) << "done: "<<(end-start)/(double)(CLOCKS_PER_SEC)<<" "<<ls << std::endl;    


    //std::string nv=d.getVal<std::string>("foo");
  }

#endif
  
  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}


int main(){
  RDLog::InitLogs();
#if 1
  Dict d;
  INT_VECT fooV;
  fooV.resize(3);
  BOOST_LOG(rdInfoLog) << "dict test" << std::endl;
  CHECK_INVARIANT(!d.hasVal("foo"),"bad init");
  int x = 1;
  d.setVal("foo", x);
  CHECK_INVARIANT(d.hasVal("foo"),"should be there");
  CHECK_INVARIANT(!d.hasVal("bar"),"bad other key");
  int v,v2;
  d.getVal("foo",v);
  CHECK_INVARIANT(v==1,"bad val");
  v2=d.getVal<int>("foo");
  CHECK_INVARIANT(v2==v,"bad val");
  d.setVal("bar",fooV);
  d.getVal("foo",v);
  CHECK_INVARIANT(v==1,"bad val");
  v2=d.getVal<int>("foo");
  CHECK_INVARIANT(v2==v,"bad val");
  INT_VECT fooV2,fooV3;
  d.getVal("bar",fooV2);
  fooV3=d.getVal<INT_VECT>("bar");
  CHECK_INVARIANT(fooV==fooV2,"bad getVal");
  CHECK_INVARIANT(fooV2==fooV3,"bad getVal");
  

  VECT_INT_VECT fooV4;
  fooV4.resize(3);
  CHECK_INVARIANT(!d.hasVal("baz"),"bad get");
  d.setVal("baz",fooV4);
  CHECK_INVARIANT(d.hasVal("baz"),"bad get");

  DictCon dc1;
  CHECK_INVARIANT(!dc1.getDict()->hasVal("foo"),"bad init");
  int y = 1;
  dc1.getDict()->setVal("foo",y);
  CHECK_INVARIANT(dc1.getDict()->hasVal("foo"),"should be there");
  CHECK_INVARIANT(!dc1.getDict()->hasVal("bar"),"bad other key");
  dc1.getDict()->setVal("bar",fooV);
  dc1.getDict()->getVal("foo",v);
  CHECK_INVARIANT(v==1,"bad val");
  dc1.getDict()->getVal("bar",fooV2);
  CHECK_INVARIANT(fooV==fooV2,"bad getVal");
  fooV4.resize(3);
  CHECK_INVARIANT(!dc1.getDict()->hasVal("baz"),"bad get");
  dc1.getDict()->setVal("baz",fooV4);
  CHECK_INVARIANT(dc1.getDict()->hasVal("baz"),"bad get");

  dc1.getDict()->reset();

  DictCon dc2=dc1;
  CHECK_INVARIANT(!dc2.getDict()->hasVal("foo"),"bad init");
  int z = 1;
  dc2.getDict()->setVal("foo",z);
  CHECK_INVARIANT(dc2.getDict()->hasVal("foo"),"should be there");
  CHECK_INVARIANT(!dc2.getDict()->hasVal("bar"),"bad other key");
  dc2.getDict()->setVal("bar",fooV);
  dc2.getDict()->getVal("foo",v);
  CHECK_INVARIANT(v==1,"bad val");
  dc2.getDict()->getVal("bar",fooV2);
  CHECK_INVARIANT(fooV==fooV2,"bad getVal");
  fooV4.resize(3);
  CHECK_INVARIANT(!dc2.getDict()->hasVal("baz"),"bad get");
  dc2.getDict()->setVal("baz",fooV4);
  CHECK_INVARIANT(dc2.getDict()->hasVal("baz"),"bad get");

  DictCon dc3(dc2);
  CHECK_INVARIANT(dc3.getDict()->hasVal("foo"),"should be there");
  dc3.getDict()->getVal("foo",v);
  CHECK_INVARIANT(v==1,"bad val");
  dc3.getDict()->getVal("bar",fooV2);
  CHECK_INVARIANT(fooV==fooV2,"bad getVal");
  fooV4.resize(3);
  CHECK_INVARIANT(dc3.getDict()->hasVal("baz"),"bad get");

  CHECK_INVARIANT(dc3.getDict()->hasVal("foo"),"should be there");
  dc3.getDict()->getVal("foo",v);
  CHECK_INVARIANT(v==1,"bad val");
  dc3.getDict()->getVal("bar",fooV2);
  CHECK_INVARIANT(fooV==fooV2,"bad getVal");
  fooV4.resize(3);
  CHECK_INVARIANT(dc3.getDict()->hasVal("baz"),"bad get");
  
  testStringVals();
  testVectToString();
#endif
  testConstReturns();

  return 0;

}
