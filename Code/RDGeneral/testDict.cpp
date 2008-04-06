//
// Copyright 2001-2008 Randal M. Henne, Greg Landrum and
//                     Rational Discovery LLC 
//
//  @@ All Rights Reserved @@
//
//

#include "types.h"
#include <RDGeneral/Invariant.h>

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

int main(){
  RDLog::InitLogs();
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
  

  return 0;

}
