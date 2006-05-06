//
// Copyright 2001-2006 Randal M. Henne, Greg Landrum and
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

  Dict d;
  INT_VECT fooV;
  fooV.resize(3);

  CHECK_INVARIANT(!d.hasVal("foo"),"bad init");
  int x = 1;
  d.setVal("foo", x);
  CHECK_INVARIANT(d.hasVal("foo"),"should be there");
  CHECK_INVARIANT(!d.hasVal("bar"),"bad other key");
  int v;
  d.getVal("foo",v);
  CHECK_INVARIANT(v==1,"bad val");
#if 1
  d.setVal("bar",fooV);
  d.getVal("foo",v);
  CHECK_INVARIANT(v==1,"bad val");
  INT_VECT fooV2;
  d.getVal("bar",fooV2);
  CHECK_INVARIANT(fooV==fooV2,"bad getVal");

  VECT_INT_VECT fooV3;
  fooV3.resize(3);
  CHECK_INVARIANT(!d.hasVal("baz"),"bad get");
  d.setVal("baz",fooV3);
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
  fooV3.resize(3);
  CHECK_INVARIANT(!dc1.getDict()->hasVal("baz"),"bad get");
  dc1.getDict()->setVal("baz",fooV3);
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
  fooV3.resize(3);
  CHECK_INVARIANT(!dc2.getDict()->hasVal("baz"),"bad get");
  dc2.getDict()->setVal("baz",fooV3);
  CHECK_INVARIANT(dc2.getDict()->hasVal("baz"),"bad get");

  DictCon dc3(dc2);
  CHECK_INVARIANT(dc3.getDict()->hasVal("foo"),"should be there");
  dc3.getDict()->getVal("foo",v);
  CHECK_INVARIANT(v==1,"bad val");
  dc3.getDict()->getVal("bar",fooV2);
  CHECK_INVARIANT(fooV==fooV2,"bad getVal");
  fooV3.resize(3);
  CHECK_INVARIANT(dc3.getDict()->hasVal("baz"),"bad get");

  
  //Atom newAt(Atom(5));
  CHECK_INVARIANT(dc3.getDict()->hasVal("foo"),"should be there");
  dc3.getDict()->getVal("foo",v);
  CHECK_INVARIANT(v==1,"bad val");
  dc3.getDict()->getVal("bar",fooV2);
  CHECK_INVARIANT(fooV==fooV2,"bad getVal");
  fooV3.resize(3);
  CHECK_INVARIANT(dc3.getDict()->hasVal("baz"),"bad get");
#endif  
  

  return 0;

}
