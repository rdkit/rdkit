// $Id$
//
// Copyright (C)  2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <iostream>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>
#include <Geometry/point.h>

#include <Features/Feature.h>

#include <boost/spirit/core.hpp>
#include <boost/spirit/actor/assign_actor.hpp>
using namespace boost::spirit;

using namespace RDKit;
using namespace RDGeom;
using namespace RDFeatures;

typedef enum {
  fooType,
  barType,
  bazType,
  grnType,
} TypeMarker;

void test1(){
  std::cerr << "-------------------------------------" << std::endl;
  std::cerr << "Basics for ExplicitFeatures." << std::endl;

  ExplicitFeature<TypeMarker> f1;
  f1.setFamily(bazType);
  TEST_ASSERT(f1.getFamily()==bazType);
  f1.setType(grnType);
  TEST_ASSERT(f1.getType()==grnType);
  TEST_ASSERT(feq(f1.getLoc().x,0.0));
  TEST_ASSERT(feq(f1.getLoc().y,0.0));
  TEST_ASSERT(feq(f1.getLoc().z,0.0));
  TEST_ASSERT(f1.getDirs().size()==0);

  f1 = ExplicitFeature<TypeMarker>(barType,fooType);
  TEST_ASSERT(f1.getFamily()==barType);
  TEST_ASSERT(f1.getType()==fooType);
  TEST_ASSERT(feq(f1.getLoc().x,0.0));
  TEST_ASSERT(feq(f1.getLoc().y,0.0));
  TEST_ASSERT(feq(f1.getLoc().z,0.0));
  TEST_ASSERT(f1.getDirs().size()==0);

  f1 = ExplicitFeature<TypeMarker>(barType,fooType,Point3D(1.0,2.0,3.0));
  TEST_ASSERT(f1.getFamily()==barType);
  TEST_ASSERT(f1.getType()==fooType);
  TEST_ASSERT(feq(f1.getLoc().x,1.0));
  TEST_ASSERT(feq(f1.getLoc().y,2.0));
  TEST_ASSERT(feq(f1.getLoc().z,3.0));
  TEST_ASSERT(f1.getDirs().size()==0);

  f1.setLoc(Point3D(-1.0,-2.0,-3.0));
  TEST_ASSERT(feq(f1.getLoc().x,-1.0));
  TEST_ASSERT(feq(f1.getLoc().y,-2.0));
  TEST_ASSERT(feq(f1.getLoc().z,-3.0));
  TEST_ASSERT(f1.getDirs().size()==0);

  
  std::cerr << "  done" << std::endl;
}

void test2(){
  std::cerr << "-------------------------------------" << std::endl;
  std::cerr << "Basics for ImplicitFeatures." << std::endl;

  ImplicitFeature<TypeMarker> f1;
  f1.setFamily(fooType);
  TEST_ASSERT(f1.getType()==fooType);
  f1.setType(grnType);
  TEST_ASSERT(f1.getType()==grnType);
  TEST_ASSERT(feq(f1.getLoc().x,0.0));
  TEST_ASSERT(feq(f1.getLoc().y,0.0));
  TEST_ASSERT(feq(f1.getLoc().z,0.0));
  TEST_ASSERT(f1.getDirs().size()==0);

  f1 = ImplicitFeature<TypeMarker>(barType,fooType);
  TEST_ASSERT(f1.getFamily()==barType);
  TEST_ASSERT(f1.getType()==barType);
  TEST_ASSERT(feq(f1.getLoc().x,0.0));
  TEST_ASSERT(feq(f1.getLoc().y,0.0));
  TEST_ASSERT(feq(f1.getLoc().z,0.0));
  TEST_ASSERT(f1.getDirs().size()==0);

  Point3D p1(0,0,0),p2(1,0,0),p3(0,1,0);
  f1.addPoint(&p1);
  f1.addPoint(&p2);
  TEST_ASSERT(feq(f1.getLoc().x,0.50));
  TEST_ASSERT(feq(f1.getLoc().y,0.0));
  TEST_ASSERT(feq(f1.getLoc().z,0.0));
  
  f1.addPoint(&p3);
  TEST_ASSERT(feq(f1.getLoc().x,0.3333));
  TEST_ASSERT(feq(f1.getLoc().y,0.3333));
  TEST_ASSERT(feq(f1.getLoc().z,0.0));

  f1.reset();
  TEST_ASSERT(feq(f1.getLoc().x,0.0));
  TEST_ASSERT(feq(f1.getLoc().y,0.0));
  TEST_ASSERT(feq(f1.getLoc().z,0.0));

  
  std::cerr << "  done" << std::endl;
}


void test3(){
  std::cerr << "-------------------------------------" << std::endl;
  std::cerr << "ExplicitFeatures 2D, string type." << std::endl;

  typedef ExplicitFeature<std::string,std::string,Point2D> LocalFeature;
  LocalFeature f1;
  f1.setType("foo");
  TEST_ASSERT(f1.getType()=="foo");
  f1.setFamily("foob");
  TEST_ASSERT(f1.getFamily()=="foob");
  TEST_ASSERT(feq(f1.getLoc().x,0.0));
  TEST_ASSERT(feq(f1.getLoc().y,0.0));
  TEST_ASSERT(f1.getDirs().size()==0);

  f1 = LocalFeature("foo","bar");
  TEST_ASSERT(f1.getFamily()=="bar");
  TEST_ASSERT(f1.getType()=="bar");
  TEST_ASSERT(feq(f1.getLoc().x,0.0));
  TEST_ASSERT(feq(f1.getLoc().y,0.0));
  TEST_ASSERT(f1.getDirs().size()==0);

  f1 = LocalFeature("grm","grn",Point2D(1.0,2.0));
  TEST_ASSERT(f1.getFamily()=="grm");
  TEST_ASSERT(f1.getType()=="grn");
  TEST_ASSERT(feq(f1.getLoc().x,1.0));
  TEST_ASSERT(feq(f1.getLoc().y,2.0));
  TEST_ASSERT(f1.getDirs().size()==0);

  f1.setLoc(Point2D(-1.0,-2.0));
  TEST_ASSERT(feq(f1.getLoc().x,-1.0));
  TEST_ASSERT(feq(f1.getLoc().y,-2.0));
  TEST_ASSERT(f1.getDirs().size()==0);
  
  std::cerr << "  done" << std::endl;
}

void testParser(){
  std::cerr << "-------------------------------------" << std::endl;
  std::cerr << "parser testing." << std::endl;

  std::string text;
  int n;
  bool ok;

  text="p1";
  ok=parse(text.c_str(),
		( ch_p('p') >> int_p[assign_a(n)] ),space_p).full;
  TEST_ASSERT(ok);
  TEST_ASSERT(n==1);

  text="p12";
  ok=parse(text.c_str(),
		( ch_p('p') >> int_p[assign_a(n)] ),space_p).full;
  TEST_ASSERT(ok);
  TEST_ASSERT(n==12);
  
  text="p2 + p3";
  ok=parse(text.c_str(),
		( ch_p('p') >> int_p[assign_a(n)] ),space_p).full;
  TEST_ASSERT(!ok);

  
  

  std::cerr << "  done" << std::endl;
}
int main(){
#if 0
  test1();
  test2();
  test3();
#endif
  testParser();
}
