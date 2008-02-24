// $Id$
//
// Copyright (c) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#include "QueryObjects.h"
#include <cmath>
#include <istream>
#include <sstream>
#include <boost/lexical_cast.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>


using namespace Queries;


void test1(){
  std::cout << "Float" << std::endl;
  EqualityQuery<double> q(1.0);

  CHECK_INVARIANT(!q.Match(0.0),"");
  CHECK_INVARIANT(q.Match(1.0),"");
  CHECK_INVARIANT(!q.Match(1.001),"");
  CHECK_INVARIANT(!q.Match(1.1),"");
  CHECK_INVARIANT(!q.Match(-2),"");


  std::cout << "With Tolerance" << std::endl;
  q.setTol(0.002);
  CHECK_INVARIANT(!q.Match(0.0),"");
  CHECK_INVARIANT(q.Match(1.0),"");
  CHECK_INVARIANT(q.Match(1.001),"");
  CHECK_INVARIANT(!q.Match(1.1),"");
  CHECK_INVARIANT(!q.Match(-2),"");

  Query<double> *newQ;
  newQ = q.copy();
  CHECK_INVARIANT(!newQ->Match(0.0),"");
  CHECK_INVARIANT(newQ->Match(1.0),"");
  CHECK_INVARIANT(newQ->Match(1.001),"");
  CHECK_INVARIANT(!newQ->Match(1.1),"");
  CHECK_INVARIANT(!newQ->Match(-2),"");
}

void test2(){
  std::cout << "Set" << std::endl;
  SetQuery<int> q;
  q.insert(1);
  q.insert(3);
  q.insert(5);
  
  CHECK_INVARIANT(!q.Match(0),"");
  CHECK_INVARIANT(q.Match(1),"");
  CHECK_INVARIANT(q.Match(3),"");
  CHECK_INVARIANT(!q.Match(-3),"");

  Query<int> *newQ;
  newQ = q.copy();
  CHECK_INVARIANT(!newQ->Match(0),"");
  CHECK_INVARIANT(newQ->Match(1),"");
  CHECK_INVARIANT(newQ->Match(3),"");
  CHECK_INVARIANT(!newQ->Match(-3),"");

}

void test3(){
  std::cout << "And" << std::endl;
  AndQuery<int> *q = new AndQuery<int>;
  LessQuery<int> *l = new LessQuery<int>;
  l->setVal(0);
  GreaterEqualQuery<int> *g = new GreaterEqualQuery<int>;
  g->setVal(3);

  q->addChild(Query<int>::CHILD_TYPE(l));
  q->addChild(Query<int>::CHILD_TYPE(g));

  CHECK_INVARIANT(!q->Match(0),"");
  CHECK_INVARIANT(q->Match(1),"");
  CHECK_INVARIANT(q->Match(3),"");
  CHECK_INVARIANT(!q->Match(-3),"");

  Query<int> *newQ;
  newQ = q->copy();
  CHECK_INVARIANT(!newQ->Match(0),"");
  CHECK_INVARIANT(newQ->Match(1),"");
  CHECK_INVARIANT(newQ->Match(3),"");
  CHECK_INVARIANT(!newQ->Match(-3),"");
  
}

void test4(){
  std::cout << "Or" << std::endl;
  OrQuery<int> *q = new OrQuery<int>;
  LessQuery<int> *l = new LessQuery<int>;
  l->setVal(0);
  GreaterEqualQuery<int> *g = new GreaterEqualQuery<int>;
  g->setVal(3);

  q->addChild(Query<int>::CHILD_TYPE(l));
  q->addChild(Query<int>::CHILD_TYPE(g));

  CHECK_INVARIANT(q->Match(0),"");
  CHECK_INVARIANT(q->Match(1),"");
  CHECK_INVARIANT(q->Match(3),"");
  CHECK_INVARIANT(q->Match(-3),"");

  Query<int> *newQ;
  newQ = q->copy();
  CHECK_INVARIANT(newQ->Match(0),"");
  CHECK_INVARIANT(newQ->Match(1),"");
  CHECK_INVARIANT(newQ->Match(3),"");
  CHECK_INVARIANT(newQ->Match(-3),"");

}

void test5(){
  std::cout << "XOr" << std::endl;
  XOrQuery<int> *q = new XOrQuery<int>;
  LessQuery<int> *l = new LessQuery<int>;
  l->setVal(0);
  GreaterEqualQuery<int> *g = new GreaterEqualQuery<int>;
  g->setVal(3);

  q->addChild(Query<int>::CHILD_TYPE(l));
  q->addChild(Query<int>::CHILD_TYPE(g));

  CHECK_INVARIANT(q->Match(-1),"");
  CHECK_INVARIANT(q->Match(0),"");
  CHECK_INVARIANT(!q->Match(1),"");
  CHECK_INVARIANT(!q->Match(3),"");
  CHECK_INVARIANT(q->Match(-3),"");

  Query<int> *newQ;
  newQ = q->copy();
  CHECK_INVARIANT(newQ->Match(-1),"");
  CHECK_INVARIANT(newQ->Match(0),"");
  CHECK_INVARIANT(!newQ->Match(1),"");
  CHECK_INVARIANT(!newQ->Match(3),"");
  CHECK_INVARIANT(newQ->Match(-3),"");
}

int foofun(double bar) { return int(floor(bar)); };

void test6(){
  std::cout << "pointer and copy foo" << std::endl;
  EqualityQuery<int,double,true> q;
  q.setDataFunc(foofun);
  q.setVal(6);
  CHECK_INVARIANT(q.Match(6.0),"");
  CHECK_INVARIANT(q.Match(6.1),"");
  CHECK_INVARIANT(!q.Match(5.0),"");

  Query<int,double,true> *newQ;
  newQ = q.copy();
  CHECK_INVARIANT(newQ->Match(6.0),"");
  CHECK_INVARIANT(newQ->Match(6.1),"");
  CHECK_INVARIANT(!newQ->Match(5.0),"");

  Query<int,double,true> *newQ2 = &q;
  CHECK_INVARIANT(newQ2->Match(6.0),"");
  CHECK_INVARIANT(newQ2->Match(6.1),"");
  CHECK_INVARIANT(!newQ2->Match(5.0),"");

  Query<int,double,true> *newQ3;
  newQ3 = newQ2->copy();
  CHECK_INVARIANT(newQ3->Match(6.0),"");
  CHECK_INVARIANT(newQ3->Match(6.1),"");
  CHECK_INVARIANT(!newQ3->Match(5.0),"");

  
}

bool matchF(int v){
  return v == 3;
}

int dataF(float v){
  return int(floor(v))*3;
}

bool cmp(int v){
  return v < 3;
}

void basics1(){
  std::cout << "Query" << std::endl;
  Query<int,float,true> q;
  q.setMatchFunc(matchF);
  q.setDataFunc(dataF);

  CHECK_INVARIANT(!q.Match(0.0),"");
  CHECK_INVARIANT(q.Match(1.0),"");
  CHECK_INVARIANT(q.Match(1.1),"");
  CHECK_INVARIANT(!q.Match(-2.0),"");
  
  std::cout << "Query2" << std::endl;
  Query<bool,int,true> q2;
  q2.setDataFunc(cmp);
  CHECK_INVARIANT(q2.Match(0),"");
  CHECK_INVARIANT(q2.Match(1),"");
  CHECK_INVARIANT(!q2.Match(3),"");
  CHECK_INVARIANT(!q2.Match(4),"");
  CHECK_INVARIANT(!q2.Match(4.0),"");
}

void basics2(){
  std::cout << "Equality" << std::endl;
  EqualityQuery<int> q2;
  q2.setVal(3);
  CHECK_INVARIANT(!q2.Match(0),"");
  CHECK_INVARIANT(!q2.Match(1),"");
  CHECK_INVARIANT(q2.Match(3),"");
  CHECK_INVARIANT(!q2.Match(-3),"");
  
  std::cout << "Greater" << std::endl;
  GreaterQuery<int> q3;
  q3.setVal(3);
  CHECK_INVARIANT(q3.Match(0),"");
  CHECK_INVARIANT(q3.Match(1),"");
  CHECK_INVARIANT(!q3.Match(3),"");
  CHECK_INVARIANT(!q3.Match(5),"");
  

  std::cout << "GreaterEqual" << std::endl;
  GreaterEqualQuery<int> q4(3);
  CHECK_INVARIANT(q4.Match(0),"");
  CHECK_INVARIANT(q4.Match(1),"");
  CHECK_INVARIANT(q4.Match(3),"");
  CHECK_INVARIANT(!q4.Match(5),"");

  std::cout << "Less" << std::endl;
  LessQuery<int> q5;
  q5.setVal(3);
  CHECK_INVARIANT(!q5.Match(0),"");
  CHECK_INVARIANT(!q5.Match(1),"");
  CHECK_INVARIANT(!q5.Match(3),"");
  CHECK_INVARIANT(q5.Match(5),"");

  std::cout << "LessEqual" << std::endl;
  LessEqualQuery<int> q6(3);

  CHECK_INVARIANT(!q6.Match(0),"");
  CHECK_INVARIANT(!q6.Match(1),"");
  CHECK_INVARIANT(q6.Match(3),"");
  CHECK_INVARIANT(q6.Match(5),"");


  std::cout << "Open Range" << std::endl;
  RangeQuery<int> q7(0,3);
  CHECK_INVARIANT(!q7.Match(0),"");
  CHECK_INVARIANT(q7.Match(1),"");
  CHECK_INVARIANT(!q7.Match(3),"");
  CHECK_INVARIANT(!q7.Match(5),"");


  std::cout << "Closed Range" << std::endl;
  q7.setEndsOpen(false,false);
  CHECK_INVARIANT(q7.Match(0),"");
  CHECK_INVARIANT(q7.Match(1),"");
  CHECK_INVARIANT(q7.Match(3),"");
  CHECK_INVARIANT(!q7.Match(5),"");

}


int convFunc(const char*arg) {return boost::lexical_cast<int>(arg);};

void test7(){
  std::cout << "Set2" << std::endl;
  SetQuery<int,const char*,true> q;
  q.setDataFunc(convFunc);
  q.insert(1);
  q.insert(3);
  q.insert(5);
  
  CHECK_INVARIANT(!q.Match("0"),"");
  CHECK_INVARIANT(q.Match("1"),"");
  CHECK_INVARIANT(q.Match("3"),"");
  CHECK_INVARIANT(!q.Match("-3"),"");

  Query<int,const char*,true> *newQ;
  newQ = q.copy();
  CHECK_INVARIANT(!newQ->Match("0"),"");
  CHECK_INVARIANT(newQ->Match("1"),"");
  CHECK_INVARIANT(newQ->Match("3"),"");
  CHECK_INVARIANT(!newQ->Match("-3"),"");

}

void test8(){
  std::cout << "test serialization" << std::endl;
  {
    Query<int> q;
    q.setMatchFunc(matchF);
  
    CHECK_INVARIANT(!q.Match(0),"");
    CHECK_INVARIANT(q.Match(3),"");

    std::stringstream ss;
    boost::archive::text_oarchive oa(ss);
    oa << (const Query<int>)(q);

    boost::archive::text_iarchive ia(ss);
    Query<int> q2;
    ia >> q2;

    q2.setMatchFunc(matchF);
    CHECK_INVARIANT(!q2.Match(0),"");
    CHECK_INVARIANT(q2.Match(3),"");
  }
  {
    typedef EqualityQuery<int> TestQuery;
    TestQuery q;
    q.setVal(3);
    CHECK_INVARIANT(!q.Match(0),"");
    CHECK_INVARIANT(!q.Match(1),"");
    CHECK_INVARIANT(q.Match(3),"");
    CHECK_INVARIANT(!q.Match(-3),"");

    std::stringstream ss;
    boost::archive::text_oarchive oa(ss);
    oa << (const TestQuery)(q);

    boost::archive::text_iarchive ia(ss);
    TestQuery q2;
    ia >> q2;

    CHECK_INVARIANT(!q2.Match(0),"");
    CHECK_INVARIANT(!q2.Match(1),"");
    CHECK_INVARIANT(q2.Match(3),"");
    CHECK_INVARIANT(!q2.Match(-3),"");
  }
  
  {
    typedef EqualityQuery<double> TestQuery;
    TestQuery q(1.0);
    q.setTol(0.002);
    CHECK_INVARIANT(!q.Match(0.0),"");
    CHECK_INVARIANT(q.Match(1.0),"");
    CHECK_INVARIANT(q.Match(1.001),"");
    CHECK_INVARIANT(!q.Match(1.1),"");
    CHECK_INVARIANT(!q.Match(-2),"");

    std::stringstream ss;
    boost::archive::text_oarchive oa(ss);
    oa << (const TestQuery)(q);

    boost::archive::text_iarchive ia(ss);
    TestQuery q2;
    ia >> q2;
    CHECK_INVARIANT(!q2.Match(0.0),"");
    CHECK_INVARIANT(q2.Match(1.0),"");
    CHECK_INVARIANT(q2.Match(1.001),"");
    CHECK_INVARIANT(!q2.Match(1.1),"");
    CHECK_INVARIANT(!q2.Match(-2),"");
  }

  {
    typedef SetQuery<int,const char*,true> TestQuery;
    TestQuery q;
    q.setDataFunc(convFunc);
    q.insert(1);
    q.insert(3);
    q.insert(5);
  
    CHECK_INVARIANT(!q.Match("0"),"");
    CHECK_INVARIANT(q.Match("1"),"");
    CHECK_INVARIANT(q.Match("3"),"");
    CHECK_INVARIANT(!q.Match("-3"),"");
    std::stringstream ss;
    boost::archive::text_oarchive oa(ss);
    oa << (const TestQuery)(q);

    boost::archive::text_iarchive ia(ss);
    TestQuery q2;
    ia >> q2;

    q2.setDataFunc(convFunc);
    CHECK_INVARIANT(!q2.Match("0"),"");
    CHECK_INVARIANT(q2.Match("1"),"");
    CHECK_INVARIANT(q2.Match("3"),"");
    CHECK_INVARIANT(!q2.Match("-3"),"");
  }

  {
    typedef RangeQuery<int> TestQuery;
    TestQuery q(0,3);
    q.setEndsOpen(false,false);
    CHECK_INVARIANT(q.Match(0),"");
    CHECK_INVARIANT(q.Match(1),"");
    CHECK_INVARIANT(q.Match(3),"");
    CHECK_INVARIANT(!q.Match(5),"");
    std::stringstream ss;
    boost::archive::text_oarchive oa(ss);
    oa << (const TestQuery)(q);

    boost::archive::text_iarchive ia(ss);
    TestQuery q2;
    ia >> q2;
    CHECK_INVARIANT(q2.Match(0),"");
    CHECK_INVARIANT(q2.Match(1),"");
    CHECK_INVARIANT(q2.Match(3),"");
    CHECK_INVARIANT(!q2.Match(5),"");
  }
}

int main()
{
  basics1();
  basics2();

  test1();
  test2();
  test3();
  test4();
  test5();
  test6();
  test7();
  test8();
  return 0;
}
