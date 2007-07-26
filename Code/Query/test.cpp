// $Id$
//
// Copyright (c) 2003-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#include <RDGeneral/RDLog.h>
#include "QueryObjects.h"
#include <iostream>
#include <math.h>
#include <boost/lexical_cast.hpp>
#include <sstream>

using namespace std;
using namespace Queries;



// we need to use the helper because the serialization fns want
// const references to things and because we need to register derived classes:
template <typename T, typename U, bool F>
void serialHelper(std::ostream &dest,Queries::Query<T,U,F> const *q){
  boost::archive::text_oarchive oa(dest);  

  oa.register_type(static_cast< Queries::Query<T,U,F> *>(NULL));
  oa.register_type(static_cast< Queries::EqualityQuery<T,U,F> *>(NULL));
  oa.register_type(static_cast< Queries::GreaterEqualQuery<T,U,F> *>(NULL));
  oa.register_type(static_cast< Queries::GreaterQuery<T,U,F> *>(NULL));
  oa.register_type(static_cast< Queries::LessEqualQuery<T,U,F> *>(NULL));
  oa.register_type(static_cast< Queries::LessQuery<T,U,F> *>(NULL));
  oa.register_type(static_cast< Queries::RangeQuery<T,U,F> *>(NULL));
  oa.register_type(static_cast< Queries::SetQuery<T,U,F> *>(NULL));
  oa.register_type(static_cast< Queries::AndQuery<T,U,F> *>(NULL));
  oa.register_type(static_cast< Queries::OrQuery<T,U,F> *>(NULL));
  oa.register_type(static_cast< Queries::XOrQuery<T,U,F> *>(NULL));

  oa << BOOST_SERIALIZATION_NVP(q);
}

// !! !! !! !! !! !! !! !! !!
// This doesn't work
template <typename T, typename U, bool F>
void deserialHelper(std::stringstream &src,Queries::Query<T,U,F> *q){
  boost::archive::text_iarchive ia(src);  

  ia.register_type(static_cast< Queries::Query<T,U,F> *>(NULL));
  ia.register_type(static_cast< Queries::EqualityQuery<T,U,F> *>(NULL));
  ia.register_type(static_cast< Queries::GreaterEqualQuery<T,U,F> *>(NULL));
  ia.register_type(static_cast< Queries::GreaterQuery<T,U,F> *>(NULL));
  ia.register_type(static_cast< Queries::LessEqualQuery<T,U,F> *>(NULL));
  ia.register_type(static_cast< Queries::LessQuery<T,U,F> *>(NULL));
  ia.register_type(static_cast< Queries::RangeQuery<T,U,F> *>(NULL));
  ia.register_type(static_cast< Queries::SetQuery<T,U,F> *>(NULL));
  ia.register_type(static_cast< Queries::AndQuery<T,U,F> *>(NULL));
  ia.register_type(static_cast< Queries::OrQuery<T,U,F> *>(NULL));
  ia.register_type(static_cast< Queries::XOrQuery<T,U,F> *>(NULL));
  

  ia >> *q;
}
// !! !! !! !! !! !! !! !! !!


#define DESERIALIZE_QUERY(_a_,_b_,T,U,F) {\
    boost::archive::text_iarchive ia((_a_));\
    ia.register_type(static_cast< Queries::Query<T,U,F> *>(NULL));\
    ia.register_type(static_cast< Queries::EqualityQuery<T,U,F> *>(NULL));\
    ia.register_type(static_cast< Queries::GreaterEqualQuery<T,U,F> *>(NULL));\
    ia.register_type(static_cast< Queries::GreaterQuery<T,U,F> *>(NULL));\
    ia.register_type(static_cast< Queries::LessEqualQuery<T,U,F> *>(NULL));\
    ia.register_type(static_cast< Queries::LessQuery<T,U,F> *>(NULL));\
    ia.register_type(static_cast< Queries::RangeQuery<T,U,F> *>(NULL));\
    ia.register_type(static_cast< Queries::SetQuery<T,U,F> *>(NULL));\
    ia.register_type(static_cast< Queries::AndQuery<T,U,F> *>(NULL));\
    ia.register_type(static_cast< Queries::OrQuery<T,U,F> *>(NULL));\
    ia.register_type(static_cast< Queries::XOrQuery<T,U,F> *>(NULL));\
    ia >> (_b_);\
  }

void testSerialize(){
  BOOST_LOG(rdInfoLog) << "Serialize" << std::endl;
  
  {
    Query<double,double,0> *p1=new Query<double,double,0>();

    p1->setDescription("foo");
    p1->setNegation(1);
    
    std::stringstream opkl;
    serialHelper(opkl,p1);
  
    delete p1;
    p1 = new Query<double,double,0>();
    p1->setDescription("bar");
    DESERIALIZE_QUERY(opkl,p1,double,double,0);
  
    TEST_ASSERT(p1->getDescription()=="foo");
    TEST_ASSERT(p1->getNegation());

  }
  {
    EqualityQuery<double> *p1=new EqualityQuery<double>(1.0);

    CHECK_INVARIANT(!p1->Match(0.0),"");
    CHECK_INVARIANT(p1->Match(1.0),"");
    std::stringstream opkl;
    serialHelper(opkl,p1);
  
    EqualityQuery<double> *p2 = new EqualityQuery<double>(2.0);
    
    std::cout << "1" << std::endl;
    DESERIALIZE_QUERY(opkl,p2,double,double,0);
    std::cout << "2" << std::endl;

    TEST_ASSERT(p2->getVal()==1.0);
    TEST_ASSERT(!p2->Match(2.0));
    TEST_ASSERT(p2->Match(1.0));
  }

  {
    GreaterEqualQuery<int> *p1 = new GreaterEqualQuery<int>(3);
    CHECK_INVARIANT(p1->Match(0),"");
    CHECK_INVARIANT(p1->Match(1),"");
    CHECK_INVARIANT(p1->Match(3),"");
    CHECK_INVARIANT(!p1->Match(5),"");
    std::stringstream pkl;
    serialHelper(pkl,p1);
    GreaterEqualQuery<double> *p12=new GreaterEqualQuery<double>(2.0);
    DESERIALIZE_QUERY(pkl,p12,double,double,0)
    CHECK_INVARIANT(p12->Match(0),"");
    CHECK_INVARIANT(p12->Match(1),"");
    CHECK_INVARIANT(p12->Match(3),"");
    CHECK_INVARIANT(!p12->Match(5),"");
  }

  {  
    BOOST_LOG(rdInfoLog) << "Greater" << std::endl;
    GreaterQuery<int> *p3=new GreaterQuery<int>();
    p3->setVal(3);
    CHECK_INVARIANT(p3->Match(0),"");
    CHECK_INVARIANT(p3->Match(1),"");
    CHECK_INVARIANT(!p3->Match(3),"");
    CHECK_INVARIANT(!p3->Match(5),"");

    std::stringstream pkl;
    serialHelper(pkl,p3);
    GreaterQuery<int> *p32=new GreaterQuery<int>();
    DESERIALIZE_QUERY(pkl,p32,int,int,0)
    
    CHECK_INVARIANT(p32->Match(0),"");
    CHECK_INVARIANT(p32->Match(1),"");
    CHECK_INVARIANT(!p32->Match(3),"");
    CHECK_INVARIANT(!p32->Match(5),"");
  } 
  
  { 
    GreaterEqualQuery<int> *p4 = new GreaterEqualQuery<int>(3);
    CHECK_INVARIANT(p4->Match(0),"");
    CHECK_INVARIANT(p4->Match(1),"");
    CHECK_INVARIANT(p4->Match(3),"");
    CHECK_INVARIANT(!p4->Match(5),"");

    GreaterEqualQuery<int> *p42 = new GreaterEqualQuery<int>(5);
    std::stringstream pkl;
    serialHelper(pkl,p4);
    DESERIALIZE_QUERY(pkl,p42,int,int,0)
    
    CHECK_INVARIANT(p42->Match(0),"");
    CHECK_INVARIANT(p42->Match(1),"");
    CHECK_INVARIANT(p42->Match(3),"");
    CHECK_INVARIANT(!p42->Match(5),"");
  }

  {  
    LessQuery<int> *q5=new LessQuery<int>();
    q5->setVal(3);
    CHECK_INVARIANT(!q5->Match(0),"");
    CHECK_INVARIANT(!q5->Match(1),"");
    CHECK_INVARIANT(!q5->Match(3),"");
    CHECK_INVARIANT(q5->Match(5),"");

    LessQuery<int> *q52=new LessQuery<int>();
    std::stringstream pkl;
    serialHelper(pkl,q5);
    DESERIALIZE_QUERY(pkl,q52,int,int,0)

    CHECK_INVARIANT(!q52->Match(0),"");
    CHECK_INVARIANT(!q52->Match(1),"");
    CHECK_INVARIANT(!q52->Match(3),"");
    CHECK_INVARIANT(q52->Match(5),"");
  }
  {
    LessEqualQuery<int> *q6=new LessEqualQuery<int>(3);
  
    CHECK_INVARIANT(!q6->Match(0),"");
    CHECK_INVARIANT(!q6->Match(1),"");
    CHECK_INVARIANT(q6->Match(3),"");
    CHECK_INVARIANT(q6->Match(5),"");

    LessEqualQuery<int> *q62=new LessEqualQuery<int>(-1);
    std::stringstream pkl;
    serialHelper(pkl,q6);
    DESERIALIZE_QUERY(pkl,q62,int,int,0)

    CHECK_INVARIANT(!q62->Match(0),"");
    CHECK_INVARIANT(!q62->Match(1),"");
    CHECK_INVARIANT(q62->Match(3),"");
    CHECK_INVARIANT(q62->Match(5),"");
  }
  {
    RangeQuery<int> *q7=new RangeQuery<int>(0,3);
    CHECK_INVARIANT(!q7->Match(0),"");
    CHECK_INVARIANT(q7->Match(1),"");
    CHECK_INVARIANT(!q7->Match(3),"");
    CHECK_INVARIANT(!q7->Match(5),"");

    RangeQuery<int> *q72= new RangeQuery<int>(9,12);
    std::stringstream pkl;
    serialHelper(pkl,q7);
    DESERIALIZE_QUERY(pkl,q72,int,int,0)

    CHECK_INVARIANT(!q72->Match(0),"");
    CHECK_INVARIANT(q72->Match(1),"");
    CHECK_INVARIANT(!q72->Match(3),"");
    CHECK_INVARIANT(!q72->Match(5),"");
  }

  {
    RangeQuery<int> *q7=new RangeQuery<int>(0,3);
    q7->setEndsOpen(false,false);
    CHECK_INVARIANT(q7->Match(0),"");
    CHECK_INVARIANT(q7->Match(1),"");
    CHECK_INVARIANT(q7->Match(3),"");
    CHECK_INVARIANT(!q7->Match(5),"");

    RangeQuery<int> *q72=new RangeQuery<int>(9,12);
    std::stringstream pkl;
    serialHelper(pkl,q7);
    DESERIALIZE_QUERY(pkl,q72,int,int,0)

    CHECK_INVARIANT(q72->Match(0),"");
    CHECK_INVARIANT(q72->Match(1),"");
    CHECK_INVARIANT(q72->Match(3),"");
    CHECK_INVARIANT(!q72->Match(5),"");
  }

  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
}



void test1(){
  BOOST_LOG(rdInfoLog) << "Float" << std::endl;
  EqualityQuery<double> q(1.0);

  CHECK_INVARIANT(!q.Match(0.0),"");
  CHECK_INVARIANT(q.Match(1.0),"");
  CHECK_INVARIANT(!q.Match(1.001),"");
  CHECK_INVARIANT(!q.Match(1.1),"");
  CHECK_INVARIANT(!q.Match(-2),"");


  BOOST_LOG(rdInfoLog) << "With Tolerance" << std::endl;
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
  BOOST_LOG(rdInfoLog) << "Set" << std::endl;
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

  {
    Query<int> *q2=new SetQuery<int>();
    
    std::stringstream pkl;
    serialHelper(pkl,newQ);
    DESERIALIZE_QUERY(pkl,q2,int,int,0)

    CHECK_INVARIANT(!q2->Match(0),"");
    CHECK_INVARIANT(q2->Match(1),"");
    CHECK_INVARIANT(q2->Match(3),"");
    CHECK_INVARIANT(!q2->Match(-3),"");
  }

}

void test3(){
  BOOST_LOG(rdInfoLog) << "And" << std::endl;
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
  
  {
    Query<int> *q2=new AndQuery<int>();
    
    std::stringstream pkl;
    serialHelper(pkl,newQ);
    DESERIALIZE_QUERY(pkl,q2,int,int,0)

    CHECK_INVARIANT(!q2->Match(0),"");
    CHECK_INVARIANT(q2->Match(1),"");
    CHECK_INVARIANT(q2->Match(3),"");
    CHECK_INVARIANT(!q2->Match(-3),"");
  }


}

void test4(){
  BOOST_LOG(rdInfoLog) << "Or" << std::endl;
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

  {
    Query<int> *q2=new OrQuery<int>();
    
    std::stringstream pkl;
    serialHelper(pkl,newQ);
    DESERIALIZE_QUERY(pkl,q2,int,int,0)

    CHECK_INVARIANT(q2->Match(0),"");
    CHECK_INVARIANT(q2->Match(1),"");
    CHECK_INVARIANT(q2->Match(3),"");
    CHECK_INVARIANT(q2->Match(-3),"");
  }

}

void test5(){
  BOOST_LOG(rdInfoLog) << "XOr" << std::endl;
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

  {
    Query<int> *q2=new XOrQuery<int>();
    
    std::stringstream pkl;
    serialHelper(pkl,newQ);
    DESERIALIZE_QUERY(pkl,q2,int,int,0)
    CHECK_INVARIANT(q2->Match(-1),"");
    CHECK_INVARIANT(q2->Match(0),"");
    CHECK_INVARIANT(!q2->Match(1),"");
    CHECK_INVARIANT(!q2->Match(3),"");
    CHECK_INVARIANT(q2->Match(-3),"");
  }


}

int foofun(double bar) { return int(floor(bar)); };

void test6(){
  BOOST_LOG(rdInfoLog) << "pointer and copy foo" << std::endl;
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
  BOOST_LOG(rdInfoLog) << "Query" << std::endl;
  Query<int,float,true> q;
  q.setMatchFunc(matchF);
  q.setDataFunc(dataF);

  CHECK_INVARIANT(!q.Match(0.0),"");
  CHECK_INVARIANT(q.Match(1.0),"");
  CHECK_INVARIANT(q.Match(1.1),"");
  CHECK_INVARIANT(!q.Match(-2.0),"");
  
  BOOST_LOG(rdInfoLog) << "Query2" << std::endl;
  Query<bool,int,true> q2;
  q2.setDataFunc(cmp);
  CHECK_INVARIANT(q2.Match(0),"");
  CHECK_INVARIANT(q2.Match(1),"");
  CHECK_INVARIANT(!q2.Match(3),"");
  CHECK_INVARIANT(!q2.Match(4),"");

}

void basics2(){
  BOOST_LOG(rdInfoLog) << "Equality" << std::endl;
  EqualityQuery<int> q2;
  q2.setVal(3);
  CHECK_INVARIANT(!q2.Match(0),"");
  CHECK_INVARIANT(!q2.Match(1),"");
  CHECK_INVARIANT(q2.Match(3),"");
  CHECK_INVARIANT(!q2.Match(-3),"");
  
  BOOST_LOG(rdInfoLog) << "Greater" << std::endl;
  GreaterQuery<int> q3;
  q3.setVal(3);
  CHECK_INVARIANT(q3.Match(0),"");
  CHECK_INVARIANT(q3.Match(1),"");
  CHECK_INVARIANT(!q3.Match(3),"");
  CHECK_INVARIANT(!q3.Match(5),"");
  

  BOOST_LOG(rdInfoLog) << "GreaterEqual" << std::endl;
  GreaterEqualQuery<int> q4(3);
  CHECK_INVARIANT(q4.Match(0),"");
  CHECK_INVARIANT(q4.Match(1),"");
  CHECK_INVARIANT(q4.Match(3),"");
  CHECK_INVARIANT(!q4.Match(5),"");

  BOOST_LOG(rdInfoLog) << "Less" << std::endl;
  LessQuery<int> q5;
  q5.setVal(3);
  CHECK_INVARIANT(!q5.Match(0),"");
  CHECK_INVARIANT(!q5.Match(1),"");
  CHECK_INVARIANT(!q5.Match(3),"");
  CHECK_INVARIANT(q5.Match(5),"");

  BOOST_LOG(rdInfoLog) << "LessEqual" << std::endl;
  LessEqualQuery<int> q6(3);

  CHECK_INVARIANT(!q6.Match(0),"");
  CHECK_INVARIANT(!q6.Match(1),"");
  CHECK_INVARIANT(q6.Match(3),"");
  CHECK_INVARIANT(q6.Match(5),"");


  BOOST_LOG(rdInfoLog) << "Open Range" << std::endl;
  RangeQuery<int> q7(0,3);
  CHECK_INVARIANT(!q7.Match(0),"");
  CHECK_INVARIANT(q7.Match(1),"");
  CHECK_INVARIANT(!q7.Match(3),"");
  CHECK_INVARIANT(!q7.Match(5),"");


  BOOST_LOG(rdInfoLog) << "Closed Range" << std::endl;
  q7.setEndsOpen(false,false);
  CHECK_INVARIANT(q7.Match(0),"");
  CHECK_INVARIANT(q7.Match(1),"");
  CHECK_INVARIANT(q7.Match(3),"");
  CHECK_INVARIANT(!q7.Match(5),"");

}


int convFunc(const char*arg) {return boost::lexical_cast<int>(arg);};

void test7(){
  BOOST_LOG(rdInfoLog) << "Set2" << std::endl;
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


int main()
{
  RDLog::InitLogs();
  basics1();
  basics2();
  testSerialize();
#if 1
  test1();
  test2();
  test3();
  test4();
  test5();
  test6();
  test7();
#endif
  return 0;
}
