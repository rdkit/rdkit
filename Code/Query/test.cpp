//
// Copyright (c) 2003-2020 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/test.h>
#include "QueryObjects.h"
#include <iostream>
#include <cmath>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace Queries;

void test1() {
  cout << "Float" << endl;
  EqualityQuery<double> q(1.0);

  TEST_ASSERT(!q.Match(0.0));
  TEST_ASSERT(q.Match(1.0));
  TEST_ASSERT(!q.Match(1.001));
  TEST_ASSERT(!q.Match(1.1));
  TEST_ASSERT(!q.Match(-2));

  cout << "With Tolerance" << endl;
  q.setTol(0.002);
  TEST_ASSERT(!q.Match(0.0));
  TEST_ASSERT(q.Match(1.0));
  TEST_ASSERT(q.Match(1.001));
  TEST_ASSERT(!q.Match(1.1));
  TEST_ASSERT(!q.Match(-2));
  TEST_ASSERT(q.getTypeLabel().empty());
  q.setTypeLabel("FloatEquality");

  Query<double> *newQ;
  newQ = q.copy();
  TEST_ASSERT(!newQ->Match(0.0));
  TEST_ASSERT(newQ->Match(1.0));
  TEST_ASSERT(newQ->Match(1.001));
  TEST_ASSERT(!newQ->Match(1.1));
  TEST_ASSERT(!newQ->Match(-2));
  TEST_ASSERT(newQ->getTypeLabel() == "FloatEquality");
  delete newQ;
}

void test2() {
  cout << "Set" << endl;
  SetQuery<int> q;
  q.insert(1);
  q.insert(3);
  q.insert(5);

  TEST_ASSERT(!q.Match(0));
  TEST_ASSERT(q.Match(1));
  TEST_ASSERT(q.Match(3));
  TEST_ASSERT(!q.Match(-3));

  Query<int> *newQ;
  newQ = q.copy();
  TEST_ASSERT(!newQ->Match(0));
  TEST_ASSERT(newQ->Match(1));
  TEST_ASSERT(newQ->Match(3));
  TEST_ASSERT(!newQ->Match(-3));
  delete newQ;
}

void test3() {
  cout << "And" << endl;
  auto *q = new AndQuery<int>;
  auto *l = new LessQuery<int>;
  l->setVal(0);
  auto *g = new GreaterEqualQuery<int>;
  g->setVal(3);

  q->addChild(Query<int>::CHILD_TYPE(l));
  q->addChild(Query<int>::CHILD_TYPE(g));

  TEST_ASSERT(!q->Match(0));
  TEST_ASSERT(q->Match(1));
  TEST_ASSERT(q->Match(3));
  TEST_ASSERT(!q->Match(-3));

  Query<int> *newQ;
  newQ = q->copy();
  TEST_ASSERT(!newQ->Match(0));
  TEST_ASSERT(newQ->Match(1));
  TEST_ASSERT(newQ->Match(3));
  TEST_ASSERT(!newQ->Match(-3));

  delete newQ;
  delete q;
}

void test4() {
  cout << "Or" << endl;
  auto *q = new OrQuery<int>;
  auto *l = new LessQuery<int>;
  l->setVal(0);
  auto *g = new GreaterEqualQuery<int>;
  g->setVal(3);

  q->addChild(Query<int>::CHILD_TYPE(l));
  q->addChild(Query<int>::CHILD_TYPE(g));

  TEST_ASSERT(q->Match(0));
  TEST_ASSERT(q->Match(1));
  TEST_ASSERT(q->Match(3));
  TEST_ASSERT(q->Match(-3));

  Query<int> *newQ;
  newQ = q->copy();
  TEST_ASSERT(newQ->Match(0));
  TEST_ASSERT(newQ->Match(1));
  TEST_ASSERT(newQ->Match(3));
  TEST_ASSERT(newQ->Match(-3));

  delete newQ;
  delete q;
}

void test5() {
  cout << "XOr" << endl;
  auto *q = new XOrQuery<int>;
  auto *l = new LessQuery<int>;
  l->setVal(0);
  auto *g = new GreaterEqualQuery<int>;
  g->setVal(3);

  q->addChild(Query<int>::CHILD_TYPE(l));
  q->addChild(Query<int>::CHILD_TYPE(g));

  TEST_ASSERT(q->Match(-1));
  TEST_ASSERT(q->Match(0));
  TEST_ASSERT(!q->Match(1));
  TEST_ASSERT(!q->Match(3));
  TEST_ASSERT(q->Match(-3));

  Query<int> *newQ;
  newQ = q->copy();
  TEST_ASSERT(newQ->Match(-1));
  TEST_ASSERT(newQ->Match(0));
  TEST_ASSERT(!newQ->Match(1));
  TEST_ASSERT(!newQ->Match(3));
  TEST_ASSERT(newQ->Match(-3));

  delete newQ;
  delete q;
}

int foofun(double bar) { return int(floor(bar)); };

void test6() {
  cout << "pointer and copy foo" << endl;
  EqualityQuery<int, double, true> q;
  q.setDataFunc(foofun);
  q.setVal(6);
  TEST_ASSERT(q.Match(6.0));
  TEST_ASSERT(q.Match(6.1));
  TEST_ASSERT(!q.Match(5.0));

  Query<int, double, true> *newQ;
  newQ = q.copy();
  TEST_ASSERT(newQ->Match(6.0));
  TEST_ASSERT(newQ->Match(6.1));
  TEST_ASSERT(!newQ->Match(5.0));

  Query<int, double, true> *newQ2 = &q;
  TEST_ASSERT(newQ2->Match(6.0));
  TEST_ASSERT(newQ2->Match(6.1));
  TEST_ASSERT(!newQ2->Match(5.0));

  Query<int, double, true> *newQ3;
  newQ3 = newQ2->copy();
  TEST_ASSERT(newQ3->Match(6.0));
  TEST_ASSERT(newQ3->Match(6.1));
  TEST_ASSERT(!newQ3->Match(5.0));

  delete newQ;
  delete newQ3;
}

bool matchF(int v) { return v == 3; }

int dataF(float v) { return int(floor(v)) * 3; }

bool cmp(int v) { return v < 3; }

void basics1() {
  cout << "Query" << endl;
  Query<int, float, true> q;
  q.setMatchFunc(matchF);
  q.setDataFunc(dataF);

  TEST_ASSERT(!q.Match(0.0));
  TEST_ASSERT(q.Match(1.0));
  TEST_ASSERT(q.Match(1.1));
  TEST_ASSERT(!q.Match(-2.0));

  TEST_ASSERT(!q.getMatchFunc()(0));
  TEST_ASSERT(q.getMatchFunc()(3));
  TEST_ASSERT(q.getDataFunc()(1.0) == dataF(1.0));

  cout << "Query2" << endl;
  Query<bool, int, true> q2;
  q2.setDataFunc(cmp);
  TEST_ASSERT(q2.Match(0));
  TEST_ASSERT(q2.Match(1));
  TEST_ASSERT(!q2.Match(3));
  TEST_ASSERT(!q2.Match(4));
  TEST_ASSERT(!q2.Match(4.0));
}

void basics2() {
  cout << "Equality" << endl;
  EqualityQuery<int> q2;
  q2.setVal(3);
  TEST_ASSERT(!q2.Match(0));
  TEST_ASSERT(!q2.Match(1));
  TEST_ASSERT(q2.Match(3));
  TEST_ASSERT(!q2.Match(-3));

  cout << "Greater" << endl;
  GreaterQuery<int> q3;
  q3.setVal(3);
  TEST_ASSERT(q3.Match(0));
  TEST_ASSERT(q3.Match(1));
  TEST_ASSERT(!q3.Match(3));
  TEST_ASSERT(!q3.Match(5));

  cout << "GreaterEqual" << endl;
  GreaterEqualQuery<int> q4(3);
  TEST_ASSERT(q4.Match(0));
  TEST_ASSERT(q4.Match(1));
  TEST_ASSERT(q4.Match(3));
  TEST_ASSERT(!q4.Match(5));

  cout << "Less" << endl;
  LessQuery<int> q5;
  q5.setVal(3);
  TEST_ASSERT(!q5.Match(0));
  TEST_ASSERT(!q5.Match(1));
  TEST_ASSERT(!q5.Match(3));
  TEST_ASSERT(q5.Match(5));

  cout << "LessEqual" << endl;
  LessEqualQuery<int> q6(3);

  TEST_ASSERT(!q6.Match(0));
  TEST_ASSERT(!q6.Match(1));
  TEST_ASSERT(q6.Match(3));
  TEST_ASSERT(q6.Match(5));

  cout << "Open Range" << endl;
  RangeQuery<int> q7(0, 3);
  TEST_ASSERT(!q7.Match(0));
  TEST_ASSERT(q7.Match(1));
  TEST_ASSERT(!q7.Match(3));
  TEST_ASSERT(!q7.Match(5));

  cout << "Closed Range" << endl;
  q7.setEndsOpen(false, false);
  TEST_ASSERT(q7.Match(0));
  TEST_ASSERT(q7.Match(1));
  TEST_ASSERT(q7.Match(3));
  TEST_ASSERT(!q7.Match(5));
}

int convFunc(const char *arg) { return boost::lexical_cast<int>(arg); };

void test7() {
  cout << "Set2" << endl;
  SetQuery<int, const char *, true> q;
  q.setDataFunc(convFunc);
  q.insert(1);
  q.insert(3);
  q.insert(5);

  TEST_ASSERT(!q.Match("0"));
  TEST_ASSERT(q.Match("1"));
  TEST_ASSERT(q.Match("3"));
  TEST_ASSERT(!q.Match("-3"));

  Query<int, const char *, true> *newQ;
  newQ = q.copy();
  TEST_ASSERT(!newQ->Match("0"));
  TEST_ASSERT(newQ->Match("1"));
  TEST_ASSERT(newQ->Match("3"));
  TEST_ASSERT(!newQ->Match("-3"));

  delete newQ;
}

int main() {
  basics1();
  basics2();

  test1();
  test2();
  test3();
  test4();
  test5();
  test6();
  test7();
}
