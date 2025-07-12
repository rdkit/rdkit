//
// Copyright (c) 2003-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>
#include "QueryObjects.h"
#include <iostream>
#include <cmath>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace Queries;

int foofun(double bar) { return int(floor(bar)); };

bool matchF(int v) { return v == 3; }

int dataF(float v) { return int(floor(v)) * 3; }

bool cmp(int v) { return v < 3; }

TEST_CASE("basics1_Query") {
  Query<int, float, true> q;
  q.setMatchFunc(matchF);
  q.setDataFunc(dataF);

  REQUIRE(!q.Match(0.0));
  REQUIRE(q.Match(1.0));
  REQUIRE(q.Match(1.1));
  REQUIRE(!q.Match(-2.0));

  REQUIRE(!q.getMatchFunc()(0));
  REQUIRE(q.getMatchFunc()(3));
  REQUIRE(q.getDataFunc()(1.0) == dataF(1.0));
}

TEST_CASE("basics1_Query2") {
  Query<bool, int, true> q2;
  q2.setDataFunc(cmp);

  REQUIRE(q2.Match(0));
  REQUIRE(q2.Match(1));
  REQUIRE(!q2.Match(3));
  REQUIRE(!q2.Match(4));
  REQUIRE(!q2.Match(4.0));
}

TEST_CASE("basics2_Equality") {
  EqualityQuery<int> q2;
  q2.setVal(3);

  REQUIRE(!q2.Match(0));
  REQUIRE(!q2.Match(1));
  REQUIRE(q2.Match(3));
  REQUIRE(!q2.Match(-3));
}

TEST_CASE("basics2_Greater") {
  GreaterQuery<int> q3;
  q3.setVal(3);

  REQUIRE(q3.Match(0));
  REQUIRE(q3.Match(1));
  REQUIRE(!q3.Match(3));
  REQUIRE(!q3.Match(5));
}

TEST_CASE("basics2_GreaterEqual") {
  GreaterEqualQuery<int> q4(3);

  REQUIRE(q4.Match(0));
  REQUIRE(q4.Match(1));
  REQUIRE(q4.Match(3));
  REQUIRE(!q4.Match(5));
}

TEST_CASE("basics2_Less") {
  LessQuery<int> q5;
  q5.setVal(3);

  REQUIRE(!q5.Match(0));
  REQUIRE(!q5.Match(1));
  REQUIRE(!q5.Match(3));
  REQUIRE(q5.Match(5));
}

TEST_CASE("basics2_LessEqual") {
  LessEqualQuery<int> q6(3);

  REQUIRE(!q6.Match(0));
  REQUIRE(!q6.Match(1));
  REQUIRE(q6.Match(3));
  REQUIRE(q6.Match(5));
}

TEST_CASE("basics2_OpenRange") {
  RangeQuery<int> q7(0, 3);

  REQUIRE(!q7.Match(0));
  REQUIRE(q7.Match(1));
  REQUIRE(!q7.Match(3));
  REQUIRE(!q7.Match(5));
}

TEST_CASE("basics2_ClosedRange") {
  RangeQuery<int> q7(0, 3);
  q7.setEndsOpen(false, false);

  REQUIRE(q7.Match(0));
  REQUIRE(q7.Match(1));
  REQUIRE(q7.Match(3));
  REQUIRE(!q7.Match(5));
}

TEST_CASE("testFloatEquality") {
  EqualityQuery<double> q(1.0);

  REQUIRE(!q.Match(0.0));
  REQUIRE(q.Match(1.0));
  REQUIRE(!q.Match(1.001));
  REQUIRE(!q.Match(1.1));
  REQUIRE(!q.Match(-2));

  q.setTol(0.002);
  REQUIRE(!q.Match(0.0));
  REQUIRE(q.Match(1.0));
  REQUIRE(q.Match(1.001));
  REQUIRE(!q.Match(1.1));
  REQUIRE(!q.Match(-2));
  REQUIRE(q.getTypeLabel().empty());
  q.setTypeLabel("FloatEquality");

  Query<double> *newQ;
  newQ = q.copy();
  REQUIRE(!newQ->Match(0.0));
  REQUIRE(newQ->Match(1.0));
  REQUIRE(newQ->Match(1.001));
  REQUIRE(!newQ->Match(1.1));
  REQUIRE(!newQ->Match(-2));
  REQUIRE(newQ->getTypeLabel() == "FloatEquality");
  delete newQ;
}

TEST_CASE("testSetQuery") {
  SetQuery<int> q;
  q.insert(1);
  q.insert(3);
  q.insert(5);

  REQUIRE(!q.Match(0));
  REQUIRE(q.Match(1));
  REQUIRE(q.Match(3));
  REQUIRE(!q.Match(-3));

  Query<int> *newQ;
  newQ = q.copy();
  REQUIRE(!newQ->Match(0));
  REQUIRE(newQ->Match(1));
  REQUIRE(newQ->Match(3));
  REQUIRE(!newQ->Match(-3));
  delete newQ;
}

TEST_CASE("testAndQuery") {
  auto *q = new AndQuery<int>;
  auto *l = new LessQuery<int>;
  l->setVal(0);
  auto *g = new GreaterEqualQuery<int>;
  g->setVal(3);

  q->addChild(Query<int>::CHILD_TYPE(l));
  q->addChild(Query<int>::CHILD_TYPE(g));

  REQUIRE(!q->Match(0));
  REQUIRE(q->Match(1));
  REQUIRE(q->Match(3));
  REQUIRE(!q->Match(-3));

  Query<int> *newQ;
  newQ = q->copy();
  REQUIRE(!newQ->Match(0));
  REQUIRE(newQ->Match(1));
  REQUIRE(newQ->Match(3));
  REQUIRE(!newQ->Match(-3));

  delete newQ;
  delete q;
}

TEST_CASE("testOrQuery") {
  auto *q = new OrQuery<int>;
  auto *l = new LessQuery<int>;
  l->setVal(0);
  auto *g = new GreaterEqualQuery<int>;
  g->setVal(3);

  q->addChild(Query<int>::CHILD_TYPE(l));
  q->addChild(Query<int>::CHILD_TYPE(g));

  REQUIRE(q->Match(0));
  REQUIRE(q->Match(1));
  REQUIRE(q->Match(3));
  REQUIRE(q->Match(-3));

  Query<int> *newQ;
  newQ = q->copy();
  REQUIRE(newQ->Match(0));
  REQUIRE(newQ->Match(1));
  REQUIRE(newQ->Match(3));
  REQUIRE(newQ->Match(-3));

  delete newQ;
  delete q;
}

TEST_CASE("testXOrQuery") {
  auto *q = new XOrQuery<int>;
  auto *l = new LessQuery<int>;
  l->setVal(0);
  auto *g = new GreaterEqualQuery<int>;
  g->setVal(3);

  q->addChild(Query<int>::CHILD_TYPE(l));
  q->addChild(Query<int>::CHILD_TYPE(g));

  REQUIRE(q->Match(-1));
  REQUIRE(q->Match(0));
  REQUIRE(!q->Match(1));
  REQUIRE(!q->Match(3));
  REQUIRE(q->Match(-3));

  Query<int> *newQ;
  newQ = q->copy();
  REQUIRE(newQ->Match(-1));
  REQUIRE(newQ->Match(0));
  REQUIRE(!newQ->Match(1));
  REQUIRE(!newQ->Match(3));
  REQUIRE(newQ->Match(-3));

  delete newQ;
  delete q;
}

TEST_CASE("testPointerAndCopyFoo") {
  EqualityQuery<int, double, true> q;
  q.setDataFunc(foofun);
  q.setVal(6);
  REQUIRE(q.Match(6.0));
  REQUIRE(q.Match(6.1));
  REQUIRE(!q.Match(5.0));

  Query<int, double, true> *newQ;
  newQ = q.copy();
  REQUIRE(newQ->Match(6.0));
  REQUIRE(newQ->Match(6.1));
  REQUIRE(!newQ->Match(5.0));

  Query<int, double, true> *newQ2 = &q;
  REQUIRE(newQ2->Match(6.0));
  REQUIRE(newQ2->Match(6.1));
  REQUIRE(!newQ2->Match(5.0));

  Query<int, double, true> *newQ3;
  newQ3 = newQ2->copy();
  REQUIRE(newQ3->Match(6.0));
  REQUIRE(newQ3->Match(6.1));
  REQUIRE(!newQ3->Match(5.0));

  delete newQ;
  delete newQ3;
}

int convFunc(const char *arg) { return boost::lexical_cast<int>(arg); };

TEST_CASE("testSetQueryWithDataFunc") {
  SetQuery<int, const char *, true> q;
  q.setDataFunc(convFunc);
  q.insert(1);
  q.insert(3);
  q.insert(5);

  REQUIRE(!q.Match("0"));
  REQUIRE(q.Match("1"));
  REQUIRE(q.Match("3"));
  REQUIRE(!q.Match("-3"));

  Query<int, const char *, true> *newQ;
  newQ = q.copy();
  REQUIRE(!newQ->Match("0"));
  REQUIRE(newQ->Match("1"));
  REQUIRE(newQ->Match("3"));
  REQUIRE(!newQ->Match("-3"));

  delete newQ;
}
