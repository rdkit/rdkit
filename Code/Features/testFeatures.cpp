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
#include <RDGeneral/test.h>
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

void testParser() {
  std::cerr << "-------------------------------------" << std::endl;
  std::cerr << "parser testing." << std::endl;

  std::string text;
  int n;
  bool ok;

  text = "p1";
  ok = parse(text.c_str(), (ch_p('p') >> int_p[assign_a(n)]), space_p).full;
  TEST_ASSERT(ok);
  TEST_ASSERT(n == 1);

  text = "p12";
  ok = parse(text.c_str(), (ch_p('p') >> int_p[assign_a(n)]), space_p).full;
  TEST_ASSERT(ok);
  TEST_ASSERT(n == 12);

  text = "p2 + p3";
  ok = parse(text.c_str(), (ch_p('p') >> int_p[assign_a(n)]), space_p).full;
  TEST_ASSERT(!ok);

  std::cerr << "  done" << std::endl;
}
int main() {
  testParser();
}
