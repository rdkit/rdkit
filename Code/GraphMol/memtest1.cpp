// $Id$
//
//  Copyright (C) 2014 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/test.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/types.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/Dict.h>

#include <iostream>
using namespace RDKit;

void testDict() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n dict" << std::endl;
  int val = 1;

  auto *d = new Dict();
  delete d;

  d = new Dict();
  d->setVal<int>("foo", val);
  delete d;

  d = new Dict();
  d->setVal<int>("foo", val);
  d->setVal<int>("bar", val);
  delete d;

  d = new Dict();
  d->setVal<int>("foo", val);
  d->setVal<int>("bar", val);
  d->setVal<int>("baz", val);
  delete d;

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}
// -------------------------------------------------------------------
int main() {
  RDLog::InitLogs();
  testDict();

  return 0;
}
