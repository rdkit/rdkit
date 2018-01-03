// $Id$
//
//  Copyright (C) 2001-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/RDKitQueries.h>
#include <RDGeneral/types.h>
#include <RDGeneral/RDLog.h>
//#include <boost/log/functions.hpp>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <sstream>
#include <iostream>
#include <boost/random/uniform_int.hpp>

using namespace std;
using namespace RDKit;

// memory tests for valgrind
void testRemoveAtomBond(RWMol &m, int atomidx, int bondidx)
{
  if (atomidx>=0) m.removeAtom(atomidx);
  if (bondidx>=0) m.removeAtom(bondidx);
}

void testRemovals(RWMol m) {
  boost::random::mt19937 rng;   
  for(int i=0;i<m.getNumAtoms()/2;++i) {
    boost::random::uniform_int_distribution<> arng(0,m.getNumAtoms()-1);
    boost::random::uniform_int_distribution<> brng(0,m.getNumBonds()-1);
    int atomidx = arng(rng);
    int bondidx = brng(rng);
    testRemoveAtomBond(m, atomidx, bondidx);
  }
}

void test1() {
  std::string smi = "CCC";
  for(int i=4; i<20; ++i) {
    smi += "C";
    RWMol *m = SmilesToMol(smi.c_str());
    testRemovals(*m);
    delete m;
  }
}

// -------------------------------------------------------------------
int main() {
  RDLog::InitLogs();
// boost::logging::enable_logs("rdApp.info");
  test1();

  return 0;
}
