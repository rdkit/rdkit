//
//  Copyright (C) 2018 Novartis Institutes Of BioMedical Research
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
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <sstream>
#include <iostream>

using namespace std;
using namespace RDKit;

// memory tests for valgrind
void testRemoveAtomBond(RWMol &m, int atomidx, int bondidx)
{
  const Bond *b = m.getBondWithIdx(bondidx);
  if (bondidx>=0) m.removeBond(b->getBeginAtomIdx(), b->getEndAtomIdx());
  if (atomidx>=0) m.removeAtom(atomidx);
}

void testRemovals(RWMol m) {
  for(unsigned int i=0;i<m.getNumAtoms()/2;i++) {
    int atomidx = (i%2==0) ? m.getNumAtoms()-1 : 0;
    int bondidx = (i%2==0) ? m.getNumBonds()-1 : 0;
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
