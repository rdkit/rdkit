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
using namespace std;
using namespace RDKit;

void testBasics() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Basic Allocations"
                       << std::endl;
  auto *a1 = new Atom(6);
  auto *b1 = new Bond();
  auto *m1 = new ROMol();
  (void)a1;
  (void)b1;
  (void)m1;
  a1 = nullptr;  // intentional leak
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testSMILES() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n SMILES Read" << std::endl;
  string smi = "CCOC";
  ROMol *m = SmilesToMol(smi);
  (void)m;  // leak al the things
  m = SmilesToMol(smi, 0, false);
  smi = "C1COC1";
  RWMol *m2 = SmilesToMol(smi);
  m2 = SmilesToMol(smi, 0, false);
  m2 = SmilesToMol(smi, 0, false);
  unsigned int failed;
  MolOps::sanitizeMol(*m2, failed,
                      MolOps::SANITIZE_CLEANUP | MolOps::SANITIZE_PROPERTIES);
  m2 = SmilesToMol(smi, 0, false);
  MolOps::sanitizeMol(*m2, failed,
                      MolOps::SANITIZE_CLEANUP | MolOps::SANITIZE_PROPERTIES);
  MolOps::symmetrizeSSSR(*m2);
  m2 = SmilesToMol(smi, 0, false);
  MolOps::sanitizeMol(*m2);

  m2 = SmilesToMol(smi, 0, false);
  MolOps::sanitizeMol(*m2);
  MolOps::assignStereochemistry(*m2, true, true, false);

  m2 = SmilesToMol(smi, 0, false);
  MolOps::sanitizeMol(*m2);
  MolOps::assignStereochemistry(*m2, true, true, true);

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testMol() {
  auto *m1 = new RWMol();
  m1->addAtom(new Atom(6), true, true);
  m1->addAtom(new Atom(6), true, true);
  m1->addAtom(new Atom(7), true, true);
  m1->addAtom(new Atom(6), true, true);
  m1->addBond(0, 1, Bond::SINGLE);
  m1->addBond(1, 2, Bond::SINGLE);
  m1->addBond(2, 3, Bond::SINGLE);
  MolOps::sanitizeMol(*m1);
}

void testMols() {
  std::string smis[] = {
      "CN1CCC[C@H]1c2cccnc2",
      "CC1(C)[C@@H](N2[C@@H](CC2=O)S1(=O)=O)C(=O)O",
      "C[C@]1(Cn2ccnn2)[C@@H](N3[C@@H](CC3=O)S1(=O)=O)C(=O)O",
      "CCN(CC)C(=O)[C@@H]1CN(C)[C@H]2Cc3c[nH]c4cccc(C2=C1)c34",
      "CCCN(CCC)[C@H]1CCc2c(O)cccc2C1",
      "CC(=O)NC[C@H]1CN(C(=O)O1)c2ccc(cc2)C(=O)C",
      "CC1(C)Oc2ccc3C=CC(=O)Oc3c2[C@@H](OC(=O)C45CCC(C)(C(=O)O4)C5(C)C)[C@H]"
      "1OC(=O)C67CCC(C)(C(=O)O6)C7(C)C",
      "CCC(C)(C)C(=O)C(=O)N1CCC[C@H]1C(=O)OCCCc2cccnc2",
      "CN1N=C(S/C/1=N/C(=O)C)S(=O)(=O)N",
      "COc1ccc(cc1)[C@@H]2Sc3ccccc3N(CCN(C)C)C(=O)[C@@H]2OC(=O)C",
      "EOS"};
  for (int i = 0; smis[i] != "EOS"; ++i) {
    SmilesToMol(smis[i], 0, false);
  }
  for (int i = 0; smis[i] != "EOS"; ++i) {
    RWMol *m = SmilesToMol(smis[i], 0, false);
    MolOps::sanitizeMol(*m);
  }
}

void testProps() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n properties" << std::endl;
  string smi = "C1COC1";
  RWMol *m = SmilesToMol(smi, 0, false);
  m = SmilesToMol(smi, 0, false);
  for (ROMol::AtomIterator ai = m->beginAtoms(); ai != m->endAtoms(); ++ai) {
    unsigned int v = 1;
    (*ai)->setProp("foo", v, true);
  }
  m = SmilesToMol(smi, 0, false);
  for (ROMol::BondIterator ai = m->beginBonds(); ai != m->endBonds(); ++ai) {
    unsigned int v = 1;
    (*ai)->setProp("foo", v, true);
  }
  m = SmilesToMol(smi, 0, false);
  for (ROMol::AtomIterator ai = m->beginAtoms(); ai != m->endAtoms(); ++ai) {
    unsigned int v = 1;
    (*ai)->setProp("foo", v, true);
    (*ai)->setProp("bar", v, true);
  }
  m = SmilesToMol(smi, 0, false);
  for (ROMol::AtomIterator ai = m->beginAtoms(); ai != m->endAtoms(); ++ai) {
    unsigned int v = 1;
    (*ai)->setProp("foo", v, true);
    (*ai)->setProp("bar", v, true);
    (*ai)->setProp("baz", v, true);
  }
  m = SmilesToMol(smi, 0, false);
  for (ROMol::AtomIterator ai = m->beginAtoms(); ai != m->endAtoms(); ++ai) {
    unsigned int v = 1;
    (*ai)->setProp("foo", v, false);
  }
  m = SmilesToMol(smi, 0, false);
  for (ROMol::AtomIterator ai = m->beginAtoms(); ai != m->endAtoms(); ++ai) {
    unsigned int v = 1;
    (*ai)->setProp("foo", v, false);
    (*ai)->setProp("bar", v, false);
  }
  m = SmilesToMol(smi, 0, false);
  for (ROMol::AtomIterator ai = m->beginAtoms(); ai != m->endAtoms(); ++ai) {
    unsigned int v = 1;
    (*ai)->setProp("foo", v, false);
    (*ai)->setProp("bar", v, false);
    (*ai)->setProp("baz", v, false);
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

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
// boost::logging::enable_logs("rdApp.info");
// test1();  // <- this doesn't seem to actually do anything
#if 1
  // testBasics();
  // testSMILES();
  // testMol();
  // testMols();
  // testProps();
  testDict();

#endif

  return 0;
}
