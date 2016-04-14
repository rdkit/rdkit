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
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <ForceField/ForceField.h>
#include <Geometry/Trajectory.h>
#include <sstream>
#include <iostream>

using namespace std;
using namespace RDKit;

// -------------------------------------------------------------------
void testBookmarks(ROMol m) {
  int i;

  // ------------------------
  // simple bookmark stuff
  Atom *a1 = m.getAtomWithIdx(1);
  m.setAtomBookmark(a1, 666);
  Atom *a2 = m.getAtomWithBookmark(666);

  TEST_ASSERT(a2->getIdx() == a1->getIdx());

  bool ok;
  m.clearAtomBookmark(666);

  boost::logging::disable_logs("rdApp.error");
  try {
    a2 = m.getAtomWithBookmark(666);
    ok = 0;
  } catch (...) {
    ok = 1;
  }
  boost::logging::enable_logs("rdApp.error");
  CHECK_INVARIANT(ok, "atom bookmark not properly cleared");

  // ------------------------
  // repeat a bookmark
  a1 = m.getAtomWithIdx(1);
  CHECK_INVARIANT(a1->getIdx() == 1, "");
  m.setAtomBookmark(a1, 666);
  m.setAtomBookmark(m.getAtomWithIdx(0), 666);
  a2 = m.getAtomWithBookmark(666);
  CHECK_INVARIANT(a2->getIdx() == 1, "");
  CHECK_INVARIANT(a2->getIdx() == a1->getIdx(), "");
  m.clearAtomBookmark(666, a2);
  a2 = m.getAtomWithBookmark(666);
  i = a2->getIdx();
  CHECK_INVARIANT(i == 0, "");
  m.clearAtomBookmark(666, a2);
  boost::logging::disable_logs("rdApp.error");
  try {
    a2 = m.getAtomWithBookmark(666);
    ok = 0;
  } catch (...) {
    ok = 1;
  }
  boost::logging::enable_logs("rdApp.error");
  CHECK_INVARIANT(ok, "atom bookmark not properly cleared");

  // make sure clearAtomBookmark doesn't barf if there's no
  // such bookmark:
  m.clearAtomBookmark(777);

  //----------------------------
  // now do bond bookmarks
  Bond *b1 = m.getBondWithIdx(0);
  m.setBondBookmark(b1, 23);
  Bond *b2 = m.getBondWithBookmark(23);
  CHECK_INVARIANT(b2->getIdx() == b1->getIdx(), "");

  m.clearBondBookmark(23);
  boost::logging::disable_logs("rdApp.error");
  try {
    b2 = m.getBondWithBookmark(23);
    ok = 0;
  } catch (...) {
    ok = 1;
  }
  boost::logging::enable_logs("rdApp.error");
  CHECK_INVARIANT(ok, "bond bookmark not properly cleared");

  m.setBondBookmark(b1, 23);
  m.setBondBookmark(m.getBondWithIdx(1), 23);
  b2 = m.getBondWithBookmark(23);
  CHECK_INVARIANT(b2->getIdx() == b1->getIdx(), "");
  m.clearBondBookmark(23, b2);
  b2 = m.getBondWithBookmark(23);
  CHECK_INVARIANT(b2->getIdx() == 1, "");
  m.clearBondBookmark(23, b2);
  boost::logging::disable_logs("rdApp.error");
  try {
    b2 = m.getBondWithBookmark(23);
    ok = 0;
  } catch (...) {
    ok = 1;
  }
  boost::logging::enable_logs("rdApp.error");
  CHECK_INVARIANT(ok, "bond bookmark not properly cleared");

  // make sure clearAtomBookmark doesn't barf if there's no
  // such bookmark:
  m.clearBondBookmark(777);
}

void testMolProps() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing Mol Property Caches" << std::endl;
  RWMol m2;
  STR_VECT propNames;

  m2.addAtom(new Atom(6));
  m2.addAtom(new Atom(6));
  m2.addBond(0, 1, Bond::TRIPLE);

  CHECK_INVARIANT(!m2.hasProp("prop1"), "");
  CHECK_INVARIANT(!m2.hasProp("prop2"), "");
  m2.setProp("prop1", 2);
  int tmpI;
  std::string tmpS;
  CHECK_INVARIANT(m2.hasProp("prop1"), "");
  m2.getProp("prop1", tmpI);
  CHECK_INVARIANT(tmpI == 2, "");
  m2.getProp("prop1", tmpS);
  CHECK_INVARIANT(tmpS == "2", "");
  m2.setProp("prop1", std::string("2"));
  CHECK_INVARIANT(m2.hasProp("prop1"), "");
  m2.getProp("prop1", tmpS);
  CHECK_INVARIANT(tmpS == "2", "");
  std::string tmpString("2");
  m2.setProp("prop1", tmpString.c_str());
  CHECK_INVARIANT(m2.hasProp("prop1"), "");
  m2.getProp("prop1", tmpS);
  CHECK_INVARIANT(tmpS == "2", "");

  tmpS = "name";
  m2.setProp(common_properties::_Name, tmpS);

  propNames = m2.getPropList(false, false);
  TEST_ASSERT(propNames.size() == 1);
  propNames = m2.getPropList(true, false);
  TEST_ASSERT(propNames.size() == 2);

  // check for computed properties
  m2.setProp("cprop1", 1, true);
  m2.setProp("cprop2", 2, true);
  STR_VECT cplst;
  m2.getProp(detail::computedPropName, cplst);
  CHECK_INVARIANT(cplst.size() == 2, "");
  CHECK_INVARIANT(cplst[0] == "cprop1", "");
  CHECK_INVARIANT(cplst[1] == "cprop2", "");

  propNames = m2.getPropList(false, false);
  TEST_ASSERT(propNames.size() == 1);
  propNames = m2.getPropList(true, false);
  TEST_ASSERT(propNames.size() == 2);
  propNames = m2.getPropList(false, true);
  TEST_ASSERT(propNames.size() == 3);
  propNames = m2.getPropList(true, true);
  TEST_ASSERT(propNames.size() == 5);
  propNames = m2.getPropList();
  TEST_ASSERT(propNames.size() == 5);

  m2.clearProp("cprop1");
  CHECK_INVARIANT(!m2.hasProp("cprop1"), "");
  m2.getProp(detail::computedPropName, cplst);
  CHECK_INVARIANT(cplst.size() == 1, "");

  m2.clearComputedProps();
  CHECK_INVARIANT(!m2.hasProp("cprop2"), "");
  m2.getProp(detail::computedPropName, cplst);
  CHECK_INVARIANT(cplst.size() == 0, "");

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testClearMol() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing RWMol.clear()"
                       << std::endl;
  RWMol m2;

  m2.addAtom(new Atom(6));
  m2.addAtom(new Atom(6));
  m2.addBond(0, 1, Bond::TRIPLE);

  TEST_ASSERT(!m2.hasProp("prop1"));
  m2.setProp("prop1", 2);
  int tmpI;
  TEST_ASSERT(m2.hasProp("prop1"));
  m2.getProp("prop1", tmpI);
  TEST_ASSERT(tmpI == 2);

  TEST_ASSERT(m2.hasProp(detail::computedPropName));

  m2.clear();
  TEST_ASSERT(!m2.hasProp("prop1"));
  TEST_ASSERT(m2.getNumAtoms() == 0);
  TEST_ASSERT(m2.getNumBonds() == 0);
  TEST_ASSERT(m2.getAtomBookmarks()->empty());
  TEST_ASSERT(m2.getBondBookmarks()->empty());

  TEST_ASSERT(m2.hasProp(detail::computedPropName));  // <- github issue 176
  TEST_ASSERT(m2.getPropList().size() == 1);

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testAtomProps() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing Atom Property Caches" << std::endl;
  RWMol m2;
  m2.addAtom(new Atom(6));
  m2.addAtom(new Atom(6));
  m2.addBond(0, 1, Bond::TRIPLE);

  Atom *a1 = m2.getAtomWithIdx(0);
  Atom *a2 = m2.getAtomWithIdx(0);
  Atom *a3 = &(*a1);
  CHECK_INVARIANT(!a1->hasProp("prop1"), "");
  CHECK_INVARIANT(!a1->hasProp("prop2"), "");
  CHECK_INVARIANT(!a2->hasProp("prop1"), "");
  CHECK_INVARIANT(!a2->hasProp("prop2"), "");
  CHECK_INVARIANT(!a3->hasProp("prop1"), "");
  CHECK_INVARIANT(!a3->hasProp("prop2"), "");
  a1->setProp("prop1", 3);
  a1->setProp("prop2", 4);
  CHECK_INVARIANT(a1->hasProp("prop1"), "");
  CHECK_INVARIANT(a1->hasProp("prop2"), "");
  CHECK_INVARIANT(a2->hasProp("prop1"), "");
  CHECK_INVARIANT(a2->hasProp("prop2"), "");
  CHECK_INVARIANT(a3->hasProp("prop1"), "");
  CHECK_INVARIANT(a3->hasProp("prop2"), "");
  CHECK_INVARIANT(!a1->hasProp("bogus"), "");
  CHECK_INVARIANT(!a2->hasProp("bogus"), "");
  CHECK_INVARIANT(!a3->hasProp("bogus"), "");

  bool ok = false;
  a1->setProp<double>("dprop", 4);
  TEST_ASSERT(a1->hasProp("dprop"));
  try {
    a1->getProp<int>("dprop");
  } catch (const boost::bad_any_cast &e) {
    ok = true;
  }
  TEST_ASSERT(ok);
  a1->setProp<int>("iprop", 4);
  TEST_ASSERT(a1->hasProp("iprop"));
  ok = false;
  try {
    a1->getProp<double>("iprop");
  } catch (const boost::bad_any_cast &e) {
    ok = true;
  }
  TEST_ASSERT(ok);

  int tmp;
  a1->getProp("prop1", tmp);
  CHECK_INVARIANT(tmp == 3, "");
  a1->getProp("prop2", tmp);
  CHECK_INVARIANT(tmp == 4, "");
  a2->getProp("prop1", tmp);
  CHECK_INVARIANT(tmp == 3, "");
  a2->getProp("prop2", tmp);
  CHECK_INVARIANT(tmp == 4, "");
  a3->getProp("prop1", tmp);
  CHECK_INVARIANT(tmp == 3, "");
  a3->getProp("prop2", tmp);
  CHECK_INVARIANT(tmp == 4, "");

  // check for computed properties
  a1->setProp("cprop1", 1, true);
  a1->setProp("cprop2", 2, true);
  STR_VECT cplst;
  a1->getProp(detail::computedPropName, cplst);
  CHECK_INVARIANT(cplst.size() == 2, "");
  CHECK_INVARIANT(cplst[0] == "cprop1", "");
  CHECK_INVARIANT(cplst[1] == "cprop2", "");

  a1->clearProp("cprop1");
  CHECK_INVARIANT(!a1->hasProp("cprop1"), "");
  a1->getProp(detail::computedPropName, cplst);
  CHECK_INVARIANT(cplst.size() == 1, "");

  a1->clearComputedProps();
  CHECK_INVARIANT(!a1->hasProp("cprop2"), "");
  a1->getProp(detail::computedPropName, cplst);
  CHECK_INVARIANT(cplst.size() == 0, "");

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testBondProps() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing Bond Property Caches" << std::endl;
  RWMol m2;
  m2.addAtom(new Atom(6));
  m2.addAtom(new Atom(6));
  m2.addBond(0, 1, Bond::TRIPLE);

  Bond *b1 = m2.getBondWithIdx(0);
  Bond *b2 = m2.getBondWithIdx(0);
  CHECK_INVARIANT(!b1->hasProp("prop1"), "");
  CHECK_INVARIANT(!b1->hasProp("prop2"), "");
  CHECK_INVARIANT(!b2->hasProp("prop1"), "");
  CHECK_INVARIANT(!b2->hasProp("prop2"), "");
  b1->setProp("prop1", 3);
  b1->setProp("prop2", 4);
  CHECK_INVARIANT(b1->hasProp("prop1"), "");
  CHECK_INVARIANT(b1->hasProp("prop2"), "");
  CHECK_INVARIANT(b2->hasProp("prop1"), "");
  CHECK_INVARIANT(b2->hasProp("prop2"), "");
  CHECK_INVARIANT(!b1->hasProp("bogus"), "");
  CHECK_INVARIANT(!b2->hasProp("bogus"), "");

  int tmp;
  b1->getProp("prop1", tmp);
  CHECK_INVARIANT(tmp == 3, "");
  b1->getProp("prop2", tmp);
  CHECK_INVARIANT(tmp == 4, "");
  b2->getProp("prop1", tmp);
  CHECK_INVARIANT(tmp == 3, "");
  b2->getProp("prop2", tmp);
  CHECK_INVARIANT(tmp == 4, "");

  // check for computed properties
  b1->setProp("cprop1", 1, true);
  b1->setProp("cprop2", 2, true);
  STR_VECT cplst;
  b1->getProp(detail::computedPropName, cplst);
  CHECK_INVARIANT(cplst.size() == 2, "");
  CHECK_INVARIANT(cplst[0] == "cprop1", "");
  CHECK_INVARIANT(cplst[1] == "cprop2", "");

  b1->clearProp("cprop1");
  CHECK_INVARIANT(!b1->hasProp("cprop1"), "");
  b1->getProp(detail::computedPropName, cplst);
  CHECK_INVARIANT(cplst.size() == 1, "");

  b1->clearComputedProps();
  CHECK_INVARIANT(!b1->hasProp("cprop2"), "");
  b1->getProp(detail::computedPropName, cplst);
  CHECK_INVARIANT(cplst.size() == 0, "");

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

// this is here because there was at one time a problem with crashes when doing
// this stuff
void testPropLeak() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing Atom and Bond Property Caches"
      << std::endl;
  RWMol m2;
  m2.addAtom(new Atom(6));
  m2.addAtom(new Atom(6));
  m2.addBond(0, 1, Bond::TRIPLE);

  Atom *a1 = m2.getAtomWithIdx(0);
  Atom *a2 = m2.getAtomWithIdx(0);
  CHECK_INVARIANT(!a1->hasProp("prop1"), "");
  CHECK_INVARIANT(!a1->hasProp("prop2"), "");
  CHECK_INVARIANT(!a2->hasProp("prop1"), "");
  CHECK_INVARIANT(!a2->hasProp("prop2"), "");
  a1->setProp("prop1", 3);
  a1->setProp("prop2", 4);
  CHECK_INVARIANT(a1->hasProp("prop1"), "");
  CHECK_INVARIANT(a1->hasProp("prop2"), "");
  CHECK_INVARIANT(a2->hasProp("prop1"), "");
  CHECK_INVARIANT(a2->hasProp("prop2"), "");
  CHECK_INVARIANT(!a1->hasProp("bogus"), "");
  CHECK_INVARIANT(!a2->hasProp("bogus"), "");

  int tmp;
  a1->getProp("prop1", tmp);
  CHECK_INVARIANT(tmp == 3, "");
  a1->getProp("prop2", tmp);
  CHECK_INVARIANT(tmp == 4, "");
  a2->getProp("prop1", tmp);
  CHECK_INVARIANT(tmp == 3, "");
  a2->getProp("prop2", tmp);
  CHECK_INVARIANT(tmp == 4, "");

  Bond *b1 = m2.getBondWithIdx(0);
  Bond *b2 = m2.getBondWithIdx(0);
  CHECK_INVARIANT(!b1->hasProp("prop1"), "");
  CHECK_INVARIANT(!b1->hasProp("prop2"), "");
  CHECK_INVARIANT(!b2->hasProp("prop1"), "");
  CHECK_INVARIANT(!b2->hasProp("prop2"), "");
  b1->setProp("prop1", 3);
  b1->setProp("prop2", 4);
  CHECK_INVARIANT(b1->hasProp("prop1"), "");
  CHECK_INVARIANT(b1->hasProp("prop2"), "");
  CHECK_INVARIANT(b2->hasProp("prop1"), "");
  CHECK_INVARIANT(b2->hasProp("prop2"), "");
  CHECK_INVARIANT(!b1->hasProp("bogus"), "");
  CHECK_INVARIANT(!b2->hasProp("bogus"), "");

  b1->getProp("prop1", tmp);
  CHECK_INVARIANT(tmp == 3, "");
  b1->getProp("prop2", tmp);
  CHECK_INVARIANT(tmp == 4, "");
  b2->getProp("prop1", tmp);
  CHECK_INVARIANT(tmp == 3, "");
  b2->getProp("prop2", tmp);
  CHECK_INVARIANT(tmp == 4, "");

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testMisc() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Misc Properties"
                       << std::endl;
  RWMol m2;
  m2.addAtom(new Atom(6));
  m2.addAtom(new Atom(6));
  m2.addAtom(new Atom(6));
  m2.addAtom(new Atom(6));
  m2.addBond(0, 1, Bond::SINGLE);
  m2.addBond(1, 2, Bond::SINGLE);
  m2.addBond(0, 2, Bond::SINGLE);
  m2.addBond(2, 3, Bond::SINGLE);

  MolOps::sanitizeMol(m2);

  Bond *bnd;
  bnd = m2.getBondBetweenAtoms(0, 1);
  CHECK_INVARIANT(bnd, "");
  bnd = m2.getBondBetweenAtoms(1, 0);
  CHECK_INVARIANT(bnd, "");
  bnd = m2.getBondBetweenAtoms(3, 0);
  CHECK_INVARIANT(!bnd, "");
  bnd = m2.getBondBetweenAtoms(0, 3);
  CHECK_INVARIANT(!bnd, "");
  const Bond *cbnd;
  cbnd = m2.getBondBetweenAtoms(0, 1);
  CHECK_INVARIANT(cbnd, "");
  cbnd = m2.getBondBetweenAtoms(1, 0);
  CHECK_INVARIANT(cbnd, "");
  cbnd = m2.getBondBetweenAtoms(0, 3);
  CHECK_INVARIANT(!cbnd, "");
  cbnd = m2.getBondBetweenAtoms(3, 0);
  CHECK_INVARIANT(!cbnd, "");

  CHECK_INVARIANT(m2.getAtomWithIdx(0)->getTotalNumHs() == 2, "");

  // we'll check atom deletion and handling of bookmarks on deletion
  // simultaneously:
  //  (The bookmark thing was the root of Issue 96)
  m2.setAtomBookmark(m2.getAtomWithIdx(0), 2342);
  m2.setBondBookmark(m2.getBondWithIdx(0), 2343);
  m2.removeAtom(static_cast<unsigned int>(0));
  CHECK_INVARIANT(!m2.hasAtomBookmark(2342), "");
  CHECK_INVARIANT(!m2.hasBondBookmark(2343), "");
  CHECK_INVARIANT(m2.getNumAtoms() == 3, "");
  CHECK_INVARIANT(m2.getNumBonds() == 2, "");
  MolOps::sanitizeMol(m2);
  CHECK_INVARIANT(m2.getAtomWithIdx(0)->getTotalNumHs() == 3, "");

  m2.addAtom(new Atom(1));
  m2.addBond(2, 3, Bond::SINGLE);
  MolOps::sanitizeMol(m2);

  CHECK_INVARIANT(m2.getAtomWithIdx(0)->getTotalNumHs() == 3, "");
  CHECK_INVARIANT(m2.getAtomWithIdx(0)->getTotalNumHs(true) == 3, "");
  CHECK_INVARIANT(m2.getAtomWithIdx(2)->getTotalNumHs() == 2, "");
  CHECK_INVARIANT(m2.getAtomWithIdx(2)->getTotalNumHs(true) == 3, "");

  Atom *other = m2.getBondWithIdx(1)->getOtherAtom(m2.getAtomWithIdx(1));
  CHECK_INVARIANT(other, "");

  const Atom *at = m2.getAtomWithIdx(1);
  ROMol::OEDGE_ITER begin, end;
  boost::tie(begin, end) = m2.getAtomBonds(at);
  while (begin != end) {
    const Atom *at2 = m2[*begin]->getOtherAtom(at);
    TEST_ASSERT(at2);
    begin++;
  }

  ROMol::VERTEX_ITER atBegin, atEnd;
  boost::tie(atBegin, atEnd) = m2.getVertices();
  TEST_ASSERT(atBegin != atEnd);
  while (atBegin != atEnd) {
    const ATOM_SPTR at2 = m2[*atBegin];
    TEST_ASSERT(at2->getIdx() == *atBegin);
    atBegin++;
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testDegree() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing degree operations"
                       << std::endl;
  RWMol *m;

  m = new RWMol();
  m->addAtom(new Atom(6));
  m->addAtom(new Atom(6));
  m->addAtom(new Atom(6));
  m->addAtom(new Atom(6));
  m->addBond(0, 1, Bond::SINGLE);
  m->addBond(1, 2, Bond::SINGLE);
  m->addBond(0, 2, Bond::SINGLE);
  m->addBond(2, 3, Bond::SINGLE);

  MolOps::sanitizeMol(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->getDegree() == 2);
  TEST_ASSERT(m->getAtomWithIdx(0)->getTotalNumHs() == 2);
  TEST_ASSERT(m->getAtomWithIdx(0)->getTotalDegree() == 4);
  TEST_ASSERT(m->getAtomWithIdx(2)->getDegree() == 3);
  TEST_ASSERT(m->getAtomWithIdx(2)->getTotalNumHs() == 1);
  TEST_ASSERT(m->getAtomWithIdx(2)->getTotalDegree() == 4);

  delete m;
  m = new RWMol();
  m->addAtom(new Atom(6));
  m->addAtom(new Atom(6));
  m->addAtom(new Atom(6));
  m->addAtom(new Atom(6));
  m->addAtom(new Atom(1));
  m->addBond(0, 1, Bond::SINGLE);
  m->addBond(1, 2, Bond::SINGLE);
  m->addBond(0, 2, Bond::SINGLE);
  m->addBond(2, 3, Bond::SINGLE);
  m->addBond(0, 4, Bond::SINGLE);

  MolOps::sanitizeMol(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->getDegree() == 3);
  TEST_ASSERT(m->getAtomWithIdx(0)->getTotalNumHs() == 1);
  TEST_ASSERT(m->getAtomWithIdx(0)->getTotalNumHs(true) == 2);
  TEST_ASSERT(m->getAtomWithIdx(0)->getTotalDegree() == 4);
  TEST_ASSERT(m->getAtomWithIdx(2)->getDegree() == 3);
  TEST_ASSERT(m->getAtomWithIdx(2)->getTotalNumHs() == 1);
  TEST_ASSERT(m->getAtomWithIdx(2)->getTotalDegree() == 4);

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testIssue1993296() {
  RWMol *m = new RWMol();
  bool ok;
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 1993296" << std::endl;

  m->addAtom(new Atom(6));
  m->addAtom(new Atom(6));
  m->addBond(0, 1, Bond::SINGLE);
  ok = false;
  try {
    m->addBond(0, 1, Bond::SINGLE);
  } catch (...) {
    ok = true;
  }
  TEST_ASSERT(ok);
  ok = false;
  try {
    m->addBond(1, 0, Bond::SINGLE);
  } catch (...) {
    ok = true;
  }
  TEST_ASSERT(ok);

  // not technically part of 1993296, but related: we also throw
  // on adding self bonds
  ok = false;
  try {
    m->addBond(1, 1, Bond::SINGLE);
  } catch (...) {
    ok = true;
  }

  Bond *newB = new Bond();
  newB->setBeginAtomIdx(0);
  newB->setEndAtomIdx(1);
  newB->setBondType(Bond::SINGLE);
  ok = false;
  try {
    m->addBond(newB);
  } catch (...) {
    ok = true;
  }
  TEST_ASSERT(ok);

  // not technically part of 1993296, but related: we also throw
  // on adding self bonds
  newB->setBeginAtomIdx(0);
  newB->setEndAtomIdx(0);
  ok = false;
  try {
    m->addBond(newB);
  } catch (...) {
    ok = true;
  }
  TEST_ASSERT(ok);
  delete newB;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIssue2381580() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 2381580" << std::endl;

  {
    RWMol *m = new RWMol();
    m->addAtom(new Atom(5));
    m->addAtom(new Atom(6));
    m->addAtom(new Atom(6));
    m->addAtom(new Atom(6));
    m->addBond(0, 1, Bond::SINGLE);
    m->addBond(0, 2, Bond::SINGLE);
    m->addBond(0, 3, Bond::SINGLE);
    MolOps::sanitizeMol(*m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getFormalCharge() == 0);
    TEST_ASSERT(m->getAtomWithIdx(0)->getExplicitValence() == 3);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumImplicitHs() == 0);
    delete m;
  }

  {
    RWMol *m = new RWMol();
    m->addAtom(new Atom(5));
    m->addAtom(new Atom(6));
    m->addAtom(new Atom(6));
    m->addAtom(new Atom(6));
    m->addAtom(new Atom(6));
    m->addBond(0, 1, Bond::SINGLE);
    m->addBond(0, 2, Bond::SINGLE);
    m->addBond(0, 3, Bond::SINGLE);
    m->addBond(0, 4, Bond::SINGLE);
    m->getAtomWithIdx(0)->setFormalCharge(-1);
    MolOps::sanitizeMol(*m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getFormalCharge() == -1);
    TEST_ASSERT(m->getAtomWithIdx(0)->getExplicitValence() == 4);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumImplicitHs() == 0);
    delete m;
  }

  {
    RWMol *m = new RWMol();
    m->addAtom(new Atom(5));
    m->addAtom(new Atom(6));
    m->addAtom(new Atom(6));
    m->addAtom(new Atom(6));
    m->addAtom(new Atom(6));
    m->addBond(0, 1, Bond::SINGLE);
    m->addBond(0, 2, Bond::SINGLE);
    m->addBond(0, 3, Bond::SINGLE);
    m->addBond(0, 4, Bond::SINGLE);
    bool ok = false;
    try {
      MolOps::sanitizeMol(*m);
    } catch (MolSanitizeException &e) {
      ok = true;
    }
    TEST_ASSERT(ok);
    delete m;
  }

  {
    RWMol *m = new RWMol();
    m->addAtom(new Atom(5));
    m->addAtom(new Atom(6));
    m->addAtom(new Atom(6));
    m->addAtom(new Atom(6));
    m->addAtom(new Atom(6));
    m->addBond(0, 1, Bond::SINGLE);
    m->addBond(0, 2, Bond::SINGLE);
    m->addBond(0, 3, Bond::SINGLE);
    m->addBond(0, 4, Bond::SINGLE);
    m->getAtomWithIdx(0)->setFormalCharge(+1);
    bool ok = false;
    try {
      MolOps::sanitizeMol(*m);
    } catch (MolSanitizeException &e) {
      ok = true;
    }
    TEST_ASSERT(ok);
    delete m;
  }

  {
    RWMol *m = new RWMol();
    m->addAtom(new Atom(5));
    m->addAtom(new Atom(6));
    m->addAtom(new Atom(6));
    m->addBond(0, 1, Bond::SINGLE);
    m->addBond(0, 2, Bond::SINGLE);
    m->getAtomWithIdx(0)->setFormalCharge(+1);
    MolOps::sanitizeMol(*m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getFormalCharge() == 1);
    TEST_ASSERT(m->getAtomWithIdx(0)->getExplicitValence() == 2);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumImplicitHs() == 0);
    delete m;
  }

  {
    RWMol *m = new RWMol();
    m->addAtom(new Atom(5));
    m->addAtom(new Atom(6));
    m->addAtom(new Atom(6));
    m->addAtom(new Atom(6));
    m->addBond(0, 1, Bond::SINGLE);
    m->addBond(0, 2, Bond::SINGLE);
    m->addBond(0, 3, Bond::SINGLE);
    m->getAtomWithIdx(0)->setFormalCharge(-1);
    MolOps::sanitizeMol(*m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getFormalCharge() == -1);
    TEST_ASSERT(m->getAtomWithIdx(0)->getExplicitValence() == 3);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumImplicitHs() == 1);
    TEST_ASSERT(m->getAtomWithIdx(0)->getExplicitValence() +
                    m->getAtomWithIdx(0)->getImplicitValence() ==
                rdcast<int>(m->getAtomWithIdx(0)->getTotalValence()));
    TEST_ASSERT(m->getAtomWithIdx(1)->getExplicitValence() +
                    m->getAtomWithIdx(1)->getImplicitValence() ==
                rdcast<int>(m->getAtomWithIdx(1)->getTotalValence()));
    TEST_ASSERT(m->getAtomWithIdx(2)->getExplicitValence() +
                    m->getAtomWithIdx(2)->getImplicitValence() ==
                rdcast<int>(m->getAtomWithIdx(2)->getTotalValence()));
    TEST_ASSERT(m->getAtomWithIdx(3)->getExplicitValence() +
                    m->getAtomWithIdx(3)->getImplicitValence() ==
                rdcast<int>(m->getAtomWithIdx(3)->getTotalValence()));
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIssue2840217() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 2840217" << std::endl;

  {
    RWMol *m = new RWMol();
    for (unsigned int i = 0; i < 200; ++i) {
      m->addAtom(new Atom(6));
      m->addAtom(new Atom(6));
      m->addAtom(new Atom(6));
      m->addAtom(new Atom(6));
      m->addAtom(new Atom(6));
      m->addAtom(new Atom(6));
      for (unsigned int j = 0; j < 5; ++j) {
        m->addBond(i * 6 + j, i * 6 + j + 1, Bond::AROMATIC);
      }
      m->addBond(i * 6, i * 6 + 5, Bond::AROMATIC);
    }
    TEST_ASSERT(m->getNumBonds() == 1200);
    MolOps::sanitizeMol(*m);
    TEST_ASSERT(m->getNumAtoms() == 1200);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test1() {
  {
    RWMol m;
    Atom *newAtom = new Atom(8);

    m.addAtom(newAtom);
    CHECK_INVARIANT(m.getAtomWithIdx(0)->getIdx() == 0, "");
    newAtom = new Atom(6);
    m.addAtom(newAtom);
    CHECK_INVARIANT(m.getAtomWithIdx(0)->getIdx() == 0, "");
    CHECK_INVARIANT(m.getAtomWithIdx(1)->getIdx() == 1, "");

    newAtom = new Atom(7);
    m.addAtom(newAtom);
    CHECK_INVARIANT(m.getAtomWithIdx(0)->getIdx() == 0, "");
    CHECK_INVARIANT(m.getAtomWithIdx(1)->getIdx() == 1, "");
    CHECK_INVARIANT(m.getAtomWithIdx(2)->getIdx() == 2, "");
    CHECK_INVARIANT(
        m.getAtomWithIdx(1)->getOwningMol().getAtomWithIdx(1)->getIdx() == 1,
        "");

    m.addBond(0, 1, Bond::SINGLE);
    m.addBond(1, 2, Bond::DOUBLE);
    CHECK_INVARIANT(m.getBondWithIdx(0)->getIdx() == 0, "");
    CHECK_INVARIANT(m.getBondWithIdx(1)->getIdx() == 1, "");

    CHECK_INVARIANT(m.getBondWithIdx(0)->getBondType() == Bond::SINGLE, "");
    CHECK_INVARIANT(m.getBondWithIdx(1)->getBondType() == Bond::DOUBLE, "");

    CHECK_INVARIANT(m.getBondWithIdx(0)->getBeginAtom()->getIdx() == 0, "");
    CHECK_INVARIANT(m.getBondWithIdx(0)->getEndAtom()->getIdx() == 1, "");
    CHECK_INVARIANT(m.getBondWithIdx(1)->getBeginAtom()->getIdx() == 1, "");
    CHECK_INVARIANT(m.getBondWithIdx(1)->getEndAtom()->getIdx() == 2, "");

    testBookmarks(m);

    // Using operator<< on a non-sanitized molecule is a test of Issue156:
    ROMol::ADJ_ITER ai1, ai2;
    boost::tie(ai1, ai2) = m.getAtomNeighbors(m.getAtomWithIdx(1));
    m.updatePropertyCache();
    boost::logging::disable_logs("rdApp.info");
    while (ai1 != ai2) {
      BOOST_LOG(rdInfoLog) << *m.getAtomWithIdx(*ai1) << endl;
      ai1++;
    }

    m.addAtom(new Atom(6));
    Bond *bsp = m.createPartialBond(2);
    m.setBondBookmark(bsp, 47);
    m.finishPartialBond(3, 47, Bond::SINGLE);
    m.clearBondBookmark(47);
    BOOST_LOG(rdInfoLog) << "partial bond added:" << endl;
    unsigned int i;
    m.updatePropertyCache();
    for (i = 0; i < m.getNumAtoms(); i++) {
      Atom *a = m.getAtomWithIdx(i);
      BOOST_LOG(rdInfoLog) << "\t" << *a << endl;
    }

    int newAtNum = m.addAtom(new Atom(6));
    m.addBond(0, newAtNum, Bond::SINGLE);

    BOOST_LOG(rdInfoLog) << "Again:" << endl;
    m.updatePropertyCache();
    for (i = 0; i < m.getNumAtoms(); i++) {
      Atom *a = m.getAtomWithIdx(i);
      BOOST_LOG(rdInfoLog) << "\t" << *a << endl;
    }

    RWMol m2;
    m2.addAtom(new Atom(6));
    m2.addAtom(new Atom(6));
    // QueryAtom *qA = new QueryAtom;
    // qA->setAtomicNum(7);
    // m2.addAtom(qA);
    m2.addAtom(new QueryAtom(7));
    m2.addBond(0, 1, Bond::TRIPLE);
    m2.addBond(1, 2, Bond::SINGLE);

    m.insertMol(m2);
    m.updatePropertyCache();
    BOOST_LOG(rdInfoLog) << "post-insert:" << endl;
    for (i = 0; i < m.getNumAtoms(); i++) {
      Atom *a = m.getAtomWithIdx(i);
      BOOST_LOG(rdInfoLog) << "\t" << *a << endl;
    }

    BOOST_LOG(rdInfoLog) << " ------------------- " << endl;
    Atom *newA = new Atom(12);
    int newIdx = m.addAtom(newA);
    m.addBond(newIdx - 1, newIdx, Bond::AROMATIC);
    // m.debugMol(cout);
    BOOST_LOG(rdInfoLog) << " trying a replace " << endl;
    Atom *repA = new Atom(22);
    m.replaceAtom(newIdx, repA);
  }
  {
    RWMol m;
    m.addAtom(new Atom(6));
    m.addAtom(new Atom(6));
    m.addBond(0, 1, Bond::SINGLE);
    Conformer *conf = new Conformer(m.getNumAtoms());
    m.addConformer(conf);
    m.getConformer().setAtomPos(0, RDGeom::Point3D(1.0, 0.0, 0.0));
    m.getConformer().setAtomPos(1, RDGeom::Point3D(0.0, 1.0, 0.0));

    RWMol m2;
    // insert molecule without a conf:
    m2.addAtom(new Atom(6));
    m.insertMol(m2);
    TEST_ASSERT(m.getConformer().getNumAtoms() == m.getNumAtoms());
    TEST_ASSERT(feq(m.getConformer().getAtomPos(2).x, 0.0));
    TEST_ASSERT(feq(m.getConformer().getAtomPos(2).y, 0.0));
    TEST_ASSERT(feq(m.getConformer().getAtomPos(2).z, 0.0));

    // insert molecule with a conf:
    conf = new Conformer(m2.getNumAtoms());
    m2.addConformer(conf);
    m2.getConformer().setAtomPos(0, RDGeom::Point3D(1.0, 1.0, 0.0));
    m.insertMol(m2);
    TEST_ASSERT(m.getConformer().getNumAtoms() == m.getNumAtoms());
    TEST_ASSERT(feq(m.getConformer().getAtomPos(2).x, 0.0));
    TEST_ASSERT(feq(m.getConformer().getAtomPos(2).y, 0.0));
    TEST_ASSERT(feq(m.getConformer().getAtomPos(2).z, 0.0));
    TEST_ASSERT(feq(m.getConformer().getAtomPos(3).x, 1.0));
    TEST_ASSERT(feq(m.getConformer().getAtomPos(3).y, 1.0));
    TEST_ASSERT(feq(m.getConformer().getAtomPos(3).z, 0.0));
  }
  {
    // start with a molecule with no conf
    RWMol m;
    m.addAtom(new Atom(6));
    m.addAtom(new Atom(6));
    m.addBond(0, 1, Bond::SINGLE);
    TEST_ASSERT(m.getNumConformers() == 0);

    RWMol m2;
    // insert molecule without a conf:
    m2.addAtom(new Atom(6));
    m.insertMol(m2);
    TEST_ASSERT(m.getNumConformers() == 0);

    // insert molecule with a conf:
    Conformer *conf = new Conformer(m2.getNumAtoms());
    m2.addConformer(conf);
    m2.getConformer().setAtomPos(0, RDGeom::Point3D(1.0, 1.0, 0.0));
    m.insertMol(m2);
    TEST_ASSERT(m.getNumConformers() == 1);
    TEST_ASSERT(m.getConformer().getNumAtoms() == m.getNumAtoms());
    TEST_ASSERT(feq(m.getConformer().getAtomPos(0).x, 0.0));
    TEST_ASSERT(feq(m.getConformer().getAtomPos(0).y, 0.0));
    TEST_ASSERT(feq(m.getConformer().getAtomPos(0).z, 0.0));
    TEST_ASSERT(feq(m.getConformer().getAtomPos(3).x, 1.0));
    TEST_ASSERT(feq(m.getConformer().getAtomPos(3).y, 1.0));
    TEST_ASSERT(feq(m.getConformer().getAtomPos(3).z, 0.0));
  }
}

void testAddConformersFromTrajectory() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n";
  BOOST_LOG(rdInfoLog) << "Testing adding conformers from a trajectory" << std::endl;
  std::string molBlock =
    "\n"
    "     RDKit          3D\n"
    "\n"
    " 71 74  0  0  0  0  0  0  0  0999 V2000\n"
    "    8.2543    3.1901   -0.3005 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    7.4558    1.9712    0.0938 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    7.3934    1.0441   -0.9483 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    6.6660   -0.0533   -0.4641 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    5.1928    0.2346   -0.4609 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    4.3713   -0.9410   -0.5770 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    3.1852   -1.0034   -1.2291 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    2.2914    0.1276   -1.6316 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.9308   -0.4468   -1.9908 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.1417   -0.7821   -0.7545 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.1848    0.3695    0.0456 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -1.5661    0.7686   -0.0745 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -2.4768   -0.0640    0.8206 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -3.8874    0.1143    0.3941 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -4.6333   -0.9984    0.0264 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -6.0127   -0.9516   -0.0400 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -6.7062    0.1599    0.3963 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -8.0408    0.4828   -0.1977 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -7.7914    1.1180   -1.5591 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -8.7622    1.4403    0.7265 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -8.8409   -0.7397   -0.4395 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -8.9121   -1.6637    0.4258 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -9.7414   -0.7636   -1.5059 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -5.9736    1.2357    0.8565 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -4.5843    1.2252    0.8530 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.6263    1.4884   -0.3942 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    2.0541    1.0258   -0.4230 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    2.9225   -2.3317   -1.2963 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    3.6061   -2.9745   -0.3180 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    3.3554   -4.1536    0.3735 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    3.7653   -4.2712    1.6948 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    4.8254   -3.4613    2.0796 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    5.1978   -2.3436    1.3419 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    4.5694   -2.0799    0.1305 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    9.3138    3.1372    0.0031 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    7.8117    4.0754    0.1798 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    8.2358    3.3535   -1.4074 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    6.4027    2.2146    0.3634 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    7.9270    1.5444    1.0040 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    7.0677   -0.2415    0.5615 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    6.9530   -0.9105   -1.1025 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    4.9578    0.7259    0.5137 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    4.9985    0.9430   -1.3033 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    2.7171    0.7264   -2.4494 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.3994    0.2339   -2.6810 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.1342   -1.4171   -2.5076 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.7632   -1.3370   -1.0391 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.7845   -1.4394   -0.1311 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.0125    0.1989    1.0673 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -1.6672    1.8215    0.2925 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -1.8705    0.7271   -1.1337 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -2.3045    0.3159    1.8590 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -2.1980   -1.1367    0.7635 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -4.1513   -1.9468   -0.2114 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -6.6138   -1.7460   -0.4718 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -7.0727    0.4399   -2.0858 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -7.3144    2.1076   -1.4482 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -8.7609    1.1720   -2.1135 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -8.3137    2.4504    0.5729 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -8.6170    1.0817    1.7580 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -9.8244    1.4444    0.4200 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -6.4629    2.0541    1.3719 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -4.0445    2.0563    1.3058 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.3329    1.8224   -1.3991 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.4920    2.3164    0.3160 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    2.2025    0.3766    0.4766 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    2.7945    1.8369   -0.3969 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    2.4404   -4.6964    0.1303 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    3.3157   -5.0055    2.3587 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    5.4272   -3.7654    2.9380 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    5.5668   -1.5069    1.9380 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  2  1  0\n"
    "  2  3  1  0\n"
    "  3  4  1  0\n"
    "  4  5  1  0\n"
    "  5  6  1  0\n"
    "  6  7  1  0\n"
    "  7  8  1  0\n"
    "  8  9  1  0\n"
    "  9 10  1  0\n"
    " 10 11  1  0\n"
    " 11 12  1  0\n"
    " 12 13  1  0\n"
    " 13 14  1  0\n"
    " 14 15  2  0\n"
    " 15 16  1  0\n"
    " 16 17  2  0\n"
    " 17 18  1  0\n"
    " 18 19  1  0\n"
    " 18 20  1  0\n"
    " 18 21  1  0\n"
    " 21 22  2  0\n"
    " 21 23  1  0\n"
    " 17 24  1  0\n"
    " 24 25  2  0\n"
    " 11 26  1  0\n"
    " 26 27  1  0\n"
    "  7 28  2  0\n"
    " 28 29  1  0\n"
    " 29 30  2  0\n"
    " 30 31  1  0\n"
    " 31 32  2  0\n"
    " 32 33  1  0\n"
    " 33 34  2  0\n"
    " 34  6  1  0\n"
    " 27  8  1  0\n"
    " 34 29  1  0\n"
    " 25 14  1  0\n"
    "  1 35  1  0\n"
    "  1 36  1  0\n"
    "  1 37  1  0\n"
    "  2 38  1  0\n"
    "  2 39  1  0\n"
    "  4 40  1  0\n"
    "  4 41  1  0\n"
    "  5 42  1  0\n"
    "  5 43  1  0\n"
    "  8 44  1  0\n"
    "  9 45  1  0\n"
    "  9 46  1  0\n"
    " 10 47  1  0\n"
    " 10 48  1  0\n"
    " 11 49  1  0\n"
    " 12 50  1  0\n"
    " 12 51  1  0\n"
    " 13 52  1  0\n"
    " 13 53  1  0\n"
    " 15 54  1  0\n"
    " 16 55  1  0\n"
    " 19 56  1  0\n"
    " 19 57  1  0\n"
    " 19 58  1  0\n"
    " 20 59  1  0\n"
    " 20 60  1  0\n"
    " 20 61  1  0\n"
    " 24 62  1  0\n"
    " 25 63  1  0\n"
    " 26 64  1  0\n"
    " 26 65  1  0\n"
    " 27 66  1  0\n"
    " 27 67  1  0\n"
    " 30 68  1  0\n"
    " 31 69  1  0\n"
    " 32 70  1  0\n"
    " 33 71  1  0\n"
    "M  CHG  2  11   1  23  -1\n"
    "M  END\n";
  RWMol *mol = MolBlockToMol(molBlock, true, false);
  const unsigned int everySteps = 10;
  const unsigned int maxIts = 1000;
  double gradTol = 0.01;
  std::string rdbase = getenv("RDBASE");
  std::string fName = rdbase + "/Code/GraphMol/test_data/bilastine_trajectory.sdf";
  SDWriter w(fName);
  ForceFields::ForceField *field = MMFF::constructForceField(*mol);
  field->initialize();
  RDGeom::Trajectory traj(3, mol->getNumAtoms());
  traj.setFreePosOnDestroy(true);
  field->minimize(everySteps, &traj, maxIts, gradTol);
  mol->removeConformer(0);
  mol->addConformersFromTrajectory(&traj);
  for (unsigned int nConf = 0;
    nConf < mol->getNumConformers(); ++nConf) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(4)
      << traj.getSnapshot(nConf).getEnergy();
    mol->setProp("ENERGY", ss.str(), false);
    w.write(*mol, nConf);
  }
  w.close();
  delete field;
  delete mol;
}

void testPeriodicTable() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n";
  BOOST_LOG(rdInfoLog) << "Testing properties from periodic table" << std::endl;

  TEST_ASSERT(PeriodicTable::getTable()->getDefaultValence(6) == 4);
  TEST_ASSERT(PeriodicTable::getTable()->getNouterElecs(6) == 4);
  TEST_ASSERT(PeriodicTable::getTable()->getMostCommonIsotope(6) == 12);
  TEST_ASSERT(PeriodicTable::getTable()->getMostCommonIsotopeMass(6) == 12.0);
  TEST_ASSERT(PeriodicTable::getTable()->getMostCommonIsotopeMass(6) == 12.0);
  TEST_ASSERT(feq(PeriodicTable::getTable()->getMostCommonIsotopeMass(4),
                  9.0122, 1e-4));
  TEST_ASSERT(feq(PeriodicTable::getTable()->getRb0(6), 0.77, 1e-2));

  TEST_ASSERT(PeriodicTable::getTable()->getDefaultValence(26) == -1);
  TEST_ASSERT(PeriodicTable::getTable()->getDefaultValence(57) == -1);

  // this was sf.net issue 269
  int anum;
  anum = PeriodicTable::getTable()->getAtomicNumber("C");
  TEST_ASSERT(anum == 6);
  try {
    anum = PeriodicTable::getTable()->getAtomicNumber("Xx");
  } catch (...) {
    anum = -1;
  }
  TEST_ASSERT(anum == -1);

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testAddAtomWithConf() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n";
  BOOST_LOG(rdInfoLog) << "Testing issue 264: adding atoms to molecules that "
                          "already have conformers" << std::endl;
  {
    RWMol m;

    m.addAtom(new Atom(6));
    m.addAtom(new Atom(6));

    Conformer *conf = new Conformer(m.getNumAtoms());
    m.addConformer(conf);

    m.addAtom(new Atom(6));
    TEST_ASSERT(m.getConformer().getNumAtoms() == m.getNumAtoms());
  }
  {
    RWMol m;

    m.addAtom(new Atom(6));
    m.addAtom(new Atom(6));

    Conformer *conf = new Conformer(m.getNumAtoms());
    m.addConformer(conf);

    m.addAtom();
    TEST_ASSERT(m.getConformer().getNumAtoms() == m.getNumAtoms());
  }
  {  // make sure things are ok even if there is no conformer
    RWMol m;

    m.addAtom(new Atom(6));
    m.addAtom(new Atom(6));
    m.addAtom(new Atom(6));
    TEST_ASSERT(m.getNumConformers() == 0);
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testIssue267() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n";
  BOOST_LOG(rdInfoLog) << "Testing issue 267: default valence of *"
                       << std::endl;
  {
    RWMol m;

    m.addAtom(new Atom(0));
    m.updatePropertyCache();

    TEST_ASSERT(m.getAtomWithIdx(0)->getImplicitValence() == 0);
  }
  {
    RWMol m;

    m.addAtom(new Atom(0));
    for (unsigned int i = 0; i < 8; ++i) {
      m.addAtom(new Atom(1));
      m.addBond(0, i + 1, Bond::SINGLE);
    }
    m.updatePropertyCache();
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testIssue284() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n";
  BOOST_LOG(rdInfoLog) << "Testing issue 284: removeBond not updating indices"
                       << std::endl;
  {
    RWMol m;

    m.addAtom(new Atom(6));
    m.addAtom(new Atom(6));
    m.addAtom(new Atom(6));
    m.addBond(0, 1, Bond::SINGLE);
    m.addBond(1, 2, Bond::SINGLE);
    m.updatePropertyCache();
    TEST_ASSERT(m.getBondBetweenAtoms(0, 1)->getIdx() == 0);
    TEST_ASSERT(m.getBondBetweenAtoms(1, 2)->getIdx() == 1);
    m.removeBond(0, 1);
    TEST_ASSERT(!m.getBondBetweenAtoms(0, 1));
    TEST_ASSERT(m.getBondBetweenAtoms(1, 2)->getIdx() == 0);
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testAtomResidues() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n";
  BOOST_LOG(rdInfoLog) << "Testing residue information handling on atoms"
                       << std::endl;
  {
    RWMol *m = new RWMol();

    m->addAtom(new Atom(6));
    m->addAtom(new Atom(6));
    m->addBond(0, 1, Bond::SINGLE);
    m->addAtom(new Atom(6));
    m->addBond(1, 2, Bond::SINGLE);
    m->addAtom(new Atom(6));
    m->addBond(2, 3, Bond::SINGLE);

    TEST_ASSERT(!(m->getAtomWithIdx(0)->getMonomerInfo()));
    TEST_ASSERT(!(m->getAtomWithIdx(1)->getMonomerInfo()));
    TEST_ASSERT(!(m->getAtomWithIdx(2)->getMonomerInfo()));
    TEST_ASSERT(!(m->getAtomWithIdx(3)->getMonomerInfo()));

    m->getAtomWithIdx(0)
        ->setMonomerInfo(new AtomMonomerInfo(AtomMonomerInfo::OTHER, "m1"));
    TEST_ASSERT((m->getAtomWithIdx(0)->getMonomerInfo()));
    TEST_ASSERT(m->getAtomWithIdx(0)->getMonomerInfo()->getName() == "m1");

    m->getAtomWithIdx(1)->setMonomerInfo(new AtomPDBResidueInfo("Ca", 3));
    TEST_ASSERT((m->getAtomWithIdx(1)->getMonomerInfo()));
    TEST_ASSERT(m->getAtomWithIdx(1)->getMonomerInfo()->getName() == "Ca");
    TEST_ASSERT(
        static_cast<const AtomPDBResidueInfo *>(
            m->getAtomWithIdx(1)->getMonomerInfo())->getSerialNumber() == 3);

    RWMol *m2 = new RWMol(*m);
    delete m;

    TEST_ASSERT((m2->getAtomWithIdx(0)->getMonomerInfo()));
    TEST_ASSERT(m2->getAtomWithIdx(0)->getMonomerInfo()->getName() == "m1");
    TEST_ASSERT((m2->getAtomWithIdx(1)->getMonomerInfo()));
    TEST_ASSERT(m2->getAtomWithIdx(1)->getMonomerInfo()->getName() == "Ca");
    TEST_ASSERT(
        static_cast<const AtomPDBResidueInfo *>(
            m2->getAtomWithIdx(1)->getMonomerInfo())->getSerialNumber() == 3);
    TEST_ASSERT(!(m2->getAtomWithIdx(2)->getMonomerInfo()));
    TEST_ASSERT(!(m2->getAtomWithIdx(3)->getMonomerInfo()));
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testNeedsUpdatePropertyCache() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n";
  BOOST_LOG(rdInfoLog) << "Testing function needsUpdatePropertyCache"
                       << std::endl;
  {
    RWMol m;

    m.addAtom(new Atom(0));
    TEST_ASSERT(m.needsUpdatePropertyCache() == true);
    m.updatePropertyCache();

    TEST_ASSERT(m.getAtomWithIdx(0)->getImplicitValence() == 0);
    TEST_ASSERT(m.needsUpdatePropertyCache() == false);
  }
  {
    RWMol m;

    m.addAtom(new Atom(6));
    for (ROMol::AtomIterator atomIt = m.beginAtoms(); atomIt != m.endAtoms();
         ++atomIt) {
      (*atomIt)->calcExplicitValence(false);
      (*atomIt)->calcImplicitValence(false);
    }
    m.addAtom(new Atom(6));
    m.addBond(0, 1, Bond::SINGLE);
    TEST_ASSERT(m.needsUpdatePropertyCache() == true);
    m.updatePropertyCache();
    TEST_ASSERT(m.needsUpdatePropertyCache() == false);
  }
  {
    RWMol m;
    m.addAtom(new Atom(6));
    m.getAtomWithIdx(0)->calcExplicitValence(false);
    TEST_ASSERT(m.getAtomWithIdx(0)->needsUpdatePropertyCache());
    m.getAtomWithIdx(0)->setNoImplicit(true);
    TEST_ASSERT(!m.getAtomWithIdx(0)->needsUpdatePropertyCache());
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

namespace {
std::string qhelper(Atom::QUERYATOM_QUERY *q, unsigned int depth = 0) {
  std::string res = "";
  if (q) {
    for (unsigned int i = 0; i < depth; ++i) res += "  ";
    res += q->getFullDescription() + "\n";
    for (Atom::QUERYATOM_QUERY::CHILD_VECT_CI ci = q->beginChildren();
         ci != q->endChildren(); ++ci) {
      res += qhelper((*ci).get(), depth + 1);
    }
  }
  return res;
}
}

const char *m_als_mol =
    "\n"
    "  Marvin  08200814552D          \n"
    "\n"
    "  9  8  0  0  0  0            999 V2000\n"
    "   -1.9152    1.6205    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -1.0902    1.6205    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.5068    2.2039    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -2.3277    0.9061    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -2.3277    2.3350    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -3.1527    2.3350    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -3.6830    2.8727    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -3.1527    0.9061    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -3.6771    0.2814    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  2  2  0  0  0  0\n"
    "  2  3  1  0  0  0  0\n"
    "  1  4  1  0  0  0  0\n"
    "  1  5  1  0  0  0  0\n"
    "  5  6  2  0  0  0  0\n"
    "  6  7  1  0  0  0  0\n"
    "  4  8  2  0  0  0  0\n"
    "  8  9  1  0  0  0  0\n"
    "M  ALS   4  2 F O   Cl  \n"
    "M  END\n";

void testAtomListLineRoundTrip() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n";
  BOOST_LOG(rdInfoLog) << "Test AtomListLine RoundTrip" << std::endl;
  std::string rdbase = getenv("RDBASE");
  std::string fName = rdbase + "/Code/GraphMol/test_data/m_als_round_trip.mol";
  const bool sanitize = false;
  const bool removeHs = true;
  const bool strictParsing = true;
  unsigned int line = 0;

  std::istringstream inStream(m_als_mol);
  RWMol *m =
      MolDataStreamToMol(inStream, line, sanitize, removeHs, strictParsing);
  std::string desc = qhelper(m->getAtomWithIdx(3)->getQuery());
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 9);

  std::string molblock = MolToMolBlock(*m);
  std::istringstream inStream2(molblock);
  RWMol *m2 =
      MolDataStreamToMol(inStream2, line, sanitize, removeHs, strictParsing);
  TEST_ASSERT(m2);
  TEST_ASSERT(desc == qhelper(m2->getAtomWithIdx(3)->getQuery()));
  Atom::ATOM_SPTR cl(new Atom(17));
  Atom::ATOM_SPTR o(new Atom(17));
  TEST_ASSERT(dynamic_cast<QueryAtom *>(m->getAtomWithIdx(3))->Match(cl));
  TEST_ASSERT(dynamic_cast<QueryAtom *>(m->getAtomWithIdx(3))->Match(o));
  TEST_ASSERT(dynamic_cast<QueryAtom *>(m2->getAtomWithIdx(3))->Match(cl));
  TEST_ASSERT(dynamic_cast<QueryAtom *>(m2->getAtomWithIdx(3))->Match(o));
  delete m;
  delete m2;
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testGithub608() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n";
  BOOST_LOG(rdInfoLog) << "Test github 608: stereo bonds wrong after insertMol"
                       << std::endl;

  {
    RWMol *m = SmilesToMol("N1NN1");
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 3);
    RWMol *f = SmilesToMol("C/C=C/C");
    TEST_ASSERT(f);
    TEST_ASSERT(f->getNumAtoms() == 4);
    TEST_ASSERT(f->getBondBetweenAtoms(1, 2)->getStereoAtoms().size() == 2);
    TEST_ASSERT(f->getBondBetweenAtoms(1, 2)->getStereoAtoms()[0] == 0);
    TEST_ASSERT(f->getBondBetweenAtoms(1, 2)->getStereoAtoms()[1] == 3);

    m->insertMol(*f);
    TEST_ASSERT(m->getNumAtoms() == 7);
    TEST_ASSERT(m->getBondBetweenAtoms(4, 5)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(4, 5)->getStereoAtoms().size() == 2);
    TEST_ASSERT(m->getBondBetweenAtoms(4, 5)->getStereoAtoms()[0] == 3);
    TEST_ASSERT(m->getBondBetweenAtoms(4, 5)->getStereoAtoms()[1] == 6);

    delete m;
    delete f;
  }

  {
    INT_VECT nAtoms;
    RWMol *m = SmilesToMol("N1NN1");
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 3);
    RWMol *f = SmilesToMol("C[C@]1(F)CC[C@](Cl)(Br)CC1");
    TEST_ASSERT(f);
    TEST_ASSERT(f->getNumAtoms() == 10);
    TEST_ASSERT(f->getAtomWithIdx(1)->getPropIfPresent(
        common_properties::_ringStereoAtoms, nAtoms));
    TEST_ASSERT(std::find(nAtoms.begin(), nAtoms.end(), 6) != nAtoms.end());
    m->insertMol(*f);
    TEST_ASSERT(m->getNumAtoms() == 13);
    TEST_ASSERT(m->getAtomWithIdx(4)->getPropIfPresent(
        common_properties::_ringStereoAtoms, nAtoms));
    TEST_ASSERT(std::find(nAtoms.begin(), nAtoms.end(), 9) != nAtoms.end());

    delete m;
    delete f;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

// -------------------------------------------------------------------
int main() {
  RDLog::InitLogs();
// boost::logging::enable_logs("rdApp.info");
#if 1
  test1();
  testAddConformersFromTrajectory();
  testPropLeak();
  testMolProps();
  testAtomProps();
  testBondProps();
  testMisc();
  testDegree();
  testIssue1993296();
  testIssue2381580();
  testIssue2840217();
#endif
  testPeriodicTable();
  testAddAtomWithConf();
  testIssue267();
  testIssue284();
  testClearMol();
  testAtomResidues();
  testNeedsUpdatePropertyCache();
  testAtomListLineRoundTrip();
  testGithub608();

  return 0;
}
