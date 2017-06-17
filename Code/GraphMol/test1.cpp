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
  m2.getProp(RDKit::detail::computedPropName, cplst);
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
  m2.getProp(RDKit::detail::computedPropName, cplst);
  CHECK_INVARIANT(cplst.size() == 1, "");

  m2.clearComputedProps();
  CHECK_INVARIANT(!m2.hasProp("cprop2"), "");
  m2.getProp(RDKit::detail::computedPropName, cplst);
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

  TEST_ASSERT(m2.hasProp(RDKit::detail::computedPropName));

  m2.clear();
  TEST_ASSERT(!m2.hasProp("prop1"));
  TEST_ASSERT(m2.getNumAtoms() == 0);
  TEST_ASSERT(m2.getNumBonds() == 0);
  TEST_ASSERT(m2.getAtomBookmarks()->empty());
  TEST_ASSERT(m2.getBondBookmarks()->empty());

  TEST_ASSERT(
      m2.hasProp(RDKit::detail::computedPropName));  // <- github issue 176
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
  a1->getProp(RDKit::detail::computedPropName, cplst);
  CHECK_INVARIANT(cplst.size() == 2, "");
  CHECK_INVARIANT(cplst[0] == "cprop1", "");
  CHECK_INVARIANT(cplst[1] == "cprop2", "");

  a1->clearProp("cprop1");
  CHECK_INVARIANT(!a1->hasProp("cprop1"), "");
  a1->getProp(RDKit::detail::computedPropName, cplst);
  CHECK_INVARIANT(cplst.size() == 1, "");

  a1->clearComputedProps();
  CHECK_INVARIANT(!a1->hasProp("cprop2"), "");
  a1->getProp(RDKit::detail::computedPropName, cplst);
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
  b1->getProp(RDKit::detail::computedPropName, cplst);
  CHECK_INVARIANT(cplst.size() == 2, "");
  CHECK_INVARIANT(cplst[0] == "cprop1", "");
  CHECK_INVARIANT(cplst[1] == "cprop2", "");

  b1->clearProp("cprop1");
  CHECK_INVARIANT(!b1->hasProp("cprop1"), "");
  b1->getProp(RDKit::detail::computedPropName, cplst);
  CHECK_INVARIANT(cplst.size() == 1, "");

  b1->clearComputedProps();
  CHECK_INVARIANT(!b1->hasProp("cprop2"), "");
  b1->getProp(RDKit::detail::computedPropName, cplst);
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
    delete repA;
    TEST_ASSERT(m.getAtomWithIdx(newIdx)->getAtomicNum() == 22);
    Bond *nbnd = new Bond(Bond::DOUBLE);
    TEST_ASSERT(m.getBondWithIdx(m.getNumBonds() - 1)->getBondType() ==
                Bond::AROMATIC);
    m.replaceBond(m.getNumBonds() - 1, nbnd);
    m.debugMol(std::cerr);
    TEST_ASSERT(m.getBondWithIdx(m.getNumBonds() - 1)->getBondType() ==
                nbnd->getBondType());
    delete nbnd;
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
                          "already have conformers"
                       << std::endl;
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

    m->getAtomWithIdx(0)->setMonomerInfo(
        new AtomMonomerInfo(AtomMonomerInfo::OTHER, "m1"));
    TEST_ASSERT((m->getAtomWithIdx(0)->getMonomerInfo()));
    TEST_ASSERT(m->getAtomWithIdx(0)->getMonomerInfo()->getName() == "m1");

    m->getAtomWithIdx(1)->setMonomerInfo(new AtomPDBResidueInfo("Ca", 3));
    TEST_ASSERT((m->getAtomWithIdx(1)->getMonomerInfo()));
    TEST_ASSERT(m->getAtomWithIdx(1)->getMonomerInfo()->getName() == "Ca");
    TEST_ASSERT(static_cast<const AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(1)->getMonomerInfo())
                    ->getSerialNumber() == 3);

    RWMol *m2 = new RWMol(*m);
    delete m;

    TEST_ASSERT((m2->getAtomWithIdx(0)->getMonomerInfo()));
    TEST_ASSERT(m2->getAtomWithIdx(0)->getMonomerInfo()->getName() == "m1");
    TEST_ASSERT((m2->getAtomWithIdx(1)->getMonomerInfo()));
    TEST_ASSERT(m2->getAtomWithIdx(1)->getMonomerInfo()->getName() == "Ca");
    TEST_ASSERT(static_cast<const AtomPDBResidueInfo *>(
                    m2->getAtomWithIdx(1)->getMonomerInfo())
                    ->getSerialNumber() == 3);
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

#ifdef RDK_TEST_MULTITHREADED
#include <RDGeneral/BoostStartInclude.h>
#include <boost/thread.hpp>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/BoostEndInclude.h>
namespace {
void runblock(std::vector<const PeriodicTable *> *pts, int idx) {
  const PeriodicTable *pt = PeriodicTable::getTable();
  (*pts)[idx] = pt;
  TEST_ASSERT(pt->getAtomicNumber("C") == 6);
  TEST_ASSERT(pt->getAtomicNumber("N") == 7);
  TEST_ASSERT(pt->getAtomicNumber("O") == 8);
};
}
void testGithub381() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "    Test github381: thread-safe initialization of the periodic table"
      << std::endl;

  boost::thread_group tg;
  unsigned int count = 32;
  std::vector<const PeriodicTable *> pts(count);
#if 1
  for (unsigned int i = 0; i < count; ++i) {
    std::cerr.flush();
    tg.add_thread(new boost::thread(runblock, &pts, i));
  }
  tg.join_all();
  TEST_ASSERT(pts[0] != NULL);
  for (unsigned int i = 1; i < count; ++i) {
    TEST_ASSERT(pts[i] == pts[0]);
  }
#endif

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}
#else
void testGithub381() {}
#endif

void testGithub1041() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test github1041: Segfault for atom with no "
                           "owner (expect some warnings)"
                        << std::endl;
  {
    Atom at(6);
    bool ok = false;
    try {
      at.getOwningMol();
    } catch (const Invar::Invariant &err) {
      ok = true;
    }
    TEST_ASSERT(ok);
  }
  {
    Bond b;
    bool ok = false;
    try {
      b.getOwningMol();
    } catch (const Invar::Invariant &err) {
      ok = true;
    }
    TEST_ASSERT(ok);
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testGithub1453() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "    Test github1453: RWMol.clear() should reset the numBonds count."
      << std::endl;

  RWMol m2;

  m2.addAtom(new Atom(6));
  m2.addAtom(new Atom(6));
  m2.addBond(0, 1, Bond::TRIPLE);
  TEST_ASSERT(m2.getNumAtoms() == 2);
  TEST_ASSERT(m2.getNumBonds() == 1);
  m2.clear();
  TEST_ASSERT(m2.getNumAtoms() == 0);
  TEST_ASSERT(m2.getNumBonds() == 0);

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

// -------------------------------------------------------------------
int main() {
  RDLog::InitLogs();
// boost::logging::enable_logs("rdApp.info");
#if 1
  test1();
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
  testGithub381();
  testGithub1041();
  testGithub1041();
  testGithub1453();

  return 0;
}
