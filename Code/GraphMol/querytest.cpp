//
//   Copyright (C) 2002-2017 Rational Discovery LLC and Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <iostream>

using namespace RDKit;
using namespace std;
typedef class RWMol Mol;

void test1() {
  BOOST_LOG(rdErrorLog) << "---------------------- Test1" << std::endl;
  Mol qM;
  Mol m;

  auto *a = new Atom(6);
  // we copy in addAtom, so this is safe
  m.addAtom(a);
  m.addAtom(a);
  delete a;
  m.addBond(0, 1, Bond::SINGLE);
  a = new Atom(8);
  m.addAtom(a);
  delete a;
  m.addBond(1, 2, Bond::DOUBLE);
  MolOps::sanitizeMol(m);

  auto *qA = new QueryAtom(6);
  CHECK_INVARIANT(qA->Match(m.getAtomWithIdx(0)), "");
  CHECK_INVARIANT(qA->Match(m.getAtomWithIdx(1)), "");
  CHECK_INVARIANT(!qA->Match(m.getAtomWithIdx(2)), "");
  qA->expandQuery(makeAtomImplicitValenceQuery(3));
  CHECK_INVARIANT(qA->Match(m.getAtomWithIdx(0)), "");
  CHECK_INVARIANT(!qA->Match(m.getAtomWithIdx(1)), "");
  CHECK_INVARIANT(!qA->Match(m.getAtomWithIdx(2)), "");

  qM.addAtom(qA);
  delete qA;
  qA = new QueryAtom(6);
  qA->expandQuery(makeAtomNumQuery(8), Queries::COMPOSITE_OR);
  qM.addAtom(qA);
  delete qA;
  qM.addAtom(new QueryAtom(8), true, true);

  QueryBond *qB;
  qB = new QueryBond(Bond::UNSPECIFIED);
  qB->setOwningMol(qM);
  qB->setBeginAtomIdx(0);
  qB->setEndAtomIdx(1);

  CHECK_INVARIANT(qB->Match(m.getBondWithIdx(0)), "");
  CHECK_INVARIANT(qB->Match(m.getBondWithIdx(1)), "");
  qM.addBond(qB, true);
  qB = new QueryBond(Bond::DOUBLE);
  qB->setOwningMol(qM);
  qB->setBeginAtomIdx(1);
  qB->setEndAtomIdx(2);
  qM.addBond(qB, true);

  CHECK_INVARIANT(qM.getAtomWithIdx(0)->Match(m.getAtomWithIdx(0)), "");
  CHECK_INVARIANT(!qM.getAtomWithIdx(0)->Match(m.getAtomWithIdx(1)), "");
  CHECK_INVARIANT(!qM.getAtomWithIdx(0)->Match(m.getAtomWithIdx(2)), "");
  CHECK_INVARIANT(qM.getAtomWithIdx(1)->Match(m.getAtomWithIdx(0)), "");
  CHECK_INVARIANT(qM.getAtomWithIdx(1)->Match(m.getAtomWithIdx(1)), "");
  CHECK_INVARIANT(qM.getAtomWithIdx(1)->Match(m.getAtomWithIdx(2)), "");
  CHECK_INVARIANT(!qM.getAtomWithIdx(2)->Match(m.getAtomWithIdx(0)), "");
  CHECK_INVARIANT(!qM.getAtomWithIdx(2)->Match(m.getAtomWithIdx(1)), "");
  CHECK_INVARIANT(qM.getAtomWithIdx(2)->Match(m.getAtomWithIdx(2)), "");

  CHECK_INVARIANT(qM.getBondWithIdx(0)->Match(m.getBondWithIdx(0)), "");
  CHECK_INVARIANT(qM.getBondWithIdx(0)->Match(m.getBondWithIdx(1)), "");
  CHECK_INVARIANT(!qM.getBondWithIdx(1)->Match(m.getBondWithIdx(0)), "");
  CHECK_INVARIANT(qM.getBondWithIdx(1)->Match(m.getBondWithIdx(1)), "");

  BOOST_LOG(rdErrorLog) << "Done!" << std::endl;
}

void test2() {
  BOOST_LOG(rdErrorLog) << "---------------------- Test2" << std::endl;
  Mol qM;
  Mol m;

  auto *a = new Atom(6);
  // we copy in addAtom, so this is safe
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  delete a;
  m.addBond(0, 1, Bond::TRIPLE);
  m.addBond(1, 2, Bond::SINGLE);
  a = new Atom(8);
  m.addAtom(a);
  delete a;
  m.addBond(2, 3, Bond::DOUBLE);
  MolOps::sanitizeMol(m);

  QueryBond *qB;
  qB = new QueryBond(Bond::SINGLE);

  CHECK_INVARIANT(!qB->Match(m.getBondWithIdx(0)), "");
  CHECK_INVARIANT(qB->Match(m.getBondWithIdx(1)), "");
  CHECK_INVARIANT(!qB->Match(m.getBondWithIdx(2)), "");

  BOND_EQUALS_QUERY *newQ = makeBondOrderEqualsQuery(Bond::DOUBLE);
  qB->expandQuery(newQ, Queries::COMPOSITE_OR);
  CHECK_INVARIANT(!qB->Match(m.getBondWithIdx(0)), "");
  CHECK_INVARIANT(qB->Match(m.getBondWithIdx(1)), "");
  CHECK_INVARIANT(qB->Match(m.getBondWithIdx(2)), "");

  delete qB;
  BOOST_LOG(rdErrorLog) << "Done!" << std::endl;
}

void test3() {
  BOOST_LOG(rdErrorLog) << "---------------------- Test3" << std::endl;
  Mol m;

  auto *a = new Atom(6);
  // we copy in addAtom, so this is safe
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  delete a;
  m.addBond(0, 1, Bond::SINGLE);
  m.addBond(1, 2, Bond::DOUBLE);
  m.addBond(2, 3, Bond::SINGLE);
  m.addBond(3, 4, Bond::DOUBLE);
  m.addBond(4, 5, Bond::SINGLE);
  m.addBond(5, 0, Bond::DOUBLE);
  m.addBond(5, 6, Bond::SINGLE);
  MolOps::sanitizeMol(m);

  ATOM_EQUALS_QUERY *aeq = makeAtomExplicitDegreeQuery(3);
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(0)), "");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(1)), "");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(5)), "");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(6)), "");
  delete aeq;
  aeq = makeAtomExplicitDegreeQuery(2);
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(0)), "");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(1)), "");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(5)), "");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(6)), "");
  delete aeq;

  aeq = makeAtomHCountQuery(1);
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(0)), "");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(1)), "");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(5)), "");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(6)), "");
  delete aeq;

  aeq = makeAtomInNRingsQuery(1);
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(0)), "");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(1)), "");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(5)), "");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(6)), "");
  delete aeq;

  aeq = makeAtomInNRingsQuery(0);
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(0)), "");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(1)), "");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(5)), "");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(6)), "");
  delete aeq;

  aeq = makeAtomAromaticQuery();
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(0)), "");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(1)), "");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(5)), "");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(6)), "");
  delete aeq;

  aeq = makeAtomAliphaticQuery();
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(0)), "");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(1)), "");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(5)), "");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(6)), "");
  delete aeq;

  aeq = makeAtomInRingQuery();
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(0)), "");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(1)), "");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(5)), "");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(6)), "");
  aeq->setNegation(true);
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(0)), "");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(1)), "");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(5)), "");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(6)), "");
  delete aeq;

  aeq = makeAtomInRingOfSizeQuery(6);
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(0)), "");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(1)), "");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(5)), "");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(6)), "");
  aeq->setNegation(true);
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(0)), "");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(1)), "");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(5)), "");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(6)), "");
  delete aeq;

  aeq = makeAtomInRingOfSizeQuery(5);
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(0)), "");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(1)), "");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(5)), "");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(6)), "");
  aeq->setNegation(true);
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(0)), "");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(1)), "");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(5)), "");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(6)), "");
  delete aeq;

  BOND_EQUALS_QUERY *beq = makeBondIsInRingQuery();
  CHECK_INVARIANT(beq->Match(m.getBondWithIdx(0)), "");
  CHECK_INVARIANT(!beq->Match(m.getBondWithIdx(6)), "");
  CHECK_INVARIANT(beq->Match(m.getBondBetweenAtoms(0, 1)), "");
  CHECK_INVARIANT(beq->Match(m.getBondBetweenAtoms(1, 0)), "");
  CHECK_INVARIANT(!beq->Match(m.getBondBetweenAtoms(5, 6)), "");
  CHECK_INVARIANT(!beq->Match(m.getBondBetweenAtoms(6, 5)), "");

  delete beq;
  BOOST_LOG(rdErrorLog) << "Done!" << std::endl;
}

void test4() {
  BOOST_LOG(rdErrorLog) << "---------------------- Test4" << std::endl;
  Mol m;

  auto *a = new Atom(6);
  // we copy in addAtom, so this is safe
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  delete a;
  m.addBond(0, 1, Bond::SINGLE);
  m.addBond(1, 2, Bond::DOUBLE);
  m.addBond(2, 3, Bond::SINGLE);
  m.addBond(3, 4, Bond::DOUBLE);
  m.addBond(4, 5, Bond::SINGLE);
  m.addBond(5, 0, Bond::DOUBLE);
  m.addBond(5, 6, Bond::SINGLE);
  MolOps::sanitizeMol(m);

  ATOM_NULL_QUERY *aeq = makeAtomNullQuery();
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(0)), "");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(1)), "");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(5)), "");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(6)), "");
  delete aeq;

  BOOST_LOG(rdErrorLog) << "Done!" << std::endl;
}

void test5() {
  BOOST_LOG(rdErrorLog) << "---------------------- Test5" << std::endl;
  Mol m;

  auto *a = new Atom(6);
  // we copy in addAtom, so this is safe
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  delete a;
  m.addBond(0, 1, Bond::SINGLE);
  m.addBond(1, 2, Bond::DOUBLE);
  m.addBond(2, 3, Bond::SINGLE);
  m.addBond(3, 4, Bond::DOUBLE);
  m.addBond(4, 5, Bond::SINGLE);
  m.addBond(5, 0, Bond::DOUBLE);
  m.addBond(5, 6, Bond::SINGLE);
  MolOps::sanitizeMol(m);

  unsigned int valenceProd;
  valenceProd = queryAtomBondProduct(m.getAtomWithIdx(0));
  TEST_ASSERT(valenceProd == 1681);  // aromatic*aromatic = 41 * 41 = 1681
  valenceProd = queryAtomBondProduct(m.getAtomWithIdx(1));
  TEST_ASSERT(valenceProd == 1681);

  valenceProd = queryAtomBondProduct(m.getAtomWithIdx(5));
  TEST_ASSERT(valenceProd == 5043);

  valenceProd = queryAtomBondProduct(m.getAtomWithIdx(6));
  TEST_ASSERT(valenceProd == 3);

  valenceProd = queryAtomAllBondProduct(m.getAtomWithIdx(6));
  TEST_ASSERT(valenceProd == 81);

  BOOST_LOG(rdErrorLog) << "Done!" << std::endl;
}

void testQueryQueryMatches() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing query--query matches" << std::endl;

  // =====================================
  // ATOMS
  // =====================================
  {
    QueryAtom a1, a2;
    a1.setQuery(makeAtomNullQuery());
    a2.setQuery(makeAtomNullQuery());
    TEST_ASSERT(a1.QueryMatch(&a2));
    TEST_ASSERT(a2.QueryMatch(&a1));
  }

  {
    QueryAtom a1, a2(6);
    a1.setQuery(makeAtomNullQuery());
    TEST_ASSERT(a1.QueryMatch(&a2));
    TEST_ASSERT(a2.QueryMatch(&a1));
  }

  {
    QueryAtom a1(6), a2(8);
    TEST_ASSERT(!a1.QueryMatch(&a2));
    TEST_ASSERT(!a2.QueryMatch(&a1));
  }

  {
    QueryAtom a1(6), a2(6);
    TEST_ASSERT(a1.QueryMatch(&a2));
    TEST_ASSERT(a2.QueryMatch(&a1));
  }

  {
    QueryAtom a1, a2;
    a1.setQuery(makeAtomAromaticQuery());
    a2.setQuery(makeAtomAromaticQuery());
    TEST_ASSERT(a1.QueryMatch(&a2));
    TEST_ASSERT(a2.QueryMatch(&a1));
  }

  {
    QueryAtom a1, a2;
    a1.setQuery(makeAtomAromaticQuery());
    a2.setQuery(makeAtomAliphaticQuery());
    TEST_ASSERT(!a1.QueryMatch(&a2));
    TEST_ASSERT(!a2.QueryMatch(&a1));
  }

  {
    QueryAtom a1(6), a2(8);
    a1.expandQuery(makeAtomNumQuery(8), Queries::COMPOSITE_OR);
    a2.expandQuery(makeAtomNumQuery(9), Queries::COMPOSITE_OR);
    TEST_ASSERT(a1.QueryMatch(&a2));
    TEST_ASSERT(a2.QueryMatch(&a1));
  }

  {
    QueryAtom a1(6), a2(8);
    a1.expandQuery(makeAtomNumQuery(7), Queries::COMPOSITE_OR);
    a2.expandQuery(makeAtomNumQuery(9), Queries::COMPOSITE_OR);
    TEST_ASSERT(!a1.QueryMatch(&a2));
    TEST_ASSERT(!a2.QueryMatch(&a1));
  }

  {
    QueryAtom a1(6), a2(6);
    a1.expandQuery(makeAtomExplicitValenceQuery(3), Queries::COMPOSITE_AND);
    a2.expandQuery(makeAtomExplicitValenceQuery(3), Queries::COMPOSITE_AND);
    TEST_ASSERT(a1.QueryMatch(&a2));
    TEST_ASSERT(a2.QueryMatch(&a1));
  }

  {
    QueryAtom a1(6), a2(6);
    a1.expandQuery(makeAtomExplicitValenceQuery(3), Queries::COMPOSITE_AND);
    a2.expandQuery(makeAtomExplicitValenceQuery(4), Queries::COMPOSITE_AND);
    TEST_ASSERT(!a1.QueryMatch(&a2));
    TEST_ASSERT(!a2.QueryMatch(&a1));
  }

  {
    QueryAtom a1(6), a2(8);
    a1.expandQuery(makeAtomExplicitValenceQuery(3), Queries::COMPOSITE_AND);
    a2.expandQuery(makeAtomExplicitValenceQuery(3), Queries::COMPOSITE_AND);
    TEST_ASSERT(!a1.QueryMatch(&a2));
    TEST_ASSERT(!a2.QueryMatch(&a1));
  }

  {
    QueryAtom a1(6), a2(8);
    a1.expandQuery(makeAtomNumQuery(8), Queries::COMPOSITE_OR);
    TEST_ASSERT(a1.QueryMatch(&a2));
    TEST_ASSERT(a2.QueryMatch(&a1));
  }

  // =====================================
  // BONDS
  // =====================================

  {
    QueryBond a1, a2;
    a1.setQuery(makeBondNullQuery());
    a2.setQuery(makeBondNullQuery());
    TEST_ASSERT(a1.Match(&a2));
    TEST_ASSERT(a2.QueryMatch(&a1));
  }

  {
    QueryBond a1, a2(Bond::SINGLE);
    a1.setQuery(makeBondNullQuery());
    TEST_ASSERT(a1.QueryMatch(&a2));
    TEST_ASSERT(a2.QueryMatch(&a1));
  }

  {
    QueryBond a1(Bond::SINGLE), a2(Bond::SINGLE);
    TEST_ASSERT(a1.QueryMatch(&a2));
    TEST_ASSERT(a2.QueryMatch(&a1));
  }

  {
    QueryBond a1(Bond::SINGLE), a2(Bond::DOUBLE);
    TEST_ASSERT(!a1.QueryMatch(&a2));
    TEST_ASSERT(!a2.QueryMatch(&a1));
  }

  {
    QueryBond a1(Bond::SINGLE), a2(Bond::DOUBLE);
    a1.expandQuery(makeBondOrderEqualsQuery(Bond::DOUBLE),
                   Queries::COMPOSITE_OR);
    a2.expandQuery(makeBondOrderEqualsQuery(Bond::AROMATIC),
                   Queries::COMPOSITE_OR);
    TEST_ASSERT(a1.QueryMatch(&a2));
    TEST_ASSERT(a2.QueryMatch(&a1));
  }

  {
    QueryBond a1(Bond::SINGLE), a2(Bond::DOUBLE);
    a1.expandQuery(makeBondOrderEqualsQuery(Bond::TRIPLE),
                   Queries::COMPOSITE_OR);
    a2.expandQuery(makeBondOrderEqualsQuery(Bond::AROMATIC),
                   Queries::COMPOSITE_OR);
    TEST_ASSERT(!a1.QueryMatch(&a2));
    TEST_ASSERT(!a2.QueryMatch(&a1));
  }

  {
    QueryBond a1(Bond::SINGLE), a2(Bond::DOUBLE);
    a1.expandQuery(makeBondMinRingSizeQuery(4), Queries::COMPOSITE_OR);
    a2.expandQuery(makeBondMinRingSizeQuery(4), Queries::COMPOSITE_OR);
    TEST_ASSERT(a1.QueryMatch(&a2));
    TEST_ASSERT(a2.QueryMatch(&a1));
  }

  {
    QueryBond a1(Bond::SINGLE), a2(Bond::DOUBLE);
    a1.expandQuery(makeBondMinRingSizeQuery(4), Queries::COMPOSITE_AND);
    a2.expandQuery(makeBondMinRingSizeQuery(4), Queries::COMPOSITE_AND);
    TEST_ASSERT(!a1.QueryMatch(&a2));
    TEST_ASSERT(!a2.QueryMatch(&a1));
  }

  {
    QueryBond a1(Bond::SINGLE), a2(Bond::SINGLE);
    a1.expandQuery(makeBondMinRingSizeQuery(5), Queries::COMPOSITE_AND);
    a2.expandQuery(makeBondMinRingSizeQuery(4), Queries::COMPOSITE_AND);
    TEST_ASSERT(!a1.QueryMatch(&a2));
    TEST_ASSERT(!a2.QueryMatch(&a1));
  }

  {
    QueryBond a1(Bond::SINGLE), a2(Bond::AROMATIC);
    a1.expandQuery(makeBondOrderEqualsQuery(Bond::AROMATIC),
                   Queries::COMPOSITE_OR);
    TEST_ASSERT(a1.QueryMatch(&a2));
    TEST_ASSERT(a2.QueryMatch(&a1));
  }

  BOOST_LOG(rdErrorLog) << "Done!" << std::endl;
}

void testIssue2892580() {
  BOOST_LOG(rdErrorLog) << "---------------------- Test issue 2892580"
                        << std::endl;
  Mol m;

  auto *a = new Atom(6);

  int massVal;
  massVal = queryAtomMass(a);
  TEST_ASSERT(massVal == static_cast<int>(RDKit::round(
                             12.011 * massIntegerConversionFactor)));

  a->setIsotope(13);
  massVal = queryAtomMass(a);
  TEST_ASSERT(massVal == static_cast<int>(RDKit::round(
                             13.003 * massIntegerConversionFactor)));

  delete a;
  BOOST_LOG(rdErrorLog) << "Done!" << std::endl;
}

void testGithub153() {
  BOOST_LOG(rdErrorLog) << "---------------------- Test github issue 53: query "
                           "molecules not matching [R]"
                        << std::endl;
  RWMol *m = SmartsToMol("[C]1-[C]-[C]1");
  MolOps::findSSSR(*m);

  RWMol *q = SmartsToMol("[R]");

  std::vector<MatchVectType> mvv;
  TEST_ASSERT(SubstructMatch(*m, *q, mvv));
  TEST_ASSERT(mvv.size() == 3);
  TEST_ASSERT(mvv[0].size() == 1);

  delete m;
  delete q;

  BOOST_LOG(rdErrorLog) << "Done!" << std::endl;
}

void testQualifiedQueries() {
  BOOST_LOG(rdErrorLog)
      << "---------------------- Test queries using qualifiers instead of =="
      << std::endl;
  RWMol *m = SmilesToMol("CNO");

  {
    QueryAtom qA;
    qA.setQuery(makeAtomNumQuery<ATOM_GREATER_QUERY>(7, "test"));
    TEST_ASSERT(qA.Match(m->getAtomWithIdx(0)));
    TEST_ASSERT(!qA.Match(m->getAtomWithIdx(1)));
    TEST_ASSERT(!qA.Match(m->getAtomWithIdx(2)));
  }
  {
    QueryAtom qA;
    qA.setQuery(makeAtomNumQuery<ATOM_GREATEREQUAL_QUERY>(7, "test"));
    TEST_ASSERT(qA.Match(m->getAtomWithIdx(0)));
    TEST_ASSERT(qA.Match(m->getAtomWithIdx(1)));
    TEST_ASSERT(!qA.Match(m->getAtomWithIdx(2)));
  }
  {
    QueryAtom qA;
    qA.setQuery(makeAtomNumQuery<ATOM_LESS_QUERY>(7, "test"));
    TEST_ASSERT(!qA.Match(m->getAtomWithIdx(0)));
    TEST_ASSERT(!qA.Match(m->getAtomWithIdx(1)));
    TEST_ASSERT(qA.Match(m->getAtomWithIdx(2)));
  }
  {
    QueryAtom qA;
    qA.setQuery(makeAtomNumQuery<ATOM_LESSEQUAL_QUERY>(7, "test"));
    TEST_ASSERT(!qA.Match(m->getAtomWithIdx(0)));
    TEST_ASSERT(qA.Match(m->getAtomWithIdx(1)));
    TEST_ASSERT(qA.Match(m->getAtomWithIdx(2)));
  }

  delete m;
  BOOST_LOG(rdErrorLog) << "Done!" << std::endl;
}

void testGithub165() {
  BOOST_LOG(rdErrorLog) << "---------------------- Test Github issue 165: "
                           "radicals not used in atom-atom matching"
                        << std::endl;

  Atom a1(6);
  Atom a2(6);

  TEST_ASSERT(a1.Match(&a2));
  TEST_ASSERT(a2.Match(&a1));
  a1.setNumRadicalElectrons(2);
  TEST_ASSERT(!a1.Match(&a2));
  TEST_ASSERT(a2.Match(&a1));
  a2.setNumRadicalElectrons(2);
  TEST_ASSERT(a1.Match(&a2));
  TEST_ASSERT(a2.Match(&a1));
  a2.setNumRadicalElectrons(3);
  TEST_ASSERT(!a1.Match(&a2));
  TEST_ASSERT(!a2.Match(&a1));

  BOOST_LOG(rdErrorLog) << "Done!" << std::endl;
}

void testHasPropMatch() {
  QueryAtom qA;
  qA.setQuery(makeHasPropQuery<Atom>("foo"));
  Atom a1(6);
  TEST_ASSERT(!qA.Match(&a1));
  a1.setProp<int>("foo", 1);
  TEST_ASSERT(qA.Match(&a1));
}

void testHasPropWithValueMatch() {
  {
    QueryAtom qA;
    qA.setQuery(makePropQuery<Atom, int>("foo", 2));
    Atom a1(6);
    TEST_ASSERT(!qA.Match(&a1));
    a1.setProp<int>("foo", 1);
    TEST_ASSERT(!qA.Match(&a1));
    a1.setProp<int>("foo", 2);
    TEST_ASSERT(qA.Match(&a1));

    a1.clearProp("foo");
    a1.setProp<double>("foo", 2);
    TEST_ASSERT(!qA.Match(&a1));
  }

  {
    QueryAtom qA;
    qA.setQuery(makePropQuery<Atom, std::string>("foo", "bar"));
    Atom a1(6);
    TEST_ASSERT(!qA.Match(&a1));
    a1.setProp<std::string>("foo", "bar");
    TEST_ASSERT(qA.Match(&a1));
  }

  {
    QueryBond qA;
    qA.setQuery(makePropQuery<Bond, int>("foo", 2));
    Bond a1;
    TEST_ASSERT(!qA.Match(&a1));
    a1.setProp<int>("foo", 1);
    TEST_ASSERT(!qA.Match(&a1));
    a1.setProp<int>("foo", 2);
    TEST_ASSERT(qA.Match(&a1));

    a1.clearProp("foo");
    a1.setProp<double>("foo", 2);
    TEST_ASSERT(!qA.Match(&a1));
  }

  {
    QueryBond qA;
    qA.setQuery(makePropQuery<Bond, std::string>("foo", "bar"));
    Bond a1;
    TEST_ASSERT(!qA.Match(&a1));
    a1.setProp<std::string>("foo", "bar");
    TEST_ASSERT(qA.Match(&a1));
  }
}

void testHasPropWithDoubleValueMatch() {
  {
    QueryAtom qA;
    qA.setQuery(makePropQuery<Atom, double>("foo", 2));
    Atom a1(6);
    TEST_ASSERT(!qA.Match(&a1));
    a1.setProp<double>("foo", 1);
    TEST_ASSERT(!qA.Match(&a1));
    a1.setProp<double>("foo", 2);
    TEST_ASSERT(qA.Match(&a1));
  }
  {
    QueryBond qA;
    qA.setQuery(makePropQuery<Bond, double>("foo", 2));
    Bond a1;
    TEST_ASSERT(!qA.Match(&a1));
    a1.setProp<double>("foo", 1);
    TEST_ASSERT(!qA.Match(&a1));
    a1.setProp<double>("foo", 2);
    TEST_ASSERT(qA.Match(&a1));
  }
}

void testExtraAtomQueries() {
  BOOST_LOG(rdErrorLog) << "---------------------- Test extra atom queries"
                        << std::endl;

  {  // radicals
    QueryAtom qA;
    qA.setQuery(makeAtomNumRadicalElectronsQuery(1));
    Atom a1(6);
    TEST_ASSERT(!qA.Match(&a1));
    a1.setNumRadicalElectrons(1);
    TEST_ASSERT(qA.Match(&a1));
    a1.setNumRadicalElectrons(2);
    TEST_ASSERT(!qA.Match(&a1));
    qA.getQuery()->setNegation(true);
    TEST_ASSERT(qA.Match(&a1));
    a1.setNumRadicalElectrons(0);
    TEST_ASSERT(qA.Match(&a1));
    a1.setNumRadicalElectrons(1);
    TEST_ASSERT(!qA.Match(&a1));
  }
  {  // chiral tags
    QueryAtom qA;
    qA.setQuery(makeAtomHasChiralTagQuery());
    Atom a1(6);
    TEST_ASSERT(!qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
    TEST_ASSERT(qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
    TEST_ASSERT(qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_OTHER);
    TEST_ASSERT(qA.Match(&a1));
    qA.getQuery()->setNegation(true);
    a1.setChiralTag(Atom::CHI_UNSPECIFIED);
    TEST_ASSERT(qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
    TEST_ASSERT(!qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
    TEST_ASSERT(!qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_OTHER);
    TEST_ASSERT(!qA.Match(&a1));
  }

  {  // missing chiral tags
    QueryAtom qA;
    qA.setQuery(makeAtomMissingChiralTagQuery());
    Atom a1(6);
    TEST_ASSERT(!qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
    TEST_ASSERT(!qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
    TEST_ASSERT(!qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_OTHER);
    TEST_ASSERT(!qA.Match(&a1));

    a1.setChiralTag(Atom::CHI_UNSPECIFIED);
    a1.setProp(common_properties::_ChiralityPossible, 1);
    TEST_ASSERT(qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
    TEST_ASSERT(!qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
    TEST_ASSERT(!qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_OTHER);
    TEST_ASSERT(!qA.Match(&a1));

    qA.getQuery()->setNegation(true);
    a1.clearProp(common_properties::_ChiralityPossible);
    a1.setChiralTag(Atom::CHI_UNSPECIFIED);
    TEST_ASSERT(qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
    TEST_ASSERT(qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
    TEST_ASSERT(qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_OTHER);
    TEST_ASSERT(qA.Match(&a1));

    a1.setChiralTag(Atom::CHI_UNSPECIFIED);
    a1.setProp(common_properties::_ChiralityPossible, 1);
    TEST_ASSERT(!qA.Match(&a1));
  }

  BOOST_LOG(rdErrorLog) << "Done!" << std::endl;
}
void testExtraBondQueries() {
  BOOST_LOG(rdErrorLog) << "---------------------- Test extra bond queries"
                        << std::endl;

  {  // stereo
    QueryBond qB;
    qB.setQuery(makeBondHasStereoQuery());
    Bond b1(Bond::DOUBLE);
    TEST_ASSERT(!qB.Match(&b1));
    b1.setStereo(Bond::STEREOE);
    TEST_ASSERT(qB.Match(&b1));
    b1.setStereo(Bond::STEREOANY);
    TEST_ASSERT(qB.Match(&b1));
    qB.getQuery()->setNegation(true);
    b1.setStereo(Bond::STEREONONE);
    TEST_ASSERT(qB.Match(&b1));
    b1.setStereo(Bond::STEREOE);
    TEST_ASSERT(!qB.Match(&b1));
    b1.setStereo(Bond::STEREOANY);
    TEST_ASSERT(!qB.Match(&b1));
  }
  BOOST_LOG(rdErrorLog) << "Done!" << std::endl;
}
void testNumHeteroatomNeighborQueries() {
  BOOST_LOG(rdErrorLog)
      << "---------------------- Test num heteroatom neighbor queries"
      << std::endl;

  RWMol *m = SmilesToMol("CCNCOO");
  {
    QueryAtom qA;
    qA.setQuery(makeAtomNumHeteroatomNbrsQuery(1));
    TEST_ASSERT(!qA.Match(m->getAtomWithIdx(0)));
    TEST_ASSERT(qA.Match(m->getAtomWithIdx(1)));
    TEST_ASSERT(!qA.Match(m->getAtomWithIdx(2)));
    TEST_ASSERT(!qA.Match(m->getAtomWithIdx(3)));
    TEST_ASSERT(qA.Match(m->getAtomWithIdx(5)));
  }
  {
    QueryAtom qA;
    qA.setQuery(
        makeAtomRangeQuery(0, 3, true, true, queryAtomNumHeteroatomNbrs));
    TEST_ASSERT(!qA.Match(m->getAtomWithIdx(0)));
    TEST_ASSERT(qA.Match(m->getAtomWithIdx(1)));
    TEST_ASSERT(!qA.Match(m->getAtomWithIdx(2)));
    TEST_ASSERT(qA.Match(m->getAtomWithIdx(3)));
    TEST_ASSERT(qA.Match(m->getAtomWithIdx(5)));
  }
  {
    QueryAtom qA;
    qA.setQuery(
        makeAtomRangeQuery(1, 2, false, false, queryAtomNumHeteroatomNbrs));
    TEST_ASSERT(!qA.Match(m->getAtomWithIdx(0)));
    TEST_ASSERT(qA.Match(m->getAtomWithIdx(1)));
    TEST_ASSERT(!qA.Match(m->getAtomWithIdx(2)));
    TEST_ASSERT(qA.Match(m->getAtomWithIdx(3)));
    TEST_ASSERT(qA.Match(m->getAtomWithIdx(5)));
  }
  {
    QueryAtom qA;
    qA.setQuery(
        makeAtomRangeQuery(0, 2, true, false, queryAtomNumHeteroatomNbrs));
    TEST_ASSERT(!qA.Match(m->getAtomWithIdx(0)));
    TEST_ASSERT(qA.Match(m->getAtomWithIdx(1)));
    TEST_ASSERT(!qA.Match(m->getAtomWithIdx(2)));
    TEST_ASSERT(qA.Match(m->getAtomWithIdx(3)));
    TEST_ASSERT(qA.Match(m->getAtomWithIdx(5)));
  }
  delete m;

  BOOST_LOG(rdErrorLog) << "Done!" << std::endl;
}


void testAtomTypeQueries() {
  BOOST_LOG(rdErrorLog)
      << "---------------------- Test atom type queries"
      << std::endl;

  RWMol *m = SmilesToMol("CCc1ccccc1");
  {
    QueryAtom qA1;
    qA1.setQuery(makeAtomTypeQuery(6,true));
    QueryAtom qA2;
    qA2.setQuery(makeAtomTypeQuery(6,false));
    QueryAtom qA3;
    qA3.setQuery(makeAtomTypeQuery(7,true));
    TEST_ASSERT(!qA1.Match(m->getAtomWithIdx(0)));
    TEST_ASSERT(qA2.Match(m->getAtomWithIdx(0)));
    TEST_ASSERT(!qA3.Match(m->getAtomWithIdx(0)));
    TEST_ASSERT(qA1.Match(m->getAtomWithIdx(2)));
    TEST_ASSERT(!qA2.Match(m->getAtomWithIdx(2)));
    TEST_ASSERT(!qA3.Match(m->getAtomWithIdx(2)));
  }

  delete m;
  BOOST_LOG(rdErrorLog) << "Done!" << std::endl;
}

int main() {
  RDLog::InitLogs();
#if 1
  test1();
  test2();
  test3();
  test4();
  test5();
  testQueryQueryMatches();
  testIssue2892580();
  testGithub153();
  testQualifiedQueries();
  testGithub165();
  testHasPropMatch();
  testHasPropWithValueMatch();
  testHasPropWithDoubleValueMatch();
#endif
  testExtraAtomQueries();
  testExtraBondQueries();
  testNumHeteroatomNeighborQueries();
  testAtomTypeQueries();
  return 0;
}
