//
//   Copyright (C) 2002-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>

#include <memory>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

using namespace RDKit;
using namespace std;
typedef class RWMol Mol;

TEST_CASE("query atom and bond matching") {
  Mol qM;
  Mol m;

  Atom a(6);
  m.addAtom(&a);
  m.addAtom(&a);
  m.addBond(0, 1, Bond::SINGLE);
  a.setAtomicNum(8);
  m.addAtom(&a);
  m.addBond(1, 2, Bond::DOUBLE);
  MolOps::sanitizeMol(m);

  QueryAtom qA(6);
  CHECK(qA.Match(m.getAtomWithIdx(0)));
  CHECK(qA.Match(m.getAtomWithIdx(1)));
  CHECK(!qA.Match(m.getAtomWithIdx(2)));
  qA.expandQuery(makeAtomImplicitValenceQuery(3));
  CHECK(qA.Match(m.getAtomWithIdx(0)));
  CHECK(!qA.Match(m.getAtomWithIdx(1)));
  CHECK(!qA.Match(m.getAtomWithIdx(2)));

  qM.addAtom(&qA);
  QueryAtom qA2(6);
  qA2.expandQuery(makeAtomNumQuery(8), Queries::COMPOSITE_OR);
  qM.addAtom(&qA2);
  QueryAtom qA3(8);
  qM.addAtom(&qA3);

  QueryBond qB(Bond::UNSPECIFIED);
  qB.setOwningMol(qM);
  qB.setBeginAtomIdx(0);
  qB.setEndAtomIdx(1);

  CHECK(qB.Match(m.getBondWithIdx(0)));
  CHECK(qB.Match(m.getBondWithIdx(1)));
  qM.addBond(&qB);
  QueryBond qB2(Bond::DOUBLE);
  qB2.setOwningMol(qM);
  qB2.setBeginAtomIdx(1);
  qB2.setEndAtomIdx(2);
  qM.addBond(&qB2);

  CHECK(qM.getAtomWithIdx(0)->Match(m.getAtomWithIdx(0)));
  CHECK(!qM.getAtomWithIdx(0)->Match(m.getAtomWithIdx(1)));
  CHECK(!qM.getAtomWithIdx(0)->Match(m.getAtomWithIdx(2)));
  CHECK(qM.getAtomWithIdx(1)->Match(m.getAtomWithIdx(0)));
  CHECK(qM.getAtomWithIdx(1)->Match(m.getAtomWithIdx(1)));
  CHECK(qM.getAtomWithIdx(1)->Match(m.getAtomWithIdx(2)));
  CHECK(!qM.getAtomWithIdx(2)->Match(m.getAtomWithIdx(0)));
  CHECK(!qM.getAtomWithIdx(2)->Match(m.getAtomWithIdx(1)));
  CHECK(qM.getAtomWithIdx(2)->Match(m.getAtomWithIdx(2)));

  CHECK(qM.getBondWithIdx(0)->Match(m.getBondWithIdx(0)));
  CHECK(qM.getBondWithIdx(0)->Match(m.getBondWithIdx(1)));
  CHECK(!qM.getBondWithIdx(1)->Match(m.getBondWithIdx(0)));
  CHECK(qM.getBondWithIdx(1)->Match(m.getBondWithIdx(1)));
}

TEST_CASE("query bond order matching") {
  Mol qM;
  Mol m;

  Atom a(6);
  m.addAtom(&a);
  m.addAtom(&a);
  m.addAtom(&a);
  m.addBond(0, 1, Bond::TRIPLE);
  m.addBond(1, 2, Bond::SINGLE);
  a.setAtomicNum(8);
  m.addAtom(&a);
  m.addBond(2, 3, Bond::DOUBLE);
  MolOps::sanitizeMol(m);

  auto qB = std::make_unique<QueryBond>(Bond::SINGLE);

  CHECK(!qB->Match(m.getBondWithIdx(0)));
  CHECK(qB->Match(m.getBondWithIdx(1)));
  CHECK(!qB->Match(m.getBondWithIdx(2)));

  auto newQ = std::unique_ptr<BOND_EQUALS_QUERY>(
      makeBondOrderEqualsQuery(Bond::DOUBLE));
  qB->expandQuery(newQ.release(), Queries::COMPOSITE_OR);
  CHECK(!qB->Match(m.getBondWithIdx(0)));
  CHECK(qB->Match(m.getBondWithIdx(1)));
  CHECK(qB->Match(m.getBondWithIdx(2)));
}

TEST_CASE("atom and bond query operations") {
  Mol m;

  Atom a(6);
  m.addAtom(&a);
  m.addAtom(&a);
  m.addAtom(&a);
  m.addAtom(&a);
  m.addAtom(&a);
  m.addAtom(&a);
  m.addAtom(&a);
  m.addBond(0, 1, Bond::SINGLE);
  m.addBond(1, 2, Bond::DOUBLE);
  m.addBond(2, 3, Bond::SINGLE);
  m.addBond(3, 4, Bond::DOUBLE);
  m.addBond(4, 5, Bond::SINGLE);
  m.addBond(5, 0, Bond::DOUBLE);
  m.addBond(5, 6, Bond::SINGLE);
  MolOps::sanitizeMol(m);

  std::unique_ptr<ATOM_EQUALS_QUERY> aeq(makeAtomExplicitDegreeQuery(3));
  CHECK(!aeq->Match(m.getAtomWithIdx(0)));
  CHECK(!aeq->Match(m.getAtomWithIdx(1)));
  CHECK(aeq->Match(m.getAtomWithIdx(5)));
  CHECK(!aeq->Match(m.getAtomWithIdx(6)));
  aeq.reset(makeAtomExplicitDegreeQuery(2));
  CHECK(aeq->Match(m.getAtomWithIdx(0)));
  CHECK(aeq->Match(m.getAtomWithIdx(1)));
  CHECK(!aeq->Match(m.getAtomWithIdx(5)));
  CHECK(!aeq->Match(m.getAtomWithIdx(6)));

  aeq.reset(makeAtomHCountQuery(1));
  CHECK(aeq->Match(m.getAtomWithIdx(0)));
  CHECK(aeq->Match(m.getAtomWithIdx(1)));
  CHECK(!aeq->Match(m.getAtomWithIdx(5)));
  CHECK(!aeq->Match(m.getAtomWithIdx(6)));

  aeq.reset(makeAtomInNRingsQuery(1));
  CHECK(aeq->Match(m.getAtomWithIdx(0)));
  CHECK(aeq->Match(m.getAtomWithIdx(1)));
  CHECK(aeq->Match(m.getAtomWithIdx(5)));
  CHECK(!aeq->Match(m.getAtomWithIdx(6)));

  aeq.reset(makeAtomInNRingsQuery(0));
  CHECK(!aeq->Match(m.getAtomWithIdx(0)));
  CHECK(!aeq->Match(m.getAtomWithIdx(1)));
  CHECK(!aeq->Match(m.getAtomWithIdx(5)));
  CHECK(aeq->Match(m.getAtomWithIdx(6)));

  aeq.reset(makeAtomAromaticQuery());
  CHECK(aeq->Match(m.getAtomWithIdx(0)));
  CHECK(aeq->Match(m.getAtomWithIdx(1)));
  CHECK(aeq->Match(m.getAtomWithIdx(5)));
  CHECK(!aeq->Match(m.getAtomWithIdx(6)));

  aeq.reset(makeAtomAliphaticQuery());
  CHECK(!aeq->Match(m.getAtomWithIdx(0)));
  CHECK(!aeq->Match(m.getAtomWithIdx(1)));
  CHECK(!aeq->Match(m.getAtomWithIdx(5)));
  CHECK(aeq->Match(m.getAtomWithIdx(6)));

  aeq.reset(makeAtomInRingQuery());
  CHECK(aeq->Match(m.getAtomWithIdx(0)));
  CHECK(aeq->Match(m.getAtomWithIdx(1)));
  CHECK(aeq->Match(m.getAtomWithIdx(5)));
  CHECK(!aeq->Match(m.getAtomWithIdx(6)));
  aeq->setNegation(true);
  CHECK(!aeq->Match(m.getAtomWithIdx(0)));
  CHECK(!aeq->Match(m.getAtomWithIdx(1)));
  CHECK(!aeq->Match(m.getAtomWithIdx(5)));
  CHECK(aeq->Match(m.getAtomWithIdx(6)));

  aeq.reset(makeAtomInRingOfSizeQuery(6));
  CHECK(aeq->Match(m.getAtomWithIdx(0)));
  CHECK(aeq->Match(m.getAtomWithIdx(1)));
  CHECK(aeq->Match(m.getAtomWithIdx(5)));
  CHECK(!aeq->Match(m.getAtomWithIdx(6)));
  aeq->setNegation(true);
  CHECK(!aeq->Match(m.getAtomWithIdx(0)));
  CHECK(!aeq->Match(m.getAtomWithIdx(1)));
  CHECK(!aeq->Match(m.getAtomWithIdx(5)));
  CHECK(aeq->Match(m.getAtomWithIdx(6)));

  aeq.reset(makeAtomInRingOfSizeQuery(5));
  CHECK(!aeq->Match(m.getAtomWithIdx(0)));
  CHECK(!aeq->Match(m.getAtomWithIdx(1)));
  CHECK(!aeq->Match(m.getAtomWithIdx(5)));
  CHECK(!aeq->Match(m.getAtomWithIdx(6)));
  aeq->setNegation(true);
  CHECK(aeq->Match(m.getAtomWithIdx(0)));
  CHECK(aeq->Match(m.getAtomWithIdx(1)));
  CHECK(aeq->Match(m.getAtomWithIdx(5)));
  CHECK(aeq->Match(m.getAtomWithIdx(6)));

  std::unique_ptr<BOND_EQUALS_QUERY> beq(makeBondIsInRingQuery());
  CHECK(beq->Match(m.getBondWithIdx(0)));
  CHECK(!beq->Match(m.getBondWithIdx(6)));
  CHECK(beq->Match(m.getBondBetweenAtoms(0, 1)));
  CHECK(beq->Match(m.getBondBetweenAtoms(1, 0)));
  CHECK(!beq->Match(m.getBondBetweenAtoms(5, 6)));
  CHECK(!beq->Match(m.getBondBetweenAtoms(6, 5)));
}

TEST_CASE("null atom queries") {
  Mol m;

  Atom a(6);
  m.addAtom(&a);
  m.addAtom(&a);
  m.addAtom(&a);
  m.addAtom(&a);
  m.addAtom(&a);
  m.addAtom(&a);
  m.addAtom(&a);
  m.addBond(0, 1, Bond::SINGLE);
  m.addBond(1, 2, Bond::DOUBLE);
  m.addBond(2, 3, Bond::SINGLE);
  m.addBond(3, 4, Bond::DOUBLE);
  m.addBond(4, 5, Bond::SINGLE);
  m.addBond(5, 0, Bond::DOUBLE);
  m.addBond(5, 6, Bond::SINGLE);
  MolOps::sanitizeMol(m);

  std::unique_ptr<ATOM_NULL_QUERY> aeq(makeAtomNullQuery());
  CHECK(aeq->Match(m.getAtomWithIdx(0)));
  CHECK(aeq->Match(m.getAtomWithIdx(1)));
  CHECK(aeq->Match(m.getAtomWithIdx(5)));
  CHECK(aeq->Match(m.getAtomWithIdx(6)));
}

TEST_CASE("query atom bond products") {
  Mol m;

  Atom a(6);
  m.addAtom(&a);
  m.addAtom(&a);
  m.addAtom(&a);
  m.addAtom(&a);
  m.addAtom(&a);
  m.addAtom(&a);
  m.addAtom(&a);
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
  CHECK(valenceProd == 1681);  // aromatic*aromatic = 41 * 41 = 1681
  valenceProd = queryAtomBondProduct(m.getAtomWithIdx(1));
  CHECK(valenceProd == 1681);

  valenceProd = queryAtomBondProduct(m.getAtomWithIdx(5));
  CHECK(valenceProd == 5043);

  valenceProd = queryAtomBondProduct(m.getAtomWithIdx(6));
  CHECK(valenceProd == 3);

  valenceProd = queryAtomAllBondProduct(m.getAtomWithIdx(6));
  CHECK(valenceProd == 81);
}

TEST_CASE("query-query matching") {
  // =====================================
  // ATOMS
  // =====================================
  {
    QueryAtom a1, a2;
    a1.setQuery(makeAtomNullQuery());
    a2.setQuery(makeAtomNullQuery());
    CHECK(a1.QueryMatch(&a2));
    CHECK(a2.QueryMatch(&a1));
  }

  {
    QueryAtom a1, a2(6);
    a1.setQuery(makeAtomNullQuery());
    CHECK(a1.QueryMatch(&a2));
    CHECK(a2.QueryMatch(&a1));
  }

  {
    QueryAtom a1(6), a2(8);
    CHECK(!a1.QueryMatch(&a2));
    CHECK(!a2.QueryMatch(&a1));
  }

  {
    QueryAtom a1(6), a2(6);
    CHECK(a1.QueryMatch(&a2));
    CHECK(a2.QueryMatch(&a1));
  }

  {
    QueryAtom a1, a2;
    a1.setQuery(makeAtomAromaticQuery());
    a2.setQuery(makeAtomAromaticQuery());
    CHECK(a1.QueryMatch(&a2));
    CHECK(a2.QueryMatch(&a1));
  }

  {
    QueryAtom a1, a2;
    a1.setQuery(makeAtomAromaticQuery());
    a2.setQuery(makeAtomAliphaticQuery());
    CHECK(!a1.QueryMatch(&a2));
    CHECK(!a2.QueryMatch(&a1));
  }

  {
    QueryAtom a1(6), a2(8);
    a1.expandQuery(makeAtomNumQuery(8), Queries::COMPOSITE_OR);
    a2.expandQuery(makeAtomNumQuery(9), Queries::COMPOSITE_OR);
    CHECK(a1.QueryMatch(&a2));
    CHECK(a2.QueryMatch(&a1));
  }

  {
    QueryAtom a1(6), a2(8);
    a1.expandQuery(makeAtomNumQuery(7), Queries::COMPOSITE_OR);
    a2.expandQuery(makeAtomNumQuery(9), Queries::COMPOSITE_OR);
    CHECK(!a1.QueryMatch(&a2));
    CHECK(!a2.QueryMatch(&a1));
  }

  {
    QueryAtom a1(6), a2(6);
    a1.expandQuery(makeAtomExplicitValenceQuery(3), Queries::COMPOSITE_AND);
    a2.expandQuery(makeAtomExplicitValenceQuery(3), Queries::COMPOSITE_AND);
    CHECK(a1.QueryMatch(&a2));
    CHECK(a2.QueryMatch(&a1));
  }

  {
    QueryAtom a1(6), a2(6);
    a1.expandQuery(makeAtomExplicitValenceQuery(3), Queries::COMPOSITE_AND);
    a2.expandQuery(makeAtomExplicitValenceQuery(4), Queries::COMPOSITE_AND);
    CHECK(!a1.QueryMatch(&a2));
    CHECK(!a2.QueryMatch(&a1));
  }

  {
    QueryAtom a1(6), a2(8);
    a1.expandQuery(makeAtomExplicitValenceQuery(3), Queries::COMPOSITE_AND);
    a2.expandQuery(makeAtomExplicitValenceQuery(3), Queries::COMPOSITE_AND);
    CHECK(!a1.QueryMatch(&a2));
    CHECK(!a2.QueryMatch(&a1));
  }

  {
    QueryAtom a1(6), a2(8);
    a1.expandQuery(makeAtomNumQuery(8), Queries::COMPOSITE_OR);
    CHECK(a1.QueryMatch(&a2));
    CHECK(a2.QueryMatch(&a1));
  }

  {
    QueryAtom a1;
    std::string excStr;
    std::unique_ptr<ATOM_EQUALS_QUERY> q(makeAtomNumQuery(7));
    try {
      a1.expandQuery(q.get(), Queries::COMPOSITE_OR);
    } catch (std::exception &e) {
      excStr = e.what();
    }
    CHECK(excStr == "Can't expand empty query");
  }

  // =====================================
  // BONDS
  // =====================================

  {
    QueryBond a1, a2;
    a1.setQuery(makeBondNullQuery());
    a2.setQuery(makeBondNullQuery());
    CHECK(a1.Match(&a2));
    CHECK(a2.QueryMatch(&a1));
  }

  {
    QueryBond a1, a2(Bond::SINGLE);
    a1.setQuery(makeBondNullQuery());
    CHECK(a1.QueryMatch(&a2));
    CHECK(a2.QueryMatch(&a1));
  }

  {
    QueryBond a1(Bond::SINGLE), a2(Bond::SINGLE);
    CHECK(a1.QueryMatch(&a2));
    CHECK(a2.QueryMatch(&a1));
  }

  {
    QueryBond a1(Bond::SINGLE), a2(Bond::DOUBLE);
    CHECK(!a1.QueryMatch(&a2));
    CHECK(!a2.QueryMatch(&a1));
  }

  {
    QueryBond a1(Bond::SINGLE), a2(Bond::DOUBLE);
    a1.expandQuery(makeBondOrderEqualsQuery(Bond::DOUBLE),
                   Queries::COMPOSITE_OR);
    a2.expandQuery(makeBondOrderEqualsQuery(Bond::AROMATIC),
                   Queries::COMPOSITE_OR);
    CHECK(a1.QueryMatch(&a2));
    CHECK(a2.QueryMatch(&a1));
  }

  {
    QueryBond a1(Bond::SINGLE), a2(Bond::DOUBLE);
    a1.expandQuery(makeBondOrderEqualsQuery(Bond::TRIPLE),
                   Queries::COMPOSITE_OR);
    a2.expandQuery(makeBondOrderEqualsQuery(Bond::AROMATIC),
                   Queries::COMPOSITE_OR);
    CHECK(!a1.QueryMatch(&a2));
    CHECK(!a2.QueryMatch(&a1));
  }

  {
    QueryBond a1(Bond::SINGLE), a2(Bond::DOUBLE);
    a1.expandQuery(makeBondMinRingSizeQuery(4), Queries::COMPOSITE_OR);
    a2.expandQuery(makeBondMinRingSizeQuery(4), Queries::COMPOSITE_OR);
    CHECK(a1.QueryMatch(&a2));
    CHECK(a2.QueryMatch(&a1));
  }

  {
    QueryBond a1(Bond::SINGLE), a2(Bond::DOUBLE);
    a1.expandQuery(makeBondMinRingSizeQuery(4), Queries::COMPOSITE_AND);
    a2.expandQuery(makeBondMinRingSizeQuery(4), Queries::COMPOSITE_AND);
    CHECK(!a1.QueryMatch(&a2));
    CHECK(!a2.QueryMatch(&a1));
  }

  {
    QueryBond a1(Bond::SINGLE), a2(Bond::SINGLE);
    a1.expandQuery(makeBondMinRingSizeQuery(5), Queries::COMPOSITE_AND);
    a2.expandQuery(makeBondMinRingSizeQuery(4), Queries::COMPOSITE_AND);
    CHECK(!a1.QueryMatch(&a2));
    CHECK(!a2.QueryMatch(&a1));
  }

  {
    QueryBond a1(Bond::SINGLE), a2(Bond::AROMATIC);
    a1.expandQuery(makeBondOrderEqualsQuery(Bond::AROMATIC),
                   Queries::COMPOSITE_OR);
    CHECK(a1.QueryMatch(&a2));
    CHECK(a2.QueryMatch(&a1));
  }
}

TEST_CASE("query atom mass") {
  Mol m;

  Atom a(6);

  int massVal;
  massVal = queryAtomMass(&a);
  CHECK(massVal ==
        static_cast<int>(std::round(12.011 * massIntegerConversionFactor)));

  a.setIsotope(13);
  massVal = queryAtomMass(&a);
  CHECK(massVal ==
        static_cast<int>(std::round(13.003 * massIntegerConversionFactor)));
}

TEST_CASE("ring query matches") {
  std::unique_ptr<RWMol> m(SmartsToMol("[C]1-[C]-[C]1"));
  REQUIRE(m);
  MolOps::findSSSR(*m);

  std::unique_ptr<RWMol> q(SmartsToMol("[R]"));
  REQUIRE(q);

  std::vector<MatchVectType> mvv;
  CHECK(SubstructMatch(*m, *q, mvv));
  CHECK(mvv.size() == 3);
  CHECK(mvv[0].size() == 1);
}

TEST_CASE("qualified atom queries") {
  std::unique_ptr<RWMol> m(SmilesToMol("CNO"));
  REQUIRE(m);

  {
    QueryAtom qA;
    qA.setQuery(makeAtomNumQuery<ATOM_GREATER_QUERY>(7, "test"));
    CHECK(qA.Match(m->getAtomWithIdx(0)));
    CHECK(!qA.Match(m->getAtomWithIdx(1)));
    CHECK(!qA.Match(m->getAtomWithIdx(2)));
  }
  {
    QueryAtom qA;
    qA.setQuery(makeAtomNumQuery<ATOM_GREATEREQUAL_QUERY>(7, "test"));
    CHECK(qA.Match(m->getAtomWithIdx(0)));
    CHECK(qA.Match(m->getAtomWithIdx(1)));
    CHECK(!qA.Match(m->getAtomWithIdx(2)));
  }
  {
    QueryAtom qA;
    qA.setQuery(makeAtomNumQuery<ATOM_LESS_QUERY>(7, "test"));
    CHECK(!qA.Match(m->getAtomWithIdx(0)));
    CHECK(!qA.Match(m->getAtomWithIdx(1)));
    CHECK(qA.Match(m->getAtomWithIdx(2)));
  }
  {
    QueryAtom qA;
    qA.setQuery(makeAtomNumQuery<ATOM_LESSEQUAL_QUERY>(7, "test"));
    CHECK(!qA.Match(m->getAtomWithIdx(0)));
    CHECK(qA.Match(m->getAtomWithIdx(1)));
    CHECK(qA.Match(m->getAtomWithIdx(2)));
  }
}

TEST_CASE("radicals in atom matching") {
  Atom a1(6);
  Atom a2(6);

  CHECK(a1.Match(&a2));
  CHECK(a2.Match(&a1));
  a1.setNumRadicalElectrons(2);
  CHECK(!a1.Match(&a2));
  CHECK(a2.Match(&a1));
  a2.setNumRadicalElectrons(2);
  CHECK(a1.Match(&a2));
  CHECK(a2.Match(&a1));
  a2.setNumRadicalElectrons(3);
  CHECK(!a1.Match(&a2));
  CHECK(!a2.Match(&a1));
}

TEST_CASE("atom has-property queries") {
  QueryAtom qA;
  qA.setQuery(makeHasPropQuery<Atom>("foo"));
  Atom a1(6);
  CHECK(!qA.Match(&a1));
  a1.setProp<int>("foo", 1);
  CHECK(qA.Match(&a1));
}

TEST_CASE("atom and bond property-value queries") {
  {
    QueryAtom qA;
    qA.setQuery(makePropQuery<Atom, int>("foo", 2));
    Atom a1(6);
    CHECK(!qA.Match(&a1));
    a1.setProp<int>("foo", 1);
    CHECK(!qA.Match(&a1));
    a1.setProp<int>("foo", 2);
    CHECK(qA.Match(&a1));

    a1.clearProp("foo");
    a1.setProp<double>("foo", 2);
    CHECK(!qA.Match(&a1));
  }

  {
    QueryAtom qA;
    qA.setQuery(makePropQuery<Atom, std::string>("foo", "bar"));
    Atom a1(6);
    CHECK(!qA.Match(&a1));
    a1.setProp<std::string>("foo", "bar");
    CHECK(qA.Match(&a1));
  }

  {
    QueryBond qA;
    qA.setQuery(makePropQuery<Bond, int>("foo", 2));
    Bond a1;
    CHECK(!qA.Match(&a1));
    a1.setProp<int>("foo", 1);
    CHECK(!qA.Match(&a1));
    a1.setProp<int>("foo", 2);
    CHECK(qA.Match(&a1));

    a1.clearProp("foo");
    a1.setProp<double>("foo", 2);
    CHECK(!qA.Match(&a1));
  }

  {
    QueryBond qA;
    qA.setQuery(makePropQuery<Bond, std::string>("foo", "bar"));
    Bond a1;
    CHECK(!qA.Match(&a1));
    a1.setProp<std::string>("foo", "bar");
    CHECK(qA.Match(&a1));
  }
}

TEST_CASE("double property-value queries") {
  {
    QueryAtom qA;
    qA.setQuery(makePropQuery<Atom, double>("foo", 2));
    Atom a1(6);
    CHECK(!qA.Match(&a1));
    a1.setProp<double>("foo", 1);
    CHECK(!qA.Match(&a1));
    a1.setProp<double>("foo", 2);
    CHECK(qA.Match(&a1));
  }
  {
    QueryBond qA;
    qA.setQuery(makePropQuery<Bond, double>("foo", 2));
    Bond a1;
    CHECK(!qA.Match(&a1));
    a1.setProp<double>("foo", 1);
    CHECK(!qA.Match(&a1));
    a1.setProp<double>("foo", 2);
    CHECK(qA.Match(&a1));
  }
}

TEST_CASE("additional atom queries") {
  {  // radicals
    QueryAtom qA;
    qA.setQuery(makeAtomNumRadicalElectronsQuery(1));
    Atom a1(6);
    CHECK(!qA.Match(&a1));
    a1.setNumRadicalElectrons(1);
    CHECK(qA.Match(&a1));
    a1.setNumRadicalElectrons(2);
    CHECK(!qA.Match(&a1));
    qA.getQuery()->setNegation(true);
    CHECK(qA.Match(&a1));
    a1.setNumRadicalElectrons(0);
    CHECK(qA.Match(&a1));
    a1.setNumRadicalElectrons(1);
    CHECK(!qA.Match(&a1));
  }
  {  // chiral tags
    QueryAtom qA;
    qA.setQuery(makeAtomHasChiralTagQuery());
    Atom a1(6);
    CHECK(!qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
    CHECK(qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
    CHECK(qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_OTHER);
    CHECK(qA.Match(&a1));
    qA.getQuery()->setNegation(true);
    a1.setChiralTag(Atom::CHI_UNSPECIFIED);
    CHECK(qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
    CHECK(!qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
    CHECK(!qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_OTHER);
    CHECK(!qA.Match(&a1));
  }

  {  // missing chiral tags
    QueryAtom qA;
    qA.setQuery(makeAtomMissingChiralTagQuery());
    Atom a1(6);
    CHECK(!qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
    CHECK(!qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
    CHECK(!qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_OTHER);
    CHECK(!qA.Match(&a1));

    a1.setChiralTag(Atom::CHI_UNSPECIFIED);
    a1.setProp(common_properties::_ChiralityPossible, 1);
    CHECK(qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
    CHECK(!qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
    CHECK(!qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_OTHER);
    CHECK(!qA.Match(&a1));

    qA.getQuery()->setNegation(true);
    a1.clearProp(common_properties::_ChiralityPossible);
    a1.setChiralTag(Atom::CHI_UNSPECIFIED);
    CHECK(qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
    CHECK(qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
    CHECK(qA.Match(&a1));
    a1.setChiralTag(Atom::CHI_OTHER);
    CHECK(qA.Match(&a1));

    a1.setChiralTag(Atom::CHI_UNSPECIFIED);
    a1.setProp(common_properties::_ChiralityPossible, 1);
    CHECK(!qA.Match(&a1));
  }
}
TEST_CASE("additional bond queries") {
  {  // stereo
    QueryBond qB;
    qB.setQuery(makeBondHasStereoQuery());
    Bond b1(Bond::DOUBLE);
    CHECK(!qB.Match(&b1));
    b1.setStereo(Bond::STEREOE);
    CHECK(qB.Match(&b1));
    b1.setStereo(Bond::STEREOANY);
    CHECK(qB.Match(&b1));
    qB.getQuery()->setNegation(true);
    b1.setStereo(Bond::STEREONONE);
    CHECK(qB.Match(&b1));
    b1.setStereo(Bond::STEREOE);
    CHECK(!qB.Match(&b1));
    b1.setStereo(Bond::STEREOANY);
    CHECK(!qB.Match(&b1));
  }
}
TEST_CASE("heteroatom neighbor queries") {
  std::unique_ptr<RWMol> m(SmilesToMol("CCNCOO"));
  REQUIRE(m);
  {
    QueryAtom qA;
    qA.setQuery(makeAtomNumHeteroatomNbrsQuery(1));
    CHECK(!qA.Match(m->getAtomWithIdx(0)));
    CHECK(qA.Match(m->getAtomWithIdx(1)));
    CHECK(!qA.Match(m->getAtomWithIdx(2)));
    CHECK(!qA.Match(m->getAtomWithIdx(3)));
    CHECK(qA.Match(m->getAtomWithIdx(5)));
  }
  {
    QueryAtom qA;
    qA.setQuery(
        makeAtomRangeQuery(0, 3, true, true, queryAtomNumHeteroatomNbrs));
    CHECK(!qA.Match(m->getAtomWithIdx(0)));
    CHECK(qA.Match(m->getAtomWithIdx(1)));
    CHECK(!qA.Match(m->getAtomWithIdx(2)));
    CHECK(qA.Match(m->getAtomWithIdx(3)));
    CHECK(qA.Match(m->getAtomWithIdx(5)));
  }
  {
    QueryAtom qA;
    qA.setQuery(
        makeAtomRangeQuery(1, 2, false, false, queryAtomNumHeteroatomNbrs));
    CHECK(!qA.Match(m->getAtomWithIdx(0)));
    CHECK(qA.Match(m->getAtomWithIdx(1)));
    CHECK(!qA.Match(m->getAtomWithIdx(2)));
    CHECK(qA.Match(m->getAtomWithIdx(3)));
    CHECK(qA.Match(m->getAtomWithIdx(5)));
  }
  {
    QueryAtom qA;
    qA.setQuery(
        makeAtomRangeQuery(0, 2, true, false, queryAtomNumHeteroatomNbrs));
    CHECK(!qA.Match(m->getAtomWithIdx(0)));
    CHECK(qA.Match(m->getAtomWithIdx(1)));
    CHECK(!qA.Match(m->getAtomWithIdx(2)));
    CHECK(qA.Match(m->getAtomWithIdx(3)));
    CHECK(qA.Match(m->getAtomWithIdx(5)));
  }
}

TEST_CASE("atom type queries") {
  std::unique_ptr<RWMol> m(SmilesToMol("CCc1ccccc1"));
  REQUIRE(m);
  {
    QueryAtom qA1;
    qA1.setQuery(makeAtomTypeQuery(6, true));
    QueryAtom qA2;
    qA2.setQuery(makeAtomTypeQuery(6, false));
    QueryAtom qA3;
    qA3.setQuery(makeAtomTypeQuery(7, true));
    CHECK(!qA1.Match(m->getAtomWithIdx(0)));
    CHECK(qA2.Match(m->getAtomWithIdx(0)));
    CHECK(!qA3.Match(m->getAtomWithIdx(0)));
    CHECK(qA1.Match(m->getAtomWithIdx(2)));
    CHECK(!qA2.Match(m->getAtomWithIdx(2)));
    CHECK(!qA3.Match(m->getAtomWithIdx(2)));
  }
}

TEST_CASE("dummy atoms in recursive queries") {
  auto m = "[*;$(CC)][*]"_smarts;
  REQUIRE(m);
  CHECK(m->getAtomWithIdx(0)->hasQuery());
  CHECK(m->getAtomWithIdx(1)->hasQuery());
  CHECK(isComplexQuery(m->getAtomWithIdx(0)));
  CHECK(!isComplexQuery(m->getAtomWithIdx(1)));
}

TEST_CASE("bond type queries") {
  {
    auto m = "C-C-@C@-C@CC-,=C!-C"_smarts;
    REQUIRE(m);
    for (const auto bond : m->bonds()) {
      bool hasQ = QueryOps::hasBondTypeQuery(*bond);
      if (bond->getIdx() != 3) {
        CHECK(hasQ);
      } else {
        CHECK(!hasQ);
      }
    }
    for (const auto bond : m->bonds()) {
      bool hasQ = QueryOps::hasComplexBondTypeQuery(*bond);
      if (bond->getIdx() < 4) {
        CHECK(!hasQ);
      } else {
        CHECK(hasQ);
      }
    }
  }
  {
    auto m = "CC-C"_smiles;
    REQUIRE(m);
    for (const auto bond : m->bonds()) {
      CHECK(!QueryOps::hasBondTypeQuery(*bond));
    }
  }
}

TEST_CASE("query-query matches with ring queries") {
  SECTION("atom basics") {
    std::vector<std::tuple<std::string, std::string, bool>> checks = {
        {"[R]", "[R]", true},
        {"[R]", "[!R]", false},
        {"[R]", "[R3]", false},
        {"[R]", "[!R3]", false}};
    for (const auto &[smarts1, smarts2, expected] : checks) {
      INFO(smarts1 << " " << smarts2);
      auto a1 = v2::SmilesParse::AtomFromSmarts(smarts1);
      REQUIRE(a1);
      CHECK(a1->hasQuery());
      auto a2 = v2::SmilesParse::AtomFromSmarts(smarts2);
      REQUIRE(a2);
      CHECK(a2->hasQuery());
      CHECK(dynamic_cast<QueryAtom *>(a1.get())->QueryMatch(
                dynamic_cast<QueryAtom *>(a2.get())) == expected);
      CHECK(dynamic_cast<QueryAtom *>(a2.get())->QueryMatch(
                dynamic_cast<QueryAtom *>(a1.get())) == expected);
    }
  }
  SECTION("atoms in molecules basics") {
    std::vector<std::tuple<std::string, std::string, size_t>> checks = {
        {"[R][R]", "[R][R]", 1},
        {"[R][R]", "[R][!R]", 0},
        {"[R][R]", "[!R][!R]", 0},
        //
        {"[!R][R]", "[R][R]", 0},
        {"[!R][R]", "[R][!R]", 1},
        {"[!R][R]", "[!R][!R]", 0},
        //
        {"[!R][!R]", "[R][R]", 0},
        {"[!R][!R]", "[R][!R]", 0},
        {"[!R][!R]", "[!R][!R]", 1}};
    for (const auto &[smarts1, smarts2, expected_size] : checks) {
      INFO(smarts1 << " " << smarts2);
      auto m1 = v2::SmilesParse::MolFromSmarts(smarts1);
      REQUIRE(m1);
      auto m2 = v2::SmilesParse::MolFromSmarts(smarts2);
      REQUIRE(m2);
      SubstructMatchParameters params;
      params.useQueryQueryMatches = true;
      CHECK(SubstructMatch(*m1, *m2, params).size() == expected_size);
    }
  }
  SECTION("bond basics") {
    std::vector<std::tuple<std::string, std::string, bool>> checks = {
        {"@", "@", true}, {"@", "!@", false}};
    for (const auto &[smarts1, smarts2, expected] : checks) {
      INFO(smarts1 << " " << smarts2);
      auto b1 = v2::SmilesParse::BondFromSmarts(smarts1);
      REQUIRE(b1);
      auto b2 = v2::SmilesParse::BondFromSmarts(smarts2);
      REQUIRE(b2);
      CHECK(dynamic_cast<QueryBond *>(b1.get())->QueryMatch(
                dynamic_cast<QueryBond *>(b2.get())) == expected);
      CHECK(dynamic_cast<QueryBond *>(b2.get())->QueryMatch(
                dynamic_cast<QueryBond *>(b1.get())) == expected);
    }
  }
  SECTION("bonds in molecules basics") {
    std::vector<std::tuple<std::string, std::string, size_t>> checks = {
        {"[R]@[R]", "[R]@[R]", 1},
        {"[R]@[R]", "[R]!@[R]", 0},
        {"[R]!@[R]", "[R]@[R]", 0},
        {"[R]!@[R]", "[R]!@[R]", 1}};
    for (const auto &[smarts1, smarts2, expected_size] : checks) {
      INFO(smarts1 << " " << smarts2);
      auto m1 = v2::SmilesParse::MolFromSmarts(smarts1);
      REQUIRE(m1);
      auto m2 = v2::SmilesParse::MolFromSmarts(smarts2);
      REQUIRE(m2);
      SubstructMatchParameters params;
      params.useQueryQueryMatches = true;
      CHECK(SubstructMatch(*m1, *m2, params).size() == expected_size);
    }
  }
}