//
//  Copyright (C) 2001-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//

#include <catch2/catch_all.hpp>

// RD bits
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/Chirality.h>
#include "SubstructMatch.h"
#include "SubstructUtils.h"

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/test_fixtures.h>

#include "vf2.hpp"

using namespace RDKit;

TEST_CASE("test1", "[substruct]") {
  bool updateLabel = true;
  bool takeOwnership = true;
  std::unique_ptr<RWMol> m = std::make_unique<RWMol>();
  m->addAtom(new Atom(8), updateLabel, takeOwnership);
  m->addAtom(new Atom(6), updateLabel, takeOwnership);
  m->addAtom(new Atom(6), updateLabel, takeOwnership);
  m->addBond(0, 1, Bond::SINGLE);
  m->addBond(1, 2, Bond::SINGLE);

  std::unique_ptr<RWMol> q1 = std::make_unique<RWMol>();
  q1->addAtom(new QueryAtom(6), updateLabel, takeOwnership);
  q1->addAtom(new QueryAtom(6), updateLabel, takeOwnership);
  q1->addBond(0, 1, Bond::SINGLE);

  std::vector<MatchVectType> matches;
  auto n = SubstructMatch(*m, *q1, matches, false);
  REQUIRE(n == 2);
  REQUIRE(matches.size() == n);
  REQUIRE(matches[0].size() == 2);
  CHECK(matches[0][0].first == 0);
  CHECK((matches[0][0].second == 1 || matches[0][0].second == 2));
  CHECK(matches[0][1].first == 1);
  CHECK(matches[0][1].second != matches[0][0].second);

  CHECK((matches[1][1].second == 1 || matches[0][1].second == 2));
  CHECK(matches[1][0].first == 0);
  CHECK((matches[1][0].second == 1 || matches[1][0].second == 2));
  CHECK(matches[1][0].second != matches[0][0].second);
  CHECK(matches[1][0].second == matches[0][1].second);
  CHECK(matches[1][1].first == 1);
  CHECK(matches[1][1].second != matches[1][0].second);
  CHECK(matches[1][1].second == matches[0][0].second);

  n = SubstructMatch(*m, *q1, matches, true);
  REQUIRE(n == 1);
  REQUIRE(matches.size() == n);
  REQUIRE(matches[0].size() == 2);
  CHECK(matches[0][0].first == 0);
  CHECK((matches[0][0].second == 1 || matches[0][0].second == 2));
  CHECK(matches[0][1].first == 1);
  CHECK(matches[0][1].second != matches[0][0].second);
  CHECK((matches[0][1].second == 1 || matches[0][1].second == 2));

  MatchVectType matchV;
  REQUIRE(SubstructMatch(*m, *q1, matchV));
  REQUIRE(matchV.size() == 2);
  // make sure we reset the match vectors.
  // build a query we won't match:
  q1->addAtom(new QueryAtom(6), updateLabel, takeOwnership);
  q1->addBond(1, 2, Bond::SINGLE);
  q1->addAtom(new QueryAtom(6), updateLabel, takeOwnership);
  q1->addBond(2, 3, Bond::SINGLE);

  CHECK(!SubstructMatch(*m, *q1, matchV));
  CHECK(matchV.size() == 0);

  n = SubstructMatch(*m, *q1, matches, false);
  CHECK(n == 0);
  CHECK(matches.size() == 0);
}

TEST_CASE("test2", "[substruct]") {
  std::unique_ptr<RWMol> m = std::make_unique<RWMol>();
  bool updateLabel = true;
  bool takeOwnership = true;
  m->addAtom(new Atom(6), updateLabel, takeOwnership);
  m->addAtom(new Atom(6), updateLabel, takeOwnership);
  m->addAtom(new Atom(8), updateLabel, takeOwnership);
  m->addBond(0, 1, Bond::SINGLE);
  m->addBond(1, 2, Bond::SINGLE);

  std::unique_ptr<RWMol> q1 = std::make_unique<RWMol>();
  q1->addAtom(new QueryAtom(6), updateLabel, takeOwnership);
  q1->addAtom(new QueryAtom(8), updateLabel, takeOwnership);
  q1->addBond(0, 1, Bond::SINGLE);

  MatchVectType matchV;
  auto n = SubstructMatch(*m, *q1, matchV);
  CHECK(n);
  CHECK(matchV.size() == 2);
  CHECK(matchV[0].first == 0);
  CHECK(matchV[0].second == 1);
  CHECK(matchV[1].first == 1);
  CHECK(matchV[1].second == 2);

  std::vector<MatchVectType> matches;
  n = SubstructMatch(*m, *q1, matches, false);
  REQUIRE(n == 1);
  REQUIRE(matches.size() == n);
  REQUIRE(matches[0].size() == 2);
  n = SubstructMatch(*m, *q1, matches, true);
  REQUIRE(n == 1);
  REQUIRE(matches.size() == n);
  REQUIRE(matches[0].size() == 2);
  REQUIRE(SubstructMatch(*m, *q1, matchV));
  REQUIRE(matchV.size() == 2);
  m.reset(new RWMol());
  m->addAtom(new Atom(6), updateLabel, takeOwnership);
  m->addAtom(new Atom(6), updateLabel, takeOwnership);
  m->addAtom(new Atom(8), updateLabel, takeOwnership);
  m->addBond(0, 1, Bond::SINGLE);
  m->addBond(1, 2, Bond::DOUBLE);

  matches.clear();
  n = SubstructMatch(*m, *q1, matches, false);
  REQUIRE(n == 0);
  REQUIRE(matches.size() == n);
  n = SubstructMatch(*m, *q1, matches, true);
  REQUIRE(n == 0);
  REQUIRE(matches.size() == n);
  REQUIRE(!SubstructMatch(*m, *q1, matchV));
}

TEST_CASE("test3", "[substruct]") {
  std::unique_ptr<RWMol> m = std::make_unique<RWMol>();
  bool updateLabel = true;
  bool takeOwnership = true;
  m->addAtom(new Atom(6), updateLabel, takeOwnership);
  m->addAtom(new Atom(6), updateLabel, takeOwnership);
  m->addAtom(new Atom(8), updateLabel, takeOwnership);
  m->addBond(0, 1, Bond::SINGLE);
  m->addBond(1, 2, Bond::SINGLE);

  std::unique_ptr<RWMol> q1 = std::make_unique<RWMol>();
  q1->addAtom(new QueryAtom(6), updateLabel, takeOwnership);
  q1->addAtom(new QueryAtom(8), updateLabel, takeOwnership);
  q1->addBond(0, 1, Bond::UNSPECIFIED);

  std::vector<MatchVectType> matches;
  auto n = SubstructMatch(*m, *q1, matches, false);
  REQUIRE(n == 1);
  REQUIRE(matches.size() == n);
  REQUIRE(matches[0].size() == 2);
  n = SubstructMatch(*m, *q1, matches, true);
  REQUIRE(n == 1);
  REQUIRE(matches.size() == n);
  REQUIRE(matches[0].size() == 2);

  MatchVectType matchV;
  REQUIRE(SubstructMatch(*m, *q1, matchV));
  REQUIRE(matchV.size() == 2);
  m = std::make_unique<RWMol>();
  m->addAtom(new Atom(6), updateLabel, takeOwnership);
  m->addAtom(new Atom(6), updateLabel, takeOwnership);
  m->addAtom(new Atom(8), updateLabel, takeOwnership);
  m->addBond(0, 1, Bond::SINGLE);
  m->addBond(1, 2, Bond::DOUBLE);

  matches.clear();
  n = SubstructMatch(*m, *q1, matches, false);
  REQUIRE(n == 1);
  REQUIRE(matches.size() == n);
  REQUIRE(matches[0].size() == 2);
  n = SubstructMatch(*m, *q1, matches, true);
  REQUIRE(n == 1);
  REQUIRE(matches.size() == n);
  REQUIRE(matches[0].size() == 2);
  REQUIRE(SubstructMatch(*m, *q1, matchV));
  REQUIRE(matchV.size() == 2);
  q1 = std::make_unique<RWMol>();
  q1->addAtom(new QueryAtom(6), updateLabel, takeOwnership);
  q1->addAtom(new QueryAtom(6), updateLabel, takeOwnership);
  q1->addBond(0, 1, Bond::UNSPECIFIED);
  n = SubstructMatch(*m, *q1, matches, false);
  CHECK(n == 2);
  CHECK(matches.size() == n);
  CHECK(matches[0].size() == 2);
  CHECK(matches[1].size() == 2);
  CHECK(matches[0][0].second != matches[1][0].second);
  CHECK(matches[0][1].second != matches[1][1].second);
  n = SubstructMatch(*m, *q1, matches, true);
  CHECK(n == 1);
  CHECK(matches.size() == n);
}

TEST_CASE("test4", "[substruct]") {
  bool updateLabel = true;
  bool takeOwnership = true;

  std::unique_ptr<Atom> a6 = std::make_unique<Atom>(6);
  std::unique_ptr<Atom> a8 = std::make_unique<Atom>(8);
  std::unique_ptr<RWMol> m = std::make_unique<RWMol>();
  m->addAtom(a6.get());
  m->addAtom(a6.get());
  m->addAtom(a8.get());
  m->addAtom(a6.get());
  m->addAtom(a6.get());
  m->addBond(1, 0, Bond::SINGLE);
  m->addBond(1, 2, Bond::SINGLE);
  m->addBond(1, 3, Bond::SINGLE);
  m->addBond(2, 4, Bond::SINGLE);

  // this will be the recursive query
  std::unique_ptr<RWMol> q1 = std::make_unique<RWMol>();
  q1->addAtom(new QueryAtom(6), updateLabel, takeOwnership);
  q1->addAtom(new QueryAtom(8), updateLabel, takeOwnership);
  q1->addBond(0, 1, Bond::UNSPECIFIED);

  // here's the main query
  std::unique_ptr<RWMol> q2 = std::make_unique<RWMol>();
  auto *rsq = new RecursiveStructureQuery(q1.release());
  auto *qA = new QueryAtom(6);
  qA->expandQuery(rsq, Queries::COMPOSITE_AND);
  // std::cout << "post expand: " << qA->getQuery() << std::endl;
  q2->addAtom(qA, true, true);
  // std::cout << "mol: " << q2->getAtomWithIdx(0)->getQuery() << std::endl;
  q2->addAtom(new QueryAtom(6), true, true);
  q2->addBond(0, 1, Bond::UNSPECIFIED);

  MatchVectType matchV;
  bool found = SubstructMatch(*m, *q2, matchV);
  REQUIRE(found);
  REQUIRE(matchV.size() == 2);
  CHECK(matchV[0].first == 0);
  CHECK(matchV[0].second == 1);
  CHECK(matchV[1].first == 1);
  CHECK((matchV[1].second == 0 || matchV[1].second == 3));
  std::vector<MatchVectType> matches;
  auto n = SubstructMatch(*m, *q2, matches, true);
  CHECK(n == 2);
  CHECK(matches.size() == (size_t)n);
  CHECK(matches[0].size() == 2);
  CHECK(matches[1].size() == 2);
  CHECK(matches[0][0].second == matches[1][0].second);
  CHECK(matches[0][1].second != matches[1][1].second);
}

TEST_CASE("test5", "[substruct]") {
  auto a6 = std::make_unique<Atom>(6);
  auto a8 = std::make_unique<Atom>(8);
  // CC(OC)C
  auto m = std::make_unique<RWMol>();
  m->addAtom(a6.get());
  m->addAtom(a6.get());
  m->addAtom(a8.get());
  m->addAtom(a6.get());
  m->addAtom(a6.get());
  m->addBond(0, 1, Bond::SINGLE);
  m->addBond(1, 2, Bond::SINGLE);
  m->addBond(1, 4, Bond::SINGLE);
  m->addBond(2, 3, Bond::SINGLE);

  // this will be the recursive query
  bool updateLabel = true;
  bool takeOwnership = true;
  auto q1 = std::make_unique<RWMol>();
  q1->addAtom(new QueryAtom(6), updateLabel, takeOwnership);
  q1->addAtom(new QueryAtom(8), updateLabel, takeOwnership);
  q1->addBond(0, 1, Bond::UNSPECIFIED);

  // here's the main query
  auto q2 = std::make_unique<RWMol>();
  auto qA = std::make_unique<QueryAtom>();
  auto *rsq = new RecursiveStructureQuery(q1.release());
  qA->setQuery(rsq);
  q2->addAtom(qA.release(), true, true);
  q2->addAtom(new QueryAtom(6), true, true);
  q2->addBond(0, 1, Bond::UNSPECIFIED);

  MatchVectType matchV;
  bool found = SubstructMatch(*m, *q2, matchV);
  REQUIRE(found);
  REQUIRE(matchV.size() == 2);
  std::vector<MatchVectType> matches;
  auto n = SubstructMatch(*m, *q2, matches, true);
  REQUIRE(n == 2);
  REQUIRE(matches[0].size() == 2);
}
TEST_CASE("test5QueryRoot", "[substruct]") {
  auto a6 = std::unique_ptr<Atom>(new Atom(6));
  auto a8 = std::unique_ptr<Atom>(new Atom(8));
  // CC(OC)C
  auto m = std::make_unique<RWMol>();
  m->addAtom(a6.get());
  m->addAtom(a6.get());
  m->addAtom(a8.get());
  m->addAtom(a6.get());
  m->addAtom(a6.get());
  m->addBond(0, 1, Bond::SINGLE);
  m->addBond(1, 2, Bond::SINGLE);
  m->addBond(1, 4, Bond::SINGLE);
  m->addBond(2, 3, Bond::SINGLE);

  // this will be the recursive query
  bool updateLabel = true;
  bool takeOwnership = true;
  auto q1 = std::make_unique<RWMol>();
  q1->addAtom(new QueryAtom(8), updateLabel, takeOwnership);
  q1->addAtom(new QueryAtom(6), updateLabel, takeOwnership);
  q1->addBond(0, 1, Bond::UNSPECIFIED);
  q1->setProp(common_properties::_queryRootAtom, 1);

  // here's the main query
  auto q2 = std::make_unique<RWMol>();
  auto qA = std::make_unique<QueryAtom>();
  auto *rsq = new RecursiveStructureQuery(q1.release());
  qA->setQuery(rsq);
  q2->addAtom(qA.release(), true, true);
  q2->addAtom(new QueryAtom(6), true, true);
  q2->addBond(0, 1, Bond::UNSPECIFIED);

  MatchVectType matchV;
  bool found = SubstructMatch(*m, *q2, matchV);
  REQUIRE(found);
  REQUIRE(matchV.size() == 2);
  std::vector<MatchVectType> matches;
  auto n = SubstructMatch(*m, *q2, matches, true);
  REQUIRE(n == 2);
  REQUIRE(matches[0].size() == 2);
}

TEST_CASE("test6", "[substruct][Issue71]") {
  auto a6 = std::make_unique<Atom>(6);

  auto m = std::make_unique<RWMol>();
  m->addAtom(a6.get());
  m->addAtom(a6.get());
  m->addAtom(a6.get());
  m->addBond(0, 1, Bond::SINGLE);
  m->addBond(1, 2, Bond::SINGLE);
  m->addBond(0, 2, Bond::SINGLE);

  auto q1 = std::make_unique<RWMol>();
  bool updateLabel = true;
  bool takeOwnership = true;
  q1->addAtom(new QueryAtom(6), updateLabel, takeOwnership);
  q1->addAtom(new QueryAtom(6), updateLabel, takeOwnership);
  q1->addAtom(new QueryAtom(6), updateLabel, takeOwnership);
  q1->addBond(0, 1, Bond::UNSPECIFIED);
  q1->addBond(1, 2, Bond::UNSPECIFIED);

  MatchVectType matchV;
  bool found = SubstructMatch(*m, *q1, matchV);
  REQUIRE(found);
  REQUIRE(matchV.size() == 3);
  std::vector<MatchVectType> matches;
  auto n = SubstructMatch(*m, *q1, matches, true);
  REQUIRE(n == 1);
  REQUIRE(matches[0].size() == 3);
  // close the loop and try again (we should still match)
  q1->addBond(0, 2, Bond::UNSPECIFIED);
  found = SubstructMatch(*m, *q1, matchV);
  REQUIRE(found);
  REQUIRE(matchV.size() == 3);
  n = SubstructMatch(*m, *q1, matches, true);
  REQUIRE(n == 1);
  REQUIRE(matches[0].size() == 3);
}

TEST_CASE("test7", "[substruct][leak]") {
  auto a6 = std::make_unique<Atom>(6);

  std::unique_ptr<RWMol> m = std::make_unique<RWMol>();
  m->addAtom(a6.get());
  m->addAtom(a6.get());
  m->addAtom(a6.get());
  m->addBond(0, 1, Bond::SINGLE);
  m->addBond(1, 2, Bond::SINGLE);
  m->addBond(0, 2, Bond::SINGLE);

  std::unique_ptr<RWMol> q1 = std::make_unique<RWMol>();
  bool updateLabel = true;
  bool takeOwnership = true;
  q1->addAtom(new QueryAtom(6), updateLabel, takeOwnership);
  q1->addAtom(new QueryAtom(6), updateLabel, takeOwnership);
  q1->addAtom(new QueryAtom(6), updateLabel, takeOwnership);
  q1->addBond(0, 1, Bond::UNSPECIFIED);
  q1->addBond(1, 2, Bond::UNSPECIFIED);

  MatchVectType matchV;
  bool found = SubstructMatch(*m, *q1, matchV);
  REQUIRE(found);
  REQUIRE(matchV.size() == 3);
  std::vector<MatchVectType> matches;
  for (int i = 0; i < 30000; i++) {
    auto n = SubstructMatch(*m, *q1, matches, true, true);
    REQUIRE(n == 1);
    REQUIRE(matches[0].size() == 3);
  }
}

TEST_CASE("test9", "[substruct][chiral]") {
  auto a6 = std::make_unique<Atom>(6);

  auto m = std::make_unique<RWMol>();
  m->addAtom(a6.get());
  bool updateLabel = true;
  bool takeOwnership = true;
  m->addAtom(new Atom(6), updateLabel, takeOwnership);
  m->addAtom(new Atom(7), updateLabel, takeOwnership);
  m->addAtom(new Atom(8), updateLabel, takeOwnership);
  m->addAtom(new Atom(9), updateLabel, takeOwnership);
  m->addBond(0, 1, Bond::SINGLE);
  m->addBond(0, 2, Bond::SINGLE);
  m->addBond(0, 3, Bond::SINGLE);
  m->addBond(0, 4, Bond::SINGLE);
  m->getAtomWithIdx(0)->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);

  auto q1 = std::make_unique<RWMol>();
  q1->addAtom(a6.get());
  q1->addAtom(new Atom(6));
  q1->addAtom(new Atom(7));
  q1->addAtom(new Atom(8));
  q1->addAtom(new Atom(9));
  q1->addBond(0, 1, Bond::SINGLE);
  q1->addBond(0, 2, Bond::SINGLE);
  q1->addBond(0, 3, Bond::SINGLE);
  q1->addBond(0, 4, Bond::SINGLE);
  q1->getAtomWithIdx(0)->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);

  MolOps::sanitizeMol(*m);
  MolOps::assignStereochemistry(*m);
  MolOps::sanitizeMol(*q1);
  MolOps::assignStereochemistry(*q1);

  // test with default options (no chirality):
  MatchVectType matchV;
  auto found = SubstructMatch(*m, *q1, matchV);
  CHECK(found);
  std::vector<MatchVectType> matches;
  auto n = SubstructMatch(*m, *q1, matches, true);
  CHECK(n == 1);

  // test with chirality
  found = SubstructMatch(*m, *q1, matchV, true, true);
  CHECK(!found);
  n = SubstructMatch(*m, *q1, matches, true, true, true);
  CHECK(n == 0);

  // self matches:
  found = SubstructMatch(*m, *m, matchV, true, true);
  CHECK(found);
  n = SubstructMatch(*m, *m, matches, true, true, true);
  CHECK(n == 1);
  found = SubstructMatch(*q1, *q1, matchV, true, true);
  CHECK(found);
  n = SubstructMatch(*q1, *q1, matches, true, true, true);
  CHECK(n == 1);
}

TEST_CASE("testRecursiveSerialNumbers", "[substruct][recursive]") {
  auto a6 = std::make_unique<Atom>(6);
  auto a8 = std::make_unique<Atom>(8);
  auto m = std::make_unique<RWMol>();
  m->addAtom(a6.get());
  m->addAtom(a6.get());
  m->addAtom(a8.get());
  m->addAtom(a6.get());
  m->addAtom(a6.get());
  m->addBond(1, 0, Bond::SINGLE);
  m->addBond(1, 2, Bond::SINGLE);
  m->addBond(1, 3, Bond::SINGLE);
  m->addBond(2, 4, Bond::SINGLE);

  {
    // this will be the recursive query
    auto q1 = std::make_unique<RWMol>();
    bool updateLabel = true;
    bool takeOwnership = true;
    q1->addAtom(new QueryAtom(6), updateLabel, takeOwnership);
    q1->addAtom(new QueryAtom(8), updateLabel, takeOwnership);
    q1->addBond(0, 1, Bond::UNSPECIFIED);

    // here's the main query
    auto q2 = std::make_unique<RWMol>();
    auto qA = std::make_unique<QueryAtom>(6);
    auto rsq = std::make_unique<RecursiveStructureQuery>(new RWMol(*q1), 1);
    qA->expandQuery(rsq.release(), Queries::COMPOSITE_AND);
    // std::cout << "post expand: " << qA->getQuery() << std::endl;
    q2->addAtom(qA.release(), true, true);
    // std::cout << "mol: " << q2->getAtomWithIdx(0)->getQuery() << std::endl;
    q2->addAtom(new QueryAtom(8), true, true);
    q2->addBond(0, 1, Bond::UNSPECIFIED);

    qA.reset(new QueryAtom(6));
    rsq.reset(new RecursiveStructureQuery(new RWMol(*q1), 1));
    qA->expandQuery(rsq.release(), Queries::COMPOSITE_AND);
    q2->addAtom(qA.release(), true, true);
    q2->addBond(1, 2, Bond::UNSPECIFIED);

    MatchVectType matchV;
    bool found = SubstructMatch(*m, *q2, matchV);
    REQUIRE(found);
    REQUIRE(matchV.size() == 3);
    std::vector<MatchVectType> matches;
    auto n = SubstructMatch(*m, *q2, matches, true);
    CHECK(n == 1);
    CHECK(matches.size() == 1);
    CHECK(matches[0].size() == 3);
  }
}

#ifdef RDK_TEST_MULTITHREADED
#include <RDGeneral/BoostStartInclude.h>
#include <thread>
#include <future>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/BoostEndInclude.h>
namespace {
void runblock(const std::vector<std::unique_ptr<ROMol>> *mols,
              const ROMol *query, const boost::dynamic_bitset<> &hits,
              unsigned int count, unsigned int idx) {
  for (unsigned int j = 0; j < 100; j++) {
    for (unsigned int i = 0; i < mols->size(); ++i) {
      if (i % count != idx) {
        continue;
      }
      const auto *mol = mols->at(i).get();
      MatchVectType matchV;
      bool found = SubstructMatch(*mol, *query, matchV);

      CHECK(found == hits[i]);
    }
  }
}
}  // namespace
TEST_CASE("testMultiThread", "[substruct][multithread]") {
  std::string fName = getenv("RDBASE");
  fName += "/Data/NCI/first_200.props.sdf";
  SDMolSupplier suppl(fName);
  std::vector<std::unique_ptr<ROMol>> mols;
  while (!suppl.atEnd() && mols.size() < 100) {
    ROMol *mol = nullptr;
    try {
      mol = suppl.next();
    } catch (...) {
      continue;
    }
    if (!mol) {
      continue;
    }
    mols.emplace_back(mol);
  }
  std::vector<std::future<void>> tg;
  auto query = v2::SmilesParse::MolFromSmarts("[#6;$([#6]([#6])[!#6])]");
  boost::dynamic_bitset<> hits(mols.size());
  for (unsigned int i = 0; i < mols.size(); ++i) {
    MatchVectType matchV;
    hits[i] = SubstructMatch(*mols[i], *query, matchV);
  }
  unsigned int count = 4;
  for (unsigned int i = 0; i < count; ++i) {
    tg.emplace_back(std::async(std::launch::async, runblock, &mols, query.get(),
                               hits, count, i));
  }
  for (auto &fut : tg) {
    fut.get();
  }
  tg.clear();

  query.reset(SmartsToMol("[#6]([#6])[!#6]"));
  for (unsigned int i = 0; i < mols.size(); ++i) {
    MatchVectType matchV;
    hits[i] = SubstructMatch(*mols[i], *query, matchV);
  }
  for (unsigned int i = 0; i < count; ++i) {
    tg.emplace_back(std::async(std::launch::async, runblock, &mols, query.get(),
                               hits, count, i));
  }
  for (auto &fut : tg) {
    fut.get();
  }
  tg.clear();

  query.reset(SmartsToMol("[$([O,S]-[!$(*=O)])]"));
  for (unsigned int i = 0; i < mols.size(); ++i) {
    MatchVectType matchV;
    hits[i] = SubstructMatch(*mols[i], *query, matchV);
  }
  for (unsigned int i = 0; i < count; ++i) {
    tg.emplace_back(std::async(std::launch::async, runblock, &mols, query.get(),
                               hits, count, i));
  }
  for (auto &fut : tg) {
    fut.get();
  }
}
#else
TEST_CASE("testMultiThread", "[substruct][multithread]") {}
#endif

TEST_CASE("testChiralMatch", "[substruct][chiral]") {
  {
    std::string qSmi = "Cl[C@](C)(F)Br";
    std::string mSmi = "Cl[C@](C)(F)Br";
    auto query = v2::SmilesParse::MolFromSmiles(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(matched);
  }
  {
    std::string qSmi = "Cl[C@](C)(F)Br";
    std::string mSmi = "Cl[C@@](C)(F)Br";
    auto query = v2::SmilesParse::MolFromSmiles(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(!matched);
  }
  {
    std::string qSmi = "Cl[C@](C)(F)Br";
    std::string mSmi = "Cl[C@@](F)(C)Br";
    auto query = v2::SmilesParse::MolFromSmiles(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(matched);
  }
  {
    std::string qSmi = "Cl[C@](C)(F)Br";
    std::string mSmi = "Cl[C@](F)(C)Br";
    auto query = v2::SmilesParse::MolFromSmiles(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(!matched);
  }
  {
    std::string qSmi = "Cl[C@](C)(F)Br";
    std::string mSmi = "Cl[C@@](Br)(C)F";
    auto query = v2::SmilesParse::MolFromSmiles(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(!matched);
  }
  {
    std::string qSmi = "Cl[C@](C)(F)Br";
    std::string mSmi = "Cl[C@](Br)(C)F";
    auto query = v2::SmilesParse::MolFromSmiles(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(matched);
  }

  {
    std::string qSmi = "C[C@](O)(F)Br";
    std::string mSmi = "O[C@](F)(Br)CC[C@](O)(F)Br";
    auto query = v2::SmilesParse::MolFromSmiles(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(matched);
  }

  {
    std::string qSmi = "C[C@](O)(F)Br";
    std::string mSmi = "O[C@](F)(Br)CC[C@@](O)(F)Br";
    auto query = v2::SmilesParse::MolFromSmiles(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(!matched);
  }

  {
    std::string qSmi = "C[C@](O)(F)Br";
    std::string mSmi = "O[C@](F)(Br)CC(C[C@](O)(F)Br)C[C@](O)(F)Br";
    auto query = v2::SmilesParse::MolFromSmiles(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    std::vector<MatchVectType> matches;
    int count = SubstructMatch(*mol, *query, matches, true, true, true);
    CHECK(count == 2);
  }

  {
    std::string qSmi = "C[C@](O)(F)Br";
    std::string mSmi = "O[C@@](F)(Br)CC(C[C@](O)(F)Br)C[C@](O)(F)Br";
    auto query = v2::SmilesParse::MolFromSmiles(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    std::vector<MatchVectType> matches;
    int count = SubstructMatch(*mol, *query, matches, true, true, true);
    CHECK(count == 3);
  }
  {
    std::string qSmi = "C[C@](O)(F)Br";
    std::string mSmi = "O[C@](F)(Br)CC(C[C@@](O)(F)Br)C[C@](O)(F)Br";
    auto query = v2::SmilesParse::MolFromSmiles(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    std::vector<MatchVectType> matches;
    int count = SubstructMatch(*mol, *query, matches, true, true, true);
    CHECK(count == 1);
  }
  {
    std::string qSmi = "C[C@](O)(F)Br";
    std::string mSmi = "O[C@](F)(Br)CC(C[C@@](O)(F)Br)C[C@@](O)(F)Br";
    auto query = v2::SmilesParse::MolFromSmiles(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    std::vector<MatchVectType> matches;
    // std::cerr<<"\n\n------------------------------------------\n"<<qSmi<<"
    // "<<mSmi<<"\n"<<std::endl;
    int count = SubstructMatch(*mol, *query, matches, true, true, true);
    // std::cerr<<"res: "<<count<<std::endl;
    CHECK(count == 0);
  }
  {
    std::string qSmi = "Cl[C@](*)(F)Br";
    std::string mSmi = "Cl[C@](C)(F)Br";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(matched);
  }
  {
    std::string qSmi = "Cl[C@@](*)(F)Br";
    std::string mSmi = "Cl[C@](C)(F)Br";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(!matched);
  }
  {
    std::string qSmi = "Cl[C@](*)(*)Br";
    std::string mSmi = "Cl[C@](C)(F)Br";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(matched);
  }
  {
    std::string qSmi = "Cl[C@@](*)(*)Br";
    std::string mSmi = "Cl[C@](C)(F)Br";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(matched);
  }
  {
    std::string qSmi = "[C@](C)(F)Br";
    std::string mSmi = "Cl[C@](C)(F)Br";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(!matched);
  }
  {
    std::string qSmi = "[C@@](C)(F)Br";
    std::string mSmi = "Cl[C@](C)(F)Br";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(matched);
  }
  {
    std::string qSmi = "C[C@](F)Br";
    std::string mSmi = "Cl[C@](C)(F)Br";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(!matched);
  }
  {
    std::string qSmi = "C[C@@](F)Br";
    std::string mSmi = "Cl[C@](C)(F)Br";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(matched);
  }
  {
    std::string qSmi = "Cl[C@](F)Br";
    std::string mSmi = "Cl[C@](C)(F)Br";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(matched);
  }
  {
    std::string qSmi = "Cl[C@](C)F";
    std::string mSmi = "Cl[C@](C)(F)Br";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(matched);
  }
}

TEST_CASE("testCisTransMatch", "[substruct][stereochemistry]") {
  {
    std::string qSmi = "CC=CC";
    std::string mSmi = "CC=CC";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(matched);
  }
  {
    std::string qSmi = "CC=CC";
    std::string mSmi = "C/C=C/C";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(matched);
  }
  {
    std::string qSmi = "C/C=C/C";
    std::string mSmi = "CC=CC";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(!matched);
  }
  {
    std::string qSmi = "C/C=C/C";
    std::string mSmi = "C/C=C\\C";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(!matched);
  }
  {
    std::string qSmi = "C/C=C/C";
    std::string mSmi = "C/C=C/C";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(matched);
  }
  {
    std::string qSmi = "C/C=C/C";
    std::string mSmi = "C/C=C(/F)C";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(!matched);
  }
  {
    std::string qSmi = "C/C=C/C";
    std::string mSmi = "C/C=C(\\F)C";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(matched);
  }
  {
    std::string qSmi = "C/C=C/C";
    std::string mSmi = "C/C(F)=C(\\F)C";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(matched);
  }
  {
    std::string qSmi = "C/C=C/C";
    std::string mSmi = "CC(/F)=C(\\F)C";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(matched);
  }
  {
    std::string qSmi = "C/C=C/C";
    std::string mSmi = "CC(\\F)=C(\\F)C";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true);
    CHECK(!matched);
  }
}

TEST_CASE("testCisTransMatch2", "[substruct][stereochemistry]") {
  {
    std::string qSmi = "CC=CC";
    std::string mSmi = "CCC=CC";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    CHECK(SubstructMatch(*mol, *query, matchV, true, true));

    CHECK(query->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
    query->getBondWithIdx(1)->getStereoAtoms().push_back(0);
    query->getBondWithIdx(1)->getStereoAtoms().push_back(3);
    query->getBondWithIdx(1)->setStereo(Bond::STEREOCIS);
    CHECK(!SubstructMatch(*mol, *query, matchV, true, true));
    CHECK(SubstructMatch(*mol, *query, matchV, true, false));

    CHECK(mol->getBondWithIdx(2)->getBondType() == Bond::DOUBLE);
    mol->getBondWithIdx(2)->getStereoAtoms().push_back(1);
    mol->getBondWithIdx(2)->getStereoAtoms().push_back(4);
    mol->getBondWithIdx(2)->setStereo(Bond::STEREOCIS);
    CHECK(SubstructMatch(*mol, *query, matchV, true, true));

    mol->getBondWithIdx(2)->setStereo(Bond::STEREOTRANS);
    CHECK(!SubstructMatch(*mol, *query, matchV, true, true));
    CHECK(SubstructMatch(*mol, *query, matchV, true, false));

    query->getBondWithIdx(1)->setStereo(Bond::STEREONONE);
    CHECK(SubstructMatch(*mol, *query, matchV, true, true));
    query->getBondWithIdx(1)->setStereo(Bond::STEREOANY);
    CHECK(SubstructMatch(*mol, *query, matchV, true, true));
  }
  {
    std::string qSmi = "CC=CC";
    std::string mSmi = "CC=C(C)F";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    CHECK(SubstructMatch(*mol, *query, matchV, true, true));

    CHECK(query->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
    query->getBondWithIdx(1)->getStereoAtoms().push_back(0);
    query->getBondWithIdx(1)->getStereoAtoms().push_back(3);
    query->getBondWithIdx(1)->setStereo(Bond::STEREOCIS);
    CHECK(!SubstructMatch(*mol, *query, matchV, true, true));
    CHECK(SubstructMatch(*mol, *query, matchV, true, false));

    CHECK(mol->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
    mol->getBondWithIdx(1)->getStereoAtoms().push_back(0);
    mol->getBondWithIdx(1)->getStereoAtoms().push_back(3);
    mol->getBondWithIdx(1)->setStereo(Bond::STEREOCIS);
    CHECK(SubstructMatch(*mol, *query, matchV, true, true));

    mol->getBondWithIdx(1)->setStereo(Bond::STEREOTRANS);
    CHECK(!SubstructMatch(*mol, *query, matchV, true, true));
    CHECK(SubstructMatch(*mol, *query, matchV, true, false));

    query->getBondWithIdx(1)->setStereo(Bond::STEREONONE);
    CHECK(SubstructMatch(*mol, *query, matchV, true, true));
    query->getBondWithIdx(1)->setStereo(Bond::STEREOANY);
    CHECK(SubstructMatch(*mol, *query, matchV, true, true));
  }
  {
    std::string qSmi = "CC=CC";
    std::string mSmi = "CCC=C(C)F";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    CHECK(SubstructMatch(*mol, *query, matchV, true, true));
    CHECK(query->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
    query->getBondWithIdx(1)->getStereoAtoms().push_back(0);
    query->getBondWithIdx(1)->getStereoAtoms().push_back(3);
    query->getBondWithIdx(1)->setStereo(Bond::STEREOCIS);
    CHECK(!SubstructMatch(*mol, *query, matchV, true, true));
    CHECK(SubstructMatch(*mol, *query, matchV, true, false));

    CHECK(mol->getBondWithIdx(2)->getBondType() == Bond::DOUBLE);
    mol->getBondWithIdx(2)->getStereoAtoms().push_back(1);
    mol->getBondWithIdx(2)->getStereoAtoms().push_back(4);
    mol->getBondWithIdx(2)->setStereo(Bond::STEREOCIS);
    CHECK(SubstructMatch(*mol, *query, matchV, true, true));

    mol->getBondWithIdx(2)->setStereo(Bond::STEREOTRANS);
    CHECK(!SubstructMatch(*mol, *query, matchV, true, true));
    CHECK(SubstructMatch(*mol, *query, matchV, true, false));

    query->getBondWithIdx(1)->setStereo(Bond::STEREONONE);
    CHECK(SubstructMatch(*mol, *query, matchV, true, true));
    query->getBondWithIdx(1)->setStereo(Bond::STEREOANY);
    CHECK(SubstructMatch(*mol, *query, matchV, true, true));
  }
  {  // now make it harder: the stereoatoms don't match, but the stereochemistry
     // does
    std::string qSmi = "CC=CC";
    std::string mSmi = "CCC=C(C)F";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    CHECK(SubstructMatch(*mol, *query, matchV, true, true));

    CHECK(query->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
    query->getBondWithIdx(1)->getStereoAtoms().push_back(0);
    query->getBondWithIdx(1)->getStereoAtoms().push_back(3);
    query->getBondWithIdx(1)->setStereo(Bond::STEREOCIS);
    CHECK(!SubstructMatch(*mol, *query, matchV, true, true));
    CHECK(SubstructMatch(*mol, *query, matchV, true, false));

    CHECK(mol->getBondWithIdx(2)->getBondType() == Bond::DOUBLE);
    mol->getBondWithIdx(2)->getStereoAtoms().push_back(1);
    mol->getBondWithIdx(2)->getStereoAtoms().push_back(5);
    mol->getBondWithIdx(2)->setStereo(Bond::STEREOTRANS);
    CHECK(SubstructMatch(*mol, *query, matchV, true, true));

    mol->getBondWithIdx(2)->setStereo(Bond::STEREOCIS);
    CHECK(!SubstructMatch(*mol, *query, matchV, true, true));
    CHECK(SubstructMatch(*mol, *query, matchV, true, false));

    query->getBondWithIdx(1)->setStereo(Bond::STEREONONE);
    CHECK(SubstructMatch(*mol, *query, matchV, true, true));
    query->getBondWithIdx(1)->setStereo(Bond::STEREOANY);
    CHECK(SubstructMatch(*mol, *query, matchV, true, true));
  }
  {  // now make it harder: the stereoatoms don't match on either end, but the
     // stereochemistry does
    std::string qSmi = "CC=CC";
    std::string mSmi = "CCC(F)=C(C)F";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    CHECK(SubstructMatch(*mol, *query, matchV, true, true));

    CHECK(query->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
    query->getBondWithIdx(1)->getStereoAtoms().push_back(0);
    query->getBondWithIdx(1)->getStereoAtoms().push_back(3);
    query->getBondWithIdx(1)->setStereo(Bond::STEREOCIS);
    CHECK(!SubstructMatch(*mol, *query, matchV, true, true));
    CHECK(SubstructMatch(*mol, *query, matchV, true, false));

    CHECK(mol->getBondWithIdx(3)->getBondType() == Bond::DOUBLE);
    mol->getBondWithIdx(3)->getStereoAtoms().push_back(3);
    mol->getBondWithIdx(3)->getStereoAtoms().push_back(6);
    mol->getBondWithIdx(3)->setStereo(Bond::STEREOCIS);
    CHECK(SubstructMatch(*mol, *query, matchV, true, true));

    mol->getBondWithIdx(3)->setStereo(Bond::STEREOTRANS);
    CHECK(!SubstructMatch(*mol, *query, matchV, true, true));
    CHECK(SubstructMatch(*mol, *query, matchV, true, false));

    query->getBondWithIdx(1)->setStereo(Bond::STEREONONE);
    CHECK(SubstructMatch(*mol, *query, matchV, true, true));
    query->getBondWithIdx(1)->setStereo(Bond::STEREOANY);
    CHECK(SubstructMatch(*mol, *query, matchV, true, true));
  }
}

TEST_CASE("testGitHubIssue15", "[substruct][github][issue15]") {
  {
    std::string qSmi = "[R2]~[R1]~[R2]";
    std::string mSmi = "CCC";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true, true);
    CHECK(!matched);
  }
  {
    std::string qSmi = "[R2]~[R1]~[R2]";
    std::string mSmi = "CCC";
    auto query = v2::SmilesParse::MolFromSmarts(qSmi);
    auto mol = v2::SmilesParse::MolFromSmiles(mSmi);
    MatchVectType matchV;
    bool matched = SubstructMatch(*mol, *query, matchV, true, true, true);
    CHECK(!matched);
  }
}

TEST_CASE("testGitHubIssue409", "[substruct][github][issue409]") {
  {
    std::string smi = "FC(F)(F)CC(F)(F)F";
    auto mol = v2::SmilesParse::MolFromSmiles(smi);
    std::vector<MatchVectType> matches;
    unsigned int matched =
        SubstructMatch(*mol, *mol, matches, false, true, false, false);
    CHECK(matched == matches.size());
    CHECK(matches.size() == 72);
    matched =
        SubstructMatch(*mol, *mol, matches, false, true, false, false, 16);
    CHECK(matches.size() == 16);
  }
}

TEST_CASE("testGitHubIssue688", "[substruct][github][issue688]") {
  {
    std::string smi = "C1CC[C@](Cl)(N)O1";
    auto mol = v2::SmilesParse::MolFromSmiles(smi);
    CHECK(mol);
    // mol->debugMol(std::cerr);
    std::string sma = "C1CC[C@](N)O1";
    auto qmol = v2::SmilesParse::MolFromSmarts(sma);
    CHECK(qmol);
    // qmol->updatePropertyCache();
    // qmol->debugMol(std::cerr);

    MatchVectType match;
    CHECK(SubstructMatch(*mol, *qmol, match, true, false));
    CHECK(match.size() == qmol->getNumAtoms());

    CHECK(SubstructMatch(*mol, *qmol, match, true, true));
    CHECK(match.size() == qmol->getNumAtoms());
  }
  {
    std::string smi = "C1CC[C@](Cl)(N)O1";
    auto mol = v2::SmilesParse::MolFromSmiles(smi);
    CHECK(mol);
    // mol->debugMol(std::cerr);
    std::string sma = "C1CC[C@@](N)O1";
    auto qmol = v2::SmilesParse::MolFromSmarts(sma);
    CHECK(qmol);
    // qmol->updatePropertyCache();
    // qmol->debugMol(std::cerr);

    MatchVectType match;
    CHECK(SubstructMatch(*mol, *qmol, match, true, false));
    CHECK(match.size() == qmol->getNumAtoms());

    CHECK(!SubstructMatch(*mol, *qmol, match, true, true));
  }
  {
    std::string smi = "N[C@]1(Cl)CCCO1";
    auto mol = v2::SmilesParse::MolFromSmiles(smi);
    CHECK(mol);
    // mol->debugMol(std::cerr);
    std::string sma = "N[C@]1CCCO1";
    auto qmol = v2::SmilesParse::MolFromSmarts(sma);
    CHECK(qmol);
    // qmol->updatePropertyCache();
    // qmol->debugMol(std::cerr);
    // std::cerr << MolToSmiles(*qmol, true) << std::endl;

    MatchVectType match;
    CHECK(SubstructMatch(*mol, *qmol, match, true, false));
    CHECK(match.size() == qmol->getNumAtoms());

    CHECK(SubstructMatch(*mol, *qmol, match, true, true));
    CHECK(match.size() == qmol->getNumAtoms());
  }
  {
    std::string smi = "N[C@]1(Cl)CCCO1";
    auto mol = v2::SmilesParse::MolFromSmiles(smi);
    CHECK(mol);
    // mol->debugMol(std::cerr);
    std::string sma = "N[C@@]1CCCO1";
    auto qmol = v2::SmilesParse::MolFromSmarts(sma);
    CHECK(qmol);
    // qmol->updatePropertyCache();
    // qmol->debugMol(std::cerr);
    // std::cerr << MolToSmiles(*qmol, true) << std::endl;

    MatchVectType match;
    CHECK(SubstructMatch(*mol, *qmol, match, true, false));
    CHECK(match.size() == qmol->getNumAtoms());

    CHECK(!SubstructMatch(*mol, *qmol, match, true, true));
  }
}

TEST_CASE("testDativeMatch", "[substruct][dative]") {
  {
    std::string smi = "[Cu]->[Fe]";
    auto mol = v2::SmilesParse::MolFromSmiles(smi);
    CHECK(mol);

    // make sure a self-match works
    MatchVectType match;
    CHECK(SubstructMatch(*mol, *mol, match));
    CHECK(match.size() == mol->getNumAtoms());

    {  // reverse the order and make sure that works
      std::string sma = "[Fe]<-[Cu]";
      auto qmol = v2::SmilesParse::MolFromSmarts(sma);
      CHECK(qmol);
      MatchVectType match;
      CHECK(SubstructMatch(*mol, *qmol, match));
      CHECK(match.size() == qmol->getNumAtoms());
    }
    {  // reverse the direction and make sure that does not work.
      std::string sma = "[Fe]->[Cu]";
      auto qmol = v2::SmilesParse::MolFromSmarts(sma);
      CHECK(qmol);
      MatchVectType match;
      CHECK(!SubstructMatch(*mol, *qmol, match));
    }
  }
}

TEST_CASE("testGithubIssue1489", "[substruct][github][issue1489]") {
  {
    std::string smi1 = "CCC[C@@H]1CN(CCC)CCN1";
    std::string smi2 = "CCC[C@H]1CN(CCC)CCN1";
    auto mol1 = v2::SmilesParse::MolFromSmiles(smi1);
    CHECK(mol1);
    auto mol2 = v2::SmilesParse::MolFromSmiles(smi2);
    CHECK(mol2);

    bool recursionPossible = true;
    bool useChirality = true;

    MatchVectType match;
    // make sure self-matches work
    CHECK(SubstructMatch(*mol1, *mol1, match, recursionPossible, useChirality));
    CHECK(SubstructMatch(*mol2, *mol2, match, recursionPossible, useChirality));

    // check matches using the molecules from smiles:
    CHECK(
        !SubstructMatch(*mol1, *mol2, match, recursionPossible, useChirality));
    CHECK(SubstructMatch(*mol1, *mol2, match, recursionPossible, false));

    {
      auto qmol1 = v2::SmilesParse::MolFromSmarts(smi1);

      qmol1->updatePropertyCache();
      CHECK(qmol1);
      CHECK(SubstructMatch(*mol1, *qmol1, match, recursionPossible, false));
      CHECK(SubstructMatch(*mol1, *qmol1, match, recursionPossible,
                           useChirality));
      CHECK(SubstructMatch(*mol2, *qmol1, match, recursionPossible, false));
      CHECK(!SubstructMatch(*mol2, *qmol1, match, recursionPossible,
                            useChirality));
    }
  }
  {
    std::string smi1 = "F([C@@H](Cl)Br)";
    auto mol1 = v2::SmilesParse::MolFromSmiles(smi1);
    CHECK(mol1);

    bool recursionPossible = true;
    bool useChirality = true;

    MatchVectType match;
    // make sure self-matches work
    CHECK(SubstructMatch(*mol1, *mol1, match, recursionPossible, useChirality));
    {
      auto qmol1 = v2::SmilesParse::MolFromSmarts(smi1);

      qmol1->updatePropertyCache();
      CHECK(qmol1);
      CHECK(SubstructMatch(*mol1, *qmol1, match, recursionPossible, false));
      CHECK(SubstructMatch(*mol1, *qmol1, match, recursionPossible,
                           useChirality));
    }
    {
      std::string smi2 = "F([C@H](Br)Cl)";
      auto qmol1 = v2::SmilesParse::MolFromSmarts(smi2);
      qmol1->updatePropertyCache();
      CHECK(qmol1);
      CHECK(SubstructMatch(*mol1, *qmol1, match, recursionPossible, false));
      CHECK(SubstructMatch(*mol1, *qmol1, match, recursionPossible,
                           useChirality));
    }
  }
}

TEST_CASE("testGithub2570", "[substruct][github][issue2570]") {
  bool uniquify = true;
  bool recursionPossible = true;
  bool useChirality = true;
  {
    const auto mol = R"(C[C@](Cl)(Br)F)"_smiles;

    {
      const auto query = R"([C@](Cl)(Br)F)"_smarts;
      std::vector<MatchVectType> matches;
      CHECK(!SubstructMatch(*mol, *query, matches, uniquify, recursionPossible,
                            useChirality));
    }
    {
      const auto query = R"([C@@](Cl)(Br)F)"_smarts;
      std::vector<MatchVectType> matches;
      CHECK(SubstructMatch(*mol, *query, matches, uniquify, recursionPossible,
                           useChirality));
    }
    {  // Swap order of a pair of atoms
      const auto query = R"([C@](Br)(Cl)F)"_smarts;
      std::vector<MatchVectType> matches;
      CHECK(SubstructMatch(*mol, *query, matches, uniquify, recursionPossible,
                           useChirality));
    }
    {
      const auto query = R"([C@@](Br)(Cl)F)"_smarts;
      std::vector<MatchVectType> matches;
      CHECK(!SubstructMatch(*mol, *query, matches, uniquify, recursionPossible,
                            useChirality));
    }
    {  // Smaller fragments should always match as long as they have have a
       // chiral tag,
      // as these don't have enough neighbors to define CW/CCW chirality
      const std::vector<std::string> smarts({"[C@](Cl)Br", "[C@@](Cl)Br",
                                             "[C@](Br)F", "[C@@](Br)F", "[C@]F",
                                             "[C@@]F", "[C@]", "[C@@]"});
      std::vector<MatchVectType> matches;
      for (const auto &sma : smarts) {
        std::unique_ptr<ROMol> query(SmartsToMol(sma));
        CHECK(SubstructMatch(*mol, *query, matches, uniquify, recursionPossible,
                             useChirality));
      }
    }
  }
  {  // Mol also starting with the chiral atom
    const auto mol = R"([C@](C)(Cl)(Br)F)"_smiles;
    {
      const auto query = R"([C@](C)(Cl)(Br)F)"_smarts;
      std::vector<MatchVectType> matches;
      CHECK(SubstructMatch(*mol, *query, matches, uniquify, recursionPossible,
                           useChirality));
    }
    {
      const auto query = R"([C@@](C)(Cl)(Br)F)"_smarts;
      std::vector<MatchVectType> matches;
      CHECK(!SubstructMatch(*mol, *query, matches, uniquify, recursionPossible,
                            useChirality));
    }
    {
      const auto query = R"([C@](C)(Cl)Br)"_smarts;
      std::vector<MatchVectType> matches;
      CHECK(SubstructMatch(*mol, *query, matches, uniquify, recursionPossible,
                           useChirality));
    }
    {
      const auto query = R"([C@@](C)(Cl)Br)"_smarts;
      std::vector<MatchVectType> matches;
      CHECK(!SubstructMatch(*mol, *query, matches, uniquify, recursionPossible,
                            useChirality));
    }
    {
      const auto query = R"([C@](Cl)(Br)F)"_smarts;
      std::vector<MatchVectType> matches;
      CHECK(!SubstructMatch(*mol, *query, matches, uniquify, recursionPossible,
                            useChirality));
    }
    {
      const auto query = R"([C@@](Cl)(Br)F)"_smarts;
      std::vector<MatchVectType> matches;
      CHECK(SubstructMatch(*mol, *query, matches, uniquify, recursionPossible,
                           useChirality));
    }
    {  // Swap order of a pair of atoms
      const auto query = R"([C@](Br)(Cl)F)"_smarts;
      std::vector<MatchVectType> matches;
      CHECK(SubstructMatch(*mol, *query, matches, uniquify, recursionPossible,
                           useChirality));
    }
    {
      const auto query = R"([C@@](Br)(Cl)F)"_smarts;
      std::vector<MatchVectType> matches;
      CHECK(!SubstructMatch(*mol, *query, matches, uniquify, recursionPossible,
                            useChirality));
    }
  }

  {  // Start from a physical H atom
    const auto mol = R"([H][C@](O)(F)Cl)"_smiles;
    const auto smarts = MolToSmarts(*mol);
    std::unique_ptr<ROMol> query(SmartsToMol(smarts));
    CHECK(smarts == R"([#8]-[#6@@H](-[#9])-[#17])");
    std::vector<MatchVectType> matches;
    CHECK(SubstructMatch(*mol, *query, matches, uniquify, recursionPossible,
                         useChirality));
  }
  {
    const auto mol = R"([H][C@](O)(F)Cl)"_smiles;
    const auto query = R"([C@H](O)(F)Cl)"_smarts;
    std::vector<MatchVectType> matches;
    CHECK(SubstructMatch(*mol, *query, matches, uniquify, recursionPossible,
                         useChirality));
  }
  {
    const auto mol = R"([H][C@](O)(F)Cl)"_smiles;
    const auto query = R"([C@@H](O)(F)Cl)"_smarts;
    std::vector<MatchVectType> matches;
    CHECK(!SubstructMatch(*mol, *query, matches, uniquify, recursionPossible,
                          useChirality));
  }
  {  // Start from an attached H atom
    const auto mol = R"([C@H](O)(F)Cl)"_smiles;
    const auto smarts = MolToSmarts(*mol);
    CHECK(smarts == R"([#8]-[#6@@H](-[#9])-[#17])");
    const auto query = v2::SmilesParse::MolFromSmarts(smarts);
    std::vector<MatchVectType> matches;
    CHECK(SubstructMatch(*mol, *query, matches, uniquify, recursionPossible,
                         useChirality));
  }
  {
    const auto mol = R"([C@H](O)(F)Cl)"_smiles;
    const auto query = R"([C@H](O)(F)Cl)"_smarts;
    std::vector<MatchVectType> matches;
    CHECK(SubstructMatch(*mol, *query, matches, uniquify, recursionPossible,
                         useChirality));
  }
  {
    const auto mol = R"([C@H](O)(F)Cl)"_smiles;
    const auto query = R"([C@@H](O)(F)Cl)"_smarts;
    std::vector<MatchVectType> matches;
    CHECK(!SubstructMatch(*mol, *query, matches, uniquify, recursionPossible,
                          useChirality));
  }
  {  // Without H
    const auto mol = R"([C@](O)(F)(Cl)C)"_smiles;
    const auto smarts = MolToSmarts(*mol);
    const auto query = v2::SmilesParse::MolFromSmarts(smarts);
    CHECK(smarts == R"([#8]-[#6@](-[#9])(-[#17])-[#6])");
    std::vector<MatchVectType> matches;
    CHECK(SubstructMatch(*mol, *query, matches, uniquify, recursionPossible,
                         useChirality));
  }
  {
    const auto mol = R"([C@](O)(F)(Cl)C)"_smiles;
    const auto query = R"([#6@](-[#8])(-[#9])(-[#17])-[#6])"_smarts;
    std::vector<MatchVectType> matches;
    CHECK(SubstructMatch(*mol, *query, matches, uniquify, recursionPossible,
                         useChirality));
  }

  {  // What about queries not coming from SMARTS?
    const std::vector<std::string> smiles(  // These are all equivalent
        {"N[C@@]([H])(C)C(=O)O", "N[C@@H](C)C(=O)O", "N[C@H](C(=O)O)C",
         "[H][C@](N)(C)C(=O)O", "[C@H](N)(C)C(=O)O"});
    for (const auto &smi1 : smiles) {
      const auto mol1 = std::unique_ptr<ROMol>(SmilesToMol(smi1));
      for (const auto &smi2 : smiles) {  // Test them in both directions
        const auto mol2 = std::unique_ptr<ROMol>(SmilesToMol(smi2));
        std::vector<MatchVectType> matches;
        CHECK(SubstructMatch(*mol1, *mol2, matches, uniquify, recursionPossible,
                             useChirality));
      };
    }
  }
}

TEST_CASE("testEZVsCisTransMatch", "[substruct][stereochemistry]") {
  UseLegacyStereoPerceptionFixture fx(true);
  const auto mol = R"(F/C(C)=C(C)/Cl)"_smiles;
  {
    const Bond *stereoBnd = mol->getBondWithIdx(2);
    CHECK(stereoBnd->getStereo() == Bond::STEREOE);
  }

  // pairs of {query, matching expectation}
  const std::vector<std::pair<std::string, bool>> checks({
      {R"(F/C(C)=C(C)/Cl)", true},   // identical
      {R"(F\C(C)=C(C)\Cl)", true},   // symmetric
      {R"(F/C(C)=C(C)\Cl)", false},  // opposite
      {R"(F\C(C)=C(C)/Cl)", false}   // symmetric opposite
  });

  // Test with same stereoatoms as mol
  for (const auto &check : checks) {
    auto query = v2::SmilesParse::MolFromSmiles(check.first);
    {
      Bond *stereoBnd = query->getBondWithIdx(2);
      auto stereo = stereoBnd->getStereo();
      CHECK((stereo == Bond::STEREOE || stereo == Bond::STEREOZ));

      stereoBnd->setStereoAtoms(0, 5);  // Same as mol
      stereo = Chirality::translateEZLabelToCisTrans(stereo);
      CHECK((stereo == Bond::STEREOCIS || stereo == Bond::STEREOTRANS));
      stereoBnd->setStereo(stereo);
    }
    MatchVectType match;
    bool recursionPossible = true;
    bool useChirality = true;
    CHECK(check.second ==
          SubstructMatch(*mol, *query, match, recursionPossible, useChirality));
  }
  // Symmetrize stereoatoms
  for (const auto &check : checks) {
    auto query = v2::SmilesParse::MolFromSmiles(check.first);
    {
      Bond *stereoBnd = query->getBondWithIdx(2);
      auto stereo = stereoBnd->getStereo();
      CHECK((stereo == Bond::STEREOE || stereo == Bond::STEREOZ));

      stereoBnd->setStereoAtoms(2, 4);  // symmetric to mol
      stereo = Chirality::translateEZLabelToCisTrans(stereo);
      CHECK((stereo == Bond::STEREOCIS || stereo == Bond::STEREOTRANS));
      stereoBnd->setStereo(stereo);
    }
    MatchVectType match;
    bool recursionPossible = true;
    bool useChirality = true;
    CHECK(check.second ==
          SubstructMatch(*mol, *query, match, recursionPossible, useChirality));
  }
  // Flip one stereoatom and the label
  for (const auto &check : checks) {
    auto query = v2::SmilesParse::MolFromSmiles(check.first);
    {
      Bond *stereoBnd = query->getBondWithIdx(2);
      auto stereo = stereoBnd->getStereo();
      CHECK((stereo == Bond::STEREOE || stereo == Bond::STEREOZ));

      stereoBnd->setStereoAtoms(0, 4);  // Reverse second stereoatom
      if (stereo == Bond::STEREOE) {
        stereo = Bond::STEREOCIS;
      } else {
        stereo = Bond::STEREOTRANS;
      }
      stereoBnd->setStereo(stereo);
    }
    MatchVectType match;
    bool recursionPossible = true;
    bool useChirality = true;
    CHECK(check.second ==
          SubstructMatch(*mol, *query, match, recursionPossible, useChirality));
  }
}

TEST_CASE("testMostSubstitutedCoreMatch", "[substruct][core]") {
  auto core = "[*:1]c1cc([*:2])ccc1[*:3]"_smarts;
  auto orthoMeta = "c1ccc(-c2ccc(-c3ccccc3)c(-c3ccccc3)c2)cc1"_smiles;
  auto ortho = "c1ccc(-c2ccccc2-c2ccccc2)cc1"_smiles;
  auto meta = "c1ccc(-c2cccc(-c3ccccc3)c2)cc1"_smiles;
  auto biphenyl = "c1ccccc1-c1ccccc1"_smiles;
  auto phenyl = "c1ccccc1"_smiles;

  struct numHsMatchingDummies {
    static unsigned int get(const ROMol &mol, const ROMol &core,
                            const MatchVectType &match) {
      return std::count_if(
          match.begin(), match.end(),
          [&mol, &core](const std::pair<int, int> &pair) {
            return (core.getAtomWithIdx(pair.first)->getAtomicNum() == 0 &&
                    mol.getAtomWithIdx(pair.second)->getAtomicNum() == 1);
          });
    }
  };

  const auto &coreRef = *core;
  for (auto &molResPair : std::vector<std::pair<RDKit::RWMol *, unsigned>>{
           std::make_pair(orthoMeta.get(), 0u), std::make_pair(ortho.get(), 1u),
           std::make_pair(meta.get(), 1u), std::make_pair(biphenyl.get(), 2u),
           std::make_pair(phenyl.get(), 3u)}) {
    auto &mol = *molResPair.first;
    const auto res = molResPair.second;
    MolOps::addHs(mol);
    auto matches = SubstructMatch(mol, coreRef);
    auto bestMatch = getMostSubstitutedCoreMatch(mol, coreRef, matches);
    CHECK(numHsMatchingDummies::get(mol, coreRef, bestMatch) == res);
    std::vector<unsigned int> ctrlCounts(matches.size());
    std::transform(matches.begin(), matches.end(), ctrlCounts.begin(),
                   [&mol, &coreRef](const MatchVectType &match) {
                     return numHsMatchingDummies::get(mol, coreRef, match);
                   });
    std::sort(ctrlCounts.begin(), ctrlCounts.end());
    std::vector<unsigned int> sortedCounts(matches.size());
    auto sortedMatches =
        sortMatchesByDegreeOfCoreSubstitution(mol, coreRef, matches);
    std::transform(sortedMatches.begin(), sortedMatches.end(),
                   sortedCounts.begin(),
                   [&mol, &coreRef](const MatchVectType &match) {
                     return numHsMatchingDummies::get(mol, coreRef, match);
                   });
    CHECK(ctrlCounts == sortedCounts);
  }
  std::vector<MatchVectType> emptyMatches;
  bool raised = false;
  try {
    getMostSubstitutedCoreMatch(*orthoMeta, coreRef, emptyMatches);
  } catch (const Invar::Invariant &) {
    raised = true;
  }
  CHECK(raised);
  raised = false;
  try {
    sortMatchesByDegreeOfCoreSubstitution(*orthoMeta, coreRef, emptyMatches);
  } catch (const Invar::Invariant &) {
    raised = true;
  }
  CHECK(raised);
}

TEST_CASE("testLongRing", "[substruct][ring]") {
  std::string mol_smiles = "c12ccc(CCCCCCCc5ccc(C2)cc5)cc1";
  std::string query_smiles = "c1cc2ccc1CCCCCCCc1ccc(cc1)C2";
  auto mol = v2::SmilesParse::MolFromSmiles(mol_smiles);
  auto query = v2::SmilesParse::MolFromSmiles(query_smiles);
  CHECK(MolToSmiles(*query) == MolToSmiles(*mol));
  MatchVectType match1;
  MatchVectType match2;
  SubstructMatchParameters params;
  CHECK(SubstructMatch(*mol, *query, match1));
  CHECK(SubstructMatch(*query, *mol, match2));
}

TEST_CASE("testIsAtomTerminalRGroupOrQueryHydrogen", "[substruct][rgroup]") {
  {
    auto mol = R"CTAB(
  MJ201100                      

  7  7  0  0  0  0  0  0  0  0999 V2000
   -0.3795    1.5839    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0939    1.1714    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0939    0.3463    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3795   -0.0661    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3349    0.3463    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3349    1.1714    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0494    1.5839    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  4  0  0  0  0
  2  3  4  0  0  0  0
  3  4  4  0  0  0  0
  4  5  4  0  0  0  0
  5  6  4  0  0  0  0
  6  1  4  0  0  0  0
  6  7  1  0  0  0  0
M  RGP  1   7   1
M  END
)CTAB"_ctab;
    const auto rAtom = mol->getAtomWithIdx(mol->getNumAtoms() - 1);
    CHECK(isAtomTerminalRGroupOrQueryHydrogen(rAtom));
  }
  {
    auto mol = R"CTAB(
  MJ201100                      

  6  6  0  0  0  0  0  0  0  0999 V2000
   -0.7589    1.4277    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4733    1.0152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4733    0.1901    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7589   -0.2223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0444    0.1901    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0444    1.0152    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
M  RGP  1   6   1
M  END
)CTAB"_ctab;
    const auto rAtom = mol->getAtomWithIdx(mol->getNumAtoms() - 1);
    CHECK(!isAtomTerminalRGroupOrQueryHydrogen(rAtom));
  }
  {
    auto mol = R"CTAB(
  MJ201100                      

  7  7  0  0  0  0  0  0  0  0999 V2000
   -0.9152    0.2893    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6296   -0.1231    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6296   -0.9482    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9152   -1.3607    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2007   -0.9482    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2007   -0.1231    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5137    0.2893    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  4  0  0  0  0
  2  3  4  0  0  0  0
  3  4  4  0  0  0  0
  4  5  4  0  0  0  0
  5  6  4  0  0  0  0
  6  1  4  0  0  0  0
  6  7  1  0  0  0  0
M  ALS   7 10 F H   C   N   O   F   P   S   Cl  Br  I   
M  END
)CTAB"_ctab;
    const auto rAtom = mol->getAtomWithIdx(mol->getNumAtoms() - 1);
    CHECK(isAtomTerminalRGroupOrQueryHydrogen(rAtom));
  }
  {
    auto mol = R"CTAB(
  MJ201100                      

  7  7  0  0  0  0  0  0  0  0999 V2000
   -0.9152    0.2893    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6296   -0.1231    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6296   -0.9482    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9152   -1.3607    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2007   -0.9482    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2007   -0.1231    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5137    0.2893    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  4  0  0  0  0
  2  3  4  0  0  0  0
  3  4  4  0  0  0  0
  4  5  4  0  0  0  0
  5  6  4  0  0  0  0
  6  1  4  0  0  0  0
  6  7  1  0  0  0  0
M  ALS   7  9 F C   N   O   F   P   S   Cl  Br  I   
M  END
)CTAB"_ctab;
    const auto rAtom = mol->getAtomWithIdx(mol->getNumAtoms() - 1);
    CHECK(!isAtomTerminalRGroupOrQueryHydrogen(rAtom));
  }
}