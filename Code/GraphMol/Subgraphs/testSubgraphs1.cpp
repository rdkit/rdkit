//
//  Copyright (C) 2003-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace std;
using namespace RDKit;

TEST_CASE("testSubgraphs", "[subgraphs]") {
  // build: CCC(C)CC
  RWMol mol;
  bool updateLabel = true;
  bool takeOwnership = true;
  mol.addAtom(new Atom(6), updateLabel, takeOwnership);
  mol.addAtom(new Atom(6), updateLabel, takeOwnership);
  mol.addAtom(new Atom(6), updateLabel, takeOwnership);
  mol.addAtom(new Atom(6), updateLabel, takeOwnership);
  mol.addAtom(new Atom(6), updateLabel, takeOwnership);
  mol.addAtom(new Atom(6), updateLabel, takeOwnership);
  mol.addBond(0, 1, Bond::SINGLE);
  mol.addBond(1, 2, Bond::SINGLE);
  mol.addBond(2, 3, Bond::SINGLE);
  mol.addBond(2, 4, Bond::SINGLE);
  mol.addBond(3, 5, Bond::SINGLE);

  PATH_LIST tmp;
  PATH_LIST::iterator i;

  int totPs = 0;
  tmp = findAllSubgraphsOfLengthN(mol, 1);
  CHECK(tmp.size() == 5);
  totPs += tmp.size();
  tmp = findAllSubgraphsOfLengthN(mol, 2);
  CHECK(tmp.size() == 5);
  totPs += tmp.size();
  tmp = findAllSubgraphsOfLengthN(mol, 3);
  CHECK(tmp.size() == 5);
  totPs += tmp.size();
  tmp = findAllSubgraphsOfLengthN(mol, 4);
  CHECK(tmp.size() == 3);
  totPs += tmp.size();
  tmp = findAllSubgraphsOfLengthN(mol, 5);
  CHECK(tmp.size() == 1);
  totPs += tmp.size();
  tmp = findAllSubgraphsOfLengthN(mol, 6);
  CHECK(tmp.empty());
  totPs += tmp.size();

  // now use the direct range function and check that we get the
  // same answer
  INT_PATH_LIST_MAP tmpm;
  tmpm = findAllSubgraphsOfLengthsMtoN(mol, 1, 6);
  int newTot, idx;
  newTot = 0;
  for (idx = 1; idx <= 6; idx++) {
    newTot += tmpm[idx].size();
  }
  CHECK(totPs == newTot);

  // add an H and make sure things don't change:
  mol.addAtom(new Atom(1), updateLabel, takeOwnership);
  mol.addBond(5, 6, Bond::SINGLE);

  tmp = findAllSubgraphsOfLengthN(mol, 1);

  CHECK(tmp.size() == 5);
  tmp = findAllSubgraphsOfLengthN(mol, 2);
  CHECK(tmp.size() == 5);
  tmp = findAllSubgraphsOfLengthN(mol, 3);
  CHECK(tmp.size() == 5);
  tmp = findAllSubgraphsOfLengthN(mol, 4);
  CHECK(tmp.size() == 3);
  tmp = findAllSubgraphsOfLengthN(mol, 5);
  CHECK(tmp.size() == 1);
  tmp = findAllSubgraphsOfLengthN(mol, 6);
  CHECK(tmp.empty());
}

TEST_CASE("testSubgraphs2a", "[subgraphs]") {
  // these have been moved here from test2.cpp
  auto mol = "C1CC2C1N2"_smiles;
  REQUIRE(mol);
  int nAll = 0;
  int nUnique = 0;
  for (unsigned int i = 1; i < 7; ++i) {
    PATH_LIST tmp;
    tmp = findAllSubgraphsOfLengthN(*mol, i);
    nAll += tmp.size();
    tmp = findUniqueSubgraphsOfLengthN(*mol, i);
    nUnique += tmp.size();
  }

  CHECK(nAll == 49);
  CHECK(nUnique == 21);
}

TEST_CASE("testSubgraphs2", "[subgraphs]") {
  // these have been moved here from test2.cpp
  auto mol = "C12C3C4C1C1C2C3N41"_smiles;
  REQUIRE(mol);

  int nAll = 0;
  int nUnique = 0;
  int i;
  for (i = 1; i < 13; i++) {
    PATH_LIST tmp;
    tmp = findAllSubgraphsOfLengthN(*mol, i);
    nAll += tmp.size();
    tmp = findUniqueSubgraphsOfLengthN(*mol, i);
    nUnique += tmp.size();
  }

  CHECK(nAll == 2433);
  CHECK(nUnique == 300);

  mol = "CCC(O)C(c1ccccc1)CC(C)N(C)C"_smiles;
  REQUIRE(mol);

  nAll = 0;
  nUnique = 0;
  for (i = 1; i < 18; i++) {
    PATH_LIST tmp;
    tmp = findAllSubgraphsOfLengthN(*mol, i);
    nAll += tmp.size();
    tmp = findUniqueSubgraphsOfLengthN(*mol, i);
    nUnique += tmp.size();
  }

  CHECK(nAll == 1990);
  CHECK(nUnique == 907);
}

void dumpVIV(PATH_LIST v) {
  PATH_LIST::iterator i;
  PATH_TYPE::iterator j;
  for (i = v.begin(); i != v.end(); i++) {
    for (j = i->begin(); j != i->end(); j++) {
      std::cout << *j << " ";
    }
    std::cout << std::endl;
  }
}

TEST_CASE("testPaths", "[subgraphs]") {
  // build: CCC(C)CC
  RWMol mol;
  bool updateLabel = true;
  bool takeOwnership = true;
  mol.addAtom(new Atom(6), updateLabel, takeOwnership);
  mol.addAtom(new Atom(6), updateLabel, takeOwnership);
  mol.addAtom(new Atom(6), updateLabel, takeOwnership);
  mol.addAtom(new Atom(6), updateLabel, takeOwnership);
  mol.addAtom(new Atom(6), updateLabel, takeOwnership);
  mol.addAtom(new Atom(6), updateLabel, takeOwnership);
  mol.addBond(0, 1, Bond::SINGLE);
  mol.addBond(1, 2, Bond::SINGLE);
  mol.addBond(2, 3, Bond::SINGLE);
  mol.addBond(2, 4, Bond::SINGLE);
  mol.addBond(3, 5, Bond::SINGLE);

  PATH_LIST tmp;

  //
  //  Retrieve using bonds
  //
  tmp = findAllPathsOfLengthN(mol, 1);
  CHECK(tmp.size() == 5);
  tmp = findAllPathsOfLengthN(mol, 2);
  CHECK(tmp.size() == 5);
  tmp = findAllPathsOfLengthN(mol, 3);
  CHECK(tmp.size() == 4);
  tmp = findAllPathsOfLengthN(mol, 4);
  CHECK(tmp.size() == 1);
  tmp = findAllPathsOfLengthN(mol, 5);
  CHECK(tmp.empty());
  tmp = findAllPathsOfLengthN(mol, 6);
  CHECK(tmp.empty());

  //
  //  Retrieve using atoms, which gives the results shifted by
  //  one (it takes two atoms to make one bond)
  //
  tmp = findAllPathsOfLengthN(mol, 1, false);
  CHECK(tmp.size() == 6);
  tmp = findAllPathsOfLengthN(mol, 2, false);
  CHECK(tmp.size() == 5);
  tmp = findAllPathsOfLengthN(mol, 3, false);
  CHECK(tmp.size() == 5);
  tmp = findAllPathsOfLengthN(mol, 4, false);
  CHECK(tmp.size() == 4);
  tmp = findAllPathsOfLengthN(mol, 5, false);
  CHECK(tmp.size() == 1);
  tmp = findAllPathsOfLengthN(mol, 6, false);
  CHECK(tmp.empty());

  //
  //  try m->n
  //
  INT_PATH_LIST_MAP pths;
  pths = findAllPathsOfLengthsMtoN(mol, 1, 6);
  CHECK(pths[1].size() == 5);
  CHECK(pths[2].size() == 5);
  CHECK(pths[3].size() == 4);
  CHECK(pths[4].size() == 1);
  CHECK(pths[5].empty());
  CHECK(pths[6].empty());

  pths = findAllPathsOfLengthsMtoN(mol, 1, 6, false);
  CHECK(pths[1].size() == 6);
  CHECK(pths[2].size() == 5);
  CHECK(pths[3].size() == 5);
  CHECK(pths[4].size() == 4);
  CHECK(pths[5].size() == 1);
  CHECK(pths[6].empty());

  //
  //  add an atom, close the ring and re-check a couple indices:
  //   (leaves us with CC1CCCCC1)
  //
  mol.addAtom(new Atom(6), updateLabel, takeOwnership);
  mol.addBond(5, 6, Bond::SINGLE);
  mol.addBond(0, 6, Bond::SINGLE);
  tmp = findAllPathsOfLengthN(mol, 4);
  CHECK(tmp.size() == 8);
  tmp = findAllPathsOfLengthN(mol, 5, false);
  CHECK(tmp.size() == 8);
}

TEST_CASE("testPaths2", "[subgraphs]") {
  // build: CCC(C)CC
  RWMol mol;
  bool updateLabel = true;
  bool takeOwnership = true;
  mol.addAtom(new Atom(6), updateLabel, takeOwnership);
  mol.addAtom(new Atom(6), updateLabel, takeOwnership);
  mol.addAtom(new Atom(6), updateLabel, takeOwnership);
  mol.addAtom(new Atom(6), updateLabel, takeOwnership);
  mol.addBond(0, 1, Bond::SINGLE);
  mol.addBond(1, 2, Bond::SINGLE);
  mol.addBond(2, 3, Bond::SINGLE);
  mol.addBond(2, 0, Bond::SINGLE);

  //
  //  Retrieve using bonds
  //
  PATH_LIST tmp = findAllPathsOfLengthN(mol, 3);
  // std::cout << "\n3:" << std::endl;
  // dumpVIV(tmp);
  CHECK(tmp.size() == 3);
}

TEST_CASE("testUniqueSubgraphs", "[subgraphs]") {
  // build: CCC(C)CC
  RWMol mol;
  bool updateLabel = true;
  bool takeOwnership = true;
  mol.addAtom(new Atom(6), updateLabel, takeOwnership);
  mol.addAtom(new Atom(6), updateLabel, takeOwnership);
  mol.addAtom(new Atom(6), updateLabel, takeOwnership);
  mol.addAtom(new Atom(6), updateLabel, takeOwnership);
  mol.addAtom(new Atom(6), updateLabel, takeOwnership);
  mol.addAtom(new Atom(6), updateLabel, takeOwnership);
  mol.addBond(0, 1, Bond::SINGLE);
  mol.addBond(1, 2, Bond::SINGLE);
  mol.addBond(2, 3, Bond::SINGLE);
  mol.addBond(2, 4, Bond::SINGLE);
  mol.addBond(3, 5, Bond::SINGLE);

  PATH_LIST tmp;
  PATH_LIST::iterator i;

  tmp = findAllSubgraphsOfLengthN(mol, 1);
  CHECK(tmp.size() == 5);

  tmp = findAllSubgraphsOfLengthN(mol, 2);
  CHECK(tmp.size() == 5);
  tmp = findUniqueSubgraphsOfLengthN(mol, 2);
  CHECK(tmp.size() == 1);

  tmp = findAllSubgraphsOfLengthN(mol, 3);
  CHECK(tmp.size() == 5);
  tmp = findUniqueSubgraphsOfLengthN(mol, 3);
  CHECK(tmp.size() == 2);

  tmp = findAllSubgraphsOfLengthN(mol, 4);
  CHECK(tmp.size() == 3);
  tmp = findUniqueSubgraphsOfLengthN(mol, 4);
  CHECK(tmp.size() == 2);

  tmp = findAllSubgraphsOfLengthN(mol, 5);
  CHECK(tmp.size() == 1);
  tmp = findUniqueSubgraphsOfLengthN(mol, 5);
  CHECK(tmp.size() == 1);

  tmp = findAllSubgraphsOfLengthN(mol, 6);
  CHECK(tmp.empty());
  tmp = findUniqueSubgraphsOfLengthN(mol, 6);
  CHECK(tmp.empty());

  // add an H and make sure things don't change:
  mol.addAtom(new Atom(1), updateLabel, takeOwnership);
  mol.addBond(5, 6, Bond::SINGLE);

  tmp = findAllSubgraphsOfLengthN(mol, 2);
  CHECK(tmp.size() == 5);
  tmp = findUniqueSubgraphsOfLengthN(mol, 2);
  CHECK(tmp.size() == 1);

  tmp = findAllSubgraphsOfLengthN(mol, 3);
  CHECK(tmp.size() == 5);
  tmp = findUniqueSubgraphsOfLengthN(mol, 3);
  CHECK(tmp.size() == 2);

  tmp = findAllSubgraphsOfLengthN(mol, 4);
  CHECK(tmp.size() == 3);
  tmp = findUniqueSubgraphsOfLengthN(mol, 4);
  CHECK(tmp.size() == 2);

  tmp = findAllSubgraphsOfLengthN(mol, 5);
  CHECK(tmp.size() == 1);
  tmp = findUniqueSubgraphsOfLengthN(mol, 5);
  CHECK(tmp.size() == 1);

  tmp = findAllSubgraphsOfLengthN(mol, 6);
  CHECK(tmp.empty());
  tmp = findUniqueSubgraphsOfLengthN(mol, 6);
  CHECK(tmp.empty());
}

TEST_CASE("testUniqueSubgraphs2", "[subgraphs]") {
  // moved here from test2.cpp
  auto mol = "O=C(O)CCCC=CC(C1C(O)CC(O)C1(C=CC(O)CCCCC))"_smiles;
  REQUIRE(mol);

  int nAll = 0;
  int nUnique = 0;
  for (int i = 1; i < 26; i++) {
    PATH_LIST tmp;
    tmp = findAllSubgraphsOfLengthN(*mol, i);
    // std::cout << i << "\t" << tmp.size();
    nAll += tmp.size();
    tmp = findUniqueSubgraphsOfLengthN(*mol, i);
    // std::cout << "\t" << tmp.size() << std::endl;;
    nUnique += tmp.size();
  }

  CHECK(nAll == 6435);
  CHECK(nUnique == 5618);
}

TEST_CASE("testLeak", "[subgraphs]") {
  // moved here from test2.cpp
  // testing for a core leak (Issue 42)
  auto mol = "O=C(O)CCCC=CC(C1C(O)CC(O)C1(C=CC(O)CCCCC))"_smiles;
  REQUIRE(mol);

  for (int rep = 0; rep < 100; rep++) {
    int nAll = 0;
    int nUnique = 0;
    for (int i = 1; i < 26; i++) {
      PATH_LIST tmp;
      tmp = findAllSubgraphsOfLengthN(*mol, i);
      // std::cout << i << "\t" << tmp.size();
      nAll += tmp.size();
      tmp = findUniqueSubgraphsOfLengthN(*mol, i);
      // std::cout << "\t" << tmp.size() << std::endl;;
      nUnique += tmp.size();
    }
    CHECK(nAll == 6435);
    CHECK(nUnique == 5618);
  }
}

TEST_CASE("testRootedSubgraphs", "[subgraphs]") {
  {
    auto mol = "CC1CC1"_smiles;
    REQUIRE(mol);

    PATH_LIST tmp;
    tmp = findAllSubgraphsOfLengthN(*mol, 1, false, 0);
    CHECK(tmp.size() == 1);
    tmp = findAllSubgraphsOfLengthN(*mol, 2, false, 0);
    CHECK(tmp.size() == 2);
    tmp = findAllSubgraphsOfLengthN(*mol, 3, false, 0);
    CHECK(tmp.size() == 3);
    tmp = findUniqueSubgraphsOfLengthN(*mol, 2, false, false, 0);
    CHECK(tmp.size() == 1);
    tmp = findUniqueSubgraphsOfLengthN(*mol, 3, false, false, 0);
    CHECK(tmp.size() == 2);

    INT_PATH_LIST_MAP tmpm;
    tmpm = findAllSubgraphsOfLengthsMtoN(*mol, 1, 3, false, 0);
    CHECK(tmpm[1].size() == 1);
    CHECK(tmpm[2].size() == 2);
    CHECK(tmpm[3].size() == 3);

    // edge case:
    tmp = findAllSubgraphsOfLengthN(*mol, 1, false, 10);
    CHECK(tmp.empty());
  }

  {  // tests for sf.net issue 250
    auto mol = "C1CC1C"_smiles;
    REQUIRE(mol);

    PATH_LIST tmp;
    tmp = findAllSubgraphsOfLengthN(*mol, 1, false, 3);
    CHECK(tmp.size() == 1);
    tmp = findAllSubgraphsOfLengthN(*mol, 2, false, 3);
    CHECK(tmp.size() == 2);
    tmp = findAllSubgraphsOfLengthN(*mol, 3, false, 3);
    CHECK(tmp.size() == 3);
    tmp = findUniqueSubgraphsOfLengthN(*mol, 2, false, false, 3);
    CHECK(tmp.size() == 1);
    tmp = findUniqueSubgraphsOfLengthN(*mol, 3, false, false, 3);
    CHECK(tmp.size() == 2);

    INT_PATH_LIST_MAP tmpm;
    tmpm = findAllSubgraphsOfLengthsMtoN(*mol, 1, 3, false, 3);
    CHECK(tmpm[1].size() == 1);
    CHECK(tmpm[2].size() == 2);
    CHECK(tmpm[3].size() == 3);

    // edge case:
    tmp = findAllSubgraphsOfLengthN(*mol, 1, false, 10);
    CHECK(tmp.empty());
  }

  {
    auto mol = "CC1CC1"_smiles;
    REQUIRE(mol);

    INT_PATH_LIST_MAP tmpm;
    tmpm = findAllSubgraphsOfLengthsMtoN(*mol, 1, 2, false, 0);
    CHECK(tmpm[1].size() == 1);
    CHECK(tmpm[2].size() == 2);
    CHECK(tmpm[3].empty());
    tmpm = findAllSubgraphsOfLengthsMtoN(*mol, 1, 3, false, 0);
    CHECK(tmpm[1].size() == 1);
    CHECK(tmpm[2].size() == 2);
    CHECK(tmpm[3].size() == 3);
  }
  {  // tests for sf.net issue 250
    auto mol = "C1CC1C"_smiles;
    REQUIRE(mol);

    INT_PATH_LIST_MAP tmpm;
    tmpm = findAllSubgraphsOfLengthsMtoN(*mol, 1, 2, false, 3);
    CHECK(tmpm[1].size() == 1);
    CHECK(tmpm[2].size() == 2);
    CHECK(tmpm[3].empty());
    tmpm = findAllSubgraphsOfLengthsMtoN(*mol, 1, 3, false, 3);
    CHECK(tmpm[1].size() == 1);
    CHECK(tmpm[2].size() == 2);
    CHECK(tmpm[3].size() == 3);
  }
}

TEST_CASE("testRootedPaths", "[subgraphs]") {
  {
    auto mol = "CC1CC1"_smiles;
    REQUIRE(mol);

    PATH_LIST tmp;

    // bond paths:
    tmp = findAllPathsOfLengthN(*mol, 1, true, false, 0);
    CHECK(tmp.size() == 1);
    tmp = findAllPathsOfLengthN(*mol, 2, true, false, 0);
    CHECK(tmp.size() == 2);
    tmp = findAllPathsOfLengthN(*mol, 3, true, false, 0);
    CHECK(tmp.size() == 2);

    // edge case:
    tmp = findAllPathsOfLengthN(*mol, 1, true, false, 10);
    CHECK(tmp.empty());

    // atom paths:
    tmp = findAllPathsOfLengthN(*mol, 1, false, false, 0);
    CHECK(tmp.size() == 1);
    tmp = findAllPathsOfLengthN(*mol, 2, false, false, 0);
    CHECK(tmp.size() == 1);
    tmp = findAllPathsOfLengthN(*mol, 3, false, false, 0);
    CHECK(tmp.size() == 2);
    tmp = findAllPathsOfLengthN(*mol, 4, false, false, 0);
    CHECK(tmp.size() == 2);
  }

  {  // tests for sf.net issue 250
    auto mol = "C1CC1C"_smiles;
    REQUIRE(mol);

    PATH_LIST tmp;

    // bond paths:
    tmp = findAllPathsOfLengthN(*mol, 1, true, false, 3);
    CHECK(tmp.size() == 1);
    tmp = findAllPathsOfLengthN(*mol, 2, true, false, 3);
    CHECK(tmp.size() == 2);
    tmp = findAllPathsOfLengthN(*mol, 3, true, false, 3);
    CHECK(tmp.size() == 2);

    // edge case:
    tmp = findAllPathsOfLengthN(*mol, 1, true, false, 10);
    CHECK(tmp.empty());

    // atom paths:
    tmp = findAllPathsOfLengthN(*mol, 1, false, false, 3);
    CHECK(tmp.size() == 1);
    tmp = findAllPathsOfLengthN(*mol, 2, false, false, 3);
    CHECK(tmp.size() == 1);
    tmp = findAllPathsOfLengthN(*mol, 3, false, false, 3);
    CHECK(tmp.size() == 2);
    tmp = findAllPathsOfLengthN(*mol, 4, false, false, 3);
    CHECK(tmp.size() == 2);
  }

  {
    auto mol = "CC1CC1"_smiles;
    REQUIRE(mol);

    INT_PATH_LIST_MAP tmpm;
    tmpm = findAllPathsOfLengthsMtoN(*mol, 1, 2, false, false, 0);
    CHECK(tmpm[1].size() == 1);
    CHECK(tmpm[2].size() == 1);
    CHECK(tmpm[3].empty());
    tmpm = findAllPathsOfLengthsMtoN(*mol, 1, 3, false, false, 0);
    CHECK(tmpm[1].size() == 1);
    CHECK(tmpm[2].size() == 1);
    CHECK(tmpm[3].size() == 2);
  }
  {  // tests for sf.net issue 250
    auto mol = "C1CC1C"_smiles;
    REQUIRE(mol);

    INT_PATH_LIST_MAP tmpm;
    tmpm = findAllPathsOfLengthsMtoN(*mol, 1, 2, false, false, 3);
    CHECK(tmpm[1].size() == 1);
    CHECK(tmpm[2].size() == 1);
    CHECK(tmpm[3].empty());
    tmpm = findAllPathsOfLengthsMtoN(*mol, 1, 3, false, false, 3);
    CHECK(tmpm[1].size() == 1);
    CHECK(tmpm[2].size() == 1);
    CHECK(tmpm[3].size() == 2);
  }
}
