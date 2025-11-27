//   Copyright (C) 2002-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

#include <algorithm>
#include <map>
#include <cmath>
#include <random>

#include <catch2/catch_all.hpp>

#include <RDGeneral/utils.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/test_fixtures.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace RDKit;
using namespace std;
RWMol _t;
typedef class ROMol Mol;

TEST_CASE("test1") {
  string smi;
  Mol *m;
  INT_VECT iv;
  unsigned int count;
  std::vector<ROMOL_SPTR> frags;
  std::vector<std::unique_ptr<ROMol>> otherFrags;

  smi = "CCCC(=O)O";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m, "smiles parse failed");
  count = MolOps::getMolFrags(*m, iv);
  CHECK_INVARIANT(count == 1, "bad frag count");
  frags = MolOps::getMolFrags(*m);
  CHECK_INVARIANT(frags.size() == 1, "bad frag count");
  REQUIRE(frags[0]->getNumAtoms() == 6);
  count = MolOps::getMolFrags(*m, otherFrags);
  CHECK_INVARIANT(count == 1, "bad frag count");
  CHECK_INVARIANT(otherFrags.size() == 1, "bad frag count");
  delete m;

  smi = "CCCC(=O)[O-].[Na+]";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m, "smiles parse failed");
  count = MolOps::getMolFrags(*m, iv);
  CHECK_INVARIANT(count == 2, "bad frag count");
  frags = MolOps::getMolFrags(*m);
  CHECK_INVARIANT(frags.size() == 2, "bad frag count");
  REQUIRE(frags[0]->getNumAtoms() == 6);
  REQUIRE(frags[1]->getNumAtoms() == 1);
  frags = MolOps::getMolFrags(*m, true, &iv);
  CHECK_INVARIANT(frags.size() == 2, "bad frag count");
  REQUIRE(frags[0]->getNumAtoms() == 6);
  REQUIRE(frags[1]->getNumAtoms() == 1);
  REQUIRE(iv.size() == 7);
  REQUIRE(iv[0] == 0);
  REQUIRE(iv[6] == 1);
  count = MolOps::getMolFrags(*m, otherFrags, true, &iv);
  CHECK_INVARIANT(count == 2, "bad frag count");
  CHECK_INVARIANT(otherFrags.size() == 2, "bad frag count");
  REQUIRE(frags[0]->getNumAtoms() == 6);
  REQUIRE(frags[1]->getNumAtoms() == 1);
  REQUIRE(iv.size() == 7);
  REQUIRE(iv[0] == 0);
  REQUIRE(iv[6] == 1);
  delete m;

  smi = "CCCC(=O)[O-].[Na+].[NH4+].[Cl-]";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m, "smiles parse failed");
  count = MolOps::getMolFrags(*m, iv);
  CHECK_INVARIANT(count == 4, "bad frag count");
  frags = MolOps::getMolFrags(*m);
  CHECK_INVARIANT(frags.size() == 4, "bad frag count");
  REQUIRE(frags[0]->getNumAtoms() == 6);
  REQUIRE(frags[1]->getNumAtoms() == 1);
  REQUIRE(frags[2]->getNumAtoms() == 1);
  REQUIRE(frags[3]->getNumAtoms() == 1);
  count = MolOps::getMolFrags(*m, otherFrags);
  CHECK_INVARIANT(count == 4, "bad frag count");
  CHECK_INVARIANT(otherFrags.size() == 4, "bad frag count");
  REQUIRE(otherFrags[0]->getNumAtoms() == 6);
  REQUIRE(otherFrags[1]->getNumAtoms() == 1);
  REQUIRE(otherFrags[2]->getNumAtoms() == 1);
  REQUIRE(otherFrags[3]->getNumAtoms() == 1);
  delete m;
};

TEST_CASE("test2") {
  string smi;
  Mol *m;
  INT_VECT iv;
  int count;
  smi = "CCCC(=O)O";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m, "smiles parse failed");
  count = MolOps::getMolFrags(*m, iv);
  CHECK_INVARIANT(count == 1, "bad frag count");
  CHECK_INVARIANT(iv[0] == 0, "bad frag membership");
  CHECK_INVARIANT(iv[1] == 0, "bad frag membership");
  CHECK_INVARIANT(iv[4] == 0, "bad frag membership");
  CHECK_INVARIANT(iv[5] == 0, "bad frag membership");
  delete m;

  smi = "CCCC(=O)[O-].[Na+]";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m, "smiles parse failed");
  count = MolOps::getMolFrags(*m, iv);
  CHECK_INVARIANT(count == 2, "bad frag count");
  CHECK_INVARIANT(iv[0] == 0, "bad frag membership");
  CHECK_INVARIANT(iv[1] == 0, "bad frag membership");
  CHECK_INVARIANT(iv[4] == 0, "bad frag membership");
  CHECK_INVARIANT(iv[5] == 0, "bad frag membership");
  CHECK_INVARIANT(iv[6] == 1, "bad frag membership");
  delete m;
};

TEST_CASE("test3") {
  string smi;
  Mol *m;
  VECT_INT_VECT sssr;
  INT_VECT rings;
  int count;

  smi = "C1CC1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getRingInfo()->numRings() == 1);
  count = MolOps::findSSSR(*m, sssr);
  REQUIRE(count == 1);
  REQUIRE(sssr[0].size() == 3);
  REQUIRE(m->getRingInfo()->numRings() == 1);

  for (unsigned int i = 0; i < m->getNumAtoms(); i++) {
    REQUIRE(m->getRingInfo()->isAtomInRingOfSize(i, 3));
    REQUIRE(!m->getRingInfo()->isAtomInRingOfSize(i, 4));
    REQUIRE(m->getRingInfo()->numAtomRings(i) == 1);
    REQUIRE(m->getRingInfo()->atomRingSizes(i) == (INT_VECT{3}));
    REQUIRE((m->getRingInfo()->atomMembers(i).size() == 1 &&
             m->getRingInfo()->atomMembers(i).at(0) == 0));
  }
  for (unsigned int i = 0; i < m->getNumBonds(); i++) {
    REQUIRE(m->getRingInfo()->isBondInRingOfSize(i, 3));
    REQUIRE(!m->getRingInfo()->isBondInRingOfSize(i, 4));
    REQUIRE(m->getRingInfo()->numBondRings(i) == 1);
    REQUIRE(m->getRingInfo()->bondRingSizes(i) == (INT_VECT{3}));
    REQUIRE((m->getRingInfo()->bondMembers(i).size() == 1 &&
             m->getRingInfo()->bondMembers(i).at(0) == 0));
  }
  REQUIRE(m->getRingInfo()->areAtomsInSameRing(0, 1));
  REQUIRE(m->getRingInfo()->areAtomsInSameRingOfSize(0, 1, 3));
  REQUIRE(!m->getRingInfo()->areAtomsInSameRingOfSize(0, 1, 4));
  REQUIRE(m->getRingInfo()->areBondsInSameRing(0, 1));
  REQUIRE(m->getRingInfo()->areBondsInSameRingOfSize(0, 1, 3));
  REQUIRE(!m->getRingInfo()->areBondsInSameRingOfSize(0, 1, 4));
  BOOST_LOG(rdInfoLog) << smi << "\n";
  delete m;

  smi = "C1CCC1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getRingInfo()->numRings() == 1);
  count = MolOps::findSSSR(*m, sssr);
  REQUIRE(count == 1);
  REQUIRE(sssr[0].size() == 4);
  REQUIRE(m->getRingInfo()->numRings() == 1);
  for (unsigned int i = 0; i < m->getNumAtoms(); i++) {
    REQUIRE(m->getRingInfo()->isAtomInRingOfSize(i, 4));
    REQUIRE(!m->getRingInfo()->isAtomInRingOfSize(i, 3));
    REQUIRE(m->getRingInfo()->numAtomRings(i) == 1);
    REQUIRE(m->getRingInfo()->atomRingSizes(i) == (INT_VECT{4}));
    REQUIRE((m->getRingInfo()->atomMembers(i).size() == 1 &&
             m->getRingInfo()->atomMembers(i).at(0) == 0));
  }
  REQUIRE(m->getRingInfo()->isBondInRingOfSize(0, 4));
  REQUIRE(m->getRingInfo()->numBondRings(0) == 1);
  REQUIRE(m->getRingInfo()->bondRingSizes(0) == (INT_VECT{4}));
  REQUIRE(m->getRingInfo()->areBondsInSameRing(0, 1));
  REQUIRE(m->getRingInfo()->areBondsInSameRingOfSize(0, 1, 4));
  REQUIRE(!m->getRingInfo()->areBondsInSameRingOfSize(0, 1, 5));

  BOOST_LOG(rdInfoLog) << smi << "\n";
  delete m;

  smi = "C1CCCCCC1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getRingInfo()->numRings() == 1);
  count = MolOps::findSSSR(*m, sssr);
  REQUIRE(count == 1);
  REQUIRE(sssr[0].size() == 7);
  REQUIRE(m->getRingInfo()->numRings() == 1);
  for (unsigned int i = 0; i < m->getNumAtoms(); i++) {
    REQUIRE(m->getRingInfo()->isAtomInRingOfSize(i, 7));
    REQUIRE(m->getRingInfo()->numAtomRings(i) == 1);
  }
  REQUIRE(m->getRingInfo()->isBondInRingOfSize(0, 7));
  REQUIRE(m->getRingInfo()->numBondRings(0) == 1);

  BOOST_LOG(rdInfoLog) << smi << "\n";
  delete m;

  smi = "C1C(CCC)CC(C(C)CCC(CC))CCC1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  count = MolOps::findSSSR(*m, sssr);
  REQUIRE(count == 1);
  REQUIRE(sssr[0].size() == 7);
  REQUIRE(m->getRingInfo()->numAtomRings(0) == 1);
  REQUIRE(
      m->getRingInfo()->numBondRings(m->getBondBetweenAtoms(0, 1)->getIdx()));
  REQUIRE(m->getRingInfo()->atomMembers(0).size() == 1);
  REQUIRE(m->getRingInfo()->atomMembers(2).empty());
  REQUIRE(m->getRingInfo()->atomRingSizes(2).empty());
  REQUIRE(m->getRingInfo()->atomRingSizes(5) == (INT_VECT{7}));
  REQUIRE(m->getRingInfo()->atomRingSizes(1) == (INT_VECT{7}));
  REQUIRE(m->getRingInfo()->atomRingSizes(99).empty());
  REQUIRE(m->getRingInfo()->areAtomsInSameRing(0, 1));
  REQUIRE(m->getRingInfo()->areAtomsInSameRingOfSize(0, 1, 7));
  REQUIRE(!m->getRingInfo()->areAtomsInSameRingOfSize(0, 1, 5));
  REQUIRE(!m->getRingInfo()->areAtomsInSameRing(1, 2));
  REQUIRE(!m->getRingInfo()->areAtomsInSameRingOfSize(1, 2, 4));
  REQUIRE(m->getRingInfo()->bondRingSizes(1).empty());
  REQUIRE(m->getRingInfo()->bondRingSizes(4) == (INT_VECT{7}));
  REQUIRE(m->getRingInfo()->bondRingSizes(5) == (INT_VECT{7}));
  REQUIRE(m->getRingInfo()->bondRingSizes(99).empty());
  REQUIRE(m->getRingInfo()->bondMembers(0).size() == 1);
  REQUIRE(m->getRingInfo()->bondMembers(1).empty());
  REQUIRE(!m->getRingInfo()->isRingFused(0));
  REQUIRE(m->getRingInfo()->numFusedBonds(0) == 0);
  REQUIRE(
      !m->getRingInfo()->numBondRings(m->getBondBetweenAtoms(1, 2)->getIdx()));
  BOOST_LOG(rdInfoLog) << smi << "\n";
  delete m;

  smi = "CC1CCC1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  count = MolOps::findSSSR(*m, sssr);
  REQUIRE(count == 1);
  REQUIRE(sssr[0].size() == 4);
  REQUIRE(!m->getBondBetweenAtoms(0, 1)->hasProp(
      common_properties::ringMembership));
  REQUIRE(
      !m->getRingInfo()->numBondRings(m->getBondBetweenAtoms(0, 1)->getIdx()));
  REQUIRE(
      m->getRingInfo()->numBondRings(m->getBondBetweenAtoms(1, 2)->getIdx()));
  BOOST_LOG(rdInfoLog) << smi << "\n";
  delete m;

  smi = "CC1C(C2)CCC2C1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  count = MolOps::findSSSR(*m, sssr);
  REQUIRE(count == 2);
  REQUIRE(sssr[0].size() == 5);
  REQUIRE(sssr[1].size() == 5);
  REQUIRE(!m->getRingInfo()->isAtomInRingOfSize(0, 5));
  REQUIRE(m->getRingInfo()->numAtomRings(0) == 0);
  REQUIRE(m->getRingInfo()->isAtomInRingOfSize(1, 5));
  REQUIRE(m->getRingInfo()->numAtomRings(1) == 1);
  REQUIRE(m->getRingInfo()->isAtomInRingOfSize(2, 5));
  REQUIRE(m->getRingInfo()->numAtomRings(2) == 2);
  REQUIRE(m->getRingInfo()->atomRingSizes(0).empty());
  REQUIRE(m->getRingInfo()->atomRingSizes(1) == (INT_VECT{5}));
  REQUIRE(m->getRingInfo()->atomRingSizes(2) == (INT_VECT{5, 5}));
  REQUIRE(m->getRingInfo()->atomMembers(2).size() == 2);
  REQUIRE(m->getRingInfo()->atomMembers(2).at(0) == 0);
  REQUIRE(m->getRingInfo()->atomMembers(2).at(1) == 1);
  REQUIRE(m->getRingInfo()->isRingFused(0));
  REQUIRE(m->getRingInfo()->isRingFused(1));
  REQUIRE(m->getRingInfo()->areRingsFused(0, 1));
  REQUIRE(m->getRingInfo()->numFusedBonds(0) == 2);
  REQUIRE(m->getRingInfo()->numFusedBonds(1) == 2);
  delete m;

  smi = "C(C1C2C3C41)(C2C35)C45";  // cubane
  // smi = "C1(C2C3C4C5C6C72)C3C4C5C6C71"; // from Figureras paper
  // smi = "C17C5C4C3C2C1C6C2C3C4C5C67";
  // we cannot use the sanitization code, because that finds *symmetric*
  // rings, which will break this case:
  m = SmilesToMol(smi, 0, 0);
  int bfs = MolOps::findSSSR(*m);
  REQUIRE(bfs == 5);
  BOOST_LOG(rdInfoLog) << "BFSR: " << bfs << "\n";
  VECT_INT_VECT bfrs;
  bfrs.resize(0);
  bfs = MolOps::symmetrizeSSSR(*m, bfrs);
  REQUIRE(bfs == 6);
  BOOST_LOG(rdInfoLog) << "BFSR: " << bfs << "\n";
  // VECT_INT_VECT_I ri;
  // for (ri == bfrs.begin(); ri != bfrs.end(); ri++) {
  for (auto bring : bfrs) {
    INT_VECT_I mi;
    BOOST_LOG(rdInfoLog) << "( ";
    // for (mi = (*ri).begin(); mi != (*ri).end(); mi++) {
    for (mi = bring.begin(); mi != bring.end(); mi++) {
      BOOST_LOG(rdInfoLog) << " " << (*mi);
    }
    BOOST_LOG(rdInfoLog) << ")\n";
  }
  BOOST_LOG(rdInfoLog) << smi << "\n";
  for (unsigned int i = 0; i < m->getRingInfo()->numRings(); ++i) {
    REQUIRE(m->getRingInfo()->isRingFused(i));
    REQUIRE(m->getRingInfo()->numFusedBonds(i) == 4);
  }

  delete m;

  // Figueras figure 4:
  //  * The third ring is bigger, and shouldn't be accessed in symmetrizeSSSR
  smi = "C12CC(CC2)CC1";
  m = SmilesToMol(smi, 0, 0);
  bfs = MolOps::findSSSR(*m);
  REQUIRE(bfs == 2);
  bfrs.resize(0);
  bfs = MolOps::symmetrizeSSSR(*m, bfrs);
  REQUIRE(bfs == 2);
  delete m;

  // Counterexamples in ring perception figure 4:
  //  * The native Figueras algorithm cannot work on this molecule, it will
  //    fail after finding one ring. Naive modified Figueras finds a 6 membered
  //    ring, which is wrong.
  smi = "C123C4C5C6(C3)C7C1C8C2C4C5C6C78";
  m = SmilesToMol(smi, 0, 0);
  bfs = MolOps::findSSSR(*m);
  REQUIRE(bfs == 7);
  bfrs.resize(0);
  bfs = MolOps::symmetrizeSSSR(*m, bfrs);
  REQUIRE(bfs == 8);
  for (auto bring : bfrs) {
    REQUIRE(bring.size() < 6);
  }
  delete m;

  smi = "C1CC2C1CCC2";
  m = SmilesToMol(smi);
  REQUIRE(m);
  count = MolOps::findSSSR(*m, sssr);
  REQUIRE(count == 2);
  REQUIRE(sssr[0].size() == 4);
  REQUIRE(sssr[1].size() == 5);
  BOOST_LOG(rdInfoLog) << smi << "\n";
  delete m;

  smi = "C12=C3C=CC=C1C=CC2=CC=C3";
  BOOST_LOG(rdInfoLog) << "\n" << smi << "\n";
  m = SmilesToMol(smi, 0, 0);
  REQUIRE(m);
  count = MolOps::findSSSR(*m, sssr);
  REQUIRE(count == 3);
  REQUIRE(sssr[0].size() == 6);
  REQUIRE(sssr[1].size() == 5);
  REQUIRE(sssr[2].size() == 6);
  BOOST_LOG(rdInfoLog) << smi << "\n";
  delete m;

  smi = "C1(O)C(O)C(O)C1O";
  m = SmilesToMol(smi, 0, 0);
  REQUIRE(m);
  count = MolOps::findSSSR(*m, sssr);
  REQUIRE(count == 1);
  REQUIRE(sssr[0].size() == 4);
  for (unsigned i = 0; i < m->getNumAtoms(); i++) {
    if (!(i % 2)) {
      REQUIRE(m->getRingInfo()->numAtomRings(i) == 1);
      REQUIRE(m->getRingInfo()->isAtomInRingOfSize(i, 4));
    } else {
      REQUIRE(m->getRingInfo()->numAtomRings(i) == 0);
    }
  }
  BOOST_LOG(rdInfoLog) << smi << "\n";
  delete m;

  // this molecule is from issue 134
  // it should come up with three rings
  smi = "SC(C3C1CC(C3)CC(C2S)(O)C1)2S";
  m = SmilesToMol(smi, 0, 0);
  REQUIRE(m);
  count = MolOps::findSSSR(*m, sssr);
  REQUIRE(count == 3);
  REQUIRE(sssr[0].size() == 5);
  REQUIRE(sssr[1].size() == 6);
  REQUIRE(sssr[2].size() == 6);
  delete m;

  // this yet another painful case
  smi = "CC1=CC=C(C=C1)S(=O)(=O)O[CH]2[CH]3CO[CH](O3)[CH]4OC(C)(C)O[CH]24";
  m = SmilesToMol(smi, 0, 0);
  REQUIRE(m);
  count = MolOps::findSSSR(*m, sssr);
  REQUIRE(count == 4);
  REQUIRE(sssr[0].size() == 6);
  REQUIRE(sssr[1].size() == 5);
  REQUIRE(sssr[2].size() == 5);
  REQUIRE(sssr[3].size() == 6);
  delete m;

  smi = "C1CC2C1C2";
  m = SmilesToMol(smi, 0, 0);
  REQUIRE(m);
  count = MolOps::findSSSR(*m, sssr);
  REQUIRE(count == 2);
  REQUIRE(sssr[0].size() == 4);
  REQUIRE(sssr[1].size() == 3);

  REQUIRE(m->getRingInfo()->numAtomRings(0) == 1);
  REQUIRE(m->getRingInfo()->isAtomInRingOfSize(0, 4));

  REQUIRE(m->getRingInfo()->numAtomRings(1) == 1);
  REQUIRE(m->getRingInfo()->isAtomInRingOfSize(1, 4));

  REQUIRE(m->getRingInfo()->numAtomRings(2) == 2);
  REQUIRE(m->getRingInfo()->isAtomInRingOfSize(2, 3));
  REQUIRE(m->getRingInfo()->isAtomInRingOfSize(2, 4));

  REQUIRE(m->getRingInfo()->numAtomRings(3) == 2);
  REQUIRE(m->getRingInfo()->isAtomInRingOfSize(3, 4));
  REQUIRE(m->getRingInfo()->isAtomInRingOfSize(3, 3));

  REQUIRE(m->getRingInfo()->numAtomRings(4) == 1);
  REQUIRE(!m->getRingInfo()->isAtomInRingOfSize(4, 4));
  REQUIRE(m->getRingInfo()->isAtomInRingOfSize(4, 3));

  REQUIRE(m->getRingInfo()->atomMembers(2).size() == 2);
  REQUIRE(m->getRingInfo()->areAtomsInSameRing(2, 3));
  REQUIRE(!m->getRingInfo()->areAtomsInSameRing(1, 4));
  REQUIRE(m->getRingInfo()->areAtomsInSameRingOfSize(2, 3, 3));
  REQUIRE(m->getRingInfo()->areAtomsInSameRingOfSize(2, 3, 4));
  REQUIRE(!m->getRingInfo()->areAtomsInSameRingOfSize(2, 3, 5));
  REQUIRE(m->getRingInfo()->bondMembers(2).size() == 2);
  REQUIRE(m->getRingInfo()->bondMembers(0).size() == 1);
  REQUIRE(m->getRingInfo()->areBondsInSameRing(1, 2));
  REQUIRE(m->getRingInfo()->areBondsInSameRing(2, 5));
  REQUIRE(!m->getRingInfo()->areBondsInSameRing(1, 3));
  REQUIRE(m->getRingInfo()->areBondsInSameRingOfSize(1, 2, 4));
  REQUIRE(m->getRingInfo()->areBondsInSameRingOfSize(2, 5, 3));
  REQUIRE(!m->getRingInfo()->areBondsInSameRingOfSize(1, 2, 3));
  REQUIRE(!m->getRingInfo()->areBondsInSameRingOfSize(1, 3, 4));
  REQUIRE(m->getRingInfo()->isRingFused(0));
  REQUIRE(m->getRingInfo()->isRingFused(1));
  REQUIRE(m->getRingInfo()->areRingsFused(0, 1));
  REQUIRE(m->getRingInfo()->numFusedBonds(0) == 1);
  REQUIRE(m->getRingInfo()->numFusedBonds(1) == 1);
  delete m;

  // This is a test of Issue 217
  smi = "C=C1C2CC1C2";
  m = SmilesToMol(smi, 0, 0);
  REQUIRE(m);
  count = MolOps::findSSSR(*m, sssr);
  REQUIRE(count == 2);
  REQUIRE(sssr[0].size() == 4);
  REQUIRE(sssr[1].size() == 4);
  count = MolOps::symmetrizeSSSR(*m, sssr);
  REQUIRE(count == 3);
  REQUIRE(sssr[0].size() == 4);
  REQUIRE(sssr[1].size() == 4);
  REQUIRE(sssr[2].size() == 4);

  REQUIRE(m->getRingInfo()->numAtomRings(0) == 0);
  REQUIRE(m->getRingInfo()->numAtomRings(1) == 2);
  REQUIRE(m->getRingInfo()->isAtomInRingOfSize(1, 4));
  REQUIRE(m->getRingInfo()->numAtomRings(2) == 3);
  REQUIRE(m->getRingInfo()->isAtomInRingOfSize(2, 4));
  REQUIRE(m->getRingInfo()->numAtomRings(3) == 2);
  REQUIRE(m->getRingInfo()->isAtomInRingOfSize(3, 4));
  REQUIRE(m->getRingInfo()->numAtomRings(4) == 3);
  REQUIRE(m->getRingInfo()->isAtomInRingOfSize(4, 4));
  delete m;
}

TEST_CASE("test4") {
  auto m = "C=C"_smiles;
  REQUIRE(m);
  double *adjMat = MolOps::getAdjacencyMatrix(*m);
  REQUIRE(adjMat);
  REQUIRE(adjMat[0] == 0);
  REQUIRE(adjMat[1] == 1);
  REQUIRE(adjMat[2] == 1);
  REQUIRE(adjMat[3] == 0);
  adjMat = MolOps::getAdjacencyMatrix(*m);
  REQUIRE(adjMat);
  REQUIRE(adjMat[0] == 0);
  REQUIRE(adjMat[1] == 1);
  REQUIRE(adjMat[2] == 1);
  REQUIRE(adjMat[3] == 0);
  bool useBO = true;
  adjMat = MolOps::getAdjacencyMatrix(*m, useBO);
  REQUIRE(adjMat);
  REQUIRE(adjMat[0] == 0);
  REQUIRE(adjMat[1] == 2.0);
  REQUIRE(adjMat[2] == 2.0);
  REQUIRE(adjMat[3] == 0);
}

TEST_CASE("test5") {
  string smi;
  Mol *m;
  VECT_INT_VECT sssr;

  int count;
  smi = "C1C4C5C3C(=O)C2C5C1C2C34";
  m = SmilesToMol(smi, 0, 0);
  count = MolOps::findSSSR(*m, sssr);
  BOOST_LOG(rdInfoLog) << "Count: " << count << "\n";
  CHECK_INVARIANT(count == 5, "");
  delete m;

  smi = "C1C(C2)CCC2C1";
  m = SmilesToMol(smi);
  count = MolOps::findSSSR(*m, sssr);
  CHECK_INVARIANT(count == 2, "");
  delete m;
}

TEST_CASE("Testing Hydrogen Ops") {
  ROMol *m, *m2, *m3;
  INT_VECT tree;

  std::string smi = "CCC";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m, "");
  CHECK_INVARIANT(m->getNumAtoms() == 3, "");

  m2 = MolOps::addHs(*m);
  CHECK_INVARIANT(m2->getNumAtoms() == 11, "");

  // addHs should not set the noImplicit flag.
  // This was Github Issue #7123
  for (auto at : m2->atoms()) {
    REQUIRE(at->getNoImplicit() == false);
  }

  delete m;
  delete m2;

  smi = "CC(=O)[OH]";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m, "");
  CHECK_INVARIANT(m->getNumAtoms() == 4, "");

  m2 = MolOps::addHs(*m, true);
  CHECK_INVARIANT(m2->getNumAtoms() == 5, "");

  m3 = MolOps::addHs(*m2, false);
  CHECK_INVARIANT(m3->getNumAtoms() == 8, "");
  delete m2;
  delete m3;

  m2 = MolOps::addHs(*m, false);
  CHECK_INVARIANT(m2->getNumAtoms() == 8, "");
  delete m;
  // remove all
  // rdInfoLog) << "5" << std::endl;
  m3 = MolOps::removeHs(*m2, false);
  CHECK_INVARIANT(m3->getNumAtoms() == 4, "");
  delete m3;

  // remove only implicit
  m3 = MolOps::removeHs(*m2, true);
  CHECK_INVARIANT(m3->getNumAtoms() == 5, "");

  MolOps::removeHs(static_cast<RWMol &>(*m3), false);
  CHECK_INVARIANT(m3->getNumAtoms() == 4, "");
  delete m2;
  delete m3;

  // this test is also done in the same order in the python tests:
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m, "");
  CHECK_INVARIANT(m->getNumAtoms() == 4, "");
  m2 = MolOps::addHs(*m, true);
  CHECK_INVARIANT(m2->getNumAtoms() == 5, "");

  m3 = MolOps::removeHs(*m2, true);
  CHECK_INVARIANT(m3->getNumAtoms() == 5, "");
  delete m3;

  m3 = MolOps::removeHs(*m2, false);
  CHECK_INVARIANT(m3->getNumAtoms() == 4, "");
  delete m2;
  delete m3;

  m2 = MolOps::addHs(*m, false);
  CHECK_INVARIANT(m2->getNumAtoms() == 8, "");

  m3 = MolOps::removeHs(*m2, true);
  CHECK_INVARIANT(m3->getNumAtoms() == 5, "");
  delete m3;

  m3 = MolOps::removeHs(*m2, false);
  CHECK_INVARIANT(m3->getNumAtoms() == 4, "");
  delete m3;
  delete m2;
  delete m;

  // related to RDTrack Issues 109 and 110:
  smi =
      "C1C=C([C@H](N)C(=O)N[C@@]2([H])[C@]3([H])SC(C)(C)[C@@H](C(=O)O)N3C(=O)2)"
      "C=CC=1";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m, "");
  CHECK_INVARIANT(m->getNumAtoms() == 24, "");
  m3 = MolOps::removeHs(*m, false);
  CHECK_INVARIANT(m3->getNumAtoms() == 24, "");
  delete m;
  delete m3;

  // RDTrack Issue 130:
  smi = "[H][N+]([H])([H])[H]";
  m = SmilesToMol(smi, false, false);
  CHECK_INVARIANT(m, "");
  CHECK_INVARIANT(m->getNumAtoms() == 5, "");

  m2 = MolOps::removeHs(*m, 0, false);
  CHECK_INVARIANT(m2->getNumAtoms() == 1, "");
  delete m;
  delete m2;

  smi = "[H][N+]([H])([H])[H]";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m, "");
  CHECK_INVARIANT(m->getNumAtoms() == 1, "");

  delete m;
  smi = "[H][H]";
  m = SmilesToMol(smi, false, false);
  CHECK_INVARIANT(m, "");
  CHECK_INVARIANT(m->getNumAtoms() == 2, "");

  m2 = MolOps::removeHs(*m, 0, false);
  CHECK_INVARIANT(m2->getNumAtoms() == 2, "");
  delete m;
  delete m2;

  std::string sma;
  smi = "C-C";
  m = SmartsToMol(smi);
  MolOps::sanitizeMol(*((RWMol *)m));
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 2);
  sma = SmartsWrite::GetAtomSmarts(
      static_cast<const QueryAtom *>(m->getAtomWithIdx(0)));
  REQUIRE(sma == "C");

  m2 = MolOps::addHs(*m);
  REQUIRE(m2->getNumAtoms() == 8);
  sma = SmartsWrite::GetAtomSmarts(
      static_cast<const QueryAtom *>(m2->getAtomWithIdx(0)));
  REQUIRE(sma == "C");
  delete m;

  m = MolOps::mergeQueryHs(*m2);
  REQUIRE(m->getNumAtoms() == 2);
  sma = SmartsWrite::GetAtomSmarts(
      static_cast<const QueryAtom *>(m->getAtomWithIdx(0)));

  // this was sf.net issue 3415204:
  REQUIRE(sma == "[C&!H0&!H1&!H2]");
  delete m;
  delete m2;

  // RDTrack Issue 1228:
  smi = "c1c[nH]cc1";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m, "");
  CHECK_INVARIANT(m->getNumAtoms() == 5, "");

  m2 = MolOps::addHs(*m, false, false);
  CHECK_INVARIANT(m2->getNumAtoms() == 10, "");

  delete m;
  m = MolOps::removeHs(*m2);
  CHECK_INVARIANT(m, "");
  CHECK_INVARIANT(m->getNumAtoms() == 5, "");
  delete m;
  delete m2;

  // labelling:
  smi = "c1cn([H])cc1";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m, "");
  CHECK_INVARIANT(m->getNumAtoms() == 5, "");

  m2 = MolOps::removeHs(*m);
  CHECK_INVARIANT(m, "");
  CHECK_INVARIANT(m->getNumAtoms() == 5, "");
  delete m;
  delete m2;

  smi = "c1cn([2H])cc1";
  m = SmilesToMol(smi);
  CHECK_INVARIANT(m, "");
  CHECK_INVARIANT(m->getNumAtoms() == 6, "");

  m2 = MolOps::removeHs(*m);
  CHECK_INVARIANT(m, "");
  CHECK_INVARIANT(m->getNumAtoms() == 6, "");
  delete m;
  delete m2;

  smi = "CC[H]";
  m = SmilesToMol(smi, false, false);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 3);
  m2 = MolOps::mergeQueryHs(*m);
  REQUIRE(m2->getNumAtoms() == 2);
  REQUIRE(!m->getAtomWithIdx(1)->hasQuery());
  REQUIRE(m2->getAtomWithIdx(1)->hasQuery());
  delete m;
  delete m2;

  // sf.net issue 3415206
  smi = "CO[H]";
  m = SmilesToMol(smi, false, false);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 3);
  m2 = MolOps::mergeQueryHs(*m);
  REQUIRE(m2->getNumAtoms() == 2);
  sma = SmartsWrite::GetAtomSmarts(
      static_cast<const QueryAtom *>(m2->getAtomWithIdx(1)));

  REQUIRE(sma == "[#8&!H0]");
  delete m;
  delete m2;

  smi = "CN([H])[H]";
  m = SmilesToMol(smi, false, false);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 4);
  m2 = MolOps::mergeQueryHs(*m);
  REQUIRE(m2->getNumAtoms() == 2);
  sma = SmartsWrite::GetAtomSmarts(
      static_cast<const QueryAtom *>(m2->getAtomWithIdx(1)));

  REQUIRE(sma == "[#7&!H0&!H1]");
  delete m;
  delete m2;

  {
    // test the onlyOnAtoms option (github #758)
    std::string smi = "CCC";
    m = SmilesToMol(smi);
    CHECK_INVARIANT(m, "");
    CHECK_INVARIANT(m->getNumAtoms() == 3, "");

    UINT_VECT onlyOn;
    onlyOn.push_back(0);
    onlyOn.push_back(2);
    m2 = MolOps::addHs(*m, false, false, &onlyOn);
    CHECK_INVARIANT(m2->getNumAtoms() == 9, "");
    CHECK_INVARIANT(m2->getAtomWithIdx(0)->getDegree() == 4, "");
    CHECK_INVARIANT(m2->getAtomWithIdx(1)->getDegree() == 2, "");
    CHECK_INVARIANT(m2->getAtomWithIdx(2)->getDegree() == 4, "");
    delete m;
    delete m2;
  }
}

TEST_CASE("Testing Distance Matrix Operations") {
  auto m = "CC=O"_smiles;
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 3);

  bool useBO = false;
  bool useAtomWts = false;
  double *dMat;
  dMat = MolOps::getDistanceMat(*m, useBO, useAtomWts);
  REQUIRE(dMat);
  REQUIRE(dMat[0] == 0.0);
  REQUIRE(dMat[1] == 1.0);
  REQUIRE(dMat[2] == 2.0);
  REQUIRE(dMat[3] == 1.0);
  REQUIRE(dMat[4] == 0.0);
  REQUIRE(dMat[5] == 1.0);
  REQUIRE(dMat[6] == 2.0);
  REQUIRE(dMat[7] == 1.0);
  REQUIRE(dMat[8] == 0.0);

  dMat = MolOps::getDistanceMat(*m, useBO, useAtomWts);
  REQUIRE(dMat);
  REQUIRE(dMat[0] == 0.0);
  REQUIRE(dMat[1] == 1.0);
  REQUIRE(dMat[2] == 2.0);
  REQUIRE(dMat[3] == 1.0);
  REQUIRE(dMat[4] == 0.0);
  REQUIRE(dMat[5] == 1.0);
  REQUIRE(dMat[6] == 2.0);
  REQUIRE(dMat[7] == 1.0);
  REQUIRE(dMat[8] == 0.0);

  // test Issue328:
  useBO = true;
  dMat = MolOps::getDistanceMat(*m, useBO, useAtomWts);
  REQUIRE(dMat);
  REQUIRE(dMat[0] == 0.0);
  REQUIRE(dMat[1] == 1.0);
  REQUIRE(dMat[2] == 1.5);
  REQUIRE(dMat[3] == 1.0);
  REQUIRE(dMat[4] == 0.0);
  REQUIRE(dMat[5] == 0.5);
  REQUIRE(dMat[6] == 1.5);
  REQUIRE(dMat[7] == 0.5);
  REQUIRE(dMat[8] == 0.0);

  useBO = false;
  dMat = MolOps::getDistanceMat(*m, useBO, useAtomWts);
  REQUIRE(dMat);
  REQUIRE(dMat[0] == 0.0);
  REQUIRE(dMat[1] == 1.0);
  REQUIRE(dMat[2] == 2.0);
  REQUIRE(dMat[3] == 1.0);
  REQUIRE(dMat[4] == 0.0);
  REQUIRE(dMat[5] == 1.0);
  REQUIRE(dMat[6] == 2.0);
  REQUIRE(dMat[7] == 1.0);
  REQUIRE(dMat[8] == 0.0);

  useBO = false;
  useAtomWts = true;
  dMat = MolOps::getDistanceMat(*m, useBO, useAtomWts);
  REQUIRE(dMat);
  for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
    for (unsigned int j = 0; j < m->getNumAtoms(); ++j) {
      std::cerr << dMat[i * m->getNumAtoms() + j] << " ";
    }
    std::cerr << std::endl;
  }
  REQUIRE(dMat[0] == 1.0);
  REQUIRE(dMat[1] == 1.0);
  REQUIRE(dMat[2] == 2.0);
  REQUIRE(dMat[3] == 1.0);
  REQUIRE(dMat[4] == 1.0);
  REQUIRE(dMat[5] == 1.0);
  REQUIRE(dMat[6] == 2.0);
  REQUIRE(dMat[7] == 1.0);
  REQUIRE(dMat[8] == 6. / 8.);

  useBO = true;
  useAtomWts = true;
  dMat = MolOps::getDistanceMat(*m, useBO, useAtomWts);
  REQUIRE(dMat);
  REQUIRE(dMat[0] == 1.0);
  REQUIRE(dMat[1] == 1.0);
  REQUIRE(dMat[2] == 1.5);
  REQUIRE(dMat[3] == 1.0);
  REQUIRE(dMat[4] == 1.0);
  REQUIRE(dMat[5] == 0.5);
  REQUIRE(dMat[6] == 1.5);
  REQUIRE(dMat[7] == 0.5);
  REQUIRE(dMat[8] == 6. / 8.);

  useBO = false;
  useAtomWts = false;
  dMat = MolOps::getDistanceMat(*m, useBO, useAtomWts);
  REQUIRE(dMat);
  REQUIRE(dMat[0] == 0.0);
  REQUIRE(dMat[1] == 1.0);
  REQUIRE(dMat[2] == 2.0);
  REQUIRE(dMat[3] == 1.0);
  REQUIRE(dMat[4] == 0.0);
  REQUIRE(dMat[5] == 1.0);
  REQUIRE(dMat[6] == 2.0);
  REQUIRE(dMat[7] == 1.0);
  REQUIRE(dMat[8] == 0.0);

  // limit participating atoms and bonds
  std::vector<int> activeAtoms = {1, 2};
  std::vector<const Bond *> activeBonds = {m->getBondWithIdx(1)};
  useBO = false;
  useAtomWts = false;
  std::unique_ptr<double[]> dMat2{
      MolOps::getDistanceMat(*m, activeAtoms, activeBonds, useBO, useAtomWts)};
  REQUIRE(dMat2);
  REQUIRE(dMat2[0] == 0.0);
  REQUIRE(dMat2[1] == 1.0);
  REQUIRE(dMat2[2] == 1.0);
  REQUIRE(dMat2[3] == 0.0);

  useBO = true;
  useAtomWts = true;
  dMat2.reset(
      MolOps::getDistanceMat(*m, activeAtoms, activeBonds, useBO, useAtomWts));
  REQUIRE(dMat2);
  REQUIRE(dMat2[0] == 1.0);
  REQUIRE(dMat2[1] == 0.5);
  REQUIRE(dMat2[2] == 0.5);
  REQUIRE(dMat2[3] == 6.0 / 8.0);
}

TEST_CASE("Testing Atom Ranking") {
  ROMol *m;
  std::string smi = "FC(Cl)(Br)C";
  m = SmilesToMol(smi);
  REQUIRE(m);

  UINT_VECT ranks;
  ranks.resize(m->getNumAtoms());
  Chirality::assignAtomCIPRanks(*m, ranks);

  unsigned int cip1, cip2;
  REQUIRE(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPRank));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPRank, cip1);
  REQUIRE(cip1 == ranks[0]);
  REQUIRE(m->getAtomWithIdx(2)->hasProp(common_properties::_CIPRank));
  m->getAtomWithIdx(2)->getProp(common_properties::_CIPRank, cip2);
  REQUIRE(cip2 == ranks[2]);
  REQUIRE(cip1 < cip2);
  REQUIRE(m->getAtomWithIdx(4)->hasProp(common_properties::_CIPRank));
  m->getAtomWithIdx(4)->getProp(common_properties::_CIPRank, cip2);
  REQUIRE(cip1 > cip2);
  m->getAtomWithIdx(2)->getProp(common_properties::_CIPRank, cip1);
  REQUIRE(m->getAtomWithIdx(3)->hasProp(common_properties::_CIPRank));
  m->getAtomWithIdx(3)->getProp(common_properties::_CIPRank, cip2);
  REQUIRE(cip1 < cip2);

  delete m;
  smi = "FC(Cl)(Br)C(F)(F)F";
  m = SmilesToMol(smi);
  REQUIRE(m);
  ranks.resize(m->getNumAtoms());
  Chirality::assignAtomCIPRanks(*m, ranks);
  for (unsigned int i = 0; i < m->getNumAtoms(); i++) {
    unsigned int cip;
    REQUIRE(m->getAtomWithIdx(i)->hasProp(common_properties::_CIPRank));
    m->getAtomWithIdx(i)->getProp(common_properties::_CIPRank, cip);
  }
  delete m;
}

TEST_CASE("Testing CIP chirality") {
  ROMol *m;
  std::string cip;
  std::string smi = "F[C@]([C@])(Cl)Br";

  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  // make sure the cleanup worked:
  REQUIRE(m->getAtomWithIdx(2)->getChiralTag() == Atom::CHI_UNSPECIFIED);

  REQUIRE(!(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode)));
  REQUIRE(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  REQUIRE(!(m->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode)));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");

  delete m;
  smi = "F[C@H](C)C";
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  REQUIRE(!(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode)));
  REQUIRE(!(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode)));
  REQUIRE(!(m->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode)));
  // test Issue 194:
  REQUIRE(m->getAtomWithIdx(1)->getNumExplicitHs() == 0);

  delete m;
  smi = "F[C@]1CC(Cl)C1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  REQUIRE(!(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode)));

  delete m;
  smi = "F[C@H]1C(Cl)CC1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  REQUIRE(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));

  delete m;
  smi = "F[C@@](C)(Cl)Br";
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(!(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode)));
  REQUIRE(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  REQUIRE(!(m->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode)));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");

  delete m;
  smi = "F[C@](Br)(C)Cl";
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");
  delete m;
  smi = "F[C@](Cl)(Br)C";
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");

  delete m;
  smi = "FC(F)(F)[C@](Br)(F)C(Cl)(Cl)Cl";
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");

  delete m;
  smi = "C[C@](C=C)(F)Br";
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");

  delete m;
  smi = "CC[C@](C=C)(F)Br";
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");
  delete m;
  std::cerr << "-------------------------------" << std::endl;
  smi = "[CH2-][C@](C)(F)Br";
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");

  delete m;
  smi = "F[C@]([H])(Cl)Br";
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");

  delete m;
  smi = "[C@H](Cl)(F)Br";
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");

  delete m;
  smi = "[C@]([H])(Cl)(F)Br";
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");

  delete m;
  smi = "F[C@H](Cl)Br";
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");

  delete m;
  smi = "CC[C@H](C=C)C";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------"
                        << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");

  delete m;
  smi = "OC[C@H](C=C)C";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------"
                        << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");

  delete m;
  smi = "CC[C@H](C=C)O";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------"
                        << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");

  delete m;
  smi = "OC[C@H](C=C)O";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------"
                        << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");

  delete m;
  smi = "C[C@H]1C[C@H](C=C1)N";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------"
                        << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");
  REQUIRE(m->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");

  // a couple random molecules from the BBB data set
  delete m;
  smi = "OC[C@H]1C[C@@H](N2C=NC3=C2N=C(N)N=C3NC4CC4)C=C1";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------"
                        << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");
  REQUIRE(m->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");

  delete m;
  smi = "N[C@H]1O[C@@H](SC1)CO";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------"
                        << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");
  REQUIRE(m->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");

  delete m;
  smi = "C1(N([C@H]2O[C@H](CO)SC2)C=CC(N)=N1)=O";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------"
                        << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");
  REQUIRE(m->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");

  // this is Issue 152:
  smi = "C1[C@H](N)C[C@H](C)C=1";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------"
                        << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  delete m;
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");
  REQUIRE(m->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");

  // -----------------------------------------------
  // these are related to Issue 397:
  smi = "C(=O)[C@@H](C)N";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------"
                        << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  delete m;
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");

  smi = "C(=O)[C@@H](CO)N";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------"
                        << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  delete m;
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");

  smi = "C(O)[C@@H](C)N";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------"
                        << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  delete m;
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");
  // -----------------------------------------------

  // NOTE: This test gives correct results according to the current
  // CIP ranking procedure, but the results are actually incorrect.
  // This arises because of the inclusion of hybridization in the
  // chiral atom invariants
  // (see the note in Chirality.cpp:buildCIPInvariants())
  smi = "[H][C@@](O)(C=C)C(C)CC";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------"
                        << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  delete m;
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");

  smi = "[H][C@@](O)(C=C)C(C)CO";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------"
                        << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  delete m;
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");

  smi = "[H][C@@]12C[C@@](NC1)(OC2)[H]";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------"
                        << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  delete m;
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");
  REQUIRE(m->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");

  smi = "[H][C@@]12C[C@@](C=C1)(CC2)[H]";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------"
                        << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  delete m;
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");
  REQUIRE(m->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");

  smi = "[H][C@@]12O[C@@](CC1)(C3C2C(NC3=O)=O)[H]";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------"
                        << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  delete m;
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");
  REQUIRE(m->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");

  smi = "[H][C@@]12O[C@@](C=C1)(C3C2C(NC3=O)=O)[H]";
#ifdef VERBOSE_CANON
  BOOST_LOG(rdDebugLog) << " ----------------- ------------- ----------------"
                        << std::endl;
  BOOST_LOG(rdDebugLog) << "\t>" << smi << std::endl;
#endif
  delete m;
  m = SmilesToMol(smi);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true);
  REQUIRE(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "R");
  REQUIRE(m->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, cip);
  REQUIRE(cip == "S");
  delete m;
  // -----------------------------------------------
}

TEST_CASE("Testing double bond stereochemistry") {
  for (const bool useLegacy : {true, false}) {
    UseLegacyStereoPerceptionFixture fx(useLegacy);

    std::string smi = "F\\C=C/Cl";
    auto m = v2::SmilesParse::MolFromSmiles(smi);
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(0)->getStereo() == Bond::STEREONONE);
    if (useLegacy) {
      REQUIRE(m->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);
    } else {
      REQUIRE(m->getBondWithIdx(1)->getStereo() == Bond::STEREOCIS);
      REQUIRE(m->getBondWithIdx(1)->getStereoAtoms()[0] == 0);
      REQUIRE(m->getBondWithIdx(1)->getStereoAtoms()[1] == 3);
    }
    REQUIRE(m->getBondWithIdx(2)->getStereo() == Bond::STEREONONE);

    smi = "F/C=CCl";
    m.reset(SmilesToMol(smi));
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(1)->getStereo() == Bond::STEREONONE);

    smi = "F/C=C/Cl";
    m.reset(SmilesToMol(smi));
    REQUIRE(m);
    if (useLegacy) {
      REQUIRE(m->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
    } else {
      REQUIRE(m->getBondWithIdx(1)->getStereo() == Bond::STEREOTRANS);
      REQUIRE(m->getBondWithIdx(1)->getStereoAtoms()[0] == 0);
      REQUIRE(m->getBondWithIdx(1)->getStereoAtoms()[1] == 3);
    }

    smi = "F/C=C(/Br)Cl";
    m.reset(SmilesToMol(smi));
    REQUIRE(m);
    if (useLegacy) {
      REQUIRE(m->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
    } else {
      REQUIRE(m->getBondWithIdx(1)->getStereo() == Bond::STEREOTRANS);
      REQUIRE(m->getBondWithIdx(1)->getStereoAtoms()[0] == 0);
      REQUIRE(m->getBondWithIdx(1)->getStereoAtoms()[1] == 3);
    }

    smi = "F/C=C(/Cl)Br";
    m.reset(SmilesToMol(smi));
    REQUIRE(m);
    if (useLegacy) {
      REQUIRE(m->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);
    } else {
      REQUIRE(m->getBondWithIdx(1)->getStereo() == Bond::STEREOTRANS);
      REQUIRE(m->getBondWithIdx(1)->getStereoAtoms()[0] == 0);
      REQUIRE(m->getBondWithIdx(1)->getStereoAtoms()[1] == 3);
    }

    smi = "F/C(Br)=C/Cl";
    m.reset(SmilesToMol(smi));
    REQUIRE(m);
    if (useLegacy) {
      REQUIRE(m->getBondWithIdx(2)->getStereo() == Bond::STEREOZ);
    } else {
      REQUIRE(m->getBondWithIdx(2)->getStereo() == Bond::STEREOTRANS);
      REQUIRE(m->getBondWithIdx(2)->getStereoAtoms()[0] == 0);
      REQUIRE(m->getBondWithIdx(2)->getStereoAtoms()[1] == 4);
    }

    smi = "F/C=C(/Cl)Cl";
    m.reset(SmilesToMol(smi));
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(1)->getStereo() == Bond::STEREONONE);

    // build a molecule from scratch to test problems
    // around Issue 180. The molecule corresponds to SMILES
    // F/C=C(/Br)C
    auto m2 = std::make_unique<RWMol>();
    m2->addAtom(new Atom(9), true, true);
    m2->addAtom(new Atom(6), true, true);
    m2->addAtom(new Atom(6), true, true);
    m2->addAtom(new Atom(35), true, true);
    m2->addAtom(new Atom(6), true, true);
    m2->addBond(0, 1, Bond::SINGLE);
    m2->addBond(1, 2, Bond::DOUBLE);
    m2->addBond(2, 3, Bond::SINGLE);
    m2->addBond(2, 4, Bond::SINGLE);
    m2->getBondWithIdx(0)->setBondDir(Bond::ENDUPRIGHT);
    m2->getBondWithIdx(2)->setBondDir(Bond::ENDUPRIGHT);
    m2->getBondWithIdx(3)->setBondDir(Bond::ENDDOWNRIGHT);
    MolOps::sanitizeMol(*m2);
    MolOps::assignStereochemistry(*m2);
    if (useLegacy) {
      REQUIRE(m2->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
    } else {
      REQUIRE(m2->getBondWithIdx(1)->getStereo() == Bond::STEREOTRANS);
      REQUIRE(m2->getBondWithIdx(1)->getStereoAtoms()[0] == 0);
      REQUIRE(m2->getBondWithIdx(1)->getStereoAtoms()[1] == 3);
    }

    m2->getBondWithIdx(0)->setBondDir(Bond::ENDDOWNRIGHT);
    MolOps::assignStereochemistry(*m2, true, true);
    if (useLegacy) {
      REQUIRE(m2->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);
    } else {
      REQUIRE(m2->getBondWithIdx(1)->getStereo() == Bond::STEREOCIS);
      REQUIRE(m2->getBondWithIdx(1)->getStereoAtoms()[0] == 0);
      REQUIRE(m2->getBondWithIdx(1)->getStereoAtoms()[1] == 3);
    }

    m2.reset(new RWMol());
    m2->addAtom(new Atom(9), true, true);
    m2->addAtom(new Atom(6), true, true);
    m2->addAtom(new Atom(6), true, true);
    m2->addAtom(new Atom(35), true, true);
    m2->addAtom(new Atom(6), true, true);
    m2->addBond(1, 0, Bond::SINGLE);
    m2->addBond(1, 2, Bond::DOUBLE);
    m2->addBond(2, 3, Bond::SINGLE);
    m2->addBond(2, 4, Bond::SINGLE);
    m2->getBondWithIdx(0)->setBondDir(Bond::ENDUPRIGHT);
    m2->getBondWithIdx(2)->setBondDir(Bond::ENDUPRIGHT);
    m2->getBondWithIdx(3)->setBondDir(Bond::ENDDOWNRIGHT);
    MolOps::sanitizeMol(*m2);
    MolOps::assignStereochemistry(*m2);
    if (useLegacy) {
      REQUIRE(m2->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);
    } else {
      REQUIRE(m2->getBondWithIdx(1)->getStereo() == Bond::STEREOCIS);
      REQUIRE(m2->getBondWithIdx(1)->getStereoAtoms()[0] == 0);
      REQUIRE(m2->getBondWithIdx(1)->getStereoAtoms()[1] == 3);
    }

    m2->getBondWithIdx(0)->setBondDir(Bond::ENDDOWNRIGHT);
    MolOps::assignStereochemistry(*m2, true, true);
    if (useLegacy) {
      REQUIRE(m2->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
    } else {
      REQUIRE(m2->getBondWithIdx(1)->getStereo() == Bond::STEREOTRANS);
      REQUIRE(m2->getBondWithIdx(1)->getStereoAtoms()[0] == 0);
      REQUIRE(m2->getBondWithIdx(1)->getStereoAtoms()[1] == 3);
    }

    // ----------------------
    // test Issue 174:
    smi = "O\\N=C\\C=N/O";
    m2.reset(SmilesToMol(smi));
    REQUIRE(m2);
    if (useLegacy) {
      REQUIRE(m2->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
      REQUIRE(m2->getBondWithIdx(3)->getStereo() == Bond::STEREOZ);
    } else {
      REQUIRE(m2->getBondWithIdx(1)->getStereo() == Bond::STEREOTRANS);

      REQUIRE(m2->getBondWithIdx(1)->getStereoAtoms()[0] == 0);
      REQUIRE(m2->getBondWithIdx(1)->getStereoAtoms()[1] == 3);
      REQUIRE(m2->getBondWithIdx(3)->getStereo() == Bond::STEREOCIS);
      REQUIRE(m2->getBondWithIdx(3)->getStereoAtoms()[0] == 2);
      REQUIRE(m2->getBondWithIdx(3)->getStereoAtoms()[1] == 5);
    }
    auto refSmi = MolToSmiles(*m2, 1);
    m.reset(SmilesToMol(refSmi));
    REQUIRE(m);
    smi = MolToSmiles(*m, 1);
    REQUIRE(refSmi == smi);
  }
}

TEST_CASE("Testing Issue 183") {
  // ----------------------
  // test "unsetting" of redundant bond directions:

  RWMol *m, *m2;
  std::string smi;
  std::string refSmi;

  smi = "Cl\\C(C)=C(\\C(F)=C(/F)C)/C(C)=C(\\F)C";
  m2 = SmilesToMol(smi);
  REQUIRE(m2);
  REQUIRE(m2->getBondWithIdx(2)->getStereo() == Bond::STEREOE);
  REQUIRE(m2->getBondWithIdx(5)->getStereo() == Bond::STEREOE);
  REQUIRE(m2->getBondWithIdx(10)->getStereo() == Bond::STEREOZ);

  // m2->debugMol(std::cerr);
  refSmi = MolToSmiles(*m2, 1);
  BOOST_LOG(rdInfoLog) << "ref: " << refSmi << std::endl;
  m = SmilesToMol(refSmi);
  REQUIRE(m);
  smi = MolToSmiles(*m, 1);
  BOOST_LOG(rdInfoLog) << "smi: " << smi << std::endl;
  REQUIRE(refSmi == smi);

  int nEs = 0, nZs = 0, nDbl = 0;
  for (RWMol::BondIterator bondIt = m->beginBonds(); bondIt != m->endBonds();
       bondIt++) {
    if ((*bondIt)->getBondType() == Bond::DOUBLE) {
      nDbl++;
      if ((*bondIt)->getStereo() == Bond::STEREOE) {
        nEs++;
      } else if ((*bondIt)->getStereo() == Bond::STEREOZ) {
        nZs++;
      }
    }
  }

  REQUIRE(nDbl == 3);
  REQUIRE(nEs == 2);
  REQUIRE(nZs == 1);
  delete m;
  delete m2;
}

TEST_CASE("Testing Issue 188: bad CIP rankings") {
  ROMol *m;
  std::string smi;
  unsigned int cip1, cip2, cip3;

  smi = "OC[C@H](C=C)C";
  m = SmilesToMol(smi);
  UINT_VECT ranks;
  ranks.resize(m->getNumAtoms());
  Chirality::assignAtomCIPRanks(*m, ranks);
  REQUIRE(m);
  REQUIRE(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPRank));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPRank, cip1);
  REQUIRE(m->getAtomWithIdx(3)->hasProp(common_properties::_CIPRank));
  m->getAtomWithIdx(3)->getProp(common_properties::_CIPRank, cip2);
  REQUIRE(cip1 > cip2);
  REQUIRE(m->getAtomWithIdx(5)->hasProp(common_properties::_CIPRank));
  m->getAtomWithIdx(5)->getProp(common_properties::_CIPRank, cip3);
  REQUIRE(cip1 > cip3);
  REQUIRE(cip2 > cip3);

  delete m;
  smi = "CC(=N\\N)/C=N/N";
  m = SmilesToMol(smi);
  REQUIRE(m);
  ranks.resize(m->getNumAtoms());
  Chirality::assignAtomCIPRanks(*m, ranks);
  REQUIRE(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPRank));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPRank, cip1);
  REQUIRE(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPRank));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPRank, cip2);
  REQUIRE(cip2 > cip1);
  REQUIRE(m->getAtomWithIdx(4)->hasProp(common_properties::_CIPRank));
  m->getAtomWithIdx(4)->getProp(common_properties::_CIPRank, cip3);
  REQUIRE(cip3 > cip1);
  REQUIRE(cip2 > cip3);
  delete m;
}

TEST_CASE("Testing Issue 189: BondDirs not getting properly cleared") {
  ROMol *m;
  std::string smi, refSmi;
  int count;

  smi = "C(=S)/N=c(/n1C)scc1";
  m = SmilesToMol(smi);
  REQUIRE(m);

  REQUIRE(m->getBondWithIdx(2)->getBondType() == Bond::DOUBLE);
  REQUIRE(m->getBondWithIdx(2)->getStereo() == Bond::STEREOZ);

  count = 0;
  for (unsigned int i = 0; i < m->getNumBonds(); i++) {
    if (m->getBondWithIdx(i)->getBondDir() != Bond::NONE) {
      count++;
    }
  }
  REQUIRE(count == 2);

  refSmi = MolToSmiles(*m, 1);
  count = 0;
  for (unsigned int i = 0; i < m->getNumBonds(); i++) {
    if (m->getBondWithIdx(i)->getBondDir() != Bond::NONE) {
      count++;
    }
  }
  REQUIRE(count == 2);

  delete m;
  m = SmilesToMol(refSmi);
  REQUIRE(m);
  count = 0;
  for (unsigned int i = 0; i < m->getNumBonds(); i++) {
    if (m->getBondWithIdx(i)->getBondDir() != Bond::NONE) {
      count++;
    }
  }
  REQUIRE(count == 2);
  smi = MolToSmiles(*m, 1);
  count = 0;
  for (unsigned int i = 0; i < m->getNumBonds(); i++) {
    if (m->getBondWithIdx(i)->getBondDir() != Bond::NONE) {
      count++;
    }
  }
  REQUIRE(count == 2);
  REQUIRE(smi == refSmi);

  delete m;
}

TEST_CASE("Testing Issue 190: BondDirs incorrectly cleared") {
  ROMol *m;
  std::string smi, refSmi;
  int count;

  smi = "O\\N=C\\NC(\\C)=N/OC";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
  REQUIRE(m->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  REQUIRE(m->getBondWithIdx(5)->getBondType() == Bond::DOUBLE);
  REQUIRE(m->getBondWithIdx(5)->getStereo() == Bond::STEREOZ);
  refSmi = MolToSmiles(*m, 1);

  count = 0;
  for (unsigned int i = 0; i < m->getNumBonds(); i++) {
    if (m->getBondWithIdx(i)->getBondDir() != Bond::NONE) {
      count++;
    }
  }
  REQUIRE(count == 4);

  delete m;
  m = SmilesToMol(refSmi);
  REQUIRE(m);
  count = 0;
  for (unsigned int i = 0; i < m->getNumBonds(); i++) {
    if (m->getBondWithIdx(i)->getBondDir() != Bond::NONE) {
      count++;
    }
  }
  REQUIRE(count == 4);
  smi = MolToSmiles(*m, 1);
  count = 0;
  for (unsigned int i = 0; i < m->getNumBonds(); i++) {
    if (m->getBondWithIdx(i)->getBondDir() != Bond::NONE) {
      count++;
    }
  }
  REQUIRE(count == 4);
  REQUIRE(smi == refSmi);

  delete m;
  smi = "O\\N=C\\CC(\\C)=N/OC";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
  REQUIRE(m->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  REQUIRE(m->getBondWithIdx(5)->getBondType() == Bond::DOUBLE);
  REQUIRE(m->getBondWithIdx(5)->getStereo() == Bond::STEREOZ);
  refSmi = MolToSmiles(*m, 1);

  count = 0;
  for (unsigned int i = 0; i < m->getNumBonds(); i++) {
    if (m->getBondWithIdx(i)->getBondDir() != Bond::NONE) {
      count++;
    }
  }
  REQUIRE(count == 4);
  delete m;
  m = SmilesToMol(refSmi);
  REQUIRE(m);
  count = 0;
  for (unsigned int i = 0; i < m->getNumBonds(); i++) {
    if (m->getBondWithIdx(i)->getBondDir() != Bond::NONE) {
      count++;
    }
  }
  REQUIRE(count == 4);
  smi = MolToSmiles(*m, 1);
  count = 0;
  for (unsigned int i = 0; i < m->getNumBonds(); i++) {
    if (m->getBondWithIdx(i)->getBondDir() != Bond::NONE) {
      count++;
    }
  }
  REQUIRE(count == 4);
  REQUIRE(smi == refSmi);

  delete m;
  smi = "O\\N=C\\C(=O)C(\\C)=N/OC";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
  REQUIRE(m->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  REQUIRE(m->getBondWithIdx(6)->getBondType() == Bond::DOUBLE);
  REQUIRE(m->getBondWithIdx(6)->getStereo() == Bond::STEREOZ);
  refSmi = MolToSmiles(*m, 1);

  count = 0;
  for (unsigned int i = 0; i < m->getNumBonds(); i++) {
    if (m->getBondWithIdx(i)->getBondDir() != Bond::NONE) {
      count++;
    }
  }
  REQUIRE(count == 4);
  delete m;
  m = SmilesToMol(refSmi);
  REQUIRE(m);
  count = 0;
  for (unsigned int i = 0; i < m->getNumBonds(); i++) {
    if (m->getBondWithIdx(i)->getBondDir() != Bond::NONE) {
      count++;
    }
  }
  REQUIRE(count == 4);
  smi = MolToSmiles(*m, 1);
  count = 0;
  for (unsigned int i = 0; i < m->getNumBonds(); i++) {
    if (m->getBondWithIdx(i)->getBondDir() != Bond::NONE) {
      count++;
    }
  }
  REQUIRE(count == 4);
  REQUIRE(smi == refSmi);

  delete m;
}

TEST_CASE("Testing shortest path code") {
  {
    std::string smi = "CC(OC1C(CCCC3)C3C(CCCC2)C2C1OC(C)=O)=O";
    ROMol *m = SmilesToMol(smi);

    INT_LIST path = MolOps::getShortestPath(*m, 1, 20);
    CHECK_INVARIANT(path.size() == 7, "");
    INT_LIST_CI pi = path.begin();
    CHECK_INVARIANT((*pi) == 1, "");
    pi++;
    CHECK_INVARIANT((*pi) == 2, "");
    pi++;
    CHECK_INVARIANT((*pi) == 3, "");
    pi++;
    CHECK_INVARIANT((*pi) == 16, "");
    pi++;
    CHECK_INVARIANT((*pi) == 17, "");
    pi++;
    CHECK_INVARIANT((*pi) == 18, "");
    pi++;
    CHECK_INVARIANT((*pi) == 20, "");
    pi++;
    delete m;
  }
  {
    // issue 2219400
    std::string smi = "CC.C";
    ROMol *m = SmilesToMol(smi);

    INT_LIST path = MolOps::getShortestPath(*m, 0, 1);
    std::cerr << "path: " << path.size() << std::endl;
    CHECK_INVARIANT(path.size() == 2, "");
    INT_LIST_CI pi = path.begin();
    CHECK_INVARIANT((*pi) == 0, "");
    pi++;
    CHECK_INVARIANT((*pi) == 1, "");

    path = MolOps::getShortestPath(*m, 1, 2);
    CHECK_INVARIANT(path.size() == 0, "");

    path = MolOps::getShortestPath(*m, 0, 2);
    CHECK_INVARIANT(path.size() == 0, "");
    delete m;
  }
  // fused ring test
  {
    std::string smi = "[H]c1nc2c(C(=O)N([H])C2([H])Cl)c([H])c1Cl";
    ROMol *m = SmilesToMol(smi);

    INT_LIST path = MolOps::getShortestPath(*m, 8, 11);
    CHECK_INVARIANT(path.size() == 7, "");
    INT_LIST_CI pi = path.begin();
    CHECK_INVARIANT((*pi) == 8, "");
    pi++;
    CHECK_INVARIANT((*pi) == 7, "");
    pi++;
    CHECK_INVARIANT((*pi) == 2, "");
    pi++;
    pi++;  // two equally long routes here
    pi++;  // two equally long routes here
    CHECK_INVARIANT((*pi) == 10, "");
    pi++;
    CHECK_INVARIANT((*pi) == 11, "");
    pi++;
    delete m;
  }
}

TEST_CASE("Testing Issue 210") {
  ROMol *m, *m2;

  std::string smi = "C1CC1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 3);
  REQUIRE(m->getRingInfo()->isInitialized());

  m2 = MolOps::addHs(*m);
  REQUIRE(m2->getNumAtoms() == 9);
  REQUIRE(m2->getRingInfo()->isInitialized());

  delete m;
  delete m2;
}

TEST_CASE("Testing Issue 211") {
  ROMol *m;

  std::string smi = "P(c1ccccc1)(c1ccccc1)c1ccccc1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 19);

  const Atom *at = m->getAtomWithIdx(0);
  REQUIRE(at->getHybridization() == Atom::SP3);
  delete m;
}

TEST_CASE("Testing Issue 212") {
  ROMol *m, *m2;
  std::string smi, mb;
  smi = "C";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 1);
  auto *conf = new Conformer(1);
  m->addConformer(conf);
  conf->setAtomPos(0, RDGeom::Point3D(0, 0, 0));
  m2 = MolOps::addHs(*m, false, true);
  REQUIRE(m2->getNumAtoms() == 5);

  try {
    mb = MolToMolBlock(*m2);
  } catch (...) {
    REQUIRE(0);  //,"MolToMolBlock() failed");
  }

  delete m;
  delete m2;
}

TEST_CASE("Testing AddHs with coordinates") {
  ROMol *m, *m2;
  RDGeom::Point3D v;
  double bondLength = PeriodicTable::getTable()->getRb0(1) +
                      PeriodicTable::getTable()->getRb0(6);
  double tetDist = 2. * sin((109.471 / 2.) * M_PI / 180) * bondLength;
  double sp2Dist = 2. * sin(60. * M_PI / 180) * bondLength;

  std::string smi;

  smi = "C";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 1);
  auto *conf = new Conformer(1);
  m->addConformer(conf);
  conf->setAtomPos(0, RDGeom::Point3D(0, 0, 0));

  m2 = MolOps::addHs(*m, false, true);
  const Conformer *conf2 = &(m2->getConformer());
  REQUIRE(m2->getNumAtoms() == 5);
  REQUIRE(
      feq((conf2->getAtomPos(0) - conf2->getAtomPos(1)).length(), bondLength));
  REQUIRE(
      feq((conf2->getAtomPos(0) - conf2->getAtomPos(2)).length(), bondLength));
  REQUIRE(
      feq((conf2->getAtomPos(0) - conf2->getAtomPos(3)).length(), bondLength));
  REQUIRE(
      feq((conf2->getAtomPos(0) - conf2->getAtomPos(4)).length(), bondLength));
  REQUIRE(feq((conf2->getAtomPos(1) - conf2->getAtomPos(2)).length(), tetDist));
  REQUIRE(feq((conf2->getAtomPos(1) - conf2->getAtomPos(3)).length(), tetDist));
  REQUIRE(feq((conf2->getAtomPos(1) - conf2->getAtomPos(4)).length(), tetDist));
  REQUIRE(feq((conf2->getAtomPos(2) - conf2->getAtomPos(3)).length(), tetDist));
  REQUIRE(feq((conf2->getAtomPos(2) - conf2->getAtomPos(4)).length(), tetDist));
  REQUIRE(feq((conf2->getAtomPos(3) - conf2->getAtomPos(4)).length(), tetDist));
  delete m;
  delete m2;

  smi = "CC";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 2);
  conf = new Conformer(2);
  m->addConformer(conf);
  conf->setAtomPos(0, RDGeom::Point3D(0, 0, 0));
  conf->setAtomPos(1, RDGeom::Point3D(1.54, 0, 0));

  m2 = MolOps::addHs(*m, false, true);
  conf2 = &(m2->getConformer());
  REQUIRE(m2->getNumAtoms() == 8);
  REQUIRE(
      feq((conf2->getAtomPos(0) - conf2->getAtomPos(2)).length(), bondLength));
  REQUIRE(
      feq((conf2->getAtomPos(0) - conf2->getAtomPos(3)).length(), bondLength));
  REQUIRE(
      feq((conf2->getAtomPos(0) - conf2->getAtomPos(4)).length(), bondLength));
  REQUIRE(feq((conf2->getAtomPos(2) - conf2->getAtomPos(3)).length(), tetDist));
  REQUIRE(feq((conf2->getAtomPos(2) - conf2->getAtomPos(4)).length(), tetDist));
  REQUIRE(feq((conf2->getAtomPos(3) - conf2->getAtomPos(4)).length(), tetDist));
  REQUIRE(
      feq((conf2->getAtomPos(1) - conf2->getAtomPos(5)).length(), bondLength));
  REQUIRE(
      feq((conf2->getAtomPos(1) - conf2->getAtomPos(6)).length(), bondLength));
  REQUIRE(
      feq((conf2->getAtomPos(1) - conf2->getAtomPos(7)).length(), bondLength));
  REQUIRE(feq((conf2->getAtomPos(5) - conf2->getAtomPos(6)).length(), tetDist));
  REQUIRE(feq((conf2->getAtomPos(5) - conf2->getAtomPos(7)).length(), tetDist));
  REQUIRE(feq((conf2->getAtomPos(6) - conf2->getAtomPos(7)).length(), tetDist));

  delete m;
  delete m2;

  smi = "C=C";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 2);
  conf = new Conformer(2);
  m->addConformer(conf);
  conf->setAtomPos(0, RDGeom::Point3D(0, 0, 0));
  conf->setAtomPos(1, RDGeom::Point3D(1.3, 0, 0));

  m2 = MolOps::addHs(*m, false, true);

  conf2 = &(m2->getConformer());

  REQUIRE(m2->getNumAtoms() == 6);

  REQUIRE(
      feq((conf2->getAtomPos(0) - conf2->getAtomPos(2)).length(), bondLength));
  REQUIRE(
      feq((conf2->getAtomPos(0) - conf2->getAtomPos(3)).length(), bondLength));
  REQUIRE(
      feq((conf2->getAtomPos(1) - conf2->getAtomPos(4)).length(), bondLength));
  REQUIRE(
      feq((conf2->getAtomPos(1) - conf2->getAtomPos(5)).length(), bondLength));
  REQUIRE(feq((conf2->getAtomPos(2) - conf2->getAtomPos(3)).length(), sp2Dist));
  REQUIRE(feq((conf2->getAtomPos(4) - conf2->getAtomPos(5)).length(), sp2Dist));
  delete m;
  delete m2;

  {
    // make sure Hs are on the xy plane if conformer is 2D
    auto m = "C=C"_smiles;
    REQUIRE(m.get());
    REQUIRE(m->getNumAtoms() == 2);
    conf = new Conformer(2);
    m->addConformer(conf);
    conf->setAtomPos(0, RDGeom::Point3D(0, 0, 0));
    conf->setAtomPos(1, RDGeom::Point3D(1.3, 0, 0));
    conf->set3D(false);

    MolOps::addHs(*m, false, true);

    conf = &(m->getConformer());

    REQUIRE(m->getNumAtoms() == 6);

    REQUIRE(feq(conf->getAtomPos(2).z, 0.0));
    REQUIRE(feq(conf->getAtomPos(3).z, 0.0));
    REQUIRE(feq(conf->getAtomPos(4).z, 0.0));
    REQUIRE(feq(conf->getAtomPos(5).z, 0.0));
  }

  {
    // make sure NHs are on the same plane as the double bond
    auto m = "NC=C"_smiles;
    REQUIRE(m.get());
    REQUIRE(m->getNumAtoms() == 3);
    unsigned int nh2Idx = 0;
    unsigned int chIdx = 1;
    unsigned int ch2Idx = 2;
    conf = new Conformer(3);
    m->addConformer(conf);
    conf->setAtomPos(nh2Idx, RDGeom::Point3D(1.759236, 0.825542, 1.347849));
    conf->setAtomPos(chIdx, RDGeom::Point3D(0.817392, 0.181048, 2.180373));
    conf->setAtomPos(ch2Idx, RDGeom::Point3D(-0.070943, 0.888262, 2.875625));

    MolOps::addHs(*m, false, true);

    conf = &(m->getConformer());

    REQUIRE(m->getNumAtoms() == 8);

    std::vector<unsigned int> nh2HIdxs{3, 4};
    unsigned int chHIdx = 5;
    std::vector<unsigned int> ch2HIdxs{6, 7};
    RDGeom::Point3D nccNormal =
        (conf->getAtomPos(nh2Idx) - conf->getAtomPos(chIdx))
            .crossProduct((conf->getAtomPos(ch2Idx) - conf->getAtomPos(chIdx)));
    nccNormal.normalize();
    RDGeom::Point3D hnhNormal =
        (conf->getAtomPos(nh2HIdxs[0]) - conf->getAtomPos(nh2Idx))
            .crossProduct(conf->getAtomPos(nh2HIdxs[1]) -
                          conf->getAtomPos(nh2Idx));
    hnhNormal.normalize();
    RDGeom::Point3D hchNormal =
        (conf->getAtomPos(ch2HIdxs[0]) - conf->getAtomPos(ch2Idx))
            .crossProduct(conf->getAtomPos(nh2HIdxs[1]) -
                          conf->getAtomPos(ch2Idx));
    hchNormal.normalize();
    RDGeom::Point3D hcnNormal =
        (conf->getAtomPos(chHIdx) - conf->getAtomPos(chIdx))
            .crossProduct((conf->getAtomPos(nh2Idx) - conf->getAtomPos(chIdx)));
    hcnNormal.normalize();
    REQUIRE(feq(fabs(nccNormal.dotProduct(hnhNormal)), 1.0));
    REQUIRE(feq(fabs(nccNormal.dotProduct(hchNormal)), 1.0));
    REQUIRE(feq(fabs(nccNormal.dotProduct(hcnNormal)), 1.0));
  }

  smi = "C#C";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 2);
  conf = new Conformer(2);
  m->addConformer(conf);
  conf->setAtomPos(0, RDGeom::Point3D(0, 0, 0));
  conf->setAtomPos(1, RDGeom::Point3D(1.2, 0, 0));

  m2 = MolOps::addHs(*m, false, true);
  conf2 = &(m2->getConformer());
  REQUIRE(m2->getNumAtoms() == 4);
  REQUIRE(
      feq((conf2->getAtomPos(0) - conf2->getAtomPos(2)).length(), bondLength));
  REQUIRE(
      feq((conf2->getAtomPos(1) - conf2->getAtomPos(3)).length(), bondLength));
  REQUIRE(feq((conf2->getAtomPos(0) - conf2->getAtomPos(3)).length(),
              bondLength + 1.2));
  REQUIRE(feq((conf2->getAtomPos(1) - conf2->getAtomPos(2)).length(),
              bondLength + 1.2));

  delete m;
  delete m2;
}

TEST_CASE("Sanitization special cases") {
  ROMol *m;
  std::string smi, pathName;

  smi = "CN(=O)=O";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 4);
  REQUIRE(m->getAtomWithIdx(1)->getFormalCharge() == 1);
  delete m;

  smi = "C[N+](=O)[O-]";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 4);
  REQUIRE(m->getAtomWithIdx(1)->getFormalCharge() == 1);
  delete m;

  smi = "Cl(=O)(=O)(=O)[O-]";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 5);
  REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == 3);
  delete m;

  smi = "Cl(=O)(=O)(=O)O";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 5);
  REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == 3);
  delete m;

  smi = "Br(=O)(=O)(=O)[O-]";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 5);
  REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == 3);
  delete m;

  smi = "Br(=O)(=O)(=O)O";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 5);
  REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == 3);
  delete m;

  smi = "I(=O)(=O)(=O)[O-]";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 5);
  REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == 3);
  delete m;

  smi = "I(=O)(=O)(=O)O";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 5);
  REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == 3);
  delete m;

  pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/test_data/";
  m = MolFileToMol(pathName + "perchlorate1.mol");
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 51);
  REQUIRE(m->getAtomWithIdx(7)->getFormalCharge() == 3);
  delete m;

  smi = "CN=N#N";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 4);
  REQUIRE(m->getAtomWithIdx(1)->getFormalCharge() == 0);
  REQUIRE(m->getAtomWithIdx(2)->getFormalCharge() == 1);
  REQUIRE(m->getAtomWithIdx(3)->getFormalCharge() == -1);
  REQUIRE(m->getBondBetweenAtoms(2, 3)->getBondType() == Bond::DOUBLE);
  delete m;

  smi = "N#N=NC";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 4);
  REQUIRE(m->getAtomWithIdx(2)->getFormalCharge() == 0);
  REQUIRE(m->getAtomWithIdx(1)->getFormalCharge() == 1);
  REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == -1);
  REQUIRE(m->getBondBetweenAtoms(0, 1)->getBondType() == Bond::DOUBLE);
  delete m;

  smi = "N#N";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 2);
  REQUIRE(m->getAtomWithIdx(1)->getFormalCharge() == 0);
  REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == 0);
  REQUIRE(m->getBondBetweenAtoms(0, 1)->getBondType() == Bond::TRIPLE);
  delete m;

  smi = "Cl(=O)(=O)O";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 4);
  REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == 2);
  delete m;

  smi = "Cl(=O)(=O)[O-]";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 4);
  REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == 2);
  delete m;

  smi = "Br(=O)(=O)O";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 4);
  REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == 2);
  delete m;

  smi = "Br(=O)(=O)[O-]";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 4);
  REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == 2);
  delete m;

  smi = "I(=O)(=O)O";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 4);
  REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == 2);
  delete m;

  smi = "I(=O)(=O)[O-]";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 4);
  REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == 2);
  delete m;

  smi = "Cl(=O)O";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 3);
  REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == 1);
  delete m;

  smi = "Cl(=O)[O-]";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 3);
  REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == 1);
  delete m;

  smi = "Br(=O)O";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 3);
  REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == 1);
  delete m;

  smi = "Br(=O)[O-]";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 3);
  REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == 1);
  delete m;

  smi = "I(=O)O";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 3);
  REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == 1);
  delete m;

  smi = "I(=O)[O-]";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 3);
  REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == 1);
  delete m;

  smi = "I(O)(O)(O)(O)O";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 6);
  REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == 0);
  delete m;

  smi = "I(O)(O)O";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 4);
  REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == 0);
  delete m;

  smi = "I(=O)(O)(O)(O)";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 5);
  REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == 1);
  delete m;

  smi = "CC(=O)O[IH2](O)OC(C)=O";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 10);
  REQUIRE(m->getAtomWithIdx(4)->getFormalCharge() == 0);
  delete m;
}

TEST_CASE("Testing Add Confomers") {
  std::string smi = "CC";
  ROMol *m = SmilesToMol(smi);
  unsigned int i;
  for (i = 0; i < 5; i++) {
    auto *conf = new Conformer(2);
    conf->setAtomPos(0, RDGeom::Point3D(0.0, 0.0, 0.0));
    conf->setAtomPos(1, RDGeom::Point3D(1.5, 0.0, 0.0));
    m->addConformer(conf, true);
  }
  CHECK_INVARIANT(m->getNumConformers() == 5, "");

  ROMol *m2 = MolOps::addHs(*m, false, true);
  CHECK_INVARIANT(m2->getNumConformers() == 5, "");
  // const ROMol::CONF_SPTR_LIST &confs = m2->getConformers();
  ROMol::ConstConformerIterator ci;
  i = 0;
  for (ci = m2->beginConformers(); ci != m2->endConformers(); ci++) {
    CHECK_INVARIANT((*ci)->getNumAtoms() == 8, "");
    CHECK_INVARIANT((*ci)->getId() == i, "");
    const ROMol *mn = &((*ci)->getOwningMol());
    CHECK_INVARIANT(mn->getNumAtoms() == 8, "");
    i++;
  }
  // std::cout << m2->getNumAtoms() << " " << m2->getNumConformers() << "\n";
  delete m;
  delete m2;
}

TEST_CASE("testIssue252") {
  // lets check if we can sanitize C60
  std::string smi =
      "C12=C3C4=C5C6=C1C7=C8C9=C1C%10=C%11C(=C29)C3=C2C3=C4C4=C5C5=C9C6=C7C6="
      "C7C8=C1C1=C8C%10=C%10C%11=C2C2=C3C3=C4C4=C5C5=C%11C%12=C(C6=C95)C7=C1C1="
      "C%12C5=C%11C4=C3C3=C5C(=C81)C%10=C23";
  ROMol *mol = SmilesToMol(smi);
  for (ROMol::BondIterator it = mol->beginBonds(); it != mol->endBonds();
       it++) {
    REQUIRE((*it)->getIsAromatic());
    REQUIRE((*it)->getBondType() == Bond::AROMATIC);
  }
  std::string asmi = MolToSmiles(*mol);
  // check if we can do it in the aromatic form
  ROMol *nmol = SmilesToMol(asmi);
  for (ROMol::BondIterator it = nmol->beginBonds(); it != nmol->endBonds();
       it++) {
    REQUIRE((*it)->getIsAromatic());
    REQUIRE((*it)->getBondType() == Bond::AROMATIC);
  }

  std::string nsmi = MolToSmiles(*nmol);
  delete mol;
  delete nmol;

  // This is a check for Issue253
  CHECK_INVARIANT(asmi == nsmi, "");
}

TEST_CASE("Issue 276") {
  std::string smi = "CP1(C)=CC=CN=C1C";
  ROMol *mol = SmilesToMol(smi);
  REQUIRE(mol);
  // as of this writing, I'm not 100% sure what the right answer is here,
  // but the hybridization definitely should *not* be SP2:
  REQUIRE(mol->getAtomWithIdx(1)->getHybridization() > Atom::SP2);
  delete mol;
}

TEST_CASE("testHsAndAromaticity: Additional Aromaticity Cases") {
  std::string smi;
  ROMol *mol;

  smi = "[CH]1-[CH]-[CH]-[CH]-[CH]-[CH]-1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  // std::cerr << mol->getAtomWithIdx(0)->getHybridization() << std::endl;
  REQUIRE(mol->getAtomWithIdx(0)->getHybridization() == Atom::SP3);
  REQUIRE(mol->getAtomWithIdx(0)->getValence(Atom::ValenceType::IMPLICIT) == 0);
  REQUIRE(mol->getAtomWithIdx(0)->getNumImplicitHs() == 0);
  REQUIRE(mol->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
  REQUIRE(!mol->getAtomWithIdx(0)->getIsAromatic());
  REQUIRE(!mol->getBondBetweenAtoms(0, 1)->getIsAromatic());
  delete mol;

  smi = "C1=CC(=C)C(=C)C=C1";
  mol = SmilesToMol(smi);
  REQUIRE(mol);
  REQUIRE(mol->getAtomWithIdx(0)->getHybridization() == Atom::SP2);
  REQUIRE(mol->getAtomWithIdx(2)->getHybridization() == Atom::SP2);
  REQUIRE(mol->getAtomWithIdx(0)->getIsAromatic());
  REQUIRE(mol->getBondBetweenAtoms(0, 1)->getIsAromatic());

  delete mol;
}

TEST_CASE(
    "Testing sf.net issue 1694023 (bad chiral smiles after removing Hs)") {
  ROMol *m;

  std::string smi;

  smi = "[C@@]([H])(F)(Cl)Br";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 4);
  REQUIRE(m->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);

  smi = "[C@@](F)([H])(Cl)Br";
  delete m;
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 4);
  REQUIRE(m->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);

  smi = "[C@@](F)(Cl)([H])Br";
  delete m;
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 4);
  REQUIRE(m->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);

  smi = "[C@@](F)(Cl)(Br)[H]";
  delete m;
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 4);
  REQUIRE(m->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);

  smi = "[H][C@@](F)(Cl)Br";
  delete m;
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 4);
  REQUIRE(m->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);

  smi = "F[C@@]([H])(Cl)Br";
  delete m;
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 4);
  REQUIRE(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);

  smi = "F[C@@](Cl)([H])Br";
  delete m;
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 4);
  REQUIRE(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);

  smi = "F[C@@](Cl)(Br)[H]";
  delete m;
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 4);
  REQUIRE(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);

  smi = "C1CO[C@@H]1Cl";
  delete m;
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 5);
  REQUIRE(m->getAtomWithIdx(3)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);

  smi = "C1CO[C@]1([H])Cl";
  delete m;
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 5);
  REQUIRE(m->getAtomWithIdx(3)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);

  smi = "C1CO[C@@]1(Cl)[H]";
  delete m;
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 5);
  REQUIRE(m->getAtomWithIdx(3)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);

  delete m;
}

TEST_CASE(
    "Testing sf.net issue 1719053 (Ring stereochemistry incorrectly removed)") {
  ROMol *m;

  std::string smi;

  smi = "C[C@@H]1CCCCC1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);

  delete m;
  smi = "C[C@@H]1CC[C@@H](C)CC1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  REQUIRE(m->getAtomWithIdx(4)->getChiralTag() != Atom::CHI_UNSPECIFIED);

  delete m;
  smi = "C[C@@H]1C(C)CCCC1C";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);

  delete m;
  smi = "C[C@@H]1[C@H](C)CCC[C@H]1C";
  m = SmilesToMol(smi);
  REQUIRE(m);
  // this is a truly symmetric case, so the stereochem should be removed:
  REQUIRE(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  REQUIRE(m->getAtomWithIdx(2)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  REQUIRE(m->getAtomWithIdx(7)->getChiralTag() != Atom::CHI_UNSPECIFIED);

  delete m;
  smi = "C[C@@H]1C=C[C@@H](C)C=C1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  REQUIRE(m->getAtomWithIdx(4)->getChiralTag() != Atom::CHI_UNSPECIFIED);

  delete m;
  smi = "C[N@@]1C=C[C@@H](C)C=C1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  REQUIRE(m->getAtomWithIdx(4)->getChiralTag() == Atom::CHI_UNSPECIFIED);

  // N in rings that aren't 3 rings is not chiral
  delete m;
  smi = "C[N@@]1CC[C@@H](C)CC1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  REQUIRE(m->getAtomWithIdx(4)->getChiralTag() == Atom::CHI_UNSPECIFIED);

  delete m;
}

TEST_CASE("Testing sf.net issue 1811276 (kekulization failing)") {
  ROMol *m;

  std::string smi;

  smi = "[O-]N1C=C[N+](=O)C=C1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  smi = MolToSmiles(*m);
  REQUIRE(smi == "O=[n+]1ccn([O-])cc1");
  delete m;

  smi = "o1ccc(=O)cc1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  smi = MolToSmiles(*m);
  REQUIRE(smi == "O=c1ccocc1");
  delete m;

  smi = "O=[n+]1ccocc1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  smi = MolToSmiles(*m);
  REQUIRE(smi == "O=[n+]1ccocc1");
  delete m;

  smi = "O=[n+]1ccn([O-])cc1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  smi = MolToSmiles(*m);
  REQUIRE(smi == "O=[n+]1ccn([O-])cc1");
  delete m;

  smi = "O=n1ccccc1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  smi = MolToSmiles(*m);
  REQUIRE(smi == "[O-][n+]1ccccc1");

  delete m;
}

TEST_CASE("Testing sf.net issue 1836576 (sanitization crash)") {
  RWMol *m;

  std::string smi;
  bool ok;

  // the original form of the test runs foul of the rules for explicit
  // valence on B:
  smi =
      "[BH]123[BH]45[BH]167[BH]289[BH]312[BH]838[BH]966[Co]74479%10%11%12[CH]"
      "633[BH]811[CH]345[BH]21[BH]1234[BH]75[BH]911[BH]226[BH]%1011[BH]227[BH]"
      "633[BH]44[BH]322[CH]%1145[CH]%12271";
  m = SmilesToMol(smi, false, false);
  REQUIRE(m);

  unsigned int opThatFailed;
  ok = false;
  try {
    MolOps::sanitizeMol(*m, opThatFailed);
  } catch (MolSanitizeException &) {
    ok = true;
  }
  REQUIRE(ok);
  REQUIRE(opThatFailed == MolOps::SANITIZE_PROPERTIES);
  delete m;
}

TEST_CASE("Testing impact of removeHs on chirality") {
  ROMol *m, *m2;

  std::string smi, code;

  smi = "F[C@]([H])(Cl)Br";
  m = SmilesToMol(smi, false, false);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true, true);
  REQUIRE(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, code);
  REQUIRE(code == "R");
  m2 = MolOps::removeHs(*m);
  REQUIRE(m2);
  MolOps::assignStereochemistry(*m2, true, true);
  REQUIRE(m2->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m2->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, code);
  REQUIRE(code == "R");
  delete m;
  delete m2;

  smi = "F[C@H](Cl)Br";
  m = SmilesToMol(smi, false, false);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true, true);
  REQUIRE(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, code);
  REQUIRE(code == "R");
  m2 = MolOps::removeHs(*m);
  REQUIRE(m2);
  MolOps::assignStereochemistry(*m2, true, true);
  REQUIRE(m2->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m2->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, code);
  REQUIRE(code == "R");
  delete m;
  delete m2;

  smi = "[C@]([H])(Cl)(F)Br";
  m = SmilesToMol(smi, false, false);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true, true);
  REQUIRE(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, code);
  REQUIRE(code == "R");
  m2 = MolOps::removeHs(*m);
  REQUIRE(m2);
  MolOps::assignStereochemistry(*m2, true, true);
  REQUIRE(m2->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m2->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, code);
  REQUIRE(code == "R");
  delete m;
  delete m2;

  smi = "[C@H](Cl)(F)Br";
  m = SmilesToMol(smi, false, false);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true, true);
  REQUIRE(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, code);
  REQUIRE(code == "R");
  m2 = MolOps::removeHs(*m);
  REQUIRE(m2);
  MolOps::assignStereochemistry(*m2, true, true);
  REQUIRE(m2->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m2->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, code);
  REQUIRE(code == "R");
  delete m;
  delete m2;

  smi = "[H]1.F[C@]1(Cl)Br";
  m = SmilesToMol(smi, false, false);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true, true);
  REQUIRE(m->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, code);
  REQUIRE(code == "R");
  m2 = MolOps::removeHs(*m);
  REQUIRE(m2);
  MolOps::assignStereochemistry(*m2, true, true);
  REQUIRE(m2->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m2->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, code);
  REQUIRE(code == "R");
  delete m;
  delete m2;

  smi = "F[C@]1(Cl)Br.[H]1";
  m = SmilesToMol(smi, false, false);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true, true);
  REQUIRE(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, code);
  REQUIRE(code == "R");
  m2 = MolOps::removeHs(*m);
  REQUIRE(m2);
  MolOps::assignStereochemistry(*m2, true, true);
  REQUIRE(m2->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m2->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, code);
  REQUIRE(code == "R");
  delete m;
  delete m2;

  smi = "[H]1.[C@]1(Cl)(F)Br";
  m = SmilesToMol(smi, false, false);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true, true);
  REQUIRE(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, code);
  REQUIRE(code == "R");
  m2 = MolOps::removeHs(*m);
  REQUIRE(m2);
  MolOps::assignStereochemistry(*m2, true, true);
  REQUIRE(m2->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m2->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, code);
  REQUIRE(code == "R");
  delete m;
  delete m2;

  smi = "[C@]1(Cl)(F)Br.[H]1";
  m = SmilesToMol(smi, false, false);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true, true);
  REQUIRE(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, code);
  REQUIRE(code == "R");
  m2 = MolOps::removeHs(*m);
  REQUIRE(m2);
  MolOps::assignStereochemistry(*m2, true, true);
  REQUIRE(m2->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m2->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, code);
  REQUIRE(code == "R");
  delete m;
  delete m2;

  smi = "Cl1.F2.Br3.[C@H]123";
  m = SmilesToMol(smi, false, false);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true, true);
  REQUIRE(m->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, code);
  REQUIRE(code == "R");
  delete m;

  smi = "[C@H]123.Cl1.F2.Br3";
  m = SmilesToMol(smi, false, false);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true, true);
  REQUIRE(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, code);
  REQUIRE(code == "R");
  delete m;

  smi = "F2.Cl1.Br3.[C@H]123";
  m = SmilesToMol(smi, false, false);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true, true);
  REQUIRE(m->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, code);
  REQUIRE(code == "R");
  delete m;

  smi = "Cl2.F1.Br3.[C@H]213";
  m = SmilesToMol(smi, false, false);
  REQUIRE(m);
  MolOps::assignStereochemistry(*m, true, true);
  REQUIRE(m->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, code);
  REQUIRE(code == "R");
  delete m;
}

TEST_CASE("Testing SFIssue1894348 (impact of removeHs on bond stereo atoms)") {
  RWMol *m, *m2;

  std::string smi;

  smi = "Cl/C([H])=C/Cl";
  m = SmilesToMol(smi, false, false);
  REQUIRE(m);
  MolOps::sanitizeMol(*m);
  MolOps::assignStereochemistry(*m);
  REQUIRE(m->getBondWithIdx(2)->getStereoAtoms().size() == 2);
  REQUIRE(m->getBondWithIdx(2)->getStereoAtoms()[0] == 0);
  REQUIRE(m->getBondWithIdx(2)->getStereoAtoms()[1] == 4);
  // we remove an H attached to a stereo bond
  m2 = static_cast<RWMol *>(MolOps::removeHs(static_cast<const ROMol &>(*m)));
  REQUIRE(m->getBondWithIdx(2)->getStereoAtoms().size() == 2);
  REQUIRE(m->getBondWithIdx(2)->getStereoAtoms()[0] == 0);
  REQUIRE(m->getBondWithIdx(2)->getStereoAtoms()[1] == 4);
  // at first the stereoatoms are gone:
  REQUIRE(m2->getBondWithIdx(2)->getStereoAtoms().size() == 0);
  // but they can be re-perceived:
  MolOps::assignStereochemistry(*m2, true, true);
  REQUIRE(m2->getBondWithIdx(1)->getStereoAtoms().size() == 2);
  REQUIRE(m2->getBondWithIdx(1)->getStereoAtoms()[0] == 0);
  REQUIRE(m2->getBondWithIdx(1)->getStereoAtoms()[1] == 3);

  delete m;
  delete m2;

  smi = "Cl/C([H])=C/Cl";
  m = SmilesToMol(smi, false, false);
  REQUIRE(m);
  MolOps::sanitizeMol(*m);
  REQUIRE(m->getBondWithIdx(2)->getStereoAtoms().size() == 0);
  m2 = static_cast<RWMol *>(MolOps::removeHs(static_cast<const ROMol &>(*m)));
  // if we don't assign stereocodes in the original we shouldn't have them here:
  REQUIRE(m2->getBondWithIdx(1)->getStereoAtoms().size() == 0);
  delete m;
  delete m2;
}

TEST_CASE("Testing some aromaticity edge cases") {
  RWMol *m;

  std::string smi;

  // ------
  // this was sf.net bug 1934360
  smi = "C1=C=NC=N1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(!m->getAtomWithIdx(0)->getIsAromatic());
  REQUIRE(!m->getBondWithIdx(0)->getIsAromatic());
  delete m;

  smi = "C1=CNC=N1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getAtomWithIdx(0)->getIsAromatic());
  REQUIRE(m->getBondWithIdx(0)->getIsAromatic());
  delete m;

  smi = "[C+]1=CNC=N1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(!m->getAtomWithIdx(1)->getIsAromatic());
  REQUIRE(!m->getBondWithIdx(1)->getIsAromatic());
  delete m;

  // ------
  // this was sf.net bug 1940646
  smi = "C1#CC=C1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(!m->getAtomWithIdx(0)->getIsAromatic());
  REQUIRE(!m->getBondWithIdx(0)->getIsAromatic());
  delete m;
  smi = "C1#CC=CC=C1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getAtomWithIdx(0)->getIsAromatic());
  REQUIRE(m->getBondWithIdx(0)->getIsAromatic());
  delete m;

  // ------
  // this was sf.net bug 2091839

  smi = "c1cccc[c]1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getAtomWithIdx(0)->getIsAromatic());
  REQUIRE(m->getBondWithIdx(0)->getIsAromatic());
  delete m;

  smi = "C1=CC=CC=[C]1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getAtomWithIdx(0)->getIsAromatic());
  REQUIRE(m->getBondWithIdx(0)->getIsAromatic());
  delete m;

  smi = "c1cccc[n+]1";  // disqualified because N has a radical
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(!m->getAtomWithIdx(0)->getIsAromatic());
  REQUIRE(!m->getBondWithIdx(0)->getIsAromatic());
  delete m;

  smi = "[N]1C=CC=C1";  // disqualified because N has a radical
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(!m->getAtomWithIdx(0)->getIsAromatic());
  REQUIRE(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
  REQUIRE(!m->getBondWithIdx(0)->getIsAromatic());
  delete m;

  smi = "[n]1ccccc1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getAtomWithIdx(0)->getIsAromatic());
  REQUIRE(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 0);
  REQUIRE(m->getBondWithIdx(0)->getIsAromatic());
  delete m;

  smi = "[H]n1cccc1";
  m = SmilesToMol(smi, 0, 0);
  REQUIRE(m);
  MolOps::sanitizeMol(*m);
  REQUIRE(m->getAtomWithIdx(1)->getIsAromatic());
  REQUIRE(m->getAtomWithIdx(1)->getNumRadicalElectrons() == 0);
  REQUIRE(!m->getAtomWithIdx(0)->getIsAromatic());
  REQUIRE(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 0);
  delete m;

  smi = "[H]";
  m = SmilesToMol(smi, 0, 0);
  REQUIRE(m);
  MolOps::sanitizeMol(*m);
  REQUIRE(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
  delete m;

  // ------
  // this was sf.net bug 2787221.
  smi = "O=C1C(=O)C=C1";
  m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getAtomWithIdx(1)->getIsAromatic());
  REQUIRE(m->getBondBetweenAtoms(1, 2)->getIsAromatic());
  delete m;
}

TEST_CASE("Testing sf.net issue 1942657") {
  RWMol *m;

  std::string smi;

  smi = "C[C](C)(C)(C)C";
  try {
    m = SmilesToMol(smi);
  } catch (MolSanitizeException &) {
    m = nullptr;
  }
  REQUIRE(!m);

  smi = "C[CH](C)(C)C";
  try {
    m = SmilesToMol(smi);
  } catch (MolSanitizeException &) {
    m = nullptr;
  }
  REQUIRE(!m);

  smi = "C[C](=C)(C)C";
  try {
    m = SmilesToMol(smi);
  } catch (MolSanitizeException &) {
    m = nullptr;
  }
  REQUIRE(!m);

  smi = "C[Si](=C)(=C)=C";
  try {
    m = SmilesToMol(smi);
  } catch (MolSanitizeException &) {
    m = nullptr;
  }
  REQUIRE(!m);
}

TEST_CASE("Testing sf.net issue 198608") {
  RWMol *m;

  std::string smi;

  smi = "C1CC1CC1CC1";
  m = SmilesToMol(smi);
  REQUIRE(m->getRingInfo()->minAtomRingSize(0) == 3);
  REQUIRE(m->getRingInfo()->minAtomRingSize(3) == 0);
  REQUIRE(m->getRingInfo()->minBondRingSize(0) == 3);
  REQUIRE(m->getRingInfo()->minBondRingSize(3) == 0);
  delete m;
}

TEST_CASE("Testing hybridization assignment") {
  {
    RWMol *m;
    std::string smi = "CCC";
    m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getHybridization() == Atom::SP3);
    REQUIRE(m->getAtomWithIdx(0)->getHybridization() == Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi = "CNC";
    m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getHybridization() == Atom::SP3);
    REQUIRE(m->getAtomWithIdx(0)->getHybridization() == Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi = "COC";
    m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getHybridization() == Atom::SP3);
    REQUIRE(m->getAtomWithIdx(0)->getHybridization() == Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi = "C[C-2]C";
    m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getHybridization() == Atom::SP3);
    REQUIRE(m->getAtomWithIdx(0)->getHybridization() == Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi = "C[CH-]C";
    m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getHybridization() == Atom::SP3);
    REQUIRE(m->getAtomWithIdx(0)->getHybridization() == Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi = "C[CH]C";
    m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getHybridization() == Atom::SP3);
    REQUIRE(m->getAtomWithIdx(0)->getHybridization() == Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi = "C[C]C";
    m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getHybridization() == Atom::SP3);
    REQUIRE(m->getAtomWithIdx(0)->getHybridization() == Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi = "C[C-]C";
    m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getHybridization() == Atom::SP3);
    REQUIRE(m->getAtomWithIdx(0)->getHybridization() == Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi = "C[CH+]C";
    m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getHybridization() == Atom::SP2);
    REQUIRE(m->getAtomWithIdx(0)->getHybridization() == Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi = "CC=C";
    m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getHybridization() == Atom::SP2);
    REQUIRE(m->getAtomWithIdx(0)->getHybridization() == Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi = "CN=C";
    m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getHybridization() == Atom::SP2);
    REQUIRE(m->getAtomWithIdx(0)->getHybridization() == Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi = "C[C-]=C";
    m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getHybridization() == Atom::SP2);
    REQUIRE(m->getAtomWithIdx(0)->getHybridization() == Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi = "C[C]=C";
    m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getHybridization() == Atom::SP2);
    REQUIRE(m->getAtomWithIdx(0)->getHybridization() == Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi = "C[N+]=C";
    m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getHybridization() == Atom::SP2);
    REQUIRE(m->getAtomWithIdx(0)->getHybridization() == Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi = "C#C";
    m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getHybridization() == Atom::SP);
    delete m;
  }

  {
    RWMol *m;
    std::string smi = "C#[C-]";
    m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getHybridization() == Atom::SP);
    delete m;
  }

  {
    RWMol *m;
    std::string smi = "C#[C]";
    m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getHybridization() == Atom::SP);
    delete m;
  }

  {
    RWMol *m;
    std::string smi = "C[O]";
    m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getHybridization() == Atom::SP3);
    delete m;
  }

  {
    RWMol *m;
    std::string smi = "C[N-]";
    m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getHybridization() == Atom::SP3);
    delete m;
  }
}

TEST_CASE("Testing sf.net issue 2196817: handling of aromatic dummies") {
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName + "dummyArom.mol");
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(0)->getAtomicNum() == 0);
    REQUIRE(m->getAtomWithIdx(0)->getIsAromatic() == true);

    MolOps::Kekulize(*m);
    REQUIRE(m->getBondBetweenAtoms(0, 1)->getBondType() == Bond::SINGLE);
    REQUIRE(m->getBondBetweenAtoms(0, 4)->getBondType() == Bond::SINGLE);

    delete m;
  }

  {
    std::string smi = "*1cncc1";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(0)->getAtomicNum() == 0);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic() == true);
    delete m;
  }

  {
    std::string smi = "*1C=NC=C1";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(0)->getAtomicNum() == 0);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic() == true);
    delete m;
  }

  {
    // case where all must be ignored:
    std::string smi = "c1*ccc1-c1*ccc1-c1*ccc1";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    delete m;
  }

  {
    std::string smi = "c1*[nH]*c1";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    smi = MolToSmiles(*m);
    REQUIRE(smi == "*1cc*[nH]1");
    delete m;
    smi = "c1***c1";
    m = SmilesToMol(smi);
    REQUIRE(m);
    smi = MolToSmiles(*m);
    REQUIRE(smi == "*1:*cc*:1");
    delete m;
    smi = "c:1:*:*:*:*1";
    m = SmilesToMol(smi);
    REQUIRE(m);
    smi = MolToSmiles(*m);
    REQUIRE(smi == "*1:*:*c*:1");

    delete m;
    // we don't kekulize rings that are all dummies, this was github #1478
    smi = "*:1:*:*:*:*:1";
    m = SmilesToMol(smi);
    REQUIRE(m);
    smi = MolToSmiles(*m);
    REQUIRE(smi == "*1:*:*:*:*:1");
    delete m;
  }

  {
    std::string smi = "c1*[nH]cc1-c1*[nH]cc1-c1*ccc1";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    delete m;
    smi = "c1*[nH]cc1-c1*ccc1-c1*[nH]cc1";
    m = SmilesToMol(smi);
    REQUIRE(m);
    delete m;
    smi = "c1*ccc1-c1*[nH]cc1-c1*[nH1]cc1";
    m = SmilesToMol(smi);
    REQUIRE(m);
    delete m;
  }

  {
    std::string smi = "c1*[nH]cc1-c1*[nH]cc1-c1*[nH]cc1";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    delete m;
  }

  {
    std::string smi = "c1ccc(C2CC(n4cc*c4=C2))cc1";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getBondBetweenAtoms(0, 1)->getIsAromatic());
    REQUIRE(m->getBondBetweenAtoms(0, 14)->getIsAromatic());
    REQUIRE(m->getBondBetweenAtoms(0, 1)->getBondType() == Bond::AROMATIC);
    REQUIRE(m->getBondBetweenAtoms(0, 14)->getBondType() == Bond::AROMATIC);
    MolOps::Kekulize(*m);
    REQUIRE(!m->getBondBetweenAtoms(0, 1)->getIsAromatic());
    REQUIRE(!m->getBondBetweenAtoms(0, 14)->getIsAromatic());
    REQUIRE((m->getBondBetweenAtoms(0, 1)->getBondType() == Bond::DOUBLE ||
             m->getBondBetweenAtoms(0, 14)->getBondType() == Bond::DOUBLE));
    MolOps::setAromaticity(*m);
    REQUIRE(m->getBondBetweenAtoms(0, 1)->getIsAromatic());
    REQUIRE(m->getBondBetweenAtoms(0, 14)->getIsAromatic());
    REQUIRE(m->getBondBetweenAtoms(0, 1)->getBondType() == Bond::AROMATIC);
    REQUIRE(m->getBondBetweenAtoms(0, 14)->getBondType() == Bond::AROMATIC);
    delete m;
  }

  {
    for (auto smi : {"c12ccccc1**CC2", "c12ccccc1C**C2", "*12ccccc1CCCC2",
                     "*12ccccc1***C2"}) {
      std::unique_ptr<RWMol> m(SmilesToMol(smi));
      REQUIRE(m);
      for (size_t i = 0; i < 6; ++i) {
        size_t j = i < 5 ? i + 1 : 0;
        REQUIRE(m->getBondBetweenAtoms(i, j)->getIsAromatic());
        REQUIRE(m->getBondBetweenAtoms(i, j)->getBondType() == Bond::AROMATIC);
      }
      for (size_t i = 5; i < 10; ++i) {
        size_t j = i < 9 ? i + 1 : 0;
        REQUIRE(!m->getBondBetweenAtoms(i, j)->getIsAromatic());
        REQUIRE(m->getBondBetweenAtoms(i, j)->getBondType() == Bond::SINGLE);
      }
      MolOps::Kekulize(*m);
      for (size_t i = 0; i < 6; ++i) {
        size_t j = i < 5 ? i + 1 : 0;
        REQUIRE(!m->getBondBetweenAtoms(i, j)->getIsAromatic());
        REQUIRE((m->getBondBetweenAtoms(i, j)->getBondType() == Bond::SINGLE ||
                 m->getBondBetweenAtoms(i, j)->getBondType() == Bond::DOUBLE));
      }
      for (size_t i = 5; i < 10; ++i) {
        size_t j = i < 9 ? i + 1 : 0;
        REQUIRE(!m->getBondBetweenAtoms(i, j)->getIsAromatic());
        REQUIRE(m->getBondBetweenAtoms(i, j)->getBondType() == Bond::SINGLE);
      }
    }
  }

  {
    for (auto smi : {"c12ccccc1****2", "c12ccccc1***2", "C12=CC=CC=C1****2",
                     "c1ccccc1N1****1", "c1ccccc1C1****1", "c1ccccc1*1****1",
                     "c1ccccc1*1*=**=*1"}) {
      std::unique_ptr<RWMol> m(SmilesToMol(smi));
      REQUIRE(m);
      for (unsigned int i = 6; i < m->getNumAtoms(); ++i) {
        REQUIRE(!m->getAtomWithIdx(i)->getIsAromatic());
      }
    }
  }

  {
    for (auto smi : {"C1CC*2ccccc21", "C1C**2ccccc21"}) {
      std::unique_ptr<RWMol> m(SmilesToMol(smi));
      REQUIRE(m);
      for (unsigned int i = 0; i < 3; ++i) {
        REQUIRE(!m->getAtomWithIdx(i)->getIsAromatic());
      }
      for (unsigned int i = 3; i < m->getNumAtoms(); ++i) {
        REQUIRE(m->getAtomWithIdx(i)->getIsAromatic());
      }
    }
  }

  {
    auto m = "*12ccccc1CCC2"_smiles;
    for (unsigned int i = 0; i < 6; ++i) {
      REQUIRE(m->getAtomWithIdx(i)->getIsAromatic());
    }
    for (unsigned int i = 6; i < m->getNumAtoms(); ++i) {
      REQUIRE(!m->getAtomWithIdx(i)->getIsAromatic());
    }
  }

  {
    for (auto smi : {"N1****1", "C1=C*2C=CC=C*2C=C1", "N1*C=CC=C1",
                     "C1=CC2=CC=C3C=CC4=CC=C5C=CN1*1*2*3*4*51"}) {
      std::unique_ptr<RWMol> m(SmilesToMol(smi));
      REQUIRE(m);
      for (const auto b : m->bonds()) {
        REQUIRE(!b->getIsAromatic());
      }
    }
  }

  {
    ROMOL_SPTR m = "C1=CC2=CC=C3C=CC4=CC=C5C=CN1*1*2*3*4N51"_smiles;
    for (const auto a : m->atoms()) {
      if (a->getIdx() == 16 || a->getIdx() == 17) {
        REQUIRE(!a->getIsAromatic());
      } else {
        REQUIRE(a->getIsAromatic());
      }
    }
    unsigned int nNonAromaticBonds = 0;
    for (const auto b : m->bonds()) {
      if (!b->getIsAromatic()) {
        ++nNonAromaticBonds;
      }
    }
    REQUIRE(nNonAromaticBonds == 5);
  }

  {
    for (auto smi :
         {"*1C=CC=C1", "N1*=**=*1", "C1=CC2=CC=C3C=CC4=CC=C5C=CN1N1N2N3N4N51",
          "C1=CC2=CC=C3C=CC4=CC=C5C=CN1*1N2*3N4N51"}) {
      std::unique_ptr<RWMol> m(SmilesToMol(smi));
      REQUIRE(m);
      for (const auto b : m->bonds()) {
        REQUIRE(b->getIsAromatic());
      }
    }
  }
}

TEST_CASE("Testing sf.net issue 2208994 : kekulization error") {
  {
    std::string smi = "Cn1ccc(=O)n1C";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getIsAromatic() == true);
    REQUIRE(m->getBondWithIdx(1)->getIsAromatic() == true);

    delete m;
  }

  {
    std::string smi = "c:1:c:c:c:c:c1";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getIsAromatic() == true);
    REQUIRE(m->getBondWithIdx(1)->getIsAromatic() == true);

    delete m;
  }

  {
    std::string smi = "c1:c:c:c:c:c:1";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getIsAromatic() == true);
    REQUIRE(m->getBondWithIdx(1)->getIsAromatic() == true);

    delete m;
  }
}

TEST_CASE("Testing sf.net issue 2313979: aromaticity assignment hangs") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/test_data/";
  SDMolSupplier suppl(pathName + "Issue2313979.sdf", false);

  while (!suppl.atEnd()) {
    ROMol *m = suppl.next();
    REQUIRE(m);
    std::string nm;
    m->getProp(common_properties::_Name, nm);
    BOOST_LOG(rdInfoLog) << "   Doing molecule: " << nm << std::endl;

    BOOST_LOG(rdInfoLog) << "     This should finish in a few seconds.  >>>"
                         << std::endl;
    MolOps::sanitizeMol(*(RWMol *)m);
    delete m;
    BOOST_LOG(rdInfoLog) << "   <<< Done." << std::endl;
  }
}

TEST_CASE("Testing sf.net issue 2316677 : canonicalization error") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/test_data/";
  RWMol *m = MolFileToMol(pathName + "Issue2316677.mol");
  REQUIRE(m);
  std::string smi = MolToSmiles(*m, true);
  std::cerr << "smi: " << smi << std::endl;
  REQUIRE(smi == "Cc1ccc(S(=O)(=O)/N=C2\\CC(=N\\C(C)(C)C)/C2=N\\C(C)(C)C)cc1");
  delete m;
}

TEST_CASE(
    "Testing sf.net issue 2830244: make sure that non-ring aromatic atoms generate errors") {
  {
    std::string smi = "c-C";

    RWMol *m = SmilesToMol(smi, 0, false);
    bool ok = false;
    try {
      MolOps::Kekulize(*m);
    } catch (MolSanitizeException &) {
      ok = true;
    }
    REQUIRE(ok);
    delete m;
  }
  {
    std::string smi = "c-C";

    RWMol *m = SmilesToMol(smi, 0, false);
    bool ok = false;
    unsigned int opThatFailed;
    try {
      MolOps::sanitizeMol(*m, opThatFailed);
    } catch (MolSanitizeException &) {
      ok = true;
    }
    REQUIRE(ok);
    REQUIRE(opThatFailed == MolOps::SANITIZE_KEKULIZE);
    delete m;
  }
}

TEST_CASE(
    "Testing sf.net issue 2951221 : hydrogens added with bad coordinates") {
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    ROMol *m = MolFileToMol(pathName + "Issue2951221.1.mol");
    REQUIRE(m);
    REQUIRE(m->getConformer().is3D());
    ROMol *m2 = MolOps::addHs(*m, false, true);
    REQUIRE(m2);
    delete m;
    REQUIRE(m2->getNumAtoms(false) == 12);
    RDGeom::Point3D coords[4];
    coords[0] = m2->getConformer().getAtomPos(2);
    coords[1] = m2->getConformer().getAtomPos(0);
    coords[2] = m2->getConformer().getAtomPos(1);
    coords[3] = m2->getConformer().getAtomPos(9);
    double dot =
        (coords[3] - coords[0])
            .dotProduct(
                (coords[1] - coords[0]).crossProduct(coords[2] - coords[0]));
    REQUIRE(dot > 1.0);
    delete m2;
  }

  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    ROMol *m = MolFileToMol(pathName + "Issue2951221.2.mol");
    REQUIRE(m);
    ROMol *m2 = MolOps::addHs(*m, false, true);
    REQUIRE(m2);
    delete m;
    REQUIRE(m2->getNumAtoms(false) == 5);
    MolOps::assignChiralTypesFrom3D(*m2);
    MolOps::assignStereochemistry(*m2, true, true);
    REQUIRE(m2->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    std::string cip;
    m2->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    REQUIRE(cip == "S");
    delete m2;
  }
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    ROMol *m = MolFileToMol(pathName + "Issue2951221.3.mol");
    REQUIRE(m);
    ROMol *m2 = MolOps::addHs(*m, false, true);
    REQUIRE(m2);
    delete m;
    REQUIRE(m2->getNumAtoms(false) == 5);
    MolOps::assignChiralTypesFrom3D(*m2);
    MolOps::assignStereochemistry(*m2, true, true);
    REQUIRE(m2->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    std::string cip;
    m2->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    REQUIRE(cip == "R");
    delete m2;
  }
}

TEST_CASE(
    "Testing sf.net issue 2952255 : bad assignment of radicals to early elements") {
  {
    std::string smi = "[C](C)(C)C";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
    delete m;
  }
  {
    std::string smi = "[C](C)C";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 2);
    delete m;
  }
  {
    std::string smi = "[CH](C)C";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
    delete m;
  }
  {
    std::string smi = "[CH+](C)C";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 0);
    delete m;
  }
  {
    std::string smi = "[C-](C)C";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
    delete m;
  }
  {
    std::string smi = "C(C)(C)(C)C";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 0);
    delete m;
  }
  {
    std::string smi = "[N](C)C";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
    delete m;
  }
  {
    std::string smi = "[N+](C)(C)C";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
    delete m;
  }
  {
    std::string smi = "[Cl]";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
    delete m;
  }
  {
    std::string smi = "[Cl-]";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 0);
    delete m;
  }
  {
    std::string smi = "[Cl]C";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 0);
    delete m;
  }
  {
    std::string smi = "[Na]";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
    delete m;
  }
  {
    std::string smi = "[Na+]";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 0);
    delete m;
  }
  {
    std::string smi = "[Na]C";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 0);
    delete m;
  }
  {
    std::string smi = "[Mg+]C";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 0);
    delete m;
  }
  {
    std::string smi = "[Mg]C";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
    delete m;
  }
  {
    std::string smi = "[Mg+]";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
    delete m;
  }
  {
    std::string smi = "[Mg+2]";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 0);
    delete m;
  }
}

TEST_CASE("Testing sf.net issue 3185548 : problems with SSSR code") {
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    BOOST_LOG(rdInfoLog) << "  Starting file read 1" << std::endl;
    RWMol *m = MolFileToMol(pathName + "Issue3185548.mol");
    BOOST_LOG(rdInfoLog) << "  finished" << std::endl;
    REQUIRE(m);
    delete m;
  }

  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    BOOST_LOG(rdInfoLog) << "  Starting file read 2" << std::endl;
    RWMol *m = MolFileToMol(pathName + "Issue3185548.2.mol");
    BOOST_LOG(rdInfoLog) << "  finished" << std::endl;
    REQUIRE(m);

    m->getRingInfo()->reset();
    unsigned int nsssr;
    VECT_INT_VECT sssrs;
    nsssr = MolOps::findSSSR(*m, sssrs);
    REQUIRE(nsssr == 48);
    nsssr = MolOps::symmetrizeSSSR(*m, sssrs);
    REQUIRE(nsssr == 56);
    delete m;
  }
}

TEST_CASE("Testing Issue 3349243") {
  {
    std::string smi = "c1cccc[n+]1";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    MolOps::Kekulize(*m);
    // just finishing is good
    REQUIRE(m->getBondWithIdx(0)->getBondType() != Bond::AROMATIC);
    delete m;
  }
}

TEST_CASE("Testing fast find rings") {
  {
    std::string smi = "CCC";
    RWMol *m = SmilesToMol(smi, 0, 0);
    REQUIRE(m);
    MolOps::fastFindRings(*m);
    REQUIRE(m->getRingInfo());
    REQUIRE(m->getRingInfo()->isInitialized());
    REQUIRE(m->getRingInfo()->numRings() == 0);
    delete m;
  }
  {
    std::string smi = "C1CC1";
    RWMol *m = SmilesToMol(smi, 0, 0);
    REQUIRE(m);
    MolOps::fastFindRings(*m);
    REQUIRE(m->getRingInfo());
    REQUIRE(m->getRingInfo()->isInitialized());
    REQUIRE(m->getRingInfo()->numRings() == 1);
    delete m;
  }

  {
    std::string smi = "CC1CC1";
    RWMol *m = SmilesToMol(smi, 0, 0);
    REQUIRE(m);
    MolOps::fastFindRings(*m);
    REQUIRE(m->getRingInfo());
    REQUIRE(m->getRingInfo()->isInitialized());
    REQUIRE(m->getRingInfo()->numRings() == 1);
    delete m;
  }

  {
    std::string smi = "C1CC1.C1CC1";
    RWMol *m = SmilesToMol(smi, 0, 0);
    REQUIRE(m);
    MolOps::fastFindRings(*m);
    REQUIRE(m->getRingInfo());
    REQUIRE(m->getRingInfo()->isInitialized());
    REQUIRE(m->getRingInfo()->numRings() == 2);
    delete m;
  }
  {
    std::string smi = "C1C(C)C1";
    RWMol *m = SmilesToMol(smi, 0, 0);
    REQUIRE(m);
    MolOps::fastFindRings(*m);
    REQUIRE(m->getRingInfo());
    REQUIRE(m->getRingInfo()->isInitialized());
    REQUIRE(m->getRingInfo()->numRings() == 1);
    delete m;
  }
  {
    std::string smi = "c1c(=O)nc2[nH]cnn2c1O";
    RWMol *m = SmilesToMol(smi, 0, 0);
    REQUIRE(m);
    MolOps::fastFindRings(*m);
    REQUIRE(m->getRingInfo());
    REQUIRE(m->getRingInfo()->isInitialized());
    REQUIRE(m->getRingInfo()->numRings() == 2);
    delete m;
  }
}

TEST_CASE("Testing Issue 3487473") {
  {
    std::string smi = "C*C";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getHybridization() == Atom::UNSPECIFIED);
    delete m;
  }

  {
    std::string smi = "C*C";
    RWMol *m = SmartsToMol(smi);
    REQUIRE(m);
    m->updatePropertyCache(false);
    MolOps::setConjugation(*m);
    REQUIRE(m->getAtomWithIdx(1)->getHybridization() == Atom::UNSPECIFIED);
    delete m;
  }
}

TEST_CASE("Testing Issue 3480481") {
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName + "Issue3480481.mol");
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(0)->getIsAromatic() == true);
    REQUIRE(m->getAtomWithIdx(0)->getValence(Atom::ValenceType::EXPLICIT) == 4);
    REQUIRE(m->getAtomWithIdx(0)->getValence(Atom::ValenceType::IMPLICIT) == 0);
    REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == -1);
    delete m;
  }
}

void aamatchtest(std::string smi1, std::string smi2, bool shouldMatch, int idx1,
                 int idx2) {
  RWMol *m1 = SmilesToMol(smi1);
  RWMol *m2 = SmilesToMol(smi2);
  REQUIRE(m1);
  REQUIRE(m2);
  // std::cerr<<"   "<<smi1<<" "<<smi2<<std::endl;
  REQUIRE(m2->getAtomWithIdx(idx2)->Match(m1->getAtomWithIdx(idx1)) ==
          shouldMatch);
  delete m1;
  delete m2;
}

TEST_CASE("Testing Atom-Atom matching behavior") {
  /* Here's what we're testing:

     | Molecule | Query   | Match |
     | CCO      | CCO     | Yes   |
     | CC[O-]   | CCO     | Yes   |
     | CCO      | CC[O-]  | No    |
     | CC[O-]   | CC[O-]  | Yes   |
     | CC[O-]   | CC[OH]  | Yes   |
     | CCOC     | CC[OH]  | Yes   |
     | CCOC     | CCO     | Yes   |
     | CCC      | CCC     | Yes   |
     | CC[14C]  | CCC     | Yes   |
     | CCC      | CC[14C] | No    |
     | CC[14C]  | CC[14C] | Yes   |
     | OCO      | C       | Yes   |
     | OCO      | [CH2]   | Yes   |
     | OCO      | [CH3]   | Yes   |
     | O[CH2]O  | C       | Yes   |
     | O[CH2]O  | [CH2]   | Yes   |
     | OCO      | [CH]    | Yes   |

     This is a large superset of issue 3495370

  */

  // note that in some cases here we have to be fairly careful about Hs on
  // the query to make sure that it doesn't have radicals (radical handling
  // added to fix github #165
  aamatchtest("CCO", "O", true, 2, 0);
  aamatchtest("CC[O-]", "O", true, 2, 0);
  aamatchtest("CCO", "[OH-]", false, 2, 0);
  aamatchtest("CC[O-]", "[OH-]", true, 2, 0);
  aamatchtest("CC[O-]", "[OH2]", true, 2, 0);
  aamatchtest("CCOC", "[OH2]", true, 2, 0);
  aamatchtest("CCOC", "O", true, 2, 0);
  aamatchtest("CCC", "C", true, 2, 0);
  aamatchtest("CC[14C]", "C", true, 2, 0);
  aamatchtest("CCC", "[14CH4]", false, 2, 0);
  aamatchtest("CC[14C]", "[14CH4]", true, 2, 0);
  aamatchtest("CC[13C]", "[14CH4]", false, 2, 0);
  aamatchtest("OCO", "C", true, 1, 0);
  aamatchtest("OCO", "[CH4]", true, 1, 0);
  aamatchtest("O[CH2]O", "C", true, 1, 0);
  aamatchtest("O[CH2]O", "[CH4]", true, 1, 0);
  aamatchtest("OCO", "[CH2]", false, 1,
              0);  // doesn't match due to radical count
  aamatchtest("O[CH2]O", "[CH2]", false, 1,
              0);  // doesn't match due to radical count
  aamatchtest("O[CH]O", "[CH3]", true, 1, 0);
  aamatchtest("O[CH]O", "[CH2]", false, 1,
              0);  // doesn't match due to radical count
  aamatchtest("CC", "*", false, 1, 0);
  aamatchtest("C*", "*", true, 1, 0);
  aamatchtest("C[1*]", "*", true, 1, 0);
  aamatchtest("C[1*]", "[1*]", true, 1, 0);
  aamatchtest("C*", "[1*]", true, 1, 0);
  aamatchtest("C[2*]", "[1*]", false, 1, 0);
}

TEST_CASE("Testing Issue 3525076") {
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName + "Issue3525076.sdf");
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(18)->getIsAromatic() == false);
    REQUIRE(m->getBondWithIdx(18)->getBondType() == Bond::SINGLE);
    MolOps::Kekulize(*m);
    REQUIRE(m->getBondWithIdx(18)->getIsAromatic() == false);
    REQUIRE(m->getBondWithIdx(18)->getBondType() == Bond::SINGLE);
    MolOps::sanitizeMol(*m);
    REQUIRE(m->getBondWithIdx(18)->getIsAromatic() == false);
    REQUIRE(m->getBondWithIdx(18)->getBondType() == Bond::SINGLE);
    delete m;
  }
}

TEST_CASE("Testing canonicalization basics") {
  // these are all cases that were problematic at one time or another during
  // the canonicalization rewrite.
  {
    std::string smi = "FC1C(=C/Cl)\\C1";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getBondBetweenAtoms(2, 3)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondBetweenAtoms(2, 3)->getStereo() == Bond::STEREOZ);

    // std::cerr<<"-------------\n";
    std::string csmi1 = MolToSmiles(*m, true);
    // std::cerr<<csmi1<<std::endl;
    RWMol *m2 = SmilesToMol(csmi1);
    REQUIRE(m2);

    MatchVectType mv;
    REQUIRE(SubstructMatch(*m, *m2, mv));
    std::map<int, int> mmap;
    for (MatchVectType::const_iterator mvit = mv.begin(); mvit != mv.end();
         ++mvit) {
      mmap[mvit->second] = mvit->first;
    }
    REQUIRE(m2->getBondBetweenAtoms(mmap[2], mmap[3])->getBondType() ==
            Bond::DOUBLE);
    REQUIRE(m2->getBondBetweenAtoms(mmap[2], mmap[3])->getStereo() ==
            Bond::STEREOZ);

    // std::cerr<<"-------------\n";
    std::string csmi2 = MolToSmiles(*m2, true);
    // std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    REQUIRE(csmi1 == csmi2);
    delete m;
    delete m2;
  }

  {
    std::string smi = "CC1(C)C2CCC1(C)C(=O)/C2=C\\C(N=N/c1ccccc1)=N/Nc1ccccc1";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getBondBetweenAtoms(10, 11)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondBetweenAtoms(10, 11)->getStereo() == Bond::STEREOZ);
    REQUIRE(m->getBondBetweenAtoms(12, 21)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondBetweenAtoms(12, 21)->getStereo() == Bond::STEREOE);
    REQUIRE(m->getBondBetweenAtoms(13, 14)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondBetweenAtoms(13, 14)->getStereo() == Bond::STEREONONE);
    // std::cerr<<"-------------\n";
    std::string csmi1 = MolToSmiles(*m, true);

    RWMol *m2 = SmilesToMol(csmi1);
    REQUIRE(m2);

    MatchVectType mv;
    REQUIRE(SubstructMatch(*m, *m2, mv));
    std::map<int, int> mmap;
    for (MatchVectType::const_iterator mvit = mv.begin(); mvit != mv.end();
         ++mvit) {
      mmap[mvit->second] = mvit->first;
    }
    REQUIRE(m2->getBondBetweenAtoms(mmap[10], mmap[11])->getBondType() ==
            Bond::DOUBLE);
    REQUIRE(m2->getBondBetweenAtoms(mmap[10], mmap[11])->getStereo() ==
            Bond::STEREOZ);
    REQUIRE(m2->getBondBetweenAtoms(mmap[12], mmap[21])->getBondType() ==
            Bond::DOUBLE);
    REQUIRE(m2->getBondBetweenAtoms(mmap[12], mmap[21])->getStereo() ==
            Bond::STEREOE);
    REQUIRE(m2->getBondBetweenAtoms(mmap[13], mmap[14])->getBondType() ==
            Bond::DOUBLE);
    REQUIRE(m2->getBondBetweenAtoms(mmap[13], mmap[14])->getStereo() ==
            Bond::STEREONONE);

    // std::cerr<<"-------------\n";
    std::string csmi2 = MolToSmiles(*m2, true);
    // std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    REQUIRE(csmi1 == csmi2);
    delete m;
    delete m2;
  }

  {
    std::string smi = "COc1ccc(OC)c2[nH]c(=O)cc(C)c21";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    // std::cerr<<"-------------\n";
    std::string csmi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(csmi1);
    REQUIRE(m);
    // std::cerr<<"-------------\n";
    std::string csmi2 = MolToSmiles(*m, true);
    // std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    REQUIRE(csmi1 == csmi2);
    delete m;
  }
  {
    std::string smi = "COc1cc(C)c(C(=O)[O-])cc1OC";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    // std::cerr<<"-------------\n";
    std::string csmi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(csmi1);
    REQUIRE(m);
    // std::cerr<<"-------------\n";
    std::string csmi2 = MolToSmiles(*m, true);
    // std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    REQUIRE(csmi1 == csmi2);
    delete m;
  }
  {
    std::string smi = "COc1ccc(C(=O)OC(c2ccc(OC)cc2)C(C)O)cc1";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    // std::cerr<<"-------------\n";
    std::string csmi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(csmi1);
    REQUIRE(m);
    // std::cerr<<"-------------\n";
    std::string csmi2 = MolToSmiles(*m, true);
    // std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    REQUIRE(csmi1 == csmi2);
    delete m;
  }
  {
    std::string smi = "CC(C)C1CCC(C)=CC1=NNC(N)=O";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    // std::cerr<<"-------------\n";
    std::string csmi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(csmi1);
    REQUIRE(m);
    // std::cerr<<"-------------\n";
    std::string csmi2 = MolToSmiles(*m, true);
    // std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    REQUIRE(csmi1 == csmi2);
    delete m;
  }
  {
    std::string smi = "COCCNC(=O)c1ccccc1N1C(=O)C2(C)c3[nH]c4ccccc4c3CCN2C1=O";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    // std::cerr<<"-------------\n";
    std::string csmi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(csmi1);
    REQUIRE(m);
    // std::cerr<<"-------------\n";
    std::string csmi2 = MolToSmiles(*m, true);
    // std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    REQUIRE(csmi1 == csmi2);
    delete m;
  }
  {
    std::string smi = "Cc1c(Br)cc(Br)cc1C(F)(F)F";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    // std::cerr<<"-------------\n";
    std::string csmi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(csmi1);
    REQUIRE(m);
    // std::cerr<<"-------------\n";
    std::string csmi2 = MolToSmiles(*m, true);
    // std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    REQUIRE(csmi1 == csmi2);
    delete m;
  }
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName + "zinc4235774a.mol");
    REQUIRE(m);
    REQUIRE(m->getBondBetweenAtoms(1, 2)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondBetweenAtoms(1, 2)->getStereo() == Bond::STEREOZ);
    REQUIRE(m->getBondBetweenAtoms(7, 8)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondBetweenAtoms(7, 8)->getStereo() == Bond::STEREOZ);
    std::string smi = MolToSmiles(*m, true);
    // std::cerr<<"SMILES: "<<smi<<std::endl;
    RWMol *m2 = SmilesToMol(smi);
    MatchVectType mv;
    REQUIRE(SubstructMatch(*m, *m2, mv));
    std::map<int, int> mmap;
    for (MatchVectType::const_iterator mvit = mv.begin(); mvit != mv.end();
         ++mvit) {
      mmap[mvit->second] = mvit->first;
    }
    REQUIRE(m2->getBondBetweenAtoms(mmap[1], mmap[2])->getBondType() ==
            Bond::DOUBLE);
    REQUIRE(m2->getBondBetweenAtoms(mmap[7], mmap[8])->getBondType() ==
            Bond::DOUBLE);
    REQUIRE(m2->getBondBetweenAtoms(mmap[1], mmap[2])->getStereo() ==
            Bond::STEREOZ);
    REQUIRE(m2->getBondBetweenAtoms(mmap[7], mmap[8])->getStereo() ==
            Bond::STEREOZ);
    delete m;
    delete m2;
  }
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName + "zinc4235774.mol");
    REQUIRE(m);
    REQUIRE(m->getBondBetweenAtoms(4, 5)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondBetweenAtoms(4, 5)->getStereo() == Bond::STEREOZ);
    REQUIRE(m->getBondBetweenAtoms(14, 15)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondBetweenAtoms(14, 15)->getStereo() == Bond::STEREOZ);
    std::string smi = MolToSmiles(*m, true);
    RWMol *m2 = SmilesToMol(smi);
    MatchVectType mv;
    REQUIRE(SubstructMatch(*m, *m2, mv));
    std::map<int, int> mmap;
    for (MatchVectType::const_iterator mvit = mv.begin(); mvit != mv.end();
         ++mvit) {
      mmap[mvit->second] = mvit->first;
    }
    REQUIRE(m2->getBondBetweenAtoms(mmap[4], mmap[5])->getBondType() ==
            Bond::DOUBLE);
    REQUIRE(m2->getBondBetweenAtoms(mmap[14], mmap[15])->getBondType() ==
            Bond::DOUBLE);
    REQUIRE(m2->getBondBetweenAtoms(mmap[4], mmap[5])->getStereo() ==
            Bond::STEREOZ);
    REQUIRE(m2->getBondBetweenAtoms(mmap[14], mmap[15])->getStereo() ==
            Bond::STEREOZ);
    delete m;
    delete m2;
  }

  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName + "zinc3850436piece.mol");
    REQUIRE(m);
    std::string csmi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(csmi1);
    // std::cerr<<"-------------\n";
    std::string csmi2 = MolToSmiles(*m, true);
    // std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    REQUIRE(csmi1 == csmi2);
    delete m;
  }
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName + "zinc13403961piece.mol");
    REQUIRE(m);
    REQUIRE(m->getBondBetweenAtoms(1, 2)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondBetweenAtoms(1, 2)->getStereo() == Bond::STEREOZ);
    REQUIRE(m->getBondBetweenAtoms(3, 7)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondBetweenAtoms(3, 7)->getStereo() == Bond::STEREOE);
    REQUIRE(m->getBondBetweenAtoms(4, 5)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondBetweenAtoms(4, 5)->getStereo() == Bond::STEREOE);
    std::string csmi1 = MolToSmiles(*m, true);
    // std::cerr<<"SMI1: "<<csmi1<<std::endl;
    RWMol *m2 = SmilesToMol(csmi1);
    MatchVectType mv;
    REQUIRE(SubstructMatch(*m, *m2, mv));
    std::map<int, int> mmap;
    for (MatchVectType::const_iterator mvit = mv.begin(); mvit != mv.end();
         ++mvit) {
      mmap[mvit->second] = mvit->first;
    }

    REQUIRE(m2->getBondBetweenAtoms(mmap[1], mmap[2])->getBondType() ==
            Bond::DOUBLE);
    REQUIRE(m2->getBondBetweenAtoms(mmap[1], mmap[2])->getStereo() ==
            Bond::STEREOZ);
    REQUIRE(m2->getBondBetweenAtoms(mmap[3], mmap[7])->getBondType() ==
            Bond::DOUBLE);
    REQUIRE(m2->getBondBetweenAtoms(mmap[3], mmap[7])->getStereo() ==
            Bond::STEREOE);
    REQUIRE(m2->getBondBetweenAtoms(mmap[4], mmap[5])->getBondType() ==
            Bond::DOUBLE);
    REQUIRE(m2->getBondBetweenAtoms(mmap[4], mmap[5])->getStereo() ==
            Bond::STEREOE);

    // std::cerr<<"-------------\n";
    std::string csmi2 = MolToSmiles(*m2, true);
    // std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    REQUIRE(csmi1 == csmi2);
    delete m;
    delete m2;
  }
  {
    std::string smi = "C\\N=c1/s/c(=N\\Cl)/c/1=N/F";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    // std::cerr<<"-------------\n";
    std::string csmi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(csmi1);
    REQUIRE(m);
    // std::cerr<<"-------------\n";
    std::string csmi2 = MolToSmiles(*m, true);
    // std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    REQUIRE(csmi1 == csmi2);
    delete m;
  }
  {
    std::string smi =
        "Cc1ccc(S(=O)(=O)/N=c2sc(=N\\C(C)(C)C)/c\\2=N/C(C)(C)C)cc1";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    // std::cerr<<"-------------\n";
    std::string csmi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(csmi1);
    REQUIRE(m);
    // std::cerr<<"-------------\n";
    std::string csmi2 = MolToSmiles(*m, true);
    // std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    REQUIRE(csmi1 == csmi2);
    delete m;
  }
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName + "zinc6624278.mol");
    REQUIRE(m);
    REQUIRE(m->getBondBetweenAtoms(21, 13)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondBetweenAtoms(21, 13)->getStereo() == Bond::STEREOE);
    REQUIRE(m->getBondBetweenAtoms(5, 12)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondBetweenAtoms(5, 12)->getStereo() == Bond::STEREOZ);

    std::string csmi1 = MolToSmiles(*m, true);
    // std::cerr<<"SMI1: "<<csmi1<<std::endl;
    RWMol *m2 = SmilesToMol(csmi1);
    MatchVectType mv;
    REQUIRE(SubstructMatch(*m, *m2, mv));
    std::map<int, int> mmap;
    for (MatchVectType::const_iterator mvit = mv.begin(); mvit != mv.end();
         ++mvit) {
      mmap[mvit->second] = mvit->first;
    }

    REQUIRE(m2->getBondBetweenAtoms(mmap[21], mmap[13])->getBondType() ==
            Bond::DOUBLE);
    REQUIRE(m2->getBondBetweenAtoms(mmap[21], mmap[13])->getStereo() ==
            Bond::STEREOE);
    REQUIRE(m2->getBondBetweenAtoms(mmap[5], mmap[12])->getBondType() ==
            Bond::DOUBLE);
    REQUIRE(m2->getBondBetweenAtoms(mmap[5], mmap[12])->getStereo() ==
            Bond::STEREOZ);

    // std::cerr<<"-------------\n";
    std::string csmi2 = MolToSmiles(*m2, true);
    // std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    REQUIRE(csmi1 == csmi2);
    delete m2;

    std::string tsmi = MolToSmiles(*m, true, false, 7, false);
    // std::cerr<<"-------------\n";
    // std::cerr<<"T:\n"<<tsmi<<"\n-------------\n"<<std::endl;
    m2 = SmilesToMol(tsmi);
    REQUIRE(SubstructMatch(*m, *m2, mv));
    mmap.clear();
    for (MatchVectType::const_iterator mvit = mv.begin(); mvit != mv.end();
         ++mvit) {
      mmap[mvit->second] = mvit->first;
    }
    REQUIRE(m2->getBondBetweenAtoms(mmap[21], mmap[13])->getBondType() ==
            Bond::DOUBLE);
    REQUIRE(m2->getBondBetweenAtoms(mmap[21], mmap[13])->getStereo() ==
            Bond::STEREOE);
    REQUIRE(m2->getBondBetweenAtoms(mmap[5], mmap[12])->getBondType() ==
            Bond::DOUBLE);
    REQUIRE(m2->getBondBetweenAtoms(mmap[5], mmap[12])->getStereo() ==
            Bond::STEREOZ);

    // std::cerr<<"-------------\n";
    csmi2 = MolToSmiles(*m2, true);
    // std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    REQUIRE(csmi1 == csmi2);
    delete m2;
    delete m;
  }
  {
    std::string smi = "F/C=C/C=C(C)/C=C/Cl";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getBondBetweenAtoms(1, 2)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondBetweenAtoms(1, 2)->getStereo() == Bond::STEREOE);
    REQUIRE(m->getBondBetweenAtoms(3, 4)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondBetweenAtoms(3, 4)->getStereo() == Bond::STEREOE);
    REQUIRE(m->getBondBetweenAtoms(6, 7)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondBetweenAtoms(6, 7)->getStereo() == Bond::STEREOE);

    std::string tsmi = MolToSmiles(*m, true, false, 3, false);
    // std::cerr<<"-------------\n";
    // std::cerr<<"T:\n"<<tsmi<<"\n-------------\n"<<std::endl;
    RWMol *m2 = SmilesToMol(tsmi);
    MatchVectType mv;
    REQUIRE(SubstructMatch(*m, *m2, mv));
    std::map<int, int> mmap;
    mmap.clear();
    for (MatchVectType::const_iterator mvit = mv.begin(); mvit != mv.end();
         ++mvit) {
      mmap[mvit->second] = mvit->first;
    }
    REQUIRE(m2->getBondBetweenAtoms(mmap[1], mmap[2])->getBondType() ==
            Bond::DOUBLE);
    REQUIRE(m2->getBondBetweenAtoms(mmap[1], mmap[2])->getStereo() ==
            Bond::STEREOE);
    REQUIRE(m2->getBondBetweenAtoms(mmap[3], mmap[4])->getBondType() ==
            Bond::DOUBLE);
    REQUIRE(m2->getBondBetweenAtoms(mmap[3], mmap[4])->getStereo() ==
            Bond::STEREOE);
    REQUIRE(m2->getBondBetweenAtoms(mmap[6], mmap[7])->getBondType() ==
            Bond::DOUBLE);
    REQUIRE(m2->getBondBetweenAtoms(mmap[6], mmap[7])->getStereo() ==
            Bond::STEREOE);

    std::string csmi1 = MolToSmiles(*m, true);
    // std::cerr<<"SMI1: "<<csmi1<<std::endl;
    // std::cerr<<"-------------\n";
    std::string csmi2 = MolToSmiles(*m2, true);
    // std::cerr<<csmi1<<"\n"<<csmi2<<"\n-------------\n"<<std::endl;
    REQUIRE(csmi1 == csmi2);

    delete m2;
    delete m;
  }

  {
    // this was issue 3528556
    std::string smi = "N12.N13.C24.C35.C46.C56";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    std::string csmi = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(csmi);
    REQUIRE(m);
    smi = MolToSmiles(*m, true);
    REQUIRE(csmi == smi);
    delete m;
  }
  {
    // this was issue 3526831
    std::string smi = "CO/N=C/C(=C(\\O)/c1ccc(Cl)cc1)/C=N\\OC";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    std::string csmi = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(csmi);
    REQUIRE(m);
    smi = MolToSmiles(*m, true);
    REQUIRE(csmi == smi);
    delete m;
  }
}

TEST_CASE("Testing sf.net issue 3549146: problems after mergeQueryHs") {
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    ROMol *m = MolFileToMol(pathName + "Issue3549146.mol", true, false);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 16);
    ROMol *m2 = MolOps::mergeQueryHs(*m);
    REQUIRE(m2);
    REQUIRE(m2->getNumAtoms() == 13);
    REQUIRE(!(m2->getRingInfo()->isInitialized()));
    delete m;
    delete m2;
  }
  {
    std::string smi = "CCC.C";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 4);
    REQUIRE((m->getRingInfo()->isInitialized()));
    m->addBond(1, 3, Bond::SINGLE);
    REQUIRE(m->getAtomWithIdx(1)->needsUpdatePropertyCache());
    REQUIRE(m->getAtomWithIdx(3)->needsUpdatePropertyCache());
    REQUIRE((m->getRingInfo()->isInitialized()));
    m->addBond(0, 2, Bond::SINGLE);
    REQUIRE((m->getRingInfo()->isInitialized()));
    REQUIRE(m->getAtomWithIdx(0)->needsUpdatePropertyCache());
    REQUIRE(m->getAtomWithIdx(2)->needsUpdatePropertyCache());
    delete m;
  }
  {
    std::string smi = "C1CC1C";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 4);
    REQUIRE((m->getRingInfo()->isInitialized()));
    m->removeBond(2, 3);
    REQUIRE(!(m->getRingInfo()->isInitialized()));
    delete m;
  }
  {
    std::string smi = "C1CC1C";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 4);
    REQUIRE((m->getRingInfo()->isInitialized()));
    m->removeAtom(3);
    REQUIRE(!(m->getRingInfo()->isInitialized()));
    delete m;
  }
}

TEST_CASE("Testing sf.net issue 249: finding rings consumes all memory") {
  {
    std::string smi =
        "Cc1cc2cc(c1)C(=O)NCc1cc-3cc(CNC(=O)c4cc(C)cc(c4)C(=O)NCc4cc(cc(CNC2=O)"
        "c4O)-c2cc4CNC(=O)c5cc(C)cc(c5)C(=O)NCc5cc-3cc(CNC(=O)c3cc(C)cc(c3)C(="
        "O)NCc(c2)c4O)c5O)c1O";
    ROMol *m = SmilesToMol(smi, 0, 0);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 88);
    m->updatePropertyCache(false);
    std::cerr << "starting ring finding" << std::endl;
    MolOps::findSSSR(*m);
    std::cerr << "done" << std::endl;

    delete m;
  }

  {
    std::string smi =
        "CCCOc1c2CNC(=O)c3cc(cc(c3)C(=O)NCc3cc4cc(CNC(=O)c5cc(C(=O)NCc1cc(c2)"
        "c1cc2CNC(=O)c6cc(cc(c6)C(=O)NCc6cc4cc(CNC(=O)c4cc(C(=O)NCc(c1)c2OCCC)"
        "cc(c4)C(=O)NC(COCCC(=O)O)(COCCC(=O)O)COCCC(=O)O)c6OCCC)C(=O)NC(COCCC(="
        "O)O)(COCCC(=O)O)COCCC(=O)O)cc(c5)C(=O)NC(COCCC(=O)O)(COCCC(=O)O)COCCC("
        "=O)O)c3OCCC)C(=O)NC(COCCC(=O)O)(COCCC(=O)O)COCCC(=O)O";
    ROMol *m = SmilesToMol(smi, 0, 0);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 196);
    m->updatePropertyCache(false);
    std::cerr << "starting ring finding" << std::endl;
    MolOps::findSSSR(*m);
    std::cerr << "done" << std::endl;

    delete m;
  }

  {
    std::string smi =
        "CCn1nnc(c1)CN=C(C1CC2C3CCC4C5C3C3C6C2C2C1C1CCC7C(C1)C1C8C9C7C(C(=O)O)"
        "C(C(=O)O)C7C9C(C9C8C8C%10C1C1C(C2C2C6C6C%11C3C(C5)C3C(C(=O)O)C5C%"
        "12CC9C9C8C8C(C%10)C%10C%13C(C%14C(C2C1C(=O)O)C6C1C2C%11C3C3C5C(C5C%"
        "12C9C(C8CC%10)CC5)C(CC3C2C(C(C1C%14CC%13C(=NCc1nnn(c1)CC)O)C(=O)O)C(="
        "O)O)C(=NCc1nnn(c1)CC)O)C(=O)O)C(=O)O)CC1C(C4C(CC71)C(=NCc1nnn(c1)CC)O)"
        "C(=O)O)O";
    ROMol *m = SmilesToMol(smi, 0, 0);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 167);
    m->updatePropertyCache(false);
    std::cerr << "starting ring finding" << std::endl;
    MolOps::findSSSR(*m);
    std::cerr << "done" << std::endl;
    delete m;
  }

  {
    std::string smi =
        "C/C=N/"
        "c1ccc(cc1)c1cc2cc(c1)c1ccc(cc1)N=Cc1ccc(cc1)c1cc(cc(c1)c1ccc(cc1)C)"
        "c1ccc(cc1)C=Nc1ccc(cc1)c1cc(cc(c1)c1ccc(cc1)/N=C/"
        "C)c1ccc(cc1)N=Cc1ccc(cc1)c1cc(c3ccc(C=Nc4ccc(c5cc6c7ccc(N=Cc8ccc(c9cc("
        "c%10ccc(C=Nc%11ccc2cc%11)cc%10)cc(c9)c2ccc(cc2)C=Nc2ccc(cc2)c2cc(cc("
        "c2)c2ccc(cc2)N=Cc2ccc(cc2)c2cc(cc(c2)c2ccc(cc2)C)c2ccc(cc2)C=Nc2ccc("
        "cc2)c2cc(c9ccc(N=Cc%10ccc(c%11cc(c%12ccc(C=Nc%13ccc(c(c6)c5)cc%13)cc%"
        "12)cc(c%11)c5ccc(cc5)C)cc%10)cc9)cc(c2)c2ccc(cc2)/N=C/C)c2ccc(cc2)/"
        "N=C/C)cc8)cc7)cc4)cc3)cc(c1)c1ccc(cc1)C";
    RWMol *m = SmilesToMol(smi, 0, 0);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 278);
    std::cerr << smi << std::endl;
    std::cerr << "starting sanitization" << std::endl;
    MolOps::sanitizeMol(*m);
    std::cerr << "done" << std::endl;
    delete m;
  }

  {
    std::string smi =
        "COc1cc2ccc1n1c(=O)n(c3ccc(cc3OC)c3ccc(c(c3)OC)n3c(=O)n(c4ccc(cc4OC)"
        "c4ccc(c(c4)OC)n4c(=O)n(c5ccc(cc5OC)c5ccc(n6c(=O)n(c7ccc(c8ccc(n9c(=O)"
        "n(c%10ccc(c%11ccc(n%12c(=O)n(c%13ccc2cc%13OC)c(=O)n(c%12=O)c2ccc("
        "cc2OC)C)c(OC)c%11)cc%10OC)c(=O)n(c9=O)c2ccc(cc2OC)C)c(OC)c8)cc7OC)c(="
        "O)n(c6=O)c2ccc(cc2OC)C)c(c5)OC)c(=O)n(c4=O)c2ccc(cc2OC)C)c(=O)n(c3=O)"
        "c2ccc(cc2OC)C)c(=O)n(c1=O)c1ccc(cc1OC)C";
    RWMol *m = SmilesToMol(smi, 0, 0);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 204);
    std::cerr << "starting sanitization" << std::endl;
    MolOps::sanitizeMol(*m);
    std::cerr << "done" << std::endl;
    delete m;
  }
}

TEST_CASE("Testing sf.net issue 256: bad atom counts") {
  {
    std::string smi = "*CC[H]";
    ROMol *m = SmilesToMol(smi, 0, 0);
    REQUIRE(m);
    m->updatePropertyCache(false);
    REQUIRE(m->getNumAtoms() == 4);
    REQUIRE(m->getNumAtoms(false) == 8);
    REQUIRE(m->getNumHeavyAtoms() == 2);
    delete m;
  }
  {
    std::string smi = "*CC[2H]";
    ROMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 4);
    REQUIRE(m->getNumAtoms(false) == 8);
    REQUIRE(m->getNumHeavyAtoms() == 2);
    delete m;
  }
}

TEST_CASE("Testing sf.net issue 266: ring finding error") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/test_data/";
  RWMol *m = MolFileToMol(pathName + "Issue266.mol", false);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 19);
  REQUIRE(m->getNumBonds() == 25);
  std::cerr << "starting ring finding" << std::endl;
  MolOps::findSSSR(*m);
  std::cerr << "done" << std::endl;
  REQUIRE(m->getRingInfo()->numRings() ==
          (m->getNumBonds() - m->getNumAtoms() + 1));
  delete m;
}

TEST_CASE("Testing sf.net issue 272: removing two-coordinate Hs") {
  {
    std::string smi = "C[H-]C";
    ROMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 3);
    delete m;
  }
  {
    std::string smi = "C[H].C";
    ROMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 2);
    delete m;
  }
}

TEST_CASE(
    "Testing Github issue 8 (impact of removeAtom on bond stereo atoms)") {
  {
    std::string smi = "Cl/C=C/Cl";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    MolOps::assignStereochemistry(*m);
    REQUIRE(m->getBondWithIdx(1)->getStereoAtoms().size() == 2);
    m->removeAtom((unsigned int)0);
    REQUIRE(m->getBondWithIdx(1)->getStereoAtoms().size() == 0);
    delete m;
  }
  {
    std::string smi = "CC/C=C/Cl";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    MolOps::assignStereochemistry(*m);
    INT_VECT &sas = m->getBondWithIdx(2)->getStereoAtoms();
    REQUIRE(sas.size() == 2);
    REQUIRE(std::find(sas.begin(), sas.end(), 1) != sas.end());
    m->removeAtom((unsigned int)0);
    REQUIRE(m->getBondWithIdx(1)->getStereoAtoms().size() == 2);
    REQUIRE(std::find(sas.begin(), sas.end(), 0) != sas.end());
    REQUIRE(std::find(sas.begin(), sas.end(), 1) == sas.end());
    delete m;
  }
  {
    std::string smi = "C/C=C/CC";
    RWMol *m = SmilesToMol(smi);
    REQUIRE(m);
    MolOps::assignStereochemistry(*m);
    INT_VECT &sas = m->getBondWithIdx(1)->getStereoAtoms();
    REQUIRE(sas.size() == 2);
    REQUIRE(std::find(sas.begin(), sas.end(), 0) != sas.end());
    REQUIRE(std::find(sas.begin(), sas.end(), 3) != sas.end());
    m->removeAtom((unsigned int)4);
    REQUIRE(m->getBondWithIdx(1)->getStereoAtoms().size() == 2);
    REQUIRE(std::find(sas.begin(), sas.end(), 0) != sas.end());
    REQUIRE(std::find(sas.begin(), sas.end(), 3) != sas.end());
    delete m;
  }
}

TEST_CASE("Testing Github issue 42 (impact of removeAtom on atom stereochem)") {
  std::string smi =
      "CCN1CCN(c2cc3[nH]c(C(=O)[C@@]4(CC)CC[C@](C)(O)CC4)nc3cc2Cl)CC1";
  RWMol *m = SmilesToMol(smi);
  REQUIRE(m);
  int indices[] = {29, 28, 27, 26, 25, 24, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1};
  for (unsigned int i = 0; indices[i] > -1; ++i) {
    m->removeAtom((unsigned int)indices[i]);
  }
  smi = MolToSmiles(*m, true);
  std::cerr << "smiles: " << smi << std::endl;
  delete m;
}

TEST_CASE(
    "Testing Github issue 65 (kekulization of boron-containing aromatic rings)") {
  std::string smi = "C[B-]1=CC=CC=C1";
  RWMol *m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getAtomWithIdx(1)->getIsAromatic());
  REQUIRE(m->getBondWithIdx(1)->getIsAromatic());

  MolOps::Kekulize(*m);

  delete m;
}

TEST_CASE(
    "Testing Github issue 72 (problems with bad benzothiazolium structure)") {
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName + "github72.mol");
    REQUIRE(m);
    REQUIRE(!m->getBondBetweenAtoms(0, 8)->getIsAromatic());
    REQUIRE(m->getBondBetweenAtoms(1, 6)->getIsAromatic());
    delete m;
  }

  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName + "github72.2.mol");
    REQUIRE(m);
    REQUIRE(m->getBondBetweenAtoms(0, 8)->getIsAromatic());
    REQUIRE(m->getBondBetweenAtoms(1, 6)->getIsAromatic());
    delete m;
  }
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName + "github72.3.mol");
    REQUIRE(m);
    REQUIRE(!m->getBondBetweenAtoms(0, 8)->getIsAromatic());
    REQUIRE(m->getBondBetweenAtoms(1, 6)->getIsAromatic());

    std::string smi = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smi);
    REQUIRE(m);

    delete m;
  }
}

namespace {
void _renumberTest(const ROMol *m) {
  PRECONDITION(m, "no molecule");
  std::vector<unsigned int> idxV(m->getNumAtoms());
  for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
    idxV[i] = i;
  }

  std::string refSmi = MolToSmiles(*m, true);
  for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
    std::vector<unsigned int> nVect(idxV);
    std::shuffle(nVect.begin(), nVect.end(), std::mt19937(0xf00d));
    // std::copy(nVect.begin(),nVect.end(),std::ostream_iterator<int>(std::cerr,",
    // "));
    // std::cerr<<std::endl;

    ROMol *nm = MolOps::renumberAtoms(*m, nVect);
    REQUIRE(nm);
    REQUIRE(nm->getNumAtoms() == m->getNumAtoms());
    REQUIRE(nm->getNumBonds() == m->getNumBonds());

    // checking the SSS is a test for Github #317
    MatchVectType mv;
    REQUIRE(SubstructMatch(*m, *nm, mv));
    REQUIRE(mv.size() == nm->getNumAtoms());

    for (unsigned int j = 0; j < m->getNumAtoms(); ++j) {
      REQUIRE(m->getAtomWithIdx(nVect[j])->getAtomicNum() ==
              nm->getAtomWithIdx(j)->getAtomicNum());
    }

    // checking the conformation is a test for Github #441
    REQUIRE(m->getNumConformers() == nm->getNumConformers());
    if (m->getNumConformers()) {
      for (unsigned int j = 0; j < m->getNumAtoms(); ++j) {
        RDGeom::Point3D po = m->getConformer().getAtomPos(nVect[j]);
        RDGeom::Point3D pn = nm->getConformer().getAtomPos(j);
        REQUIRE(po.x == pn.x);
        REQUIRE(po.y == pn.y);
        REQUIRE(po.z == pn.z);
      }
      // checking conformer dimensionality is a test for Github #584
      REQUIRE(m->getConformer().is3D() == nm->getConformer().is3D());
    }

    std::string nSmi = MolToSmiles(*nm, true);
    if (nSmi != refSmi) {
      std::cerr << refSmi << std::endl << nSmi << std::endl;
    }
    REQUIRE(nSmi == refSmi);
    delete nm;
  }
}
}  // namespace

TEST_CASE("Testing renumbering atoms") {
  {
    std::string smiles = "CC1CCCC(C)C1C";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    _renumberTest(m);
    delete m;
  }
  {
    std::string smiles = "C[C@H]1C[C@H](F)C1";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    _renumberTest(m);
    delete m;
  }
  {
    std::string smiles = "C[C@H]1CC[C@H](C/C=C/[C@H](F)Cl)CC1";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    _renumberTest(m);
    delete m;
  }
  {  // github issue #441 and #584
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName + "Issue266.mol");  // no significance to
                                                         // choice of files, we
                                                         // just need something
                                                         // with 2D coords
    REQUIRE(m);
    _renumberTest(m);
    delete m;
  }

  {  // github issue 1735 renumber empty molecules
    auto *m = new ROMol;
    REQUIRE(m);
    std::vector<unsigned int> nVect;
    auto *m1 = MolOps::renumberAtoms(*m, nVect);
    delete m;
    delete m1;
  }
}
TEST_CASE(
    "Testing github issue 141: Kekulization of molecule with aromatic N leaves the explicit H there") {
  std::string smiles = "N1C=CC=C1";
  RWMol *m = SmilesToMol(smiles);
  REQUIRE(m);
  MolOps::Kekulize(*m, true);
  m->updatePropertyCache(true);
  REQUIRE(!m->getAtomWithIdx(0)->getIsAromatic());
  REQUIRE(m->getAtomWithIdx(0)->getNumImplicitHs() == 1);
  REQUIRE(m->getAtomWithIdx(0)->getNumExplicitHs() == 0);

  delete m;
}

TEST_CASE("Testing ZBO basics") {
  auto *m = new RWMol();

  m->addAtom(new Atom(26), true, true);
  m->addAtom(new Atom(6), true, true);
  m->addAtom(new Atom(6), true, true);
  m->addAtom(new Atom(6), true, true);
  m->addAtom(new Atom(6), true, true);
  m->addAtom(new Atom(6), true, true);
  m->addAtom(new Atom(6), true, true);
  m->addBond(1, 2, Bond::AROMATIC);
  m->addBond(2, 3, Bond::AROMATIC);
  m->addBond(3, 4, Bond::AROMATIC);
  m->addBond(4, 5, Bond::AROMATIC);
  m->addBond(5, 6, Bond::AROMATIC);
  m->addBond(1, 6, Bond::AROMATIC);

  m->addBond(1, 0, Bond::ZERO);
  m->addBond(2, 0, Bond::ZERO);
  m->addBond(3, 0, Bond::ZERO);
  m->addBond(4, 0, Bond::ZERO);
  m->addBond(5, 0, Bond::ZERO);
  m->addBond(6, 0, Bond::ZERO);

  MolOps::sanitizeMol(*m);

  REQUIRE(m->getRingInfo()->numAtomRings(0) == 0);
  REQUIRE(m->getRingInfo()->numAtomRings(1) == 1);

  REQUIRE(m->getAtomWithIdx(1)->getHybridization() == Atom::SP2);
  REQUIRE(m->getAtomWithIdx(2)->getHybridization() == Atom::SP2);
  REQUIRE(m->getAtomWithIdx(3)->getHybridization() == Atom::SP2);
  REQUIRE(m->getAtomWithIdx(4)->getHybridization() == Atom::SP2);
  REQUIRE(m->getAtomWithIdx(5)->getHybridization() == Atom::SP2);
  REQUIRE(m->getAtomWithIdx(6)->getHybridization() == Atom::SP2);

  REQUIRE(m->getBondWithIdx(0)->getIsAromatic());
  REQUIRE(m->getBondWithIdx(1)->getIsAromatic());
  REQUIRE(m->getBondWithIdx(2)->getIsAromatic());
  REQUIRE(m->getBondWithIdx(3)->getIsAromatic());
  REQUIRE(m->getBondWithIdx(4)->getIsAromatic());
  REQUIRE(m->getBondWithIdx(5)->getIsAromatic());
  delete m;
}

TEST_CASE("Testing operator= on molecules") {
  std::string smi =
      "CCN1CCN(c2cc3[nH]c(C(=O)[C@@]4(CC)CC[C@](C)(O)CC4)nc3cc2Cl)CC1";
  RWMol *m = SmilesToMol(smi);
  REQUIRE(m);
  std::string csmi = MolToSmiles(*m, true);

  RWMol m2 = *m;
  std::string nsmi = MolToSmiles(m2, true);
  REQUIRE(nsmi == csmi);

  RWMol *m3 = SmilesToMol("C2CC2[C@H](F)Cl");
  REQUIRE(m3);
  *m3 = *m;
  nsmi = MolToSmiles(*m3, true);
  REQUIRE(nsmi == csmi);
  delete m3;
  delete m;
}

TEST_CASE("Testing github issue 190: Don't merge Hs onto dummy atoms") {
  std::string smiles = "*[H]";
  RWMol *m = SmilesToMol(smiles);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 2);
  delete m;
}

namespace {
inline int getAtNum(const ROMol &, const Atom *at) {
  return at->getAtomicNum();
}
}  // namespace
TEST_CASE("Testing getMolFragsWithQuery()") {
  {
    std::string smiles = "C1CCC1ONNC";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 8);
    std::map<int, boost::shared_ptr<ROMol>> res =
        MolOps::getMolFragsWithQuery(*m, getAtNum);
    REQUIRE(res.size() == 3);
    REQUIRE(res.find(6) != res.end());
    REQUIRE(res.find(7) != res.end());
    REQUIRE(res.find(8) != res.end());
    REQUIRE(res.find(5) == res.end());
    REQUIRE(res[6]->getNumAtoms() == 5);
    REQUIRE(res[6]->getNumBonds() == 4);
    REQUIRE(res[7]->getNumAtoms() == 2);
    REQUIRE(res[7]->getNumBonds() == 1);
    REQUIRE(res[8]->getNumAtoms() == 1);
    REQUIRE(res[8]->getNumBonds() == 0);

    std::map<int, std::unique_ptr<ROMol>> otherRes;
    MolOps::getMolFragsWithQuery(*m, getAtNum, otherRes);
    REQUIRE(otherRes.size() == 3);
    REQUIRE(otherRes.find(6) != otherRes.end());
    REQUIRE(otherRes.find(7) != otherRes.end());
    REQUIRE(otherRes.find(8) != otherRes.end());
    REQUIRE(otherRes.find(5) == otherRes.end());
    REQUIRE(otherRes[6]->getNumAtoms() == 5);
    REQUIRE(otherRes[6]->getNumBonds() == 4);
    REQUIRE(otherRes[7]->getNumAtoms() == 2);
    REQUIRE(otherRes[7]->getNumBonds() == 1);
    REQUIRE(otherRes[8]->getNumAtoms() == 1);
    REQUIRE(otherRes[8]->getNumBonds() == 0);
    delete m;
  }
  {
    std::string smiles = "C1CCC1ONNC";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 8);
    std::vector<int> keep;
    keep.push_back(6);
    keep.push_back(8);
    std::map<int, boost::shared_ptr<ROMol>> res =
        MolOps::getMolFragsWithQuery(*m, getAtNum, true, &keep);
    REQUIRE(res.size() == 2);
    REQUIRE(res.find(6) != res.end());
    REQUIRE(res.find(7) == res.end());
    REQUIRE(res.find(8) != res.end());
    REQUIRE(res[6]->getNumAtoms() == 5);
    REQUIRE(res[6]->getNumBonds() == 4);
    REQUIRE(res[8]->getNumAtoms() == 1);
    REQUIRE(res[8]->getNumBonds() == 0);
    std::map<int, std::unique_ptr<ROMol>> otherRes;
    MolOps::getMolFragsWithQuery(*m, getAtNum, otherRes, true, &keep);
    REQUIRE(otherRes.size() == 2);
    REQUIRE(otherRes.find(6) != otherRes.end());
    REQUIRE(otherRes.find(7) == otherRes.end());
    REQUIRE(otherRes.find(8) != otherRes.end());
    REQUIRE(otherRes[6]->getNumAtoms() == 5);
    REQUIRE(otherRes[6]->getNumBonds() == 4);
    REQUIRE(otherRes[8]->getNumAtoms() == 1);
    REQUIRE(otherRes[8]->getNumBonds() == 0);

    delete m;
  }
  {
    std::string smiles = "C1CCC1ONNC";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 8);
    std::vector<int> keep;
    keep.push_back(6);
    keep.push_back(8);
    std::map<int, boost::shared_ptr<ROMol>> res =
        MolOps::getMolFragsWithQuery(*m, getAtNum, true, &keep, true);
    REQUIRE(res.size() == 1);
    REQUIRE(res.find(6) == res.end());
    REQUIRE(res.find(7) != res.end());
    REQUIRE(res.find(8) == res.end());
    REQUIRE(res[7]->getNumAtoms() == 2);
    REQUIRE(res[7]->getNumBonds() == 1);

    std::map<int, std::unique_ptr<ROMol>> otherRes;
    MolOps::getMolFragsWithQuery(*m, getAtNum, otherRes, true, &keep, true);
    REQUIRE(otherRes.size() == 1);
    REQUIRE(otherRes.find(6) == otherRes.end());
    REQUIRE(otherRes.find(7) != otherRes.end());
    REQUIRE(otherRes.find(8) == otherRes.end());
    REQUIRE(otherRes[7]->getNumAtoms() == 2);
    REQUIRE(otherRes[7]->getNumBonds() == 1);

    delete m;
  }
}

TEST_CASE("Testing github issue 418: removeHs not updating H count") {
  {
    auto *m2 = new RWMol();
    m2->addAtom(new Atom(7), true, true);
    m2->addAtom(new Atom(1), true, true);
    m2->addAtom(new Atom(1), true, true);
    m2->addAtom(new Atom(1), true, true);
    m2->addAtom(new Atom(1), true, true);
    m2->addBond(0, 1, Bond::SINGLE);
    m2->addBond(0, 2, Bond::SINGLE);
    m2->addBond(0, 3, Bond::SINGLE);
    m2->addBond(0, 4, Bond::SINGLE);
    MolOps::removeHs(*m2, false, true, false);
    REQUIRE(m2->getAtomWithIdx(0)->getNumExplicitHs() == 4);
    delete m2;
  }
  {
    std::string smiles = "[H][N+]([H])([H])[H]";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 1);
    REQUIRE(m->getAtomWithIdx(0)->getNumExplicitHs() == 4);
    delete m;
  }
  {
    std::string smiles = "[H]N([H])([H])[H]";
    bool ok = false;
    try {
      SmilesToMol(smiles);
    } catch (MolSanitizeException &) {
      ok = true;
    }
    REQUIRE(ok);
  }
}

TEST_CASE(
    "Testing github issue 432: problems caused by aromatic Ns with radical electrons") {
  {
    std::string smiles = "C1=NN=N[N]1";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(4)->getNumRadicalElectrons() == 1);
    REQUIRE(!m->getAtomWithIdx(4)->getIsAromatic());
    REQUIRE(!m->getAtomWithIdx(0)->getIsAromatic());
    REQUIRE(!m->getBondWithIdx(0)->getIsAromatic());
    delete m;
  }
  {  // test round-tripping:
    std::string smiles = "C1=NN=N[N]1";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    smiles = MolToSmiles(*m);
    delete m;
    m = SmilesToMol(smiles);
    REQUIRE(m);
    delete m;
  }
  {  // test round-tripping:
    std::string smiles = "OC(=O)C(=O)Nc1cccc(c1)C2=NN=N[N]2";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    smiles = MolToSmiles(*m);
    delete m;
    m = SmilesToMol(smiles);
    REQUIRE(m);
    delete m;
  }
  {
    std::string smiles = "C1=C[N]C=C1";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(2)->getNumRadicalElectrons() == 1);
    REQUIRE(!m->getAtomWithIdx(2)->getIsAromatic());
    REQUIRE(!m->getAtomWithIdx(0)->getIsAromatic());
    REQUIRE(!m->getBondWithIdx(0)->getIsAromatic());
    delete m;
  }
}

TEST_CASE(
    "Testing github issue 443: kekulization problems caused by any bonds") {
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName + "github443.min.mol");
    REQUIRE(m);
    REQUIRE(!m->getAtomWithIdx(0)->getIsAromatic());
    REQUIRE(m->getBondBetweenAtoms(0, 3));
    MolOps::Kekulize(*m);
    delete m;
  }

  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName + "github443.mol");
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(19)->getIsAromatic());
    REQUIRE(m->getAtomWithIdx(12)->getIsAromatic());
    // we might normally expect these to be aromatic because the outer porphyrin
    // ring
    // is 4n+2 aromatic. However, the current fused ring aromaticity perception
    // uses
    // the symmetrized SSSR rings and only works if all atoms are aromatic. This
    // cannot
    // happen when the Mg is involved
    // REQUIRE(m->getAtomWithIdx(13)->getIsAromatic());
    // REQUIRE(m->getAtomWithIdx(11)->getIsAromatic());
    REQUIRE(m->getBondBetweenAtoms(13, 20));
    REQUIRE(m->getBondBetweenAtoms(19, 20));
    REQUIRE(m->getBondBetweenAtoms(11, 20));
    REQUIRE(m->getBondBetweenAtoms(12, 20));
    MolOps::Kekulize(*m);

    delete m;
  }
}

TEST_CASE(
    "Testing github issue 447: Radicals are not correctly assigned when reading from SMILES") {
  {
    std::string smiles = "C[S]";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getNoImplicit());
    REQUIRE(m->getAtomWithIdx(1)->getNumRadicalElectrons() == 1);
    delete m;
  }
  {
    std::string smiles = "C[SH]C";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getNoImplicit());
    REQUIRE(m->getAtomWithIdx(1)->getNumRadicalElectrons() == 1);
    delete m;
  }
  {
    std::string smiles = "C[SH3]C";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getNoImplicit());
    REQUIRE(m->getAtomWithIdx(1)->getNumRadicalElectrons() == 1);
    delete m;
  }
  {
    std::string smiles = "C[SH4]C";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getNoImplicit());
    REQUIRE(m->getAtomWithIdx(1)->getNumRadicalElectrons() == 0);
    delete m;
  }
  {
    std::string smiles = "C[SH3]C";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getNoImplicit());
    REQUIRE(m->getAtomWithIdx(1)->getNumRadicalElectrons() == 1);
    delete m;
  }
  {
    std::string smiles = "C[SH2+]C";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getNoImplicit());
    REQUIRE(m->getAtomWithIdx(1)->getNumRadicalElectrons() == 1);
    delete m;
  }

  {
    std::string smiles = "C[P]C";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getNoImplicit());
    REQUIRE(m->getAtomWithIdx(1)->getNumRadicalElectrons() == 1);
    delete m;
  }
  {
    std::string smiles = "C[PH2]C";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->getNoImplicit());
    REQUIRE(m->getAtomWithIdx(1)->getNumRadicalElectrons() == 1);
    delete m;
  }
}

TEST_CASE("Testing generation of new molecules from molecule fragments") {
  {
    std::string smiles = "c1ccccc1.O.CCC(=O)O";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);

    INT_VECT fragsMapping;
    VECT_INT_VECT fragsMolAtomMapping;
    std::vector<ROMOL_SPTR> frags =
        MolOps::getMolFrags(*m, false, &fragsMapping, &fragsMolAtomMapping);

    REQUIRE(frags.size() == 3);
    REQUIRE(fragsMapping.size() == m->getNumAtoms());

    REQUIRE(fragsMapping[2] == 0);
    REQUIRE(fragsMapping[6] == 1);
    REQUIRE(fragsMapping[8] == 2);
    REQUIRE(fragsMolAtomMapping[0].size() == frags[0]->getNumAtoms());
    REQUIRE(fragsMolAtomMapping[1].size() == frags[1]->getNumAtoms());
    REQUIRE(fragsMolAtomMapping[2].size() == frags[2]->getNumAtoms());
    REQUIRE(fragsMolAtomMapping[0][1] == 1);
    REQUIRE(fragsMolAtomMapping[1][0] == 6);
    REQUIRE(fragsMolAtomMapping[2][1] == 8);

    REQUIRE(MolToSmiles(*frags[0], true) == "c1ccccc1");
    REQUIRE(MolToSmiles(*frags[1], true) == "O");
    REQUIRE(MolToSmiles(*frags[2], true) == "CCC(=O)O");
    delete m;
  }
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    RWMol *m = MolFileToMol(pathName + "chembl1203199.mol");
    REQUIRE(m);
    std::string smi = "C[C@H](NC(=O)[C@H]1Cc2c(sc3ccccc23)CN1)c1ccccc1.Cl";
    REQUIRE(MolToSmiles(*m, true) == smi);

    INT_VECT fragsMapping;
    VECT_INT_VECT fragsMolAtomMapping;
    std::vector<ROMOL_SPTR> frags = MolOps::getMolFrags(
        *m, false, &fragsMapping, &fragsMolAtomMapping, true);

    REQUIRE(frags.size() == 2);
    REQUIRE(fragsMapping.size() == m->getNumAtoms());
    REQUIRE(fragsMapping[2] == 0);
    REQUIRE(fragsMapping[24] == 1);
    REQUIRE(fragsMolAtomMapping[0].size() == frags[0]->getNumAtoms());
    REQUIRE(fragsMolAtomMapping[1].size() == frags[1]->getNumAtoms());
    REQUIRE(fragsMolAtomMapping[0][1] == 1);
    REQUIRE(fragsMolAtomMapping[1][0] == 24);

    REQUIRE(frags[0]->getNumConformers() == 1);
    REQUIRE(frags[1]->getNumConformers() == 1);

    REQUIRE(frags[0]->getConformer(0).getAtomPos(0).x ==
            m->getConformer(0).getAtomPos(0).x);
    REQUIRE(frags[0]->getConformer(0).getAtomPos(0).y ==
            m->getConformer(0).getAtomPos(0).y);
    REQUIRE(frags[0]->getConformer(0).getAtomPos(0).z ==
            m->getConformer(0).getAtomPos(0).z);

    REQUIRE(frags[0]->getConformer(0).getAtomPos(3).x ==
            m->getConformer(0).getAtomPos(3).x);
    REQUIRE(frags[0]->getConformer(0).getAtomPos(3).y ==
            m->getConformer(0).getAtomPos(3).y);
    REQUIRE(frags[0]->getConformer(0).getAtomPos(3).z ==
            m->getConformer(0).getAtomPos(3).z);

    REQUIRE(frags[1]->getConformer(0).getAtomPos(0).x ==
            m->getConformer(0).getAtomPos(24).x);
    REQUIRE(frags[1]->getConformer(0).getAtomPos(0).y ==
            m->getConformer(0).getAtomPos(24).y);
    REQUIRE(frags[1]->getConformer(0).getAtomPos(0).z ==
            m->getConformer(0).getAtomPos(24).z);
    delete m;
  }
  {  // confirm bond-only stereogroups are not removed during GetmolFrags
    std::string smiles =
        "Cc1cccc(Cl)c1-c1c(C)cccc1I.Cc1cccc(F)c1-c1c(C)cccc1Cl |wD:8.15,wU:23.23,o1:23,&1:8|";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);

    INT_VECT fragsMapping;
    VECT_INT_VECT fragsMolAtomMapping;
    std::vector<ROMOL_SPTR> frags =
        MolOps::getMolFrags(*m, false, &fragsMapping, &fragsMolAtomMapping);
    REQUIRE(frags.size() == 2);
    REQUIRE(fragsMapping.size() == m->getNumAtoms());

    for (const auto &frag : frags) {
      REQUIRE(frag->getNumAtoms() == 16);
      REQUIRE(frag->getNumBonds() == 17);
      REQUIRE(frag->getStereoGroups().size() == 1);
      REQUIRE(frag->getStereoGroups()[0].getBonds().size() == 1);
    }

    RDKit::SmilesWriteParams sps;
    REQUIRE(MolToCXSmiles(*frags[0], sps,
                          SmilesWrite::CXSmilesFields::CX_ALL_BUT_COORDS) ==
            "Cc1cccc(Cl)c1-c1c(C)cccc1I |wU:7.6,&1:7|");
    REQUIRE(MolToCXSmiles(*frags[1], sps,
                          SmilesWrite::CXSmilesFields::CX_ALL_BUT_COORDS) ==
            "Cc1cccc(F)c1-c1c(C)cccc1Cl |wU:7.6,o1:7|");
    delete m;
  }
}

namespace {
void hypervalent_check(const char *smiles) {
  RWMol *m = SmilesToMol(smiles);
  REQUIRE(m);
  REQUIRE(m->getBondBetweenAtoms(0, 1)->getBondType() == Bond::SINGLE);
  REQUIRE(m->getAtomWithIdx(1)->getFormalCharge() == -1);
  delete m;
}
}  // namespace

TEST_CASE("Testing github issue 510: Hexafluorophosphate cannot be handled") {
  hypervalent_check("F[P-](F)(F)(F)(F)F");
  // test #1668 too, it's the same thing but with As, Sb, and Bi
  hypervalent_check("F[As-](F)(F)(F)(F)F");
  hypervalent_check("F[Sb-](F)(F)(F)(F)F");
  hypervalent_check("F[Bi-](F)(F)(F)(F)F");

  hypervalent_check("F[Sb-](F)(F)(F)(F)F");
  hypervalent_check("F[Bi-](F)(F)(F)(F)F");

  // we also added a valence of 5 for Bi:
  hypervalent_check("F[Bi-](F)(F)F");
}

TEST_CASE(
    "Testing github issue 526: Bad ring finding in a complex fused ring") {
  {
    std::string smiles = "N1C2[C@@H]3N[C@H]4[C@@H]5N[C@@H]([C@@H]1C35)C24";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getRingInfo()->numRings() == 6);
    delete m;
  }
  {
    std::string smiles =
        "NN1C2C3[C@@H]4[C@@H]1C1[C@H]2N([C@H]3[C@@H]1N4N1C(=O)C2=C(C=CC=C2)C1="
        "O)N1C(=O)C2=C(C=CC=C2)C1=O";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getRingInfo()->numRings() == 10);
    delete m;
  }
}

TEST_CASE(
    "Testing github issue 539: Lack of conjugation in allyl cations, lack of aromaticity perception/ability to kekulize aromatic carbocations such as cyclopropenyl and tropylium") {
  {
    std::string smiles = "C=C-[CH2+]";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    bool allConjugated = true;
    for (unsigned int i = 0; allConjugated && i < m->getNumBonds(); ++i) {
      allConjugated = m->getBondWithIdx(i)->getIsConjugated();
    }
    REQUIRE(allConjugated);
    delete m;
  }
  {
    std::vector<std::string> smilesVec;
    smilesVec.emplace_back("C1=C[CH+]1");
    smilesVec.emplace_back("C1=CC=C[CH+]C=C1");
    smilesVec.emplace_back("c1c[cH+]1");
    smilesVec.emplace_back("c1ccc[cH+]cc1");
    for (std::vector<std::string>::const_iterator smiles = smilesVec.begin();
         smiles != smilesVec.end(); ++smiles) {
      RWMol *m = SmilesToMol(*smiles);
      REQUIRE(m);
      bool allConjugated = true;
      for (unsigned int i = 0; allConjugated && i < m->getNumBonds(); ++i) {
        allConjugated = m->getBondWithIdx(i)->getIsConjugated();
      }
      REQUIRE(allConjugated);
      bool allAromatic = true;
      for (unsigned int i = 0; allAromatic && i < m->getNumBonds(); ++i) {
        allAromatic = m->getBondWithIdx(i)->getIsAromatic();
      }
      REQUIRE(allAromatic);
      delete m;
    }
  }

  {
    std::string smiles = "C=C-C";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(2)->getHybridization() == Atom::SP3);
    REQUIRE(m->getBondWithIdx(1)->getBondType() == Bond::SINGLE);
    REQUIRE(!m->getBondWithIdx(1)->getIsConjugated());
    delete m;
  }

  {
    std::string smiles = "C=C-O";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(2)->getHybridization() == Atom::SP2);
    REQUIRE(m->getBondWithIdx(1)->getBondType() == Bond::SINGLE);
    REQUIRE(m->getBondWithIdx(1)->getIsConjugated());
    delete m;
  }

  {
    std::string smiles = "C=C-N";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(2)->getHybridization() == Atom::SP2);
    REQUIRE(m->getBondWithIdx(1)->getBondType() == Bond::SINGLE);
    REQUIRE(m->getBondWithIdx(1)->getIsConjugated());
    delete m;
  }

  {
    std::string smiles = "C=C-[NH3+]";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(2)->getHybridization() == Atom::SP3);
    REQUIRE(m->getBondWithIdx(1)->getBondType() == Bond::SINGLE);
    REQUIRE(!m->getBondWithIdx(1)->getIsConjugated());
    delete m;
  }
  {
    std::string smiles = "Cc1ccccc1";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(0)->getBondType() == Bond::SINGLE);
    // the bond to the CH3 should not be conjugated, but the others are
    REQUIRE(!m->getBondWithIdx(0)->getIsConjugated());
    for (unsigned int i = 1; i < m->getNumBonds(); ++i) {
      REQUIRE(m->getBondWithIdx(i)->getIsConjugated());
    }
    delete m;
  }
  {
    std::string smiles = "Fc1c[nH]c(=O)[nH]c1=O";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(0)->getBondType() == Bond::SINGLE);
    // the bond to the F should not be conjugated, but the others are
    REQUIRE(!m->getBondWithIdx(0)->getIsConjugated());
    for (unsigned int i = 1; i < m->getNumBonds(); ++i) {
      REQUIRE(m->getBondWithIdx(i)->getIsConjugated());
    }
    delete m;
  }
}

TEST_CASE("Testing adjustQueryProperties()") {
  {  // basics from SMILES
    std::string smiles = "C1CCC1C";
    ROMol *qm = SmilesToMol(smiles);
    REQUIRE(qm);
    REQUIRE(qm->getNumAtoms() == 5);
    ROMol *aqm = MolOps::adjustQueryProperties(*qm);
    REQUIRE(aqm);
    REQUIRE(aqm->getNumAtoms() == 5);
    {
      smiles = "C1CCC1CC";
      ROMol *m = SmilesToMol(smiles);
      REQUIRE(m);
      MatchVectType match;
      REQUIRE(SubstructMatch(*m, *qm, match));
      REQUIRE(SubstructMatch(*m, *aqm, match));
      delete m;
    }
    {
      smiles = "C1C(C)CC1CC";
      ROMol *m = SmilesToMol(smiles);
      REQUIRE(m);
      MatchVectType match;
      REQUIRE(SubstructMatch(*m, *qm, match));
      REQUIRE(!SubstructMatch(*m, *aqm, match));
      delete m;
    }

    delete qm;
    delete aqm;
  }
  {  // basics from SMARTS
    std::string smiles = "C1CCC1*";
    ROMol *qm = SmartsToMol(smiles);
    REQUIRE(qm);
    REQUIRE(qm->getNumAtoms() == 5);
    ROMol *aqm = MolOps::adjustQueryProperties(*qm);
    REQUIRE(aqm);
    REQUIRE(aqm->getNumAtoms() == 5);
    {
      smiles = "C1CCC1CC";
      ROMol *m = SmilesToMol(smiles);
      REQUIRE(m);
      MatchVectType match;
      REQUIRE(SubstructMatch(*m, *qm, match));
      REQUIRE(SubstructMatch(*m, *aqm, match));
      delete m;
    }
    {
      smiles = "C1C(C)CC1CC";
      ROMol *m = SmilesToMol(smiles);
      REQUIRE(m);
      MatchVectType match;
      REQUIRE(SubstructMatch(*m, *qm, match));
      REQUIRE(!SubstructMatch(*m, *aqm, match));
      delete m;
    }

    delete qm;
    delete aqm;
  }

  {
    std::string smiles = "C1CC(*)C1*";
    ROMol *qm = SmartsToMol(smiles);
    REQUIRE(qm);
    REQUIRE(qm->getNumAtoms() == 6);
    ROMol *aqm = MolOps::adjustQueryProperties(*qm);
    REQUIRE(aqm);
    REQUIRE(aqm->getNumAtoms() == 6);
    {
      smiles = "C1CC2C1CC2";
      ROMol *m = SmilesToMol(smiles);
      REQUIRE(m);
      MatchVectType match;
      REQUIRE(SubstructMatch(*m, *qm, match));
      REQUIRE(SubstructMatch(*m, *aqm, match));

      MolOps::AdjustQueryParameters aqp;

      delete aqm;
      aqp.adjustDegree = false;
      aqp.adjustRingCount = false;
      aqm = MolOps::adjustQueryProperties(*qm, &aqp);
      REQUIRE(aqm);
      REQUIRE(aqm->getNumAtoms() == 6);
      REQUIRE(SubstructMatch(*m, *aqm, match));

      delete aqm;
      aqp.adjustDegree = false;
      aqp.adjustRingCount = true;
      aqm = MolOps::adjustQueryProperties(*qm, &aqp);
      REQUIRE(aqm);
      REQUIRE(aqm->getNumAtoms() == 6);
      REQUIRE(!SubstructMatch(*m, *aqm, match));

      delete aqm;
      aqp.adjustDegree = true;
      aqp.adjustRingCount = false;
      aqp.adjustDegreeFlags = MolOps::ADJUST_IGNORENONE;
      aqm = MolOps::adjustQueryProperties(*qm, &aqp);
      REQUIRE(aqm);
      REQUIRE(aqm->getNumAtoms() == 6);
      REQUIRE(!SubstructMatch(*m, *aqm, match));

      delete m;
    }

    {
      smiles = "C1CC(C2CC2)C1C2CC2";
      ROMol *m = SmilesToMol(smiles);
      REQUIRE(m);
      MatchVectType match;
      delete aqm;
      aqm = MolOps::adjustQueryProperties(*qm);
      REQUIRE(SubstructMatch(*m, *qm, match));
      REQUIRE(SubstructMatch(*m, *aqm, match));

      MolOps::AdjustQueryParameters aqp;

      delete aqm;
      aqp.adjustRingCount = true;
      aqm = MolOps::adjustQueryProperties(*qm, &aqp);
      REQUIRE(aqm);
      REQUIRE(aqm->getNumAtoms() == 6);
      REQUIRE(SubstructMatch(*m, *aqm, match));

      delete aqm;
      aqp.adjustRingCountFlags =
          MolOps::ADJUST_IGNORENONE;  // neither "not dummy"
                                      // nor "in ring"
                                      // restrictions
      aqm = MolOps::adjustQueryProperties(*qm, &aqp);
      REQUIRE(aqm);
      REQUIRE(aqm->getNumAtoms() == 6);
      REQUIRE(!SubstructMatch(*m, *aqm, match));

      delete aqm;
      aqp.adjustRingCountFlags =
          MolOps::ADJUST_IGNOREDUMMIES;  // no "in ring" restrictions
      aqm = MolOps::adjustQueryProperties(*qm, &aqp);
      REQUIRE(aqm);
      REQUIRE(aqm->getNumAtoms() == 6);
      REQUIRE(SubstructMatch(*m, *aqm, match));

      delete m;
    }

    delete qm;
    delete aqm;
  }

  {  // dummies from SMILES
    std::string smiles = "C1CCC1*";
    ROMol *qm = SmilesToMol(smiles);
    REQUIRE(qm);
    REQUIRE(qm->getNumAtoms() == 5);
    ROMol *aqm = MolOps::adjustQueryProperties(*qm);
    REQUIRE(aqm);
    REQUIRE(aqm->getNumAtoms() == 5);

    smiles = "C1CCC1CC";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    MatchVectType match;
    REQUIRE(!SubstructMatch(*m, *qm, match));
    REQUIRE(SubstructMatch(*m, *aqm, match));

    delete aqm;
    MolOps::AdjustQueryParameters aqp;
    aqp.makeDummiesQueries = false;
    aqm = MolOps::adjustQueryProperties(*qm, &aqp);
    REQUIRE(!SubstructMatch(*m, *aqm, match));

    delete m;
    delete qm;
    delete aqm;
  }
  {  // dummies from SMILES 2
    std::string smiles = "C1CCC1[1*]";
    ROMol *qm = SmilesToMol(smiles);
    REQUIRE(qm);
    REQUIRE(qm->getNumAtoms() == 5);
    ROMol *aqm = MolOps::adjustQueryProperties(*qm);
    REQUIRE(aqm);
    REQUIRE(aqm->getNumAtoms() == 5);
    {
      smiles = "C1CCC1CC";
      ROMol *m = SmilesToMol(smiles);
      REQUIRE(m);
      MatchVectType match;
      REQUIRE(!SubstructMatch(*m, *qm, match));
      REQUIRE(!SubstructMatch(*m, *aqm, match));
      delete m;
    }
    delete qm;
    delete aqm;
  }
  {  // dummies from SMILES 2
    std::string smiles = "C1CCC1[*:1]";
    ROMol *qm = SmilesToMol(smiles);
    qm->getAtomWithIdx(4)->setProp<int>("foo", 2);

    REQUIRE(qm);
    REQUIRE(qm->getNumAtoms() == 5);
    ROMol *aqm = MolOps::adjustQueryProperties(*qm);
    REQUIRE(aqm);
    REQUIRE(aqm->getNumAtoms() == 5);
    REQUIRE(aqm->getAtomWithIdx(4)->getProp<int>("foo") == 2);
    REQUIRE(aqm->getAtomWithIdx(4)->getAtomMapNum() == 1);
    {
      smiles = "C1CCC1CC";
      ROMol *m = SmilesToMol(smiles);
      REQUIRE(m);
      MatchVectType match;
      REQUIRE(!SubstructMatch(*m, *qm, match));
      REQUIRE(SubstructMatch(*m, *aqm, match));
      delete m;
    }
    delete qm;
    delete aqm;
  }
  {  // CTAB
    //  -- only match rgroups
    std::string mb =
        "adjust.mol\n"
        "  ChemDraw06271617272D\n"
        "\n"
        "  7  7  0  0  0  0  0  0  0  0999 V2000\n"
        "   -1.0717    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "   -1.0717   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "   -0.3572   -0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "    0.3572   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "    0.3572    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "   -0.3572    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "    1.0717    0.8250    0.0000 R   0  0  0  0  0  0  0  0  0  1  0  "
        "0\n"
        "  1  2  1  0      \n"
        "  2  3  2  0      \n"
        "  3  4  1  0      \n"
        "  4  5  2  0      \n"
        "  5  6  1  0      \n"
        "  6  1  2  0      \n"
        "  5  7  1  0      \n"
        "M  END\n";
    MolOps::AdjustQueryParameters params;
    params.aromatizeIfPossible = false;
    params.makeDummiesQueries = true;
    params.adjustDegreeFlags =
        (MolOps::ADJUST_IGNOREDUMMIES | MolOps::ADJUST_IGNORECHAINS |
         MolOps::ADJUST_IGNOREMAPPED);

    RWMol *m = MolBlockToMol(mb, false, false);
    MolOps::adjustQueryProperties(*m, &params);
    MatchVectType match;
    ROMol *t = SmilesToMol("c1ccccc1Cl");
    // shouldn't match (aromaticity):
    REQUIRE(!SubstructMatch(*t, *m, match));
    // adjust aromaticity and then it should match:
    params.aromatizeIfPossible = true;
    MolOps::adjustQueryProperties(*m, &params);
    REQUIRE(SubstructMatch(*t, *m, match));

    delete t;
    // shouldn't match (explicit degree)
    t = SmilesToMol("c1ccc(Cl)cc1Cl");
    REQUIRE(!SubstructMatch(*t, *m, match));
    delete m;
    delete t;
  }

  {  // CTAB
    //  -- match non rgroups if mapped
    std::string mb =
        "adjust.mol\n"
        "  ChemDraw06271617272D\n"
        "\n"
        "  7  7  0  0  0  0  0  0  0  0999 V2000\n"
        "   -1.0717    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  "
        "0\n"
        "   -1.0717   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  3  0  "
        "0\n"
        "   -0.3572   -0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  4  0  "
        "0\n"
        "    0.3572   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  5  0  "
        "0\n"
        "    0.3572    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  6  0  "
        "0\n"
        "   -0.3572    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  7  0  "
        "0\n"
        "    1.0717    0.8250    0.0000 R   0  0  0  0  0  0  0  0  0  1  0  "
        "0\n"
        "  1  2  1  0      \n"
        "  2  3  2  0      \n"
        "  3  4  1  0      \n"
        "  4  5  2  0      \n"
        "  5  6  1  0      \n"
        "  6  1  2  0      \n"
        "  5  7  1  0      \n"
        "M  END\n";
    MolOps::AdjustQueryParameters params;
    params.aromatizeIfPossible = true;
    params.makeDummiesQueries = true;
    params.adjustDegreeFlags =
        (MolOps::ADJUST_IGNOREDUMMIES | MolOps::ADJUST_IGNORECHAINS |
         MolOps::ADJUST_IGNOREMAPPED);

    RWMol *m = MolBlockToMol(mb);
    MolOps::adjustQueryProperties(*m, &params);
    MatchVectType match;
    ROMol *t = SmilesToMol("c1ccccc1Cl");
    REQUIRE(SubstructMatch(*t, *m, match));
    delete t;
    // should match (mapped!)
    t = SmilesToMol("c1c(Cl)cc(Cl)cc1Cl");
    REQUIRE(SubstructMatch(*t, *m, match));
    delete m;
    delete t;
  }
  {  // make atoms generic
    std::string smiles = "C1CC1CC";
    ROMol *qm = SmilesToMol(smiles);
    REQUIRE(qm);
    REQUIRE(qm->getNumAtoms() == 5);

    {
      MolOps::AdjustQueryParameters params;
      params.makeAtomsGeneric = true;

      ROMol *aqm = MolOps::adjustQueryProperties(*qm, &params);
      // std::cerr << MolToSmarts(*aqm) << std::endl;
      REQUIRE(aqm);
      REQUIRE(aqm->getNumAtoms() == qm->getNumAtoms());
      {
        MatchVectType match;
        REQUIRE(SubstructMatch(*qm, *aqm, match));
        std::string smiles = "O1CN1NN";
        ROMol *m = SmilesToMol(smiles);
        // std::cerr << MolToSmiles(*m) << std::endl;

        REQUIRE(m);
        REQUIRE(!SubstructMatch(*m, *qm, match));
        REQUIRE(SubstructMatch(*m, *aqm, match));
        delete m;
      }
      delete aqm;
    }
    {
      MolOps::AdjustQueryParameters params;
      params.makeAtomsGeneric = true;
      params.makeAtomsGenericFlags = MolOps::ADJUST_IGNORECHAINS;

      ROMol *aqm = MolOps::adjustQueryProperties(*qm, &params);
      REQUIRE(aqm);
      REQUIRE(aqm->getNumAtoms() == qm->getNumAtoms());
      {
        MatchVectType match;
        REQUIRE(SubstructMatch(*qm, *aqm, match));
        std::string smiles = "O1CN1NN";
        ROMol *m = SmilesToMol(smiles);
        REQUIRE(m);
        REQUIRE(!SubstructMatch(*m, *qm, match));
        REQUIRE(!SubstructMatch(*m, *aqm, match));
        delete m;
        smiles = "O1CN1CC";
        m = SmilesToMol(smiles);
        REQUIRE(m);
        REQUIRE(!SubstructMatch(*m, *qm, match));
        REQUIRE(SubstructMatch(*m, *aqm, match));
        delete m;
      }
      delete aqm;
    }
    {
      MolOps::AdjustQueryParameters params;
      params.makeAtomsGeneric = true;
      params.makeAtomsGenericFlags = MolOps::ADJUST_IGNORERINGS;

      ROMol *aqm = MolOps::adjustQueryProperties(*qm, &params);
      REQUIRE(aqm);
      REQUIRE(aqm->getNumAtoms() == qm->getNumAtoms());
      {
        MatchVectType match;
        REQUIRE(SubstructMatch(*qm, *aqm, match));
        std::string smiles = "O1CN1NN";
        ROMol *m = SmilesToMol(smiles);
        REQUIRE(m);
        REQUIRE(!SubstructMatch(*m, *qm, match));
        REQUIRE(!SubstructMatch(*m, *aqm, match));
        delete m;
        smiles = "C1CC1NN";
        m = SmilesToMol(smiles);
        REQUIRE(m);
        REQUIRE(!SubstructMatch(*m, *qm, match));
        REQUIRE(SubstructMatch(*m, *aqm, match));
        delete m;
      }
      delete aqm;
    }

    delete qm;
  }

  {  // make bonds generic
    std::string smiles = "N1C=C1C=C";
    ROMol *qm = SmilesToMol(smiles);
    REQUIRE(qm);
    REQUIRE(qm->getNumAtoms() == 5);

    {
      MolOps::AdjustQueryParameters params;
      params.makeBondsGeneric = true;

      ROMol *aqm = MolOps::adjustQueryProperties(*qm, &params);
      REQUIRE(aqm);
      REQUIRE(aqm->getNumAtoms() == qm->getNumAtoms());
      {
        MatchVectType match;
        REQUIRE(SubstructMatch(*qm, *aqm, match));
        std::string smiles = "N1=CC1=CC";
        ROMol *m = SmilesToMol(smiles);
        REQUIRE(m);
        REQUIRE(!SubstructMatch(*m, *qm, match));
        REQUIRE(SubstructMatch(*m, *aqm, match));
        delete m;
      }
      delete aqm;
    }
    {
      MolOps::AdjustQueryParameters params;
      params.makeBondsGeneric = true;
      params.makeBondsGenericFlags = MolOps::ADJUST_IGNORECHAINS;

      ROMol *aqm = MolOps::adjustQueryProperties(*qm, &params);
      REQUIRE(aqm);
      REQUIRE(aqm->getNumAtoms() == qm->getNumAtoms());
      {
        MatchVectType match;
        REQUIRE(SubstructMatch(*qm, *aqm, match));
        std::string smiles = "N1=CC1=C=C";
        ROMol *m = SmilesToMol(smiles);
        REQUIRE(m);
        REQUIRE(!SubstructMatch(*m, *qm, match));
        REQUIRE(!SubstructMatch(*m, *aqm, match));
        delete m;
        smiles = "N1=CC1C=C";
        m = SmilesToMol(smiles);
        REQUIRE(m);
        REQUIRE(!SubstructMatch(*m, *qm, match));
        REQUIRE(SubstructMatch(*m, *aqm, match));
        delete m;
      }
      delete aqm;
    }
    {
      MolOps::AdjustQueryParameters params;
      params.makeBondsGeneric = true;
      params.makeBondsGenericFlags = MolOps::ADJUST_IGNORERINGS;

      ROMol *aqm = MolOps::adjustQueryProperties(*qm, &params);
      REQUIRE(aqm);
      REQUIRE(aqm->getNumAtoms() == qm->getNumAtoms());
      {
        MatchVectType match;
        REQUIRE(SubstructMatch(*qm, *aqm, match));
        std::string smiles = "N1=CC1=C=C";
        ROMol *m = SmilesToMol(smiles);
        REQUIRE(m);
        REQUIRE(!SubstructMatch(*m, *qm, match));
        REQUIRE(!SubstructMatch(*m, *aqm, match));
        delete m;
        smiles = "N1C=C1CC#C";
        m = SmilesToMol(smiles);

        REQUIRE(m);
        REQUIRE(!SubstructMatch(*m, *qm, match));
        REQUIRE(SubstructMatch(*m, *aqm, match));
        delete m;
      }
      delete aqm;
    }

    delete qm;
  }

  {  // heavy atom degree
    std::string smiles = "C1CC(*)C1*";
    ROMol *qm = SmartsToMol(smiles);
    REQUIRE(qm);
    REQUIRE(qm->getNumAtoms() == 6);
    MolOps::AdjustQueryParameters params;
    params.adjustDegree = false;
    params.adjustHeavyDegree = true;
    ROMol *aqm = MolOps::adjustQueryProperties(*qm, &params);
    REQUIRE(aqm);
    REQUIRE(aqm->getNumAtoms() == 6);
    {
      smiles = "C1CC(C)C1(C)";
      ROMol *m = SmilesToMol(smiles);
      REQUIRE(m);
      MatchVectType match;
      REQUIRE(SubstructMatch(*m, *qm, match));
      REQUIRE(SubstructMatch(*m, *aqm, match));
      delete m;
    }
    {
      smiles = "C1CC([2H])C1(C)";
      ROMol *m = SmilesToMol(smiles);
      REQUIRE(m);
      MatchVectType match;
      REQUIRE(SubstructMatch(*m, *qm, match));
      REQUIRE(!SubstructMatch(*m, *aqm, match));
      delete m;
    }
    delete qm;
    delete aqm;
  }

  {  // ring-chain membership
    std::string smiles = "CC1CCC1";
    ROMol *qm = SmartsToMol(smiles);
    REQUIRE(qm);
    REQUIRE(qm->getNumAtoms() == 5);
    MolOps::AdjustQueryParameters params;
    params.adjustRingChain = true;
    params.adjustDegree = false;
    ROMol *aqm = MolOps::adjustQueryProperties(*qm, &params);
    REQUIRE(aqm);
    REQUIRE(aqm->getNumAtoms() == 5);
    {
      smiles = "C1CCC12CCC2";
      ROMol *m = SmilesToMol(smiles);
      REQUIRE(m);
      MatchVectType match;
      REQUIRE(SubstructMatch(*m, *qm, match));
      REQUIRE(!SubstructMatch(*m, *aqm, match));
      delete m;
    }
    {
      smiles = "C1CCC1CCC";
      ROMol *m = SmilesToMol(smiles);
      REQUIRE(m);
      MatchVectType match;
      REQUIRE(SubstructMatch(*m, *qm, match));
      REQUIRE(SubstructMatch(*m, *aqm, match));
      delete m;
    }
    delete qm;
    delete aqm;
  }
}

TEST_CASE(
    "Testing github issue 678: failure in AddHs when addCoords is true and coords are all zero") {
  {
    std::string smiles = "CC";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    auto *conf = new Conformer(2);
    m->addConformer(conf);
    MolOps::addHs(*m, false, true);
    REQUIRE(m->getNumAtoms() == 8);
    delete m;
  }

  {  // single connected atom with degenerate coords
    std::string mb =
        "example\n"
        "  Mrv0541 12171503572D          \n"
        "\n"
        "  7  8  0  0  0  0            999 V2000\n"
        "   -3.5063    2.5339    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "   -3.5063    1.7089    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "   -2.6813    2.5339    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "   -2.6813    1.7089    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "   -1.8563    2.5339    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "   -1.8563    1.7089    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "   -1.8563    2.5339    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "  1  2  1  0  0  0  0\n"
        "  1  3  1  0  0  0  0\n"
        "  3  4  1  0  0  0  0\n"
        "  2  4  1  0  0  0  0\n"
        "  3  5  1  0  0  0  0\n"
        "  5  6  1  0  0  0  0\n"
        "  4  6  1  0  0  0  0\n"
        "  5  7  1  0  0  0  0\n"
        "M  END\n";
    RWMol *m = MolBlockToMol(mb);
    REQUIRE(m);
    MolOps::addHs(*m, false, true);
    REQUIRE(m->getNumAtoms() == 19);
    delete m;
  }
  {  // doubly connected atom(s) with degenerate coords
    std::string mb =
        "example\n"
        "  Mrv0541 12171503572D          \n"
        "\n"
        "  7  8  0  0  0  0            999 V2000\n"
        "   -3.5063    2.5339    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "   -3.5063    2.5339    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "   -2.6813    2.5339    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "   -2.6813    1.7089    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "   -1.8563    2.5339    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "   -1.8563    1.7089    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "   -1.2729    3.1173    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "  1  2  1  0  0  0  0\n"
        "  1  3  1  0  0  0  0\n"
        "  3  4  1  0  0  0  0\n"
        "  2  4  1  0  0  0  0\n"
        "  3  5  1  0  0  0  0\n"
        "  5  6  1  0  0  0  0\n"
        "  4  6  1  0  0  0  0\n"
        "  5  7  1  0  0  0  0\n"
        "M  END\n";
    RWMol *m = MolBlockToMol(mb);
    REQUIRE(m);
    MolOps::addHs(*m, false, true);
    REQUIRE(m->getNumAtoms() == 19);
    delete m;
  }

  {  // triply connected atom(s) with degenerate coords
    std::string mb =
        "example\n"
        "  Mrv0541 12171503572D          \n"
        "\n"
        "  7  8  0  0  0  0            999 V2000\n"
        "   -3.5063    2.5339    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "   -3.5063    1.7089    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "   -3.5063    2.5339    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "   -3.5063    2.5339    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "   -1.8563    2.5339    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "   -1.8563    1.7089    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "   -1.2729    3.1173    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "  1  2  1  0  0  0  0\n"
        "  1  3  1  0  0  0  0\n"
        "  3  4  1  0  0  0  0\n"
        "  2  4  1  0  0  0  0\n"
        "  3  5  1  0  0  0  0\n"
        "  5  6  1  0  0  0  0\n"
        "  4  6  1  0  0  0  0\n"
        "  5  7  1  0  0  0  0\n"
        "M  END\n";
    RWMol *m = MolBlockToMol(mb);
    REQUIRE(m);
    MolOps::addHs(*m, false, true);
    REQUIRE(m->getNumAtoms() == 19);
    delete m;
  }
}

TEST_CASE(
    "Testing github issue 717: AddHs cip rank is declared <int> should be unsigned int") {
  // single connected atom with degenerate coords
  std::string mb =
      "mol\n"
      "  Mrv1561 01051606293D\n"
      "\n"
      "  4  3  0  0  0  0            999 V2000\n"
      "   -0.0080   -0.0000   -0.0004 C   0  0  0  0  0  0  0  0  0  0  0  "
      "0\n"
      "    1.5343   -0.0050    0.0032 C   0  0  1  0  0  0  0  0  0  0  0  "
      "0\n"
      "    2.1517   -0.8276    1.4332 Cl  0  0  0  0  0  0  0  0  0  0  0  "
      "0\n"
      "    2.0123   -0.6447   -1.1142 F   0  0  0  0  0  0  0  0  0  0  0  "
      "0\n"
      "  2  3  1  0  0  0  0\n"
      "  2  4  1  0  0  0  0\n"
      "  2  1  1  0  0  0  0\n"
      "M  END\n";
  RWMol *m = MolBlockToMol(mb);
  REQUIRE(m);
  MolOps::assignChiralTypesFrom3D(*m);
  MolOps::assignStereochemistry(*m, true, true);
  MolOps::addHs(*m, false, true);
  REQUIRE(m->getNumAtoms() == 8);
  delete m;
}

TEST_CASE("Testing findPotentialStereoBonds") {
  {  // starting point: full sanitization
    auto m1 = "BrC(=NN=c1nn[nH][nH]1)c1ccncc1"_smiles;
    REQUIRE(m1);
    std::string smiles =
        "Br/C(=N\\N=c1/nn[nH][nH]1)c1ccncc1";  // possible problem reported by
                                               // Steve Roughley
    ROMol *m = SmilesToMol(smiles);
    // m->updatePropertyCache(false);
    // m->debugMol(std::cerr);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 15);
    REQUIRE(m->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(1)->getStereoAtoms().size() == 2);
    REQUIRE(m->getBondWithIdx(3)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(3)->getStereoAtoms().size() == 2);
    delete m;

    // partial sanitization:
    m = SmilesToMol(smiles, false, false);
    REQUIRE(m);
    m->updatePropertyCache(true);
    MolOps::findSSSR(*m);
    MolOps::findPotentialStereoBonds(*m, false);
    REQUIRE(m->getNumAtoms() == 15);
    REQUIRE(m->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(1)->getStereoAtoms().size() == 2);
    REQUIRE(m->getBondWithIdx(3)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(3)->getStereoAtoms().size() == 2);
    delete m;
  }

  // this next block is for github1230: FindPotentialStereoBonds() failure
  {  // simple
    std::string smiles = "CC=CC";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 4);
    REQUIRE(m->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(1)->getStereo() == Bond::STEREONONE);

    MolOps::findPotentialStereoBonds(*m, true);
    REQUIRE(m->getBondWithIdx(1)->getStereo() == Bond::STEREOANY);

    delete m;
  }
  {  // simple2
    std::string smiles = "CC=C(C)C";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 5);
    REQUIRE(m->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(1)->getStereo() == Bond::STEREONONE);

    MolOps::findPotentialStereoBonds(*m, true);
    REQUIRE(m->getBondWithIdx(1)->getStereo() == Bond::STEREONONE);

    delete m;
  }
  {  // the real problem
    std::string smiles = "CC/C=C/C(\\C=C/CC)=C(CC)CO";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 14);
    REQUIRE(m->getBondWithIdx(2)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(2)->getStereo() == Bond::STEREOE);
    REQUIRE(m->getBondWithIdx(5)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(5)->getStereo() == Bond::STEREOZ);
    REQUIRE(m->getBondWithIdx(8)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(8)->getStereo() == Bond::STEREONONE);
    MolOps::findPotentialStereoBonds(*m, true);
    REQUIRE(m->getBondWithIdx(2)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(2)->getStereo() == Bond::STEREOANY);
    REQUIRE(m->getBondWithIdx(5)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(5)->getStereo() == Bond::STEREOANY);
    REQUIRE(m->getBondWithIdx(8)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(8)->getStereo() == Bond::STEREOANY);

    delete m;
  }
  {  // repeat the real problem, but set the cleanIt argument to false
    std::string smiles = "CC/C=C/C(\\C=C/CC)=C(CC)CO";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 14);
    REQUIRE(m->getBondWithIdx(2)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(2)->getStereo() == Bond::STEREOE);
    REQUIRE(m->getBondWithIdx(5)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(5)->getStereo() == Bond::STEREOZ);
    REQUIRE(m->getBondWithIdx(8)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(8)->getStereo() == Bond::STEREONONE);
    MolOps::findPotentialStereoBonds(*m, false);
    REQUIRE(m->getBondWithIdx(2)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(2)->getStereo() == Bond::STEREOE);
    REQUIRE(m->getBondWithIdx(5)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(5)->getStereo() == Bond::STEREOZ);
    REQUIRE(m->getBondWithIdx(8)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(8)->getStereo() == Bond::STEREOANY);

    delete m;
  }

  {  // just do document that we still don't do this one, which is much harder
    std::string smiles = "CC/C=C/C(/C=C/CC)=C(CC)CO";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 14);
    REQUIRE(m->getBondWithIdx(2)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(2)->getStereo() == Bond::STEREOE);
    REQUIRE(m->getBondWithIdx(5)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(5)->getStereo() == Bond::STEREOE);
    REQUIRE(m->getBondWithIdx(8)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(8)->getStereo() == Bond::STEREONONE);
    MolOps::findPotentialStereoBonds(*m, true);
    REQUIRE(m->getBondWithIdx(2)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(2)->getStereo() == Bond::STEREOANY);
    REQUIRE(m->getBondWithIdx(5)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(5)->getStereo() == Bond::STEREOANY);
    REQUIRE(m->getBondWithIdx(8)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(8)->getStereo() == Bond::STEREONONE);

    delete m;
  }
}

TEST_CASE("Testing Bond::setStereo(Bond::STEREOCIS / Bond::STEREOTRANS)") {
  // tests to make sure neighboring bond stereo is handled properly
  {
    const char *smiles[] = {"CC=CC",         "CC=C/C=C/C",    "CC=C/C=C\\C",
                            "CC=C\\C=C/C",   "CC=C\\C=C\\C",  "C(C)=CC",
                            "C(C)=C/C=C/C",  "C(C)=C/C=C\\C", "C(C)=C\\C=C/C",
                            "C(C)=C\\C=C\\C"};
    const Bond::BondStereo stereos[] = {Bond::STEREOCIS, Bond::STEREOTRANS};
    const Bond::BondStereo ezstros[] = {Bond::STEREOZ, Bond::STEREOE};

    for (auto &smile : smiles) {
      ROMol *m = SmilesToMol(smile);
      MolOps::findPotentialStereoBonds(*m);
      Bond *bond = m->getBondWithIdx(1);

      for (size_t j = 0; j < 2; ++j) {
        Bond::BondStereo desired_stereo = stereos[j];
        bond->setStereo(desired_stereo);

        bool doIsomericSmiles = true;
        bool doKekule = false;
        int rootedAtAtom = -1;
        bool canonical = false;
        std::string isosmi = MolToSmiles(*m, doIsomericSmiles, doKekule,
                                         rootedAtAtom, canonical);

        ROMol *isomol = SmilesToMol(isosmi);
        Bond *isobond = isomol->getBondWithIdx(1);
        const Bond::BondStereo expected_ez_stereo = ezstros[j];
        REQUIRE(isobond->getStereo() == expected_ez_stereo);

        std::string round_trip_isosmi = MolToSmiles(
            *m, doIsomericSmiles, doKekule, rootedAtAtom, canonical);
        REQUIRE(isosmi == round_trip_isosmi);

        BOOST_LOG(rdInfoLog) << isosmi << " == " << round_trip_isosmi << " "
                             << desired_stereo << std::endl;

        delete isomol;
      }
      delete m;
    }
  }

  // tests enumerating all possible smiles with halogens still yield
  // the same isomeric canonical smiles strings.
  {
    const char *smiles[] = {"ClC=CF", "FC=CCl", "C(Cl)=CF", "C(F)=CCl"};
    const Bond::BondStereo stereos[] = {Bond::STEREOCIS, Bond::STEREOTRANS};

    for (auto desired_stereo : stereos) {
      std::string refSmiles;
      for (auto &smile : smiles) {
        ROMol *m = SmilesToMol(smile);
        MolOps::findPotentialStereoBonds(*m);
        REQUIRE(m->getNumAtoms() == 4);

        Bond *doubleBond = m->getBondWithIdx(1);
        doubleBond->setStereo(desired_stereo);

        bool doIsomericSmiles = true;
        std::string isocansmi = MolToSmiles(*m, doIsomericSmiles);

        if (refSmiles.empty()) {
          refSmiles = isocansmi;
        }
        BOOST_LOG(rdInfoLog) << refSmiles << " == " << isocansmi << " "
                             << desired_stereo << std::endl;
        REQUIRE(refSmiles == isocansmi);

        delete m;
      }
    }
  }
}

TEST_CASE("Testing Bond::setStereoAtoms(...)") {
  // tests to make sure setStereoAtoms works as expected

  bool doIsomericSmiles = true;
  std::string unspec_smiles = "FC(Cl)=C(Br)I";

  ROMol *m = SmilesToMol(unspec_smiles);

  Bond *doubleBond = m->getBondWithIdx(2);
  REQUIRE(doubleBond->getBondType() == 2);

  doubleBond->setStereoAtoms(0, 4);
  doubleBond->setStereo(Bond::STEREOCIS);
  BOOST_LOG(rdInfoLog) << MolToSmiles(*m, doIsomericSmiles) << std::endl;
  REQUIRE(MolToSmiles(*m, doIsomericSmiles) == "F/C(Cl)=C(\\Br)I");
  // this should be the same as the previous
  doubleBond->setStereoAtoms(0, 5);
  doubleBond->setStereo(Bond::STEREOTRANS);
  BOOST_LOG(rdInfoLog) << MolToSmiles(*m, doIsomericSmiles) << std::endl;
  REQUIRE(MolToSmiles(*m, doIsomericSmiles) == "F/C(Cl)=C(\\Br)I");

  doubleBond->setStereoAtoms(0, 4);
  doubleBond->setStereo(Bond::STEREOTRANS);
  BOOST_LOG(rdInfoLog) << MolToSmiles(*m, doIsomericSmiles) << std::endl;
  REQUIRE(MolToSmiles(*m, doIsomericSmiles) == "F/C(Cl)=C(/Br)I");
  // this should be the same as the previous
  doubleBond->setStereoAtoms(0, 5);
  doubleBond->setStereo(Bond::STEREOCIS);
  BOOST_LOG(rdInfoLog) << MolToSmiles(*m, doIsomericSmiles) << std::endl;
  REQUIRE(MolToSmiles(*m, doIsomericSmiles) == "F/C(Cl)=C(/Br)I");

  doubleBond->setStereoAtoms(3, 4);
  doubleBond->setStereo(Bond::STEREOTRANS);
  BOOST_LOG(rdInfoLog) << MolToSmiles(*m, doIsomericSmiles) << std::endl;
  REQUIRE(MolToSmiles(*m, doIsomericSmiles) == "F/C(Cl)=C(\\Br)I");
  // this should be the same as the previous
  doubleBond->setStereoAtoms(3, 5);
  doubleBond->setStereo(Bond::STEREOCIS);
  BOOST_LOG(rdInfoLog) << MolToSmiles(*m, doIsomericSmiles) << std::endl;
  REQUIRE(MolToSmiles(*m, doIsomericSmiles) == "F/C(Cl)=C(\\Br)I");

  doubleBond->setStereoAtoms(3, 4);
  doubleBond->setStereo(Bond::STEREOCIS);
  BOOST_LOG(rdInfoLog) << MolToSmiles(*m, doIsomericSmiles) << std::endl;
  REQUIRE(MolToSmiles(*m, doIsomericSmiles) == "F/C(Cl)=C(/Br)I");
  // this should be the same as the previous
  doubleBond->setStereoAtoms(3, 5);
  doubleBond->setStereo(Bond::STEREOTRANS);
  BOOST_LOG(rdInfoLog) << MolToSmiles(*m, doIsomericSmiles) << std::endl;
  REQUIRE(MolToSmiles(*m, doIsomericSmiles) == "F/C(Cl)=C(/Br)I");
  delete m;
}

TEST_CASE("Testing github #754 : loss of double bond geometry with removeHs") {
  {  // starting point: full sanitization
    std::string smiles =
        "[H]C([H])([H])/C([H])=C(/[H])C([H])([H])[H]";  // possible problem
                                                        // reported by
                                                        // Steve Roughley
    RWMol *m = SmilesToMol(smiles, false, false);
    REQUIRE(m);
    MolOps::sanitizeMol(*m);
    MolOps::assignStereochemistry(*m, true, true);
    REQUIRE(m->getNumAtoms() == 12);
    REQUIRE(m->getBondWithIdx(5)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(5)->getStereo() == Bond::STEREOZ);
    delete m;

    m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 4);
    REQUIRE(m->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);
    delete m;
  }
  {  // another basic test
    std::string smiles = "[H]/C(C)=C/C";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 4);
    REQUIRE(m->getBondBetweenAtoms(0, 2)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondBetweenAtoms(0, 2)->getStereo() == Bond::STEREOZ);
    delete m;
  }
  {  // H following the C:
    std::string smiles = "CC(\\[H])=C/C";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 4);
    REQUIRE(m->getBondBetweenAtoms(1, 2)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondBetweenAtoms(1, 2)->getStereo() == Bond::STEREOZ);
    delete m;
  }
  {  // bond dir already set :
    std::string smiles = "[H]/C(/C)=C\\C";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 4);
    REQUIRE(m->getBondBetweenAtoms(0, 2)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondBetweenAtoms(0, 2)->getStereo() == Bond::STEREOE);
    delete m;
  }

  {  // chained bonds :
    std::string smiles = "[H]/C(C=C/C)=C\\C";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 6);
    REQUIRE(m->getBondBetweenAtoms(0, 4)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondBetweenAtoms(0, 4)->getStereo() == Bond::STEREOE);
    REQUIRE(m->getBondBetweenAtoms(1, 2)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondBetweenAtoms(1, 2)->getStereo() == Bond::STEREOE);
    delete m;
  }
}

TEST_CASE("Testing github #805 : Pre-condition Violation: bad bond type") {
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    ROMol *m = MolFileToMol(pathName + "pubchem_87396055.sdf");

    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 20);
    REQUIRE(m->getBondBetweenAtoms(2, 6)->getBondType() == Bond::SINGLE);
    REQUIRE(m->getAtomWithIdx(2)->getFormalCharge() == 1);
    REQUIRE(m->getAtomWithIdx(6)->getFormalCharge() == -1);
    REQUIRE(m->getBondBetweenAtoms(2, 9)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondBetweenAtoms(2, 9)->getStereo() != Bond::STEREONONE);
    REQUIRE(m->getBondBetweenAtoms(3, 10)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondBetweenAtoms(3, 10)->getStereo() != Bond::STEREONONE);
    std::string smi = MolToSmiles(*m, true);
    REQUIRE(smi == "CCO/[P+]([O-])=C1\\CSC(c2cccs2)\\C1=[P+](\\[O-])OCC");
    delete m;
  }
  {
    std::string smi = "O=P(/O)=C/C";
    ROMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 5);
    REQUIRE(m->getBondBetweenAtoms(0, 1)->getBondType() == Bond::SINGLE);
    REQUIRE(m->getAtomWithIdx(1)->getFormalCharge() == 1);
    REQUIRE(m->getAtomWithIdx(0)->getFormalCharge() == -1);
    REQUIRE(m->getBondBetweenAtoms(1, 3)->getBondType() == Bond::DOUBLE);
    REQUIRE(m->getBondBetweenAtoms(1, 3)->getStereo() != Bond::STEREONONE);
    smi = MolToSmiles(*m, true);
    REQUIRE(smi == "C/C=[P+](/[O-])O");
    delete m;
  }
}

TEST_CASE(
    "Testing github #518 : Rings containing all dummy atoms with single bonds are flagged as aromatic") {
  {
    std::string smi = "*-1-*-*-*-1";
    ROMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 4);
    REQUIRE(!m->getBondWithIdx(0)->getIsAromatic());
    REQUIRE(m->getBondWithIdx(0)->getBondType() == Bond::SINGLE);

    REQUIRE(!m->getAtomWithIdx(0)->getIsAromatic());
    delete m;
  }
  {  // in this case we leave it aromatic since it's all dummies
    std::string smi = "*:1:*:*:*:1";
    ROMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 4);
    REQUIRE(!m->getBondWithIdx(0)->getIsAromatic());
    REQUIRE(m->getBondWithIdx(0)->getBondType() == Bond::AROMATIC);
    REQUIRE(!m->getAtomWithIdx(0)->getIsAromatic());
    delete m;
  }
  {
    std::string smi = "*-1-*-C-*-*-*-1";
    ROMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 6);
    REQUIRE(!m->getBondWithIdx(0)->getIsAromatic());
    REQUIRE(m->getBondWithIdx(0)->getBondType() == Bond::SINGLE);

    REQUIRE(!m->getAtomWithIdx(0)->getIsAromatic());
    delete m;
  }
  {
    std::string smi = "C1=CC=*2*(=C1)*1=CC=CC=*1*1=CC=CC=*21";
    ROMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 18);
    REQUIRE(!m->getBondBetweenAtoms(4, 6)->getIsAromatic());
    REQUIRE(m->getBondBetweenAtoms(4, 6)->getBondType() == Bond::SINGLE);
    REQUIRE(!m->getBondBetweenAtoms(11, 12)->getIsAromatic());
    REQUIRE(m->getBondBetweenAtoms(11, 12)->getBondType() == Bond::SINGLE);
    REQUIRE(!m->getBondBetweenAtoms(3, 17)->getIsAromatic());
    REQUIRE(m->getBondBetweenAtoms(3, 17)->getBondType() == Bond::SINGLE);

    delete m;
  }
}

TEST_CASE("Testing simple aromaticity") {
  {
    std::string smiles = "c1ccccc1";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic() == true);
    REQUIRE(m->getAtomWithIdx(0)->getIsAromatic() == true);
    MolOps::Kekulize(*m, true);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic() == false);
    REQUIRE(m->getAtomWithIdx(0)->getIsAromatic() == false);
    MolOps::setAromaticity(*m, MolOps::AROMATICITY_SIMPLE);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic() == true);
    REQUIRE(m->getAtomWithIdx(0)->getIsAromatic() == true);
    delete m;
  }
  {
    std::string smiles = "c1[nH]ccc1";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic() == true);
    REQUIRE(m->getAtomWithIdx(0)->getIsAromatic() == true);
    MolOps::Kekulize(*m, true);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic() == false);
    REQUIRE(m->getAtomWithIdx(0)->getIsAromatic() == false);
    MolOps::setAromaticity(*m, MolOps::AROMATICITY_SIMPLE);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic() == true);
    REQUIRE(m->getAtomWithIdx(0)->getIsAromatic() == true);
    delete m;
  }
  {  // ring size constraints
    std::string smiles = "c1cccoocc1";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic() == true);
    REQUIRE(m->getAtomWithIdx(0)->getIsAromatic() == true);
    MolOps::Kekulize(*m, true);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic() == false);
    REQUIRE(m->getAtomWithIdx(0)->getIsAromatic() == false);
    MolOps::setAromaticity(*m, MolOps::AROMATICITY_SIMPLE);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic() == false);
    REQUIRE(m->getAtomWithIdx(0)->getIsAromatic() == false);
    delete m;
  }
  {  // ring size constraints
    std::string smiles = "c1coo1";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic() == true);
    REQUIRE(m->getAtomWithIdx(0)->getIsAromatic() == true);
    MolOps::Kekulize(*m, true);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic() == false);
    REQUIRE(m->getAtomWithIdx(0)->getIsAromatic() == false);
    MolOps::setAromaticity(*m, MolOps::AROMATICITY_SIMPLE);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic() == false);
    REQUIRE(m->getAtomWithIdx(0)->getIsAromatic() == false);
    delete m;
  }
  {  // fused rings are not considered
    std::string smiles = "C1=CC2=CC=CC=CC2=C1";  // azulene
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic() == true);
    REQUIRE(m->getAtomWithIdx(0)->getIsAromatic() == true);
    MolOps::Kekulize(*m, true);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic() == false);
    REQUIRE(m->getAtomWithIdx(0)->getIsAromatic() == false);
    MolOps::setAromaticity(*m, MolOps::AROMATICITY_SIMPLE);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic() == false);
    REQUIRE(m->getAtomWithIdx(0)->getIsAromatic() == false);
    delete m;
  }
}

TEST_CASE("Testing MMFF94 aromaticity") {
  // test one known difference between RDKit and MMFF94 aromaticity models:
  // the latter does not recognize azulene as aromatic

  std::string smiles = "C1=CC=C2C=CC=C2C=C1";
  RWMol *m = SmilesToMol(smiles);
  MolOps::Kekulize(*m, true);

  MolOps::setAromaticity(*m, MolOps::AROMATICITY_RDKIT);
  int arombondcount = 0;
  for (auto b : m->bonds()) {
    if (b->getIsAromatic()) arombondcount++;
  }
  // all bonds, except the fused one, should be aromatic
  REQUIRE(arombondcount == 10);
  REQUIRE(m->getBondBetweenAtoms(3, 7)->getIsAromatic() == false);

  MolOps::setAromaticity(*m, MolOps::AROMATICITY_MMFF94);
  arombondcount = 0;
  for (auto b : m->bonds()) {
    if (b->getIsAromatic()) arombondcount++;
  }
  // no aromatics here
  REQUIRE(arombondcount == 0);
  delete m;
}

//! really dumb aromaticity: any conjugated ring bond is aromatic
int customAromaticity(RWMol &m) {
  m.updatePropertyCache();
  MolOps::setConjugation(m);
  MolOps::fastFindRings(m);
  int res = 0;
  for (ROMol::BondIterator bIt = m.beginBonds(); bIt != m.endBonds(); ++bIt) {
    if ((*bIt)->getIsConjugated() && queryIsBondInRing(*bIt)) {
      (*bIt)->setIsAromatic(true);
      (*bIt)->getBeginAtom()->setIsAromatic(true);
      (*bIt)->getEndAtom()->setIsAromatic(true);
      ++res;
    }
  }
  return res;
}

TEST_CASE("Testing custom aromaticity") {
  {
    std::string smiles = "C1=CC=CC=C1";
    RWMol *m = SmilesToMol(smiles, 0, false);
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic() == false);
    REQUIRE(m->getAtomWithIdx(0)->getIsAromatic() == false);
    MolOps::setAromaticity(*m, MolOps::AROMATICITY_CUSTOM, customAromaticity);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic() == true);
    REQUIRE(m->getAtomWithIdx(0)->getIsAromatic() == true);
    delete m;
  }
  {
    std::string smiles = "C1CC=CC=C1";
    RWMol *m = SmilesToMol(smiles, 0, false);
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic() == false);
    REQUIRE(m->getAtomWithIdx(0)->getIsAromatic() == false);
    MolOps::setAromaticity(*m, MolOps::AROMATICITY_CUSTOM, customAromaticity);
    REQUIRE(m->getBondWithIdx(0)->getIsAromatic() == false);
    REQUIRE(m->getAtomWithIdx(0)->getIsAromatic() == false);
    REQUIRE(m->getBondWithIdx(2)->getIsAromatic() == true);
    REQUIRE(m->getAtomWithIdx(2)->getIsAromatic() == true);
    delete m;
  }
}

TEST_CASE(
    "Testing github issue #1730: setAromaticity() should work even if there are aromatic atoms present") {
  std::string smiles = "C1=CC=CC=C1-c2ccccc2";
  RWMol *m = SmilesToMol(smiles, 0, false);
  m->updatePropertyCache();
  REQUIRE(m);
  REQUIRE(m->getBondWithIdx(0)->getIsAromatic() == false);
  REQUIRE(m->getAtomWithIdx(0)->getIsAromatic() == false);
  REQUIRE(m->getBondWithIdx(6)->getIsAromatic() == true);
  REQUIRE(m->getAtomWithIdx(6)->getIsAromatic() == true);
  MolOps::setAromaticity(*m);
  REQUIRE(m->getBondWithIdx(0)->getIsAromatic() == true);
  REQUIRE(m->getAtomWithIdx(0)->getIsAromatic() == true);
  delete m;
}

TEST_CASE("Testing error reporting for kekulization") {
  std::stringstream sstrm;
  rdErrorLog->SetTee(sstrm);
  {
    sstrm.str("");
    std::string smi = "c1ccccc1";
    ROMol *m = SmilesToMol(smi);
    REQUIRE(m);
    REQUIRE(sstrm.str() == "");
    delete m;
  }
  {
    sstrm.str("");
    std::string smi = "c1cccc1";
    ROMol *m;
    try {
      m = SmilesToMol(smi);
    } catch (MolSanitizeException &) {
      m = nullptr;
    }
    REQUIRE(m == nullptr);
    REQUIRE(sstrm.str().find("0 1 2 3 4") != std::string::npos);
    delete m;
  }
  {
    sstrm.str("");
    std::string smi = "c1ccccc1.c1cccc1";
    ROMol *m;
    try {
      m = SmilesToMol(smi);
    } catch (MolSanitizeException &) {
      m = nullptr;
    }
    REQUIRE(m == nullptr);
    REQUIRE(sstrm.str().find("6 7 8 9 10") != std::string::npos);
    delete m;
  }
  {
    sstrm.str("");
    std::string smi = "c1cccc1.c1cccc1";
    ROMol *m;
    try {
      m = SmilesToMol(smi);
    } catch (MolSanitizeException &) {
      m = nullptr;
    }
    REQUIRE(m == nullptr);
    REQUIRE(sstrm.str().find("0 1 2 3 4") != std::string::npos);
    delete m;
  }
  rdErrorLog->ClearTee();
}

TEST_CASE(
    "Testing github issue #868: inappropriate warning from MergeQueryHs") {
  std::stringstream sstrm;
  rdWarningLog->SetTee(sstrm);
  {
    sstrm.str("");

    std::string sma = "[SX3](=O)[O-,#1]";
    RWMol *m = SmartsToMol(sma);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 3);
    MolOps::mergeQueryHs(*m);
    REQUIRE(
        sstrm.str().find(
            "merging explicit H queries involved in ORs is not supported") !=
        std::string::npos);
    REQUIRE(sstrm.str().find("This query will not be merged") !=
            std::string::npos);
    delete m;
  }
  {
    sstrm.str("");

    std::string sma = "[SX3](=O)[O-,H1]";
    RWMol *m = SmartsToMol(sma);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 3);
    MolOps::mergeQueryHs(*m);
    REQUIRE(sstrm.str().find("merging explicit H queries involved in "
                             "ORs is not supported") == std::string::npos);
    REQUIRE(sstrm.str().find("This query will not be merged") ==
            std::string::npos);
    delete m;
  }
  {
    sstrm.str("");

    std::string sma = "[SX3](=O)[O-,H]";
    RWMol *m = SmartsToMol(sma);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 3);
    MolOps::mergeQueryHs(*m);
    REQUIRE(sstrm.str().find("merging explicit H queries involved in "
                             "ORs is not supported") == std::string::npos);
    REQUIRE(sstrm.str().find("This query will not be merged") ==
            std::string::npos);
    delete m;
  }

  // test atom type query merging
  for (int aromatic = 0; aromatic < 2; ++aromatic) {
    sstrm.str("");
    REQUIRE(sstrm.str() == "");
    RWMol m;
    QueryAtom *qa = new QueryAtom();
    qa->setQuery(makeAtomTypeQuery(1, aromatic));
    qa->expandQuery(makeAtomNumQuery(6),
                    Queries::CompositeQueryType::COMPOSITE_OR);
    m.addAtom(qa, true, true);
    MolOps::mergeQueryHs(m);
    REQUIRE(sstrm.str().find("merging explicit H queries involved in "
                             "ORs is not supported") != std::string::npos);
    REQUIRE(sstrm.str().find("This query will not be merged") !=
            std::string::npos);
  }

  {
    sstrm.str("");
    REQUIRE(sstrm.str() == "");
    // github 7687 - merge with more than one option in or
    std::string sma = "[#6]-[#1,#6,#7]";
    RWMol *m = SmartsToMol(sma);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 2);
    MolOps::mergeQueryHs(*m);
    std::cerr << sstrm.str() << std::endl;
    REQUIRE(sstrm.str().find("merging explicit H queries involved in "
                             "ORs is not supported") != std::string::npos);
    REQUIRE(sstrm.str().find("This query will not be merged") !=
            std::string::npos);
    delete m;
  }

  rdWarningLog->ClearTee();
}

TEST_CASE(
    "Testing github issue 908: AddHs() using 3D coordinates with 2D conformations") {
  std::string mb =
      "\n     RDKit          2D\n\n  4  3  0  0  0  0  0  0  0  0999 "
      "V2000\n "
      "  -0.0000   -1.5000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  "
      "0\n   -0.0000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  "
      "0 "
      " 0\n    1.2990    0.7500    0.0000 F   0  0  0  0  0  0  0  0  0  0 "
      " "
      "0  0\n   -1.2990    0.7500    0.0000 Cl  0  0  0  0  0  0  0  0  0  "
      "0 "
      " 0  0\n  2  1  1  1\n  2  3  1  0\n  2  4  1  0\nM  END\n";
  RWMol *m = MolBlockToMol(mb);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 4);
  MolOps::addHs(*m, false, true);
  REQUIRE(m->getNumAtoms() == 5);
  REQUIRE(feq(m->getConformer().getAtomPos(4).z, 0.0));
  delete m;
}

TEST_CASE(
    "Testing github issue 962: Kekulization issues post successful smiles parsing") {
  {
    std::string smi = "C2*c1ccccc1C2";
    RWMol *m = SmilesToMol(smi, 0, false);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 9);
    REQUIRE(m->getBondBetweenAtoms(0, 1));
    REQUIRE(m->getBondBetweenAtoms(2, 1));
    m->updatePropertyCache();
    MolOps::Kekulize(*m);
    REQUIRE(m->getBondBetweenAtoms(0, 1)->getBondType() == Bond::SINGLE);
    REQUIRE(m->getBondBetweenAtoms(2, 1)->getBondType() == Bond::SINGLE);

    delete m;
  }
  {  // this one did not cause problems before, but verify!
    std::string smi = "*2Cc1ccccc1C2";
    RWMol *m = SmilesToMol(smi, 0, false);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 9);
    REQUIRE(m->getBondBetweenAtoms(0, 1));
    REQUIRE(m->getBondBetweenAtoms(0, 8));
    m->updatePropertyCache();
    MolOps::Kekulize(*m);
    REQUIRE(m->getBondBetweenAtoms(0, 1)->getBondType() == Bond::SINGLE);
    REQUIRE(m->getBondBetweenAtoms(0, 8)->getBondType() == Bond::SINGLE);

    delete m;
  }
}

TEST_CASE(
    "Testing github issue 1021: AssignStereochemistry() giving incorrect results after  FastFindRings()") {
  std::string smi = "C[C@H]1CC2CCCC(C1)[C@H]2N";
  RWMol *m = SmilesToMol(smi);
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 11);
  REQUIRE(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  REQUIRE(m->getAtomWithIdx(9)->getChiralTag() != Atom::CHI_UNSPECIFIED);

  m->clearComputedProps();
  bool cleanit = true, force = true;
  MolOps::assignStereochemistry(*m, cleanit, force);
  REQUIRE(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  REQUIRE(m->getAtomWithIdx(9)->getChiralTag() != Atom::CHI_UNSPECIFIED);

  m->clearComputedProps();
  MolOps::fastFindRings(*m);
  MolOps::assignStereochemistry(*m, cleanit, force);
  REQUIRE(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  REQUIRE(m->getAtomWithIdx(9)->getChiralTag() != Atom::CHI_UNSPECIFIED);

  delete m;
}

TEST_CASE(
    "Testing github issue 607: AssignAtomChiralTagsFromStructure() not recognizing chiral S") {
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    ROMol *m = MolFileToMol(pathName + "1a9u.zwitterion.sdf");
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 27);
    MolOps::assignChiralTypesFrom3D(*m);

    REQUIRE(m->getAtomWithIdx(26)->getAtomicNum() == 16);
    REQUIRE(m->getAtomWithIdx(26)->getChiralTag() != Atom::CHI_UNSPECIFIED);

    delete m;
  }
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    ROMol *m = MolFileToMol(pathName + "1a9u.sdf");
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 27);
    MolOps::assignChiralTypesFrom3D(*m);

    REQUIRE(m->getAtomWithIdx(26)->getAtomicNum() == 16);
    REQUIRE(m->getAtomWithIdx(26)->getChiralTag() != Atom::CHI_UNSPECIFIED);

    delete m;
  }
  {  // convert S -> Se and test again
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    ROMol *m = MolFileToMol(pathName + "1a9u.zwitterion.sdf");
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 27);
    m->getAtomWithIdx(26)->setAtomicNum(34);
    MolOps::assignChiralTypesFrom3D(*m);
    REQUIRE(m->getAtomWithIdx(26)->getChiralTag() != Atom::CHI_UNSPECIFIED);

    delete m;
  }
  {  // convert S -> Se and test again
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    ROMol *m = MolFileToMol(pathName + "1a9u.sdf");
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 27);
    m->getAtomWithIdx(26)->setAtomicNum(34);
    MolOps::assignChiralTypesFrom3D(*m);
    REQUIRE(m->getAtomWithIdx(26)->getChiralTag() != Atom::CHI_UNSPECIFIED);

    delete m;
  }
}

TEST_CASE("Testing github issue 1204: Support tetravalent and hexavalent Te") {
  {
    std::string smiles = "F[Te](F)(F)(F)(F)F";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    delete m;
  }
  {
    std::string smiles = "F[Te](F)(F)(F)";
    RWMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    delete m;
  }
}

TEST_CASE(
    "Testing github issue 1478: Aromatic rings composed solely of dummy atoms should not be kekulized") {
  {  // basics
    std::string smiles = "*:1:*:*:*:*:*:1";
    RWMol *m = SmilesToMol(smiles, false);
    REQUIRE(m);
    m->updatePropertyCache();
    MolOps::Kekulize(*m);
    for (unsigned int i = 0; i < m->getNumBonds(); ++i) {
      REQUIRE(m->getBondWithIdx(i)->getBondType() == Bond::AROMATIC);
    }
    delete m;
  }

  {  // fused rings where one is kekulized
    std::string smiles = "*:1:*:*:*:*:2:*:1cccc2";
    RWMol *m = SmilesToMol(smiles, false);
    REQUIRE(m);
    m->updatePropertyCache();
    MolOps::Kekulize(*m);
    REQUIRE(m->getBondBetweenAtoms(0, 1)->getBondType() == Bond::AROMATIC);
    REQUIRE(m->getBondBetweenAtoms(6, 7)->getBondType() != Bond::AROMATIC);
    delete m;
  }
}

TEST_CASE(
    "Testing github issue 1439: RemoveHs() removes H atom attached to dummy if it came from AddHs()") {
  std::string smiles = "F";
  RWMol *m = SmilesToMol(smiles);
  REQUIRE(m);
  MolOps::addHs(*m);
  REQUIRE(m->getNumAtoms() == 2);
  m->getAtomWithIdx(0)->setAtomicNum(0);
  MolOps::removeHs(*m);
  REQUIRE(m->getNumAtoms() == 2);
  delete m;
}

TEST_CASE(
    "Testing github issue 1281: RDKit gets stuck on PubChem CID 102128817") {
  {  // basics
    std::string smiles =
        "COC1=CC=C(C=C1)C2C3=C(C=CC4=CC=CC=C43)OC5=CC6=C(C=C5)C7=NC8=C9C=CC1="
        "CC9=C(N8)N=C3C4=C5C=CC(=C4)OC4=C(C(C8=C(C=CC9=CC=CC=C98)OC8=CC9=C(C="
        "C8)C8=NC9=NC9=C%10C=C(C=CC%10=C(N9)N=C9C%10=C(C=C(C=C%10)OC%10=C2C2="
        "CC=CC=C2C=C%10)C(=N9)NC2=NC(=N8)C8=C2C=C(C=C8)OC2=C(C(C8=C(C=CC9=CC="
        "CC=C98)OC8=CC9=C(C=C8)C(=NC5=N3)N=C9NC6=N7)C3=CC=C(C=C3)OC)C3=CC=CC="
        "C3C=C2)OC2=C(C(C3=C(O1)C=CC1=CC=CC=C13)C1=CC=C(C=C1)OC)C1=CC=CC=C1C="
        "C2)C1=CC=C(C=C1)OC)C1=CC=CC=C1C=C4";
    {
      RWMol *m = SmilesToMol(smiles, 0, false);
      REQUIRE(m);
      REQUIRE(m->getNumAtoms() == 204);
      REQUIRE(m->getNumBonds() == 244);
      bool ok = false;
      try {
        MolOps::findSSSR(*m);
      } catch (const ValueErrorException &) {
        ok = true;
      }
      REQUIRE(ok);
      delete m;
    }
    {
      bool ok = false;
      try {
        RWMol *m = SmilesToMol(smiles);
        // Can never get here:
        delete m;
      } catch (const ValueErrorException &) {
        ok = true;
      }
      REQUIRE(ok);
    }
  }
}

TEST_CASE(
    "Testing Github issue 1605: Inappropriate bad valence exception during partial sanitization") {
  std::string smiles = "C1=CC=CC=C1N(=O)=O";
  {  // easy to test; we shouldn't throw an exception. :-)
    RWMol *m = SmilesToMol(smiles, 0, false);
    REQUIRE(m);
    unsigned int failed;
    MolOps::sanitizeMol(
        *m, failed,
        MolOps::SANITIZE_SETAROMATICITY | MolOps::SANITIZE_ADJUSTHS);
    REQUIRE(!failed);
    delete m;
  }
}

TEST_CASE("Testing Github issue 1622: add MDL aromaticity perception") {
  {
    // rings that should be aromatic
    string aromaticSmis[] = {"C1=CC=CC=C1",  // benzene, of course
                                             // heterocyclics
                             "N1=CC=CC=C1",  // pyridine
                             "N1=CC=CC=N1",  // pyridazine
                             "N1=CC=CN=C1",  // pyrimidine
                             "N1=CC=NC=C1",  // pyrazine
                             "N1=CN=CN=C1",  // 1,3,5-triazine
                             // polycyclic aromatics
                             "C1=CC2=CC=CC=CC2=C1",            // azulene
                             "C1=CC=CC2=CC=CC=C12",            // 6-6 fused
                             "C1=CC2=CC=CC=CC=C12",            // 4-8 fused
                             "C1=CC=C2C(=C1)N=CC=N2",          // 6-6 with Ns
                             "C1=CN=CC2C=CC=CC1=2",            // 6-6
                             "C1=CC=C2C(=C1)N=C3C=CC=CC3=N2",  // 6-6-6
                             "C1=CN=NC2C=CC=CC1=2",            // 6-6 with Ns
                             // macrocycle aromatics
                             "C1=CC=CC=CC=CC=C1",              // 10 atoms
                             "C1=CC=CC=CC=CC=CC=CC=CC=CC=C1",  // 18 atoms
                             "N1=CN=NC=CC=CC=CC=CC=CC=CC=CC=CC=CC=CC=CC=CC=C1",
                             "EOS"};
    unsigned int i = 0;
    while (aromaticSmis[i] != "EOS") {
      string smi = aromaticSmis[i];
      // std::cerr << smi << std::endl;
      int debugParse = 0;
      bool sanitize = false;
      RWMol *mol = SmilesToMol(smi, debugParse, sanitize);
      REQUIRE(mol);
      unsigned int whatFailed = 0;
      unsigned int sanitFlags =
          MolOps::SANITIZE_ALL ^ MolOps::SANITIZE_SETAROMATICITY;
      MolOps::sanitizeMol(*mol, whatFailed, sanitFlags);
      MolOps::setAromaticity(*mol, MolOps::AROMATICITY_MDL);
      REQUIRE(mol->getAtomWithIdx(0)->getIsAromatic());
      delete mol;
      ++i;
    }
  }
  {
    // rings that should not be aromatic
    string nonaromaticSmis[] = {
        "C1=C[N]C=C1",  // radicals are not two electron donors
        // exocyclic double bonds disqualify us
        "C1(=O)C=CNC=C1", "C1(=C)C=CC(=C)C=C1", "C1(=O)C=CC(=O)C=C1",

        "C1#CC=CC=C1",  // not benzyne
        // five-membered heterocycles
        "C1=COC=C1",      // furan
        "C1=CSC=C1",      // thiophene
        "C1=CNC=C1",      // pyrrole
        "C1=COC=N1",      // oxazole
        "C1=CSC=N1",      // thiazole
        "C1=CNC=N1",      // imidzole
        "C1=CNN=C1",      // pyrazole
        "C1=CON=C1",      // isoxazole
        "C1=CSN=C1",      // isothiazole
        "C1=CON=N1",      // 1,2,3-oxadiazole
        "C1=CNN=N1",      // 1,2,3-triazole
        "N1=CSC=N1",      // 1,3,4-thiadiazole
        "C1=CS(=O)C=C1",  // not sure how to classify this example from the
                          // OEChem docs
        //  outside the second rows
        "C1=CC=C[Si]=C1", "C1=CC=CC=P1",
        // 5-membered heterocycles outside the second row
        "C1=C[Se]C=C1", "C1=C[Te]C=C1",

        "EOS"};
    unsigned int i = 0;
    while (nonaromaticSmis[i] != "EOS") {
      string smi = nonaromaticSmis[i];
      // std::cerr << smi << std::endl;
      int debugParse = 0;
      bool sanitize = false;
      RWMol *mol = SmilesToMol(smi, debugParse, sanitize);
      REQUIRE(mol);
      unsigned int whatFailed = 0;
      unsigned int sanitFlags =
          MolOps::SANITIZE_ALL ^ MolOps::SANITIZE_SETAROMATICITY;
      MolOps::sanitizeMol(*mol, whatFailed, sanitFlags);
      MolOps::setAromaticity(*mol, MolOps::AROMATICITY_MDL);
      REQUIRE(!(mol->getAtomWithIdx(0)->getIsAromatic()));
      delete mol;
      ++i;
    }
  }

  {
    // ring systems where part is aromatic, part not
    string mixedaromaticSmis[] = {
        "O1C=CC2=CC=CC=C12",         "S1C=CC2=CC=CC=C12",
        "N1C2=CC=CC=C2C2=CC=CC=C12", "N1C=CC2=CC=CC=C12",
        "N1C=NC2=CC=CC=C12",         "N1C=NC2=CN=CN=C12",
        "C1CCCC2=CC3=CCCCC3=CC2=1",  "EOS"};
    unsigned int i = 0;
    while (mixedaromaticSmis[i] != "EOS") {
      string smi = mixedaromaticSmis[i];
      // std::cerr << smi << std::endl;
      int debugParse = 0;
      bool sanitize = false;
      RWMol *mol = SmilesToMol(smi, debugParse, sanitize);
      REQUIRE(mol);
      unsigned int whatFailed = 0;
      unsigned int sanitFlags =
          MolOps::SANITIZE_ALL ^ MolOps::SANITIZE_SETAROMATICITY;
      MolOps::sanitizeMol(*mol, whatFailed, sanitFlags);
      MolOps::setAromaticity(*mol, MolOps::AROMATICITY_MDL);
      REQUIRE(!(mol->getAtomWithIdx(0)->getIsAromatic()));
      REQUIRE((mol->getAtomWithIdx(mol->getNumAtoms() - 1)->getIsAromatic()));
      delete mol;
      ++i;
    }
  }
}

TEST_CASE(
    "Testing Github issue 1703: Dative bonds interfere with kekulization and the perception of aromaticity") {
  {  // start with zero-order bonds
    SmilesParserParams ps;
    ps.sanitize = false;
    std::unique_ptr<RWMol> mol(SmilesToMol("C1=CC=NC=N1.[Fe]", ps));
    REQUIRE(mol);
    mol->addBond(5, 6, Bond::ZERO);
    MolOps::sanitizeMol(*mol);
    REQUIRE(mol->getBondBetweenAtoms(0, 1)->getIsAromatic());
    REQUIRE(mol->getAtomWithIdx(5)->getIsAromatic());
    MolOps::Kekulize(*mol);
  }
  {  // and dative bonds:
    SmilesParserParams ps;
    ps.sanitize = false;
    std::unique_ptr<RWMol> mol(SmilesToMol("C1=CC=NC=N1->[Fe]", ps));
    REQUIRE(mol);
    MolOps::sanitizeMol(*mol);
    REQUIRE(mol->getBondBetweenAtoms(0, 1)->getIsAromatic());
    REQUIRE(mol->getAtomWithIdx(5)->getIsAromatic());
    MolOps::Kekulize(*mol);
  }
}

TEST_CASE(
    "Testing github issue 1614: AssignStereochemistry incorrectly removing CIS/TRANS bond stereo") {
  {
    RWMol m;
    m.addAtom(new Atom(9), true, true);
    m.addAtom(new Atom(6), true, true);
    m.addAtom(new Atom(6), true, true);
    m.addAtom(new Atom(17), true, true);
    m.addBond(0, 1, Bond::SINGLE);
    m.addBond(2, 3, Bond::SINGLE);
    m.addBond(1, 2, Bond::DOUBLE);
    m.getBondBetweenAtoms(1, 2)->setStereoAtoms(0, 3);
    m.getBondBetweenAtoms(1, 2)->setStereo(Bond::STEREOTRANS);
    m.updatePropertyCache();

    {
      RWMol nm(m);
      bool force = true, cleanIt = true;
      MolOps::setDoubleBondNeighborDirections(nm);

      MolOps::assignStereochemistry(nm, cleanIt, force);
      // nm.debugMol(std::cerr);
      REQUIRE(nm.getBondBetweenAtoms(1, 2)->getStereo() > Bond::STEREOANY);
      std::string smi = MolToSmiles(nm, true);
      std::cerr << smi << std::endl;
      REQUIRE(smi == "F/C=C/Cl");
    }
    {
      RWMol nm(m);
      MolOps::addHs(nm);
      bool force = true, cleanIt = true;
      MolOps::setDoubleBondNeighborDirections(nm);
      MolOps::assignStereochemistry(nm, cleanIt, force);
      // nm.debugMol(std::cerr);
      REQUIRE(nm.getBondBetweenAtoms(1, 2)->getStereo() > Bond::STEREOANY);
      std::string smi = MolToSmiles(nm, true);
      std::cerr << smi << std::endl;
      REQUIRE(smi == "[H]/C(F)=C(/[H])Cl");
    }
  }

  {
    RWMol m;
    m.addAtom(new Atom(9), true, true);
    m.addAtom(new Atom(6), true, true);
    m.addAtom(new Atom(6), true, true);
    m.addAtom(new Atom(17), true, true);
    m.addBond(0, 1, Bond::SINGLE);
    m.addBond(3, 2, Bond::SINGLE);
    m.addBond(1, 2, Bond::DOUBLE);
    m.getBondBetweenAtoms(1, 2)->setStereoAtoms(0, 3);
    m.getBondBetweenAtoms(1, 2)->setStereo(Bond::STEREOTRANS);
    m.updatePropertyCache();

    {
      RWMol nm(m);
      bool force = true, cleanIt = true;
      MolOps::setDoubleBondNeighborDirections(nm);

      MolOps::assignStereochemistry(nm, cleanIt, force);
      // nm.debugMol(std::cerr);
      REQUIRE(nm.getBondBetweenAtoms(1, 2)->getStereo() > Bond::STEREOANY);
      std::string smi = MolToSmiles(nm, true);
      std::cerr << smi << std::endl;
      REQUIRE(smi == "F/C=C/Cl");
    }
  }
  {
    RWMol m;
    m.addAtom(new Atom(9), true, true);
    m.addAtom(new Atom(6), true, true);
    m.addAtom(new Atom(6), true, true);
    m.addAtom(new Atom(17), true, true);
    m.addBond(1, 0, Bond::SINGLE);
    m.addBond(2, 3, Bond::SINGLE);
    m.addBond(1, 2, Bond::DOUBLE);
    m.getBondBetweenAtoms(1, 2)->setStereoAtoms(0, 3);
    m.getBondBetweenAtoms(1, 2)->setStereo(Bond::STEREOTRANS);
    m.updatePropertyCache();

    {
      RWMol nm(m);
      bool force = true, cleanIt = true;
      MolOps::setDoubleBondNeighborDirections(nm);

      MolOps::assignStereochemistry(nm, cleanIt, force);
      // nm.debugMol(std::cerr);
      REQUIRE(nm.getBondBetweenAtoms(1, 2)->getStereo() > Bond::STEREOANY);
      std::string smi = MolToSmiles(nm, true);
      std::cerr << smi << std::endl;
      REQUIRE(smi == "F/C=C/Cl");
    }
  }

  {
    RWMol m;
    m.addAtom(new Atom(9), true, true);
    m.addAtom(new Atom(6), true, true);
    m.addAtom(new Atom(6), true, true);
    m.addAtom(new Atom(17), true, true);
    m.addBond(0, 1, Bond::SINGLE);
    m.addBond(2, 3, Bond::SINGLE);
    m.addBond(1, 2, Bond::DOUBLE);
    m.addAtom(new Atom(6), true, true);
    m.addAtom(new Atom(6), true, true);
    m.addBond(1, 5, Bond::SINGLE);
    m.addBond(2, 4, Bond::SINGLE);

    m.getBondBetweenAtoms(1, 2)->setStereoAtoms(0, 3);
    m.getBondBetweenAtoms(1, 2)->setStereo(Bond::STEREOTRANS);
    m.updatePropertyCache();

    {
      RWMol nm(m);
      bool force = true, cleanIt = true;
      MolOps::setDoubleBondNeighborDirections(nm);

      MolOps::assignStereochemistry(nm, cleanIt, force);
      // nm.debugMol(std::cerr);
      REQUIRE(nm.getBondBetweenAtoms(1, 2)->getStereo() > Bond::STEREOANY);
      std::string smi = MolToSmiles(nm, true);
      std::cerr << smi << std::endl;
      REQUIRE(smi == "C/C(F)=C(/C)Cl");
    }
  }
  {
    RWMol *m = SmilesToMol("F/C=C(\\C/C=C/C)C/C=C\\F", false, false);
    REQUIRE(m);
    MolOps::sanitizeMol(*m);

    {
      RWMol nm(*m);
      MolOps::setDoubleBondNeighborDirections(nm);
      // nm.debugMol(std::cerr);
      bool force = true, cleanIt = true;
      MolOps::assignStereochemistry(nm, cleanIt, force);
      // nm.debugMol(std::cerr);
      REQUIRE(nm.getBondBetweenAtoms(4, 5)->getStereo() == Bond::STEREOE);
      REQUIRE(nm.getBondBetweenAtoms(8, 9)->getStereo() == Bond::STEREOZ);
      REQUIRE(nm.getBondBetweenAtoms(1, 2)->getStereo() == Bond::STEREOE);
    }
    delete m;
  }

  {
    RWMol *m = SmilesToMol("FC=C(C/C=C/C)C/C=C\\F", false, false);
    REQUIRE(m);
    MolOps::sanitizeMol(*m);

    REQUIRE(m->getBondBetweenAtoms(1, 2));
    m->getBondBetweenAtoms(1, 2)->setStereoAtoms(0, 3);
    m->getBondBetweenAtoms(1, 2)->setStereo(Bond::STEREOCIS);

    {
      RWMol nm(*m);
      MolOps::setDoubleBondNeighborDirections(nm);
      // nm.debugMol(std::cerr);
      bool force = true, cleanIt = true;
      MolOps::assignStereochemistry(nm, cleanIt, force);
      REQUIRE(nm.getBondBetweenAtoms(4, 5)->getStereo() == Bond::STEREOE);
      REQUIRE(nm.getBondBetweenAtoms(8, 9)->getStereo() == Bond::STEREOZ);
      REQUIRE(nm.getBondBetweenAtoms(1, 2)->getStereo() == Bond::STEREOE);
    }
    delete m;
  }

  {
    RWMol *m = SmilesToMol("F/C=C(\\C/C=C/C)C/C=C\\C", false, false);
    REQUIRE(m);
    MolOps::sanitizeMol(*m);

    {
      RWMol nm(*m);
      MolOps::setDoubleBondNeighborDirections(nm);
      bool force = true, cleanIt = true;
      MolOps::assignStereochemistry(nm, cleanIt, force);
      REQUIRE(nm.getBondBetweenAtoms(4, 5)->getStereo() == Bond::STEREOE);
      REQUIRE(nm.getBondBetweenAtoms(8, 9)->getStereo() == Bond::STEREOZ);
      REQUIRE(nm.getBondBetweenAtoms(1, 2)->getStereo() == Bond::STEREOE);
    }
    delete m;
  }
  {
    RWMol *m = SmilesToMol("FC=C(C/C=C/C)C/C=C\\C", false, false);
    REQUIRE(m);
    MolOps::sanitizeMol(*m);

    REQUIRE(m->getBondBetweenAtoms(1, 2));
    m->getBondBetweenAtoms(1, 2)->setStereoAtoms(0, 3);
    m->getBondBetweenAtoms(1, 2)->setStereo(Bond::STEREOCIS);

    {
      RWMol nm(*m);
      MolOps::setDoubleBondNeighborDirections(nm);
      bool force = true, cleanIt = true;
      MolOps::assignStereochemistry(nm, cleanIt, force);
      REQUIRE(nm.getBondBetweenAtoms(4, 5)->getStereo() == Bond::STEREOE);
      REQUIRE(nm.getBondBetweenAtoms(8, 9)->getStereo() == Bond::STEREOZ);
      REQUIRE(nm.getBondBetweenAtoms(1, 2)->getStereo() == Bond::STEREOE);
    }
    delete m;
  }
}

TEST_CASE(
    "Testing Github issue 1810: removeHs() should not remove H atoms that are contributing to the definition of a stereo bond") {
  {
    std::unique_ptr<RWMol> mol(SmilesToMol("F/C=C/[H]"));
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 4);
    REQUIRE(mol->getBondBetweenAtoms(1, 2)->getStereo() == Bond::STEREOE);
  }
  {
    std::unique_ptr<RWMol> mol(SmilesToMol("F/C=C(/F)[H]"));
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 4);
    REQUIRE(mol->getBondBetweenAtoms(1, 2)->getStereo() == Bond::STEREOE);
  }
  {
    std::unique_ptr<RWMol> mol(SmilesToMol("F/C=C(/[H])F"));
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 4);
    REQUIRE(mol->getBondBetweenAtoms(1, 2)->getStereo() == Bond::STEREOZ);
  }
  {
    std::unique_ptr<RWMol> mol(SmilesToMol("FC=C(F)[H]", false, false));
    REQUIRE(mol);
    MolOps::sanitizeMol(*mol);
    REQUIRE(mol->getNumAtoms() == 5);
    mol->getBondBetweenAtoms(1, 2)->setStereoAtoms(0, 4);
    mol->getBondBetweenAtoms(1, 2)->setStereo(Bond::STEREOTRANS);
    MolOps::removeHs(*mol);
    REQUIRE(mol->getNumAtoms() == 4);
    REQUIRE(mol->getBondBetweenAtoms(1, 2)->getStereoAtoms()[0] == 0);
    REQUIRE(mol->getBondBetweenAtoms(1, 2)->getStereoAtoms()[1] == 3);
  }
}

TEST_CASE(
    "Testing Github issue 1936: Bad aromaticity for rings with radical carbocations") {
  {
    std::unique_ptr<RWMol> mol(SmilesToMol("C1=CC=C[C+]C=C1"));
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 7);
    REQUIRE(mol->getAtomWithIdx(4)->getNumRadicalElectrons() == 1);
    REQUIRE(!mol->getAtomWithIdx(0)->getIsAromatic());
  }

  {  // the original report
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    std::unique_ptr<RWMol> mol(MolFileToMol(pathName + "github1936.mol"));
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 7);
    REQUIRE(mol->getAtomWithIdx(4)->getNumRadicalElectrons() == 1);
    REQUIRE(!mol->getAtomWithIdx(0)->getIsAromatic());
  }
}

TEST_CASE(
    "Testing Github issue 1928: incorrect aromatic SMILES generated for structure") {
  {
    std::unique_ptr<ROMol> mol(
        SmilesToMol("N1C2=CC3=CC=CC=C3OC1=CC1=C(O2)C=CC=C1"));
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 19);
    REQUIRE(mol->getBondBetweenAtoms(0, 1)->getBondType() == Bond::SINGLE);
    REQUIRE(!mol->getBondBetweenAtoms(0, 1)->getIsAromatic());
    REQUIRE(!mol->getAtomWithIdx(0)->getIsAromatic());
  }
  {  // the original report
    std::unique_ptr<ROMol> mol(SmilesToMol(
        "C12=C3C=CC=C1CCC(=O)C2=C4OC5=CC=CC6=C5C(=C(N4)O3)C(=O)CC6"));
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 27);
    REQUIRE(mol->getBondBetweenAtoms(20, 21)->getBondType() == Bond::SINGLE);
    REQUIRE(!mol->getBondBetweenAtoms(20, 21)->getIsAromatic());
    REQUIRE(!mol->getAtomWithIdx(21)->getIsAromatic());
  }
}

TEST_CASE("Testing Github issue 1990: removeHs screws up bond stereo") {
  {
    std::unique_ptr<RWMol> mol(SmilesToMol("F/C=C/F"));
    REQUIRE(mol);
    MolOps::addHs(*mol);
    MolOps::removeHs(*mol);
    REQUIRE(mol->getNumAtoms() == 4);
    REQUIRE(mol->getBondWithIdx(1)->getStereoAtoms().size() == 2);
  }
  {  // make sure that stereo is not removed when it comes from Hs:
    std::unique_ptr<RWMol> mol(SmilesToMol("F/C=C/F"));
    REQUIRE(mol);
    MolOps::addHs(*mol);
    REQUIRE(mol->getBondWithIdx(1)->getStereoAtoms().size() == 2);
    mol->getBondWithIdx(1)->getStereoAtoms()[0] = 4;
    MolOps::removeHs(*mol);
    REQUIRE(mol->getBondWithIdx(1)->getStereoAtoms().size() == 2);
    mol->getBondWithIdx(1)->getStereoAtoms()[0] = 0;
  }
  {  // make sure that stereo is not removed when it comes from Hs:
    std::unique_ptr<RWMol> mol(SmilesToMol("F/C=C/F"));
    REQUIRE(mol);
    MolOps::addHs(*mol);
    REQUIRE(mol->getBondWithIdx(1)->getStereoAtoms().size() == 2);
    mol->getBondWithIdx(1)->getStereoAtoms()[1] = 5;
    MolOps::removeHs(*mol);
    REQUIRE(mol->getBondWithIdx(1)->getStereoAtoms().size() == 2);
    mol->getBondWithIdx(1)->getStereoAtoms()[1] = 3;
  }
}

TEST_CASE("Testing removeAndTrackIsotopes parameter") {
  struct IsotopicHsCount {
    IsotopicHsCount(const ROMol &mol) {
      for (auto b : mol.bonds()) {
        const auto ba = b->getBeginAtom();
        const auto ea = b->getEndAtom();
        if (ba->getAtomicNum() == 1 && ba->getIsotope()) {
          ++d_map[ea->getIdx()];
        } else if (ea->getAtomicNum() == 1 && ea->getIsotope()) {
          ++d_map[ba->getIdx()];
        }
      }
    }
    // get the number of isotopic Hs attached to atom idx
    unsigned int at(unsigned int idx) {
      auto it = d_map.find(idx);
      return (it == d_map.end() ? 0 : it->second);
    }
    // count total number of isotopes in the molecule
    unsigned int total() {
      return std::accumulate(
          d_map.begin(), d_map.end(), 0U,
          [](unsigned int s, const std::pair<unsigned int, unsigned int> &p) {
            return s + p.second;
          });
    }
    static void countExplicitImplicitHs(const ROMol &m, unsigned int &expl,
                                        unsigned int &impl) {
      expl = 0;
      impl = 0;
      for (auto a : m.atoms()) {
        expl += a->getNumExplicitHs();
        impl += a->getNumImplicitHs();
      }
    }
    std::map<unsigned int, unsigned int> d_map;
  };
  // CHEMBL2024142
  auto m =
      "[2H]C1=C(C(=C2C(=C1[2H])C(=O)C(=C(C2=O)C([2H])([2H])[2H])C/C=C(\\C)/CC([2H])([2H])/C=C(/CC/C=C(\\C)/CCC=C(C)C)\\C([2H])([2H])[2H])[2H])[2H]"_smiles;
  REQUIRE(m.get());
  std::unique_ptr<IsotopicHsCount> m_isotopicHsPerHeavy(
      new IsotopicHsCount(*m));
  unsigned int m_numExplicitHs;
  unsigned int m_numImplicitHs;
  IsotopicHsCount::countExplicitImplicitHs(*m, m_numExplicitHs,
                                           m_numImplicitHs);
  REQUIRE(m_numExplicitHs == 0);
  REQUIRE(m_numImplicitHs == 28);
  REQUIRE(m_isotopicHsPerHeavy->total() == 12);
  MolOps::RemoveHsParameters ps;
  ps.removeAndTrackIsotopes = true;
  std::unique_ptr<ROMol> mNoH(removeHs(*static_cast<ROMol *>(m.get()), ps));
  REQUIRE(mNoH->getAtomWithIdx(0)->getAtomicNum() == 6);
  REQUIRE(mNoH->getAtomWithIdx(0)->hasProp(common_properties::_isotopicHs));
  std::vector<unsigned int> isoHs;
  REQUIRE(mNoH->getAtomWithIdx(0)->getPropIfPresent(
      common_properties::_isotopicHs, isoHs));
  REQUIRE(isoHs.size() == 1);
  REQUIRE(isoHs.front() == 2);
  REQUIRE(mNoH->getAtomWithIdx(30)->getAtomicNum() == 6);
  REQUIRE(!mNoH->getAtomWithIdx(30)->hasProp(common_properties::_isotopicHs));

  IsotopicHsCount mNoH_isotopicHsPerHeavy(*mNoH);
  unsigned int mNoH_numExplicitHs;
  unsigned int mNoH_numImplicitHs;
  IsotopicHsCount::countExplicitImplicitHs(*mNoH, mNoH_numExplicitHs,
                                           mNoH_numImplicitHs);
  REQUIRE(mNoH_numExplicitHs == 0);
  REQUIRE(mNoH_numImplicitHs == 40);
  REQUIRE(mNoH_isotopicHsPerHeavy.total() == 0);
  std::unique_ptr<ROMol> mH(MolOps::addHs(*mNoH));
  std::unique_ptr<IsotopicHsCount> mH_isotopicHsPerHeavy(
      new IsotopicHsCount(*mH));
  unsigned int mH_numExplicitHs;
  unsigned int mH_numImplicitHs;
  IsotopicHsCount::countExplicitImplicitHs(*mH, mH_numExplicitHs,
                                           mH_numImplicitHs);
  REQUIRE(mH_numExplicitHs == 0);
  REQUIRE(mH_numImplicitHs == 0);
  MatchVectType match;
  REQUIRE(SubstructMatch(*mH, *m, match));
  REQUIRE(match.size() == m->getNumAtoms());
  REQUIRE(mH_isotopicHsPerHeavy->total() == 12);
  std::unique_ptr<ROMol> mH2(MolOps::removeHs(*mH));
  REQUIRE(m->getNumAtoms() == mH2->getNumAtoms());
  std::unique_ptr<IsotopicHsCount> mH2_isotopicHsPerHeavy(
      new IsotopicHsCount(*mH2));
  unsigned int mH2_numExplicitHs;
  unsigned int mH2_numImplicitHs;
  IsotopicHsCount::countExplicitImplicitHs(*mH2, mH2_numExplicitHs,
                                           mH2_numImplicitHs);
  MatchVectType matchH2;
  REQUIRE(SubstructMatch(*m, *mH2, matchH2));
  REQUIRE(matchH2.size() == m->getNumAtoms());
  REQUIRE(mH2_isotopicHsPerHeavy->total() == 12);
  REQUIRE(mH2_numExplicitHs == 0);
  REQUIRE(mH2_numImplicitHs == 28);
  for (auto p : matchH2) {
    REQUIRE(mH2_isotopicHsPerHeavy->at(p.first) ==
            m_isotopicHsPerHeavy->at(p.second));
  }

  // shuffle atoms before adding Hs; result should not change
  std::vector<unsigned int> randomOrder(mNoH->getNumAtoms());
  std::iota(randomOrder.begin(), randomOrder.end(), 0U);
  std::shuffle(randomOrder.begin(), randomOrder.end(),
               std::default_random_engine());
  std::unique_ptr<ROMol> mNoHRen(MolOps::renumberAtoms(*mNoH, randomOrder));
  mH.reset(MolOps::addHs(*mNoHRen));
  mH_isotopicHsPerHeavy.reset(new IsotopicHsCount(*mH));
  IsotopicHsCount::countExplicitImplicitHs(*mH, mH_numExplicitHs,
                                           mH_numImplicitHs);
  REQUIRE(mH_numExplicitHs == 0);
  REQUIRE(mH_numImplicitHs == 0);
  MatchVectType matchRen;
  REQUIRE(SubstructMatch(*mH, *m, matchRen));
  REQUIRE(match != matchRen);
  REQUIRE(match.size() == matchRen.size());
  REQUIRE(mH_isotopicHsPerHeavy->total() == 12);
  mH2.reset(MolOps::removeHs(*mH));
  REQUIRE(m->getNumAtoms() == mH2->getNumAtoms());
  mH2_isotopicHsPerHeavy.reset(new IsotopicHsCount(*mH2));
  IsotopicHsCount::countExplicitImplicitHs(*mH2, mH2_numExplicitHs,
                                           mH2_numImplicitHs);
  MatchVectType matchH2Ren;
  REQUIRE(SubstructMatch(*m, *mH2, matchH2Ren));
  REQUIRE(matchH2 != matchH2Ren);
  REQUIRE(matchH2.size() == matchH2Ren.size());
  REQUIRE(mH2_isotopicHsPerHeavy->total() == 12);
  REQUIRE(mH2_numExplicitHs == 0);
  REQUIRE(mH2_numImplicitHs == 28);
  for (auto p : matchH2Ren) {
    REQUIRE(mH2_isotopicHsPerHeavy->at(p.first) ==
            m_isotopicHsPerHeavy->at(p.second));
  }

  // Add isotopes incrementally only on some atoms at a time
  // This should add 4 isotopes
  UINT_VECT onlyOnAtoms{0, 12};
  mH.reset(MolOps::addHs(*mNoH, false, false, &onlyOnAtoms));
  mH_isotopicHsPerHeavy.reset(new IsotopicHsCount(*mH));
  REQUIRE(mH_isotopicHsPerHeavy->total() == 4);
  REQUIRE(mH_isotopicHsPerHeavy->at(0) == 1);
  REQUIRE(mH_isotopicHsPerHeavy->at(12) == 3);
  // This should add 4 more isotopes
  onlyOnAtoms = UINT_VECT{1, 2, 18};
  mH.reset(MolOps::addHs(*mH, false, false, &onlyOnAtoms));
  mH_isotopicHsPerHeavy.reset(new IsotopicHsCount(*mH));
  REQUIRE(mH_isotopicHsPerHeavy->total() == 8);
  REQUIRE(mH_isotopicHsPerHeavy->at(0) == 1);
  REQUIRE(mH_isotopicHsPerHeavy->at(1) == 1);
  REQUIRE(mH_isotopicHsPerHeavy->at(2) == 1);
  REQUIRE(mH_isotopicHsPerHeavy->at(12) == 3);
  REQUIRE(mH_isotopicHsPerHeavy->at(18) == 2);
  // This should add the last 4 isotopes
  onlyOnAtoms = UINT_VECT{5, 32};
  mH.reset(MolOps::addHs(*mH, false, false, &onlyOnAtoms));
  mH_isotopicHsPerHeavy.reset(new IsotopicHsCount(*mH));
  REQUIRE(mH_isotopicHsPerHeavy->total() == 12);
  REQUIRE(mH_isotopicHsPerHeavy->at(0) == 1);
  REQUIRE(mH_isotopicHsPerHeavy->at(1) == 1);
  REQUIRE(mH_isotopicHsPerHeavy->at(2) == 1);
  REQUIRE(mH_isotopicHsPerHeavy->at(5) == 1);
  REQUIRE(mH_isotopicHsPerHeavy->at(12) == 3);
  REQUIRE(mH_isotopicHsPerHeavy->at(18) == 2);
  REQUIRE(mH_isotopicHsPerHeavy->at(32) == 3);
  match.clear();
  REQUIRE(SubstructMatch(*mH, *m, match));
  REQUIRE(match.size() == mH->getNumAtoms());
  for (auto p : match) {
    auto m_nIso = m_isotopicHsPerHeavy->at(p.first);
    if (!m_nIso) {
      continue;
    }
    auto mH_nIso = mH_isotopicHsPerHeavy->at(p.second);
    REQUIRE(m_nIso == mH_nIso);
  }

  // Check that chirality on centers which bear both non-isotopic
  // and isotopic Hs is preserved after...
  std::set<Atom::ChiralType> chiralTypeSet;
  std::set<Atom::ChiralType> chiralTypeSetAfterAddHs;
  std::set<Atom::ChiralType> chiralTypeSetAfterRemoveAllHsAddHs;
  std::set<Atom::ChiralType> chiralTypeSetAfterRemoveAllHsAddHsRemoveHs;
  for (unsigned int i : {0, 1}) {
    unsigned int hIdx = i + 24;
    std::string expectedCipCode(1, 'R' + i);
    std::unique_ptr<ROMol> mChiral(new ROMol(*m));
    mChiral->getAtomWithIdx(23)->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
    mChiral->getAtomWithIdx(hIdx)->setIsotope(0);
    MolOps::assignStereochemistry(*mChiral, true, true);
    REQUIRE(mChiral->getAtomWithIdx(23)->getProp<std::string>(
                common_properties::_CIPCode) == expectedCipCode);
    chiralTypeSet.insert(mChiral->getAtomWithIdx(23)->getChiralTag());
    // 1) ...Adding Hs
    mH.reset(MolOps::addHs(*mChiral));
    MolOps::assignStereochemistry(*mH, true, true);
    match.clear();
    REQUIRE(mH->getAtomWithIdx(23)->getProp<std::string>(
                common_properties::_CIPCode) == expectedCipCode);
    chiralTypeSetAfterAddHs.insert(mH->getAtomWithIdx(23)->getChiralTag());
    // 2) ...Removing all Hs including isotopes and then putting them back
    mNoH.reset(MolOps::removeHs(*static_cast<ROMol *>(mChiral.get()), ps));
    mH.reset(MolOps::addHs(*mNoH));
    MolOps::assignStereochemistry(*mH, true, true);
    match.clear();
    REQUIRE(SubstructMatch(*mH, *mChiral, match));
    REQUIRE(match.size() == mChiral->getNumAtoms());
    REQUIRE(mH->getAtomWithIdx(match[23].second)
                ->getProp<std::string>(common_properties::_CIPCode) ==
            expectedCipCode);
    chiralTypeSetAfterRemoveAllHsAddHs.insert(
        mH->getAtomWithIdx(match[23].second)->getChiralTag());
    // 3) ...Removing non-isotopic Hs
    mNoH.reset(MolOps::removeHs(*mH));
    MolOps::assignStereochemistry(*mNoH, true, true);
    REQUIRE(mNoH->getAtomWithIdx(match[23].second)
                ->getProp<std::string>(common_properties::_CIPCode) ==
            expectedCipCode);
    chiralTypeSetAfterRemoveAllHsAddHsRemoveHs.insert(
        mNoH->getAtomWithIdx(match[23].second)->getChiralTag());
  }
  // CIP chirality is preserved because when all Hs are removed
  // the parity is inverted on one of the enantiomers such that
  // chirality is preserved also when Hs are implicit.
  // So we must find a single parity before calling removeHs,
  // and two afterwards.
  REQUIRE((chiralTypeSet.size() == 1 &&
           *chiralTypeSet.begin() == Atom::CHI_TETRAHEDRAL_CW));
  REQUIRE((chiralTypeSetAfterAddHs.size() == 1 &&
           *chiralTypeSetAfterAddHs.begin() == Atom::CHI_TETRAHEDRAL_CW));
  REQUIRE(chiralTypeSetAfterRemoveAllHsAddHs.size() == 2);
  REQUIRE(chiralTypeSetAfterRemoveAllHsAddHsRemoveHs.size() == 2);

  // Check that chirality on centers which bear different
  // H isotopes is preserved after...
  chiralTypeSet.clear();
  chiralTypeSetAfterAddHs.clear();
  chiralTypeSetAfterRemoveAllHsAddHs.clear();
  chiralTypeSetAfterRemoveAllHsAddHsRemoveHs.clear();
  for (unsigned int i : {0, 1}) {
    unsigned int hIdx = i + 24;
    std::string expectedCipCode(1, 'R' + i);
    std::unique_ptr<ROMol> mChiral(new ROMol(*m));
    mChiral->getAtomWithIdx(23)->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
    mChiral->getAtomWithIdx(hIdx)->setIsotope(3);
    MolOps::assignStereochemistry(*mChiral, true, true);
    REQUIRE(mChiral->getAtomWithIdx(23)->getProp<std::string>(
                common_properties::_CIPCode) == expectedCipCode);
    chiralTypeSet.insert(mChiral->getAtomWithIdx(23)->getChiralTag());
    // 1) ...Adding Hs
    mH.reset(MolOps::addHs(*mChiral));
    MolOps::assignStereochemistry(*mH, true, true);
    match.clear();
    REQUIRE(mH->getAtomWithIdx(23)->getProp<std::string>(
                common_properties::_CIPCode) == expectedCipCode);
    chiralTypeSetAfterAddHs.insert(mH->getAtomWithIdx(23)->getChiralTag());
    // 2) ...Removing all Hs including isotopes and then putting them back
    mNoH.reset(MolOps::removeHs(*static_cast<ROMol *>(mChiral.get()), ps));
    mH.reset(MolOps::addHs(*mNoH));

    MolOps::assignStereochemistry(*mH, true, true);
    match.clear();
    REQUIRE(SubstructMatch(*mH, *mChiral, match));
    REQUIRE(match.size() == mChiral->getNumAtoms());
    REQUIRE(mH->getAtomWithIdx(match[23].second)
                ->getProp<std::string>(common_properties::_CIPCode) ==
            expectedCipCode);
    chiralTypeSetAfterRemoveAllHsAddHs.insert(
        mH->getAtomWithIdx(match[23].second)->getChiralTag());
    // 3) ...Removing non-isotopic Hs
    mNoH.reset(MolOps::removeHs(*mH));
    MolOps::assignStereochemistry(*mNoH, true, true);
    REQUIRE(mNoH->getAtomWithIdx(match[23].second)
                ->getProp<std::string>(common_properties::_CIPCode) ==
            expectedCipCode);
    chiralTypeSetAfterRemoveAllHsAddHsRemoveHs.insert(
        mNoH->getAtomWithIdx(match[23].second)->getChiralTag());
  }
  // In this case we should find a single parity throughout
  // as inverting the positions of 2H and 3H will trigger
  // an inversion in CIP chirality without need for parity change
  REQUIRE((chiralTypeSet.size() == 1 &&
           *chiralTypeSet.begin() == Atom::CHI_TETRAHEDRAL_CCW));
  REQUIRE((chiralTypeSetAfterAddHs.size() == 1 &&
           *chiralTypeSetAfterAddHs.begin() == Atom::CHI_TETRAHEDRAL_CCW));
  REQUIRE((chiralTypeSetAfterRemoveAllHsAddHs.size() == 1 &&
           *chiralTypeSetAfterRemoveAllHsAddHs.begin() ==
               Atom::CHI_TETRAHEDRAL_CCW));
  REQUIRE((chiralTypeSetAfterRemoveAllHsAddHsRemoveHs.size() == 1 &&
           *chiralTypeSetAfterRemoveAllHsAddHsRemoveHs.begin() ==
               Atom::CHI_TETRAHEDRAL_CCW));
}

TEST_CASE(
    "Testing github issue 3854: AddHs creates H atom with nan coordinates on edge case 2D structure") {
  std::string molb = R"CTAB(
     RDKit          2D

  7  8  0  0  1  0  0  0  0  0999 V2000
    5.0014    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1764    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3514    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9389    1.1271    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3514    1.8412    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1764    1.8412    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5889    1.1271    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  1
  7  1  1  1
  2  3  1  0
  2  7  1  0
  3  4  1  0
  4  5  1  0
  5  6  1  0
  7  6  1  0
M  END)CTAB";
  bool sanitize = true;
  bool removeHs = false;
  std::unique_ptr<ROMol> m(MolBlockToMol(molb, sanitize, removeHs));
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 7);

  bool explicitOnly = false;
  bool addCoords = true;
  std::vector<unsigned> onlyOnAtoms = {1, 6};
  std::unique_ptr<ROMol> m2(
      MolOps::addHs(*m, explicitOnly, addCoords, &onlyOnAtoms));
  REQUIRE(m2);
  REQUIRE(m2->getNumAtoms() == 9);

  auto conf = m2->getConformer();
  for (auto i = 7; i < 9; ++i) {
    auto atom_pos = conf.getAtomPos(i);
    REQUIRE((!isnan(atom_pos.x) && !isnan(atom_pos.y) && !isnan(atom_pos.z)));
  }

  // check that we bisect the correct angle and point outside the rings
  auto v71 = conf.getAtomPos(7) - conf.getAtomPos(1);
  auto v21 = conf.getAtomPos(2) - conf.getAtomPos(1);
  auto v01 = conf.getAtomPos(0) - conf.getAtomPos(1);
  auto v61 = conf.getAtomPos(6) - conf.getAtomPos(1);
  REQUIRE(fabs(fabs(v71.dotProduct(v01)) - fabs(v71.dotProduct(v21))) < 1e-3);
  REQUIRE(v71.dotProduct(v61) < -1e-4);

  auto v86 = conf.getAtomPos(8) - conf.getAtomPos(6);
  auto v06 = conf.getAtomPos(0) - conf.getAtomPos(6);
  auto v56 = conf.getAtomPos(5) - conf.getAtomPos(6);
  auto v16 = conf.getAtomPos(1) - conf.getAtomPos(6);
  REQUIRE(fabs(fabs(v86.dotProduct(v56)) - fabs(v86.dotProduct(v06))) < 1e-3);
  REQUIRE(v86.dotProduct(v16) < -1e-4);
}

#ifdef RDK_USE_URF
TEST_CASE("Testing ring family calculation") {
  {
    std::string smiles = "C(C1C2C3C41)(C2C35)C45";  // cubane
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 8);
    REQUIRE(!m->getRingInfo()->areRingFamiliesInitialized());
    MolOps::findRingFamilies(*m);
    REQUIRE(m->getRingInfo()->isInitialized());
    REQUIRE(m->getRingInfo()->areRingFamiliesInitialized());
    int numURF = RDL_getNofURF(m->getRingInfo()->dp_urfData.get());
    int numRC = RDL_getNofRC(m->getRingInfo()->dp_urfData.get());
    REQUIRE(numRC == 6);
    REQUIRE(numURF == 6);

    int numRings = m->getRingInfo()->numRingFamilies();
    REQUIRE(numRings == 6);
    numRings = m->getRingInfo()->numRings();
    REQUIRE(numRings == 6);

    delete m;
  }
  {
    std::string smiles = "C1CC2CCC1CC1CCC(CC1)CC1CCC(CC1)CC1CCC(CC1)C2";
    ROMol *m = SmilesToMol(smiles);
    REQUIRE(m);
    REQUIRE(m->getNumAtoms() == 28);
    REQUIRE(!m->getRingInfo()->areRingFamiliesInitialized());
    MolOps::findRingFamilies(*m);
    REQUIRE(m->getRingInfo()->isInitialized());
    REQUIRE(m->getRingInfo()->areRingFamiliesInitialized());
    int numURF = RDL_getNofURF(m->getRingInfo()->dp_urfData.get());
    int numRC = RDL_getNofRC(m->getRingInfo()->dp_urfData.get());
    // std::cerr << " URF, RC " << numURF << " " << numRC << std::endl;
    REQUIRE(numURF == 5);
    REQUIRE(numRC == 20);
    int numRings = m->getRingInfo()->numRings();
    // std::cerr << "num rings: " << numRings << std::endl;
    REQUIRE(numRings == 14);
    REQUIRE(m->getRingInfo()->numRingFamilies() == 5);
    delete m;
  }
}
#endif

TEST_CASE("Testing adding coordinates to a terminal atom") {
  auto mol = R"CTAB(
     RDKit          2D

  6  6  0  0  0  0  0  0  0  0999 V2000
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
M  END)CTAB"_ctab;
  auto atom = new Atom(0);
  auto idx = mol->addAtom(atom);
  delete atom;
  mol->addBond(idx, 0);
  MolOps::setTerminalAtomCoords(static_cast<ROMol &>(*mol), idx, 0);
  auto &coord = mol->getConformer().getAtomPos(idx);
  REQUIRE((coord.x > 2.499 && coord.x < 2.501));
  REQUIRE((coord.y > -0.001 && coord.y < 0.001));
  REQUIRE((coord.z > -0.001 && coord.z < 0.001));
}

TEST_CASE("testing get3DDistanceMat()") {
  auto mol = R"CTAB(bogus example
     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.1000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2000    0.0000    0.1000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5000    0.0000    0.1000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
M  END)CTAB"_ctab;
  REQUIRE(mol);

  double *dm = MolOps::get3DDistanceMat(*mol);
  REQUIRE(dm);
  REQUIRE(dm[0] == 0.0);
  REQUIRE(dm[1] == 1.2);
  REQUIRE(dm[2] == 2.5);
  REQUIRE(dm[3] == 1.2);
  REQUIRE(dm[4] == 0.0);
  REQUIRE(dm[5] == 1.3);
  REQUIRE(dm[6] == 2.5);
  REQUIRE(dm[7] == 1.3);
  REQUIRE(dm[8] == 0.0);
  // this will use a cached version:
  double *dm2 = MolOps::get3DDistanceMat(*mol);
  REQUIRE(dm == dm2);

  int confId = -1;
  bool useAtomWts = true;
  dm = MolOps::get3DDistanceMat(*mol, confId, useAtomWts);
  REQUIRE(dm);
  REQUIRE(dm[0] == 1.0);
  REQUIRE(dm[1] == 1.2);
  REQUIRE(dm[2] == 2.5);
  REQUIRE(dm[3] == 1.2);
  REQUIRE(dm[4] == 1.0);
  REQUIRE(dm[5] == 1.3);
  REQUIRE(dm[6] == 2.5);
  REQUIRE(dm[7] == 1.3);
  REQUIRE(dm[8] == 6.0 / 8.0);
}

TEST_CASE(
    "Testing github issue 5099: Removing H preserving only wedged ones strips all H") {
  std::string smi{"FC([H])(O)Cl"};
  std::unique_ptr<RWMol> m{SmilesToMol(smi, false, false)};
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 5);
  m->getBondBetweenAtoms(1, 2)->setBondDir(Bond::BondDir::BEGINWEDGE);

  auto ps = MolOps::RemoveHsParameters();
  ps.showWarnings = false;

  // Remove all but Wedged/dashed H
  ps.removeWithWedgedBond = false;
  ps.removeDefiningBondStereo = true;
  ps.removeDegreeZero = true;
  ps.removeDummyNeighbors = true;
  ps.removeHigherDegrees = true;
  ps.removeHydrides = true;
  ps.removeInSGroups = true;
  ps.removeIsotopes = true;
  ps.removeMapped = true;
  ps.removeNonimplicit = true;
  ps.removeOnlyHNeighbors = true;
  ps.removeWithQuery = true;

  removeHs(*m, ps);

  // H shouldn't be removed
  REQUIRE(m->getNumAtoms() == 5);
}

TEST_CASE("Testing hasQueryHs") {
  const auto has_no_query_hs = std::make_pair(false, false);
  const auto has_only_query_hs = std::make_pair(true, false);
  const auto has_unmergeable_hs = std::make_pair(true, true);

  auto m0 = "CCCC"_smarts;
  REQUIRE(RDKit::MolOps::hasQueryHs(*m0) == has_no_query_hs);

  auto m = "[#1]"_smarts;
  REQUIRE(RDKit::MolOps::hasQueryHs(*m) == has_only_query_hs);

  auto m2 = "[#1,N]"_smarts;
  REQUIRE(RDKit::MolOps::hasQueryHs(*m2) == has_unmergeable_hs);

  // remove the negation
  auto recursive = "[$(C-[H])]"_smarts;
  REQUIRE(RDKit::MolOps::hasQueryHs(*recursive) == has_only_query_hs);

  auto recursive_or = "[$([C,#1])]"_smarts;
  REQUIRE(RDKit::MolOps::hasQueryHs(*recursive_or) == has_unmergeable_hs);

  // from rd_filters for something bigger
  auto keto_def_heterocycle =
      "[$(c([C;!R;!$(C-[N,O,S]);!$(C-[H])](=O))1naaaa1),$(c([C;!R;!$(C-[N,O,S]);!$(C-[H])](=O))1naa[n,s,o]1)]"_smarts;
  REQUIRE(RDKit::MolOps::hasQueryHs(*keto_def_heterocycle) ==
          has_only_query_hs);

  auto github7687 = "[#1,#6,#7]"_smarts;
  REQUIRE(RDKit::MolOps::hasQueryHs(*github7687) == has_unmergeable_hs);
  auto github7687b = "[1;#7,#1,#6]"_smarts;
  REQUIRE(RDKit::MolOps::hasQueryHs(*github7687b) == has_unmergeable_hs);
  auto github7687c = "[1&#7,#1,#6]"_smarts;
  REQUIRE(RDKit::MolOps::hasQueryHs(*github7687c) == has_unmergeable_hs);
}

TEST_CASE("Testing isRingFused") {
  auto molOrig = "C1C(C2CC3CCCCC3C12)C1CCCCC1"_smiles;
  {
    RWMol mol(*molOrig);
    auto ri = mol.getRingInfo();
    REQUIRE(ri->numRings() == 4);
    boost::dynamic_bitset<> fusedRings(ri->numRings());
    for (size_t i = 0; i < ri->numRings(); ++i) {
      fusedRings.set(i, ri->isRingFused(i));
    }
    REQUIRE(fusedRings.count() == 3);
    REQUIRE(fusedRings.size() - fusedRings.count() == 1);
    auto query = "[$(C1CCC1)]-@[$(C1CCCCC1)]"_smarts;
    MatchVectType matchVect;
    SubstructMatch(mol, *query, matchVect);
    REQUIRE(matchVect.size() == 2);
    mol.removeBond(matchVect.at(0).second, matchVect.at(1).second);
    MolOps::sanitizeMol(mol);
    REQUIRE(MolToSmiles(mol) == "C1CCC(CC2CCC2C2CCCCC2)CC1");
    REQUIRE(ri->numRings() == 3);
    fusedRings.resize(ri->numRings());
    for (size_t i = 0; i < ri->numRings(); ++i) {
      fusedRings.set(i, ri->isRingFused(i));
    }
    REQUIRE(fusedRings.count() == 0);
    REQUIRE(fusedRings.size() - fusedRings.count() == 3);
  }
  {
    RWMol mol(*molOrig);
    auto ri = mol.getRingInfo();
    REQUIRE(ri->numRings() == 4);
    boost::dynamic_bitset<> fusedRings(ri->numRings());
    for (size_t i = 0; i < ri->numRings(); ++i) {
      fusedRings.set(i, ri->isRingFused(i));
    }
    REQUIRE(fusedRings.count() == 3);
    REQUIRE(fusedRings.size() - fusedRings.count() == 1);
    std::vector<unsigned int> fusedBonds(ri->numRings());
    for (size_t i = 0; i < ri->numRings(); ++i) {
      fusedBonds[i] = ri->numFusedBonds(i);
    }
    REQUIRE(std::count(fusedBonds.begin(), fusedBonds.end(), 0) == 1);
    REQUIRE(std::count(fusedBonds.begin(), fusedBonds.end(), 1) == 2);
    REQUIRE(std::count(fusedBonds.begin(), fusedBonds.end(), 2) == 1);
    auto query =
        "[$(C1CCCCC1-!@[CX4;R1;r4])].[$(C1C(-!@[CX4;R1;r6])CC1)]"_smarts;
    MatchVectType matchVect;
    SubstructMatch(mol, *query, matchVect);
    REQUIRE(matchVect.size() == 2);
    mol.addBond(matchVect.at(0).second, matchVect.at(1).second, Bond::SINGLE);
    MolOps::sanitizeMol(mol);
    REQUIRE(MolToSmiles(mol) == "C1CCC2C(C1)CC1C2C2C3CCCCC3C12");
    REQUIRE(ri->numRings() == 5);
    fusedRings.resize(ri->numRings());
    for (size_t i = 0; i < ri->numRings(); ++i) {
      fusedRings.set(i, ri->isRingFused(i));
    }
    REQUIRE(fusedRings.count() == 5);
    REQUIRE(fusedRings.size() - fusedRings.count() == 0);
    fusedBonds.resize(ri->numRings());
    for (size_t i = 0; i < ri->numRings(); ++i) {
      fusedBonds[i] = ri->numFusedBonds(i);
    }
    REQUIRE(std::count(fusedBonds.begin(), fusedBonds.end(), 0) == 0);
    REQUIRE(std::count(fusedBonds.begin(), fusedBonds.end(), 1) == 2);
    REQUIRE(std::count(fusedBonds.begin(), fusedBonds.end(), 2) == 3);
  }
}
