//
//  Copyright (C) 2003-2015 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <iostream>
using namespace std;
using namespace RDKit;

void test1() {
  auto *a1 = new Atom(6);
  auto *a2 = new Atom(*a1);
  Atom *a3 = a1->copy();
  delete a1;
  delete a2;
  delete a3;
}

void test2() {
  auto *mol = new RWMol();

  mol->addAtom(new Atom(6), false, true);
  mol->addAtom(new Atom(6), false, true);
  mol->addAtom(new Atom(8), false, true);
  mol->addBond(0, 1, Bond::SINGLE);
  mol->addBond(1, 2, Bond::SINGLE);
  mol->setAtomBookmark(mol->getAtomWithIdx(1), 1);
  mol->setBondBookmark(mol->getBondWithIdx(0), 2);
  CHECK_INVARIANT(mol->hasAtomBookmark(1), "");
  CHECK_INVARIANT(mol->getAtomWithBookmark(1)->getIdx() == 1, "");
  CHECK_INVARIANT(mol->hasBondBookmark(2), "");
  CHECK_INVARIANT(mol->getBondWithBookmark(2)->getIdx() == 0, "");

  auto *mol2 = new RWMol(*mol);
  CHECK_INVARIANT(mol2->hasAtomBookmark(1), "");
  CHECK_INVARIANT(mol2->getAtomWithBookmark(1)->getIdx() == 1, "");
  CHECK_INVARIANT(mol2->hasBondBookmark(2), "");
  CHECK_INVARIANT(mol2->getBondWithBookmark(2)->getIdx() == 0, "");
  delete mol;
  delete mol2;
}

void testQueryCopying() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing behavior of copied queries"
      << std::endl;

  {
    std::string smi =
        "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$"
        "(N-*=!@[O,N,P,S])]),$([nH0,o,s;+0])]";
    auto *m = static_cast<ROMol *>(SmartsToMol(smi));
    TEST_ASSERT(m);
    std::cerr << "\n\n\nCOPY" << std::endl;
    auto *m2 = new ROMol(*m, true);

    delete m;
    delete m2;
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testConformerCopying() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing copying conformers"
                       << std::endl;

  {
    std::string smi = "CCC";
    auto *m = static_cast<ROMol *>(SmilesToMol(smi));
    TEST_ASSERT(m);

    auto *conf = new Conformer(m->getNumAtoms());
    conf->setId(1);
    m->addConformer(conf, false);
    conf = new Conformer(m->getNumAtoms());
    conf->setId(2);
    m->addConformer(conf, false);
    {
      auto *m2 = new ROMol(*m);
      TEST_ASSERT(m2->getNumConformers() == 2);
      delete m2;
    }
    {
      auto *m2 = new ROMol(*m, false, 1);
      TEST_ASSERT(m2->getNumConformers() == 1);
      TEST_ASSERT(m2->getConformer().getId() == 1);
      delete m2;
    }
    {
      auto *m2 = new ROMol(*m, false, 2);
      TEST_ASSERT(m2->getNumConformers() == 1);
      TEST_ASSERT(m2->getConformer().getId() == 2);
      delete m2;
    }
    {
      auto *m2 = new ROMol(*m, false, 3);
      TEST_ASSERT(m2->getNumConformers() == 0);
      delete m2;
    }

    delete m;
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testGithub2144() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Github #2144: no "
                          "RWMol(RWMol) copy constructor"
                       << std::endl;

  {
    auto mol1 = "C1CCCCC1"_smiles;
    TEST_ASSERT(mol1);
    RWMol mol2(*mol1);
    TEST_ASSERT(mol2.getNumBonds() == mol1->getNumBonds())
    TEST_ASSERT(mol2.getNumAtoms() == mol1->getNumAtoms())
    RWMol mol3(mol2);
    TEST_ASSERT(mol2.getNumBonds() == mol3.getNumBonds())
    TEST_ASSERT(mol2.getNumAtoms() == mol3.getNumAtoms())
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

// -------------------------------------------------------------------
int main() {
  RDLog::InitLogs();

#if 1
  test1();
  test2();
#endif
  testQueryCopying();
  testConformerCopying();
  testGithub2144();
  return 0;
}
