//
//  Copyright (C) 2017 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include <iostream>
#include <string>
#include <map>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolBundle.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Resonance.h>

using namespace RDKit;

void testBaseFunctionality() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "    testBaseFunctionality" << std::endl;
  ROMOL_SPTR mol(SmilesToMol("CC[C@H](C)F"));

  {
    MolBundle bundle;
    TEST_ASSERT(bundle.size() == 0);
    TEST_ASSERT(bundle.addMol(mol) == 1);
    TEST_ASSERT(bundle.size() == 1);
    TEST_ASSERT(bundle.addMol(ROMOL_SPTR(SmilesToMol("CC[C@@H](C)F"))) == 2);
    TEST_ASSERT(bundle.size() == 2);

    MolBundle bundle2(bundle);
    TEST_ASSERT(bundle2.size() == 2);

    TEST_ASSERT(bundle.getMol(0)->getNumAtoms() == 5);
    TEST_ASSERT(bundle[0]->getNumAtoms() == 5);
  }

  {
    FixedMolSizeMolBundle bundle;
    TEST_ASSERT(bundle.size() == 0);
    TEST_ASSERT(bundle.addMol(mol) == 1);
    TEST_ASSERT(bundle.size() == 1);
    TEST_ASSERT(bundle.addMol(ROMOL_SPTR(SmilesToMol("CC[C@@H](C)F"))) == 2);
    TEST_ASSERT(bundle.size() == 2);

    FixedMolSizeMolBundle bundle2(bundle);
    TEST_ASSERT(bundle2.size() == 2);
  }

  BOOST_LOG(rdInfoLog) << "  Done." << std::endl;
}

void testSubstructFunctionality() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "    testSubstructFunctionality" << std::endl;
  {
    MolBundle bundle;
    TEST_ASSERT(bundle.addMol(ROMOL_SPTR(SmilesToMol("CC[C@H](Cl)F"))) == 1);
    TEST_ASSERT(bundle.addMol(ROMOL_SPTR(SmilesToMol("CC[C@@H](Cl)F"))) == 2);

    ROMOL_SPTR qry(SmilesToMol("C[C@@H](Cl)F"));
    MatchVectType match;
    TEST_ASSERT(SubstructMatch(bundle, *qry, match));
    TEST_ASSERT(SubstructMatch(bundle, *qry, match, true, true));
    ROMOL_SPTR qry2(SmilesToMol("C[C@@H](Br)F"));
    TEST_ASSERT(!SubstructMatch(bundle, *qry2, match));
    TEST_ASSERT(!SubstructMatch(bundle, *qry2, match, true, true));

    std::vector<MatchVectType> matches;
    TEST_ASSERT(SubstructMatch(bundle, *qry, matches) == 1);
    TEST_ASSERT(SubstructMatch(bundle, *qry, matches, false, true, true) == 1);
    TEST_ASSERT(SubstructMatch(bundle, *qry2, matches) == 0);
    TEST_ASSERT(SubstructMatch(bundle, *qry2, matches, false, true, true) == 0);
  }
  {
    MolBundle bundle;
    TEST_ASSERT(bundle.addMol(ROMOL_SPTR(SmilesToMol("C[C@H](Cl)F"))) == 1);
    TEST_ASSERT(bundle.addMol(ROMOL_SPTR(SmilesToMol("C[C@@H](Cl)F"))) == 2);

    ROMOL_SPTR mol(SmilesToMol("CC[C@@H](Cl)F"));
    MatchVectType match;
    TEST_ASSERT(SubstructMatch(*mol, bundle, match));
    TEST_ASSERT(SubstructMatch(*mol, bundle, match, true, true));
    ROMOL_SPTR mol2(SmilesToMol("CC[C@@H](Br)F"));
    TEST_ASSERT(!SubstructMatch(*mol2, bundle, match));
    TEST_ASSERT(!SubstructMatch(*mol2, bundle, match, true, true));

    std::vector<MatchVectType> matches;
    TEST_ASSERT(SubstructMatch(*mol, bundle, matches) == 1);
    TEST_ASSERT(SubstructMatch(*mol, bundle, matches, false, true, true) == 1);
    TEST_ASSERT(SubstructMatch(*mol2, bundle, matches) == 0);
    TEST_ASSERT(SubstructMatch(*mol2, bundle, matches, false, true, true) == 0);
  }
  {
    MolBundle bundle;
    TEST_ASSERT(bundle.addMol(ROMOL_SPTR(SmilesToMol("CC[C@H](Cl)F"))) == 1);
    TEST_ASSERT(bundle.addMol(ROMOL_SPTR(SmilesToMol("CC[C@@H](Cl)F"))) == 2);

    MolBundle qbundle;
    TEST_ASSERT(qbundle.addMol(ROMOL_SPTR(SmilesToMol("C[13C@H](Cl)F"))) == 1);
    TEST_ASSERT(qbundle.addMol(ROMOL_SPTR(SmilesToMol("C[C@@H](Cl)F"))) == 2);

    MatchVectType match;
    TEST_ASSERT(SubstructMatch(bundle, qbundle, match));
    TEST_ASSERT(SubstructMatch(bundle, qbundle, match, true, true));

    std::vector<MatchVectType> matches;
    TEST_ASSERT(SubstructMatch(bundle, qbundle, matches) == 1);
    TEST_ASSERT(SubstructMatch(bundle, qbundle, matches, false, true, true) ==
                1);
  }
  BOOST_LOG(rdInfoLog) << "  Done." << std::endl;
}

void testExceptions() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "    testExceptions" << std::endl;
  ROMOL_SPTR mol(SmilesToMol("C1CCC1"));

  {
    // FixedMolSizeMolBundle requires all molecules to have the same size
    FixedMolSizeMolBundle bundle;
    TEST_ASSERT(bundle.size() == 0);
    TEST_ASSERT(bundle.addMol(mol) == 1);
    TEST_ASSERT(bundle.size() == 1);
    {
      bool ok = false;
      try {
        bundle.addMol(ROMOL_SPTR(SmilesToMol("C1CC1")));
      } catch (const ValueErrorException &) {
        ok = true;
      }
      TEST_ASSERT(ok);
    }
    {
      bool ok = false;
      try {
        bundle.addMol(ROMOL_SPTR(SmilesToMol("CCCC")));
      } catch (const ValueErrorException &) {
        ok = true;
      }
      TEST_ASSERT(ok);
    }
    {
      bool ok = false;
      try {
        bundle.getMol(1);
      } catch (const IndexErrorException &) {
        ok = true;
      }
      TEST_ASSERT(ok);
    }
    {
      bool ok = false;
      try {
        bundle[1];
      } catch (const IndexErrorException &) {
        ok = true;
      }
      TEST_ASSERT(ok);
    }
  }
  {
    // MolBundle requires supports different size molecules
    MolBundle bundle;
    TEST_ASSERT(bundle.size() == 0);
    TEST_ASSERT(bundle.addMol(mol) == 1);
    TEST_ASSERT(bundle.size() == 1);
    bundle.addMol(ROMOL_SPTR(SmilesToMol("C1CC1")));
    bundle.addMol(ROMOL_SPTR(SmilesToMol("CCCC")));
    TEST_ASSERT(bundle.size() == 3);
  }
  BOOST_LOG(rdInfoLog) << "  Done." << std::endl;
}

int main() {
  RDLog::InitLogs();
#if 1
  testBaseFunctionality();
  testSubstructFunctionality();
  testExceptions();
#endif
  return 0;
}
