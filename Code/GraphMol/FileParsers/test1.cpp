//
//  Copyright (C) 2002-2021 Greg Landrum and other RDKit contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/test.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Canon.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/MolPickler.h>
#include "FileParsers.h"
#include "SequenceParsers.h"
#include "SequenceWriters.h"
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/FileParsers/ProximityBonds.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/LocaleSwitcher.h>
#include <clocale>
#include <cstdlib>

#include <regex>
#include <string>
#include <fstream>
#include <boost/lexical_cast.hpp>

using namespace RDKit;

void test1() {
  BOOST_LOG(rdInfoLog) << "testing atom query parsing" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase + "/Code/GraphMol/FileParsers/test_data/list-query.mol";
  RWMol *m = MolFileToMol(fName, false);
  // MolOps::sanitizeMol(*m);
  TEST_ASSERT(m);

  TEST_ASSERT(m->getNumAtoms() == 6);

  std::string sma = MolToSmarts(*m);
  TEST_ASSERT(sma == "[#6]1=[#6]-[#6]=[#6]-[#6]=[#6,#7,#15]-1");

  std::string smi = MolToSmiles(*m);
  TEST_ASSERT(smi == "*1=CC=CC=C1");

  m->updatePropertyCache();
  smi = MolToSmarts(*m);
  TEST_ASSERT(smi == "[#6]1=[#6]-[#6]=[#6]-[#6]=[#6,#7,#15]-1");

  smi = "C1=CC=CC=C1";
  RWMol *m2 = SmilesToMol(smi, false, false);
  MatchVectType mv;
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 6);
  // sanitize it, which will aromatize the bonds... we will not match:
  MolOps::sanitizeMol(*m2);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 0);

  delete m2;
  smi = "N1=CC=CC=C1";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 6);
  delete m2;
  smi = "S1=CC=CC=C1";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 0);
  delete m2;
  smi = "P1=CC=CC=C1";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 6);
  delete m2;

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/not-list-query.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);

  smi = "CC(=C)C";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 4);
  delete m2;

  smi = "CC(=O)C";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;

  smi = "CC(=N)C";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;

  smi = "CC(=O)C(=C)C";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 4);
  delete m2;

  smi = "C(=C)C";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;

  // make sure new-style atom lists override old-style atom lists:
  delete m;
  fName = rdbase +
          "/Code/GraphMol/FileParsers/test_data/conflicting-list-query.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);

  smi = "CC(=C)C";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 4);
  delete m2;

  smi = "CC(=O)C";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;

  // longer list queries, this was issue 2413431:
  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/list-query-long.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(14)->hasQuery());

  smi = "C1COC2=CC3=CC4=C(C=CC=C4)C=C3C=C2C1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  delete m2;
  smi = "C1C[Se]C2=CC3=CC4=C(C=CC=C4)C=C3C=C2C1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  delete m2;
  smi = "C1C[Te]C2=CC3=CC4=C(C=CC=C4)C=C3C=C2C1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  delete m2;
  smi = "C1C[As]C2=CC3=CC4=C(C=CC=C4)C=C3C=C2C1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;

  delete m;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void test2() {
  BOOST_LOG(rdInfoLog) << "testing bond query parsing" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase + "/Code/GraphMol/FileParsers/test_data/bond-query.mol";

  RWMol *m = MolFileToMol(fName);

  TEST_ASSERT(m->getNumAtoms() == 5);
  std::string smi = MolToSmiles(*m);
  TEST_ASSERT(smi == "C=CC~CC");
  smi = "C1=CC=CC=C1";
  RWMol *m2 = SmilesToMol(smi, false, false);
  MatchVectType mv;
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 5);
  // sanitize it (making bonds aromatic) ... we will not match:
  MolOps::sanitizeMol(*m2);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 0);

  delete m2;
  smi = "C=CC=CC=C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 5);

  delete m2;
  smi = "C=CCCC=C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 5);

  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/RingBondQuery.mol";
  delete m;
  m = MolFileToMol(fName);
  TEST_ASSERT(m->getNumAtoms() == 5);
  delete m2;
  smi = "C1CCC1C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;
  smi = "C1CC2C1C2";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 5);

  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ChainBondQuery.mol";
  delete m;
  m = MolFileToMol(fName);
  TEST_ASSERT(m->getNumAtoms() == 5);
  delete m2;
  smi = "C1CCC1C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 5);
  delete m2;
  smi = "C1CC2C1C2";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));

  // - - - - - - - - - - - - - - - - - - - - - - - -
  // this was github issue #269
  // - - - - - - - - - - - - - - - - - - - - - - - -
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/bond-query4.mol";
  delete m;
  m = MolFileToMol(fName);
  TEST_ASSERT(m->getNumAtoms() == 5);
  delete m2;
  smi = "C1CCC1C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;
  smi = "C1CCC1=C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 5);
  delete m2;
  smi = "C1C=CC1=C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 5);
  delete m2;
  smi = "C1C#CC1=C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;
  smi = "CCCC=C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;
  smi = "CC=CC=C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;
  smi = "CC#CC=C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));

  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/bond-query5.mol";
  delete m;
  m = MolFileToMol(fName);
  TEST_ASSERT(m->getNumAtoms() == 5);
  delete m2;
  smi = "C1CCC1C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;
  smi = "C1CCC1=C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;
  smi = "C1C=CC1=C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;
  smi = "C1C#CC1=C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;
  smi = "CCCC=C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 5);
  delete m2;
  smi = "CC=CC=C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 5);
  delete m2;
  smi = "CC#CC=C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));

  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/bond-query6.mol";
  delete m;
  m = MolFileToMol(fName);
  TEST_ASSERT(m->getNumAtoms() == 5);
  delete m2;
  smi = "C1CCC1C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;
  smi = "C1CCC1=C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;
  smi = "C1C=CC1=C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;
  smi = "C1C#CC1=C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;
  smi = "CCCC=C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 5);
  delete m2;
  smi = "CC=CC=C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 5);
  delete m2;
  smi = "CC#CC=C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));

  delete m2;
  delete m;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void test4() {
  // basic writing test
  BOOST_LOG(rdInfoLog) << " ----------> Test4 " << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/";
  std::string fName = rdbase + "test_data/mol1.mol";

  RWMol *m = MolFileToMol(fName);
  TEST_ASSERT(m);
  std::string smi = MolToSmiles(*m);
  CHECK_INVARIANT(smi == "c1cc[cH-]c1", smi);

  TEST_ASSERT(m->getConformer().is3D() == false);

  std::string molBlock = MolToMolBlock(*m);
  delete m;
  m = MolBlockToMol(molBlock);
  TEST_ASSERT(m);
  smi = MolToSmiles(*m);
  CHECK_INVARIANT(smi == "c1cc[cH-]c1", smi);
  TEST_ASSERT(m->getConformer().is3D() == false);
  delete m;

  BOOST_LOG(rdInfoLog) << " Finished <---------- " << std::endl;
}

void test5() {
  // formerly problematic molecules
  BOOST_LOG(rdInfoLog) << " ----------> Test5 " << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/";
  std::string fName = rdbase + "test_data/issue123.mol";

  RWMol *m = MolFileToMol(fName);
  CHECK_INVARIANT(m, "");
  TEST_ASSERT(m->getNumAtoms() == 23);
  TEST_ASSERT(m->getConformer().is3D() == true);
  std::string molBlock = MolToMolBlock(*m);
  delete m;
  m = MolBlockToMol(molBlock);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 23);
  TEST_ASSERT(m->getConformer().is3D() == true);

  delete m;
  // now try without removing the Hs:
  m = MolFileToMol(fName, true, false);
  CHECK_INVARIANT(m, "");
  TEST_ASSERT(m->getNumAtoms() == 39);
  delete m;
}

void test6() {
  BOOST_LOG(rdInfoLog) << "testing chirality parsing" << std::endl;

  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/";
  std::string fName = rdbase + "test_data/chiral1.mol";

  RWMol *m;
  std::string smi, cip;

  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 5);
  smi = MolToSmiles(*m, true);
  TEST_ASSERT(smi == "C[C@](F)(Cl)Br");
  delete m;

  fName = rdbase + "test_data/chiral1a.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 5);
  smi = MolToSmiles(*m, true);
  TEST_ASSERT(smi == "C[C@](F)(Cl)Br");
  delete m;

  fName = rdbase + "test_data/chiral2.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 5);
  smi = MolToSmiles(*m, true);
  TEST_ASSERT(smi == "C[C@@](F)(Cl)Br");
  delete m;

  fName = rdbase + "test_data/chiral2a.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 5);
  smi = MolToSmiles(*m, true);
  TEST_ASSERT(smi == "C[C@@](F)(Cl)Br");
  delete m;

  fName = rdbase + "test_data/chiral3.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
#if 1
  smi = MolToSmiles(*m, true);
  // BOOST_LOG(rdInfoLog) << " smi: " << smi << std::endl;
  TEST_ASSERT(smi == "C[C@H](F)Cl");
#endif
  delete m;

  fName = rdbase + "test_data/chiral3a.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
#if 1
  smi = MolToSmiles(*m, true);
  TEST_ASSERT(smi == "C[C@H](F)Cl");
#endif
  delete m;

  fName = rdbase + "test_data/chiral4.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");
#if 1
  smi = MolToSmiles(*m, true);
  TEST_ASSERT(smi == "C[C@@H](F)Cl");
#endif
  delete m;

  fName = rdbase + "test_data/chiral4a.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");
#if 1
  smi = MolToSmiles(*m, true);
  TEST_ASSERT(smi == "C[C@@H](F)Cl");
#endif
  delete m;

  fName = rdbase + "test_data/chiral5.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 5);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(!m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
#if 1
  smi = MolToSmiles(*m, true);
  TEST_ASSERT(smi == "CC(C)(Cl)Br");
#endif
  delete m;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void test7() {
  BOOST_LOG(rdInfoLog) << "testing roundtrip chirality parsing" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/";

  RWMol *m, *m2;
  std::string fName;
  std::string smi, molBlock, smi2, cip;

  fName = rdbase + "test_data/chiral1.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 5);
  smi = MolToSmiles(*m, true);
  TEST_ASSERT(smi == "C[C@](F)(Cl)Br");
  molBlock = MolToMolBlock(*m);
  m2 = MolBlockToMol(molBlock);
  TEST_ASSERT(m2)
  smi2 = MolToSmiles(*m2, true);
  TEST_ASSERT(smi == smi2);
  delete m;
  delete m2;

  fName = rdbase + "test_data/chiral2.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 5);
  smi = MolToSmiles(*m, true);
  TEST_ASSERT(smi == "C[C@@](F)(Cl)Br");
  molBlock = MolToMolBlock(*m);
  m2 = MolBlockToMol(molBlock);
  TEST_ASSERT(m2)
  smi2 = MolToSmiles(*m2, true);
  TEST_ASSERT(smi == smi2);
  delete m;
  delete m2;
  fName = rdbase + "test_data/chiral3.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
#if 1
  smi = MolToSmiles(*m, true);
  TEST_ASSERT(smi == "C[C@H](F)Cl");
  molBlock = MolToMolBlock(*m);
  m2 = MolBlockToMol(molBlock);
  TEST_ASSERT(m2)
  MolOps::assignStereochemistry(*m2);
  TEST_ASSERT(m2->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m2->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  smi2 = MolToSmiles(*m2, true);
  TEST_ASSERT(smi == smi2);
  delete m2;
#endif
  delete m;

  fName = rdbase + "test_data/chiral4.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");
#if 1
  smi = MolToSmiles(*m, true);
  TEST_ASSERT(smi == "C[C@@H](F)Cl");
  molBlock = MolToMolBlock(*m);
  // BOOST_LOG(rdInfoLog) << molBlock << std::endl;
  m2 = MolBlockToMol(molBlock);
  TEST_ASSERT(m2)
  MolOps::assignStereochemistry(*m2);
  TEST_ASSERT(m2->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m2->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");
  // smi2 = MolToSmiles(*m2,true);
  // TEST_ASSERT(smi==smi2);
  delete m2;
#endif
  delete m;

  fName = rdbase + "test_data/Issue142d.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
#if 1
  smi = MolToSmiles(*m, true);
  m2 = SmilesToMol(smi);
  smi2 = MolToSmiles(*m2, true);
  if (smi != smi2) {
    BOOST_LOG(rdInfoLog) << "\n " << smi << "\n !=\n " << smi2 << std::endl;
  }
  TEST_ASSERT(smi == smi2);
  delete m2;
  BOOST_LOG(rdInfoLog) << "SMI: " << smi << std::endl;
  std::cout << "***************************************" << std::endl;
  molBlock = MolToMolBlock(*m);
  std::cout << "***************************************" << std::endl;
  BOOST_LOG(rdInfoLog) << molBlock << std::endl;
  m2 = MolBlockToMol(molBlock);
  TEST_ASSERT(m2)
  smi2 = MolToSmiles(*m2, true);
  if (smi != smi2) {
    BOOST_LOG(rdInfoLog) << "\n " << smi << "\n !=\n " << smi2 << std::endl;
  }
  TEST_ASSERT(smi == smi2);
  delete m2;
#endif

  delete m;
  fName = rdbase + "test_data/Issue142b.mol";
  m = MolFileToMol(fName);
  // BOOST_LOG(rdInfoLog) << m->getNumAtoms() << "\n";
  // BOOST_LOG(rdInfoLog) << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" <<
  // std::endl;
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 9);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  TEST_ASSERT(m->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
#if 1
  smi = MolToSmiles(*m, true);
  m2 = SmilesToMol(smi);
  smi2 = MolToSmiles(*m2, true);
  if (smi != smi2) {
    BOOST_LOG(rdInfoLog) << "\n " << smi << "\n !=\n " << smi2 << std::endl;
  }
  TEST_ASSERT(smi == smi2);
  delete m2;
  // BOOST_LOG(rdInfoLog) << "SMI: "<< smi << std::endl;
  BOOST_LOG(rdInfoLog) << m->getNumAtoms() << " "
                       << m->getConformer().getNumAtoms() << "\n";
  molBlock = MolToMolBlock(*m);
  // BOOST_LOG(rdInfoLog) << molBlock << std::endl;
  m2 = MolBlockToMol(molBlock);
  TEST_ASSERT(m2)
  smi2 = MolToSmiles(*m2, true);
  if (smi != smi2) {
    BOOST_LOG(rdInfoLog) << "\n " << smi << "\n !=\n " << smi2 << std::endl;
  }
  TEST_ASSERT(smi == smi2);
  delete m2;
#endif
  delete m;

  fName = rdbase + "test_data/issue142a.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 28);
#if 1
  smi = MolToSmiles(*m, true);
  m2 = SmilesToMol(smi);
  smi2 = MolToSmiles(*m2, true);
  if (smi != smi2) {
    BOOST_LOG(rdInfoLog) << "\n " << smi << "\n !=\n " << smi2 << std::endl;
  }
  TEST_ASSERT(smi == smi2);
  delete m2;

  molBlock = MolToMolBlock(*m);
  // BOOST_LOG(rdInfoLog) << molBlock << std::endl;
  m2 = MolBlockToMol(molBlock);
  TEST_ASSERT(m2)
  smi2 = MolToSmiles(*m2, true);
  if (smi != smi2) {
    BOOST_LOG(rdInfoLog) << "\n " << smi << "\n !=\n " << smi2 << std::endl;
  }
  TEST_ASSERT(smi == smi2);
  delete m2;
#endif
  delete m;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void test8() {
  BOOST_LOG(rdInfoLog) << "testing reading without sanitization" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/";

  RWMol *m;
  std::string fName;
  std::string smi, molBlock, smi2;

  // in this case the test means to not remove Hs:
  fName = rdbase + "test_data/unsanitary.mol";
  m = MolFileToMol(fName, false);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 6);

  delete m;
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 3);

  delete m;
  fName = rdbase + "test_data/unsanitary2.mol";
  m = MolFileToMol(fName, false);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 9);
  delete m;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIssue145() {
  BOOST_LOG(rdInfoLog) << "testing Issue145:\n Mol parsing: molecule yields "
                          "non-canonical smiles from mol block"
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/";

  RWMol *m, *m2;
  std::string fName;
  std::string smi, molBlock, smi2;

  fName = rdbase + "test_data/issue145.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 19);
  smi = MolToSmiles(*m, true);
  m2 = SmilesToMol(smi);
  smi2 = MolToSmiles(*m2, true);
  if (smi != smi2) {
    BOOST_LOG(rdInfoLog) << "\n " << smi << "\n !=\n " << smi2 << std::endl;
  }
  TEST_ASSERT(smi == smi2);
  delete m2;

  molBlock = MolToMolBlock(*m);
  m2 = MolBlockToMol(molBlock);
  TEST_ASSERT(m2)
  smi2 = MolToSmiles(*m2, true);
  if (smi != smi2) {
    BOOST_LOG(rdInfoLog) << "\n " << smi << "\n !=\n " << smi2 << std::endl;
  }
  TEST_ASSERT(smi == smi2);
  delete m;
  delete m2;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIssue148() {
  BOOST_LOG(rdInfoLog) << "testing Issue148:\n Mol files containing mis-drawn "
                          "nitro groups not properly parsed"
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/";

  RWMol *m;
  std::string fName;
  std::string smi, molBlock, smi2;

  fName = rdbase + "test_data/issue148.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 9);
  TEST_ASSERT(m->getAtomWithIdx(0)->getFormalCharge() == 1);
  delete m;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIssue180() {
  BOOST_LOG(rdInfoLog) << "testing Issue180: bad Z/E assignments" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/";

  RWMol *m;
  std::string fName;
  std::string code;
  Bond *bond;

  fName = rdbase + "test_data/Issue180.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  bond = m->getBondWithIdx(2);
  TEST_ASSERT(bond->getBondType() == Bond::DOUBLE);
  TEST_ASSERT(bond->getStereo() == Bond::STEREOZ);

  bond = m->getBondWithIdx(5);
  TEST_ASSERT(bond->getBondType() == Bond::DOUBLE);
  TEST_ASSERT(bond->getStereo() == Bond::STEREOE);
  delete m;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIssue264() {
  BOOST_LOG(rdInfoLog) << "testing Issue264: bad stereochemistry from mol files"
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/";

  RWMol *m1, *m2;
  std::string smi1, smi2;
  std::string fName;

  fName = rdbase + "test_data/Issue264-1.mol";
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);

  fName = rdbase + "test_data/Issue264-2.mol";
  m2 = MolFileToMol(fName);
  TEST_ASSERT(m2);

  smi1 = MolToSmiles(*m1, false);
  smi2 = MolToSmiles(*m2, false);
  TEST_ASSERT(smi1 == smi2);

  smi1 = MolToSmiles(*m1, true);
  smi2 = MolToSmiles(*m2, true);
  BOOST_LOG(rdInfoLog) << smi1 << std::endl;
  BOOST_LOG(rdInfoLog) << smi2 << std::endl;
  TEST_ASSERT(smi1 != smi2);

  delete m1;
  delete m2;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIssue399() {
  BOOST_LOG(rdInfoLog) << "testing Issue399: bond wedging cleanup" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  RWMol *m1;
  std::string smi1, smi2;
  std::string fName;

  fName = rdbase + "Issue399a.mol";
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);
#if 1
  smi1 = MolToSmiles(*m1, true);
  TEST_ASSERT(smi1 == "C[C@H]1CO1");
#endif
  MolOps::assignStereochemistry(*m1);
  TEST_ASSERT(m1->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m1->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, smi2);
  TEST_ASSERT(smi2 == "S");
#if 1
  Chirality::wedgeMolBonds(*m1, &m1->getConformer());
  TEST_ASSERT(m1->getBondWithIdx(0)->getBondDir() == Bond::BEGINWEDGE);
  TEST_ASSERT(m1->getBondWithIdx(1)->getBondDir() == Bond::NONE);
  TEST_ASSERT(m1->getBondWithIdx(2)->getBondDir() == Bond::NONE);
  TEST_ASSERT(m1->getBondWithIdx(3)->getBondDir() == Bond::NONE);
#endif

  delete m1;

  // make sure we prefer wedging bonds to Hs:
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);
  MolOps::addHs(*m1, false, true);
  TEST_ASSERT(m1->getAtomWithIdx(7)->getAtomicNum() == 1);
  TEST_ASSERT(m1->getBondBetweenAtoms(1, 7));
  TEST_ASSERT(m1->getBondBetweenAtoms(1, 7)->getBondType() == Bond::SINGLE);
  TEST_ASSERT(m1->getBondBetweenAtoms(1, 7)->getBondDir() == Bond::NONE);

  Chirality::wedgeMolBonds(*m1, &m1->getConformer());
  TEST_ASSERT(m1->getBondBetweenAtoms(1, 7)->getBondDir() == Bond::BEGINDASH);

  delete m1;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMolFileChgLines() {
  BOOST_LOG(rdInfoLog) << "testing handling of charge lines" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  // SF.Net Issue1603923: problems with multiple chg lines
  {
    RWMol *m1;
    std::string fName;
    fName = rdbase + "MolFileChgBug.mol";
    m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getAtomWithIdx(24)->getFormalCharge() == -1);
    TEST_ASSERT(m1->getAtomWithIdx(25)->getFormalCharge() == -1);
    delete m1;
  }

  // many charges in one molecule:
  {
    RWMol *m1;
    std::string fName;
    fName = rdbase + "manycharges.mol";
    m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getFormalCharge() == -1);
    TEST_ASSERT(m1->getAtomWithIdx(13)->getFormalCharge() == -1);

    std::string molBlock = MolToMolBlock(*m1);
    // std::cerr<<molBlock<<std::endl;
    delete m1;
    m1 = MolBlockToMol(molBlock);
    // m1->debugMol(std::cerr);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getFormalCharge() == -1);
    TEST_ASSERT(m1->getAtomWithIdx(13)->getFormalCharge() == -1);
    delete m1;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testDblBondStereochem() {
  BOOST_LOG(rdInfoLog) << "testing basic double bond stereochemistry"
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  {
    RWMol *m1;
    std::string fName = rdbase + "simple_z.mol";
    m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getBondWithIdx(0)->getStereo() == Bond::STEREOZ);
    delete m1;
  }

  {
    RWMol *m1;
    std::string fName = rdbase + "simple_e.mol";
    m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getBondWithIdx(0)->getStereo() == Bond::STEREOE);
    delete m1;
  }

  {
    RWMol *m1;
    std::string fName = rdbase + "simple_either.mol";
    m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getBondWithIdx(0)->getStereo() == Bond::STEREOANY);
    delete m1;
  }

  // the next group for sf.net issue 3009836
  BOOST_LOG(rdInfoLog) << "  sub-test for issue 3099836" << std::endl;
  {
    RWMol *m1;
    std::string fName = rdbase + "Issue3009836.1.mol";
    m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getBondBetweenAtoms(3, 4)->getStereo() == Bond::STEREOZ);

    delete m1;
  }

  {
    RWMol *m1;
    std::string fName = rdbase + "Issue3009836.2.mol";
    m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getBondBetweenAtoms(3, 4)->getStereo() == Bond::STEREOZ);

    delete m1;
  }

  {
    RWMol *m1;
    std::string fName = rdbase + "Issue3009836.3.mol";
    m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getBondBetweenAtoms(6, 7)->getStereo() == Bond::STEREOE);
    TEST_ASSERT(m1->getBondBetweenAtoms(10, 11)->getStereo() == Bond::STEREOZ);

    delete m1;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testSymmetricDblBondStereochem() {
  // this was sf.net issue 1718794:
  // http://sourceforge.net/tracker/index.php?func=detail&aid=1718794&group_id=160139&atid=814650)
  BOOST_LOG(rdInfoLog) << "testing double bonds with symmetric substituents"
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  RWMol *m1;
  std::string fName, smi;

  fName = rdbase + "cistrans.1a.mol";
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);
  TEST_ASSERT(m1->getBondWithIdx(0)->getStereo() == Bond::STEREOE);
  smi = MolToSmiles(*m1, true);
  TEST_ASSERT(smi == "C/C=C/Cl");

  fName = rdbase + "cistrans.2a.mol";
  delete m1;
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);
  TEST_ASSERT(m1->getBondWithIdx(0)->getStereo() == Bond::STEREOZ);

  smi = MolToSmiles(*m1, true);
  TEST_ASSERT(smi == "C/C=C\\Cl");

  fName = rdbase + "cistrans.1.mol";
  delete m1;
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);
  TEST_ASSERT(m1->getBondWithIdx(0)->getStereo() == Bond::STEREOE);

  smi = MolToSmiles(*m1, true);
  TEST_ASSERT(smi == "C/C=C/C");

  fName = rdbase + "cistrans.2.mol";
  delete m1;
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);
  TEST_ASSERT(m1->getBondWithIdx(0)->getStereo() == Bond::STEREOZ);

  smi = MolToSmiles(*m1, true);
  TEST_ASSERT(smi == "C/C=C\\C");

  fName = rdbase + "cistrans.3.mol";
  delete m1;
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);
  TEST_ASSERT(m1->getBondWithIdx(0)->getStereo() == Bond::STEREOANY);

  smi = MolToSmiles(*m1, true);
  TEST_ASSERT(smi == "CC=CC");

  delete m1;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testRingDblBondStereochem() {
  // this was sf.net issue 1725068:
  // http://sourceforge.net/tracker/index.php?func=detail&aid=1725068&group_id=160139&atid=814650
  BOOST_LOG(rdInfoLog)
      << "testing double bonds in rings with stereochem specifications"
      << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  RWMol *m1;
  std::string fName, smi;

  fName = rdbase + "badringstereochem3.mol";
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);

  smi = MolToSmiles(*m1, true);
  TEST_ASSERT(smi.find("/", 0) == std::string::npos);
  TEST_ASSERT(smi.find("\\", 0) == std::string::npos);
  delete m1;

  fName = rdbase + "badringstereochem2.mol";
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);

  smi = MolToSmiles(*m1, true);
  TEST_ASSERT(smi.find("/", 0) == std::string::npos);
  TEST_ASSERT(smi.find("\\", 0) == std::string::npos);
  delete m1;

  fName = rdbase + "badringstereochem.mol";
  m1 = MolFileToMol(fName);
  TEST_ASSERT(m1);

  smi = MolToSmiles(*m1, true);
  TEST_ASSERT(smi.find("/", 0) == std::string::npos);
  TEST_ASSERT(smi.find("\\", 0) == std::string::npos);
  delete m1;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMolFileRGroups() {
  BOOST_LOG(rdInfoLog) << "testing mol file R-group parsing" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase + "/Code/GraphMol/FileParsers/test_data/rgroups1.mol";
  RWMol *m = MolFileToMol(fName);
  TEST_ASSERT(m);
  unsigned int idx;
  std::string label;

  TEST_ASSERT(m->getAtomWithIdx(3)->hasProp(common_properties::_MolFileRLabel));
  m->getAtomWithIdx(3)->getProp(common_properties::_MolFileRLabel, idx);
  TEST_ASSERT(idx == 2);
  TEST_ASSERT(m->getAtomWithIdx(3)->getAtomicNum() == 0);
  TEST_ASSERT(feq(m->getAtomWithIdx(3)->getIsotope(), 2));

  TEST_ASSERT(m->getAtomWithIdx(4)->hasProp(common_properties::_MolFileRLabel));
  m->getAtomWithIdx(4)->getProp(common_properties::_MolFileRLabel, idx);
  TEST_ASSERT(idx == 1);
  TEST_ASSERT(m->getAtomWithIdx(4)->getAtomicNum() == 0);
  TEST_ASSERT(feq(m->getAtomWithIdx(4)->getIsotope(), 1));

  //  test sf.net issue 3316600:
  TEST_ASSERT(m->getAtomWithIdx(3)->hasProp(common_properties::dummyLabel));
  m->getAtomWithIdx(3)->getProp(common_properties::dummyLabel, label);
  TEST_ASSERT(label == "R2");
  TEST_ASSERT(m->getAtomWithIdx(3)->getSymbol() == "R2");

  TEST_ASSERT(m->getAtomWithIdx(4)->hasProp(common_properties::dummyLabel));
  m->getAtomWithIdx(4)->getProp(common_properties::dummyLabel, label);
  TEST_ASSERT(label == "R1");
  TEST_ASSERT(m->getAtomWithIdx(4)->getSymbol() == "R1");

  RWMol *m2;
  MatchVectType mv;

  std::string smi;
  smi = "C1C(O)C1C";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 5);
  delete m2;

  smi = "C1CC(O)C1C";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;

  smi = "C1C(CO)C1CC";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 5);
  delete m2;

  smi = "CC(=O)CC";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/rgroups2.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(3)->hasProp(common_properties::_MolFileRLabel));
  m->getAtomWithIdx(3)->getProp(common_properties::_MolFileRLabel, idx);
  TEST_ASSERT(idx == 1);
  TEST_ASSERT(m->getAtomWithIdx(3)->getAtomicNum() == 0);
  TEST_ASSERT(feq(m->getAtomWithIdx(3)->getIsotope(), 1));

  TEST_ASSERT(m->getAtomWithIdx(4)->hasProp(common_properties::_MolFileRLabel));
  m->getAtomWithIdx(4)->getProp(common_properties::_MolFileRLabel, idx);
  TEST_ASSERT(idx == 1);
  TEST_ASSERT(m->getAtomWithIdx(4)->getAtomicNum() == 0);
  TEST_ASSERT(feq(m->getAtomWithIdx(4)->getIsotope(), 1));

  smi = "C1C(O)C1C";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 5);
  delete m2;

  smi = "C1CC(O)C1C";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;

  smi = "C1C(CO)C1CC";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 5);
  delete m2;

  smi = "CC(=O)CC";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/rgroups3.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(3)->hasProp(common_properties::_MolFileRLabel));
  m->getAtomWithIdx(3)->getProp(common_properties::_MolFileRLabel, idx);
  TEST_ASSERT(idx == 11);
  TEST_ASSERT(m->getAtomWithIdx(3)->getAtomicNum() == 0);
  TEST_ASSERT(feq(m->getAtomWithIdx(3)->getIsotope(), 11));
  TEST_ASSERT(m->getAtomWithIdx(4)->hasProp(common_properties::_MolFileRLabel));
  m->getAtomWithIdx(4)->getProp(common_properties::_MolFileRLabel, idx);
  TEST_ASSERT(idx == 503);
  TEST_ASSERT(m->getAtomWithIdx(4)->getAtomicNum() == 0);
  TEST_ASSERT(feq(m->getAtomWithIdx(4)->getIsotope(), 503));

  delete m;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMolFileDegreeQueries() {
  BOOST_LOG(rdInfoLog) << "testing mol file degree queries" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase + "/Code/GraphMol/FileParsers/test_data/subst1.mol";
  RWMol *m = MolFileToMol(fName);
  TEST_ASSERT(m);

  RWMol *m2;
  MatchVectType mv;
  std::string smi;

  smi = "CC(=O)O";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 4);
  delete m2;

  smi = "CC(=O)C";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 4);
  delete m2;

  smi = "CC(=O)";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/subst2.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);

  smi = "CC(=O)O";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 4);
  delete m2;

  smi = "CC(=O)C";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 4);
  delete m2;

  smi = "CC(O)C";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 4);
  delete m2;

  smi = "CC(O)";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;

  smi = "CC(O)(C)C";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/subst3.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);

  smi = "CC(=O)O";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 3);
  delete m2;

  smi = "CC(=O)C";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 3);
  delete m2;

  smi = "CC(O)C";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;

  smi = "CC(O)";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;

  smi = "CC(O)(C)C";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/subst4.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);

  smi = "CC(=O)O";
  m2 = SmilesToMol(smi, false, false);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  delete m2;

  {
    delete m;
    fName = rdbase + "/Code/GraphMol/FileParsers/test_data/combined.mol";
    m = MolFileToMol(fName);
    TEST_ASSERT(m);

    smi = "CC(=O)[CH-]C";
    m2 = SmilesToMol(smi, false, false);
    TEST_ASSERT(SubstructMatch(*m2, *m, mv));
    TEST_ASSERT(mv.size() == 4);
    delete m2;
    smi = "CC(=O)[C-](C)C";
    m2 = SmilesToMol(smi, false, false);
    TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
    delete m2;
    smi = "CC(=O)CC";
    m2 = SmilesToMol(smi, false, false);
    TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
    delete m2;
  }

  delete m;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMolFileRBCQueries() {
  BOOST_LOG(rdInfoLog) << "testing mol file ring-bond count queries"
                       << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName;
  RWMol *m;
  RWMol *m2;
  MatchVectType mv;
  std::string smi;

  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_2.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);

  smi = "CC";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 0);
  delete m2;
  smi = "C1CCC1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 2);
  delete m2;
  smi = "C12C3C4C1C5C2C3C45";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 0);
  delete m2;

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_3.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);

  smi = "CC";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 0);
  delete m2;
  smi = "C1CCC1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 0);
  delete m2;
  smi = "C12C3C4C1C5C2C3C45";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 2);
  delete m2;

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_0.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);

  smi = "CC";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 2);
  delete m2;
  smi = "C1CCC1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 0);
  delete m2;

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_4.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  smi = "CC";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 0);
  delete m2;
  smi = "C1CCC1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 0);
  delete m2;
  smi = "C12C3C4C1C5C2C3C45";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 0);
  delete m2;

  smi = "C1CS234C5CC2CC13CC4C5";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 2);
  delete m2;
  smi = "C1C2CC3CC4CC1S234";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 2);

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_star.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);

  delete m2;
  smi = "CC";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 2);
  delete m2;
  smi = "C1CCC1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 0);
  delete m2;

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_star2.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);

  smi = "CC";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 0);
  delete m2;
  smi = "C1CCC1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 4);
  delete m2;
  smi = "C1CC2C1CC2";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 4);
  delete m2;
  smi = "C12C3C4C1C5C2C3C45";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 0);
  delete m2;

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_star3.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);

  smi = "CC";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 0);
  delete m2;
  smi = "C1CCC1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 4);
  delete m2;
  smi = "C1CC2C1CC2";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 0);
  delete m2;
  smi = "C12C3C4C1C5C2C3C45";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 0);
  delete m2;

  delete m;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMolFileUnsaturationQueries() {
  BOOST_LOG(rdInfoLog) << "testing mol file unsaturation queries" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName;
  RWMol *m;
  RWMol *m2;
  MatchVectType mv;
  std::string smi;

  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/unsaturation.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);

  smi = "CO";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 0);
  delete m2;

  smi = "C=O";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 2);
  delete m2;

  smi = "CCO";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 0);
  delete m2;

  smi = "C=CO";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 2);
  delete m2;

  smi = "C#CO";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(SubstructMatch(*m2, *m, mv));
  TEST_ASSERT(mv.size() == 2);
  delete m2;

  delete m;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMolFileQueryToSmarts() {
  BOOST_LOG(rdInfoLog) << "testing mol file queries -> SMARTS " << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName;
  RWMol *m;
  std::string sma;

  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_2.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  sma = MolToSmarts(*m, true);
  TEST_ASSERT(sma == "[#6&x2]-[#6]")

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_3.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  sma = MolToSmarts(*m, true);
  TEST_ASSERT(sma == "[#6&x3]-[#6]")

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_0.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  sma = MolToSmarts(*m, true);
  TEST_ASSERT(sma == "[#6&x0]-[#6]")

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_4.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  sma = MolToSmarts(*m, true);
  TEST_ASSERT(sma == "[#16&x4]-[#6]")

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_star.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  sma = MolToSmarts(*m, true);
  TEST_ASSERT(sma == "[#6&x0]-[#6]")

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/ringcount_star2.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  sma = MolToSmarts(*m, true);
  TEST_ASSERT(sma.find("[#6&x2]") != std::string::npos);

  delete m;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/unsaturation.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  sma = MolToSmarts(*m, true);
  TEST_ASSERT(sma == "[#6&$(*=,:,#*)]~[#8]")

  delete m;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMissingFiles() {
  BOOST_LOG(rdInfoLog) << "testing handling of missing files" << std::endl;

  std::string fName;
  bool ok;
  RWMol *m;
  (void)m;
  fName = "bogus_file.mol";
  ok = false;
  try {
    m = MolFileToMol(fName);
  } catch (BadFileException &) {
    ok = true;
  }
  TEST_ASSERT(ok);

  ok = false;
  try {
    m = TPLFileToMol(fName);
  } catch (BadFileException &) {
    ok = true;
  }
  TEST_ASSERT(ok);

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIssue1965035() {
  BOOST_LOG(rdInfoLog)
      << "testing issue Issue1965035: problems with WedgeMolBonds "
      << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName;
  RWMol *m;
  std::string sma;

  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/Issue1965035.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);

  // the mol file parser removes bond wedging info:
  TEST_ASSERT(m->getBondWithIdx(4)->getBondDir() == Bond::NONE);
  // but a chiral tag is assigned:
  TEST_ASSERT(m->getAtomWithIdx(2)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);

  WedgeMolBonds(*m, &m->getConformer());
  TEST_ASSERT(m->getBondWithIdx(4)->getBondDir() == Bond::BEGINDASH);
  delete m;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testRadicals() {
  BOOST_LOG(rdInfoLog) << "testing handling of radicals " << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName;
  RWMol *m;
  std::string smiles;

  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/radical.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 0);
  TEST_ASSERT(m->getAtomWithIdx(1)->getNumRadicalElectrons() == 1);

  std::string molBlock = MolToMolBlock(*m);
  delete m;
  m = MolBlockToMol(molBlock);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 0);
  TEST_ASSERT(m->getAtomWithIdx(1)->getNumRadicalElectrons() == 1);
  delete m;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testBadBondOrders() {
  BOOST_LOG(rdInfoLog)
      << "testing handling of bogus bond orders (issue 2337369)" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName;
  RWMol *m;

  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/bondorder0.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getBondBetweenAtoms(0, 1)->getBondType() == Bond::UNSPECIFIED);
  TEST_ASSERT(!m->getBondBetweenAtoms(0, 1)->hasQuery());
  delete m;

  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/bondorder9.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(!m->getBondBetweenAtoms(0, 1)->hasQuery());
  TEST_ASSERT(m->getBondBetweenAtoms(0, 1)->getBondType() == Bond::DATIVE);

  delete m;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testAtomParity() {
  BOOST_LOG(rdInfoLog) << "testing handling of atom stereo parity flags"
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  {
    int parity;
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/parity.simple1.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(!m->getAtomWithIdx(0)->hasProp(common_properties::molParity));
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::molParity));
    m->getAtomWithIdx(1)->getProp(common_properties::molParity, parity);
    TEST_ASSERT(parity == 1);

    // look for the parity flags on output:
    auto molBlock = MolToMolBlock(*m);
    auto m2 = MolBlockToMol(molBlock);
    TEST_ASSERT(m2);
    TEST_ASSERT(!m2->getAtomWithIdx(0)->hasProp(common_properties::molParity));
    TEST_ASSERT(m2->getAtomWithIdx(1)->hasProp(common_properties::molParity));
    m2->getAtomWithIdx(1)->getProp(common_properties::molParity, parity);
    TEST_ASSERT(parity == 1);
    delete m2;

    // if we clear the stereo, no parity flags end up in the output:
    m->getAtomWithIdx(1)->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
    molBlock = MolToMolBlock(*m);
    m2 = MolBlockToMol(molBlock);
    TEST_ASSERT(m2);
    TEST_ASSERT(!m2->getAtomWithIdx(0)->hasProp(common_properties::molParity));
    TEST_ASSERT(!m2->getAtomWithIdx(1)->hasProp(common_properties::molParity));
    delete m2;

    delete m;
  }

  {
    int parity;
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/parity.simple2.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(!m->getAtomWithIdx(0)->hasProp(common_properties::molParity));
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::molParity));
    m->getAtomWithIdx(1)->getProp(common_properties::molParity, parity);
    TEST_ASSERT(parity == 2);

    MolOps::assignChiralTypesFrom3D(*m);
    std::string molBlock = MolToMolBlock(*m);
    RWMol *m2 = MolBlockToMol(molBlock);
    TEST_ASSERT(m2);
    TEST_ASSERT(!m2->getAtomWithIdx(0)->hasProp(common_properties::molParity));
    TEST_ASSERT(m2->getAtomWithIdx(1)->hasProp(common_properties::molParity));
    m2->getAtomWithIdx(1)->getProp(common_properties::molParity, parity);
    TEST_ASSERT(parity == 2);
    delete m2;

    delete m;
  }

  {
    // a case with an H on the chiral center:
    int parity;
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/parity.simpleH1.mol";
    RWMol *m = MolFileToMol(fName, true, false);
    TEST_ASSERT(m);
    TEST_ASSERT(!m->getAtomWithIdx(0)->hasProp(common_properties::molParity));
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::molParity));
    m->getAtomWithIdx(1)->getProp(common_properties::molParity, parity);
    TEST_ASSERT(parity == 1);

    MolOps::assignChiralTypesFrom3D(*m);
    std::string molBlock = MolToMolBlock(*m);
    RWMol *m2 = MolBlockToMol(molBlock, true, false);
    TEST_ASSERT(m2);
    TEST_ASSERT(!m2->getAtomWithIdx(0)->hasProp(common_properties::molParity));
    TEST_ASSERT(m2->getAtomWithIdx(1)->hasProp(common_properties::molParity));
    m2->getAtomWithIdx(1)->getProp(common_properties::molParity, parity);
    TEST_ASSERT(parity == 1);
    delete m2;

    // if we remove the H and write things out, we should
    // still get the right answer back:
    m2 = (RWMol *)MolOps::removeHs(*((ROMol *)m));
    molBlock = MolToMolBlock(*m2);
    delete m2;
    m2 = MolBlockToMol(molBlock);
    TEST_ASSERT(m2);
    TEST_ASSERT(!m2->getAtomWithIdx(0)->hasProp(common_properties::molParity));
    TEST_ASSERT(m2->getAtomWithIdx(1)->hasProp(common_properties::molParity));
    m2->getAtomWithIdx(1)->getProp(common_properties::molParity, parity);
    TEST_ASSERT(parity == 1);
    delete m2;
    delete m;
  }

  {
    // a case with an H on the chiral center:
    int parity;
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/parity.simpleH2.mol";
    RWMol *m = MolFileToMol(fName, true, false);
    TEST_ASSERT(m);
    TEST_ASSERT(!m->getAtomWithIdx(0)->hasProp(common_properties::molParity));
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::molParity));
    m->getAtomWithIdx(1)->getProp(common_properties::molParity, parity);
    TEST_ASSERT(parity == 2);

    MolOps::assignChiralTypesFrom3D(*m);
    std::string molBlock = MolToMolBlock(*m);
    RWMol *m2 = MolBlockToMol(molBlock, true, false);
    TEST_ASSERT(m2);
    TEST_ASSERT(!m2->getAtomWithIdx(0)->hasProp(common_properties::molParity));
    TEST_ASSERT(m2->getAtomWithIdx(1)->hasProp(common_properties::molParity));
    m2->getAtomWithIdx(1)->getProp(common_properties::molParity, parity);
    TEST_ASSERT(parity == 2);
    delete m2;

    m2 = (RWMol *)MolOps::removeHs(*((ROMol *)m));
    molBlock = MolToMolBlock(*m2);
    delete m2;
    m2 = MolBlockToMol(molBlock);
    TEST_ASSERT(m2);
    TEST_ASSERT(!m2->getAtomWithIdx(0)->hasProp(common_properties::molParity));
    TEST_ASSERT(m2->getAtomWithIdx(1)->hasProp(common_properties::molParity));
    m2->getAtomWithIdx(1)->getProp(common_properties::molParity, parity);
    TEST_ASSERT(parity == 2);
    delete m2;
    delete m;
  }

  {
    // a case with an N as the "chiral" center
    int parity;
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/parity.nitrogen.mol";
    RWMol *m = MolFileToMol(fName, true, false);
    TEST_ASSERT(m);
    TEST_ASSERT(!m->getAtomWithIdx(0)->hasProp(common_properties::molParity));
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::molParity));
    m->getAtomWithIdx(1)->getProp(common_properties::molParity, parity);
    TEST_ASSERT(parity == 1);

    MolOps::assignChiralTypesFrom3D(*m);
    std::string molBlock = MolToMolBlock(*m);
    RWMol *m2 = MolBlockToMol(molBlock, true, false);
    TEST_ASSERT(m2);
    TEST_ASSERT(!m2->getAtomWithIdx(0)->hasProp(common_properties::molParity));
    TEST_ASSERT(!m2->getAtomWithIdx(1)->hasProp(common_properties::molParity));
    delete m2;
    delete m;
  }

  {
    // a case with two Hs on the chiral center:
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/parity.twoHs.mol";
    RWMol *m = MolFileToMol(fName, true, false);
    TEST_ASSERT(m);
    TEST_ASSERT(!m->getAtomWithIdx(0)->hasProp(common_properties::molParity));
    TEST_ASSERT(!m->getAtomWithIdx(1)->hasProp(common_properties::molParity));

    // add a bogus chiral spec:
    m->getAtomWithIdx(0)->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
    std::string molBlock = MolToMolBlock(*m);
    auto *m2 = (RWMol *)MolOps::removeHs(*((ROMol *)m));
    molBlock = MolToMolBlock(*m2);
    delete m2;
    m2 = MolBlockToMol(molBlock, true, false);
    TEST_ASSERT(m2);
    TEST_ASSERT(!m2->getAtomWithIdx(0)->hasProp(common_properties::molParity));
    TEST_ASSERT(!m2->getAtomWithIdx(1)->hasProp(common_properties::molParity));
    delete m2;
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIssue2692246() {
  // basic writing test
  BOOST_LOG(rdInfoLog) << " Testing issue 2692246 " << std::endl;
  std::string smiles(120, 'C');
  smiles += "[CH2+]";
  RWMol *m = SmilesToMol(smiles);
  TEST_ASSERT(m);
  std::string molBlock = MolToMolBlock(*m);
  delete m;
  m = MolBlockToMol(molBlock);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 121);
  TEST_ASSERT(m->getAtomWithIdx(120)->getFormalCharge() == 1);
  delete m;

  BOOST_LOG(rdInfoLog) << " done" << std::endl;
}

void testKekulizationSkip() {
  // basic writing test
  BOOST_LOG(rdInfoLog) << " Testing mol blocks without kekulization "
                       << std::endl;
  std::string smiles("c1ccccc1");
  RWMol *m = SmilesToMol(smiles);
  TEST_ASSERT(m);
  std::string molBlock = MolToMolBlock(*m, true, -1, false);
  TEST_ASSERT(molBlock.find("1  2  4") != std::string::npos);
  TEST_ASSERT(molBlock.find("2  3  4") != std::string::npos);
  TEST_ASSERT(molBlock.find("3  4  4") != std::string::npos);

  molBlock = MolToMolBlock(*m);
  TEST_ASSERT(molBlock.find("1  2  4") == std::string::npos);
  TEST_ASSERT(molBlock.find("2  3  4") == std::string::npos);
  TEST_ASSERT(molBlock.find("3  4  4") == std::string::npos);

  delete m;
  BOOST_LOG(rdInfoLog) << " done" << std::endl;
}

void testMolFileAtomValues() {
  BOOST_LOG(rdInfoLog) << "testing atom values in mol files" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/";

  {
    RWMol *m;
    std::string fName, val;

    fName = rdbase + "test_data/AtomProps1.mol";
    m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(
        !m->getAtomWithIdx(0)->hasProp(common_properties::molFileValue));
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::molFileValue));
    m->getAtomWithIdx(1)->getProp(common_properties::molFileValue, val);
    TEST_ASSERT(val == "acidchloride");
    TEST_ASSERT(getAtomValue(m->getAtomWithIdx(1)) == "acidchloride")

    TEST_ASSERT(
        m->getAtomWithIdx(0)->hasProp(common_properties::molAtomMapNumber));
    TEST_ASSERT(
        m->getAtomWithIdx(1)->hasProp(common_properties::molAtomMapNumber));
    TEST_ASSERT(
        m->getAtomWithIdx(2)->hasProp(common_properties::molAtomMapNumber));
    TEST_ASSERT(
        !m->getAtomWithIdx(3)->hasProp(common_properties::molAtomMapNumber));

    TEST_ASSERT(m->getAtomWithIdx(0)->getAtomMapNum() == 1);
    TEST_ASSERT(m->getAtomWithIdx(1)->getAtomMapNum() == 2);
    TEST_ASSERT(m->getAtomWithIdx(2)->getAtomMapNum() == 3);
    TEST_ASSERT(m->getAtomWithIdx(3)->getAtomMapNum() == 0);

    // round trip
    m->getAtomWithIdx(0)->setAtomMapNum(4);
    setAtomRLabel(m->getAtomWithIdx(3), 1);
    setAtomAlias(m->getAtomWithIdx(0), "acidchloride");
    setAtomValue(m->getAtomWithIdx(0), "foobar");
    RWMol *m2 = MolBlockToMol(MolToMolBlock(*m));
    TEST_ASSERT(m2);
    TEST_ASSERT(m->getAtomWithIdx(0)->getAtomMapNum() == 4);
    TEST_ASSERT(m->getAtomWithIdx(1)->getAtomMapNum() == 2);
    TEST_ASSERT(m->getAtomWithIdx(2)->getAtomMapNum() == 3);
    TEST_ASSERT(m->getAtomWithIdx(3)->getAtomMapNum() == 0);
    TEST_ASSERT(getAtomRLabel(m->getAtomWithIdx(3)) == 1);
    TEST_ASSERT(getAtomAlias(m->getAtomWithIdx(0)) == "acidchloride");
    TEST_ASSERT(getAtomValue(m->getAtomWithIdx(0)) == "foobar");

    delete m;
    delete m2;
  }

  {
    RWMol *m;
    std::string fName, val;

    fName = rdbase + "test_data/AtomProps2.mol";
    m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(
        !m->getAtomWithIdx(0)->hasProp(common_properties::molFileValue));
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::molFileValue));
    m->getAtomWithIdx(1)->getProp(common_properties::molFileValue, val);
    TEST_ASSERT(val == "acidchloride");

    TEST_ASSERT(m->getAtomWithIdx(2)->hasProp(common_properties::molFileValue));
    m->getAtomWithIdx(2)->getProp(common_properties::molFileValue, val);
    TEST_ASSERT(val == "testing");

    TEST_ASSERT(m->getAtomWithIdx(3)->getFormalCharge() == -1);

    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMolFileAtomQueries() {
  BOOST_LOG(rdInfoLog) << "testing handling of A, Q, and * in mol files"
                       << std::endl;

  std::string rdbase = getenv("RDBASE");

  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/query_star.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);

    RWMol *m2;
    MatchVectType mv;
    std::string smi;

    smi = "[H]c1ccccc1";
    m2 = SmilesToMol(smi, false, false);
    MolOps::sanitizeMol(*m2);
    TEST_ASSERT(SubstructMatch(*m2, *m, mv));
    TEST_ASSERT(mv.size() == 7);
    delete m2;

    smi = "Cc1ccccc1";
    m2 = SmilesToMol(smi);
    TEST_ASSERT(SubstructMatch(*m2, *m, mv));
    TEST_ASSERT(mv.size() == 7);
    delete m2;

    smi = "Clc1ccccc1";
    m2 = SmilesToMol(smi);
    TEST_ASSERT(SubstructMatch(*m2, *m, mv));
    TEST_ASSERT(mv.size() == 7);
    delete m2;

    delete m;
  }
  {
    std::vector<std::string> fNames = {
        rdbase + "/Code/GraphMol/FileParsers/test_data/query_A.mol",
        rdbase + "/Code/GraphMol/FileParsers/test_data/query_A.v3k.mol",
    };
    for (const auto &fName : fNames) {
      std::unique_ptr<RWMol> m(MolFileToMol(fName));
      TEST_ASSERT(m);
      TEST_ASSERT(m->getAtomWithIdx(6)->hasQuery());
      TEST_ASSERT(m->getAtomWithIdx(6)->getQuery()->getTypeLabel() == "A");

      RWMol *m2;
      MatchVectType mv;
      std::string smi;

      smi = "[H]c1ccccc1";
      m2 = SmilesToMol(smi, false, false);
      MolOps::sanitizeMol(*m2);
      TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
      delete m2;

      smi = "Cc1ccccc1";
      m2 = SmilesToMol(smi);
      TEST_ASSERT(SubstructMatch(*m2, *m, mv));
      TEST_ASSERT(mv.size() == 7);
      delete m2;

      smi = "Clc1ccccc1";
      m2 = SmilesToMol(smi);
      TEST_ASSERT(SubstructMatch(*m2, *m, mv));
      TEST_ASSERT(mv.size() == 7);
      delete m2;
    }
  }
  {
    std::vector<std::string> fNames = {
        rdbase + "/Code/GraphMol/FileParsers/test_data/query_Q.mol",
        rdbase + "/Code/GraphMol/FileParsers/test_data/query_Q.v3k.mol",
    };
    for (const auto &fName : fNames) {
      std::unique_ptr<RWMol> m(MolFileToMol(fName));
      TEST_ASSERT(m);

      TEST_ASSERT(m->getAtomWithIdx(6)->hasQuery());
      TEST_ASSERT(m->getAtomWithIdx(6)->getQuery()->getTypeLabel() == "Q");
      RWMol *m2;
      MatchVectType mv;
      std::string smi;

      smi = "[H]c1ccccc1";
      m2 = SmilesToMol(smi, false, false);
      MolOps::sanitizeMol(*m2);
      TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
      delete m2;

      smi = "Cc1ccccc1";
      m2 = SmilesToMol(smi);
      TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
      delete m2;

      smi = "Clc1ccccc1";
      m2 = SmilesToMol(smi);
      TEST_ASSERT(SubstructMatch(*m2, *m, mv));
      TEST_ASSERT(mv.size() == 7);
      delete m2;
    }
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testListsAndValues() {
  BOOST_LOG(rdInfoLog)
      << "testing handling of mol files with atom lists and values"
      << std::endl;

  std::string rdbase = getenv("RDBASE");

  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/lists_plus_values.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);

    std::string value;

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::molFileValue));
    m->getAtomWithIdx(1)->getProp(common_properties::molFileValue, value);
    TEST_ASSERT(value == "halogen");

    TEST_ASSERT(m->getAtomWithIdx(1)->hasQuery());
    TEST_ASSERT(m->getAtomWithIdx(1)->getQuery()->getDescription() == "AtomOr");

    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void test1V3K() {
  BOOST_LOG(rdInfoLog) << "testing basic handling of v3000 mol files"
                       << std::endl;

  std::string rdbase = getenv("RDBASE");

  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/v3k.1.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 8);
    TEST_ASSERT(m->getNumBonds() == 8);

    delete m;
  }
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/v3k.3.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 9);
    TEST_ASSERT(m->getNumBonds() == 9);
    TEST_ASSERT(m->getAtomWithIdx(4)->getFormalCharge() == -1);
    TEST_ASSERT(m->getAtomWithIdx(4)->getIsotope() == 17);
    // m->debugMol(std::cerr);
    // TEST_ASSERT(m->getBondWithIdx(8)->getBondDir()==Bond::BEGINWEDGE);
    delete m;
  }
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/v3k.5a.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 5);
    TEST_ASSERT(m->getNumBonds() == 4);
    TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    delete m;
  }
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/v3k.5b.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 5);
    TEST_ASSERT(m->getNumBonds() == 4);
    TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    delete m;
  }
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/v3k.6a.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 4);
    TEST_ASSERT(m->getNumBonds() == 3);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo() == Bond::STEREOE);
    delete m;
  }
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/v3k.6b.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 4);
    TEST_ASSERT(m->getNumBonds() == 3);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo() == Bond::STEREOANY);
    delete m;
  }

  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/v3k.crash1.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 7);
    TEST_ASSERT(m->getNumBonds() == 7);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void test2V3K() {
  BOOST_LOG(rdInfoLog) << "testing more queries from v3000 mol files"
                       << std::endl;

  std::string rdbase = getenv("RDBASE");

  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/v3k.2.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 8);
    TEST_ASSERT(m->getNumBonds() == 8);

    std::string smiles = "O=C(O)C1OCCC1";
    RWMol *m2 = SmilesToMol(smiles);
    MatchVectType mv;
    TEST_ASSERT(SubstructMatch(*m2, *m, mv));

    delete m2;
    smiles = "O=C(O)C1OCCC1";
    m2 = SmilesToMol(smiles);
    TEST_ASSERT(SubstructMatch(*m2, *m, mv));

    delete m2;
    smiles = "O=C(O)C1SCCS1";
    m2 = SmilesToMol(smiles);
    TEST_ASSERT(SubstructMatch(*m2, *m, mv));

    delete m2;
    smiles = "O=C(O)C1OCCN1";
    m2 = SmilesToMol(smiles);
    TEST_ASSERT(SubstructMatch(*m2, *m, mv));

    delete m2;
    smiles = "O=C(O)C1OCCO1";
    m2 = SmilesToMol(smiles);
    TEST_ASSERT(!SubstructMatch(*m2, *m, mv));

    delete m2;
    delete m;
  }
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/v3k.4a.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 2);
    TEST_ASSERT(m->getNumBonds() == 1);

    std::string smiles = "OC1OCC1";
    RWMol *m2 = SmilesToMol(smiles);
    MatchVectType mv;
    TEST_ASSERT(SubstructMatch(*m2, *m, mv));

    delete m2;
    smiles = "C1OCC1";
    m2 = SmilesToMol(smiles);
    TEST_ASSERT(SubstructMatch(*m2, *m, mv));

    delete m2;
    smiles = "COCC";
    m2 = SmilesToMol(smiles);
    TEST_ASSERT(!SubstructMatch(*m2, *m, mv));

    delete m2;
    delete m;
  }
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/v3k.4b.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 2);
    TEST_ASSERT(m->getNumBonds() == 1);

    std::string smiles = "OC1OCC1";
    RWMol *m2 = SmilesToMol(smiles);
    MatchVectType mv;
    TEST_ASSERT(SubstructMatch(*m2, *m, mv));

    delete m2;
    smiles = "C1OCC1";
    m2 = SmilesToMol(smiles);
    TEST_ASSERT(!SubstructMatch(*m2, *m, mv));

    delete m2;
    smiles = "COCC";
    m2 = SmilesToMol(smiles);
    TEST_ASSERT(SubstructMatch(*m2, *m, mv));

    delete m2;
    delete m;
  }
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/v3k.rbc.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 3);
    TEST_ASSERT(m->getNumBonds() == 2);

    std::string smiles = "C1CC1";
    RWMol *m2 = SmilesToMol(smiles);
    TEST_ASSERT(m2);
    MatchVectType mv;
    TEST_ASSERT(SubstructMatch(*m2, *m, mv));

    delete m2;
    smiles = "CCC";
    m2 = SmilesToMol(smiles);
    TEST_ASSERT(m2);
    TEST_ASSERT(!SubstructMatch(*m2, *m, mv));

    delete m2;
    smiles = "N1NC2NNC12";
    m2 = SmilesToMol(smiles);
    TEST_ASSERT(m2);
    TEST_ASSERT(!SubstructMatch(*m2, *m, mv));

    delete m2;
    delete m;
  }

  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/chebi_15469.v3k.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 53);
    TEST_ASSERT(m->getNumBonds() == 55);
    TEST_ASSERT(m->getAtomWithIdx(52)->getAtomicNum() == 0);
    TEST_ASSERT(!m->getAtomWithIdx(52)->hasQuery());
    delete m;
  }

  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/chebi_57262.v3k.2.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 22);
    TEST_ASSERT(m->getNumBonds() == 21);
    TEST_ASSERT(m->getAtomWithIdx(18)->getAtomicNum() == 0);
    TEST_ASSERT(!m->getAtomWithIdx(18)->hasQuery());
    TEST_ASSERT(m->getAtomWithIdx(18)->getIsotope() == 1);
    TEST_ASSERT(m->getAtomWithIdx(21)->getAtomicNum() == 0);
    TEST_ASSERT(!m->getAtomWithIdx(21)->hasQuery());
    TEST_ASSERT(m->getAtomWithIdx(21)->getIsotope() == 2);
    delete m;
  }
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/chebi_57262.v3k.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 22);
    TEST_ASSERT(m->getNumBonds() == 21);
    TEST_ASSERT(m->getAtomWithIdx(18)->getAtomicNum() == 0);
    TEST_ASSERT(!m->getAtomWithIdx(18)->hasQuery());
    TEST_ASSERT(m->getAtomWithIdx(18)->getIsotope() == 1);
    TEST_ASSERT(m->getAtomWithIdx(21)->getAtomicNum() == 0);
    TEST_ASSERT(!m->getAtomWithIdx(21)->hasQuery());
    TEST_ASSERT(m->getAtomWithIdx(21)->getIsotope() == 2);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}
void test3V3K() {
  BOOST_LOG(rdInfoLog) << "testing basic writing of v3000 mol files"
                       << std::endl;

  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  std::string fName;

  {
    // charges
    fName = rdbase + "issue148.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 9);
    TEST_ASSERT(m->getAtomWithIdx(0)->getFormalCharge() == 1);

    std::string mb = MolToMolBlock(*m, true, -1, true, true);
    delete m;

    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 9);
    TEST_ASSERT(m->getAtomWithIdx(0)->getFormalCharge() == 1);

    delete m;
  }

  {
    // multiple charge lines
    fName = rdbase + "MolFileChgBug.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(24)->getFormalCharge() == -1);
    TEST_ASSERT(m->getAtomWithIdx(25)->getFormalCharge() == -1);

    std::string mb = MolToMolBlock(*m, true, -1, true, true);
    delete m;

    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(24)->getFormalCharge() == -1);
    TEST_ASSERT(m->getAtomWithIdx(25)->getFormalCharge() == -1);

    delete m;
  }
  {
    // radicals
    fName = rdbase + "radical.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 0);
    TEST_ASSERT(m->getAtomWithIdx(1)->getNumRadicalElectrons() == 1);

    std::string mb = MolToMolBlock(*m, true, -1, true, true);
    delete m;

    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 0);
    TEST_ASSERT(m->getAtomWithIdx(1)->getNumRadicalElectrons() == 1);

    delete m;
  }

  {
    // radical and valence
    fName = rdbase + "CH.v3k.mol";
    RWMol *m = MolFileToMol(fName);  // cxsmiles = "[CH] |(0,0,),^5:0|"
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 1);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNoImplicit());
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumExplicitHs() == 1);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 3);

    std::string mb = MolToMolBlock(*m, true, -1, true, true);
    delete m;

    // no bonds in this one, make sure there's no bond block:
    TEST_ASSERT(mb.find("BEGIN ATOM") != std::string::npos);
    TEST_ASSERT(mb.find("BEGIN BOND") == std::string::npos);

    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 1);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNoImplicit());
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumExplicitHs() == 1);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 3);

    delete m;
  }

  {
    // R Groups
    fName = rdbase + "rgroups1.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    unsigned int idx;
    std::string label;

    TEST_ASSERT(
        m->getAtomWithIdx(3)->hasProp(common_properties::_MolFileRLabel));
    m->getAtomWithIdx(3)->getProp(common_properties::_MolFileRLabel, idx);
    TEST_ASSERT(idx == 2);
    TEST_ASSERT(m->getAtomWithIdx(3)->getAtomicNum() == 0);
    TEST_ASSERT(feq(m->getAtomWithIdx(3)->getIsotope(), 2));
    TEST_ASSERT(
        m->getAtomWithIdx(4)->hasProp(common_properties::_MolFileRLabel));
    m->getAtomWithIdx(4)->getProp(common_properties::_MolFileRLabel, idx);
    TEST_ASSERT(idx == 1);
    TEST_ASSERT(m->getAtomWithIdx(4)->getAtomicNum() == 0);
    TEST_ASSERT(feq(m->getAtomWithIdx(4)->getIsotope(), 1));
    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp(common_properties::dummyLabel));
    m->getAtomWithIdx(3)->getProp(common_properties::dummyLabel, label);
    TEST_ASSERT(label == "R2");
    TEST_ASSERT(m->getAtomWithIdx(3)->getSymbol() == "R2");
    TEST_ASSERT(m->getAtomWithIdx(4)->hasProp(common_properties::dummyLabel));
    m->getAtomWithIdx(4)->getProp(common_properties::dummyLabel, label);
    TEST_ASSERT(label == "R1");
    TEST_ASSERT(m->getAtomWithIdx(4)->getSymbol() == "R1");

    std::string mb = MolToMolBlock(*m, true, -1, true, true);
    delete m;
    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(
        m->getAtomWithIdx(3)->hasProp(common_properties::_MolFileRLabel));
    m->getAtomWithIdx(3)->getProp(common_properties::_MolFileRLabel, idx);
    TEST_ASSERT(idx == 2);
    TEST_ASSERT(m->getAtomWithIdx(3)->getAtomicNum() == 0);
    TEST_ASSERT(feq(m->getAtomWithIdx(3)->getIsotope(), 2));
    TEST_ASSERT(
        m->getAtomWithIdx(4)->hasProp(common_properties::_MolFileRLabel));
    m->getAtomWithIdx(4)->getProp(common_properties::_MolFileRLabel, idx);
    TEST_ASSERT(idx == 1);
    TEST_ASSERT(m->getAtomWithIdx(4)->getAtomicNum() == 0);
    TEST_ASSERT(feq(m->getAtomWithIdx(4)->getIsotope(), 1));
    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp(common_properties::dummyLabel));
    m->getAtomWithIdx(3)->getProp(common_properties::dummyLabel, label);
    TEST_ASSERT(label == "R2");
    TEST_ASSERT(m->getAtomWithIdx(3)->getSymbol() == "R2");
    TEST_ASSERT(m->getAtomWithIdx(4)->hasProp(common_properties::dummyLabel));
    m->getAtomWithIdx(4)->getProp(common_properties::dummyLabel, label);
    TEST_ASSERT(label == "R1");
    TEST_ASSERT(m->getAtomWithIdx(4)->getSymbol() == "R1");

    delete m;
  }

  {
    // automatic cut over to v3k
    std::string smiles(1024, 'C');
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 1024);
    auto *conf = new Conformer(m->getNumAtoms());
    m->addConformer(conf, true);
    std::string mb = MolToMolBlock(*m);
    TEST_ASSERT(mb.find("V2000") == std::string::npos);
    TEST_ASSERT(mb.find("V3000") != std::string::npos);
    delete m;
  }

  {
    // D in CTAB
    fName = rdbase + "D_in_CTAB.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 3);
    TEST_ASSERT(m->getAtomWithIdx(2)->getAtomicNum() == 1);
    TEST_ASSERT(m->getAtomWithIdx(2)->getIsotope() == 2);

    std::string mb = MolToMolBlock(*m, true, -1, true, true);
    delete m;

    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 3);
    TEST_ASSERT(m->getAtomWithIdx(2)->getAtomicNum() == 1);
    TEST_ASSERT(m->getAtomWithIdx(2)->getIsotope() == 2);

    delete m;
  }
  {
    // T in CTAB
    fName = rdbase + "T_in_CTAB.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 3);
    TEST_ASSERT(m->getAtomWithIdx(2)->getAtomicNum() == 1);
    TEST_ASSERT(m->getAtomWithIdx(2)->getIsotope() == 3);

    std::string mb = MolToMolBlock(*m, true, -1, true, true);
    delete m;

    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 3);
    TEST_ASSERT(m->getAtomWithIdx(2)->getAtomicNum() == 1);
    TEST_ASSERT(m->getAtomWithIdx(2)->getIsotope() == 3);

    delete m;
  }

  {
    // atom list
    fName = rdbase + "list-query.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);

    TEST_ASSERT(m->getNumAtoms() == 6);
    std::string sma = MolToSmarts(*m);
    TEST_ASSERT(sma == "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6,#7,#15]:1");

    std::string mb = MolToMolBlock(*m, true, -1, true, true);
    delete m;
    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 6);
    sma = MolToSmarts(*m);
    TEST_ASSERT(sma == "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6,#7,#15]:1");
    delete m;
  }

  {
    // not atom list
    fName = rdbase + "not-list-query.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);

    TEST_ASSERT(m->getNumAtoms() == 4);
    std::string sma = MolToSmarts(*m);
    TEST_ASSERT(sma == "[#6]-[#6](-[#6])=[!#7&!#8]");

    std::string mb = MolToMolBlock(*m, true, -1, true, true);
    delete m;

    m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 4);
    sma = MolToSmarts(*m);
    TEST_ASSERT(sma == "[#6]-[#6](-[#6])=[!#7&!#8]");
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIssue2963522() {
  BOOST_LOG(rdInfoLog) << " Testing issue 2963522 " << std::endl;
  {
    std::string smiles = "CC=CC";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(1)->getStereo() == Bond::STEREONONE);

    auto *conf = new Conformer(m->getNumAtoms());
    conf->setAtomPos(0, RDGeom::Point3D(-1, 1, 0));
    conf->setAtomPos(1, RDGeom::Point3D(0, 1, 0));
    conf->setAtomPos(2, RDGeom::Point3D(0, -1, 0));
    conf->setAtomPos(3, RDGeom::Point3D(1, -1, 0));
    m->addConformer(conf, true);

    std::string molBlock = MolToMolBlock(*m);
    delete m;
    m = MolBlockToMol(molBlock);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(1)->getStereo() == Bond::STEREOANY);

    delete m;
  }

  {
    std::string smiles = "C/C=C\\C";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);

    auto *conf = new Conformer(m->getNumAtoms());
    conf->setAtomPos(0, RDGeom::Point3D(-1, 1, 0));
    conf->setAtomPos(1, RDGeom::Point3D(0, 1, 0));
    conf->setAtomPos(2, RDGeom::Point3D(0, -1, 0));
    conf->setAtomPos(3, RDGeom::Point3D(-1, -1, 0));
    m->addConformer(conf, true);

    std::string molBlock = MolToMolBlock(*m);
    delete m;
    m = MolBlockToMol(molBlock);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);

    delete m;
  }

  {
    std::string smiles = "C/C=C/C";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(1)->getStereo() == Bond::STEREOE);

    auto *conf = new Conformer(m->getNumAtoms());
    conf->setAtomPos(0, RDGeom::Point3D(-1, 1, 0));
    conf->setAtomPos(1, RDGeom::Point3D(0, 1, 0));
    conf->setAtomPos(2, RDGeom::Point3D(0, -1, 0));
    conf->setAtomPos(3, RDGeom::Point3D(1, -1, 0));
    m->addConformer(conf, true);

    std::string molBlock = MolToMolBlock(*m);
    delete m;
    m = MolBlockToMol(molBlock);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(1)->getStereo() == Bond::STEREOE);

    delete m;
  }

  {
    std::string smiles = "C1C=CC1";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(1)->getStereo() == Bond::STEREONONE);

    auto *conf = new Conformer(m->getNumAtoms());
    conf->setAtomPos(0, RDGeom::Point3D(-1, 1, 0));
    conf->setAtomPos(1, RDGeom::Point3D(0, 1, 0));
    conf->setAtomPos(2, RDGeom::Point3D(0, -1, 0));
    conf->setAtomPos(3, RDGeom::Point3D(-1, -1, 0));
    m->addConformer(conf, true);

    std::string molBlock = MolToMolBlock(*m);
    delete m;
    m = MolBlockToMol(molBlock);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(1)->getStereo() == Bond::STEREONONE);

    delete m;
  }

  {
    // this was issue 3009756:
    std::string smiles = "CC(=O)C";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(1)->getStereo() == Bond::STEREONONE);

    auto *conf = new Conformer(m->getNumAtoms());
    conf->setAtomPos(0, RDGeom::Point3D(-1, 0, 0));
    conf->setAtomPos(1, RDGeom::Point3D(0, 0, 0));
    conf->setAtomPos(2, RDGeom::Point3D(0, 1, 0));
    conf->setAtomPos(3, RDGeom::Point3D(1, 0, 0));
    m->addConformer(conf, true);

    std::string molBlock = MolToMolBlock(*m);
    delete m;
    m = MolBlockToMol(molBlock);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(1)->getStereo() == Bond::STEREONONE);

    delete m;
  }
  {
    // this was issue 3009756:
    std::string smiles = "CC(=C)Cl";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(1)->getStereo() == Bond::STEREONONE);

    auto *conf = new Conformer(m->getNumAtoms());
    conf->setAtomPos(0, RDGeom::Point3D(-1, 0, 0));
    conf->setAtomPos(1, RDGeom::Point3D(0, 0, 0));
    conf->setAtomPos(2, RDGeom::Point3D(0, 1, 0));
    conf->setAtomPos(3, RDGeom::Point3D(1, 0, 0));
    m->addConformer(conf, true);

    std::string molBlock = MolToMolBlock(*m);
    delete m;
    m = MolBlockToMol(molBlock);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(1)->getStereo() == Bond::STEREONONE);

    delete m;
  }

  BOOST_LOG(rdInfoLog) << " done" << std::endl;
}

void testIssue3073163() {
  BOOST_LOG(rdInfoLog) << " Testing issue 3073163 " << std::endl;
  {
    std::string smiles = "C[2H]";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    smiles = "[2#1]";
    RWMol *p = SmartsToMol(smiles);
    TEST_ASSERT(p);

    MatchVectType mv;
    TEST_ASSERT(SubstructMatch(*m, *p, mv));
    std::string mb = MolToMolBlock(*m);
    // std::cerr<<"mb:\n"<<mb<<"----\n";
    RWMol *m2 = MolBlockToMol(mb);
    TEST_ASSERT(m2);
    // std::cerr<<"  mol: "<<MolToSmiles(*m,true)<<std::endl;
    // std::cerr<<"  mol2: "<<MolToSmiles(*m2,true)<<std::endl;
    TEST_ASSERT(SubstructMatch(*m2, *p, mv));

    delete m2;

    delete m;
    delete p;
  }

  BOOST_LOG(rdInfoLog) << " done" << std::endl;
}

void testIssue3154208() {
  BOOST_LOG(rdInfoLog) << " Testing Issue3154208 (a large mol failure)"
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/largemol.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 476);
    TEST_ASSERT(m->getNumBonds() == 531);
    std::cerr << "generating smiles" << std::endl;
    std::string smiles = MolToSmiles(*m, false, false, -1, false);
    std::cerr << "smiles: " << smiles << std::endl;

    std::cerr << "converting back" << std::endl;
    RWMol *m2 = SmilesToMol(smiles);
    TEST_ASSERT(m2);
    TEST_ASSERT(m2->getNumAtoms() == 476);
    TEST_ASSERT(m2->getNumBonds() == 531);
    MatchVectType mv;
    std::cerr << "check isomorphism" << std::endl;
    TEST_ASSERT(SubstructMatch(*m, *m2, mv));
    TEST_ASSERT(SubstructMatch(*m2, *m, mv));
    delete m2;

#if 1
    BOOST_LOG(rdInfoLog) << "Large molecule canonical smiles test" << std::endl;
    std::string csmiles = MolToSmiles(*m);
    for (unsigned int i = 0; i < 50; ++i) {
      if (!(i % 10)) {
        BOOST_LOG(rdInfoLog) << "Iteration: " << i + 1 << " of 50" << std::endl;
      }
      std::string nsmiles = MolToSmiles(*m, false, false, 2 * i, false);
      RWMol *nm = SmilesToMol(nsmiles);
      TEST_ASSERT(nm);
      TEST_ASSERT(nm->getNumAtoms() == 476);
      TEST_ASSERT(nm->getNumBonds() == 531);
      nsmiles = MolToSmiles(*m);
      if (nsmiles != csmiles) {
        std::cerr << "MISMATCH:\n" << nsmiles << "\n" << csmiles << "\n";
      }
      TEST_ASSERT(nsmiles == csmiles);
      delete nm;
    }
#endif
    delete m;
  }

  BOOST_LOG(rdInfoLog) << " done" << std::endl;
}

void testIssue3228150() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Issue 3228150: round-trip stereochemistry failure"
                       << std::endl;

  std::string rdbase = getenv("RDBASE");
  {
    std::string fName;
    RWMol *m;
    fName = rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3228150.sdf";
    m = MolFileToMol(fName);
    TEST_ASSERT(m);

    TEST_ASSERT(m->getBondWithIdx(0)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo() == Bond::STEREOZ);
    TEST_ASSERT(m->getBondWithIdx(2)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(2)->getStereo() == Bond::STEREOZ);
    std::string smi1 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi1 << std::endl;
    m->clearComputedProps();
    m->updatePropertyCache();
    std::string smi2 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);

    delete m;
  }
  {
    std::string fName;
    RWMol *m;
    fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3228150.full.sdf";
    m = MolFileToMol(fName);
    TEST_ASSERT(m);

    TEST_ASSERT(m->getBondWithIdx(2)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(2)->getStereo() == Bond::STEREOZ);
    TEST_ASSERT(m->getBondWithIdx(4)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(4)->getStereo() == Bond::STEREOZ);
    TEST_ASSERT(m->getBondWithIdx(6)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(6)->getStereo() == Bond::STEREOZ);
    TEST_ASSERT(m->getBondWithIdx(8)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(8)->getStereo() == Bond::STEREOZ);

    std::string smi1 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi1 << std::endl;
    MolOps::assignStereochemistry(*m, true, true);
    TEST_ASSERT(m->getBondWithIdx(2)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(2)->getStereo() == Bond::STEREOZ);
    TEST_ASSERT(m->getBondWithIdx(4)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(4)->getStereo() == Bond::STEREOZ);
    TEST_ASSERT(m->getBondWithIdx(6)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(6)->getStereo() == Bond::STEREOZ);
    TEST_ASSERT(m->getBondWithIdx(8)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(8)->getStereo() == Bond::STEREOZ);

    std::string smi2 = MolToSmiles(*m, true);
    smi2 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);

    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIssue3313540() {
  BOOST_LOG(rdInfoLog) << "testing writing mol file R-groups (issue 3313540)"
                       << std::endl;

  std::string rdbase = getenv("RDBASE");
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/rgroups1.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    unsigned int idx;

    TEST_ASSERT(
        m->getAtomWithIdx(3)->hasProp(common_properties::_MolFileRLabel));
    m->getAtomWithIdx(3)->getProp(common_properties::_MolFileRLabel, idx);
    TEST_ASSERT(idx == 2);

    TEST_ASSERT(
        m->getAtomWithIdx(4)->hasProp(common_properties::_MolFileRLabel));
    m->getAtomWithIdx(4)->getProp(common_properties::_MolFileRLabel, idx);
    TEST_ASSERT(idx == 1);

    std::string mb = MolToMolBlock(*m);
    RWMol *m2 = MolBlockToMol(mb);
    TEST_ASSERT(m2);

    TEST_ASSERT(
        m2->getAtomWithIdx(3)->hasProp(common_properties::_MolFileRLabel));
    m2->getAtomWithIdx(3)->getProp(common_properties::_MolFileRLabel, idx);
    TEST_ASSERT(idx == 2);

    TEST_ASSERT(
        m2->getAtomWithIdx(4)->hasProp(common_properties::_MolFileRLabel));
    m2->getAtomWithIdx(4)->getProp(common_properties::_MolFileRLabel, idx);
    TEST_ASSERT(idx == 1);

    delete m;
    delete m2;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIssue3359739() {
  // basic writing test
  BOOST_LOG(rdInfoLog) << " ----------> Test issue 3359739 " << std::endl;

  std::string smi = "[C]C";
  RWMol *m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 3);

  std::string molBlock = MolToMolBlock(*m);
  delete m;
  m = MolBlockToMol(molBlock);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 3);
  delete m;

  BOOST_LOG(rdInfoLog) << " Finished <---------- " << std::endl;
}

void testIssue3374639() {
  // basic writing test
  BOOST_LOG(rdInfoLog) << " ----------> Test issue 3374639 " << std::endl;

  std::string rdbase = getenv("RDBASE");
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3374639.2.mol";
    RWMol *m = MolFileToMol(fName);

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    std::string cip;
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    delete m;
  }
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3374639.1.mol";
    RWMol *m = MolFileToMol(fName);

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    std::string cip;
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    delete m;
  }
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3374639.full.mol";
    RWMol *m = MolFileToMol(fName);

    TEST_ASSERT(m->getAtomWithIdx(16)->hasProp(common_properties::_CIPCode));
    std::string cip;
    m->getAtomWithIdx(16)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    delete m;
  }

  BOOST_LOG(rdInfoLog) << " Finished <---------- " << std::endl;
}

void testThreeCoordinateChirality() {
  // basic writing test
  BOOST_LOG(rdInfoLog) << " ----------> Test three-coordinate chirality "
                       << std::endl;

  std::string rdbase = getenv("RDBASE");
  {
    std::string fName =
        rdbase +
        "/Code/GraphMol/FileParsers/test_data/three_coordinate_chirality.1.mol";
    RWMol *m = MolFileToMol(fName);

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    std::string cip;
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    delete m;
  }

  {
    std::string fName =
        rdbase +
        "/Code/GraphMol/FileParsers/test_data/three_coordinate_chirality.2.mol";
    RWMol *m = MolFileToMol(fName);

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    std::string cip;
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    delete m;
  }

  {
    std::string fName =
        rdbase +
        "/Code/GraphMol/FileParsers/test_data/three_coordinate_chirality.3.mol";
    RWMol *m = MolFileToMol(fName);

    TEST_ASSERT(!m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    delete m;
  }

  {
    std::string fName =
        rdbase +
        "/Code/GraphMol/FileParsers/test_data/three_coordinate_chirality.4.mol";
    RWMol *m = MolFileToMol(fName);

    TEST_ASSERT(!m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    delete m;
  }

  {
    std::string fName =
        rdbase +
        "/Code/GraphMol/FileParsers/test_data/three_coordinate_chirality.5.mol";
    RWMol *m = MolFileToMol(fName);

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    std::string cip;
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    delete m;
  }
  {
    std::string fName =
        rdbase +
        "/Code/GraphMol/FileParsers/test_data/three_coordinate_chirality.6.mol";
    RWMol *m = MolFileToMol(fName);

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    std::string cip;
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    delete m;
  }

  BOOST_LOG(rdInfoLog) << " Finished <---------- " << std::endl;
}

void testIssue3375647() {
  // basic writing test
  BOOST_LOG(rdInfoLog) << " ----------> Test issue 3375647 " << std::endl;

  std::string rdbase = getenv("RDBASE");
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3375647.1.mol";
    RWMol *m = MolFileToMol(fName);

    TEST_ASSERT(m->getBondBetweenAtoms(2, 11)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(2, 11)->getBondDir() !=
                Bond::EITHERDOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(2, 11)->getStereo() == Bond::STEREOZ);
    delete m;
  }
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3375647.2.mol";
    RWMol *m = MolFileToMol(fName);

    TEST_ASSERT(m->getBondBetweenAtoms(2, 11)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(2, 11)->getBondDir() !=
                Bond::EITHERDOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(2, 11)->getStereo() == Bond::STEREOE);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << " Finished <---------- " << std::endl;
}

void testIssue3375684() {
  // basic writing test
  BOOST_LOG(rdInfoLog) << " ----------> Test issue 3375684 " << std::endl;

  std::string rdbase = getenv("RDBASE");
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3375684.1.mol";
    RWMol *m = MolFileToMol(fName);

    TEST_ASSERT(m->getBondBetweenAtoms(6, 7)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(6, 7)->getStereo() == Bond::STEREOANY);
    delete m;
  }
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3375684.2.mol";
    RWMol *m = MolFileToMol(fName);

    TEST_ASSERT(m->getBondBetweenAtoms(3, 9)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(3, 9)->getBondDir() !=
                Bond::EITHERDOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(3, 9)->getStereo() == Bond::STEREOE);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << " Finished <---------- " << std::endl;
}

void testChiralPhosphorous() {
  // basic writing test
  BOOST_LOG(rdInfoLog) << " ----------> Test handling of chiral phosphorous "
                       << std::endl;

  std::string rdbase = getenv("RDBASE");
  {
    std::string fName =
        rdbase +
        "/Code/GraphMol/FileParsers/test_data/chiral_phosphorous.1.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);

    TEST_ASSERT(m->getAtomWithIdx(5)->hasProp(common_properties::_CIPCode));
    std::string cip;
    m->getAtomWithIdx(5)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");
    delete m;
  }
  {
    std::string fName =
        rdbase +
        "/Code/GraphMol/FileParsers/test_data/chiral_phosphorous.2.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);

    TEST_ASSERT(m->getAtomWithIdx(5)->hasProp(common_properties::_CIPCode));
    std::string cip;
    m->getAtomWithIdx(5)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    delete m;
  }
  {
    std::string fName =
        rdbase +
        "/Code/GraphMol/FileParsers/test_data/chiral_phosphorous.3.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);

    TEST_ASSERT(!m->getAtomWithIdx(5)->hasProp(common_properties::_CIPCode));
    delete m;
  }
  {
    std::string fName =
        rdbase +
        "/Code/GraphMol/FileParsers/test_data/chiral_phosphorous.4.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);

    TEST_ASSERT(!m->getAtomWithIdx(5)->hasProp(common_properties::_CIPCode));
    delete m;
  }
  BOOST_LOG(rdInfoLog) << " Finished <---------- " << std::endl;
}

void testIssue3392107() {
  // basic writing test
  BOOST_LOG(rdInfoLog) << " ----------> Test issue 3392107 " << std::endl;

  std::string rdbase = getenv("RDBASE");
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3392107.1.mol";
    RWMol *m = MolFileToMol(fName);

    std::string smi;
    MatchVectType mv;
    RWMol *m2;

    smi = "C1CCCCC1";
    m2 = SmilesToMol(smi);
    TEST_ASSERT(m2);
    TEST_ASSERT(SubstructMatch(*m2, *m, mv));
    TEST_ASSERT(mv.size() == 6);
    delete m2;

    smi = "C1CCCCN1";
    m2 = SmilesToMol(smi);
    TEST_ASSERT(m2);
    TEST_ASSERT(SubstructMatch(*m2, *m, mv));
    TEST_ASSERT(mv.size() == 6);
    delete m2;

    smi = "C1CCNCN1";
    m2 = SmilesToMol(smi);
    TEST_ASSERT(m2);
    TEST_ASSERT(SubstructMatch(*m2, *m, mv));
    TEST_ASSERT(mv.size() == 6);
    delete m2;

    smi = "C1NCNCN1";
    m2 = SmilesToMol(smi);
    TEST_ASSERT(m2);
    TEST_ASSERT(!SubstructMatch(*m2, *m, mv));
    delete m2;

    delete m;
  }
  BOOST_LOG(rdInfoLog) << " Finished <---------- " << std::endl;
}

void testIssue3432136() {
  BOOST_LOG(rdInfoLog) << " ----------> Test issue 3432136 " << std::endl;

  std::string rdbase = getenv("RDBASE");
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3432136_1.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(!m);
  }
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3432136_1.v3k.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(!m);
  }
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3432136_2.v3k.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    delete m;
  }
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3432136_2.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << " Finished <---------- " << std::endl;
}

void testIssue3477283() {
  BOOST_LOG(rdInfoLog) << " ----------> Test issue 3477283 " << std::endl;

  std::string rdbase = getenv("RDBASE");
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3477283.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << " Finished <---------- " << std::endl;
}

void testIssue3484552() {
  BOOST_LOG(rdInfoLog) << " ----------> Test issue 3484552 " << std::endl;

  {
    std::string smi = "C[13CH3]";
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getMass() > 12.999);
    std::string molBlock = MolToMolBlock(*m);
    TEST_ASSERT(molBlock.find("M  ISO") != std::string::npos);
    delete m;
    m = MolBlockToMol(molBlock);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getMass() > 12.999);

    delete m;
  }
  BOOST_LOG(rdInfoLog) << " Finished <---------- " << std::endl;
}

void testIssue3514824() {
  BOOST_LOG(rdInfoLog) << " ----------> Test issue 3514824 " << std::endl;

  std::string rdbase = getenv("RDBASE");
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3514824.2.mol";
    RWMol *m = MolFileToMol(fName, false);
    TEST_ASSERT(m);
    m->updatePropertyCache();
    MolOps::findSSSR(*m);
    TEST_ASSERT(m->getRingInfo());
    TEST_ASSERT(m->getRingInfo()->isInitialized());
    TEST_ASSERT(m->getRingInfo()->numRings() == 6);
    delete m;
  }
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3514824.mol";
    RWMol *m = MolFileToMol(fName, false);
    TEST_ASSERT(m);
    m->updatePropertyCache();
    MolOps::findSSSR(*m);
    TEST_ASSERT(m->getRingInfo());
    TEST_ASSERT(m->getRingInfo()->isInitialized());
    TEST_ASSERT(m->getRingInfo()->numRings() == 8);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << " Finished <---------- " << std::endl;
}

void testIssue3525799() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Issue 3525799: bad smiles for r groups" << std::endl;

  std::string rdbase = getenv("RDBASE");
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3525799.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    std::string smiles = MolToSmiles(*m, true);
    std::cerr << "smiles: " << smiles << std::endl;
    TEST_ASSERT(smiles ==
                "[1*]c1c([2*])c([3*])c([4*])c(-c2c([9*])oc3c([8*])c([7*])c([6*]"
                ")c([5*])c3c2=O)c1[10*]");
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIssue3557675() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Bad issue 3557676: handling of D and T in CTABs"
                       << std::endl;

  std::string rdbase = getenv("RDBASE");
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/D_in_CTAB.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 3);
    TEST_ASSERT(m->getAtomWithIdx(2)->getAtomicNum() == 1);
    TEST_ASSERT(m->getAtomWithIdx(2)->getIsotope() == 2);
    delete m;
  }
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/T_in_CTAB.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 3);
    TEST_ASSERT(m->getAtomWithIdx(2)->getAtomicNum() == 1);
    TEST_ASSERT(m->getAtomWithIdx(2)->getIsotope() == 3);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testSkipLines() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing skip lines in CTABs" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase + "/Code/GraphMol/FileParsers/test_data/SkipLines.sdf";
  RWMol *m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 1);
  delete m;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIssue269() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Bad issue 269: handling of bad atom symbols in CTABs"
                       << std::endl;

  std::string rdbase = getenv("RDBASE");
  {
    // since the new elements were added, the original version of this no longer
    // fails. The test input file has been updated to still have an atomic
    // symbol that is not recognized. We'll be ok until Mv is an element. :-)
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue269.mol";
    RWMol *m = nullptr;
    try {
      m = MolFileToMol(fName);
    } catch (...) {
    }
    TEST_ASSERT(!m);
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMolFileChiralFlag() {
  BOOST_LOG(rdInfoLog) << "testing handling of chiral flags" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  // SF.Net Issue1603923: problems with multiple chg lines
  {
    RWMol *m1;
    std::string fName;
    fName = rdbase + "chiral_flag.mol";
    m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->hasProp(common_properties::_MolFileChiralFlag));
    unsigned int cflag;
    m1->getProp(common_properties::_MolFileChiralFlag, cflag);
    TEST_ASSERT(cflag == 1);
    delete m1;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMolFileTotalValence() {
  BOOST_LOG(rdInfoLog) << "testing handling of mol file valence flags"
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  {
    RWMol *m1;
    std::string fName;
    fName = rdbase + "Na.mol";
    m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getNumAtoms() == 1);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNoImplicit());
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumExplicitHs() == 0);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
    delete m1;
  }
  {
    RWMol *m1;
    std::string fName;
    fName = rdbase + "CH.mol";
    m1 = MolFileToMol(fName);  // cxsmi = "[CH] |(-1.2375,2.1509,),^5:0|"
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getNumAtoms() == 1);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNoImplicit());
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumExplicitHs() == 1);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumRadicalElectrons() == 3);

    delete m1;
  }
  {
    RWMol *m1;
    std::string fName;
    fName = rdbase + "CH2.mol";
    m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getNumAtoms() == 1);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNoImplicit());
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumExplicitHs() == 2);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumRadicalElectrons() == 2);
    delete m1;
  }
  {
    RWMol *m1;
    std::string fName;
    fName = rdbase + "CH3.mol";
    m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getNumAtoms() == 1);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNoImplicit());
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumExplicitHs() == 3);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
    delete m1;
  }
  {
    // make sure we get it for v3k mol blocks too:
    RWMol *m1;
    std::string fName;
    fName = rdbase + "CH.v3k.mol";
    m1 = MolFileToMol(fName);  // cxsmiles = "[CH] |(0,0,),^5:0|"
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getNumAtoms() == 1);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNoImplicit());
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumExplicitHs() == 1);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getNumRadicalElectrons() == 3);

    delete m1;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub88() {
  BOOST_LOG(rdInfoLog)
      << "testing github issue 88: M  END not being read from V3K ctabs"
      << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  {
    std::string fName;
    fName = rdbase + "github88.v3k.mol";
    bool ok = false;
    try {
      MolFileToMol(fName);
    } catch (FileParseException &) {
      ok = true;
    }
    TEST_ASSERT(ok);
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub82() {
  BOOST_LOG(rdInfoLog) << "testing github issue 82: stereochemistry only "
                          "perceived if sanitization is done"
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  {
    std::string fName;
    fName = rdbase + "github82.1.mol";
    ROMol *m;
    m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(2)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    TEST_ASSERT(m->getAtomWithIdx(3)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    TEST_ASSERT(m->getAtomWithIdx(4)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    delete m;

    m = MolFileToMol(fName, true, false);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(2)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    TEST_ASSERT(m->getAtomWithIdx(3)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    TEST_ASSERT(m->getAtomWithIdx(4)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    delete m;

    m = MolFileToMol(fName, false, false);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(2)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    TEST_ASSERT(m->getAtomWithIdx(3)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    TEST_ASSERT(m->getAtomWithIdx(4)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMolFileWithHs() {
  BOOST_LOG(rdInfoLog) << "testing impact of Hs in mol files on stereochemistry"
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  {
    std::string fName;
    fName = rdbase + "chiral_3h.mol";
    ROMol *m;
    m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(3)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    delete m;

    m = MolFileToMol(fName, true, false);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(3)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    delete m;

    m = MolFileToMol(fName, false, false);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(3)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMolFileWithRxn() {
  BOOST_LOG(rdInfoLog) << "testing reading reactions in mol files" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  {
    std::string fName;
    fName = rdbase + "rxn1.mol";
    ROMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 18);
    TEST_ASSERT(m->getNumBonds() == 16);

    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::molRxnRole));
    TEST_ASSERT(
        m->getAtomWithIdx(0)->getProp<int>(common_properties::molRxnRole) == 1);
    TEST_ASSERT(
        m->getAtomWithIdx(0)->hasProp(common_properties::molRxnComponent));
    TEST_ASSERT(m->getAtomWithIdx(0)->getProp<int>(
                    common_properties::molRxnComponent) == 1);

    TEST_ASSERT(m->getAtomWithIdx(17)->hasProp(common_properties::molRxnRole));
    TEST_ASSERT(m->getAtomWithIdx(17)->getProp<int>(
                    common_properties::molRxnRole) == 2);
    TEST_ASSERT(
        m->getAtomWithIdx(17)->hasProp(common_properties::molRxnComponent));
    TEST_ASSERT(m->getAtomWithIdx(17)->getProp<int>(
                    common_properties::molRxnComponent) == 3);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testPDBFile() {
  BOOST_LOG(rdInfoLog) << "testing reading pdb files" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  {
    std::string fName;
    fName = rdbase + "1CRN.pdb";
    ROMol *m = PDBFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 327);
    TEST_ASSERT(m->getNumBonds() == 337);
    TEST_ASSERT(m->getAtomWithIdx(0)->getMonomerInfo());
    TEST_ASSERT(m->getAtomWithIdx(0)->getMonomerInfo()->getMonomerType() ==
                AtomMonomerInfo::PDBRESIDUE);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(0)->getMonomerInfo())
                    ->getSerialNumber() == 1);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(0)->getMonomerInfo())
                    ->getResidueNumber() == 1);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(9)->getMonomerInfo())
                    ->getSerialNumber() == 10);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(9)->getMonomerInfo())
                    ->getResidueNumber() == 2);

    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(0)->getMonomerInfo())
                    ->getName() == " N  ");
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(0)->getMonomerInfo())
                    ->getResidueName() == "THR");
    TEST_ASSERT(feq(static_cast<AtomPDBResidueInfo *>(
                        m->getAtomWithIdx(0)->getMonomerInfo())
                        ->getTempFactor(),
                    13.79));
    TEST_ASSERT(m->getNumConformers() == 1);
    TEST_ASSERT(feq(m->getConformer().getAtomPos(0).x, 17.047));
    TEST_ASSERT(feq(m->getConformer().getAtomPos(0).y, 14.099));
    TEST_ASSERT(feq(m->getConformer().getAtomPos(0).z, 3.625));

    std::string pkl;
    MolPickler::pickleMol(*m, pkl);
    RWMol m2(pkl);
    for (const auto atom : m->atoms()) {
      const auto atom2 = m2.getAtomWithIdx(atom->getIdx());
      auto info1 = static_cast<AtomPDBResidueInfo *>(atom->getMonomerInfo());
      auto info2 = static_cast<AtomPDBResidueInfo *>(atom2->getMonomerInfo());
      // this is awkward because operator== isn't defined for AtomPDBResidueInfo
      // yet
      TEST_ASSERT(info1->getName() == info2->getName());
      TEST_ASSERT(info1->getChainId() == info2->getChainId());
      TEST_ASSERT(info1->getInsertionCode() == info2->getInsertionCode());
      TEST_ASSERT(info1->getTempFactor() == info2->getTempFactor());
      TEST_ASSERT(info1->getMonomerType() == info2->getMonomerType());
      TEST_ASSERT(info1->getResidueName() == info2->getResidueName());
      TEST_ASSERT(info1->getResidueNumber() == info2->getResidueNumber());
      TEST_ASSERT(info1->getSerialNumber() == info2->getSerialNumber());
      TEST_ASSERT(info1->getOccupancy() == info2->getOccupancy());
    }

    std::string mb = MolToPDBBlock(*m);
    delete m;
    m = PDBBlockToMol(mb);
    TEST_ASSERT(m->getNumAtoms() == 327);
    TEST_ASSERT(m->getNumBonds() == 337);
    TEST_ASSERT(m->getAtomWithIdx(0)->getMonomerInfo());
    TEST_ASSERT(m->getAtomWithIdx(0)->getMonomerInfo()->getMonomerType() ==
                AtomMonomerInfo::PDBRESIDUE);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(0)->getMonomerInfo())
                    ->getSerialNumber() == 1);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(0)->getMonomerInfo())
                    ->getResidueNumber() == 1);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(9)->getMonomerInfo())
                    ->getSerialNumber() == 10);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(9)->getMonomerInfo())
                    ->getResidueNumber() == 2);

    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(0)->getMonomerInfo())
                    ->getName() == " N  ");
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(0)->getMonomerInfo())
                    ->getResidueName() == "THR");
    TEST_ASSERT(feq(static_cast<AtomPDBResidueInfo *>(
                        m->getAtomWithIdx(0)->getMonomerInfo())
                        ->getTempFactor(),
                    13.79));
    TEST_ASSERT(m->getNumConformers() == 1);
    TEST_ASSERT(feq(m->getConformer().getAtomPos(0).x, 17.047));
    TEST_ASSERT(feq(m->getConformer().getAtomPos(0).y, 14.099));
    TEST_ASSERT(feq(m->getConformer().getAtomPos(0).z, 3.625));

    // test adding hydrogens
    ROMol *nm = MolOps::addHs(*m, false, false, nullptr, true);
    auto *info =
        (AtomPDBResidueInfo *)(nm->getAtomWithIdx(nm->getNumAtoms() - 1)
                                   ->getMonomerInfo());
    TEST_ASSERT(info->getMonomerType() == AtomMonomerInfo::PDBRESIDUE);
    TEST_ASSERT(info->getName() == " H7 ");
    TEST_ASSERT(info->getResidueName() == "ASN");
    delete nm;
    delete m;
  }

  {
    std::string fName;
    fName = rdbase + "2FVD.pdb";
    ROMol *m = PDBFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 2501);
    TEST_ASSERT(m->getNumBonds() == 2383);
    TEST_ASSERT(m->getAtomWithIdx(0)->getMonomerInfo());
    TEST_ASSERT(m->getAtomWithIdx(0)->getMonomerInfo()->getMonomerType() ==
                AtomMonomerInfo::PDBRESIDUE);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(0)->getMonomerInfo())
                    ->getSerialNumber() == 1);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(0)->getMonomerInfo())
                    ->getResidueNumber() == 1);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(0)->getMonomerInfo())
                    ->getIsHeteroAtom() == 0);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(2292)->getMonomerInfo())
                    ->getSerialNumber() == 2294);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(2292)->getMonomerInfo())
                    ->getResidueNumber() == 299);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(2292)->getMonomerInfo())
                    ->getIsHeteroAtom() == 1);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(2292)->getMonomerInfo())
                    ->getChainId() == "A");

    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(1)->getMonomerInfo())
                    ->getName() == " CA ");
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(1)->getMonomerInfo())
                    ->getResidueName() == "MET");
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(2292)->getMonomerInfo())
                    ->getName() == " N1 ");
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(2292)->getMonomerInfo())
                    ->getResidueName() == "LIA");

    std::string mb = MolToPDBBlock(*m, -1, 32);
    delete m;
    m = PDBBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 2501);
    TEST_ASSERT(m->getNumBonds() == 2383);
    TEST_ASSERT(m->getAtomWithIdx(0)->getMonomerInfo());
    TEST_ASSERT(m->getAtomWithIdx(0)->getMonomerInfo()->getMonomerType() ==
                AtomMonomerInfo::PDBRESIDUE);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(0)->getMonomerInfo())
                    ->getSerialNumber() == 1);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(0)->getMonomerInfo())
                    ->getResidueNumber() == 1);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(0)->getMonomerInfo())
                    ->getIsHeteroAtom() == 0);
    // FIX:
    // TEST_ASSERT(static_cast<AtomPDBResidueInfo
    // *>(m->getAtomWithIdx(2292)->getMonomerInfo())->getSerialNumber()==2294);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(2292)->getMonomerInfo())
                    ->getResidueNumber() == 299);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(2292)->getMonomerInfo())
                    ->getIsHeteroAtom() == 1);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(2292)->getMonomerInfo())
                    ->getChainId() == "A");

    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(1)->getMonomerInfo())
                    ->getName() == " CA ");
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(1)->getMonomerInfo())
                    ->getResidueName() == "MET");
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(2292)->getMonomerInfo())
                    ->getName() == " N1 ");
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(2292)->getMonomerInfo())
                    ->getResidueName() == "LIA");
    delete m;
  }

  {  // DNA
    std::string fName;
    fName = rdbase + "4BNA.pdb";
    ROMol *m = PDBFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumHeavyAtoms() == 602);
    TEST_ASSERT(m->getAtomWithIdx(0)->getMonomerInfo());
    TEST_ASSERT(m->getAtomWithIdx(0)->getMonomerInfo()->getMonomerType() ==
                AtomMonomerInfo::PDBRESIDUE);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(0)->getMonomerInfo())
                    ->getSerialNumber() == 1);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(1)->getMonomerInfo())
                    ->getResidueName() == " DC");
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(57)->getMonomerInfo())
                    ->getResidueName() == " DG");
    std::string mb = MolToPDBBlock(*m);
    delete m;
    m = PDBBlockToMol(mb);
    TEST_ASSERT(m->getNumHeavyAtoms() == 602);
    TEST_ASSERT(m->getAtomWithIdx(0)->getMonomerInfo());
    TEST_ASSERT(m->getAtomWithIdx(0)->getMonomerInfo()->getMonomerType() ==
                AtomMonomerInfo::PDBRESIDUE);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(1)->getMonomerInfo())
                    ->getResidueName() == " DC");
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(57)->getMonomerInfo())
                    ->getResidueName() == " DG");
    TEST_ASSERT(m->getBondBetweenAtoms(31, 33)->getBondType() ==
                Bond::AROMATIC);
    TEST_ASSERT(m->getBondBetweenAtoms(120, 114)->getBondType() ==
                Bond::AROMATIC);
    TEST_ASSERT(m->getBondBetweenAtoms(39, 38)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(138, 140)->getBondType() ==
                Bond::AROMATIC);
    delete m;
  }
  {  // RNA
    std::string fName;
    fName = rdbase + "4TNA.pdb";
    ROMol *m = PDBFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumHeavyAtoms() == 1656);
    TEST_ASSERT(m->getAtomWithIdx(0)->getMonomerInfo());
    TEST_ASSERT(m->getAtomWithIdx(0)->getMonomerInfo()->getMonomerType() ==
                AtomMonomerInfo::PDBRESIDUE);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(0)->getMonomerInfo())
                    ->getSerialNumber() == 1);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(1)->getMonomerInfo())
                    ->getResidueName() == "  G");
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(90)->getMonomerInfo())
                    ->getResidueName() == "  A");
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(197)->getMonomerInfo())
                    ->getResidueName() == "2MG");
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(197)->getMonomerInfo())
                    ->getIsHeteroAtom());
    TEST_ASSERT(m->getBondBetweenAtoms(104, 103)->getBondType() ==
                Bond::AROMATIC);
    TEST_ASSERT(m->getBondBetweenAtoms(60, 61)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(38, 37)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(148, 149)->getBondType() ==
                Bond::DOUBLE);

    std::string mb = MolToPDBBlock(*m);
    delete m;
    m = PDBBlockToMol(mb);
    TEST_ASSERT(m->getNumHeavyAtoms() == 1656);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(0)->getMonomerInfo())
                    ->getSerialNumber() == 1);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(1)->getMonomerInfo())
                    ->getResidueName() == "  G");
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(90)->getMonomerInfo())
                    ->getResidueName() == "  A");
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(197)->getMonomerInfo())
                    ->getResidueName() == "2MG");
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(197)->getMonomerInfo())
                    ->getIsHeteroAtom());

    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testSequences() {
  BOOST_LOG(rdInfoLog) << "testing reading sequences" << std::endl;
  {
    std::string seq = "CGCGAATTACCGCG";  // made up
    int flavor = 6;                      // DNA
    ROMol *m = SequenceToMol(seq, true, flavor);
    TEST_ASSERT(m);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(0)->getMonomerInfo())
                    ->getSerialNumber() == 1);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(1)->getMonomerInfo())
                    ->getResidueName() == " DC");
    seq = MolToSequence(*m);
    TEST_ASSERT(seq == "CGCGAATTACCGCG");
    seq = MolToHELM(*m);
    // std::cerr << seq << std::endl;
    TEST_ASSERT(seq ==
                "RNA1{[dR](C)P.[dR](G)P.[dR](C)P.[dR](G)P.[dR](A)P.[dR](A)P.["
                "dR](T)P.[dR](T)P.[dR](A)P.[dR](C)P.[dR](C)P.[dR](G)P.[dR](C)P."
                "[dR](G)}$$$$");
    {
      std::string lseq = MolToHELM(*m);
      TEST_ASSERT(lseq == seq);

      ROMol *m2 = HELMToMol(seq);
      TEST_ASSERT(m2)
      lseq = MolToSequence(*m2);
      TEST_ASSERT(lseq == "CGCGAATTACCGCG");
      lseq = MolToHELM(*m2);
      TEST_ASSERT(lseq == seq);
      delete m2;
    }

    {
      ROMol *nm = MolOps::addHs(*m);
      TEST_ASSERT(nm);
      std::string pdb = MolToPDBBlock(*nm);
      delete nm;
      nm = PDBBlockToMol(pdb);
      TEST_ASSERT(nm);
      std::string lseq = MolToSequence(*nm);
      TEST_ASSERT(lseq == "CGCGAATTACCGCG");
      delete nm;
    }

    delete m;
  }
  {
    std::string seq = "CGCGAAUUACCGCG";  // made up
    int flavor = 2;                      // RNA
    ROMol *m = SequenceToMol(seq, true, flavor);
    TEST_ASSERT(m);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(0)->getMonomerInfo())
                    ->getSerialNumber() == 1);
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(1)->getMonomerInfo())
                    ->getResidueName() == "  C");
    seq = MolToSequence(*m);
    TEST_ASSERT(seq == "CGCGAAUUACCGCG");
    seq = MolToHELM(*m);
    TEST_ASSERT(seq ==
                "RNA1{R(C)P.R(G)P.R(C)P.R(G)P.R(A)P.R(A)P.R(U)P.R(U)P.R(A)P.R("
                "C)P.R(C)P.R(G)P.R(C)P.R(G)}$$$$");
    {
      std::string lseq = MolToHELM(*m);
      TEST_ASSERT(lseq == seq);

      ROMol *m2 = HELMToMol(seq);
      TEST_ASSERT(m2)
      lseq = MolToSequence(*m2);
      TEST_ASSERT(lseq == "CGCGAAUUACCGCG");
      lseq = MolToHELM(*m2);
      TEST_ASSERT(lseq == seq);
      delete m2;
    }
    {
      ROMol *nm = MolOps::addHs(*m);
      TEST_ASSERT(nm);
      std::string pdb = MolToPDBBlock(*nm);
      delete nm;
      nm = PDBBlockToMol(pdb);
      TEST_ASSERT(nm);
      std::string lseq = MolToSequence(*nm);
      TEST_ASSERT(lseq == "CGCGAAUUACCGCG");
      delete nm;
    }

    delete m;
  }
}

void testGithub1023() {
  BOOST_LOG(rdInfoLog) << "GetSSSR interrupted by segmentation fault"
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  {
    std::string fName;
    fName = rdbase + "github1023.pdb";
    bool sanitize = false;
    ROMol *m = PDBFileToMol(fName, sanitize);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 119);
    TEST_ASSERT(m->getNumBonds() == 399);
    MolOps::findSSSR(*m);  // this was seg faulting
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub166() {
  BOOST_LOG(rdInfoLog)
      << "testing Github 166: skipping sanitization on reading pdb files"
      << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  {
    std::string fName;
    fName = rdbase + "1CRN.pdb";
    ROMol *m = PDBFileToMol(fName, false, false);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 327);
    TEST_ASSERT(m->getNumBonds() == 337);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testZBO() {
  BOOST_LOG(rdInfoLog) << "testing ZBO parsing" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  {
    std::string fName;
    fName = rdbase + "FeCO5.mol";
    ROMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 11);
    TEST_ASSERT(m->getNumBonds() == 10);
    TEST_ASSERT(m->getBondWithIdx(0)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(1)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(2)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(6)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(7)->getBondType() == Bond::ZERO);
    delete m;
  }
  {
    std::string fName;
    fName = rdbase + "CrBz.mol";
    ROMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 7);
    TEST_ASSERT(m->getNumBonds() == 12);
    TEST_ASSERT(m->getBondWithIdx(6)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(7)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(8)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(9)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(10)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(11)->getBondType() == Bond::ZERO);

    // make sure we don't screw up aromaticity:
    TEST_ASSERT(m->getBondWithIdx(0)->getIsAromatic());
    delete m;
  }
  {
    std::string fName;
    fName = rdbase + "CrBz2.mol";
    ROMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 13);
    TEST_ASSERT(m->getNumBonds() == 24);
    TEST_ASSERT(m->getBondWithIdx(6)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(7)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(8)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(9)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(10)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(11)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(18)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(19)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(20)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(21)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(22)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getBondWithIdx(23)->getBondType() == Bond::ZERO);
    delete m;
  }
  {
    std::string fName;
    fName = rdbase + "H3BNH3.mol";
    ROMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 2);
    TEST_ASSERT(m->getNumBonds() == 1);
    TEST_ASSERT(m->getBondWithIdx(0)->getBondType() == Bond::ZERO);
    TEST_ASSERT(m->getAtomWithIdx(0)->getFormalCharge() == 0);
    TEST_ASSERT(m->getAtomWithIdx(1)->getFormalCharge() == 0);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumExplicitHs() == 3);
    TEST_ASSERT(m->getAtomWithIdx(1)->getNumExplicitHs() == 0);
    TEST_ASSERT(m->getAtomWithIdx(0)->getTotalNumHs() == 3);
    TEST_ASSERT(m->getAtomWithIdx(1)->getTotalNumHs() == 3);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub164() {
  BOOST_LOG(rdInfoLog) << "testing Github 164: problems with Xe from mol files"
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  {
    std::string fName;
    fName = rdbase + "Github164.mol";
    ROMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 3);
    TEST_ASSERT(m->getNumBonds() == 2);
    TEST_ASSERT(m->getAtomWithIdx(0)->getExplicitValence() == 2);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub194() {
  BOOST_LOG(rdInfoLog) << "testing github issue 194: bad bond types from pdb"
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  {
    std::string fName;
    fName = rdbase + "1CRN.pdb";
    ROMol *m = PDBFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 327);
    TEST_ASSERT(m->getNumBonds() == 337);
    TEST_ASSERT(m->getAtomWithIdx(0)->getMonomerInfo());
    TEST_ASSERT(m->getAtomWithIdx(0)->getMonomerInfo()->getMonomerType() ==
                AtomMonomerInfo::PDBRESIDUE);
    // the root cause: problems in SamePDBResidue:
    TEST_ASSERT(SamePDBResidue(static_cast<AtomPDBResidueInfo *>(
                                   m->getAtomWithIdx(0)->getMonomerInfo()),
                               static_cast<AtomPDBResidueInfo *>(
                                   m->getAtomWithIdx(1)->getMonomerInfo())));
    TEST_ASSERT(SamePDBResidue(static_cast<AtomPDBResidueInfo *>(
                                   m->getAtomWithIdx(0)->getMonomerInfo()),
                               static_cast<AtomPDBResidueInfo *>(
                                   m->getAtomWithIdx(2)->getMonomerInfo())));
    TEST_ASSERT(!SamePDBResidue(static_cast<AtomPDBResidueInfo *>(
                                    m->getAtomWithIdx(0)->getMonomerInfo()),
                                static_cast<AtomPDBResidueInfo *>(
                                    m->getAtomWithIdx(9)->getMonomerInfo())));

    // the symptom, bond orders:
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(2)->getMonomerInfo())
                    ->getName() == " C  ");
    TEST_ASSERT(static_cast<AtomPDBResidueInfo *>(
                    m->getAtomWithIdx(3)->getMonomerInfo())
                    ->getName() == " O  ");
    TEST_ASSERT(m->getBondBetweenAtoms(2, 3));
    TEST_ASSERT(m->getBondBetweenAtoms(2, 3)->getBondType() == Bond::DOUBLE);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub196() {
  BOOST_LOG(rdInfoLog)
      << "testing github issue 196: left justified bond topology" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  {
    std::string fName;
    fName = rdbase + "github196.mol";
    ROMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 19);
    TEST_ASSERT(m->getNumBonds() == 20);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub191() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing github issue 191: "
                          "wavy bonds to Hs should affect attached double bond "
                          "stereochemistry."
                       << std::endl;
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/FileParsers/test_data/";
    RWMol *m = MolFileToMol(pathName + "github191.1.mol");
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo() == Bond::STEREOE);
    std::string smi = MolToSmiles(*m, true);
    TEST_ASSERT(smi == "C/C=C/C");
    delete m;
  }
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/FileParsers/test_data/";
    RWMol *m = MolFileToMol(pathName + "github191.2.mol");
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo() == Bond::STEREOANY);
    std::string smi = MolToSmiles(*m, true);
    TEST_ASSERT(smi == "CC=CC");
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testGithub210() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github 210: flag possible stereocenters "
                          "when calling assignStereochemistry()"
                       << std::endl;
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/FileParsers/test_data/";
    RWMol *m = MolFileToMol(pathName + "github210.mol");
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    TEST_ASSERT(
        m->getAtomWithIdx(4)->hasProp(common_properties::_ChiralityPossible));
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

namespace {
std::string getResidue(const ROMol &, const Atom *at) {
  if (at->getMonomerInfo()->getMonomerType() != AtomMonomerInfo::PDBRESIDUE) {
    return "";
  }
  return static_cast<const AtomPDBResidueInfo *>(at->getMonomerInfo())
      ->getResidueName();
}
}  // namespace
void testPDBResidues() {
  BOOST_LOG(rdInfoLog) << "testing splitting on PDB residues" << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  {
    std::string fName;
    fName = rdbase + "2NW4.pdb";
    ROMol *m = PDBFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getMonomerInfo()->getMonomerType() ==
                AtomMonomerInfo::PDBRESIDUE);
    std::map<std::string, boost::shared_ptr<ROMol>> res =
        MolOps::getMolFragsWithQuery(*m, getResidue, false);

    TEST_ASSERT(res.size() == 22);
    TEST_ASSERT(res.find(std::string("8NH")) != res.end());
    TEST_ASSERT(res.find(std::string("ALA")) != res.end());
    TEST_ASSERT(res[std::string("8NH")]->getNumAtoms() == 21);

    const ROMol *lig = res[std::string("8NH")].get();
    TEST_ASSERT(lig->getNumConformers() == 1);
    TEST_ASSERT(feq(lig->getConformer().getAtomPos(0).x, 23.517));
    TEST_ASSERT(feq(lig->getConformer().getAtomPos(0).y, 5.263));
    TEST_ASSERT(feq(lig->getConformer().getAtomPos(0).z, 4.399));
    TEST_ASSERT(feq(lig->getConformer().getAtomPos(11).x, 27.589));
    TEST_ASSERT(feq(lig->getConformer().getAtomPos(11).y, -0.311));
    TEST_ASSERT(feq(lig->getConformer().getAtomPos(11).z, 3.743));
    delete m;
  }
  {
    std::string fName;
    fName = rdbase + "2NW4.pdb";
    ROMol *m = PDBFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getMonomerInfo()->getMonomerType() ==
                AtomMonomerInfo::PDBRESIDUE);
    std::vector<std::string> keep;
    keep.emplace_back("8NH");
    std::map<std::string, boost::shared_ptr<ROMol>> res =
        MolOps::getMolFragsWithQuery(*m, getResidue, false, &keep);

    TEST_ASSERT(res.size() == 1);
    TEST_ASSERT(res.find(std::string("8NH")) != res.end());
    TEST_ASSERT(res[std::string("8NH")]->getNumAtoms() == 21);

    const ROMol *lig = res[std::string("8NH")].get();
    TEST_ASSERT(lig->getNumConformers() == 1);
    TEST_ASSERT(feq(lig->getConformer().getAtomPos(0).x, 23.517));
    TEST_ASSERT(feq(lig->getConformer().getAtomPos(0).y, 5.263));
    TEST_ASSERT(feq(lig->getConformer().getAtomPos(0).z, 4.399));
    TEST_ASSERT(feq(lig->getConformer().getAtomPos(11).x, 27.589));
    TEST_ASSERT(feq(lig->getConformer().getAtomPos(11).y, -0.311));
    TEST_ASSERT(feq(lig->getConformer().getAtomPos(11).z, 3.743));
    delete m;
  }
  {
    std::string fName;
    fName = rdbase + "2NW4.pdb";
    ROMol *m = PDBFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getMonomerInfo()->getMonomerType() ==
                AtomMonomerInfo::PDBRESIDUE);
    std::vector<std::string> keep;
    keep.emplace_back("8NH");
    std::map<std::string, boost::shared_ptr<ROMol>> res =
        MolOps::getMolFragsWithQuery(*m, getResidue, false, &keep, true);

    TEST_ASSERT(res.size() == 21);
    TEST_ASSERT(res.find(std::string("8NH")) == res.end());
    TEST_ASSERT(res.find(std::string("ALA")) != res.end());
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub337() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github 337: No double bond stereo "
                          "perception from CTABs when sanitization is turned "
                          "off"
                       << std::endl;
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/FileParsers/test_data/";
    RWMol *m = MolFileToMol(pathName + "unsanitized_stereo.mol", false);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo() == Bond::STEREONONE);
    std::string molBlock = MolToMolBlock(*m);
    TEST_ASSERT(molBlock.find("  1  2  2  0") != std::string::npos);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testGithub360() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github 360: Computed props on non-sanitized "
                          "molecule interfering with substructure matching"
                       << std::endl;
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/FileParsers/test_data/";
    RWMol *dbm = SmilesToMol("C1Cc2ccccc2CN1");
    TEST_ASSERT(dbm);
    RWMol *tmpl = MolFileToMol(pathName + "github360.mol", false);
    TEST_ASSERT(tmpl);

    MatchVectType mv;
    TEST_ASSERT(SubstructMatch(*dbm, *tmpl, mv));

    delete dbm;
    delete tmpl;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testGithub741() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github 741: Support CTABs where the second "
                          "letter in atom symbols is capitalized"
                       << std::endl;
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/FileParsers/test_data/";
    RWMol *mol = MolFileToMol(pathName + "github741.mol");
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getAtomWithIdx(1)->getSymbol() == "Br");
    TEST_ASSERT(mol->getAtomWithIdx(2)->getSymbol() == "Br");
    delete mol;
  }
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/FileParsers/test_data/";
    RWMol *mol = MolFileToMol(pathName + "github741.v3k.mol");
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getAtomWithIdx(1)->getSymbol() == "Br");
    TEST_ASSERT(mol->getAtomWithIdx(2)->getSymbol() == "Br");
    delete mol;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}
void testGithub188() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github 188: Bad E/Z assignment from CTAB"
                       << std::endl;
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/FileParsers/test_data/";
    RWMol *mol = MolFileToMol(pathName + "github188.mol");
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getBondBetweenAtoms(16, 17));
    TEST_ASSERT(mol->getBondBetweenAtoms(16, 17)->getStereo() == Bond::STEREOZ);
    TEST_ASSERT(mol->getBondBetweenAtoms(13, 14));
    TEST_ASSERT(mol->getBondBetweenAtoms(13, 14)->getStereo() == Bond::STEREOZ);

    delete mol;
  }
  {
    std::cerr << "----------------------------------" << std::endl;
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/FileParsers/test_data/";
    RWMol *mol = MolFileToMol(pathName + "github188_2.mol");
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getBondBetweenAtoms(16, 17));
    TEST_ASSERT(mol->getBondBetweenAtoms(16, 17)->getStereo() == Bond::STEREOZ);
    TEST_ASSERT(mol->getBondBetweenAtoms(13, 14));
    TEST_ASSERT(mol->getBondBetweenAtoms(13, 14)->getStereo() == Bond::STEREOZ);

    delete mol;
  }
  {
    std::cerr << "----------------------------------" << std::endl;
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/FileParsers/test_data/";
    RWMol *mol = MolFileToMol(pathName + "github192.mol");
    TEST_ASSERT(mol);
    for (unsigned int i = 0; i < mol->getNumBonds(); ++i) {
      const Bond *bnd = mol->getBondWithIdx(i);
      if (bnd->getBondType() == Bond::DOUBLE &&
          i != mol->getBondBetweenAtoms(26, 27)->getIdx() &&
          i != mol->getBondBetweenAtoms(6, 8)->getIdx()) {
        TEST_ASSERT(bnd->getStereo() == Bond::STEREOE);
      }
    }
    delete mol;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}
void testRCSBSdf() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing right-aligned elements in RCSB SDF files"
                       << std::endl;
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/FileParsers/test_data/";
  RWMol *mol = MolFileToMol(pathName + "s58_rcsb.mol");
  TEST_ASSERT(mol);
  delete mol;
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testParseCHG() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing PDB charge parsing" << std::endl;
  // BAD PDB Ligand with CHG line too long (>8) and right and mid-justified
  // symbols
  const std::string molblock_chg =
      "2D1G_DVT_B_2001\n"
      "  RCSB PDB01151500373D\n"
      "Coordinates from PDB:2D1G:B:2001 Model:1 without hydrogens\n"
      " 38 60  0  0  0  0            999 V2000\n"
      "   19.0320   93.5880   16.2640   O 0  5  0  0  0  0  0  0  0  0  0  0\n"
      "   19.6400   94.8240   15.4350  V  0  5  0  0  0  0  0  0  0  0  0  0\n"
      "   18.1700   95.2990   14.2790 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   19.6000   96.2720   16.5920   O 0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   20.3140   97.8870   15.9210   V 0  5  0  0  0  0  0  0  0  0  0  0\n"
      "   20.1880   98.9370   17.1150   O 0  5  0  0  0  0  0  0  0  0  0  0\n"
      "   21.5140   94.4720   15.5850   O 0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   22.5590   95.9150   15.0470   V 0  5  0  0  0  0  0  0  0  0  0  0\n"
      "   20.6050   96.4040   14.1390   O 0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   22.1700   97.4330   15.9940   O 0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   23.1150   97.3160   13.7160   O 0  3  0  0  0  0  0  0  0  0  0  0\n"
      "   23.9960   95.5550   15.6610   O 0  5  0  0  0  0  0  0  0  0  0  0\n"
      "   22.6470   94.8660   13.3320   O 0  3  0  0  0  0  0  0  0  0  0  0\n"
      "   23.7710   96.1210   12.1910   V 0  5  0  0  0  0  0  0  0  0  0  0\n"
      "   24.0950   97.6940   11.3470   O 0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   25.2350   95.8110   12.7760   O 0  5  0  0  0  0  0  0  0  0  0  0\n"
      "   23.5040   94.7850   10.9890   O 0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   22.0160   95.1170    9.8340   V 0  5  0  0  0  0  0  0  0  0  0  0\n"
      "   20.1580   95.5900    9.6670   O 0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   22.2300   94.0400    8.6680   O 0  5  0  0  0  0  0  0  0  0  0  0\n"
      "   22.6710   96.7240    9.1340   O 0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   20.8370   94.9790   12.6510   V 0  5  0  0  0  0  0  0  0  0  0  0\n"
      "   21.1160   93.9980   11.2840   O 0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   20.1410   93.8410   13.6860   O 0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   21.6740   96.5770   11.6070   O 0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   22.5990   98.1610   10.2810   V 0  5  0  0  0  0  0  0  0  0  0  0\n"
      "   23.1590   99.4050    9.4060   O 0  5  0  0  0  0  0  0  0  0  0  0\n"
      "   22.1350   99.1410   12.0500   O 0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   21.4440   98.0120   13.0890   V 0  5  0  0  0  0  0  0  0  0  0  0\n"
      "   21.1330   99.0140   14.4480   O 0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   19.1670   95.7180   12.0250   O 0  3  0  0  0  0  0  0  0  0  0  0\n"
      "   19.7240   97.0880   10.6480   V 0  5  0  0  0  0  0  0  0  0  0  0\n"
      "   20.7360   98.5460   10.1630   O 0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   18.2530   97.4600   10.0390   O 0  5  0  0  0  0  0  0  0  0  0  0\n"
      "   19.6430   98.1340   12.3760   O 0  3  0  0  0  0  0  0  0  0  0  0\n"
      "   18.5370   96.9390   13.5530   V 0  5  0  0  0  0  0  0  0  0  0  0\n"
      "   18.8160   98.2640   14.7770   O 0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   17.0530   97.2880   13.0100   O 0  5  0  0  0  0  0  0  0  0  0  0\n"
      "  1  2  1  0  0  0  0\n"
      "  2  3  1  0  0  0  0\n"
      "  2  4  1  0  0  0  0\n"
      "  2  7  1  0  0  0  0\n"
      "  2  9  1  0  0  0  0\n"
      "  2 24  1  0  0  0  0\n"
      "  3 36  1  0  0  0  0\n"
      "  4  5  1  0  0  0  0\n"
      "  5  6  1  0  0  0  0\n"
      "  5  9  1  0  0  0  0\n"
      "  5 10  1  0  0  0  0\n"
      "  5 30  1  0  0  0  0\n"
      "  5 37  1  0  0  0  0\n"
      "  7  8  1  0  0  0  0\n"
      "  8  9  1  0  0  0  0\n"
      "  8 10  1  0  0  0  0\n"
      "  8 11  1  0  0  0  0\n"
      "  8 12  1  0  0  0  0\n"
      "  8 13  1  0  0  0  0\n"
      "  9 22  1  0  0  0  0\n"
      "  9 29  1  0  0  0  0\n"
      "  9 36  1  0  0  0  0\n"
      " 11 14  1  0  0  0  0\n"
      " 11 29  1  0  0  0  0\n"
      " 13 14  1  0  0  0  0\n"
      " 13 22  1  0  0  0  0\n"
      " 14 15  1  0  0  0  0\n"
      " 14 16  1  0  0  0  0\n"
      " 14 17  1  0  0  0  0\n"
      " 14 28  1  0  0  0  0\n"
      " 15 26  1  0  0  0  0\n"
      " 17 18  1  0  0  0  0\n"
      " 18 19  1  0  0  0  0\n"
      " 18 20  1  0  0  0  0\n"
      " 18 21  1  0  0  0  0\n"
      " 18 23  1  0  0  0  0\n"
      " 18 28  1  0  0  0  0\n"
      " 19 32  1  0  0  0  0\n"
      " 21 26  1  0  0  0  0\n"
      " 22 23  1  0  0  0  0\n"
      " 22 24  1  0  0  0  0\n"
      " 22 28  1  0  0  0  0\n"
      " 22 31  1  0  0  0  0\n"
      " 25 26  1  0  0  0  0\n"
      " 25 29  1  0  0  0  0\n"
      " 26 27  1  0  0  0  0\n"
      " 26 28  1  0  0  0  0\n"
      " 26 33  1  0  0  0  0\n"
      " 28 29  1  0  0  0  0\n"
      " 28 32  1  0  0  0  0\n"
      " 29 30  1  0  0  0  0\n"
      " 29 35  1  0  0  0  0\n"
      " 31 32  1  0  0  0  0\n"
      " 31 36  1  0  0  0  0\n"
      " 32 33  1  0  0  0  0\n"
      " 32 34  1  0  0  0  0\n"
      " 32 35  1  0  0  0  0\n"
      " 35 36  1  0  0  0  0\n"
      " 36 37  1  0  0  0  0\n"
      " 36 38  1  0  0  0  0\n"
      "M  CHG 24   1  -1   2  -1   5  -1   6  -1   8  -1   9   4  11   1  12  "
      "-1  13   1  14  -1  16  -1  18  -1  20  -1  22  -1  26  -1  27  -1  28  "
      " 4  29  -1  31   1  32  -1  34  -1  35   1  36  -1  38  -1\n"
      "M  END\n";
  const int charges[] = {1,  -1, 2,  -1, 5,  -1, 6,  -1, 8,  -1, 9,  4,  11,
                         1,  12, -1, 13, 1,  14, -1, 16, -1, 18, -1, 20, -1,
                         22, -1, 26, -1, 27, -1, 28, 4,  29, -1, 31, 1,  32,
                         -1, 34, -1, 35, 1,  36, -1, 38, -1, 0,  0};
  // Shouldn't seg fault, throw exception or have a null mol
  bool sanitize=false;
  std::unique_ptr<RWMol> m(MolBlockToMol(molblock_chg, sanitize));
  size_t i = 0;
  while (1) {
    if (charges[i] == 0) {
      break;
    }

    TEST_ASSERT(
        m->getAtomWithIdx((unsigned int)charges[i] - 1)->getFormalCharge() ==
        charges[i + 1]);
    i += 2;
  }

  TEST_ASSERT(m);
  // Write it out in V3000 format, which makes counting the different charges
  // easier
  bool forceV3000(true);
  std::string out = MolToMolBlock(*m, true, -1, true, forceV3000);
  std::regex chg_all("CHG="), chg_m1("CHG=-1"), chg_p1("CHG=1"),
      chg_p4("CHG=4");
  TEST_ASSERT(
      std::distance(std::sregex_iterator(out.begin(), out.end(), chg_all),
                    std::sregex_iterator()) == 24);
  TEST_ASSERT(
      std::distance(std::sregex_iterator(out.begin(), out.end(), chg_m1),
                    std::sregex_iterator()) == 18);
  TEST_ASSERT(
      std::distance(std::sregex_iterator(out.begin(), out.end(), chg_p1),
                    std::sregex_iterator()) == 4);
  TEST_ASSERT(
      std::distance(std::sregex_iterator(out.begin(), out.end(), chg_p4),
                    std::sregex_iterator()) == 2);

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testMDLAtomProps() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing MDL atom properties" << std::endl;
  std::string smi = "CC";
  ROMOL_SPTR mol(SmilesToMol(smi, false, false));
  setAtomAlias(mol->getAtomWithIdx(0), "foo");
  setAtomValue(mol->getAtomWithIdx(0), "bar");
  setAtomRLabel(mol->getAtomWithIdx(0), 1);
  mol.reset(MolBlockToMol(MolToMolBlock(*mol.get())));
  TEST_ASSERT(getAtomAlias(mol->getAtomWithIdx(0)) == "foo");
  TEST_ASSERT(getAtomValue(mol->getAtomWithIdx(0)) == "bar");
  TEST_ASSERT(getAtomRLabel(mol->getAtomWithIdx(0)) == 1);
  try {
    setAtomRLabel(mol->getAtomWithIdx(0), 100);
    TEST_ASSERT(0);
  } catch (...) {
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testSupplementalSmilesLabel() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing supplemental SMILES labels" << std::endl;
  std::string smi = "C";
  ROMOL_SPTR mol(SmilesToMol(smi, false, false));
  setSupplementalSmilesLabel(mol->getAtomWithIdx(0), "xxx");
  smi = MolToSmiles(*mol.get());
  TEST_ASSERT(smi == "Cxxx");
  TEST_ASSERT(getSupplementalSmilesLabel(mol->getAtomWithIdx(0)) == "xxx");
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testGithub1034() {
  BOOST_LOG(rdInfoLog)
      << "Test github 1034: Squiggle bonds from CTABs lost post-parsing"
      << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  {  // double bond
    std::string fName;
    fName = rdbase + "github1034.1.mol";
    bool sanitize = true;
    RWMol *m = MolFileToMol(fName, sanitize);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 4);
    TEST_ASSERT(m->getBondWithIdx(0)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo() == Bond::STEREOANY);
    TEST_ASSERT(m->getBondWithIdx(1)->getBondType() == Bond::SINGLE);
    TEST_ASSERT(m->getBondWithIdx(2)->getBondType() == Bond::SINGLE);
    TEST_ASSERT(m->getBondWithIdx(2)->getBondDir() != Bond::UNKNOWN);
    TEST_ASSERT(
        m->getBondWithIdx(2)->hasProp(common_properties::_UnknownStereo));
    delete m;
  }
  {  // double bond
    std::string fName;
    fName = rdbase + "github1034.1.mol";
    bool sanitize = false;
    RWMol *m = MolFileToMol(fName, sanitize);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 4);
    TEST_ASSERT(m->getBondWithIdx(2)->getBondType() == Bond::SINGLE);
    TEST_ASSERT(m->getBondWithIdx(2)->getBondDir() != Bond::UNKNOWN);
    int explicit_unknown_stereo;
    TEST_ASSERT(
        m->getBondWithIdx(2)->getPropIfPresent<int>(
            common_properties::_UnknownStereo, explicit_unknown_stereo) &&
        explicit_unknown_stereo)

    TEST_ASSERT(m->getBondWithIdx(0)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo() ==
                Bond::STEREOANY);  // now set in detectStereochemistry
    ROMol(*m).getBondWithIdx(0)->setStereo(
        Bond::STEREOCIS);  // ensure we have stereo atoms
    MolOps::sanitizeMol(*m);
    TEST_ASSERT(m->getBondWithIdx(0)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo() == Bond::STEREOANY);
    TEST_ASSERT(m->getBondWithIdx(1)->getBondType() == Bond::SINGLE);
    TEST_ASSERT(m->getBondWithIdx(1)->getBondDir() == Bond::NONE);
    TEST_ASSERT(m->getBondWithIdx(2)->getBondType() == Bond::SINGLE);
    TEST_ASSERT(m->getBondWithIdx(2)->getBondDir() != Bond::UNKNOWN);
    MolOps::assignStereochemistry(*m, true, true);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo() == Bond::STEREOANY);
    delete m;
  }
  {  // chiral center
    std::string fName;
    fName = rdbase + "github1034.2.mol";
    bool sanitize = true;
    RWMol *m = MolFileToMol(fName, sanitize);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 5);
    TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    TEST_ASSERT(m->getBondWithIdx(3)->getBondType() == Bond::SINGLE);
    TEST_ASSERT(m->getBondWithIdx(3)->getBondDir() == Bond::NONE);
    TEST_ASSERT(
        m->getBondWithIdx(3)->hasProp(common_properties::_UnknownStereo));
    delete m;
  }
  {  // chiral center
    std::string fName;
    fName = rdbase + "github1034.2.mol";
    bool sanitize = false;
    RWMol *m = MolFileToMol(fName, sanitize);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 5);
    TEST_ASSERT(m->getBondWithIdx(3)->getBondType() == Bond::SINGLE);
    TEST_ASSERT(m->getBondWithIdx(3)->getBondDir() != Bond::UNKNOWN);
    int explicit_unknown_stereo;
    TEST_ASSERT(
        m->getBondWithIdx(3)->getPropIfPresent<int>(
            common_properties::_UnknownStereo, explicit_unknown_stereo) &&
        explicit_unknown_stereo)
    MolOps::sanitizeMol(*m);
    TEST_ASSERT(m->getBondWithIdx(3)->getBondType() == Bond::SINGLE);
    TEST_ASSERT(m->getBondWithIdx(3)->getBondDir() != Bond::UNKNOWN);
    TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub1049() {
  BOOST_LOG(rdInfoLog) << "Test github 1049: MolOps::cleanUp() being called by "
                          "CTAB parser even when sanitization isn't on"
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  {  // no stereo (this one worked before)
    std::string fName;
    fName = rdbase + "github1049.1.mol";
    bool sanitize = false;
    RWMol *m = MolFileToMol(fName, sanitize);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 7);
    TEST_ASSERT(m->getAtomWithIdx(1)->getAtomicNum() == 7);
    TEST_ASSERT(m->getAtomWithIdx(1)->getFormalCharge() == 0);
    TEST_ASSERT(m->getAtomWithIdx(2)->getAtomicNum() == 8);
    TEST_ASSERT(m->getAtomWithIdx(2)->getFormalCharge() == 0);
    TEST_ASSERT(m->getAtomWithIdx(3)->getAtomicNum() == 8);
    TEST_ASSERT(m->getAtomWithIdx(3)->getFormalCharge() == 0);

    MolOps::sanitizeMol(*m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getFormalCharge() == 1);
    TEST_ASSERT((m->getAtomWithIdx(2)->getFormalCharge() == -1 &&
                 m->getAtomWithIdx(3)->getFormalCharge() == 0) ||
                (m->getAtomWithIdx(2)->getFormalCharge() == 0 &&
                 m->getAtomWithIdx(3)->getFormalCharge() == -1));

    delete m;
  }
  {  // with stereo (this one did not work)
    std::string fName;
    fName = rdbase + "github1049.2.mol";
    bool sanitize = false;
    RWMol *m = MolFileToMol(fName, sanitize);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 7);
    TEST_ASSERT(m->getAtomWithIdx(1)->getAtomicNum() == 7);
    TEST_ASSERT(m->getAtomWithIdx(1)->getFormalCharge() == 0);
    TEST_ASSERT(m->getAtomWithIdx(2)->getAtomicNum() == 8);
    TEST_ASSERT(m->getAtomWithIdx(2)->getFormalCharge() == 0);
    TEST_ASSERT(m->getAtomWithIdx(3)->getAtomicNum() == 8);
    TEST_ASSERT(m->getAtomWithIdx(3)->getFormalCharge() == 0);

    MolOps::sanitizeMol(*m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getFormalCharge() == 1);
    TEST_ASSERT((m->getAtomWithIdx(2)->getFormalCharge() == -1 &&
                 m->getAtomWithIdx(3)->getFormalCharge() == 0) ||
                (m->getAtomWithIdx(2)->getFormalCharge() == 0 &&
                 m->getAtomWithIdx(3)->getFormalCharge() == -1));

    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMolFileDativeBonds() {
  BOOST_LOG(rdInfoLog) << "Test MDL molfiles with dative bonds (V3000 only)"
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  // Read molfiles with dative bonds.
  {
    std::string fName = rdbase + "dative_bonds_one.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumBonds() == 5);
    TEST_ASSERT(m->getBondWithIdx(4)->getBondType() == Bond::DATIVE);

    std::string smiles = MolToSmiles(*m);
    TEST_ASSERT(smiles == "CCC(=O)O->[Cu]");

    delete m;
  }

  {
    std::string fName = rdbase + "dative_bonds_two.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumBonds() == 10);
    TEST_ASSERT(m->getBondWithIdx(8)->getBondType() == Bond::DATIVE);
    TEST_ASSERT(m->getBondWithIdx(9)->getBondType() == Bond::DATIVE);

    std::string smiles = MolToSmiles(*m);
    TEST_ASSERT(smiles == "CCC(=O)O->[Cu]<-OC(O)CC");

    delete m;
  }

  {
    std::string fName = rdbase + "DativeBond2000.mol";
    RWMol *m = MolFileToMol(fName, false);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumBonds() == 5);
    TEST_ASSERT(m->getBondWithIdx(4)->getBondType() == Bond::DATIVE);

    std::string smiles = MolToSmiles(*m);
    TEST_ASSERT(smiles == "CC(C)->[Mg](Cl)Cl");

    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub1251() {
  BOOST_LOG(rdInfoLog)
      << "Test github 1251: MolFromMolBlock sanitizing when it should not be"
      << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";

  {
    std::string fName = rdbase + "github1251.mol";
    RWMol *m = MolFileToMol(fName, false);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 39);
    TEST_ASSERT(m->getNumBonds() == 44);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub1029() {
  BOOST_LOG(rdInfoLog)
      << "Test github 1029: PDB reader fails for arginine explicit hydrogens "
      << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  {  // the original bug report
    std::string fName = rdbase + "github1029.1.pdb";
    bool sanitize = true, removeHs = false;
    ROMol *m = PDBFileToMol(fName, sanitize, removeHs);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 3268);
    TEST_ASSERT(m->getNumBonds() == 3302);
    TEST_ASSERT(m->getAtomWithIdx(121)->getExplicitValence() == 4);
    TEST_ASSERT(m->getAtomWithIdx(121)->getFormalCharge() == 1);

    delete m;
  }
  {  // a second report that came in
    std::string fName = rdbase + "github1029.1jld_chaina.pdb";
    bool sanitize = false, removeHs = false;
    ROMol *m = PDBFileToMol(fName, sanitize, removeHs);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 1533);
    TEST_ASSERT(m->getNumBonds() == 1545);
    TEST_ASSERT(m->getAtomWithIdx(123)->getExplicitValence() == 4);
    TEST_ASSERT(m->getAtomWithIdx(123)->getFormalCharge() == 1);

    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub1340() {
  BOOST_LOG(rdInfoLog) << "Test github 1340: PDB parser creating H-H bonds "
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  {  // a second report that came in
    std::string fName = rdbase + "github1340.1jld_snip.pdb";
    bool sanitize = true, removeHs = false;
    ROMol *m = PDBFileToMol(fName, sanitize, removeHs);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 37);
    TEST_ASSERT(m->getNumBonds() == 35);
    TEST_ASSERT(m->getAtomWithIdx(10)->getAtomicNum() == 1);
    TEST_ASSERT(m->getAtomWithIdx(34)->getAtomicNum() == 1);
    RDGeom::Point3D p10 = m->getConformer().getAtomPos(10);
    RDGeom::Point3D p34 = m->getConformer().getAtomPos(34);
    TEST_ASSERT((p34 - p10).length() < 1.0);
    TEST_ASSERT(!m->getBondBetweenAtoms(10, 34));

    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMolBlockChirality() {
  BOOST_LOG(rdInfoLog)
      << "Test automatic generation of coordinates for mol block chirality "
      << std::endl;
  {
    std::string smi = "C[C@H](Cl)Br";
    auto mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    auto mb = MolToMolBlock(*mol, true);
    auto mol2 = MolBlockToMol(mb);
    TEST_ASSERT(mol2);
    auto csmi1 = MolToSmiles(*mol, true);
    auto csmi2 = MolToSmiles(*mol2, true);
    TEST_ASSERT(csmi1.find("@") != std::string::npos);
    TEST_ASSERT(csmi1 == csmi2);
    delete mol;
    delete mol2;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMolBlock3DStereochem() {
  BOOST_LOG(rdInfoLog)
      << "Test automatic perception of stereochem from 3D structure "
      << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  {
    std::string fName = rdbase + "stereo3d_1.mol";
    auto m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    std::string smi = "F[C@](Cl)(Br)I";
    auto m2 = SmilesToMol(smi);
    TEST_ASSERT(m2);
    auto csmi1 = MolToSmiles(*m1, true);
    auto csmi2 = MolToSmiles(*m2, true);
    TEST_ASSERT(csmi1.find("@") != std::string::npos);
    TEST_ASSERT(csmi1 == csmi2);
    delete m1;
    delete m2;
  }
  {
    std::string fName = rdbase + "stereo3d_2.mol";
    auto m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    std::string smi = "F[C@@](Cl)(Br)I";
    auto m2 = SmilesToMol(smi);
    TEST_ASSERT(m2);
    auto csmi1 = MolToSmiles(*m1, true);
    auto csmi2 = MolToSmiles(*m2, true);
    TEST_ASSERT(csmi1.find("@") != std::string::npos);
    TEST_ASSERT(csmi1 == csmi2);
    delete m1;
    delete m2;
  }
  {
    std::string fName = rdbase + "stereo3d_conflict.mol";
    auto m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    std::string smi = "F[C@](Cl)(Br)I";
    auto m2 = SmilesToMol(smi);
    TEST_ASSERT(m2);
    auto csmi1 = MolToSmiles(*m1, true);
    auto csmi2 = MolToSmiles(*m2, true);
    TEST_ASSERT(csmi1.find("@") != std::string::npos);
    TEST_ASSERT(csmi1 == csmi2);
    delete m1;
    delete m2;
  }
  {
    std::string fName = rdbase + "stereo3d_unknown.mol";
    auto m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    std::string smi = "F[C](Cl)(Br)I";
    auto m2 = SmilesToMol(smi);
    TEST_ASSERT(m2);
    auto csmi1 = MolToSmiles(*m1, true);
    auto csmi2 = MolToSmiles(*m2, true);
    TEST_ASSERT(csmi1.find("@") == std::string::npos);
    TEST_ASSERT(csmi1 == csmi2);
    delete m1;
    delete m2;
  }
  {
    std::string fName = rdbase + "stereo3d_trans.mol";
    auto m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    std::string smi = "F/C=C/F";
    auto m2 = SmilesToMol(smi);
    TEST_ASSERT(m2);
    auto csmi1 = MolToSmiles(*m1, true);
    auto csmi2 = MolToSmiles(*m2, true);
    TEST_ASSERT(csmi1.find("/") != std::string::npos);
    TEST_ASSERT(csmi1 == csmi2);
    delete m1;
    delete m2;
  }
  {
    std::string fName = rdbase + "stereo3d_cis.mol";
    auto m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    std::string smi = "F/C=C\\F";
    auto m2 = SmilesToMol(smi);
    TEST_ASSERT(m2);
    auto csmi1 = MolToSmiles(*m1, true);
    auto csmi2 = MolToSmiles(*m2, true);
    TEST_ASSERT(csmi1.find("/") != std::string::npos);
    TEST_ASSERT(csmi1.find("\\") != std::string::npos);
    TEST_ASSERT(csmi1 == csmi2);
    delete m1;
    delete m2;
  }

  {
    std::string fName = rdbase + "stereo3d_dblunknown.mol";
    auto m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    std::string smi = "FC=CF";
    auto m2 = SmilesToMol(smi);
    TEST_ASSERT(m2);
    auto csmi1 = MolToSmiles(*m1, true);
    auto csmi2 = MolToSmiles(*m2, true);
    TEST_ASSERT(csmi1.find("/") == std::string::npos);
    TEST_ASSERT(csmi1.find("\\") == std::string::npos);
    TEST_ASSERT(csmi1 == csmi2);
    delete m1;
    delete m2;
  }
}

void testMarvinSMATag() {
  BOOST_LOG(rdInfoLog) << "Test Marvin MRV SMA tag " << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  {
    std::string fName = rdbase + "mrv-sma.mol";
    RWMol *m = MolFileToMol(fName, false);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 4);
    TEST_ASSERT(m->getAtomWithIdx(0)->getProp<std::string>(
                    common_properties::MRV_SMA) == "[#6;r6]");
    TEST_ASSERT(m->getAtomWithIdx(1)->getProp<std::string>(
                    common_properties::MRV_SMA) == "[#16;H1]");
    TEST_ASSERT(m->getAtomWithIdx(2)->getProp<std::string>(
                    common_properties::MRV_SMA) == "[#6;r6]");
    TEST_ASSERT(m->getAtomWithIdx(3)->getProp<std::string>(
                    common_properties::MRV_SMA) == "[#7;H2A]");
    // this should be similar to [#7;AH2:4][c;r6:3]:[c;r6:1]-[#16H1:2]
    //  RDKit makes these recursive smarts, not "ANDED" smarts which are a
    //  simpler case
    std::string sma = MolToSmarts(*m);
    TEST_ASSERT(sma ==
                "[#6&$([#6&r6]):1](-[#16&$([#16&H1]):2]):[#6&$([#6&r6]):3]-[#7&"
                "$([#7&H2&A]):4]");
    delete m;
  }

  {
    std::string fName = rdbase + "mrv-sma-bad.mol";
    bool ok = false;
    try {
      MolFileToMol(fName, false);
    } catch (FileParseException &e) {
      ok = true;
      TEST_ASSERT(
          std::string("Cannot parse smarts: 'MyDogHasFleas' on line 12") ==
          e.what());
    }
    TEST_ASSERT(ok);
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub1689() {
  BOOST_LOG(rdInfoLog) << "Test github 1689: Play nice with naughty MOL blocks"
                       << std::endl;

  std::string molb =
      "rdkit_blank_line_before_M_END_test.sdf\n"
      "  ChemDraw12181709392D\n"
      "\n"
      "  2  1  0  0  0  0  0  0  0  0999 V2000\n"
      "   -0.3572   -0.2062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "    0.3572    0.2062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "  1  2  1  0      \n"
      "\n"
      "M  END\n";
  {
    bool sanitize = true, removeHs = true, strictParsing = false;
    ROMol *m = MolBlockToMol(molb, sanitize, removeHs, strictParsing);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 2);
    TEST_ASSERT(m->getNumBonds() == 1);
    delete m;
  }
  {
    bool sanitize = true, removeHs = true, strictParsing = true;
    bool ok = false;
    try {
      MolBlockToMol(molb, sanitize, removeHs, strictParsing);
    } catch (FileParseException &) {
      ok = true;
    }
    TEST_ASSERT(ok);
  }
}

void testWedgeBondToDoublebond() {
  BOOST_LOG(rdInfoLog) << "Test wedged bonds to double bonds " << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  {  // a second report that came in
    std::string fName = rdbase + "wedged_single_to_double_bond.mol";
    bool sanitize = false, removeHs = false;
    ROMol *m = MolFileToMol(fName, sanitize, removeHs);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 49);
    TEST_ASSERT(m->getNumBonds() == 51);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub1615() {
  BOOST_LOG(rdInfoLog) << "testing github #1615: add function WedgeBond()"
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  {
    RWMol *m1;
    std::string fName = rdbase + "Issue399a.mol";
    m1 = MolFileToMol(fName);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    TEST_ASSERT(m1->getAtomWithIdx(1)->getProp<std::string>(
                    common_properties::_CIPCode) == "S");
    TEST_ASSERT(m1->getBondWithIdx(0)->getBondDir() == Bond::NONE);
    Chirality::wedgeBond(m1->getBondWithIdx(0), 1, &m1->getConformer());
    TEST_ASSERT(m1->getBondWithIdx(0)->getBondDir() == Bond::BEGINWEDGE);

    delete m1;
  }

  {
    std::string mb =
        "example\n"
        "  ChemDraw04050615582D\n"
        "\n"
        "  4  4  0  0  0  0  0  0  0  0999 V2000\n"
        "   -0.7697    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "    0.0553    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "    0.7697    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "    0.7697   -0.4125    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  "
        "0\n"
        "  2  1  1  0\n"
        "  2  3  1  0\n"
        "  3  4  1  0\n"
        "  2  4  1  0\n"
        "M  END\n";
    RWMol *m1 = MolBlockToMol(mb);
    TEST_ASSERT(m1);
    TEST_ASSERT(!m1->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    TEST_ASSERT(m1->getBondWithIdx(0)->getBondDir() == Bond::NONE);
    TEST_ASSERT(m1->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    m1->getAtomWithIdx(1)->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
    MolOps::assignStereochemistry(*m1, true, true);
    TEST_ASSERT(m1->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    TEST_ASSERT(m1->getAtomWithIdx(1)->getProp<std::string>(
                    common_properties::_CIPCode) == "S");
    TEST_ASSERT(m1->getBondWithIdx(0)->getBondDir() == Bond::NONE);

    Chirality::wedgeBond(m1->getBondWithIdx(0), 1, &m1->getConformer());
    TEST_ASSERT(m1->getBondWithIdx(0)->getBondDir() == Bond::BEGINWEDGE);

    delete m1;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub2000() {
  BOOST_LOG(rdInfoLog)
      << "testing github #2000: Error while parsing empty atom list"
      << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  {
    std::string fName = rdbase + "github2000.sdf";
    std::unique_ptr<RWMol> m1(MolFileToMol(fName));
    TEST_ASSERT(m1);
    TEST_ASSERT(!m1->getAtomWithIdx(0)->hasQuery());
  }
  {
    std::string fName = rdbase + "github2000.2.sdf";
    bool ok = false;
    try {
      std::unique_ptr<RWMol> m1(MolFileToMol(fName));
    } catch (FileParseException &) {
      ok = true;
    }
    TEST_ASSERT(ok);
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testAtomQueries() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing atom queries contains HCOUNT, RBCOUNT, UNSAT"
                       << std::endl;
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/FileParsers/test_data/";
    std::unique_ptr<RWMol> testQuery(
        MolFileToMol(pathName + "AtomQuery1.mol", true, false, true));
    TEST_ASSERT(testQuery);

    std::unique_ptr<RWMol> testMol(
        SmilesToMol("[H]C([H])([H])C1=CN=C(C=C)C(C2CCCCC2)=C1C"));
    TEST_ASSERT(testMol);

    MatchVectType mv;
    SubstructMatch(*testMol, *testQuery, mv);
    TEST_ASSERT(mv.size() > 0);

    testMol.reset(SmilesToMol("[H]CC1=CN=C(C=C)C(C2CCCCC2)=C1C"));
    TEST_ASSERT(testMol);

    SubstructMatch(*testMol, *testQuery, mv);
    TEST_ASSERT(mv.size() > 0);

    testMol.reset(SmilesToMol("CCC1=CN=C(C=C)C(C2CCCCC2)=C1C"));
    TEST_ASSERT(testMol);

    SubstructMatch(*testMol, *testQuery, mv);
    TEST_ASSERT(mv.size() == 0);  // search fails
  }
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/FileParsers/test_data/";
    std::unique_ptr<RWMol> testQuery(
        MolFileToMol(pathName + "AtomQuery2.mol", true, false, true));
    TEST_ASSERT(testQuery);

    std::unique_ptr<RWMol> testMol(SmilesToMol("C1(C)C(C)=CN=C(C)C=1C1CCCCC1"));
    TEST_ASSERT(testMol);

    MatchVectType mv;
    SubstructMatch(*testMol, *testQuery, mv);
    TEST_ASSERT(mv.size() > 0);

    testMol.reset(SmilesToMol("C1(C)C(C)=CN=C(C)C=1C12CCC(CC1)CC2"));
    TEST_ASSERT(testMol);

    SubstructMatch(*testMol, *testQuery, mv);
    TEST_ASSERT(mv.size() == 0);  // search fails
  }
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/FileParsers/test_data/";
    std::unique_ptr<RWMol> testQuery(
        MolFileToMol(pathName + "AtomQuery3.mol", true, false, true));
    TEST_ASSERT(testQuery);

    std::unique_ptr<RWMol> testMol(
        SmilesToMol("CCC1=CN=C(C=C)C(C2CCCCC2)=C1C"));
    TEST_ASSERT(testMol);

    MatchVectType mv;
    SubstructMatch(*testMol, *testQuery, mv);
    TEST_ASSERT(mv.size() > 0);

    testMol.reset(SmilesToMol("CCC1=CN=C(C3=CC=CC=C3)C(C2CCCCC2)=C1C"));
    TEST_ASSERT(testMol);

    SubstructMatch(*testMol, *testQuery, mv);
    TEST_ASSERT(mv.size() > 0);

    testMol.reset(SmilesToMol("CCC1=CN=C(C3CCCCC3)C(C2CCCCC2)=C1C"));
    TEST_ASSERT(testMol);

    SubstructMatch(*testMol, *testQuery, mv);
    TEST_ASSERT(mv.size() == 0);
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void RunTests() {
#if 1
  test1();
  test2();
  test4();
  test5();
  test6();
  testIssue145();
  testIssue148();
  test7();
  test8();
  testIssue180();
  testIssue264();
  testIssue399();
  testMolFileChgLines();
  testDblBondStereochem();
  testSymmetricDblBondStereochem();
  testRingDblBondStereochem();
  testMolFileRGroups();
  testMolFileDegreeQueries();
  testMolFileRBCQueries();
  testMolFileUnsaturationQueries();
  testMolFileQueryToSmarts();
  testMissingFiles();
  testIssue1965035();
  testRadicals();
  testBadBondOrders();
  testAtomParity();
  testIssue2692246();
  testKekulizationSkip();
  testMolFileAtomValues();
  testMolFileAtomQueries();
  testListsAndValues();
  // testCrash();
  test1V3K();
  testIssue2963522();
  testIssue3073163();
  testIssue3154208();
  testIssue3228150();
  testIssue3313540();
  testIssue3359739();
  testIssue3374639();
  testThreeCoordinateChirality();
  testIssue3375647();
  testIssue3375684();
  testChiralPhosphorous();
  testIssue3392107();
  testIssue3477283();
  testIssue3484552();
  testIssue3514824();
  testIssue3525799();
  testSkipLines();
  testIssue269();
  testMolFileChiralFlag();
  testMolFileTotalValence();
  testGithub88();
  testGithub82();
  testMolFileWithHs();
  testMolFileWithRxn();
  testGithub166();
  testZBO();

  testGithub164();
  testGithub194();
  testGithub196();
  testIssue3557675();
  test3V3K();
  test2V3K();
  testGithub191();
  testGithub210();
  testPDBResidues();
  testGithub337();
  testGithub360();
  testGithub741();
  testGithub188();
  testRCSBSdf();
  testParseCHG();
  testMDLAtomProps();
  testSupplementalSmilesLabel();
  testGithub1023();
  testGithub1049();
  testPDBFile();
  testSequences();

  // testSequenceReaders();

  testMolFileDativeBonds();
  testGithub1251();
  testMarvinSMATag();
  testGithub1029();
  testGithub1340();
  testGithub1034();
  testMolBlockChirality();
  testMolBlock3DStereochem();
  testGithub1689();
  testWedgeBondToDoublebond();
  testGithub1615();
  testGithub2000();
  testAtomQueries();
#endif
}

// must be in German Locale for test...
void testLocaleSwitcher() {
  float d = -1.0;
  char buffer[1024];
  sprintf(buffer, "%0.2f", d);
  if (std::string(buffer) != "-1,00") {
    BOOST_LOG(rdInfoLog) << " ---- no German locale support (skipping) ---- "
                         << std::endl;
    return;
  }

  {
    RDKit::Utils::LocaleSwitcher ls;
    sprintf(buffer, "%0.2f", d);
    CHECK_INVARIANT(std::string(buffer) == "-1.00", "Locale Switcher Fail");
    // test locale switcher recursion
    {
      RDKit::Utils::LocaleSwitcher ls;
      sprintf(buffer, "%0.2f", d);
      CHECK_INVARIANT(std::string(buffer) == "-1.00", "Locale Switcher Fail");
    }
    // should still be in the "C" variant
    sprintf(buffer, "%0.2f", d);
    CHECK_INVARIANT(std::string(buffer) == "-1.00", "Locale Switcher Fail");
  }

  // Should be back in German Locale
  sprintf(buffer, "%0.2f", d);
  CHECK_INVARIANT(std::string(buffer) == "-1,00", "Locale Switcher Fail");
}

#ifdef RDK_TEST_MULTITHREADED

#include <thread>
#include <future>

namespace {
void runblock() {
  std::setlocale(LC_ALL, "de_DE.UTF-8");
  testLocaleSwitcher();
}
}  // namespace
void testMultiThreadedSwitcher() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test multithreading Locale Switching"
                        << std::endl;

  std::vector<std::future<void>> tg;
  unsigned int count = 100;
  for (unsigned int i = 0; i < count; ++i) {
    tg.emplace_back(std::async(std::launch::async, runblock));
  }
  for (auto &fut : tg) {
    fut.get();
  }
  BOOST_LOG(rdErrorLog) << "    Test multithreading (Done)" << std::endl;
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
}

#else

void testMultiThreadedSwitcher() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << " ---- Multithreaded tests disabled ---- "
                       << std::endl;
}
#endif

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;
  //  std::locale::global(std::locale("de_DE.UTF-8"));

  RDLog::InitLogs();
  BOOST_LOG(rdInfoLog) << " ---- Running with POSIX locale ----- " << std::endl;
  RunTests();  // run with C locale

  BOOST_LOG(rdInfoLog) << " ---- Running with German locale ----- "
                       << std::endl;
  std::setlocale(LC_ALL, "de_DE.UTF-8");
  std::cout << std::setlocale(LC_ALL, nullptr) << std::endl;
  testLocaleSwitcher();
  testMultiThreadedSwitcher();
  RunTests();

  return 0;
}
