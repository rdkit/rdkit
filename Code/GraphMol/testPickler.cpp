//
//  Copyright (C) 2004-2018 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include <RDGeneral/utils.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <RDGeneral/RDLog.h>

#include <stdlib.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>

using namespace RDKit;

void test1(bool doLong = 0) {
  std::string fName = getenv("RDBASE");
  fName += "/Code/GraphMol/test_data/canonSmiles.smi";

  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing String Pickle Roundtrips." << std::endl;

  SmilesMolSupplier suppl(fName, "\t", 0, 1, false);

  int count = 0;
  while (!suppl.atEnd()) {
    ROMol *m = suppl.next();
    TEST_ASSERT(m);

    std::string pickle;
    MolPickler::pickleMol(*m, pickle);
    TEST_ASSERT(pickle.size());

    ROMol m2;
    MolPickler::molFromPickle(pickle, m2);
    std::string smi1 = MolToSmiles(*m);
    // BOOST_LOG(rdInfoLog) << smi << " " << smi1 << std::endl;
    std::string smi2 = MolToSmiles(m2);
    if (smi1 != smi2) {
      BOOST_LOG(rdInfoLog) << "Line: " << count << "\n  " << smi1
                           << "\n != \n  " << smi2 << std::endl;
    }
    TEST_ASSERT(smi1 == smi2);
    delete m;
    count++;
    if (!doLong && count >= 100) break;
  }
  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void _createPickleFile() {
  std::string smiName = getenv("RDBASE");
  smiName += "/Code/GraphMol/test_data/canonSmiles.smi";
  std::string pklName = getenv("RDBASE");
#ifdef OLD_PICKLE
  pklName += "/Code/GraphMol/test_data/canonSmiles.v1.pkl";
#else
  pklName += "/Code/GraphMol/test_data/canonSmiles.v2.pkl";
#endif
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "creating pickle file." << std::endl;

  SmilesMolSupplier suppl(smiName, "\t", 0, 1, false);
  std::ofstream outStream(pklName.c_str(), std::ios_base::binary);
  int count = 0;
  while (!suppl.atEnd()) {
    ROMol *m = suppl.next();
    TEST_ASSERT(m);

    std::string pickle;
    MolPickler::pickleMol(*m, outStream);
    delete m;
    count++;
  }
  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void test2(bool doLong = 0) {
  std::string smiName = getenv("RDBASE");
  smiName += "/Code/GraphMol/test_data/canonSmiles.smi";
  std::string pklName = getenv("RDBASE");
  pklName += "/Code/GraphMol/test_data/canonSmiles.v1.pkl";

  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing reading existing pickle file (v1)."
                        << std::endl;

  SmilesMolSupplier suppl(smiName, "\t", 0, 1, false);
  std::ifstream inStream(pklName.c_str(), std::ios_base::binary);
  int count = 0;
  while (!suppl.atEnd()) {
    ROMol *m1 = suppl.next();
    TEST_ASSERT(m1);
    ROMol m2;
    MolPickler::molFromPickle(inStream, m2);

    std::string smi1 = MolToSmiles(*m1);
    std::string smi2 = MolToSmiles(m2);

    if (smi1 != smi2) {
      BOOST_LOG(rdInfoLog) << "Line: " << count << "\n  " << smi1
                           << "\n != \n  " << smi2 << std::endl;
    }
    TEST_ASSERT(smi1 == smi2);
    delete m1;
    count++;
    if (!doLong && count >= 100) break;
  }
  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void test3(bool doLong = 0) {
  std::string smiName = getenv("RDBASE");
  smiName += "/Code/GraphMol/test_data/canonSmiles.smi";
  std::string pklName = getenv("RDBASE");
  pklName += "/Code/GraphMol/test_data/canonSmiles.v2.pkl";

  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing reading existing pickle file (v2)."
                        << std::endl;

  SmilesMolSupplier suppl(smiName, "\t", 0, 1, false);
  std::ifstream inStream(pklName.c_str(), std::ios_base::binary);
  int count = 0;
  while (!suppl.atEnd()) {
    ROMol *m1 = suppl.next();
    TEST_ASSERT(m1);
    ROMol m2;
    MolPickler::molFromPickle(inStream, m2);

    std::string smi1 = MolToSmiles(*m1, 1);
    std::string smi2 = MolToSmiles(m2, 1);

    if (smi1 != smi2) {
      BOOST_LOG(rdInfoLog) << "Line: " << count << "\n  " << smi1
                           << "\n != \n  " << smi2 << std::endl;
    }
    TEST_ASSERT(smi1 == smi2);
    delete m1;
    count++;
    if (!doLong && count >= 100) break;
  }
  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void timeTest(bool doLong = 0) {
  time_t t1, t2;
  std::string smiName = getenv("RDBASE");
  smiName += "/Code/GraphMol/test_data/canonSmiles.smi";

  std::string pklName = getenv("RDBASE");
  pklName += "/Code/GraphMol/test_data/canonSmiles.v2.pkl";

  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Timing reads." << std::endl;

  t1 = std::time(nullptr);
  SmilesMolSupplier suppl(smiName, "\t", 0, 1, false);
  int count = 0;
  while (!suppl.atEnd()) {
    ROMol *m1 = suppl.next();
    TEST_ASSERT(m1);
    count++;
    if (!doLong && count >= 100) break;
    delete m1;
  }
  t2 = std::time(nullptr);
  BOOST_LOG(rdInfoLog) << " Smiles time: " << std::difftime(t2, t1)
                       << std::endl;
  ;

  std::ifstream inStream(pklName.c_str(), std::ios_base::binary);
  t1 = std::time(nullptr);
  while (count > 0) {
    ROMol m2;
    MolPickler::molFromPickle(inStream, m2);
    count--;
    if (!doLong && count >= 100) break;
  }
  t2 = std::time(nullptr);
  BOOST_LOG(rdInfoLog) << " Pickle time: " << std::difftime(t2, t1)
                       << std::endl;
  ;

  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void test4() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing combining molecules using pickles."
                        << std::endl;

  ROMol *m1 = SmilesToMol("C1CC1");
  TEST_ASSERT(m1);
  ROMol *m2 = SmilesToMol("C1COC1");
  TEST_ASSERT(m2);
  ROMol m3;
  std::string pickle;
  MolPickler::pickleMol(*m1, pickle);
  MolPickler::molFromPickle(pickle, m3);
  MolPickler::pickleMol(*m2, pickle);
  MolPickler::molFromPickle(pickle, m3);

  m3.debugMol(std::cout);
  TEST_ASSERT(m3.getNumAtoms() == 7);
  TEST_ASSERT(m3.getNumBonds() == 7);
  TEST_ASSERT(m3.getRingInfo()->numRings() == 2);
  TEST_ASSERT(m3.getRingInfo()->isAtomInRingOfSize(0, 3));
  TEST_ASSERT(!m3.getRingInfo()->isAtomInRingOfSize(0, 4));
  TEST_ASSERT(!m3.getRingInfo()->isAtomInRingOfSize(4, 3));
  TEST_ASSERT(m3.getRingInfo()->isAtomInRingOfSize(4, 4));

  TEST_ASSERT(m3.getRingInfo()->isBondInRingOfSize(0, 3));
  TEST_ASSERT(!m3.getRingInfo()->isBondInRingOfSize(0, 4));
  TEST_ASSERT(!m3.getRingInfo()->isBondInRingOfSize(4, 3));
  TEST_ASSERT(m3.getRingInfo()->isBondInRingOfSize(4, 4));
  delete m1;
  delete m2;

  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void testIssue164() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing Issue164: Roundtrip Pickle Failure."
                        << std::endl;

  std::string smi =
      "NCCCCC(NC(C(C(C)C)NC(=O)C1CCCN1C(C(N)CCCNC(=N)N)=O)=O)C(NC(C(C)C)C(NC("
      "Cc2ccc(O)cc2)C(=O)N6C(C(NC(CC(N)=O)C(NCC(NC(C)C(NC(CCC(O)=O)C(NC(CC(O)="
      "O)C(NC(C(NC(CO)C(NC(C)C(NC(CCC(O)=O)C(NC(C)C(NC(Cc3ccccc3)C(=O)N5C(C(NC("
      "CC(C)C)C(NC(C(NC(C(O)=O)Cc4ccccc4)=O)CCC(O)=O)=O)=O)CCC5)=O)=O)=O)=O)=O)"
      "CCC(O)=O)=O)=O)=O)=O)=O)=O)CCC6)=O)=O";
  ROMol *m1 = SmilesToMol(smi);
  TEST_ASSERT(m1);
  std::string pickle;
  ROMol m2;
  MolPickler::pickleMol(*m1, pickle);
  MolPickler::molFromPickle(pickle, m2);

  TEST_ASSERT(m1->getNumAtoms() == m2.getNumAtoms());

  // the issue had to do with number of atoms, so let's make an enormous
  // molecule and try again:
  RWMol *m3 = static_cast<RWMol *>(SmilesToMol(smi));
  m3->insertMol(*m1);
  m3->insertMol(*m1);
  MolOps::sanitizeMol(*m3);

  MolPickler::pickleMol(*m3, pickle);
  ROMol m4;
  MolPickler::molFromPickle(pickle, m4);

  TEST_ASSERT(m3->getNumAtoms() == m4.getNumAtoms());
  TEST_ASSERT(m3->getNumBonds() == m4.getNumBonds());
  delete m1;
  delete m3;

  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void testIssue219() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "Testing Issue219: Pickle Failure with Conformations." << std::endl;

  std::string smi = "CC";
  ROMol *m1 = SmilesToMol(smi);
  TEST_ASSERT(m1);
  std::string pickle;
  ROMol *m2;
  MolPickler::pickleMol(*m1, pickle);
  m2 = new ROMol();
  MolPickler::molFromPickle(pickle, *m2);
  TEST_ASSERT(m1->getNumAtoms() == m2->getNumAtoms());

  auto *conf = new Conformer(2);
  conf->setId(23);
  m1->addConformer(conf);
  MolPickler::pickleMol(*m1, pickle);
  delete m2;
  m2 = new ROMol();
  MolPickler::molFromPickle(pickle, *m2);
  TEST_ASSERT(m1->getNumAtoms() == m2->getNumAtoms());
  TEST_ASSERT(m2->getConformer().getId() == 23);
  delete m1;
  delete m2;

  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void testIssue220() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "Testing Issue220: Pickle problems with bond stereochemistry."
      << std::endl;

  std::string smi = "N/N=C/C";
  ROMol *m1 = SmilesToMol(smi);
  TEST_ASSERT(m1);
  TEST_ASSERT(m1->getBondWithIdx(1)->getStereo() != Bond::STEREONONE);
  std::string pickle;
  ROMol *m2;
  MolPickler::pickleMol(*m1, pickle);
  m2 = new ROMol();
  MolPickler::molFromPickle(pickle, *m2);
  TEST_ASSERT(m1->getNumAtoms() == m2->getNumAtoms());
  TEST_ASSERT(m1->getNumBonds() == m2->getNumBonds());
  TEST_ASSERT(m2->getBondWithIdx(1)->getStereo() != Bond::STEREONONE);
  for (unsigned int i = 0; i < m1->getNumBonds(); ++i) {
    TEST_ASSERT(m1->getBondWithIdx(i)->getStereo() ==
                m2->getBondWithIdx(i)->getStereo());
    TEST_ASSERT(m1->getBondWithIdx(i)->getStereoAtoms() ==
                m2->getBondWithIdx(i)->getStereoAtoms());
  }

  delete m1;
  delete m2;
  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void testQueries() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing query serialization." << std::endl;
  std::string smi, pickle, smi2;
  int tmpInt;
  ROMol *m1, *m2;
  MatchVectType matchV;

  // start simple : atom map numbers
  smi = "C";
  m1 = SmilesToMol(smi);
  TEST_ASSERT(m1);
  m1->getAtomWithIdx(0)->setProp(common_properties::molAtomMapNumber, 1);
  MolPickler::pickleMol(*m1, pickle);
  delete m1;
  m1 = new ROMol();
  MolPickler::molFromPickle(pickle, *m1);
  TEST_ASSERT(m1->getNumAtoms() == 1);
  TEST_ASSERT(
      m1->getAtomWithIdx(0)->hasProp(common_properties::molAtomMapNumber));
  m1->getAtomWithIdx(0)->getProp(common_properties::molAtomMapNumber, tmpInt);
  TEST_ASSERT(tmpInt == 1);
  delete m1;

  // now a basic query:
  smi = "C";
  m1 = SmartsToMol(smi);
  TEST_ASSERT(m1);
  MolPickler::pickleMol(*m1, pickle);
  delete m1;
  m1 = new ROMol();
  MolPickler::molFromPickle(pickle, *m1);
  TEST_ASSERT(m1->getNumAtoms() == 1);
  TEST_ASSERT(m1->getAtomWithIdx(0)->hasQuery());
  TEST_ASSERT(m1->getAtomWithIdx(0)->getQuery()->getDescription() ==
              "AtomType");
  // query should be for aliphatic C:
  smi = "C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(m1->getAtomWithIdx(0)->Match(m2->getAtomWithIdx(0)));
  delete m2;
  smi = "O";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(!m1->getAtomWithIdx(0)->Match(m2->getAtomWithIdx(0)));
  delete m2;
  smi = "c1ccccc1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(!m1->getAtomWithIdx(0)->Match(m2->getAtomWithIdx(0)));
  delete m1;
  delete m2;

  // atom lists:
  smi = "[C,N,O]";
  m1 = SmartsToMol(smi);
  TEST_ASSERT(m1);
  MolPickler::pickleMol(*m1, pickle);
  delete m1;
  m1 = new ROMol();
  MolPickler::molFromPickle(pickle, *m1);
  TEST_ASSERT(m1->getNumAtoms() == 1);
  TEST_ASSERT(m1->getAtomWithIdx(0)->hasQuery());
  TEST_ASSERT(m1->getAtomWithIdx(0)->getQuery()->getDescription() == "AtomOr");
  smi = "C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(m1->getAtomWithIdx(0)->Match(m2->getAtomWithIdx(0)));
  delete m2;
  smi = "O";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(m1->getAtomWithIdx(0)->Match(m2->getAtomWithIdx(0)));
  delete m2;
  smi = "N";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(m1->getAtomWithIdx(0)->Match(m2->getAtomWithIdx(0)));
  delete m2;
  smi = "F";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(!m1->getAtomWithIdx(0)->Match(m2->getAtomWithIdx(0)));
  delete m2;
  smi = "c1nocc1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(!m1->getAtomWithIdx(0)->Match(m2->getAtomWithIdx(0)));
  TEST_ASSERT(!m1->getAtomWithIdx(0)->Match(m2->getAtomWithIdx(1)));
  TEST_ASSERT(!m1->getAtomWithIdx(0)->Match(m2->getAtomWithIdx(2)));
  delete m1;
  delete m2;

  // more complex atom queries:
  smi = "[C&D2,N&D1]";
  m1 = SmartsToMol(smi);
  TEST_ASSERT(m1);
  MolPickler::pickleMol(*m1, pickle);
  delete m1;
  m1 = new ROMol();
  MolPickler::molFromPickle(pickle, *m1);
  TEST_ASSERT(m1->getNumAtoms() == 1);
  TEST_ASSERT(m1->getAtomWithIdx(0)->hasQuery());
  TEST_ASSERT(m1->getAtomWithIdx(0)->getQuery()->getDescription() == "AtomOr");
  smi = "C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(!m1->getAtomWithIdx(0)->Match(m2->getAtomWithIdx(0)));
  delete m2;
  smi = "CC";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(!m1->getAtomWithIdx(0)->Match(m2->getAtomWithIdx(0)));
  delete m2;
  smi = "C(C)C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(m1->getAtomWithIdx(0)->Match(m2->getAtomWithIdx(0)));
  delete m2;
  smi = "C(C)(C)C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(!m1->getAtomWithIdx(0)->Match(m2->getAtomWithIdx(0)));
  delete m2;
  smi = "N";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(!m1->getAtomWithIdx(0)->Match(m2->getAtomWithIdx(0)));
  delete m2;
  smi = "NC";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(m1->getAtomWithIdx(0)->Match(m2->getAtomWithIdx(0)));
  delete m2;
  smi = "N(C)C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(!m1->getAtomWithIdx(0)->Match(m2->getAtomWithIdx(0)));
  delete m1;
  delete m2;

  // basic recursive queries:
  smi = "[$([#6])]";
  m1 = SmartsToMol(smi);
  TEST_ASSERT(m1);
  MolPickler::pickleMol(*m1, pickle);
  delete m1;
  m1 = new ROMol();
  MolPickler::molFromPickle(pickle, *m1);
  TEST_ASSERT(m1->getNumAtoms() == 1);
  TEST_ASSERT(m1->getAtomWithIdx(0)->hasQuery());
  TEST_ASSERT(m1->getAtomWithIdx(0)->getQuery()->getDescription() ==
              "RecursiveStructure");
  smi = "C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(SubstructMatch(*m2, *m1, matchV));
  delete m2;
  smi = "N";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(!SubstructMatch(*m2, *m1, matchV));
  delete m2;
  smi = "c1ccccc1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(SubstructMatch(*m2, *m1, matchV));
  delete m1;
  delete m2;

  smi = "[R2]";
  m1 = SmartsToMol(smi);
  TEST_ASSERT(m1);
  MolPickler::pickleMol(*m1, pickle);
  delete m1;
  m1 = new ROMol();
  MolPickler::molFromPickle(pickle, *m1);
  TEST_ASSERT(m1->getNumAtoms() == 1);
  TEST_ASSERT(m1->getAtomWithIdx(0)->hasQuery());
  TEST_ASSERT(m1->getAtomWithIdx(0)->getQuery()->getDescription() ==
              "AtomInNRings");
  smi = "C1CCC1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(!m1->getAtomWithIdx(0)->Match(m2->getAtomWithIdx(0)));
  delete m2;
  smi = "C12(CCC2)CCC1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(m1->getAtomWithIdx(0)->Match(m2->getAtomWithIdx(0)));
  TEST_ASSERT(!m1->getAtomWithIdx(0)->Match(m2->getAtomWithIdx(1)));
  delete m1;
  delete m2;

  // basic bond queries:
  smi = "[#6][#6]";
  m1 = SmartsToMol(smi);
  TEST_ASSERT(m1);
  MolPickler::pickleMol(*m1, pickle);
  delete m1;
  m1 = new ROMol();
  MolPickler::molFromPickle(pickle, *m1);
  TEST_ASSERT(m1->getNumAtoms() == 2);
  TEST_ASSERT(m1->getAtomWithIdx(0)->hasQuery());
  TEST_ASSERT(m1->getAtomWithIdx(1)->hasQuery());
  TEST_ASSERT(m1->getBondWithIdx(0)->hasQuery());
  TEST_ASSERT(m1->getBondWithIdx(0)->getQuery()->getDescription() ==
              "SingleOrAromaticBond");
  smi = "CC";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(SubstructMatch(*m2, *m1, matchV));
  delete m2;
  smi = "C=C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(!SubstructMatch(*m2, *m1, matchV));
  delete m2;
  smi = "CN";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(!SubstructMatch(*m2, *m1, matchV));
  delete m2;
  smi = "c1ccccc1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(SubstructMatch(*m2, *m1, matchV));
  delete m1;
  delete m2;

  smi = "[#6]-[#6]";
  m1 = SmartsToMol(smi);
  TEST_ASSERT(m1);
  MolPickler::pickleMol(*m1, pickle);
  delete m1;
  m1 = new ROMol();
  MolPickler::molFromPickle(pickle, *m1);
  TEST_ASSERT(m1->getNumAtoms() == 2);
  TEST_ASSERT(m1->getAtomWithIdx(0)->hasQuery());
  TEST_ASSERT(m1->getAtomWithIdx(1)->hasQuery());
  TEST_ASSERT(m1->getBondWithIdx(0)->hasQuery());
  TEST_ASSERT(m1->getBondWithIdx(0)->getQuery()->getDescription() ==
              "BondOrder");
  smi = "CC";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(SubstructMatch(*m2, *m1, matchV));
  delete m2;
  smi = "C=C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(!SubstructMatch(*m2, *m1, matchV));
  delete m2;
  smi = "CN";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(!SubstructMatch(*m2, *m1, matchV));
  delete m2;
  smi = "c1ccccc1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(!SubstructMatch(*m2, *m1, matchV));
  delete m1;
  delete m2;

  smi = "C=C";
  m1 = SmartsToMol(smi);
  TEST_ASSERT(m1);
  MolPickler::pickleMol(*m1, pickle);
  delete m1;
  m1 = new ROMol();
  MolPickler::molFromPickle(pickle, *m1);
  TEST_ASSERT(m1->getNumAtoms() == 2);
  smi = MolToSmarts(*m1);
  TEST_ASSERT(smi == "C=C");
  delete m1;

  smi = "C=[$(C=O)]";
  m1 = SmartsToMol(smi);
  TEST_ASSERT(m1);
  MolPickler::pickleMol(*m1, pickle);
  delete m1;
  m1 = new ROMol();
  MolPickler::molFromPickle(pickle, *m1);
  TEST_ASSERT(m1->getNumAtoms() == 2);
  smi = MolToSmarts(*m1);
  BOOST_LOG(rdErrorLog) << smi << std::endl;
  TEST_ASSERT(smi == "C=[$(C=O)]");
  delete m1;

  smi = "S(=O)(=O)-[C;H2]([F,Br,I,Cl])";
  m1 = SmartsToMol(smi);
  TEST_ASSERT(m1);
  smi = MolToSmarts(*m1);
  MolPickler::pickleMol(*m1, pickle);
  delete m1;
  m1 = new ROMol();
  MolPickler::molFromPickle(pickle, *m1);
  smi2 = MolToSmarts(*m1);
  BOOST_LOG(rdErrorLog) << smi << " " << smi2 << std::endl;
  TEST_ASSERT(smi == smi2);
  delete m1;

  // more complex recursive queries:
  smi = "[$(C)]";
  m1 = SmartsToMol(smi);
  TEST_ASSERT(m1);
  MolPickler::pickleMol(*m1, pickle);
  delete m1;
  m1 = new ROMol();
  MolPickler::molFromPickle(pickle, *m1);
  TEST_ASSERT(m1->getNumAtoms() == 1);
  TEST_ASSERT(m1->getAtomWithIdx(0)->hasQuery());
  TEST_ASSERT(m1->getAtomWithIdx(0)->getQuery()->getDescription() ==
              "RecursiveStructure");
  smi = "C";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(SubstructMatch(*m2, *m1, matchV));
  delete m2;
  smi = "N";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(!SubstructMatch(*m2, *m1, matchV));
  delete m2;
  smi = "c1ccccc1";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(!SubstructMatch(*m2, *m1, matchV));
  delete m1;
  delete m2;

  // more complex recursive queries:
  smi = "[$(C=O);$(C-N)]";
  m1 = SmartsToMol(smi);
  TEST_ASSERT(m1);
  MolPickler::pickleMol(*m1, pickle);
  delete m1;
  m1 = new ROMol();
  MolPickler::molFromPickle(pickle, *m1);
  TEST_ASSERT(m1->getNumAtoms() == 1);
  TEST_ASSERT(m1->getAtomWithIdx(0)->hasQuery());
  TEST_ASSERT(m1->getAtomWithIdx(0)->getQuery()->getDescription() == "AtomAnd");
  smi = "NCO";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(!SubstructMatch(*m2, *m1, matchV));
  delete m2;
  smi = "CNC(=O)";
  m2 = SmilesToMol(smi);
  TEST_ASSERT(m2);
  TEST_ASSERT(SubstructMatch(*m2, *m1, matchV));
  delete m1;
  delete m2;

  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void testRadicals() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing radical handling in pickles." << std::endl;

  std::string smi = "C[CH]C";
  ROMol *m1 = SmilesToMol(smi);
  TEST_ASSERT(m1);
  TEST_ASSERT(m1->getAtomWithIdx(1)->getNumRadicalElectrons() == 1);
  std::string pickle;
  ROMol *m2;
  MolPickler::pickleMol(*m1, pickle);
  m2 = new ROMol();
  MolPickler::molFromPickle(pickle, *m2);
  TEST_ASSERT(m1->getNumAtoms() == m2->getNumAtoms());
  TEST_ASSERT(m1->getNumBonds() == m2->getNumBonds());
  TEST_ASSERT(m2->getAtomWithIdx(1)->getNumRadicalElectrons() == 1);
  delete m1;
  delete m2;

  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void testIssue2788233(bool doLong = 0) {
  (void)doLong;
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing sf.net issue 2788233." << std::endl;

  {
    std::string fName = getenv("RDBASE");
    fName += "/Code/GraphMol/test_data/Issue2788233.mol";

    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 2);
    TEST_ASSERT(m->getNumBonds() == 0);

    std::string pickle;
    MolPickler::pickleMol(*m, pickle);
    auto *m2 = new RWMol();
    MolPickler::molFromPickle(pickle, *m2);
    TEST_ASSERT(m2->getNumAtoms() == 2);
    TEST_ASSERT(m2->getNumBonds() == 0);

    delete m;
    delete m2;
  }
  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void testIssue3202580() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing sf.net issue 3202580." << std::endl;

  {
    ROMol *m1 = SmilesToMol("C");
    TEST_ASSERT(m1);
    TEST_ASSERT(feq(m1->getAtomWithIdx(0)->getMass(), 12.011, .001));
    std::string pickle;
    MolPickler::pickleMol(*m1, pickle);
    auto *m2 = new RWMol();
    MolPickler::molFromPickle(pickle, *m2);
    TEST_ASSERT(feq(m2->getAtomWithIdx(0)->getMass(), 12.011, .001));
    delete m1;
    delete m2;
  }
  {
    ROMol *m1 = SmilesToMol("[12CH4]");
    TEST_ASSERT(m1);
    TEST_ASSERT(feq(m1->getAtomWithIdx(0)->getMass(), 12.000, .001));
    std::string pickle;
    MolPickler::pickleMol(*m1, pickle);
    auto *m2 = new RWMol();
    MolPickler::molFromPickle(pickle, *m2);
    TEST_ASSERT(feq(m2->getAtomWithIdx(0)->getMass(), 12.000, .001));
    delete m1;
    delete m2;
  }

  {
    ROMol *m1 = SmilesToMol("[13CH4]");
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getAtomWithIdx(0)->getIsotope() == 13);
    std::cerr << "  mass! " << m1->getAtomWithIdx(0)->getMass() << std::endl;
    TEST_ASSERT(feq(m1->getAtomWithIdx(0)->getMass(), 13.003, .001));
    std::string pickle;
    MolPickler::pickleMol(*m1, pickle);
    auto *m2 = new RWMol();
    MolPickler::molFromPickle(pickle, *m2);
    TEST_ASSERT(feq(m2->getAtomWithIdx(0)->getMass(), 13.003, .001));
    delete m1;
    delete m2;
  }

  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void testIssue3316407() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing sf.net issue 3316407." << std::endl;

  {
    std::string fName = getenv("RDBASE");
    fName += "/Code/GraphMol/FileParsers/test_data/rgroups1.mol";

    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 5);

    std::string pickle;
    MolPickler::pickleMol(*m, pickle);
    auto *m2 = new RWMol();
    MolPickler::molFromPickle(pickle, *m2);
    TEST_ASSERT(m2->getNumAtoms() == 5);
    for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
      TEST_ASSERT(m->getAtomWithIdx(i)->getExplicitValence() ==
                  m2->getAtomWithIdx(i)->getExplicitValence());
      TEST_ASSERT(m->getAtomWithIdx(i)->getImplicitValence() ==
                  m2->getAtomWithIdx(i)->getImplicitValence());
    }
    // part of sf.net issue 285:
    TEST_ASSERT(m2->getAtomWithIdx(3)->hasProp(common_properties::dummyLabel));
    TEST_ASSERT(m2->getAtomWithIdx(4)->hasProp(common_properties::dummyLabel));

    delete m;
    delete m2;
  }

  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void testIssue3496759() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing sf.net issue 3496759." << std::endl;

  {
    ROMol *m1 = SmartsToMol("c1ncncn1");
    TEST_ASSERT(m1);
    std::string smi1 = MolToSmiles(*m1, 1);
    TEST_ASSERT(smi1 == "c1ncncn1");

    std::string pickle;
    MolPickler::pickleMol(*m1, pickle);
    RWMol *m2 = new RWMol(pickle);
    TEST_ASSERT(m2);

    std::string smi2 = MolToSmiles(*m2, 1);
    TEST_ASSERT(smi2 == "c1ncncn1");

    delete m1;
    delete m2;
  }

  {
    std::string fName = getenv("RDBASE");
    fName += "/Data/SmartsLib/RLewis_smarts.txt";
    std::ifstream inf(fName.c_str());
    while (!inf.eof()) {
      std::string inl;
      std::getline(inf, inl);
      if (inl[0] == '#' || inl.size() < 2) continue;
      std::vector<std::string> tokens;
      boost::split(tokens, inl, boost::is_any_of(" \t"));
      // std::cerr << "smarts: " << tokens[0] << std::endl;
      ROMol *m1 = SmartsToMol(tokens[0]);
      TEST_ASSERT(m1);
      std::string smi1 = MolToSmiles(*m1, 1);
      std::string sma1 = MolToSmarts(*m1);

      std::string pickle;
      MolPickler::pickleMol(*m1, pickle);
      RWMol *m2 = new RWMol(pickle);
      TEST_ASSERT(m2);

      std::string smi2 = MolToSmiles(*m2, 1);
      std::string sma2 = MolToSmarts(*m2);

      std::cerr << "smi match: " << smi1 << " " << smi2 << std::endl;
      TEST_ASSERT(smi1 == smi2);
      // std::cerr<<"sma match: "<<sma1<<" "<<sma2<<std::endl;
      TEST_ASSERT(sma1 == sma2);

      delete m1;
      delete m2;
    }
  }

  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void testIssue280() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing sf.net issue 280." << std::endl;

  {
    ROMol *m1 = SmilesToMol("CCC");
    TEST_ASSERT(m1);
    std::string v = "1";
    m1->getAtomWithIdx(0)->setProp(common_properties::molAtomMapNumber, v);
    // ints still work
    m1->getAtomWithIdx(1)->setProp(common_properties::molAtomMapNumber, 2);

    std::string pickle;
    MolPickler::pickleMol(*m1, pickle);
    RWMol *m2 = new RWMol(pickle);
    TEST_ASSERT(m2);
    TEST_ASSERT(
        m2->getAtomWithIdx(0)->hasProp(common_properties::molAtomMapNumber));
    TEST_ASSERT(
        m2->getAtomWithIdx(1)->hasProp(common_properties::molAtomMapNumber));
    int iv;
    m2->getAtomWithIdx(0)->getProp(common_properties::molAtomMapNumber, iv);
    TEST_ASSERT(iv == 1);
    m2->getAtomWithIdx(1)->getProp(common_properties::molAtomMapNumber, iv);
    TEST_ASSERT(iv == 2);
    delete m1;
    delete m2;
  }

  {  // but bogus values are not propagated:
    ROMol *m1 = SmilesToMol("CC");
    TEST_ASSERT(m1);
    std::string v = "foo";
    m1->getAtomWithIdx(0)->setProp(common_properties::molAtomMapNumber, v);

    std::string pickle;
    MolPickler::pickleMol(*m1, pickle);
    RWMol *m2 = new RWMol(pickle);
    TEST_ASSERT(m2);
    TEST_ASSERT(
        !m2->getAtomWithIdx(0)->hasProp(common_properties::molAtomMapNumber));
    delete m1;
    delete m2;
  }

  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void testIssue285() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing sf.net issue 285." << std::endl;

  {
    ROMol *m1 = SmilesToMol("*C");
    TEST_ASSERT(m1);
    std::string v = "R";
    m1->getAtomWithIdx(0)->setProp(common_properties::dummyLabel, v);

    std::string pickle;
    MolPickler::pickleMol(*m1, pickle);
    RWMol *m2 = new RWMol(pickle);
    TEST_ASSERT(m2);
    TEST_ASSERT(m2->getAtomWithIdx(0)->hasProp(common_properties::dummyLabel));
    m2->getAtomWithIdx(0)->getProp(common_properties::dummyLabel, v);
    TEST_ASSERT(v == "R");
    delete m1;
    delete m2;
  }

  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void testAtomResidues() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n";
  BOOST_LOG(rdInfoLog) << "Testing residue information handling on atoms"
                       << std::endl;
  {
    auto *m = new RWMol();

    bool updateLabel = true;
    bool takeOwnership = true;
    m->addAtom(new Atom(6), updateLabel, takeOwnership);
    m->addAtom(new Atom(6), updateLabel, takeOwnership);
    m->addBond(0, 1, Bond::SINGLE);
    m->addAtom(new Atom(6), updateLabel, takeOwnership);
    m->addBond(1, 2, Bond::SINGLE);
    m->addAtom(new Atom(6), updateLabel, takeOwnership);
    m->addBond(2, 3, Bond::SINGLE);

    m->getAtomWithIdx(0)->setMonomerInfo(
        new AtomMonomerInfo(AtomMonomerInfo::OTHER, "m1"));
    m->getAtomWithIdx(1)->setMonomerInfo(new AtomPDBResidueInfo("Ca", 3));
    MolOps::sanitizeMol(*m);

    std::string pkl;
    MolPickler::pickleMol(*m, pkl);
    delete m;
    RWMol *m2 = new RWMol(pkl);
    TEST_ASSERT(m2);
    TEST_ASSERT((m2->getAtomWithIdx(0)->getMonomerInfo()));
    TEST_ASSERT(m2->getAtomWithIdx(0)->getMonomerInfo()->getName() == "m1");
    TEST_ASSERT((m2->getAtomWithIdx(1)->getMonomerInfo()));
    TEST_ASSERT(m2->getAtomWithIdx(1)->getMonomerInfo()->getName() == "Ca");
    TEST_ASSERT(static_cast<const AtomPDBResidueInfo *>(
                    m2->getAtomWithIdx(1)->getMonomerInfo())
                    ->getSerialNumber() == 3);
    TEST_ASSERT(!(m2->getAtomWithIdx(2)->getMonomerInfo()));
    TEST_ASSERT(!(m2->getAtomWithIdx(3)->getMonomerInfo()));
    delete m2;
  }
  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void testGithub149() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n";
  BOOST_LOG(rdInfoLog)
      << "Testing Github issue 149: cannot pickle unsanitized molecules"
      << std::endl;
  {
    ROMol *m = SmilesToMol("CCO", 0, 0);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 3);

    std::string pkl;
    MolPickler::pickleMol(*m, pkl);
    delete m;
    m = new RWMol(pkl);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 3);
    delete m;
  }
  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void testGithub713() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n";
  BOOST_LOG(rdInfoLog) << "Testing Github issue 713: Molecule serialization "
                          "doesn't read/write atomic numbers above 128"
                       << std::endl;
  {
    ROMol *m = SmilesToMol("CC");
    TEST_ASSERT(m);
    m->getAtomWithIdx(0)->setAtomicNum(128);
    m->getAtomWithIdx(1)->setAtomicNum(177);

    std::string pkl;
    MolPickler::pickleMol(*m, pkl);
    delete m;
    m = new RWMol(pkl);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getAtomicNum() == 128);
    TEST_ASSERT(m->getAtomWithIdx(1)->getAtomicNum() == 177);

    delete m;
  }
  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void testPickleProps() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n";
  BOOST_LOG(rdInfoLog) << "Testing pickling of properties" << std::endl;

  std::vector<double> v;
  v.push_back(1234.);
  v.push_back(444.);
  v.push_back(1123.);

  ROMol *m = SmilesToMol("CC");
  m->setProp("double", 1.0);
  m->setProp("int", 100);
  m->setProp("bool", true);
  m->setProp("boolfalse", false);
  m->setProp("dvec", v);

  Atom *a = m->getAtomWithIdx(0);
  a->setProp("double", 1.0);
  a->setProp("int", 100);
  a->setProp("bool", true);
  a->setProp("boolfalse", false);
  a->setProp("dvec", v);
  a->setProp("_private", true);

  Bond *b = m->getBondWithIdx(0);
  b->setProp("double", 1.0);
  b->setProp("int", 100);
  b->setProp("bool", true);
  b->setProp("boolfalse", false);
  b->setProp("dvec", v);
  b->setProp("_private", true);

  std::string pkl;
  {
    MolPickler::pickleMol(*m, pkl, PicklerOps::AllProps);
    RWMol *m2 = new RWMol(pkl);
    TEST_ASSERT(m2);
    TEST_ASSERT(m2->getProp<double>("double") == 1.0);
    TEST_ASSERT(m2->getProp<int>("int") == 100);
    TEST_ASSERT(m2->getProp<bool>("bool") == true);
    TEST_ASSERT(m2->getProp<bool>("boolfalse") == false);

    a = m2->getAtomWithIdx(0);
    TEST_ASSERT(a->getProp<double>("double") == 1.0);
    TEST_ASSERT(a->getProp<int>("int") == 100);
    TEST_ASSERT(a->getProp<bool>("bool") == true);
    TEST_ASSERT(a->getProp<bool>("boolfalse") == false);
    TEST_ASSERT(a->getProp<bool>("_private") == true);

    b = m2->getBondWithIdx(0);
    TEST_ASSERT(b->getProp<double>("double") == 1.0);
    TEST_ASSERT(b->getProp<int>("int") == 100);
    TEST_ASSERT(b->getProp<bool>("bool") == true);
    TEST_ASSERT(b->getProp<bool>("boolfalse") == false);
    TEST_ASSERT(b->getProp<bool>("_private") == true);
    delete m2;
  }

  {
    MolPickler::pickleMol(*m, pkl, PicklerOps::MolProps);
    RWMol *m2 = new RWMol(pkl);
    TEST_ASSERT(m2);
    TEST_ASSERT(m2->getProp<double>("double") == 1.0);
    TEST_ASSERT(m2->getProp<int>("int") == 100);
    TEST_ASSERT(m2->getProp<bool>("bool") == true);
    TEST_ASSERT(m2->getProp<bool>("boolfalse") == false);

    a = m2->getAtomWithIdx(0);
    TEST_ASSERT(!a->hasProp("double"));
    TEST_ASSERT(!a->hasProp("int"));
    TEST_ASSERT(!a->hasProp("bool"));
    TEST_ASSERT(!a->hasProp("boolfalse"));
    TEST_ASSERT(!a->hasProp("_private"));

    b = m2->getBondWithIdx(0);
    TEST_ASSERT(!b->hasProp("double"));
    TEST_ASSERT(!b->hasProp("int"));
    TEST_ASSERT(!b->hasProp("bool"));
    TEST_ASSERT(!b->hasProp("boolfalse"));
    TEST_ASSERT(!b->hasProp("_private"));

    delete m2;
  }

  {
    MolPickler::pickleMol(*m, pkl, PicklerOps::AtomProps);
    RWMol *m2 = new RWMol(pkl);
    TEST_ASSERT(m2);
    TEST_ASSERT(!m2->hasProp("double"));
    TEST_ASSERT(!m2->hasProp("int"));
    TEST_ASSERT(!m2->hasProp("bool"));
    TEST_ASSERT(!m2->hasProp("boolfalse"));

    a = m2->getAtomWithIdx(0);
    TEST_ASSERT(a->getProp<double>("double") == 1.0);
    TEST_ASSERT(a->getProp<int>("int") == 100);
    TEST_ASSERT(a->getProp<bool>("bool") == true);
    TEST_ASSERT(a->getProp<bool>("boolfalse") == false);
    TEST_ASSERT(!a->hasProp("_private"));

    b = m2->getBondWithIdx(0);
    TEST_ASSERT(!b->hasProp("double"));
    TEST_ASSERT(!b->hasProp("int"));
    TEST_ASSERT(!b->hasProp("bool"));
    TEST_ASSERT(!b->hasProp("boolfalse"));
    TEST_ASSERT(!b->hasProp("_private"));
    delete m2;
  }

  {
    MolPickler::pickleMol(*m, pkl,
                          PicklerOps::AtomProps | PicklerOps::PrivateProps);
    RWMol *m2 = new RWMol(pkl);
    TEST_ASSERT(m2);
    TEST_ASSERT(!m2->hasProp("double"));
    TEST_ASSERT(!m2->hasProp("int"));
    TEST_ASSERT(!m2->hasProp("bool"));
    TEST_ASSERT(!m2->hasProp("boolfalse"));

    a = m2->getAtomWithIdx(0);
    TEST_ASSERT(a->getProp<double>("double") == 1.0);
    TEST_ASSERT(a->getProp<int>("int") == 100);
    TEST_ASSERT(a->getProp<bool>("bool") == true);
    TEST_ASSERT(a->getProp<bool>("boolfalse") == false);
    TEST_ASSERT(a->getProp<bool>("_private") == true);

    b = m2->getBondWithIdx(0);
    TEST_ASSERT(!b->hasProp("double"));
    TEST_ASSERT(!b->hasProp("int"));
    TEST_ASSERT(!b->hasProp("bool"));
    TEST_ASSERT(!b->hasProp("boolfalse"));
    TEST_ASSERT(!b->hasProp("_private"));

    delete m2;
  }

  {
    MolPickler::pickleMol(*m, pkl, PicklerOps::BondProps);
    RWMol *m2 = new RWMol(pkl);
    TEST_ASSERT(m2);
    TEST_ASSERT(!m2->hasProp("double"));
    TEST_ASSERT(!m2->hasProp("int"));
    TEST_ASSERT(!m2->hasProp("bool"));
    TEST_ASSERT(!m2->hasProp("boolfalse"));

    a = m2->getAtomWithIdx(0);
    TEST_ASSERT(!a->hasProp("double"));
    TEST_ASSERT(!a->hasProp("int"));
    TEST_ASSERT(!a->hasProp("bool"));
    TEST_ASSERT(!a->hasProp("boolfalse"));
    TEST_ASSERT(!a->hasProp("_private"));

    b = m2->getBondWithIdx(0);
    TEST_ASSERT(b->getProp<double>("double") == 1.0);
    TEST_ASSERT(b->getProp<int>("int") == 100);
    TEST_ASSERT(b->getProp<bool>("bool") == true);
    TEST_ASSERT(b->getProp<bool>("boolfalse") == false);
    TEST_ASSERT(!b->hasProp("_private"));
    delete m2;
  }

  {
    MolPickler::pickleMol(*m, pkl,
                          PicklerOps::BondProps | PicklerOps::PrivateProps);
    RWMol *m2 = new RWMol(pkl);
    TEST_ASSERT(m2);
    TEST_ASSERT(!m2->hasProp("double"));
    TEST_ASSERT(!m2->hasProp("int"));
    TEST_ASSERT(!m2->hasProp("bool"));
    TEST_ASSERT(!m2->hasProp("boolfalse"));

    a = m2->getAtomWithIdx(0);
    TEST_ASSERT(!a->hasProp("double"));
    TEST_ASSERT(!a->hasProp("int"));
    TEST_ASSERT(!a->hasProp("bool"));
    TEST_ASSERT(!a->hasProp("boolfalse"));
    TEST_ASSERT(!a->hasProp("_private"));

    b = m2->getBondWithIdx(0);
    TEST_ASSERT(b->getProp<double>("double") == 1.0);
    TEST_ASSERT(b->getProp<int>("int") == 100);
    TEST_ASSERT(b->getProp<bool>("bool") == true);
    TEST_ASSERT(b->getProp<bool>("boolfalse") == false);
    TEST_ASSERT(b->getProp<bool>("_private") == true);
    delete m2;
  }

  delete m;

  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void testGithub1563() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n";
  BOOST_LOG(rdInfoLog)
      << "Testing Github 1563: Add a canned Atom query for heavy atom degree"
      << std::endl;
  RWMol m;
  auto *qa = new QueryAtom();
  auto *query = makeAtomHeavyAtomDegreeQuery(3);
  qa->setQuery(query);
  m.addAtom(qa);
  TEST_ASSERT(m.getAtomWithIdx(0)->hasQuery());
  TEST_ASSERT(m.getAtomWithIdx(0)->getQuery()->getDescription() ==
              "AtomHeavyAtomDegree");
  RWMol nm(m);
  TEST_ASSERT(nm.getAtomWithIdx(0)->hasQuery());
  TEST_ASSERT(nm.getAtomWithIdx(0)->getQuery()->getDescription() ==
              "AtomHeavyAtomDegree");
  std::string pkl;
  MolPickler::pickleMol(m, pkl);
  RWMol nm2(pkl);
  TEST_ASSERT(nm2.getAtomWithIdx(0)->hasQuery());
  TEST_ASSERT(nm2.getAtomWithIdx(0)->getQuery()->getDescription() ==
              "AtomHeavyAtomDegree");
  delete qa;
  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void testGithub1710() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n";
  BOOST_LOG(rdInfoLog) << "Testing Github 1710: bonds that are STEREOCIS or "
                          "STEREOTRANS cannot be depickled"
                       << std::endl;

  RWMol m;
  m.addAtom(new Atom(6), false, true);
  m.addAtom(new Atom(6), false, true);
  m.addAtom(new Atom(6), false, true);
  m.addAtom(new Atom(6), false, true);
  m.addBond(0, 1, Bond::SINGLE);
  m.addBond(1, 2, Bond::DOUBLE);
  m.addBond(2, 3, Bond::SINGLE);
  m.getBondWithIdx(1)->setStereoAtoms(0, 3);

  {
    m.getBondWithIdx(1)->setStereo(Bond::STEREOCIS);
    std::string pkl;
    MolPickler::pickleMol(m, pkl);
    RWMol *m2 = new RWMol(pkl);
    TEST_ASSERT(m2->getBondWithIdx(1)->getStereo() == Bond::STEREOCIS);
    TEST_ASSERT(m2->getBondWithIdx(1)->getStereoAtoms()[0] == 0);
    TEST_ASSERT(m2->getBondWithIdx(1)->getStereoAtoms()[1] == 3);
    delete m2;
  }
  {
    m.getBondWithIdx(1)->setStereo(Bond::STEREOTRANS);
    std::string pkl;
    MolPickler::pickleMol(m, pkl);
    RWMol *m2 = new RWMol(pkl);
    TEST_ASSERT(m2->getBondWithIdx(1)->getStereo() == Bond::STEREOTRANS);
    TEST_ASSERT(m2->getBondWithIdx(1)->getStereoAtoms()[0] == 0);
    TEST_ASSERT(m2->getBondWithIdx(1)->getStereoAtoms()[1] == 3);
    delete m2;
  }
  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void testGithub1999() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n";
  BOOST_LOG(rdInfoLog) << "Testing Github 1999: pickling atom map numbers >127 "
                          "gives negative values"
                       << std::endl;

  {
    RWMol m;
    m.addAtom(new Atom(6), false, true);
    m.getAtomWithIdx(0)->setAtomMapNum(777);
    std::string pkl;
    MolPickler::pickleMol(m, pkl);
    std::unique_ptr<RWMol> m2(new RWMol(pkl));
    TEST_ASSERT(m2->getAtomWithIdx(0)->getAtomMapNum() == 777);
  }
  {
    RWMol m;
    m.addAtom(new Atom(6), false, true);
    m.getAtomWithIdx(0)->setAtomMapNum(128);
    std::string pkl;
    MolPickler::pickleMol(m, pkl);
    std::unique_ptr<RWMol> m2(new RWMol(pkl));
    TEST_ASSERT(m2->getAtomWithIdx(0)->getAtomMapNum() == 128);
  }
  BOOST_LOG(rdErrorLog) << "\tdone" << std::endl;
}

void testEnhancedStereoChemistry() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n";
  BOOST_LOG(rdInfoLog) << "Testing that enhanced stereochemistry is pickled"
                       << std::endl;
  RWMol m;
  for (unsigned i = 0u; i < 6; ++i) {
    m.addAtom(new Atom(6), false, true);
  }

  // add 3 stereo groups
  {
    std::vector<StereoGroup> groups;
    std::vector<Atom *> atoms0 = {{m.getAtomWithIdx(0), m.getAtomWithIdx(1)}};
    groups.emplace_back(RDKit::StereoGroupType::STEREO_ABSOLUTE,
                        std::move(atoms0));
    std::vector<Atom *> atoms1 = {{m.getAtomWithIdx(2), m.getAtomWithIdx(3)}};
    groups.emplace_back(RDKit::StereoGroupType::STEREO_OR, std::move(atoms1));
    std::vector<Atom *> atoms2 = {{m.getAtomWithIdx(4), m.getAtomWithIdx(5)}};
    groups.emplace_back(RDKit::StereoGroupType::STEREO_AND, std::move(atoms2));
    m.setStereoGroups(std::move(groups));
  }

  std::string pkl;
  MolPickler::pickleMol(m, pkl);

  std::unique_ptr<RWMol> roundTripped(new RWMol(pkl));

  auto &ref_groups = m.getStereoGroups();
  auto &new_groups = roundTripped->getStereoGroups();
  TEST_ASSERT(ref_groups.size() == new_groups.size());
  for (unsigned i = 0u; i < 3; ++i) {
    TEST_ASSERT(ref_groups[i].getGroupType() == new_groups[i].getGroupType());
  }
}

int main(int argc, char *argv[]) {
  RDLog::InitLogs();
  bool doLong = false;
  if (argc > 1) {
    if (!strncmp(argv[1], "-l", 2)) {
      doLong = true;
    }
  }

//_createPickleFile();
#if 1
  test1(doLong);
  // test2(doLong);
  test3(doLong);
  test4();
  testIssue164();
  testIssue219();
  testIssue220();
  // timeTest(doLong);
  testQueries();
  testRadicals();
  testPickleProps();
#endif
  testIssue2788233();
  testIssue3202580();
  testIssue3316407();
  testIssue3496759();
  testIssue280();
  testIssue285();
  testAtomResidues();
  testGithub149();
  testGithub713();
  testGithub1563();
  testGithub1710();
  testGithub1999();
  testEnhancedStereoChemistry();
}
