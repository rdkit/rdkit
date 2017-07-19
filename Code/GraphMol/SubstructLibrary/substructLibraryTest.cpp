//
//  Copyright (C) 2017 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//

// std bits
#include <iostream>

// RD bits
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SubstructLibrary/SubstructLibrary.h>

#include <GraphMol/Substruct/SubstructMatch.h>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>

using namespace RDKit;

namespace {
void runTest(SubstructLibrary &ssslib, const ROMol &pattern, int nThreads) {
  std::vector<unsigned int> libMatches = ssslib.getMatches(pattern, nThreads);
  boost::dynamic_bitset<> hasMatch(ssslib.size());
  BOOST_FOREACH (unsigned int idx, libMatches) { hasMatch[idx] = 1; }

  for (unsigned int i = 0; i < ssslib.size(); ++i) {
    MatchVectType match;
    bool matched = SubstructMatch(*ssslib.getMol(i), pattern, match);
    // std::cerr << MolToSmiles(*ssslib.getMol(i), true) << " " << hasMatch[i]
    //           << " " << matched << std::endl;
    TEST_ASSERT(hasMatch[i] == matched);
  }
};
}

void test1() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test1" << std::endl;

  std::string fName = getenv("RDBASE");
  fName += "/Data/NCI/first_200.props.sdf";
  SDMolSupplier suppl(fName);
  SubstructLibrary ssslib;
  while (!suppl.atEnd()) {
    ROMol *mol = 0;
    try {
      mol = suppl.next();
    } catch (...) {
      continue;
    }
    if (!mol) continue;
    ssslib.addMol(*mol);
  }

  {
    ROMol *query = SmartsToMol("[#6;$([#6]([#6])[!#6])]");
    runTest(ssslib, *query, 1);
#ifdef RDK_TEST_MULTITHREADED
    runTest(ssslib, *query, -1);
#endif
    delete query;
  }

  {
    ROMol *query = SmartsToMol("[$([O,S]-[!$(*=O)])]");
    runTest(ssslib, *query, 1);
#ifdef RDK_TEST_MULTITHREADED
    runTest(ssslib, *query, -1);
#endif
    delete query;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void test2() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test2" << std::endl;

  std::string fName = getenv("RDBASE");
  fName += "/Data/NCI/first_200.props.sdf";
  SDMolSupplier suppl(fName);
  MolHolder *mols = new MolHolder();
  PatternHolder *fps = new PatternHolder();
  boost::shared_ptr<MolHolder> mols_ptr(mols);
  boost::shared_ptr<PatternHolder> fps_ptr(fps);

  SubstructLibrary ssslib(mols_ptr, fps_ptr);
  while (!suppl.atEnd()) {
    ROMol *mol = 0;
    try {
      mol = suppl.next();
    } catch (...) {
      continue;
    }
    if (!mol) continue;
    ssslib.addMol(*mol);
  }

  {
    ROMol *query = SmartsToMol("[#6]([#6])[!#6]");
    runTest(ssslib, *query, 1);
#ifdef RDK_TEST_MULTITHREADED
    runTest(ssslib, *query, -1);
#endif
    delete query;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void test3() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test3 (stereo options)" << std::endl;

  SubstructLibrary ssslib(boost::make_shared<MolHolder>());
  for(int i=0;i<10;++i) {
    ROMol *m1 = SmilesToMol("C1CCO[C@@](N)(O)1");
    ROMol *m2 = SmilesToMol("C1CCO[C@](N)(O)1");
    ROMol *m3 = SmilesToMol("C1CCO[C@@](O)(N)1");
    ROMol *m4 = SmilesToMol("C1CCO[C@](O)(N)1");
    ssslib.addMol( *m1 );
    ssslib.addMol( *m2 );
    ssslib.addMol( *m3 );
    ssslib.addMol( *m4 );
    delete m1;
    delete m2;
    delete m3;
    delete m4;
  }

  ROMol *query = SmartsToMol("C-1-C-C-O-C(-[O])(-[N])1");
  std::vector<unsigned int> res = ssslib.getMatches(*query, true, false);
  TEST_ASSERT(res.size() == 40);

  delete query;
  query = SmartsToMol("C-1-C-C-O-[C@@](-[O])(-[N])1")          ;

  res = ssslib.getMatches(*query, true, true);
  TEST_ASSERT(res.size() == 20);

  res = ssslib.getMatches(*query, true, false);
  TEST_ASSERT(res.size() == 40);

  delete query;
  BOOST_LOG(rdErrorLog) << "    Done (stereo options)" << std::endl;

}

void test4() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test4 (trusted smiles)" << std::endl;

  boost::shared_ptr<CachedSmilesMolHolder> holder = boost::make_shared<CachedSmilesMolHolder>();
  SubstructLibrary ssslib(holder);
  
  for(int i=0;i<10;++i) {
    holder->addSmiles("C1CCO[C@@](N)(O)1");
    holder->addSmiles("C1CCO[C@](N)(O)1");
    holder->addSmiles("C1CCO[C@@](O)(N)1");
    holder->addSmiles("C1CCO[C@](O)(N)1");
  }

  ROMol *query = SmartsToMol("C-1-C-C-O-C(-[O])(-[N])1");
  std::vector<unsigned int> res = ssslib.getMatches(*query, true, false);
  TEST_ASSERT(res.size() == 40);

  delete query;
  query = SmartsToMol("C-1-C-C-O-[C@@](-[O])(-[N])1")          ;

  res = ssslib.getMatches(*query, true, true);
  TEST_ASSERT(res.size() == 20);

  res = ssslib.getMatches(*query, true, false);
  TEST_ASSERT(res.size() == 40);

  delete query;
  BOOST_LOG(rdErrorLog) << "    Done (stereo options)" << std::endl;

}

int main(int argc, char *argv[]) {
#if 1
  test1();
  test2();
  test3();
  test4();
#endif
  return 0;
}
