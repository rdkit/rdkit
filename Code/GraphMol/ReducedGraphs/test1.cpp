// $Id$
//
//  Copyright (C) 2013 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/test.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/ReducedGraphs/ReducedGraphs.h>
#include <DataStructs/ExplicitBitVect.h>

#include <RDGeneral/RDLog.h>
#include <string>

using namespace RDKit;

void test1() {
  BOOST_LOG(rdInfoLog) << "testing basics" << std::endl;
  {
    std::string smi = "c1ccccc1CCO";
    ROMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);

    RDNumeric::DoubleVector *fp = ReducedGraphs::getErGFingerprint(*m1);
    std::cerr << *fp << std::endl;

    delete fp;
    delete m1;
  }

  {
    std::string smi = "c1cnccc1CCO";
    ROMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);

    RDNumeric::DoubleVector *fp = ReducedGraphs::getErGFingerprint(*m1);
    std::cerr << *fp << std::endl;

    delete fp;
    delete m1;
  }

  {
    std::string smi = "OCCC1=CC2=C(C=CC=C2)C=C1";
    ROMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);

    RDNumeric::DoubleVector *fp = ReducedGraphs::getErGFingerprint(*m1);
    std::cerr << *fp << std::endl;

    delete fp;
    delete m1;
  }

  {
    std::string smi = "OCCC1=CC2=C(C=CC=C2)N=C1";
    ROMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);

    RDNumeric::DoubleVector *fp = ReducedGraphs::getErGFingerprint(*m1);
    std::cerr << *fp << std::endl;

    delete fp;
    delete m1;
  }
  {
    std::string smi = "OCCC1=CC2=C(CCCC2)N=C1";
    ROMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);

    RDNumeric::DoubleVector *fp = ReducedGraphs::getErGFingerprint(*m1);
    std::cerr << *fp << std::endl;

    delete fp;
    delete m1;
  }
  {
    std::string smi = "OCCc1ccc(CC(C)C)cc1";
    ROMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);

    RDNumeric::DoubleVector *fp = ReducedGraphs::getErGFingerprint(*m1);
    std::cerr << *fp << std::endl;

    delete fp;
    delete m1;
  }
  {
    std::string smi = "OCCC1CCC(CC(C)C)CC1";
    ROMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);

    RDNumeric::DoubleVector *fp = ReducedGraphs::getErGFingerprint(*m1);
    std::cerr << *fp << std::endl;

    delete fp;
    delete m1;
  }
  {
    std::string smi = "OCCC1=CCC(CC(C)C)CC1";
    ROMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);

    RDNumeric::DoubleVector *fp = ReducedGraphs::getErGFingerprint(*m1);
    std::cerr << *fp << std::endl;

    delete fp;
    delete m1;
  }
  {
    std::string smi = "OCCC1=CCC(CC(C)C)C=C1";
    ROMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);

    RDNumeric::DoubleVector *fp = ReducedGraphs::getErGFingerprint(*m1);
    std::cerr << *fp << std::endl;

    delete fp;
    delete m1;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void test2() {
  BOOST_LOG(rdInfoLog) << "testing basics2" << std::endl;
  {
    std::string smi1 = "c1ccccc1CCO";
    std::string smi2 = "C1=CC=CCC1CCO";
    ROMol *m1 = SmilesToMol(smi1);
    TEST_ASSERT(m1);
    ROMol *m2 = SmilesToMol(smi2);
    TEST_ASSERT(m2);

    RDNumeric::DoubleVector *fp1 = ReducedGraphs::getErGFingerprint(*m1);
    RDNumeric::DoubleVector *fp2 = ReducedGraphs::getErGFingerprint(*m2);

    TEST_ASSERT(feq(RDNumeric::TanimotoSimilarity(*fp1, *fp2), 1.0, 0.001));

    delete fp1;
    delete fp2;
    delete m1;
    delete m2;
  }
  {
    std::string smi1 = "c1ccccc1CCO";
    std::string smi2 = "C1CC=CCC1CCO";
    ROMol *m1 = SmilesToMol(smi1);
    TEST_ASSERT(m1);
    ROMol *m2 = SmilesToMol(smi2);
    TEST_ASSERT(m2);

    RDNumeric::DoubleVector *fp1 = ReducedGraphs::getErGFingerprint(*m1);
    RDNumeric::DoubleVector *fp2 = ReducedGraphs::getErGFingerprint(*m2);

    TEST_ASSERT(!feq(RDNumeric::TanimotoSimilarity(*fp1, *fp2), 1.0, 0.001));

    delete fp1;
    delete fp2;
    delete m1;
    delete m2;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testCanRetrieveProp() {
  BOOST_LOG(rdInfoLog) << "testing retrieving _ErGAtomTypes from property"
                       << std::endl;
  auto m = "OCCc1ccccc1"_smiles;
  std::vector<std::vector<int>> expected{{0, 1}, {}, {}, {}, {5}};
  std::vector<std::vector<int>> res;
  std::unique_ptr<ROMol> mrg(
      ReducedGraphs::generateMolExtendedReducedGraph(*m));
  for (const auto atom : mrg->atoms()) {
    std::vector<int> atomTypes;
    TEST_ASSERT(atom->getPropIfPresent("_ErGAtomTypes", atomTypes));
    res.push_back(atomTypes);
  }
  TEST_ASSERT(res == expected);
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;
  RDLog::InitLogs();
  test1();
  test2();
  testCanRetrieveProp();
  return 0;
}
