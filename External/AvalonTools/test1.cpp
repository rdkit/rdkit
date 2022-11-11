//
//  Copyright (C) 2008-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
//  Created by Greg Landrum, July 2008
//

//
//  Expected test results here correspond to v1.0 of the open-source
//  avalontoolkit
//

#include <RDGeneral/test.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <RDGeneral/Invariant.h>
#include <DataStructs/ExplicitBitVect.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include "AvalonTools.h"

#include <string>

using namespace RDKit;

void test1() {
  BOOST_LOG(rdInfoLog) << "testing canonical smiles generation" << std::endl;

  {
    ROMol *m = static_cast<ROMol *>(SmilesToMol("c1ccccc1"));
    TEST_ASSERT(m);
    std::string smi = AvalonTools::getCanonSmiles(*m);
    TEST_ASSERT(smi == "c1ccccc1");
    delete m;
  }
  {
    ROMol *m = static_cast<ROMol *>(SmilesToMol("c1cccnc1"));
    TEST_ASSERT(m);
    std::string smi = AvalonTools::getCanonSmiles(*m);
    TEST_ASSERT(smi == "c1ccncc1");
    delete m;
  }
  {
    ROMol *m = static_cast<ROMol *>(SmilesToMol("n1ccccc1"));
    TEST_ASSERT(m);
    std::string smi = AvalonTools::getCanonSmiles(*m);
    TEST_ASSERT(smi == "c1ccncc1");
    delete m;
  }
  {
    std::string smi = AvalonTools::getCanonSmiles("n1ccccc1", true);
    TEST_ASSERT(smi == "c1ccncc1");
  }
  {
    std::string smi = AvalonTools::getCanonSmiles("c1cccnc1", true);
    TEST_ASSERT(smi == "c1ccncc1");
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}
void test2() {
  BOOST_LOG(rdInfoLog) << "testing coordinate generation" << std::endl;

#if 1
  {
    RWMol *m = SmilesToMol("c1cccnc1");
    TEST_ASSERT(m);
    unsigned int confId = AvalonTools::set2DCoords(*m);
    TEST_ASSERT(m->getNumConformers() == 1);
    TEST_ASSERT(confId == 0);
    delete m;
  }
#endif
  {
    std::string molb = AvalonTools::set2DCoords("c1cccnc1", true);
    TEST_ASSERT(molb != "");
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void test3() {
  BOOST_LOG(rdInfoLog) << "testing fingerprint generation" << std::endl;

  {
    ROMol *m = static_cast<ROMol *>(SmilesToMol("c1ccccn1"));
    TEST_ASSERT(m);
    ExplicitBitVect bv(512);
    AvalonTools::getAvalonFP(*m, bv, 512, false, true, 0x00006FFF);
    BOOST_LOG(rdInfoLog) << "c1ccccn1 " << bv.getNumOnBits() << std::endl;
    TEST_ASSERT(bv.getNumOnBits() == 18);
    delete m;
  }
  {
    ROMol *m = static_cast<ROMol *>(SmilesToMol("c1ccccc1"));
    TEST_ASSERT(m);
    ExplicitBitVect bv(512);
    AvalonTools::getAvalonFP(*m, bv, 512, false, true, 0x006FFF);
    BOOST_LOG(rdInfoLog) << "c1ccccn1 " << bv.getNumOnBits() << std::endl;
    TEST_ASSERT(bv.getNumOnBits() == 6);
    delete m;
  }
  {
    ROMol *m = static_cast<ROMol *>(SmilesToMol("c1nnccc1"));
    TEST_ASSERT(m);
    ExplicitBitVect bv(512);
    AvalonTools::getAvalonFP(*m, bv, 512, false, true, 0x006FFF);
    BOOST_LOG(rdInfoLog) << "c1nnccc1 " << bv.getNumOnBits() << std::endl;
    TEST_ASSERT(bv.getNumOnBits() == 28);
    delete m;
  }
  {
    ROMol *m = static_cast<ROMol *>(SmilesToMol("c1ncncc1"));
    TEST_ASSERT(m);
    ExplicitBitVect bv(512);
    AvalonTools::getAvalonFP(*m, bv, 512, false, true, 0x006FFF);
    BOOST_LOG(rdInfoLog) << "c1ncncc1 " << bv.getNumOnBits() << std::endl;
    TEST_ASSERT(bv.getNumOnBits() == 25);
    delete m;
  }
  {
    ExplicitBitVect bv(512);
    AvalonTools::getAvalonFP("c1cccnc1", true, bv, 512, false, true, 0x006FFF);
    BOOST_LOG(rdInfoLog) << "c1cccnc1 " << bv.getNumOnBits() << std::endl;
    TEST_ASSERT(bv.getNumOnBits() == 18);
  }
  {
    ExplicitBitVect bv(512);
    AvalonTools::getAvalonFP("c1ccccc1", true, bv, 512, false, true, 0x006FFF);
    BOOST_LOG(rdInfoLog) << "c1ccccc1 " << bv.getNumOnBits() << std::endl;
    TEST_ASSERT(bv.getNumOnBits() == 6);
  }

  {
    ROMol *m = static_cast<ROMol *>(SmilesToMol("c1cccnc1"));
    TEST_ASSERT(m);
    ExplicitBitVect bv(1024);
    AvalonTools::getAvalonFP(*m, bv, 1024, false, true, 0x006FFF);
    BOOST_LOG(rdInfoLog) << "c1cccnc1 " << bv.getNumOnBits() << std::endl;
    TEST_ASSERT(bv.getNumOnBits() == 19);
    delete m;
  }
  {
    ExplicitBitVect bv(2048);
    AvalonTools::getAvalonFP("c1cocc1", true, bv, 2048, false, true, 0x006FFF);
    BOOST_LOG(rdInfoLog) << "c1cocc1 " << bv.getNumOnBits() << std::endl;
    TEST_ASSERT(bv.getNumOnBits() == 53);
  }
  {
    ExplicitBitVect bv(2048);
    AvalonTools::getAvalonFP("C1=COC=C1", true, bv, 2048, false, true,
                             0x006FFF);
    BOOST_LOG(rdInfoLog) << "C1=COC=C1 " << bv.getNumOnBits() << std::endl;
    TEST_ASSERT(bv.getNumOnBits() == 53);
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGitHub1062() {
  BOOST_LOG(rdInfoLog) << "testing GitHub issue #1062:  pyAvalonTools not "
                          "preserving the stereochemistry of double bonds when "
                          "computing 2D coordinates"
                       << std::endl;
  {
    std::string s0 = "C/C=C\\C";
    ROMol *m1 = SmilesToMol(s0);
    auto s1 = MolToSmiles(*m1);
    AvalonTools::set2DCoords(*m1);
    std::string mb = MolToMolBlock(*m1);
    ROMol *m2 = MolBlockToMol(mb);
    std::string s2 = MolToSmiles(*m2);
    delete m1;
    delete m2;
    TEST_ASSERT(s1 == s2);
  }
  {
    // repeat the test with an input smiles that is not canonical
    // to verify that the implementation is not sensitive to the
    // ordering of atoms
    std::string s0 = "C/C=C(F)\\C";
    ROMol *m1 = SmilesToMol(s0);
    auto s1 = MolToSmiles(*m1);
    TEST_ASSERT(s1 != s0);
    AvalonTools::set2DCoords(*m1);
    std::string mb = MolToMolBlock(*m1);
    ROMol *m2 = MolBlockToMol(mb);
    std::string s2 = MolToSmiles(*m2);
    delete m1;
    delete m2;
    TEST_ASSERT(s1 == s2);
  }
}

void testRDK151() {
  BOOST_LOG(rdInfoLog) << "testing Jira issue RDK-151:  pyAvalonTools not "
                          "generating chiral smiles from molecules"
                       << std::endl;

  {
    std::string tSmi = "C[C@H](F)Cl";
    auto *m = static_cast<ROMol *>(SmilesToMol(tSmi));
    TEST_ASSERT(m);
    std::string smi = AvalonTools::getCanonSmiles(tSmi, true);
    CHECK_INVARIANT(smi == tSmi, smi + "!=" + tSmi);
    smi = AvalonTools::getCanonSmiles(*m);
    CHECK_INVARIANT(smi == tSmi, smi + "!=" + tSmi);

    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testSmilesFailures() {
  BOOST_LOG(rdInfoLog) << "testing handling of bad smiles strings" << std::endl;

  {
    std::string tSmi = "C1C";
    std::string smi = AvalonTools::getCanonSmiles(tSmi, true);
    CHECK_INVARIANT(smi == "", smi);
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testSubstructFps() {
  BOOST_LOG(rdInfoLog) << "testing substructure fingerprints " << std::endl;
  {
    ExplicitBitVect bv1(512), bv2(512);
    AvalonTools::getAvalonFP("c1ccccc1", true, bv1, 512, true, true,
                             AvalonTools::avalonSSSBits);
    AvalonTools::getAvalonFP("c1ccccc1C(F)(F)F", true, bv2, 512);
    TEST_ASSERT((bv1 & bv2) == bv1);
    AvalonTools::getAvalonFP("c1ccccc1C(F)(F)F", true, bv1, 512);
    TEST_ASSERT((bv1 & bv2) == bv1);
    AvalonTools::getAvalonFP("c1cccc(C)c1C(F)(F)F", true, bv2, 512);
    TEST_ASSERT((bv1 & bv2) == bv1);
  }
  {
    ExplicitBitVect bv1(512), bv2(512);
    AvalonTools::getAvalonFP("c1ccccc1O", true, bv1, 512, true, true,
                             AvalonTools::avalonSSSBits);
    AvalonTools::getAvalonFP("c1ccccc1OC", true, bv2, 512);
    TEST_ASSERT((bv1 & bv2) == bv1);
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testStruChk() {
  BOOST_LOG(rdInfoLog) << "testing structure checking " << std::endl;
  {
    int errs = 0;
    RDKit::ROMOL_SPTR m = AvalonTools::checkMol(errs, "c1ccccc1", true);
    TEST_ASSERT(errs == 0);
    m = AvalonTools::checkMol(errs, "c1c(R)cccc1C1(CC-C(C)C1)C", true);
    TEST_ASSERT(errs != 0);
  }
  {
    int errs = 0;
    std::string res;
    boost::tie(res, errs) = AvalonTools::checkMolString("c1ccccc1", true);
    TEST_ASSERT(errs == 0);
    TEST_ASSERT(res != "");
    boost::tie(res, errs) =
        AvalonTools::checkMolString("c1c(R)cccc1C1(CC-C(C)C1)C", true);
    TEST_ASSERT(errs == 1);
    TEST_ASSERT(res == "");
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testBadMolfile() {
  BOOST_LOG(rdInfoLog) << "testing handling bad molecules " << std::endl;
  // some tests around dealing with bad mol blocks
  {
    std::string molb =
        "SNAP007157A\n\
  MACCS-II3194121345\n\
\n\
0    0  0  0  0";
    std::string smi = AvalonTools::getCanonSmiles(molb, false);
    CHECK_INVARIANT(smi == "", smi);

    ExplicitBitVect bv(1024);
    AvalonTools::getAvalonFP(molb, false, bv, 1024);
    TEST_ASSERT(bv.getNumOnBits() == 0);

    std::string oMolb;
    AvalonTools::set2DCoords(molb, false);
    CHECK_INVARIANT(oMolb == "", oMolb);
  }
}

void testSmilesSegFault() {
  BOOST_LOG(rdInfoLog)
      << "testing a canonical smiles case that led to seg faults " << std::endl;
  // some tests around dealing with bad mol blocks
  {
    std::string inSmi(1024, 'C');
    std::string smi = AvalonTools::getCanonSmiles(inSmi, true);
    TEST_ASSERT(smi == inSmi);
  }
  {
    std::string inSmi(1534, 'C');
    std::string smi = AvalonTools::getCanonSmiles(inSmi, true);
    TEST_ASSERT(smi == inSmi);
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub336() {
  BOOST_LOG(rdInfoLog) << "testing github issue 336: bad canonical smiles for "
                          "conjugated double bonds"
                       << std::endl;
  // some tests around dealing with bad mol blocks
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/External/AvalonTools/test_data/";
    std::ifstream ins((pathName + "EZ_test.2.sdf").c_str());
    std::string mb((std::istreambuf_iterator<char>(ins)),
                   std::istreambuf_iterator<char>());
    ROMol *m = MolBlockToMol(mb);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 17);

    std::string smi1 = AvalonTools::getCanonSmiles(mb, false);
    std::string smi2 = AvalonTools::getCanonSmiles(*m);
    std::cerr << "smi1: " << smi1 << std::endl;
    std::cerr << "smi2: " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);
    delete m;
  }

  {
    std::string pathName = getenv("RDBASE");
    pathName += "/External/AvalonTools/test_data/";
    std::ifstream ins((pathName + "heterocycle.mol").c_str());
    std::string mb((std::istreambuf_iterator<char>(ins)),
                   std::istreambuf_iterator<char>());
    RWMol *m = MolBlockToMol(mb, false);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 6);
    m->updatePropertyCache();
    MolOps::cleanUp(*m);
    MolOps::setAromaticity(*m);

    std::string smi1 = AvalonTools::getCanonSmiles(mb, false);
    std::string smi2 = AvalonTools::getCanonSmiles(*m);
    std::cerr << "smi1: " << smi1 << std::endl;
    std::cerr << "smi2: " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);
    TEST_ASSERT(smi1 == "CC1C=NNC=1");
    delete m;
  }

  {
    std::string pathName = getenv("RDBASE");
    pathName += "/External/AvalonTools/test_data/";
    std::ifstream ins((pathName + "heterocycle2.mol").c_str());
    std::string mb((std::istreambuf_iterator<char>(ins)),
                   std::istreambuf_iterator<char>());
    RWMol *m = MolBlockToMol(mb, false);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 11);
    m->updatePropertyCache();
    MolOps::cleanUp(*m);
    MolOps::setAromaticity(*m);

    std::string smi1 = AvalonTools::getCanonSmiles(mb, false);
    std::string smi2 = AvalonTools::getCanonSmiles(*m);
    std::cerr << "smi1: " << smi1 << std::endl;
    std::cerr << "smi2: " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);
    TEST_ASSERT(smi1 == "CN2C=CC1=CC(=O)NC=C12");
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testCountFps() {
  BOOST_LOG(rdInfoLog) << "testing substructure fingerprints " << std::endl;
  {
    SparseIntVect<boost::uint32_t> cv1(5000), cv2(5000);
    AvalonTools::getAvalonCountFP("c1ccccc1", true, cv1, 5000);
    AvalonTools::getAvalonCountFP("c1ccccc1.c1ccccc1", true, cv2, 5000);
    for (unsigned int i = 0; i < cv1.size(); ++i) {
      if (cv1[i] && (cv2[i] != 2 * cv1[i])) {
        std::cerr << "  mismatch: " << i << " " << cv1[i] << " " << cv2[i]
                  << std::endl;
      }
    }
    for (unsigned int i = 0; i < cv1.size(); ++i) {
      TEST_ASSERT(!cv1[i] || (cv2[i] == 2 * cv1[i]));
    }
  }
  {
    ROMol *m1 = static_cast<ROMol *>(SmilesToMol("c1ccccc1"));
    TEST_ASSERT(m1);
    ROMol *m2 = static_cast<ROMol *>(SmilesToMol("c1ccccc1.c1ccccc1"));
    TEST_ASSERT(m2);

    SparseIntVect<boost::uint32_t> cv1(5000), cv2(5000);
    AvalonTools::getAvalonCountFP(*m1, cv1, 5000);
    AvalonTools::getAvalonCountFP(*m2, cv2, 5000);
    for (unsigned int i = 0; i < cv1.size(); ++i) {
      if (cv1[i] && (cv2[i] != 2 * cv1[i])) {
        std::cerr << "  mismatch: " << i << " " << cv1[i] << " " << cv2[i]
                  << std::endl;
      }
    }
    for (unsigned int i = 0; i < cv1.size(); ++i) {
      TEST_ASSERT(!cv1[i] || (cv2[i] == 2 * cv1[i]));
    }
    delete m1;
    delete m2;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testInitStruChk() {
  BOOST_LOG(rdInfoLog) << "testing init struchk " << std::endl;
  {
    std::string pathName = getenv("RDBASE");
    pathName += "/Data/struchk/";
    std::string struchk_init =
        "-tm\n"
        "-ta " +
        pathName + std::string("checkfgs.trn\n") +
        "-tm\n"
        "-or\n"
        "-ca " +
        pathName + std::string("checkfgs.chk\n") +
        "-cc\n"
        "-cl 3\n"
        "-cs\n"
        "-cn 999\n"
        "-l " +
#ifdef _WIN32
        std::getenv("TEMP") +
#endif
        std::string(std::tmpnam(nullptr)) + std::string("\n");
    int errs = AvalonTools::initCheckMol(struchk_init);
    TEST_ASSERT(!errs);
    RDKit::ROMOL_SPTR m = AvalonTools::checkMol(errs, "c1ccccc1", true);
    AvalonTools::closeCheckMolFiles();
    TEST_ASSERT(errs == 0);
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub4075() {
  BOOST_LOG(rdInfoLog) << "testing Github #4075: AvalonTools.Generate2DCoords "
                          "results in an assert violation"
                       << std::endl;
  {
    auto mol = "Cl[C@H](C)CCC(F)Cl"_smiles;
    TEST_ASSERT(mol);
    AvalonTools::set2DCoords(*mol);
    TEST_ASSERT(mol->getNumConformers() == 1);
    std::unique_ptr<ROMol> newMol(MolBlockToMol(MolToMolBlock(*mol)));
    TEST_ASSERT(newMol);
    TEST_ASSERT(MolToSmiles(*newMol) == MolToSmiles(*mol));
  }
  {
    auto mol = "Cl[C@H](C)C([C@H](O)I)C[C@H](F)Cl"_smiles;
    TEST_ASSERT(mol);
    AvalonTools::set2DCoords(*mol);
    TEST_ASSERT(mol->getNumConformers() == 1);
    std::unique_ptr<ROMol> newMol(MolBlockToMol(MolToMolBlock(*mol)));
    TEST_ASSERT(newMol);
    TEST_ASSERT(MolToSmiles(*newMol) == MolToSmiles(*mol));
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub4330() {
  BOOST_LOG(rdInfoLog) << "testing Github #4330: AvalonTools::set2DCoords "
                          "failing for molecules which have Hs"
                       << std::endl;
  {
    auto mol = "F[C@H](OOF)Br"_smiles;
    TEST_ASSERT(mol);
    AvalonTools::set2DCoords(*mol);
    TEST_ASSERT(mol->getNumConformers() == 1);
    MolOps::addHs(*mol);
    AvalonTools::set2DCoords(*mol);
    TEST_ASSERT(mol->getNumConformers() == 1);
  }
  {  // example with the explicit Hs already in the graph
    SmilesParserParams params;
    params.removeHs = false;
    std::unique_ptr<RWMol> mol(SmilesToMol("F[C@]([H])(OOF)Br", params));
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 7);
    AvalonTools::set2DCoords(*mol);
    TEST_ASSERT(mol->getNumConformers() == 1);
  }
  {
    auto mol = "CC"_smiles;
    TEST_ASSERT(mol);
    AvalonTools::set2DCoords(*mol);
    TEST_ASSERT(mol->getNumConformers() == 1);
    MolOps::addHs(*mol);
    AvalonTools::set2DCoords(*mol);
    TEST_ASSERT(mol->getNumConformers() == 1);
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testNoAtomCTAB() {
  BOOST_LOG(rdInfoLog) << "testing CTAB with no atoms"
                       << std::endl;
  std::string data = "NullMol\r\n  ChemDraw04291110502D\r\n\r\n  0  0  0  0  0  0  0  0  0  0999 V2000\r\nM  END";
  bool isSmiles = false;
  int flags = 0x01 | 0x02 | 0x20;
  auto res = AvalonTools::getCanonSmiles(data, isSmiles, flags);
  TEST_ASSERT(res.empty());
}

int main() {
  RDLog::InitLogs();
#if 1
  test1();
  test2();
  test3();
  testGitHub1062();
  testRDK151();
  testSmilesFailures();
  testSubstructFps();
  testStruChk();
  testBadMolfile();
  testSmilesSegFault();
  testGithub336();
  testCountFps();
  testInitStruChk();
#endif
  testGithub4075();
  testGithub4330();
  testNoAtomCTAB();

  return 0;
}
