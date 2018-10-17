//
//  Copyright (c) 2014, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#ifdef _WIN32
#include <RDGeneral/test.h>
#include <Windows.h>
#else
#include <unistd.h>
#include <fcntl.h>
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <string>
#include <iostream>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>
#include "../RDKitBase.h"
#include "../FileParsers/FileParsers.h"  //MOL single molecule !
#include "../FileParsers/MolSupplier.h"  //SDF
#include "../SmilesParse/SmilesParse.h"
#include "../SmilesParse/SmilesWrite.h"
#include "../SmilesParse/SmartsWrite.h"
#include "FMCS.h"
#include "DebugTrace.h"  //#ifdef VERBOSE_STATISTICS_ON

#include "../Substruct/SubstructMatch.h"

using namespace RDKit;

unsigned long long T0;
unsigned long long t0;

void printTime() {
  unsigned long long t1 = nanoClock();
  double sec = double(t1 - t0) / 1000000.;
  printf("Time elapsed %.3lf seconds\n", sec);
  t0 = nanoClock();
}

std::string getSmilesOnly(
    const char* smiles,
    std::string* id = nullptr) {  // remove label, because RDKit parse FAILED
  const char* sp = strchr(smiles, ' ');
  unsigned n = (sp ? sp - smiles + 1 : strlen(smiles));
  if (id) *id = std::string(smiles + n);
  return std::string(smiles, n);
}

// UNIT Test Set:
//=========================================================================

void test1Basics() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "FMCS test1Basics()" << std::endl;

  std::vector<ROMOL_SPTR> mols;
  const char* smi[] = {
      "CC1CCC(N)CC1", "CC1CC(C)CC(C)C1",  // OK test.sdf
  };

  for (auto& i : smi) {
    std::string id;
    mols.push_back(ROMOL_SPTR(SmilesToMol(getSmilesOnly(i, &id))));
  }
  t0 = nanoClock();
  MCSResult res = findMCS(mols);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 7 && res.NumBonds == 7);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test32() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS test32" << std::endl;
  std::cout << "\n test32()\n";
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] = {
      //#  Using CHEMBL1515359 CHEMBL1590658 CHEMBL1447567 CHEMBL1384017
      // CHEMBL1456416 CHEMBL1308819 CHEMBL1703007 CHEMBL1707819 CHEMBL1500793
      // CHEMBL1334715
      // 32 . 1 31 33 0.82
      // S(-N1-C-C-O-C-C-1)(-c1:c:c:c(-N(-C-C)-C-C):c(-N-C(-C=C-c2:c:c:c:c:c:2)=O):c:1)(=O)=O
      "O=C(Nc1cc(S(N2CCOCC2)(=O)=O)ccc1N1CCOCC1)C=Cc1ccc(Cl)cc1  CHEMBL1515359",
      "c1ccc(C=CC(Nc2cc(S(N3CCOCC3)(=O)=O)ccc2N2CCOCC2)=O)cc1  CHEMBL1590658",
      "Cc1ccc(C=CC(=O)Nc2cc(S(N3CCOCC3)(=O)=O)ccc2N2CCOCC2)cc1  CHEMBL1447567",
      "c1ccc(C=CC(Nc2cc(S(N3CCOCC3)(=O)=O)ccc2N2CCCC2)=O)cc1  CHEMBL1384017",
      "O=C(C=Cc1ccc(F)cc1)Nc1cc(S(N2CCOCC2)(=O)=O)ccc1N1CCCC1  CHEMBL1456416",
      "c1cc(F)cc(C=CC(=O)Nc2c(N3CCCC3)ccc(S(N3CCOCC3)(=O)=O)c2)c1  "
      "CHEMBL1308819",
      "CCN1CCN(c2ccc(S(N3CCOCC3)(=O)=O)cc2NC(=O)C=Cc2ccc(C)cc2)CC1  "
      "CHEMBL1703007",
      "c1cc(C=CC(=O)Nc2cc(S(N3CCOCC3)(=O)=O)ccc2N2CCOCC2)c([N+]([O-])=O)cc1  "
      "CHEMBL1707819",
      "N#CC(=Cc1ccccc1)C(=O)Nc1cc(S(N2CCOCC2)(=O)=O)ccc1N1CCCC1  CHEMBL1500793",
      "C(=Cc1ccc2c(c1)OCO2)C(Nc1cc(S(=O)(=O)N2CCOCC2)ccc1N1CCOCC1)=O  "
      "CHEMBL1334715",
      // 31 33 0.35 sec MCS: CCN(CC)c1ccc(cc1NC(=O)C=Cc1ccccc1)S(=O)(=O)N1CCOCC1
  };
  for (auto& i : smi) mols.push_back(ROMOL_SPTR(SmilesToMol(getSmilesOnly(i))));
  t0 = nanoClock();
  MCSResult res = findMCS(mols);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 31 && res.NumBonds == 33);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test190() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS test190" << std::endl;
  std::cout << "\n test190()\n";
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] = {
      // # 190
      "COc1cc2nc(-c3cc(NC(=O)CSc4ccc(Cl)cc4)ccc3)oc2cc1  CHEMBL1479679",
      "COc1cc2nc(-c3cc(NC(=O)CSc4ccc(Cl)cc4)c(C)cc3)oc2cc1  CHEMBL1333382",
      "Cc1cc2oc(-c3cc(NC(=O)CSc4ccc(Cl)cc4)ccc3)nc2cc1  CHEMBL1437584",
      "COc1c(NC(=O)CSc2ccc(Cl)cc2)cc(-c2nc3ccccc3o2)cc1  CHEMBL1601350",
      "Cc1cc2nc(-c3cccc(NC(=O)CSc4ccc(Cl)cc4)c3)oc2cc1C  CHEMBL1398008",
      "Cc1cc2oc(-c3cc(NC(=O)CSc4ccc(Cl)cc4)c(C)cc3)nc2cc1  CHEMBL1612903",
      "COc1cc2nc(-c3cc(NC(=O)Cc4ccc(Cl)cc4)c(C)cc3)oc2cc1  CHEMBL1316483",
      "Cc1c(NC(=O)CSc2ccc(Cl)cc2)cccc1-c1nc2cc(Cl)ccc2o1  CHEMBL1568754",
      "COc1ccc2oc(-c3ccc(C)c(NC(=O)COc4cc(C)cc(C)c4)c3)nc2c1  CHEMBL1436972",
      "Cc1ccc(SCC(=O)Nc2cc(-c3nc4cc(C)ccc4o3)c(O)cc2)cc1  CHEMBL1611932",
      //# 19 21 1.37 sec MCS: CC(=O)Nc1cccc(c1)-c1nc2ccccc2o1
      //  19 21 2.36 sec MCS: CC(=O)Nc1cccc(c1)-c1nc2ccccc2o1 19 atoms, 21 bonds
  };
  for (auto& i : smi) mols.push_back(ROMOL_SPTR(SmilesToMol(getSmilesOnly(i))));
  t0 = nanoClock();
  MCSResult res = findMCS(mols);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 19 && res.NumBonds == 21);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test45() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS test45" << std::endl;
  std::cout << "\n test45()\n";
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] = {
      // SLOW
      // test 45 #  Using CHEMBL551656 CHEMBL563796 CHEMBL561978 CHEMBL559467
      // CHEMBL550503 CHEMBL562866 CHEMBL552190 CHEMBL181547 CHEMBL359567
      // CHEMBL373316
      // 45 . 1 30 32 27.01
      // n12-C-c:c(-c:2:c:c2-C(-O)(-C-C)-C(-O-C-c:2:c:1=O)=O):n:c(:c:c:c-O):c(:c):c-C-C-C
      "CCC1(O)c2cc3n(c(=O)c2COC1=O)Cc1c-3nc2ccc(OC)cc2c1C1CCCCC1 CHEMBL551656",
      "CCC1(O)c2cc3n(c(=O)c2COC1=O)Cc1c-3nc2ccc(OC)cc2c1C1CCCC1 CHEMBL563796",  // Q
      "CCC1(O)C(=O)OCc2c1cc1n(c2=O)Cc2c-1nc1ccc(OC)cc1c2C1CCCCCC1 CHEMBL561978",
      "CCC1(O)C(=O)OCc2c1cc1n(c2=O)Cc2c-1nc1ccc(OC)cc1c2C1CCC1 CHEMBL559467",
      "CCC1(O)C(=O)OCc2c1cc1n(c2=O)Cc2c-1nc1ccc(O)cc1c2C1CCCC1 CHEMBL550503",
      "CCC1(O)c2cc3n(c(=O)c2COC1=O)Cc1c-3nc2ccc(O)cc2c1C1CCCCCC1 CHEMBL562866",
      "CCC1(O)C(=O)OCc2c1cc1n(c2=O)Cc2c-1nc1ccc(O)cc1c2C1CCCCC1 CHEMBL552190",
      "CCC1(O)c2c(c(=O)n3c(c2)-c2nc4cc5c(cc4c(C4CCCC4)c2C3)OCO5)COC1=O "
      "CHEMBL181547",
      "CCC1(O)c2c(c(=O)n3c(c2)-c2nc4c(c(C5CCCCC5)c2C3)cc2c(c4)OCO2)COC1=O "
      "CHEMBL359567",
      "CCCc1c(OC)ccc2nc3c(c(CC)c21)Cn1c-3cc2c(c1=O)COC(=O)C2(O)CC CHEMBL373316",
  };
  for (auto& i : smi) mols.push_back(ROMOL_SPTR(SmilesToMol(getSmilesOnly(i))));
  t0 = nanoClock();
  MCSResult res = findMCS(mols);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 31 && res.NumBonds == 33);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test3() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS test3" << std::endl;
  std::cout << "\n test3()\n";
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] = {
      // TEST 3
      "CN(C)c1ccc(CC(=O)NCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL153934",
      "CN(C)c1ccc(CC(=O)NCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL152361",
      "CN(C)c1ccc(CC(=O)NCCCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL157336",
      "CN(C)c1ccc(CC(=O)NCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL157429",
      "CN(C)c1ccc(CC(=O)NCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL357551",
      "CN(C)c1ccc(CC(=O)NCCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL421974",
      "CN(C)c1ccc(CC(NCCCCCC(NO)=O)=O)cc1 CHEMBL484488",
      "CC(C)Cc1ccc(C(C)C(=O)NC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL564780",
      "c1cc([N+]([O-])=O)ccc1CC(=O)NC1CCCCCC1 CHEMBL1553142",
      "CC1(C)NC(C)(C)CC(NC(=O)Cc2ccccc2)C1 CHEMBL1703640",
      //# 3 . 1 14 14 0.08 sec MCS: CCCCNC(=O)Cc1ccccc1
  };
  for (auto& i : smi) mols.push_back(ROMOL_SPTR(SmilesToMol(getSmilesOnly(i))));
  t0 = nanoClock();
  MCSResult res = findMCS(mols);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 14 && res.NumBonds == 14);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testRing1() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testRing1" << std::endl;
  std::cout << "\ntestRing1()\n";
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] = {
      "COCc1c(ncc2[nH]c3cccc(Oc4ccc(Cl)cc4)c3c12)C(=O)OC(C)C",
      //      "COCc1cnc(C(=O)OC(C)C)c2[nH]c3cc(Oc4ccc(Cl)cc4)ccc3c12", //
      //      original molecule
      "COCc1cnc(C(=O)OC(C)C)c2[nH]ccc(Oc4ccc(Cl)cc4)cccc12",  // ring 3 removed
  };
  for (auto& i : smi)
    mols.push_back(
        ROMOL_SPTR(SmilesToMol(getSmilesOnly(i))));  // with RING INFO

  MCSParameters p;
  p.BondCompareParameters.RingMatchesRingOnly = true;
  p.BondCompareParameters.CompleteRingsOnly = true;
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 12 && res.NumBonds == 12);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test504() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS test504" << std::endl;
  std::cout << "\ntest504()\n";
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] = {
      // TEST 504
      "C(CCNC(C1CC1[c:1]1[c:2]c(Cl)c(Cl)c[c:3]1)=O)CCN1CCC(NC(Nc2ccc(Cl)cc2)=O)"
      "C1 CHEMBL545864",  // - QUERY

      "FC(F)(F)c1cc(NC(N2CCCN(CCCCCNC(C3CC3c3ccc(Cl)c(Cl)c3)=O)CC2)=O)ccc1Cl "
      "CHEMBL528228",
      "FC(F)(F)c1cc(NC(NC2CCN(CCCCCNC(C3CC3c3ccc(Cl)c(Cl)c3)=O)C2)=O)ccc1Cl "
      "CHEMBL525875",
      "Fc1ccc(NC(N2CCCN(CCCCCNC(C3CC3c3ccc(Cl)c(Cl)c3)=O)CC2)=O)cc1C(F)(F)F "
      "CHEMBL527277",
      "FC(F)(F)c1cc(NC(NC2CCN(CCCCCNC(C3CC3c3ccc(Cl)c(Cl)c3)=O)CC2)=O)ccc1Cl "
      "CHEMBL537333",
      "Fc1ccc(NC(NC2CCN(CCCCCNC(C3CC3c3ccc(Cl)c(Cl)c3)=O)C2)=O)cc1C(F)(F)F "
      "CHEMBL588077",
      "FC(F)(F)c1ccc(NC(NC2CCN(CCCCCNC(C3CC3c3cc(Cl)c(Cl)cc3)=O)C2)=O)cc1 "
      "CHEMBL525307",
      "Fc1ccc(NC(NC2CCN(CCCCCNC(C3CC3c3ccc(Cl)c(Cl)c3)=O)CC2)=O)cc1C(F)(F)F "
      "CHEMBL581847",
      "FC(F)(F)c1ccc(NC(NC2CCN(CCCCCNC(C3CC3c3cc(Cl)c(Cl)cc3)=O)CC2)=O)cc1 "
      "CHEMBL579547",
      "N#Cc1cccc(NC(NC2CCN(CCCCCNC(C3CC3c3ccc(Cl)c(Cl)c3)=O)CC2)=O)c1 "
      "CHEMBL529994",
  };
  RWMol* qm = SmilesToMol(getSmilesOnly(smi[0]));
  unsigned nq = qm->getNumAtoms();
  for (size_t ai = 0; ai < nq; ai++) {
    Atom* atom = qm->getAtomWithIdx(ai);
    atom->setProp("molAtomMapNumber", (int)ai);
  }
  std::cout << "Query +MAP " << MolToSmiles(*qm) << "\n";
  mols.push_back(ROMOL_SPTR(qm));  // with RING INFO
  for (size_t i = 1; i < sizeof(smi) / sizeof(smi[0]); i++)
    mols.push_back(
        ROMOL_SPTR(SmilesToMol(getSmilesOnly(smi[i]))));  // with RING INFO
  t0 = nanoClock();
  MCSResult res = findMCS(mols);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 34 && res.NumBonds == 36);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test18() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS test18" << std::endl;
  std::cout << "\ntest18()\n";
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] = {
      // TEST 18
      "Cc1nc(CN(C(C)c2ncccc2)CCCCN)ccc1 CHEMBL1682991",  //-- QUERY
      "Cc1ccc(CN(C(C)c2ccccn2)CCCCN)nc1 CHEMBL1682990",
      "Cc1cccnc1CN(C(C)c1ccccn1)CCCCN CHEMBL1682998",
      "CC(N(CCCCN)Cc1c(N)cccn1)c1ccccn1 CHEMBL1682987",
      "Cc1cc(C)c(CN(C(C)c2ccccn2)CCCCN)nc1 CHEMBL1682992",
      "Cc1cc(C(C)N(CCCCN)Cc2c(C)cccn2)ncc1 CHEMBL1682993",
      "Cc1nc(C(C)N(CCCCN)Cc2nc3c([nH]2)cccc3)ccc1 CHEMBL1682878",
      "CC(c1ncccc1)N(CCCCN)Cc1nc2c([nH]1)cccc2 CHEMBL1682867",
      "CC(N(CCCCN)Cc1c(C(C)(C)C)cccn1)c1ccccn1 CHEMBL1682989",
      "CC(N(CCCCN)Cc1c(C(F)(F)F)cccn1)c1ccccn1 CHEMBL1682988",
      //# 18 .  20 20 0.04 sec. Python MCS: CC(c1ccccn1)N(CCCCN)Ccnccc
  };
  RWMol* qm = SmilesToMol(getSmilesOnly(smi[0]));
  unsigned nq = qm->getNumAtoms();
  for (size_t ai = 0; ai < nq; ai++) {
    Atom* atom = qm->getAtomWithIdx(ai);
    atom->setProp("molAtomMapNumber", (int)ai);
  }
  std::cout << "Query +MAP " << MolToSmiles(*qm) << "\n";
  mols.push_back(ROMOL_SPTR(qm));  // with RING INFO
  for (size_t i = 1; i < sizeof(smi) / sizeof(smi[0]); i++)
    mols.push_back(
        ROMOL_SPTR(SmilesToMol(getSmilesOnly(smi[i]))));  // with RING INFO
  t0 = nanoClock();
  MCSResult res = findMCS(mols);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 21 && res.NumBonds == 21);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testThreshold() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testThreshold" << std::endl;
  std::cout << "\ntestThreshold()\n";
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] = {
      "CCC", "CCCO", "CCCN", "CC",
      //        "CCC", "CC", //th=0.5
  };
  for (auto& i : smi) mols.push_back(ROMOL_SPTR(SmilesToMol(getSmilesOnly(i))));
  findMCS(mols);
  MCSParameters p;
  p.Threshold = 0.7;
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  p.Threshold = 1.0;  // restore default value
  TEST_ASSERT(res.NumAtoms == 3 && res.NumBonds == 2);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test330() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS test330" << std::endl;
  std::cout << "\ntest330()\n";
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] = {
      // TEST 330  40 sec
      "CCC(C)C(NC(=O)C(NC(C(CCC(O)=O)NC(=O)C(NC(=O)C(NC(C(CC(O)=O)NC(C(CC(C)C)"
      "NC(C(Cc1ccccc1)NC(CN)=O)=O)=O)=O)C(C)CC)C(C)CC)=O)CCCCN)C(NC(C)C(NC("
      "CCCCN)C(NC(CO)C(NC(Cc1c[nH]c2c1cccc2)C(O)=O)=O)=O)=O)=O CHEMBL1240742",
      "CCC(C)C(NC(=O)C(NC(C(CCCCN)NC(=O)C(NC(=O)C(NC(C(CC(O)=O)NC(C(Cc1ccccc1)"
      "NC(C(CC(C)C)NC(CN)=O)=O)=O)=O)C(C)CC)C(C)CC)=O)CCCCN)C(NC(C)C(NC(CCC(O)="
      "O)C(NC(CO)C(NC(Cc1c[nH]c2c1cccc2)C(O)=O)=O)=O)=O)=O CHEMBL1240736",
      "CCC(C)C(NC(CN)=O)C(NC(C(NC(CC(O)=O)C(NC(C(NC(C)C(NC(CCCCN)C(NC(C(NC(CC("
      "C)C)C(NC(Cc1ccccc1)C(NC(CCC(O)=O)C(NC(CO)C(NC(Cc1c[nH]c2c1cccc2)C(O)=O)="
      "O)=O)=O)=O)=O)CCCCN)=O)=O)=O)C(C)CC)=O)=O)C(C)CC)=O CHEMBL1240738",
      "CCC(C)C(NC(CN)=O)C(NC(Cc1ccccc1)C(NC(CC(O)=O)C(NC(CCCCN)C(NC(CC(C)C)C("
      "NC(C)C(NC(CCCCN)C(NC(CCC(O)=O)C(NC(C(NC(CO)C(NC(C(NC(Cc1c[nH]c2c1cccc2)"
      "C(O)=O)=O)C(C)CC)=O)=O)C(C)CC)=O)=O)=O)=O)=O)=O)=O)=O CHEMBL1240740",
      "CCC(C)C(NC(CN)=O)C(NC(Cc1c[nH]c2c1cccc2)C(NC(CO)C(NC(CC(O)=O)C(NC(CC(C)"
      "C)C(NC(C)C(NC(CCC(O)=O)C(NC(C(NC(C(NC(CCCCN)C(NC(CCCCN)C(NC(Cc1ccccc1)C("
      "O)=O)=O)=O)=O)C(C)CC)=O)C(C)CC)=O)=O)=O)=O)=O)=O)=O CHEMBL1240741",
      "CCC(C)C(NC(=O)C(NC(=O)C(CCCCN)NC(C(CC(C)C)NC(C(Cc1c[nH]c2c1cccc2)NC(CN)="
      "O)=O)=O)CCCCN)C(NC(CCC(O)=O)C(NC(CO)C(=O)NC(C(NC(C(NC(CC(O)=O)C(NC(C)C("
      "NC(Cc1ccccc1)C(O)=O)=O)=O)=O)C(C)CC)=O)C(C)CC)=O)=O CHEMBL1240743",
      "CCC(C)C(NC(C(CCC(O)=O)NC(C(CC(O)=O)NC(C(CC(C)C)NC(C(Cc1ccccc1)NC(C)=O)="
      "O)=O)=O)=O)C(NC(Cc1c[nH]c2ccccc12)C(O)=O)=O CHEMBL431874",
      "CCC(C)C(NC(C(CC(O)=O)NC(C(CC(C)C)NC(C(Cc1ccccc1)NC(C)=O)=O)=O)=O)C(NC("
      "CCC(O)=O)C(NC(Cc1c[nH]c2ccccc12)C(O)=O)=O)=O CHEMBL262166",
      "CCC(C)C(NC(C(CC(O)=O)NC(C(CC(C)C)NC(C(Cc1ccccc1)NC(C)=O)=O)=O)=O)C(NC("
      "CCCCN)C(NC(Cc1c[nH]c2c1cccc2)C(O)=O)=O)=O CHEMBL313122",
      "CCC(C)C(NC(C(CCCCN)NC(C(CC(O)=O)NC(C(CC(C)C)NC(C(Cc1ccccc1)NC(C)=O)=O)="
      "O)=O)=O)C(NC(Cc1c[nH]c2c1cccc2)C(O)=O)=O CHEMBL314239",
      //# 330 F  42 41 30.93 sec MCS:
      //[#6]-[#6](-[#7]-[#6](-[#6](-[#6])-[#7]-[#6](-[#6](-[#6])-[#7]-[#6](-[#6](-[#6]-[#6]-[#6])-[#7]-[#6](-[#6](-[#6])-[#7]-[#6](-[#6])=[#8])=[#8])=[#8])=[#8])=[#8])-[#6](-[#7]-[#6](-[#6]-[#6](:[#6]):[#6]:[#6]:[#6]:[#6])-[#6](-[#8])=[#8])=[#8]
  };
  for (auto& i : smi) mols.push_back(ROMOL_SPTR(SmilesToMol(getSmilesOnly(i))));
  MCSParameters p;
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 42 && res.NumBonds == 41);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testTarget_no_10188_30149() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testTarget_no_10188_30149" << std::endl;
  std::cout << "\ntestTarget_no_10188_30149()\n";
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] = {
      // Target_no_10188_30149.txt // VERY SLOWER than Python
      "CN(C)CCNC(=O)c1ccc(-c2n[nH]c3cc(Nc4ccccc4Cl)ccc32)cc1 CHEMBL399167",
      "O=C(O)c1cccc(-c2[nH]nc3cc(Nc4ccccc4Cl)ccc32)c1 CHEMBL197613",
      "c1ccc(Nc2ccc3c(c2)[nH]nc3-c2ccccc2)cc1 CHEMBL383177",  /// == QUERY
      "NC(=O)c1cccc(-c2[nH]nc3cc(Nc4ccccc4Cl)ccc32)c1 CHEMBL199136",
      "Clc1ccccc1Nc1ccc2c(c1)n[nH]c2-c1ccccc1 CHEMBL440566",
      "O=C(NCCCN1CCOCC1)c1cccc(-c2[nH]nc3cc(Nc4ccccc4Cl)ccc32)c1 CHEMBL198687",
      "O=C(O)c1ccc(-c2[nH]nc3cc(Nc4ccccc4Cl)ccc32)cc1 CHEMBL197698",
      "O=C(NC1CCNCC1)c1cccc(-c2n[nH]c3cc(Nc4ccccc4Cl)ccc32)c1 CHEMBL194806",
      "COc1ccccc1Nc1ccc2c(c1)[nH]nc2-c1ccccc1 CHEMBL254443",
      "CN(C)CCNC(=O)c1cccc(-c2[nH]nc3cc(Nc4ccccc4Cl)ccc32)c1 CHEMBL198821",
  };
  for (auto& i : smi) mols.push_back(ROMOL_SPTR(SmilesToMol(getSmilesOnly(i))));
  MCSParameters p;
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 15 && res.NumBonds == 14);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testTarget_no_10188_49064() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testTarget_no_10188_49064" << std::endl;
  std::cout << "\ntestTarget_no_10188_49064()\n";
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] = {
      "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3cccc(O)c3)nc21",
      "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(F)cc3)nc21",
      "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(OC4OC(CO)C(O)C(O)C4O)cc3)nc21",
      "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(NC(=O)CCl)cc3)nc21",
      "Cn1c2nc(Nc3ccc(NC(=O)CCN)cc3)ncc2cc(-c2c(Cl)cccc2Cl)c1=O",
      "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3cccc(COCC(O)CO)c3)nc21",
      "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(O)cc3)nc21",
      "CC(=O)Nc1ccc(Nc2ncc3cc(-c4c(Cl)cccc4Cl)c(=O)n(C)c3n2)cc1",
      "Cn1c2nc(Nc3ccc(N)cc3)ncc2cc(-c2c(Cl)cccc2Cl)c1=O",
      "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(NC(=O)CCNC(=O)OC(C)(C)C)cc3)"
      "nc21",
      "Cc1ccc(Nc2ncc3cc(-c4c(Cl)cccc4Cl)c(=O)n(C)c3n2)cc1",
      "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3cccc(CO)c3)nc21",
      "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(NCC(O)CO)cc3)nc21",
      "CCc1cccc(Nc2ncc3cc(-c4c(Cl)cccc4Cl)c(=O)n(C)c3n2)c1",
      "Cn1c2nc(Nc3cccc(N)c3)ncc2cc(-c2c(Cl)cccc2Cl)c1=O",
      "CC(=O)Nc1cccc(Nc2ncc3cc(-c4c(Cl)cccc4Cl)c(=O)n(C)c3n2)c1",
      "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(CCO)cc3)nc21",
      "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(I)cc3)nc21",
      "CN1CCN(C(=O)c2ccc(Nc3ncc4cc(-c5c(Cl)cccc5Cl)c(=O)n(C)c4n3)cc2)CC1",
  };
  for (auto& i : smi) mols.push_back(ROMOL_SPTR(SmilesToMol(getSmilesOnly(i))));
  MCSParameters p;
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 15 && res.NumBonds == 14);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testSegFault() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testSegFault" << std::endl;
  std::cout << "\ntestSegFault()\n";
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] = {
      "CN(CCCN(C)CCc1ccccc1)CCOC(c1ccccc1)c1ccccc1",
      "CN(CCCc1ccccc1)CCCN(C)CCOC(c1ccccc1)c1ccccc1",
      "Fc1ccc(C(OCCNCCCNCCc2ccccc2)c2ccc(F)cc2)cc1",
      "O=C(Cc1ccccc1)NCCCCNCCOC(c1ccc(F)cc1)c1ccc(F)cc1",
      "O=C(Cc1ccccc1)NCCNCCOC(c1ccc(F)cc1)c1ccc(F)cc1",
      "O=C(Cc1ccc(Br)cc1)NC=CNCCOC(c1ccc(F)cc1)c1ccc(F)cc1",
      "O=C(Cc1ccc(F)cc1)NCCNCCOC(c1ccc(F)cc1)c1ccc(F)cc1",
      "O=C(Cc1ccccc1)NCCCNCCOC(c1ccc(F)cc1)c1ccc(F)cc1",
      "CN(CCOC(c1ccc(F)cc1)c1ccc(F)cc1)CCN(C)CCOC(c1ccc(F)cc1)c1ccc(F)cc1",
      "COC(=O)C1C2CCC(CC1C(=O)Oc1ccccc1)N2C",
      "O=C1CN(CCc2ccccc2)CCN1CCOC(c1ccc(F)cc1)c1ccc(F)cc1",
      "CN(CCOC(c1ccccc1)c1ccccc1)CCN(C)CCc1ccc(F)cc1",
  };
  for (auto& i : smi) mols.push_back(ROMOL_SPTR(SmilesToMol(getSmilesOnly(i))));
  MCSParameters p;
  p.BondCompareParameters.RingMatchesRingOnly = true;
  p.BondCompareParameters.CompleteRingsOnly = true;
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 6 && res.NumBonds == 6);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testAtomCompareIsotopes() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testAtomCompareIsotopes" << std::endl;
  std::cout << "\ntestAtomCompareIsotopes()\n";
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] = {
      "CC[13NH2]",
      "CC[13CH3]",
  };
  for (auto& i : smi) mols.push_back(ROMOL_SPTR(SmilesToMol(getSmilesOnly(i))));
  MCSParameters p;
  p.AtomTyper = MCSAtomCompareIsotopes;
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 3 && res.NumBonds == 2);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testAtomCompareAnyAtom() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testAtomCompareAnyAtom" << std::endl;
  std::cout << "\ntestAtomCompareAnyAtom()\n";
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] = {
      "c1ccccc1C", "c1ccccc1O", "c1ccccc1Cl",
      "c1ccccc1F",  // opt
      "c1ccccc1N",  // opt
  };
  for (auto& i : smi) mols.push_back(ROMOL_SPTR(SmilesToMol(getSmilesOnly(i))));
  MCSParameters p;
  p.AtomTyper = MCSAtomCompareAny;
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 7 && res.NumBonds == 7);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testAtomCompareAnyAtomBond() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testAtomCompareAnyAtomBond"
                       << std::endl;
  std::cout << "\ntestAtomCompareAnyAtom()\n";
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] = {
      "C1CCCCC1=C", "c1ccccc1O", "c1ccccc1Cl",
      "c1ccccc1F",  // opt
      "c1ccccc1N",  // opt
  };
  for (auto& i : smi) mols.push_back(ROMOL_SPTR(SmilesToMol(getSmilesOnly(i))));
  t0 = nanoClock();
  MCSParameters p;
  p.AtomTyper = MCSAtomCompareAny;
  p.BondTyper = MCSBondCompareAny;
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 7 && res.NumBonds == 7);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testSimple() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testSimple" << std::endl;
  std::cout << "\ntestSimple()\n";
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] = {
      // LONG TIME TEST for performance analisis (about 30 sec)
      "CC(C)CC(NC(=O)C(Cc1ccc(NC(C)=O)cc1)NC(=O)C(Cc1ccc(NC(C)=O)cc1)NC(C("
      "CO)NC(C(NC(c1ccncc1)=O)NC(=O)C(Cc1ccc(Cl)cc1)NC=O)=O)=O)C(NC("
      "CCCCNC(C)C)C(N1C(C(=O)NC(C)C(N)=O)CCC1)=O)=O CHEMBL439258 modified "
      "QUERY",  // CHEMBL439258
      "CC(C)CC(NC(=O)C(Cc1ccc(NC(C)=O)cc1)NC(=O)C(Cc1ccccc1)NC(C(CO)NC(C("
      "NC(c1ccncc1)=O)NC(=O)C(Cc1ccc(Cl)cc1)NC(C(NC(C)=O)Cc1cc2ccccc2cc1)="
      "O)=O)=O)C(NC(CCCCNC(C)C)C(N1C(C(=O)NC(C)C(N)=O)CCC1)=O)=O "
      "CHEMBL439258",  // CHEMBL439258
      "CC(C)CC(NC(=O)CNC(=O)C(Cc1ccc(NC(C)=O)cc1)NC(C(CO)NC(C(NC(c1ccncc1)"
      "=O)NC(=O)C(Cc1ccc(Cl)cc1)NC(C(NC(C)=O)Cc1cc2ccccc2cc1)=O)=O)=O)C("
      "NC(CCCCNC(C)C)C(N1C(C(=O)NC(C)C(N)=O)CCC1)=O)=O CHEMBL439258 "
      "modified",  // CHEMBL439258

      "CCCCC(NC(C(CCC(O)=O)NC(C(CC(C)C)NC(C(C(C)C)NC(=O)C(CCC(O)=O)NC(C("
      "CCCN=C(N)N)NC(C(NC(=O)C(NC(C(NC(C1CCCNC(=O)CCC(N)C(=O)NC(CC(C)C)C(="
      "O)NC(C(C)O)C(=O)N1)=O)Cc1c[nH]cn1)=O)CC(C)C)CC(C)C)=O)=O)=O)=O)=O)"
      "C(NC(C)C(NC(CCCN=C(N)N)C(NC(C)C(NC(CCC(O)=O)C(NC(CCC(N)=O)C(NC(CC("
      "C)C)C(NC(C)C(NC(CCC(N)=O)C(NC(CCC(N)=O)C(NC(C)C(NC(Cc1c[nH]cn1)C("
      "NC(CO)C(NC(CC(N)=O)C(NC(CCCN=C(N)N)C(NC(CCCCN)C(NC(CC(C)C)C(NC("
      "CCCC)C(NC(C(NC(C(C)CC)C(NC(C(N)=O)C(C)CC)=O)=O)CCC(O)=O)=O)=O)=O)="
      "O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O CHEMBL438567",
      "CCC(C)C(NC(CNC(=O)C(C)NC(=O)C(C)NC(C(Cc1nc[nH]c1)NC(C(CC(N)=O)NC("
      "CNC(C(CO)NC(=O)C(C)NC(=O)C(CCC(N)=O)NC(C(NC(=O)C(NC(C(CCCN=C(N)N)"
      "NC(C(CCC(N)=O)NC(C(NC(C(CCCN=C(N)N)NC(CNC(C(CCC(N)=O)NC(C(CC(C)C)"
      "NC(C(C)N)=O)=O)=O)=O)=O)CC(C)C)=O)=O)=O)CC(C)C)CC(C)C)=O)=O)=O)=O)="
      "O)=O)C(NC(CC(C)C)C(NC(C(O)C)C(NC(CCSC)C(O)=O)=O)=O)=O CHEMBL429374",
      "CC(C)CC1NC(=O)C(CCCCN)NC(=O)C(Cc2ccc(O)cc2)NC(=O)CNC(=O)C2NC(=O)C("
      "NC(C(C(C)C)NC(CNC(C3NC(=O)CC3)=O)=O)=O)CSSCC(C(O)=O)NC(=O)C3N("
      "CCC3O)C(=O)C(Cc3ccccc3)NC(=O)C(CSSC2)NC1=O CHEMBL1076370",
  };
  for (auto& i : smi) mols.push_back(ROMOL_SPTR(SmilesToMol(getSmilesOnly(i))));
  MCSParameters p;
  p.BondCompareParameters.RingMatchesRingOnly = true;
  p.BondCompareParameters.CompleteRingsOnly = true;
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 15 && res.NumBonds == 14);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testSimpleFast() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testSimpleFast" << std::endl;
  std::cout << "\ntestSimpleFast()\n";
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] = {
      // SHORT TEST for 26 bonds.
      // Python MCS = 26 bonds : COCc1cncc(c1):n:c1cccc(Oc2ccc(Cl)cc2)c1
      // MCS 26: COCc1c-ncc(c1)nc1cccc(c1)Oc1ccc(Cl)cc1 24 atoms, 26 bonds
      ///            "COCC1=C(N=CC2=C1C1=C(OC3=CC=C(Cl)C=C3)C=CC=C1N2)C(=O)OC(C)C",
      ///            "COCC1=CN=C(C(=O)OC(C)C)C2=C1C1=CC=C(OC3=CC=C(Cl)C=C3)C=C1N2",
      // The SAME, but pre AROMATIZATED (else PRECONDITION Exception with
      // Implicit Hs / 16 bonds only)
      "COCc1c(ncc2[nH]c3cccc(Oc4ccc(Cl)cc4)c3c12)C(=O)OC(C)C",
      "COCc1cnc(C(=O)OC(C)C)c2[nH]c3cc(Oc4ccc(Cl)cc4)ccc3c12",
  };
  for (auto& i : smi) mols.push_back(ROMOL_SPTR(SmilesToMol(getSmilesOnly(i))));
  MCSParameters p;
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 24 && res.NumBonds == 26);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void compareChirality(const char* target, const char* query,
                      bool useChirality) {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  std::cout << "test target = " << target << " query = " << query
            << " useChirality = " << useChirality << std::endl;
  ROMOL_SPTR target_ptr(SmilesToMol(target));
  ROMOL_SPTR query_ptr(SmartsToMol(query));

  std::vector<std::pair<int, int>> vect;
  bool sub_res = SubstructMatch(*target_ptr.get(), *query_ptr.get(), vect, true,
                                useChirality);

  std::cout << "res sub = " << sub_res << " size=" << vect.size() << std::endl;

  std::vector<ROMOL_SPTR> mols;
  mols.push_back(target_ptr);
  mols.push_back(query_ptr);

  MCSParameters p;
  p.AtomCompareParameters.MatchChiralTag = useChirality;
  p.BondCompareParameters.MatchStereo = useChirality;

  MCSResult mcs_res = findMCS(mols, &p);
  std::cout << "MCS: " << mcs_res.SmartsString << " " << mcs_res.NumAtoms
            << " atoms, " << mcs_res.NumBonds << " bonds\n";

  bool mcs_resb = query_ptr->getNumAtoms() == mcs_res.NumAtoms &&
                  query_ptr->getNumBonds() == mcs_res.NumBonds;
  TEST_ASSERT(sub_res == mcs_resb);
  if (sub_res != mcs_resb) {  // || vect.size() != mcs_res.NumAtoms) {
    BOOST_LOG(rdInfoLog) << "mcs_resb=" << mcs_resb
                         << "\t*** TEST FAILED ***\n";  // exit(1);
  }
}

void testSubMcsChirality(const char* target, const char* query) {
  compareChirality(target, query, false);
  compareChirality(target, query, true);
}

void testChirality() {
  BOOST_LOG(rdInfoLog) << "\n-------------------------------------"
                       << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testChirality" << std::endl;

  testSubMcsChirality("O[C@H](F)CCl", "C(F)C");  // MCS = CCF
  testSubMcsChirality("CC[C@H](F)Cl", "CCC");    // MCS = CCC
  testSubMcsChirality("CCC(F)Cl", "CC");         // MCS = CC

  // GREG's tests:
  testSubMcsChirality("C[C@H](F)CCl", "C[C@H](F)C");
  testSubMcsChirality("C[C@H](F)CCl", "C[C@@H](F)C");
  testSubMcsChirality(
      "O[C@H](F)CCl",
      "CC(F)C");  // DISMATCH. but actually it has non chiral MCS = CCF
  testSubMcsChirality("O[C@H](F)CCl", "OC(F)C");
  testSubMcsChirality("O[C@H](F)CCl", "O[C@H](F)C");
  testSubMcsChirality("O[C@H](F)CCl", "O[C@@H](F)C");

  testSubMcsChirality("OC(F)CCl", "OC(F)C");
  testSubMcsChirality("OC(F)CCl", "O[C@H](F)C");
  testSubMcsChirality("OC(F)CCl", "O[C@@H](F)C");

  testSubMcsChirality(
      "O[C@H](F)CCl",
      "O[C@H]C");  // Is OCC = MCS ? Degree of [C@H] is different
  testSubMcsChirality(
      "O[C@H](F)CCl",
      "O[C@@H]C");  // Is OCC = MCS ? Degree of [C@H] is different

  // ADD-IN TESTS:
  std::cout
      << "\n0. <<<<<<<<< actual MCS: [#6]-[#6]-[#6] 3 atoms, 2 bonds >>>>>>>\n";
  testSubMcsChirality(
      "CC[C@H](F)Cl",
      "CC[C@H]");  //  actual MCS is CCC always. We lost last [C@H]

  /* //TEMP (asymmetric implementation of match algorithm, that is incorrect for
  FMCS task):
     std::cout << "\nTEMP query/targ==MCS[asymmetric implementation of match
  algorithm]\n";
     std::cout << "1. <<<<<<<<< MCS: [#8]-[#6]-[#9] 3 atoms, 2 bonds >>>>>>>\n";
     testSubMcsChirality ("O[C@H](F)", "OC(F)"    ); // PASSED
  / /   testSubMcsChirality ("OC(F)"    , "O[C@H](F)"); // FAILED !!!
     std::cout << "\n2.1<<<<<<<<< MCS: [#8]-[#6]-[#6] 3 atoms, 2 bonds
  >>>>>>>\n";
     testSubMcsChirality ("OC[C@H](F)C", "OC(F)C"    ); // PASSED but incorrect
     std::cout << "\n2.2<<<<<<<<< MCS: [#8]-[#6](-[#9])-[#6] 4 atoms, 3 bonds
  >>>>>>>\n";
     testSubMcsChirality ("OC(F)C"     , "O[C@H](F)C"); // FAILED !!!
  */
  std::cout << "\tdone" << std::endl;
}

void testJSONParameters() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testJSONParameters" << std::endl;
  MCSParameters pj;

  parseMCSParametersJSON(nullptr, &pj);
  parseMCSParametersJSON("{}", nullptr);

  pj = MCSParameters();  // parsing of empty string keeps default values
  parseMCSParametersJSON("", &pj);
  TEST_ASSERT(pj.MaximizeBonds == true && pj.Threshold == 1.0 &&
              pj.Timeout == (unsigned int)-1 &&
              pj.AtomCompareParameters.MatchValences == false &&
              pj.AtomCompareParameters.MatchChiralTag == false &&
              pj.BondCompareParameters.MatchStereo == false &&
              pj.BondCompareParameters.RingMatchesRingOnly == false &&
              pj.BondCompareParameters.CompleteRingsOnly == false);

  pj = MCSParameters();
  const char json[] =
      "{\"MaximizeBonds\": false, \"Threshold\": 0.7, \"Timeout\": 3"
      ", \"MatchValences\": true, \"MatchChiralTag\": true"
      ", \"MatchStereo\": true, \"RingMatchesRingOnly\": true, "
      "\"CompleteRingsOnly\": true"
      ", \"InitialSeed\": \"CNC\""
      "}";
  parseMCSParametersJSON(json, &pj);
  TEST_ASSERT(pj.MaximizeBonds == false && pj.Threshold == 0.7 &&
              pj.Timeout == 3 &&
              pj.AtomCompareParameters.MatchValences == true &&
              pj.AtomCompareParameters.MatchChiralTag == true &&
              pj.BondCompareParameters.MatchStereo == true &&
              pj.BondCompareParameters.RingMatchesRingOnly == true &&
              pj.BondCompareParameters.CompleteRingsOnly == true &&
              0 == strcmp(pj.InitialSeed.c_str(), "CNC"));
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGithubIssue481() {
  BOOST_LOG(rdInfoLog) << "\n-------------------------------------"
                       << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing github #481 : order dependence in FMCS with chirality"
      << std::endl;

  {
    std::string s1 = "O[C@H](C)Cl";
    std::string s2 = "O[C@@H](C)Cl";
    ROMOL_SPTR ptr1(SmilesToMol(s1));
    ROMOL_SPTR ptr2(SmilesToMol(s2));
    std::vector<ROMOL_SPTR> mols;
    mols.push_back(ptr1);
    mols.push_back(ptr2);

    BOOST_LOG(rdInfoLog) << "**** mols:" << s1 << "   " << s2 << "\n";
    //#688:
    std::vector<std::pair<int, int>> vect;
    ROMOL_SPTR mcs_mol(SmilesToMol("OCCl"));
    bool sub_res = SubstructMatch(*ptr2, *mcs_mol, vect, true, false);
    std::cout << "query=OCCl "
              << "SubstructMatch(useChirality=false) res =" << sub_res
              << " size=" << vect.size() << std::endl;
    sub_res = SubstructMatch(*ptr2, *mcs_mol, vect, true, true);
    std::cout << "query=OCCl "
              << "SubstructMatch(useChirality=true ) res =" << sub_res
              << " size=" << vect.size() << std::endl;

    {
      MCSParameters p;
      p.AtomCompareParameters.MatchChiralTag = false;
      p.BondCompareParameters.MatchStereo = false;
      MCSResult mcs_res = findMCS(mols, &p);
      BOOST_LOG(rdInfoLog) << "MCS: " << mcs_res.SmartsString << " "
                           << mcs_res.NumAtoms << " atoms, " << mcs_res.NumBonds
                           << " bonds\n";
      TEST_ASSERT(mcs_res.NumAtoms == 4);
      TEST_ASSERT(mcs_res.NumBonds == 3);

      p.AtomCompareParameters.MatchChiralTag = true;
      p.BondCompareParameters.MatchStereo = true;
      mcs_res = findMCS(mols, &p);
      BOOST_LOG(rdInfoLog) << "MCS: " << mcs_res.SmartsString << " "
                           << mcs_res.NumAtoms << " atoms, " << mcs_res.NumBonds
                           << " bonds\n";
      TEST_ASSERT(mcs_res.NumAtoms == 3);  //#688 was 2);
      TEST_ASSERT(mcs_res.NumBonds == 2);  //#688 was 1);
    }

    BOOST_LOG(rdInfoLog) << "------ REVERSE mols -------- \n";

    mols.clear();
    mols.push_back(ptr2);
    mols.push_back(ptr1);
    {
      MCSParameters p;
      p.AtomCompareParameters.MatchChiralTag = false;
      p.BondCompareParameters.MatchStereo = false;
      MCSResult mcs_res = findMCS(mols, &p);
      BOOST_LOG(rdInfoLog) << "MCS: " << mcs_res.SmartsString << " "
                           << mcs_res.NumAtoms << " atoms, " << mcs_res.NumBonds
                           << " bonds\n";
      TEST_ASSERT(mcs_res.NumAtoms == 4);
      TEST_ASSERT(mcs_res.NumBonds == 3);

      p.AtomCompareParameters.MatchChiralTag = true;
      p.BondCompareParameters.MatchStereo = true;
      mcs_res = findMCS(mols, &p);
      BOOST_LOG(rdInfoLog) << "MCS: " << mcs_res.SmartsString << " "
                           << mcs_res.NumAtoms << " atoms, " << mcs_res.NumBonds
                           << " bonds\n";
      TEST_ASSERT(mcs_res.NumAtoms == 3);  //#688 was 2);
      TEST_ASSERT(mcs_res.NumBonds == 2);  //#688 was 1);
    }
  }

  {
    std::string s1 = "OC(C)Cl";
    std::string s2 = "O[C@H](C)Cl";
    ROMOL_SPTR ptr1(SmilesToMol(s1));
    ROMOL_SPTR ptr2(SmilesToMol(s2));
    std::vector<ROMOL_SPTR> mols;
    mols.push_back(ptr1);
    mols.push_back(ptr2);

    BOOST_LOG(rdInfoLog) << "**** mols:" << s1 << "   " << s2 << "\n";

    {
      MCSParameters p;
      p.AtomCompareParameters.MatchChiralTag = false;
      p.BondCompareParameters.MatchStereo = false;
      MCSResult mcs_res = findMCS(mols, &p);
      BOOST_LOG(rdInfoLog) << "MCS: " << mcs_res.SmartsString << " "
                           << mcs_res.NumAtoms << " atoms, " << mcs_res.NumBonds
                           << " bonds\n";
      TEST_ASSERT(mcs_res.NumAtoms == 4);
      TEST_ASSERT(mcs_res.NumBonds == 3);

      p.AtomCompareParameters.MatchChiralTag = true;
      p.BondCompareParameters.MatchStereo = true;
      mcs_res = findMCS(mols, &p);
      BOOST_LOG(rdInfoLog) << "MCS: " << mcs_res.SmartsString << " "
                           << mcs_res.NumAtoms << " atoms, " << mcs_res.NumBonds
                           << " bonds\n";
      std::vector<std::pair<int, int>> vect;
      bool sub_res =
          SubstructMatch(*mols[1].get(), *mols[0].get(), vect, true, true);
      if (sub_res == false) {  // actualy == true & 4, 3 !!!
        TEST_ASSERT(mcs_res.NumAtoms == 0);
        TEST_ASSERT(mcs_res.NumBonds == 0);
      }
    }

    BOOST_LOG(rdInfoLog) << "------ REVERSE mols -------- \n";

    mols.clear();
    mols.push_back(ptr2);
    mols.push_back(ptr1);
    {
      MCSParameters p;
      p.AtomCompareParameters.MatchChiralTag = false;
      p.BondCompareParameters.MatchStereo = false;
      MCSResult mcs_res = findMCS(mols, &p);
      BOOST_LOG(rdInfoLog) << "MCS: " << mcs_res.SmartsString << " "
                           << mcs_res.NumAtoms << " atoms, " << mcs_res.NumBonds
                           << " bonds\n";
      TEST_ASSERT(mcs_res.NumAtoms == 4);
      TEST_ASSERT(mcs_res.NumBonds == 3);

      p.AtomCompareParameters.MatchChiralTag = true;
      p.BondCompareParameters.MatchStereo = true;
      mcs_res = findMCS(mols, &p);
      BOOST_LOG(rdInfoLog) << "MCS: " << mcs_res.SmartsString << " "
                           << mcs_res.NumAtoms << " atoms, " << mcs_res.NumBonds
                           << " bonds\n";
      std::vector<std::pair<int, int>> vect;
      bool sub_res =
          SubstructMatch(*mols[1].get(), *mols[0].get(), vect, true, true);
      if (sub_res == false) {
        TEST_ASSERT(mcs_res.NumAtoms == 0);
        TEST_ASSERT(mcs_res.NumBonds == 0);
      }
    }
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testInitialSeed() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "FMCS testInitialSeed()" << std::endl;

  std::vector<ROMOL_SPTR> mols;
  const char* smi[] = {
      "CC1CCC(N)CC1", "CC1CC(C)CC(C)C1",  // OK test.sdf
  };

  for (auto& i : smi) {
    std::string id;
    mols.push_back(ROMOL_SPTR(SmilesToMol(getSmilesOnly(i, &id))));
  }
  MCSParameters p;
  p.InitialSeed = "CC";
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 7 && res.NumBonds == 7);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testInitialSeed2() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "FMCS testInitialSeed2()" << std::endl;

  std::vector<ROMOL_SPTR> mols;
  const char* smi[] = {
      "Cc1c(F)c(N2CCNC(C)C2)cc2c1c(=O)c(C(=O)O)cn2C1CC1",
      "COc1c(N2CCNC(C)C2)c(F)cc2c(=O)c(C(=O)O)cn(C3CC3)c12",
  };
  const char* initial_smarts = "CCNCCNcccccccnC1CC1";
  BOOST_LOG(rdInfoLog) << "initial_smarts: " << initial_smarts << std::endl;

  for (auto& i : smi) {
    std::string id;
    mols.push_back(ROMOL_SPTR(SmilesToMol(getSmilesOnly(i, &id))));
    std::auto_ptr<ROMol> seed(SmartsToMol(initial_smarts));
    MatchVectType match;
    bool matched = SubstructMatch(*mols.back(), *seed, match);
    BOOST_LOG(rdInfoLog) << (matched ? "RDKit MATCHED " : "RDKit DISmatched ")
                         << i << std::endl;
  }
  MCSParameters p;
  p.InitialSeed = initial_smarts;
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms != 0);

  // Make Initial Seed from MCS
  p.Verbose = true;
  p.InitialSeed =
      "[#6]1-[#6]-[#7]-[#6](-[#6]-[#7]-1-[#6]1:[#6](:[#6]:[#6]2:[#6](:[#6](:[#"
      "6]:[#7](-[#6]3-[#6]-[#6]-3):[#6]:2:[#6]:1)-[#6](=[#8])-[#8])=[#8])-[#9])"
      "-[#6]";  // 25 atoms, 28 bonds
  BOOST_LOG(rdInfoLog)
      << "\n\nFound MCS as the only initial seed (25 atoms, 28 bonds): \n"
      << p.InitialSeed << std::endl;
  t0 = nanoClock();
  res = findMCS(mols, &p);
  BOOST_LOG(rdInfoLog) << "MCS: " << res.SmartsString << " " << res.NumAtoms
                       << " atoms, " << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms != 0);

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGithub631() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing github issue 631: FindMCS "
                          "matchChiralTag=True does not match self"
                       << std::endl;
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] = {
      // all examples derived from the bug report
      "CN(C)[C@@H]1CCCNC1",  // == MCS
      "C1C=CCN1[C@@H]1CCCNC1",
      "Cc1cc2c(cc1C)C(=O)N([C@@H]1CCC(=O)NC1=O)C2=O",
  };

  for (int pass = 0; pass < 2; ++pass, mols.clear())
    for (auto& i : smi) {
      RWMol* m = SmilesToMol(getSmilesOnly(i));
      TEST_ASSERT(m);

      if (0 == pass)
        mols.clear();  // use a pair of the same molecules only. On the second
                       // pass use all.

      mols.push_back(ROMOL_SPTR(m));
      mols.push_back(ROMOL_SPTR(new ROMol(*m)));
      {
        MCSParameters p;
        p.AtomCompareParameters.MatchChiralTag = false;
        //          p.Verbose = true;
        MCSResult res = findMCS(mols, &p);
        BOOST_LOG(rdInfoLog)
            << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n"
            << std::endl;
        ;

        TEST_ASSERT(res.NumAtoms == mols[0]->getNumAtoms());
        TEST_ASSERT(res.NumBonds == mols[0]->getNumBonds());
      }

      {
        MCSParameters p;
        p.AtomCompareParameters.MatchChiralTag = true;
        // p.BondCompareParameters.MatchStereo = useChirality;
        //        p.Verbose = true;
        MCSResult res = findMCS(mols, &p);
        BOOST_LOG(rdInfoLog)
            << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n"
            << std::endl;
        ;

        TEST_ASSERT(res.NumAtoms == mols[0]->getNumAtoms());
        TEST_ASSERT(res.NumBonds == mols[0]->getNumBonds());
      }
      BOOST_LOG(rdInfoLog) << "============================================"
                           << std::endl;
    }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testFormalChargeMatch() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing including formal charge in matches "
                       << std::endl;
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] = {"CCNC", "CCN(C)C", "CC[N+](C)C"};

  for (auto& i : smi) {
    RWMol* m = SmilesToMol(getSmilesOnly(i));
    TEST_ASSERT(m);

    mols.push_back(ROMOL_SPTR(m));
  }
  {
    // by default charge is not used.
    MCSParameters p;
    MCSResult res = findMCS(mols, &p);
    BOOST_LOG(rdInfoLog) << "MCS: " << res.SmartsString << " " << res.NumAtoms
                         << " atoms, " << res.NumBonds << " bonds\n"
                         << std::endl;
    ;

    TEST_ASSERT(res.NumAtoms == 4);
    TEST_ASSERT(res.NumBonds == 3);
  }
  {
    // check charge
    MCSParameters p;
    p.AtomCompareParameters.MatchFormalCharge = true;
    MCSResult res = findMCS(mols, &p);
    BOOST_LOG(rdInfoLog) << "MCS: " << res.SmartsString << " " << res.NumAtoms
                         << " atoms, " << res.NumBonds << " bonds\n"
                         << std::endl;
    ;

    TEST_ASSERT(res.NumAtoms == 2);
    TEST_ASSERT(res.NumBonds == 1);
  }

  BOOST_LOG(rdInfoLog) << "============================================"
                       << std::endl;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGithub2034() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github #2034: add test for ring--non-ring "
                          "matches at the atom level"
                       << std::endl;
  std::vector<ROMOL_SPTR> mols;
  const char* smi[] = {"C1CC1N2CC2", "C1CC1N"};

  for (auto& i : smi) {
    auto m = SmilesToMol(getSmilesOnly(i));
    TEST_ASSERT(m);

    mols.push_back(ROMOL_SPTR(m));
  }
  {
    // by default we're not doing ringMatchesRingOnly.
    MCSResult res = findMCS(mols);
    BOOST_LOG(rdInfoLog) << "MCS: " << res.SmartsString << " " << res.NumAtoms
                         << " atoms, " << res.NumBonds << " bonds\n"
                         << std::endl;

    TEST_ASSERT(res.NumAtoms == 4);
    TEST_ASSERT(res.NumBonds == 4);
  }

  {
    MCSParameters p;
    p.AtomCompareParameters.RingMatchesRingOnly = true;
    //          p.Verbose = true;
    MCSResult res = findMCS(mols, &p);
    BOOST_LOG(rdInfoLog) << "MCS: " << res.SmartsString << " " << res.NumAtoms
                         << " atoms, " << res.NumBonds << " bonds\n"
                         << std::endl;

    TEST_ASSERT(res.NumAtoms == 3);
    TEST_ASSERT(res.NumBonds == 3);
  }
  {
    // set it:
    bool maximizeBonds = true, verbose = false, matchValences = false,
         ringMatchesRingOnly = true;
    double threshold = 1.0;
    unsigned timeout = 3000;
    MCSResult res = findMCS(mols, maximizeBonds, threshold, timeout, verbose,
                            matchValences, ringMatchesRingOnly);
    BOOST_LOG(rdInfoLog) << "MCS: " << res.SmartsString << " " << res.NumAtoms
                         << " atoms, " << res.NumBonds << " bonds\n"
                         << std::endl;

    TEST_ASSERT(res.NumAtoms == 3);
    TEST_ASSERT(res.NumBonds == 3);
  }

  BOOST_LOG(rdInfoLog) << "============================================"
                       << std::endl;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

//====================================================================================================
//====================================================================================================

int main(int argc, const char* argv[]) {
  (void)argc;
  (void)argv;
  // p.Verbose = true;
  BOOST_LOG(rdInfoLog)
      << "*******************************************************\n";
  BOOST_LOG(rdInfoLog) << "FMCS Unit Test \n";

// use maximum CPU resoures to increase time measuring accuracy and stability in
// multi process environment
#ifdef WIN32
  //    SetPriorityClass (GetCurrentProcess(), REALTIME_PRIORITY_CLASS );
  SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_HIGHEST);
#else
  setpriority(PRIO_PROCESS, getpid(), -20);
#endif

  T0 = nanoClock();
  t0 = nanoClock();

  testJSONParameters();

  test1Basics();

  test32();
  test190();
  test3();

  testSimpleFast();
  testSimple();
  testSegFault();

  testThreshold();
  testRing1();

  testAtomCompareIsotopes();
  testAtomCompareAnyAtom();
  testAtomCompareAnyAtomBond();

  test18();
  test504();
  // very SLOW optional tests:
  //    test330();  // SLOW test
  //    test45();   // SLOW test

  testInitialSeed();
  testInitialSeed2();

  // chirality check:
  testGithubIssue481();
  testChirality();
  testGithub631();
  //---

  testFormalChargeMatch();
  testGithub2034();

  unsigned long long t1 = nanoClock();
  double sec = double(t1 - T0) / 1000000.;
  printf("TOTAL Time elapsed %.2lf seconds\n", sec);

  BOOST_LOG(rdInfoLog)
      << "*******************************************************\n";
  return 0;
}
