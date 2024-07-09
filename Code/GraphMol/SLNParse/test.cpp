// $Id$
//
//  Copyright (c) 2008, Novartis Institutes for BioMedical Research Inc.
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
//       products derived from this software without specific prior
//       written permission.
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
//  Created by Greg Landrum September, 2006
//
#include <RDGeneral/test.h>
#include <iostream>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SLNParse/SLNParse.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <RDGeneral/RDLog.h>
using namespace std;

void test1() {
  RDKit::RWMol *mol;
  std::string sln;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Test1 " << std::endl;

  sln = "CH4";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 1);

  delete mol;
  sln = "CH3CH3";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 2);
  TEST_ASSERT(mol->getRingInfo()->numRings() == 0);

  delete mol;
  sln = "C[1]H2CH2CH2@1";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 3);
  TEST_ASSERT(mol->getRingInfo()->numRings() == 1);

  delete mol;
  sln = "C[2]H2CH2CH2@2";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 3);
  TEST_ASSERT(mol->getRingInfo()->numRings() == 1);

  delete mol;
  sln = "C[200]H2CH2CH2@200";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 3);
  TEST_ASSERT(mol->getRingInfo()->numRings() == 1);

  delete mol;
  sln = "C[1:foo]H2CH2CH2@1";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 3);
  TEST_ASSERT(mol->getRingInfo()->numRings() == 1);

  delete mol;
  sln = "C[foo;bar=1;baz=bletch]H4";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 1);
  TEST_ASSERT(mol->getRingInfo()->numRings() == 0);

  delete mol;
  sln = "CH3CH(CH3)OCH3";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 5);
  TEST_ASSERT(mol->getRingInfo()->numRings() == 0);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getDegree() == 1);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getDegree() == 3);
  TEST_ASSERT(mol->getAtomWithIdx(2)->getDegree() == 1);
  TEST_ASSERT(mol->getAtomWithIdx(3)->getDegree() == 2);

  delete mol;
  sln = "CH3CH(CH3)OCH3.Cl";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 6);
  TEST_ASSERT(mol->getRingInfo()->numRings() == 0);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getDegree() == 1);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getDegree() == 3);
  TEST_ASSERT(mol->getAtomWithIdx(2)->getDegree() == 1);
  TEST_ASSERT(mol->getAtomWithIdx(3)->getDegree() == 2);
  TEST_ASSERT(mol->getAtomWithIdx(4)->getDegree() == 1);
  TEST_ASSERT(mol->getAtomWithIdx(5)->getDegree() == 0);

  delete mol;
  sln = "H-O-H";
  mol = RDKit::SLNToMol(sln, false);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 3);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getAtomicNum() == 1);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getAtomicNum() == 8);
  TEST_ASSERT(mol->getAtomWithIdx(2)->getAtomicNum() == 1);

  delete mol;
  sln = "HC(F)(Cl)Br";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 4);

  delete mol;
  sln = "CH3C(=O)H";
  mol = RDKit::SLNToMol(sln, false);  //,1);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 7);

  delete mol;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test2() {
  RDKit::RWMol *mol;
  std::string sln;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Test2: bond orders " << std::endl;

  sln = "CH3CH3";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 2);
  TEST_ASSERT(mol->getBondBetweenAtoms(0, 1)->getBondType() ==
              RDKit::Bond::SINGLE);

  delete mol;
  sln = "CH2-CH2";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 2);
  TEST_ASSERT(mol->getBondBetweenAtoms(0, 1)->getBondType() ==
              RDKit::Bond::SINGLE);

  delete mol;
  sln = "CH2=CH2";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 2);
  TEST_ASSERT(mol->getBondBetweenAtoms(0, 1)->getBondType() ==
              RDKit::Bond::DOUBLE);

  delete mol;
  sln = "CH#CH";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 2);
  TEST_ASSERT(mol->getBondBetweenAtoms(0, 1)->getBondType() ==
              RDKit::Bond::TRIPLE);
  delete mol;

  sln = "C[1]H-CH2-CH2-CH2-CH2-CH=@1";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 6);
  TEST_ASSERT(mol->getRingInfo()->numRings() == 1);
  TEST_ASSERT(mol->getBondBetweenAtoms(0, 1)->getBondType() ==
              RDKit::Bond::SINGLE);
  TEST_ASSERT(mol->getBondBetweenAtoms(1, 2)->getBondType() ==
              RDKit::Bond::SINGLE);
  TEST_ASSERT(mol->getBondBetweenAtoms(2, 3)->getBondType() ==
              RDKit::Bond::SINGLE);
  TEST_ASSERT(mol->getBondBetweenAtoms(3, 4)->getBondType() ==
              RDKit::Bond::SINGLE);
  TEST_ASSERT(mol->getBondBetweenAtoms(4, 5)->getBondType() ==
              RDKit::Bond::SINGLE);
  TEST_ASSERT(mol->getBondBetweenAtoms(5, 0)->getBondType() ==
              RDKit::Bond::DOUBLE);

  delete mol;
  sln = "C[1]H:CH:CH:CH:CH:CH:@1";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 6);
  TEST_ASSERT(mol->getRingInfo()->numRings() == 1);
  TEST_ASSERT(mol->getBondBetweenAtoms(0, 1)->getBondType() ==
              RDKit::Bond::AROMATIC);
  TEST_ASSERT(mol->getBondBetweenAtoms(1, 2)->getBondType() ==
              RDKit::Bond::AROMATIC);
  TEST_ASSERT(mol->getBondBetweenAtoms(2, 3)->getBondType() ==
              RDKit::Bond::AROMATIC);
  TEST_ASSERT(mol->getBondBetweenAtoms(3, 4)->getBondType() ==
              RDKit::Bond::AROMATIC);
  TEST_ASSERT(mol->getBondBetweenAtoms(4, 5)->getBondType() ==
              RDKit::Bond::AROMATIC);
  TEST_ASSERT(mol->getBondBetweenAtoms(5, 0)->getBondType() ==
              RDKit::Bond::AROMATIC);

  delete mol;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test3() {
  std::string pval;
  RDKit::RWMol *mol;
  std::string sln;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Test3: atom properties " << std::endl;

  sln = "C[-]H3";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 1);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getFormalCharge() == -1);

  delete mol;
  sln = "C[charge=-]H3";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 1);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getFormalCharge() == -1);

  delete mol;
  sln = "C[charge=-1]H3";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 1);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getFormalCharge() == -1);

  delete mol;
  sln = "C[CHARGE=-1]H3";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 1);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getFormalCharge() == -1);

  delete mol;
  sln = "C[chARgE=-1]H3";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 1);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getFormalCharge() == -1);

  delete mol;
  sln = "C[1:-2]H2";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 1);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getFormalCharge() == -2);

  delete mol;
  sln = "C[1:+]H5";
  bool sanitize = false;
  mol = RDKit::SLNToMol(sln, sanitize);
  TEST_ASSERT(mol);
  mol->debugMol(std::cerr);
  TEST_ASSERT(mol->getNumAtoms() == 6);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getFormalCharge() == 1);

  delete mol;
  sln = "C[1:+2]H6";
  mol = RDKit::SLNToMol(sln, sanitize);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 7);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getFormalCharge() == 2);

  delete mol;
  sln = "C[1:foo;bar=baz]H4";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 1);
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp("foo"));
  mol->getAtomWithIdx(0)->getProp("foo", pval);
  TEST_ASSERT(pval == "");
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp("bar"));
  mol->getAtomWithIdx(0)->getProp("bar", pval);
  TEST_ASSERT(pval == "baz");
  TEST_ASSERT(!mol->getAtomWithIdx(0)->hasProp("baz"));

  delete mol;
  sln = "H[I=2]-C(-H[I=2])(F)F";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 5);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getIsotope() == 2);
  TEST_ASSERT(mol->getAtomWithIdx(2)->getIsotope() == 2);

  delete mol;
  sln = "C[*]H3";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 1);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getNumExplicitHs() == 3);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getNumImplicitHs() == 0);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getNoImplicit());
  TEST_ASSERT(mol->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);

#if 1
  // FIX: this should be accepted

  delete mol;
  sln = "CH[I=2]";
  mol = RDKit::SLNToMol(sln);  //,true,1);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 2);
#endif

  // but this should not be accepted:
  delete mol;
  sln = "CH4[I=13]";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(!mol);
  delete mol;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test4() {
  RDKit::RWMol *mol;
  std::string sln;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Test4: hydrogen handling " << std::endl;

  sln = "CH4";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 1);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getNumExplicitHs() == 4);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getNumImplicitHs() == 0);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getNoImplicit());

  delete mol;
  sln = "C[-]H3";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 1);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getNumExplicitHs() == 3);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getNumImplicitHs() == 0);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getNoImplicit());

  delete mol;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test5() {
  std::string pval;
  RDKit::RWMol *patt, *mol;
  std::vector<RDKit::MatchVectType> mV;
  std::string sln, smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Test5: basic queries " << std::endl;

#if 1
  sln = "C[charge=+1]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  TEST_ASSERT(patt->getAtomWithIdx(0)->getNumExplicitHs() == 0);
  TEST_ASSERT(patt->getAtomWithIdx(0)->getNoImplicit());
  TEST_ASSERT(patt->getAtomWithIdx(0)->getNumImplicitHs() == 0);

  smi = "C[CH2+](C)C";
  bool sanitize = false;
  bool debugParse = false;
  mol = RDKit::SmilesToMol(smi, debugParse, sanitize);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);

  delete mol;
  smi = "C[CH](C)C";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(!RDKit::SubstructMatch(*mol, *patt, mV));

  delete mol;
  smi = "C[CH3+2](C)C";
  mol = RDKit::SmilesToMol(smi, debugParse, sanitize);
  TEST_ASSERT(mol);
  TEST_ASSERT(!RDKit::SubstructMatch(*mol, *patt, mV));

  delete mol;
  smi = "C(=O)OC";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);

  delete patt;
  sln = "AnyAny";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 2);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 2);

  delete patt;
  sln = "Any-Any";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 2);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 2);

  delete patt;
  sln = "Any=Any";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 2);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);

  delete patt;
  sln = "Any#Any";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 2);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 0);

  delete patt;
  sln = "Any:Any";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 2);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 0);

  delete patt;
  sln = "Any~Any";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 2);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 3);

  delete patt;
  sln = "C[charge=+1|charge=+2]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  delete mol;
  smi = "C[CH2+](C)[CH+2]";
  mol = RDKit::SmilesToMol(smi, debugParse, sanitize);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 2);
  delete mol;
  smi = "C[CH2+](C)[N+2]";
  mol = RDKit::SmilesToMol(smi, debugParse, sanitize);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  delete mol;
  smi = "C[N+](C)[CH+2]";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);

  delete patt;
  sln = "C[charge=+1;HC=2]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  delete mol;
  smi = "C[CH2+](CC)[NH+]";
  mol = RDKit::SmilesToMol(smi, debugParse, sanitize);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  TEST_ASSERT(mV[0][0].second == 1);

  delete patt;
  sln = "Any[charge=+1;HC=2]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  delete mol;
  smi = "C[CH2+](CC)[NH+]";
  mol = RDKit::SmilesToMol(smi, debugParse, sanitize);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  TEST_ASSERT(mV[0][0].second == 1);

  delete patt;
  sln = "Any[charge=+1;!HC=2]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  delete mol;
  smi = "C[CH2+](CC)[NH+]";
  mol = RDKit::SmilesToMol(smi, debugParse, sanitize);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  TEST_ASSERT(mV[0][0].second == 4);

  delete mol;
  smi = "[H]CC[H]";
  mol = RDKit::SmilesToMol(smi, 0, false);
  TEST_ASSERT(mol);
  RDKit::MolOps::sanitizeMol(*mol);
  TEST_ASSERT(mol->getNumAtoms() == 4);

  delete patt;
  sln = "AnyAny";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 2);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 3);
  delete patt;
  sln = "HevHev";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 2);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  delete patt;
  sln = "AnyHev";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 2);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 3);
  delete patt;

  delete mol;
  smi = "FC(Cl)(Br)CI";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  RDKit::MolOps::sanitizeMol(*mol);
  TEST_ASSERT(mol->getNumAtoms() == 6);

  sln = "HalC";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 2);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 4);
  delete patt;
#endif

  delete mol;
  smi = "CO[2H]";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  RDKit::MolOps::sanitizeMol(*mol);
  TEST_ASSERT(mol->getNumAtoms() == 3);

  sln = "H";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  delete patt;

  sln = "H[i=1]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 0);
  delete patt;

  sln = "H[i=2]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  delete patt;

  sln = "Het";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  delete patt;

  delete mol;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test6() {
  std::string pval;
  RDKit::RWMol *patt, *mol;
  std::vector<RDKit::MatchVectType> mV;
  std::string sln, smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Test6: ring queries " << std::endl;

  sln = "C[rbc=2]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  TEST_ASSERT(!patt->getAtomWithIdx(0)->hasProp("rbc"));

  smi = "C1CC1";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 3);

  delete mol;
  smi = "C1C2C1C2";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 2);

  delete mol;
  smi = "CCC";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(!RDKit::SubstructMatch(*mol, *patt, mV));

  delete patt;
  sln = "C[rbc=1]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);

  delete mol;
  smi = "C1CC1";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(!RDKit::SubstructMatch(*mol, *patt, mV));

  delete mol;
  smi = "CCC";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(!RDKit::SubstructMatch(*mol, *patt, mV));

  delete patt;
  sln = "C[rbc=3]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);

  delete mol;
  smi = "C1CC1";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(!RDKit::SubstructMatch(*mol, *patt, mV));

  delete mol;
  smi = "C1C2C1C2";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 2);

  // test case sensitivity, just because we should:
  delete patt;
  sln = "C[rBC=2]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  TEST_ASSERT(!patt->getAtomWithIdx(0)->hasProp("rbc"));

  delete mol;
  smi = "C1CC1C";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 3);

  delete patt;
  sln = "C[r]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);

  delete mol;
  smi = "C1CC1C";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 3);

  delete patt;
  sln = "C[1:rbc=f]CCC@1";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 4);

  delete mol;
  smi = "C1CCC1";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);

  delete mol;
  smi = "C1CC2C1CC2";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 2);

  delete patt;
  sln = "C[1:rbc=f]C[rbc=f]CC@1";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 4);

  delete mol;
  smi = "C1CCC1";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV));

  delete mol;
  smi = "C1CC2C1CC2";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 2);

  delete patt;
  sln = "C[1:rbc=f]C[rbc=f]C[rbc=f]C@1";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 4);

  delete mol;
  smi = "C1CCC1";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);

  delete mol;
  smi = "C1CC2C1CC2";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(!RDKit::SubstructMatch(*mol, *patt, mV));

  delete mol;
  delete patt;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test7() {
  std::string pval;
  RDKit::RWMol *patt, *mol;
  std::vector<RDKit::MatchVectType> mV;
  std::string sln, smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Test7: qualified queries " << std::endl;

  sln = "C[charge=+1]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);

  smi = "C[CH2+](C)C";
  bool debugParse = false;
  bool sanitize = false;
  mol = RDKit::SmilesToMol(smi, debugParse, sanitize);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);

  delete patt;
  sln = "C[charge<1]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);

  delete mol;
  smi = "C[CH2+](C)C";
  mol = RDKit::SmilesToMol(smi, debugParse, sanitize);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 3);

  delete mol;
  smi = "C[CH](C)C";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 4);

  delete mol;
  smi = "C[CH3+2](C)C";
  mol = RDKit::SmilesToMol(smi, debugParse, sanitize);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 3);

  delete patt;
  sln = "C[charge<=1]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);

  delete mol;
  smi = "C[CH2+](C)C";
  mol = RDKit::SmilesToMol(smi, debugParse, sanitize);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 4);

  delete mol;
  smi = "C[CH](C)C";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 4);

  delete mol;
  smi = "C[CH3+2](C)C";
  mol = RDKit::SmilesToMol(smi, debugParse, sanitize);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 3);

  delete patt;
  sln = "C[charge>=1]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);

  delete mol;
  smi = "C[CH2+](C)C";
  mol = RDKit::SmilesToMol(smi, debugParse, sanitize);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);

  delete mol;
  smi = "C[CH](C)C";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 0);

  delete mol;
  smi = "C[CH3+2](C)C";
  mol = RDKit::SmilesToMol(smi, debugParse, sanitize);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);

  delete patt;
  sln = "C[charge>1]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);

  delete mol;
  smi = "C[CH2+](C)C";
  mol = RDKit::SmilesToMol(smi, debugParse, sanitize);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 0);

  delete mol;
  smi = "C[CH](C)C";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 0);

  delete mol;
  smi = "C[CH3+2](C)C";
  mol = RDKit::SmilesToMol(smi, debugParse, sanitize);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);

  delete patt;
  sln = "C[charge!=0]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);

  delete mol;
  smi = "C[CH2+](C)C";
  mol = RDKit::SmilesToMol(smi, debugParse, sanitize);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);

  delete mol;
  smi = "C[CH](C)C";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 0);

  delete mol;
  smi = "C[CH3+2](C)C";
  mol = RDKit::SmilesToMol(smi, debugParse, sanitize);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);

  delete patt;
  sln = "C[charge!=1]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);

  delete mol;
  smi = "C[CH2+](C)C";
  mol = RDKit::SmilesToMol(smi, debugParse, sanitize);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 3);

  delete mol;
  smi = "C[CH](C)C";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 4);

  delete mol;
  smi = "C[CH3+2](C)C";
  mol = RDKit::SmilesToMol(smi, debugParse, sanitize);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 4);

  delete mol;
  delete patt;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test8() {
  std::string pval;
  RDKit::RWMol *patt, *mol;
  std::vector<RDKit::MatchVectType> mV;
  std::string sln, smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Test8: more complex atom properties " << std::endl;

  sln = "Any";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);

  smi = "C(=O)OC";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 4);
  TEST_ASSERT(mV[0][0].second == 0);
  TEST_ASSERT(mV[1][0].second == 1);
  TEST_ASSERT(mV[2][0].second == 2);
  TEST_ASSERT(mV[3][0].second == 3);

  delete patt;
  sln = "CO[F]C";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 3);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  TEST_ASSERT(mV.size() == 1);
  TEST_ASSERT(mV[0][0].second == 0);
  TEST_ASSERT(mV[0][1].second == 2);

  delete patt;
  sln = "Any[TBO=4]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 2);
  // CHECK: should TBO really match both Cs here or only atom 3 (i.e. should it
  // be degree or valence)?
  TEST_ASSERT(mV.size() == 2);
  TEST_ASSERT(mV[0][0].second == 0);
  TEST_ASSERT(mV[1][0].second == 3);

  delete patt;
  sln = "Any[TAC=2]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  TEST_ASSERT(mV.size() == 1);
  TEST_ASSERT(mV[0][0].second == 2);

  delete patt;
  sln = "Any[TAC>2]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 2);
  TEST_ASSERT(mV.size() == 2);
  TEST_ASSERT(mV[0][0].second == 0);
  TEST_ASSERT(mV[1][0].second == 3);

  delete patt;
  sln = "Any[HC=1]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  TEST_ASSERT(mV.size() == 1);
  TEST_ASSERT(mV[0][0].second == 0);

  delete patt;
  sln = "Any[HC=0]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 2);
  TEST_ASSERT(mV.size() == 2);
  TEST_ASSERT(mV[0][0].second == 1);
  TEST_ASSERT(mV[1][0].second == 2);

  delete patt;
  sln = "Any[HC=F]H3";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  TEST_ASSERT(mV.size() == 1);
  TEST_ASSERT(mV[0][0].second == 3);

  delete patt;
  sln = "Any[HAC>1]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 2);
  TEST_ASSERT(mV.size() == 2);
  TEST_ASSERT(mV[0][0].second == 0);
  TEST_ASSERT(mV[1][0].second == 2);

  delete patt;
  sln = "Any[HAC=1]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 2);
  TEST_ASSERT(mV.size() == 2);
  TEST_ASSERT(mV[0][0].second == 1);
  TEST_ASSERT(mV[1][0].second == 3);

  delete patt;
  sln = "Any[HAC=F]~Any";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 2);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 2);
  TEST_ASSERT(mV.size() == 2);
  TEST_ASSERT(mV[0][0].second == 1);
  TEST_ASSERT(mV[1][0].second == 3);

  delete patt;
  sln = "AnyAny[TAC=F]Any";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 3);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  TEST_ASSERT(mV.size() == 1);
  TEST_ASSERT(mV[0][1].second == 2);

  delete mol;
  smi = "CCC=C";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 0);

  delete patt;
  sln = "CC";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 2);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 2);

  delete patt;
  sln = "CC[F]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 2);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 0);

  delete patt;
  sln = "CC[F]H3";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 2);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);

  delete patt;
  delete mol;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test9() {
  std::string pval;
  RDKit::RWMol *patt, *mol;
  std::vector<RDKit::MatchVectType> mV;
  std::string sln, smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Test9: recursive SLNs " << std::endl;

  sln = "Any[is=Cl,Br,I]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);

  smi = "C(=O)Cl";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  TEST_ASSERT(mV.size() == 1);
  TEST_ASSERT(mV[0].size() == 1);
  TEST_ASSERT(mV[0][0].second == 2);

  delete mol;
  smi = "C(=O)Br";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  TEST_ASSERT(mV.size() == 1);
  TEST_ASSERT(mV[0].size() == 1);
  TEST_ASSERT(mV[0][0].second == 2);

  delete mol;
  smi = "C(=O)I";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  TEST_ASSERT(mV.size() == 1);
  TEST_ASSERT(mV[0].size() == 1);
  TEST_ASSERT(mV[0][0].second == 2);

  delete mol;
  smi = "C(=O)C";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 0);

  delete patt;
  sln = "C[is=C=O]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);

  delete mol;
  smi = "C(=O)OC";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  TEST_ASSERT(mV.size() == 1);
  TEST_ASSERT(mV[0].size() == 1);
  TEST_ASSERT(mV[0][0].second == 0);

  delete patt;
  sln = "Any[is=C,O;charge=-1]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  delete mol;
  smi = "C(=O)OC";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 0);
  delete mol;
  smi = "C(=O)[O-]";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  TEST_ASSERT(mV.size() == 1);
  TEST_ASSERT(mV[0].size() == 1);
  TEST_ASSERT(mV[0][0].second == 2);

  delete patt;
  // make sure we aren't case sensitive:
  sln = "Any[Is=C,O;charge=-1]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  delete mol;
  smi = "[C-](=O)O";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  TEST_ASSERT(mV.size() == 1);
  TEST_ASSERT(mV[0].size() == 1);
  TEST_ASSERT(mV[0][0].second == 0);

  delete patt;
  // check handling of 'not':
  sln = "Any[not=N,O]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  delete mol;
  smi = "C(=O)OCN";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 2);
  TEST_ASSERT(mV.size() == 2);
  TEST_ASSERT(mV[0].size() == 1);
  TEST_ASSERT(mV[0][0].second == 0);
  TEST_ASSERT(mV[1].size() == 1);
  TEST_ASSERT(mV[1][0].second == 3);

  delete patt;
  delete mol;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test10() {
  std::string pval;
  RDKit::RWMol *patt, *mol;
  std::vector<RDKit::MatchVectType> mV;
  std::string sln, smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Test10: nested recursive SLNs " << std::endl;

  // do recursions in the 'is':
  sln = "Any[is=C[is=Any=O],O]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  smi = "C(=O)OC";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 3);
  TEST_ASSERT(mV.size() == 3);
  TEST_ASSERT(mV[0].size() == 1);
  TEST_ASSERT(mV[0][0].second == 0);
  TEST_ASSERT(mV[1].size() == 1);
  TEST_ASSERT(mV[1][0].second == 1);
  TEST_ASSERT(mV[2].size() == 1);
  TEST_ASSERT(mV[2][0].second == 2);

  // do recursions in the 'not':
  delete patt;
  sln = "Any[not=C[is=Any=O],O]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  delete mol;
  smi = "C(=O)OC";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  TEST_ASSERT(mV.size() == 1);
  TEST_ASSERT(mV[0].size() == 1);
  TEST_ASSERT(mV[0][0].second == 3);
  delete patt;

  // move the anchor:
  sln = "C[is=C=O]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  delete mol;
  smi = "C(=O)OC";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  TEST_ASSERT(mV.size() == 1);
  TEST_ASSERT(mV[0].size() == 1);
  TEST_ASSERT(mV[0][0].second == 0);
  delete patt;

  sln = "C[is=O=C*]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  TEST_ASSERT(mV.size() == 1);
  TEST_ASSERT(mV[0].size() == 1);
  TEST_ASSERT(mV[0][0].second == 0);
  delete patt;

  sln = "C[is=O=C]";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 1);
  TEST_ASSERT(!RDKit::SubstructMatch(*mol, *patt, mV));
  delete patt;

  delete mol;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test11() {
  std::string pval, cip;
  RDKit::RWMol *mol;
  std::vector<RDKit::MatchVectType> mV;
  std::string sln, smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Test11: chiral SLNs " << std::endl;

  sln = "CH(Cl)(F)Br";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 4);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag() ==
              RDKit::Atom::CHI_UNSPECIFIED);

  delete mol;
  sln = "C[s=N]H(Cl)(F)Br";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 4);
  std::cerr << sln << " -> " << MolToSmiles(*mol, true) << " "
            << mol->getNumAtoms() << std::endl;
  smi = MolToSmiles(*mol, true);
  // TEST_ASSERT(smi=="F[C@@H](Cl)Br");

  delete mol;
  sln = "ClC[s=i]H(F)Br";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 4);
// mol->debugMol(std::cerr);
#if 0
  RDKit::MolOps::assignAtomChiralCodes(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(RDKit::common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(RDKit::common_properties::_CIPCode,cip);
  TEST_ASSERT(cip=="R");
#endif
  std::cerr << sln << " -> " << MolToSmiles(*mol, true) << " "
            << mol->getNumAtoms() << std::endl;
  smi = MolToSmiles(*mol, true);
  // TEST_ASSERT(smi=="F[C@@H](Cl)Br");

  delete mol;
  sln = "FC[s=N]H(Cl)Br";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  std::cerr << sln << " -> " << MolToSmiles(*mol, true) << " "
            << mol->getNumAtoms() << std::endl;
  smi = MolToSmiles(*mol, true);
  // TEST_ASSERT(smi=="F[C@@H](Cl)Br");

  delete mol;
  sln = "FC[s=N]H(Br)Cl";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 4);
  std::cerr << sln << " -> " << MolToSmiles(*mol, true) << " "
            << mol->getNumAtoms() << std::endl;
  smi = MolToSmiles(*mol, true);
  // TEST_ASSERT(smi=="F[C@H](Cl)Br");

  delete mol;
  sln = "HC[s=N](Cl)(F)Br";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 4);
  std::cerr << sln << " -> " << MolToSmiles(*mol, true) << " "
            << mol->getNumAtoms() << std::endl;
  smi = MolToSmiles(*mol, true);
  // TEST_ASSERT(smi=="F[C@@H](Cl)Br");

  delete mol;
  sln = "C[s=i]H(Cl)(F)Br";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 4);
  smi = MolToSmiles(*mol, true);
  std::cerr << sln << " -> " << MolToSmiles(*mol, true) << " "
            << mol->getNumAtoms() << std::endl;
  // TEST_ASSERT(smi=="F[C@H](Cl)Br");

  delete mol;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test12() {
  RDKit::RWMol *patt, *mol;
  std::vector<RDKit::MatchVectType> mV;
  std::string sln, smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Test12: bond queries and properties " << std::endl;

  sln = "Any-Any";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 2);

  smi = "C=CN";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  TEST_ASSERT(mV.size() == 1);
  TEST_ASSERT(mV[0].size() == 2);
  TEST_ASSERT(mV[0][0].second == 1);
  TEST_ASSERT(mV[0][1].second == 2);

  delete patt;
  sln = "Any~Any";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 2);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 2);
  TEST_ASSERT(mV[0].size() == 2);
  TEST_ASSERT(mV[0][0].second == 0);
  TEST_ASSERT(mV[0][1].second == 1);
  TEST_ASSERT(mV[1].size() == 2);
  TEST_ASSERT(mV[1][0].second == 1);
  TEST_ASSERT(mV[1][1].second == 2);

  delete patt;
  sln = "Any-=Any";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 2);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 2);
  TEST_ASSERT(mV[0].size() == 2);
  TEST_ASSERT(mV[0][0].second == 0);
  TEST_ASSERT(mV[0][1].second == 1);
  TEST_ASSERT(mV[1].size() == 2);
  TEST_ASSERT(mV[1][0].second == 1);
  TEST_ASSERT(mV[1][1].second == 2);

  delete patt;
  sln = "Any=-Any";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 2);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 2);
  TEST_ASSERT(mV[0].size() == 2);
  TEST_ASSERT(mV[0][0].second == 0);
  TEST_ASSERT(mV[0][1].second == 1);
  TEST_ASSERT(mV[1].size() == 2);
  TEST_ASSERT(mV[1][0].second == 1);
  TEST_ASSERT(mV[1][1].second == 2);

  delete patt;
  sln = "Any-:Any";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 2);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  TEST_ASSERT(mV[0].size() == 2);
  TEST_ASSERT(mV[0][0].second == 1);
  TEST_ASSERT(mV[0][1].second == 2);

  sln = "Any-[type=2]Any";
  delete patt;
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 2);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  TEST_ASSERT(mV[0].size() == 2);
  TEST_ASSERT(mV[0][0].second == 0);
  TEST_ASSERT(mV[0][1].second == 1);
  delete patt;

  sln = "Any-[type=2|type=1]Any";
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 2);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 2);
  TEST_ASSERT(mV[0].size() == 2);
  TEST_ASSERT(mV[0][0].second == 0);
  TEST_ASSERT(mV[0][1].second == 1);
  TEST_ASSERT(mV[1].size() == 2);
  TEST_ASSERT(mV[1][0].second == 1);
  TEST_ASSERT(mV[1][1].second == 2);

  sln = "O~[r]C";
  delete patt;
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 2);

  delete mol;
  smi = "O=CC1COC1";
  mol = RDKit::SmilesToMol(smi);
  TEST_ASSERT(mol);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 2);
  TEST_ASSERT(mV[0].size() == 2);
  TEST_ASSERT(mV[0][0].second == 4);
  TEST_ASSERT(mV[1][0].second == 4);

  sln = "O~[!r]C";
  delete patt;
  patt = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(patt);
  TEST_ASSERT(patt->getNumAtoms() == 2);
  TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
  TEST_ASSERT(mV[0].size() == 2);
  TEST_ASSERT(mV[0][0].second == 0);

  delete mol;
  delete patt;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test13() {
  RDKit::RWMol *mol;
  std::string sln;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Test13: ring closure details " << std::endl;

  sln = "C[1]H2CH2CH2CH2CH2CH2@1";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 6);
  TEST_ASSERT(mol->getRingInfo()->numRings() == 1);
  TEST_ASSERT(mol->getRingInfo()->atomRings()[0].size() == 6);

  delete mol;
  sln = "C[1]H2CH2(CH2CH2CH2CH2@1)";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 6);
  TEST_ASSERT(mol->getRingInfo()->numRings() == 1);
  TEST_ASSERT(mol->getRingInfo()->atomRings()[0].size() == 6);

  delete mol;
  sln = "CH2(C[1]H2)CH2(CH2CH2CH2@1)";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 6);
  TEST_ASSERT(mol->getRingInfo()->numRings() == 1);
  TEST_ASSERT(mol->getRingInfo()->atomRings()[0].size() == 6);

  delete mol;
  sln = "C[1]H2CH2CH2CH2C[2]HCH@1CH2CH2@2";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 8);
  TEST_ASSERT(mol->getRingInfo()->numRings() == 2);
  TEST_ASSERT(mol->getRingInfo()->atomRings()[0].size() == 6);
  TEST_ASSERT(mol->getRingInfo()->atomRings()[1].size() == 4);

  delete mol;
  sln = "C[1]H2(CH2(CH2(CH2(C[2]H(CH@1(CH2(CH2@2)))))))";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 8);
  TEST_ASSERT(mol->getRingInfo()->numRings() == 2);
  TEST_ASSERT(mol->getRingInfo()->atomRings()[0].size() == 6);
  TEST_ASSERT(mol->getRingInfo()->atomRings()[1].size() == 4);

  delete mol;
  sln = "C[1](CH2CH2CH(CH2CH2@1)CH2CH2@1)Cl";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 9);
  TEST_ASSERT(mol->getRingInfo()->numRings() == 3);
  TEST_ASSERT(mol->getRingInfo()->atomRings()[0].size() == 6);
  TEST_ASSERT(mol->getRingInfo()->atomRings()[1].size() == 6);
  TEST_ASSERT(mol->getRingInfo()->atomRings()[2].size() == 6);

  delete mol;
  sln = "CH2(CH2@1)CH2(C[1]H2)";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(!mol);

  delete mol;
  sln = "CH2(CH2@1)";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(!mol);

  delete mol;
  sln = "C[1]H-CH2-C(@1)=O";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 4);
  TEST_ASSERT(mol->getRingInfo()->numRings() == 1);
  TEST_ASSERT(mol->getBondBetweenAtoms(0, 1)->getBondType() ==
              RDKit::Bond::SINGLE);
  TEST_ASSERT(mol->getBondBetweenAtoms(1, 2)->getBondType() ==
              RDKit::Bond::SINGLE);
  TEST_ASSERT(mol->getBondBetweenAtoms(2, 0)->getBondType() ==
              RDKit::Bond::SINGLE);
  TEST_ASSERT(mol->getBondBetweenAtoms(2, 3)->getBondType() ==
              RDKit::Bond::DOUBLE);

  delete mol;
  sln = "C[1]H-CH2-CH(@1)O";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 4);
  TEST_ASSERT(mol->getRingInfo()->numRings() == 1);
  TEST_ASSERT(mol->getBondBetweenAtoms(0, 1)->getBondType() ==
              RDKit::Bond::SINGLE);
  TEST_ASSERT(mol->getBondBetweenAtoms(1, 2)->getBondType() ==
              RDKit::Bond::SINGLE);
  TEST_ASSERT(mol->getBondBetweenAtoms(2, 0)->getBondType() ==
              RDKit::Bond::SINGLE);
  TEST_ASSERT(mol->getBondBetweenAtoms(2, 3)->getBondType() ==
              RDKit::Bond::SINGLE);

  delete mol;
  sln = "C[1]H-CH2-C(=@1)O";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 4);
  TEST_ASSERT(mol->getRingInfo()->numRings() == 1);
  TEST_ASSERT(mol->getBondBetweenAtoms(0, 1)->getBondType() ==
              RDKit::Bond::SINGLE);
  TEST_ASSERT(mol->getBondBetweenAtoms(1, 2)->getBondType() ==
              RDKit::Bond::SINGLE);
  TEST_ASSERT(mol->getBondBetweenAtoms(2, 0)->getBondType() ==
              RDKit::Bond::DOUBLE);
  TEST_ASSERT(mol->getBondBetweenAtoms(2, 3)->getBondType() ==
              RDKit::Bond::SINGLE);

  delete mol;
  sln = "C[1]H-CH2-C(O)=@1";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 4);
  TEST_ASSERT(mol->getRingInfo()->numRings() == 1);
  TEST_ASSERT(mol->getBondBetweenAtoms(0, 1)->getBondType() ==
              RDKit::Bond::SINGLE);
  TEST_ASSERT(mol->getBondBetweenAtoms(1, 2)->getBondType() ==
              RDKit::Bond::SINGLE);
  TEST_ASSERT(mol->getBondBetweenAtoms(2, 0)->getBondType() ==
              RDKit::Bond::DOUBLE);
  TEST_ASSERT(mol->getBondBetweenAtoms(2, 3)->getBondType() ==
              RDKit::Bond::SINGLE);

  delete mol;
  sln = "C[1]H-CH2-C(O)(=@1)";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 4);
  TEST_ASSERT(mol->getRingInfo()->numRings() == 1);
  TEST_ASSERT(mol->getBondBetweenAtoms(0, 1)->getBondType() ==
              RDKit::Bond::SINGLE);
  TEST_ASSERT(mol->getBondBetweenAtoms(1, 2)->getBondType() ==
              RDKit::Bond::SINGLE);
  TEST_ASSERT(mol->getBondBetweenAtoms(2, 0)->getBondType() ==
              RDKit::Bond::DOUBLE);
  TEST_ASSERT(mol->getBondBetweenAtoms(2, 3)->getBondType() ==
              RDKit::Bond::SINGLE);

  delete mol;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test14() {
  RDKit::RWMol *mol;
  std::string sln;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Test14: error catching " << std::endl;

  sln = "CH2(C@1H2)CH2(CH2CH2C[1]H2)";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(!mol);

  delete mol;
  sln = "CH2(CH2[1])CH2(CH2CH2CH2@1)";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(!mol);

  delete mol;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test15() {
  RDKit::RWMol *mol;
  std::string sln;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Test15: CTAB properties " << std::endl;

  {
    sln = "CH4<blah>";
    mol = RDKit::SLNToMol(sln);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->hasProp("blah"));
    delete mol;
  }

  {
    sln = "CH4<name=methane>";
    mol = RDKit::SLNToMol(sln);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->hasProp(RDKit::common_properties::_Name));
    std::string sval;
    mol->getProp(RDKit::common_properties::_Name, sval);
    TEST_ASSERT(sval == "methane");
    delete mol;
  }

  {
    sln = "CH4<blah;foo=\"1\">";
    mol = RDKit::SLNToMol(sln);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->hasProp("blah"));
    TEST_ASSERT(mol->hasProp("foo"));
    std::string sval;
    mol->getProp("foo", sval);
    TEST_ASSERT(sval == "1");
    delete mol;
  }
  {
    sln = "CH4<blah;foo=1>";
    mol = RDKit::SLNToMol(sln);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->hasProp("blah"));
    TEST_ASSERT(mol->hasProp("foo"));
    std::string sval;
    mol->getProp("foo", sval);
    TEST_ASSERT(sval == "1");
    delete mol;
  }
  {
    sln =
        "CH4<name=\"methane\";blah;coord2d=(1,0);too.small;test=lots and-lots "
        "of special,characters all at once.>";
    mol = RDKit::SLNToMol(sln);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->hasProp("blah"));
    std::string sval;
    mol->getProp(RDKit::common_properties::_Name, sval);
    TEST_ASSERT(sval == "methane");
    mol->getProp("coord2d", sval);
    TEST_ASSERT(sval == "(1,0)");
    TEST_ASSERT(mol->hasProp("too.small"));
    TEST_ASSERT(mol->hasProp("test"));
    mol->getProp("test", sval);
    TEST_ASSERT(sval == "lots and-lots of special,characters all at once.");
    delete mol;
  }

  {
    // Though it isn't part of the spec,
    // sometimes we have multiple ctab property blocks:
    sln = "CH4<foo=1><bar=2>";
    mol = RDKit::SLNToMol(sln);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->hasProp("foo"));
    TEST_ASSERT(mol->hasProp("bar"));
    std::string sval;
    mol->getProp("foo", sval);
    TEST_ASSERT(sval == "1");
    mol->getProp("bar", sval);
    TEST_ASSERT(sval == "2");
    delete mol;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test16() {
  RDKit::RWMol *mol;
  std::string sln;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Test16: simple macro atoms " << std::endl;

  {
    sln = "CH3ZCH3{Z:O}";
    mol = RDKit::SLNToMol(sln);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 3);
    TEST_ASSERT(mol->getAtomWithIdx(1)->getAtomicNum() == 8);
    delete mol;
  }
  {
    sln = "CH3ZGCH3{Z:O}{G:N}";
    mol = RDKit::SLNToMol(sln);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 4);
    TEST_ASSERT(mol->getAtomWithIdx(1)->getAtomicNum() == 8);
    TEST_ASSERT(mol->getAtomWithIdx(2)->getAtomicNum() == 7);
    delete mol;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testIssue278() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Test issue278: handling of 'not' and 'is' queries on Any"
      << std::endl;

  {
    std::string sln = "Any[IS=C(=O)]";
    RDKit::RWMol *mol = RDKit::SLNQueryToMol(sln);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 1);
    TEST_ASSERT(mol->getAtomWithIdx(0)->hasQuery());
    delete mol;
  }
  {
    std::string sln = "Any[NOT=C(=O)]";
    RDKit::RWMol *mol = RDKit::SLNQueryToMol(sln);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 1);
    TEST_ASSERT(mol->getAtomWithIdx(0)->hasQuery());
    delete mol;
  }
  {
    std::string sln = "Any[IS=C(=O)]";
    RDKit::RWMol *mol = RDKit::SLNToMol(sln);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 1);
    TEST_ASSERT(mol->getAtomWithIdx(0)->hasQuery());
    delete mol;
  }
  {
    std::string sln = "Any[NOT=C(=O)]";
    RDKit::RWMol *mol = RDKit::SLNToMol(sln);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 1);
    TEST_ASSERT(mol->getAtomWithIdx(0)->hasQuery());
    delete mol;
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}
void testIssue277() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Test issue277: parse error with & " << std::endl;

  {
    std::string sln = "Any[NOT=C,IS=O]";
    RDKit::RWMol *mol = RDKit::SLNQueryToMol(sln);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 1);
    TEST_ASSERT(mol->getAtomWithIdx(0)->hasQuery());
    delete mol;
  }
  {
    std::string sln = "Any[NOT=C;IS=O]";
    RDKit::RWMol *mol = RDKit::SLNQueryToMol(sln);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 1);
    TEST_ASSERT(mol->getAtomWithIdx(0)->hasQuery());
    delete mol;
  }
  {
    std::string sln = "Any[NOT=C&IS=O]";
    RDKit::RWMol *mol = RDKit::SLNQueryToMol(sln);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 1);
    TEST_ASSERT(mol->getAtomWithIdx(0)->hasQuery());
    delete mol;
  }
  {
    std::string sln = "Any[NOT=C,N&IS=O,S]";
    RDKit::RWMol *mol = RDKit::SLNQueryToMol(sln);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 1);
    TEST_ASSERT(mol->getAtomWithIdx(0)->hasQuery());
    delete mol;
  }
  {
    std::string sln = "Hev[!r;NOT=N]";
    RDKit::RWMol *patt = RDKit::SLNQueryToMol(sln);
    TEST_ASSERT(patt);
    TEST_ASSERT(patt->getNumAtoms() == 1);
    TEST_ASSERT(patt->getAtomWithIdx(0)->hasQuery());

    std::string smi = "CC1CC1N";
    RDKit::RWMol *mol = RDKit::SmilesToMol(smi);
    std::vector<RDKit::MatchVectType> mV;
    TEST_ASSERT(mol);
    TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
    TEST_ASSERT(mV.size() == 1);
    TEST_ASSERT(mV[0].size() == 1);
    TEST_ASSERT(mV[0][0].second == 0);

    delete mol;
    delete patt;
  }
  {
    std::string sln = "Hev[!r&NOT=N]";
    RDKit::RWMol *patt = RDKit::SLNQueryToMol(sln);
    TEST_ASSERT(patt);
    TEST_ASSERT(patt->getNumAtoms() == 1);
    TEST_ASSERT(patt->getAtomWithIdx(0)->hasQuery());

    std::string smi = "CC1CC1N";
    RDKit::RWMol *mol = RDKit::SmilesToMol(smi);
    std::vector<RDKit::MatchVectType> mV;
    TEST_ASSERT(mol);
    TEST_ASSERT(RDKit::SubstructMatch(*mol, *patt, mV) == 1);
    TEST_ASSERT(mV.size() == 1);
    TEST_ASSERT(mV[0].size() == 1);
    TEST_ASSERT(mV[0][0].second == 0);

    delete mol;
    delete patt;
  }

  {
    // examples from Pat Walters, taken from
    // http://pubs.acs.org/doi/abs/10.1021/ci300461a
    std::string slns[] = {
        "HetC(=O)O-[!R]C[8]:Hev(Any[NOT=O,S[TAC=2],C[TAC=4],N[TAC=3]]):Hev:Hev("
        "Any[IS=Hal,C#N,C(F)(F)F,S(=O)=O,C(=O)&NOT=C(=O)OH]):Hev:Hev(Any[NOT=O,"
        "S[TAC=2],C[TAC=4],N[TAC=3]]):@8",
        "HetC(=O)O-[!R]C[8]:Hev(Any[IS=Hal,C#N,C(F)(F)F,S(=O)=O,C(=O)&NOT=C(=O)"
        "OH]):Hev:Hev(Any[NOT=O,S[TAC=2],C[TAC=4],N[TAC=3]]):Hev:Hev(Any[NOT=O,"
        "S[TAC=2],C[TAC=4],N[TAC=3]]):@8",
        "CCH=[!R]C(Any[IS=H,C])Any[IS=C#N,C(=O)&NOT=C(=O)Any[IS=N,O]]",
        "Any[IS=C,O]CH=[!R]C(Any[IS=C(=O),S(=O),C#N,Hal,C(Hal)(Hal)Hal&NOT=C(="
        "O)OH])Any[IS=C(=O),S(=O),C#N&NOT=C(=O)OH]",
        "C[1](C(=O)OC[5]:C:C:C:C:C:@5CH=@1)Any[IS=C(=O),C(=S),S(=O),C#N,Hal,C("
        "Hal)(Hal)Hal&NOT=C(=O)OH]",
        "C[1](=CHOC[5]:C:C:C:C:C:@5C@1=O)Any[IS=C(=O),C(=S),S(=O),C#N,Hal,C("
        "Hal)(Hal)Hal&NOT=C(=O)OH]",
        "C[1]:N:Any(Any[IS=Hal,S(=O)(=O)C]):Any:Any(Any[IS=Hal,C(C)=NO,C#N,C(="
        "O),C(F)(F)F,S(=O)=O&NOT=C(=O)OH]):Any:@1",
        "C[1]:N:Any(Any[IS=Hal,C(C)=NO,C#N,C(=O),C(F)(F)F,S(=O)=O&NOT=C(=O)OH])"
        ":Any:Any(Any[IS=Hal,S(=O)(=O)C]):Any:@1",
        "C[1]:N:Any(Any[IS=Hal,S(=O)(=O)C]):Any(Any[IS=Hal,C(C)=NO,C#N,C(=O),C("
        "F)(F)F,S(=O)=O&NOT=C(=O)OH]):Any:Any:@1",
        "C[1]:N:Any(Any[IS=Hal,S(=O)(=O)C]):Any:Any:Any(Any[IS=Hal,C(C)=NO,C#N,"
        "C(=O),C(F)(F)F,S(=O)=O&NOT=C(=O)OH]):@1",
        "C[1](Any[IS=Hal,C(C)=NO,C#N,C(=O),C(F)(F)F,S(=O)=O&NOT=C(=O)OH]):N:"
        "Any(Any[IS=Hal,S(=O)(=O)C]):Any(Any[NOT=N]):Any(Any[NOT=N]):Any:@1",
        "C[1]:N:Any:Any(Any[IS=Hal,C(C)=NO,C#N,C(=O),C(F)(F)F,S(=O)=O&NOT=C(=O)"
        "OH]):Any(Any[IS=Hal,S(=O)(=O)C]):Any:@1",
        "C[1](Any[NOT=N,O]):C(Any[IS=Hal,C#N,C(=O),C(F)(F)F,S(=O)=O&NOT=C(=O)"
        "OH]):C(Any[IS=Hal,S(=O)(=O)C,C[r](=O)NC]):C(Any[NOT=N,O]):C(Any[NOT=N,"
        "O]):C(Any[IS=Hal,C#N,C(=O),C(F)(F)F,S(=O)=O&NOT=C(=O)OH]):@1",
        "C[1](Any[IS=Hal,C#N,C(=O),C(F)(F)F,S(=O)=O&NOT=C(=O)OH]):C(Any[IS=Hal,"
        "S(=O)(=O)C,C[r](=O)NC]):C(Any[IS=Hal,C#N,C(=O),C(F)(F)F,S(=O)=O&NOT=C("
        "=O)OH]):C(Any[NOT=N,O]):C(Any[NOT=N,O]):C(Any[NOT=N,O]):@1",
        "N(Any[IS=H,C[TAC=4]&NOT=C[TAC=4]-[R]C[TAC=4]N])(Any[IS=H,C[TAC=4]&NOT="
        "C[TAC=4]-[R]C[TAC=4]N])(Any[IS=H,C[TAC=4]&NOT=C[TAC=4]-[R]C[TAC=4]N])<"
        "max=1>",
        "Any[IS=H,C&NOT=C=O]N[!r](Any[IS=H,C&NOT=C=O])C(=O)C<max=2>",
        "Hev[!r&NOT=NC(=O)NC(=O)]Hev[!r&NOT=NC(=O)NC(=O)]Hev[!r&NOT=NC(=O)NC(="
        "O)]Hev[!r&NOT=NC(=O)NC(=O)]Hev[!r&NOT=NC(=O)NC(=O)]Hev[!r&NOT=NC(=O)"
        "NC(=O)]Hev[!r&NOT=NC(=O)NC(=O)]Hev[!r&NOT=NC(=O)NC(=O)]",
        "Hev[!r&NOT=C=O,S(=O)(=O),N*S(=O),N*(C=O)]Hev[!r&NOT=C=O,S(=O)(=O),N*S("
        "=O),N*(C=O)]Hev[!r&NOT=C=O,S(=O)(=O),N*S(=O),N*(C=O)]Hev[!r&NOT=C=O,S("
        "=O)(=O),N*S(=O),N*(C=O)]Hev[!r&NOT=C=O,S(=O)(=O)]Any[IS=CH3,OH,NH2,N("
        "CH3)CH3]",
        "EOF"};
    unsigned int i = 0;
    while (slns[i] != "EOF") {
      RDKit::RWMol *mol = RDKit::SLNQueryToMol(slns[i++]);
      TEST_ASSERT(mol);
      delete mol;
    }
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test17() {
  RDKit::RWMol *mol;
  std::string sln;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Test1 " << std::endl;

  // test whitespace at end
  sln = "CH4 \t";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 1);

  delete mol;
  sln = "CH4\t";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(mol);
  TEST_ASSERT(mol->getNumAtoms() == 1);

  delete mol;
  sln = "CH4\tfff";
  mol = RDKit::SLNToMol(sln);
  TEST_ASSERT(!mol);

  delete mol;
  sln = "C[charge=+1] \t";
  mol = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(mol);

  delete mol;
  sln = "C[charge=+1] \tfoo";
  mol = RDKit::SLNQueryToMol(sln);
  TEST_ASSERT(!mol);
  delete mol;
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;
  RDLog::InitLogs();

// FIX: need a test for handling Hs in the SLN itself. This should be done for
// both normal and query SLNs and must be done after the SLN parser handles
// that case (errr, duh)
#if 1
  test1();
  test2();
  test3();
  test4();
  test5();
  test6();
  test7();
  test8();
  test9();
  test10();
  // test11();
  test12();
  test13();
  test14();
  test15();
  test16();
  test17();
#endif
  testIssue277();
  testIssue278();
}
