//
//  Copyright (c) 2007-2018, Novartis Institutes for BioMedical Research Inc.
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
#include <RDGeneral/test.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>
#include <GraphMol/RDKitBase.h>
#include <string>
#include <iostream>
#include <Geometry/point.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/ReactionPickler.h>
#include <GraphMol/ChemReactions/ReactionRunner.h>
#include <GraphMol/ChemReactions/ReactionUtils.h>
#include <GraphMol/ChemReactions/SanitizeRxn.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Atom.h>
#include <fstream>

using namespace RDKit;

void test1Basics() {
  ROMol *mol = nullptr;
  ChemicalReaction rxn;
  MOL_SPTR_VECT reacts;
  std::vector<MOL_SPTR_VECT> prods;
  std::string smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing ChemicalReaction infrastructure"
                       << std::endl;

  smi = "[C:1](=[O:2])O";
  mol = SmartsToMol(smi);
  TEST_ASSERT(mol);
  rxn.addReactantTemplate(ROMOL_SPTR(mol));
  TEST_ASSERT(rxn.getNumReactantTemplates() == 1);

  smi = "[N:3][C:4]";
  mol = SmartsToMol(smi);
  TEST_ASSERT(mol);
  rxn.addReactantTemplate(ROMOL_SPTR(mol));
  TEST_ASSERT(rxn.getNumReactantTemplates() == 2);

  smi = "[C:1](=[O:2])[N:3][C:4]";
  mol = SmartsToMol(smi);
  TEST_ASSERT(mol);
  rxn.addProductTemplate(ROMOL_SPTR(mol));
  TEST_ASSERT(rxn.getNumProductTemplates() == 1);

  smi = "C(=O)O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  rxn.initReactantMatchers();
  prods = rxn.runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);

  reacts.clear();
  smi = "CC(C(=O)O)C(=O)O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  prods = rxn.runReactants(reacts);
  TEST_ASSERT(prods.size() == 2);

  reacts.clear();
  smi = "CC(C(=O)O)C(=O)O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "NCN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  prods = rxn.runReactants(reacts);
  TEST_ASSERT(prods.size() == 4);

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test2SimpleReactions() {
  ROMol *mol = nullptr;
  ChemicalReaction rxn;
  MOL_SPTR_VECT reacts;
  std::vector<MOL_SPTR_VECT> prods;
  std::string smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing simple reactions" << std::endl;

  smi = "[C:1](=[O:2])O";
  mol = SmartsToMol(smi);
  TEST_ASSERT(mol);
  rxn.addReactantTemplate(ROMOL_SPTR(mol));
  TEST_ASSERT(rxn.getNumReactantTemplates() == 1);

  smi = "[N:3][C:4]";
  mol = SmartsToMol(smi);
  TEST_ASSERT(mol);
  rxn.addReactantTemplate(ROMOL_SPTR(mol));
  TEST_ASSERT(rxn.getNumReactantTemplates() == 2);

  smi = "[C:1](=[O:2])[N:3][C:4]";
  mol = SmartsToMol(smi);
  TEST_ASSERT(mol);
  rxn.addProductTemplate(ROMOL_SPTR(mol));
  TEST_ASSERT(rxn.getNumProductTemplates() == 1);

  smi = "C(=O)O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  rxn.initReactantMatchers();
  prods = rxn.runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 4);
  TEST_ASSERT(prods[0][0]->getNumBonds() == 3);

  reacts.clear();
  smi = "CC(=O)O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  prods = rxn.runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
  TEST_ASSERT(prods[0][0]->getNumBonds() == 4);

  reacts.clear();
  smi = "CC(C(=O)O)C(=O)O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  prods = rxn.runReactants(reacts);
  TEST_ASSERT(prods.size() == 2);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[1].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 9);
  TEST_ASSERT(prods[1][0]->getNumAtoms() == 9);
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[1][0]));

  reacts.clear();
  smi = "CC(C(=O)O)C(=O)O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "NCN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));
  prods = rxn.runReactants(reacts);
  TEST_ASSERT(prods.size() == 4);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[1].size() == 1);
  TEST_ASSERT(prods[2].size() == 1);
  TEST_ASSERT(prods[3].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 10);
  TEST_ASSERT(prods[1][0]->getNumAtoms() == 10);
  TEST_ASSERT(prods[2][0]->getNumAtoms() == 10);
  TEST_ASSERT(prods[3][0]->getNumAtoms() == 10);
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[1][0]));
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[2][0]));
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[3][0]));

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test3RingFormation() {
  ROMol *mol = nullptr;
  ChemicalReaction rxn;
  MOL_SPTR_VECT reacts;
  std::vector<MOL_SPTR_VECT> prods;
  std::string smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing ring formation reactions" << std::endl;

  smi = "[C:1]=[C:2]";
  mol = SmartsToMol(smi);
  TEST_ASSERT(mol);
  rxn.addReactantTemplate(ROMOL_SPTR(mol));
  TEST_ASSERT(rxn.getNumReactantTemplates() == 1);

  smi = "[C:3]=[C:4][C:5]=[C:6]";
  mol = SmartsToMol(smi);
  TEST_ASSERT(mol);
  rxn.addReactantTemplate(ROMOL_SPTR(mol));
  TEST_ASSERT(rxn.getNumReactantTemplates() == 2);

  smi = "[C:1]1[C:2][C:3][C:4]=[C:5][C:6]1";
  mol = SmartsToMol(smi);
  TEST_ASSERT(mol);
  rxn.addProductTemplate(ROMOL_SPTR(mol));
  TEST_ASSERT(rxn.getNumProductTemplates() == 1);

  smi = "C=C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "C=CC=C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  rxn.initReactantMatchers();
  prods = rxn.runReactants(reacts);
  TEST_ASSERT(prods.size() == 4);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 6);

  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[1][0]));
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[2][0]));
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[3][0]));

  reacts.clear();
  smi = "CC=C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "C=CC=C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  prods = rxn.runReactants(reacts);
  TEST_ASSERT(prods.size() == 4);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 7);
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[1][0]));
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[2][0]));
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[3][0]));

  reacts.clear();
  smi = "CC=C[Cl]";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "[F]C=CC=C[Br]";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  prods = rxn.runReactants(reacts);
  TEST_ASSERT(prods.size() == 4);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 10);
  for (unsigned int i = 0; i < prods.size(); i++) {
    unsigned int nDifferent = 0;
    for (unsigned int j = 0; j < prods.size(); j++) {
      if (MolToSmiles(*prods[i][0]) != MolToSmiles(*prods[j][0])) {
        nDifferent += 1;
      }
    }
    TEST_ASSERT(nDifferent == 2);
  }

  reacts.clear();
  smi = "C1C=CCCC1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "C=CC=C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  prods = rxn.runReactants(reacts);
  TEST_ASSERT(prods.size() == 4);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 10);
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[1][0]));
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[2][0]));
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[3][0]));

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test4MultipleProducts() {
  ROMol *mol = nullptr;
  ChemicalReaction rxn;
  MOL_SPTR_VECT reacts;
  std::vector<MOL_SPTR_VECT> prods;
  std::string smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing reactions forming multiple products"
                       << std::endl;

  // yes, I know this is bogus... it's a test!
  smi = "[N:1][C:2][C:3](=[O:4])[O:5]";
  mol = SmartsToMol(smi);
  TEST_ASSERT(mol);
  rxn.addReactantTemplate(ROMOL_SPTR(mol));
  TEST_ASSERT(rxn.getNumReactantTemplates() == 1);

  smi = "[N:6][C:7][C:8](=[O:9])[O:10]";
  mol = SmartsToMol(smi);
  TEST_ASSERT(mol);
  rxn.addReactantTemplate(ROMOL_SPTR(mol));
  TEST_ASSERT(rxn.getNumReactantTemplates() == 2);

  smi = "[N:1]1[C:2][C:3](=[O:4])[N:6][C:7][C:8]1=[O:9]";
  mol = SmartsToMol(smi);
  TEST_ASSERT(mol);
  rxn.addProductTemplate(ROMOL_SPTR(mol));
  TEST_ASSERT(rxn.getNumProductTemplates() == 1);

  smi = "[O:5][O:10]";
  mol = SmartsToMol(smi);
  TEST_ASSERT(mol);
  rxn.addProductTemplate(ROMOL_SPTR(mol));
  TEST_ASSERT(rxn.getNumProductTemplates() == 2);

  smi = "OC(=O)CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "OC(=O)CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  rxn.initReactantMatchers();
  prods = rxn.runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 2);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 8);
  TEST_ASSERT(prods[0][1]->getNumAtoms() == 2);
  smi = MolToSmiles(*prods[0][0]);
  std::cerr << "smi: " << smi << std::endl;
  TEST_ASSERT(smi == "O=C1CNC(=O)CN1");
  TEST_ASSERT(MolToSmiles(*prods[0][1]) == "OO");

  reacts.clear();
  smi = "COC(=O)CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "COC(=O)CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  prods = rxn.runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 2);
  MolOps::sanitizeMol(*(static_cast<RWMol *>(prods[0][0].get())));
  MolOps::sanitizeMol(*(static_cast<RWMol *>(prods[0][1].get())));
  std::cerr << "1: " << MolToSmiles(*prods[0][0]) << std::endl;
  std::cerr << "2: " << MolToSmiles(*prods[0][1]) << std::endl;
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 8);
  TEST_ASSERT(prods[0][1]->getNumAtoms() == 4);
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == "O=C1CNC(=O)CN1");
  TEST_ASSERT(MolToSmiles(*prods[0][1]) == "COOC");

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test5Salts() {
  ROMol *mol = nullptr;
  ChemicalReaction rxn;
  MOL_SPTR_VECT reacts;
  std::vector<MOL_SPTR_VECT> prods;
  std::string smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing salt handling" << std::endl;

  smi = "[C:1](=[O:2])O";
  mol = SmartsToMol(smi);
  TEST_ASSERT(mol);
  rxn.addReactantTemplate(ROMOL_SPTR(mol));
  TEST_ASSERT(rxn.getNumReactantTemplates() == 1);

  smi = "[N:3][C:4]";
  mol = SmartsToMol(smi);
  TEST_ASSERT(mol);
  rxn.addReactantTemplate(ROMOL_SPTR(mol));
  TEST_ASSERT(rxn.getNumReactantTemplates() == 2);

  smi = "[C:1](=[O:2])[N:3][C:4]";
  mol = SmartsToMol(smi);
  TEST_ASSERT(mol);
  rxn.addProductTemplate(ROMOL_SPTR(mol));
  TEST_ASSERT(rxn.getNumProductTemplates() == 1);

  smi = "C(=O)O.[ClH]";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "CN.C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  rxn.initReactantMatchers();
  prods = rxn.runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 4);
  TEST_ASSERT(prods[0][0]->getNumBonds() == 3);

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test6DaylightParser() {
  ROMol *mol = nullptr;
  ChemicalReaction *rxn;
  MOL_SPTR_VECT reacts;
  std::vector<MOL_SPTR_VECT> prods;
  std::string smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Daylight parser" << std::endl;

  smi = "[C:1](=[O:2])O.[N:3][C:4]>>[C:1](=[O:2])[N:3][C:4]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);

  smi = "C(=O)O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 4);
  TEST_ASSERT(prods[0][0]->getNumBonds() == 3);

  delete rxn;
  reacts.clear();
  smi =
      "[N:1][C:2][C:3](=[O:4])[O:5].[N:6][C:7][C:8](=[O:9])[O:10]>>[N:1]1[C:2]["
      "C:3](=[O:4])[N:6][C:7][C:8]1=[O:9].[O:5][O:10]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);

  TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
  TEST_ASSERT(rxn->getNumProductTemplates() == 2);

  smi = "OC(=O)CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "OC(=O)CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 2);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 8);
  TEST_ASSERT(prods[0][1]->getNumAtoms() == 2);
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == "O=C1CNC(=O)CN1");
  TEST_ASSERT(MolToSmiles(*prods[0][1]) == "OO");

  delete rxn;
  reacts.clear();
  smi =
      "[N:1][C:2][C:3](=[O:4])[O:5].[N:6][C:7][C:8](=[O:9])[O:10].[N:1]1[C:2]["
      "C:3](=[O:4])[N:6][C:7][C:8]1=[O:9].[O:5][O:10]";
  try {
    rxn = RxnSmartsToChemicalReaction(smi);
  } catch (ChemicalReactionParserException &) {
  }

  smi =
      "[N:1][C:2][C:3](=[O:4])[O:5].[N:6][C:7][C:8](=[O:9])[O:10]>>[N:1]1[C:2]["
      "C:3](=[O:4])[N:6][C:7][C:8]1=[O:9]>>[O:5][O:10]";
  try {
    rxn = RxnSmartsToChemicalReaction(smi);
  } catch (ChemicalReactionParserException &) {
  }

  smi =
      "[Q:1][C:2][C:3](=[O:4])[O:5].[N:6][C:7][C:8](=[O:9])[O:10]>>[N:1]1[C:2]["
      "C:3](=[O:4])[N:6][C:7][C:8]1=[O:9].[O:5][O:10]";
  try {
    rxn = RxnSmartsToChemicalReaction(smi);
  } catch (ChemicalReactionParserException &) {
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test7MDLParser() {
  ROMol *mol = nullptr;
  ChemicalReaction *rxn;
  MOL_SPTR_VECT reacts;
  std::vector<MOL_SPTR_VECT> prods;
  std::string smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing MDL parser" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName;

  fName = rdbase + "/Code/GraphMol/ChemReactions/testData/AmideBond.rxn";

  rxn = RxnFileToChemicalReaction(fName);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);

  smi = "C(=O)O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 4);
  TEST_ASSERT(prods[0][0]->getNumBonds() == 3);

  delete rxn;
  reacts.clear();
  fName = rdbase + "/Code/GraphMol/ChemReactions/testData/cyclization1.rxn";
  rxn = RxnFileToChemicalReaction(fName);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);

  smi = "OC(=O)CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "OC(=O)CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 8);
  std::cerr << MolToSmiles(*prods[0][0]) << std::endl;
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == "O=C1CNC(=O)CN1");

  delete rxn;
  reacts.clear();
  fName = rdbase + "/Code/GraphMol/ChemReactions/testData/cyclization1.rxn";
  rxn = RxnFileToChemicalReaction(fName);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);

  smi = "OC(=O)CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "OC(=O)CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 8);
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == "O=C1CNC(=O)CN1");

  delete rxn;
  reacts.clear();
  fName = rdbase + "/Code/GraphMol/ChemReactions/testData/cyclization2.rxn";
  rxn = RxnFileToChemicalReaction(fName);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);

  smi = "CC=C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "C=CC=C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 4);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 7);
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[1][0]));
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[2][0]));
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[3][0]));

  reacts.clear();
  smi = "CC=C[Cl]";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "[F]C=CC=C[Br]";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 4);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 10);
  for (unsigned int i = 0; i < prods.size(); i++) {
    unsigned int nDifferent = 0;
    for (unsigned int j = 0; j < prods.size(); j++) {
      if (MolToSmiles(*prods[i][0]) != MolToSmiles(*prods[j][0])) {
        nDifferent += 1;
      }
    }
    TEST_ASSERT(nDifferent == 2);
  }

  reacts.clear();
  smi = "C1C=CCCC1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "C=CC=C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 4);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 10);
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[1][0]));
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[2][0]));
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[3][0]));

  delete rxn;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test8Validation() {
  ChemicalReaction *rxn;
  unsigned int nWarn, nError;

  std::string smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Reaction Validation." << std::endl;

  smi = "[C:1](=[O:2])O.[N:3][C:4]>>[C:1](=[O:2])[N:3][C:4]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->validate(nWarn, nError, true));
  TEST_ASSERT(nWarn == 0);
  TEST_ASSERT(nError == 0);
  delete rxn;

  smi = "[C:1](=[O:2])[O:5].[N:3][C:4]>>[C:1](=[O:2])[N:3][C:4]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->validate(nWarn, nError, true));
  TEST_ASSERT(nWarn == 1);
  TEST_ASSERT(nError == 0);
  delete rxn;

  smi = "[C:1](=[O:2])O.[N:1][C:4]>>[C:1](=[O:2])[N:3][C:4]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(!rxn->validate(nWarn, nError, true));
  TEST_ASSERT(nWarn == 1);
  TEST_ASSERT(nError == 1);
  delete rxn;

  smi = "[C:1](=[O:2])O.[N:3][C:4]>>[C:1](=[O:2])[N:3][C:5]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->validate(nWarn, nError, true));
  TEST_ASSERT(nWarn == 2);
  TEST_ASSERT(nError == 0);
  delete rxn;

  smi = "[C:1](=[O:2])O.[N:3][C:4]>>[C:1](=[O:2])[N:3][C:1]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->validate(nWarn, nError, false));
  TEST_ASSERT(nWarn == 2);
  TEST_ASSERT(nError == 0);

  delete rxn;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test9ProductQueries() {
  ROMol *mol = nullptr;
  ChemicalReaction *rxn;
  MOL_SPTR_VECT reacts;
  std::vector<MOL_SPTR_VECT> prods;
  std::string smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing handling of query atoms in the products"
                       << std::endl;

  smi = "[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[*:3]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);

  smi = "C(=O)O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 4);
  TEST_ASSERT(prods[0][0]->getNumBonds() == 3);

  TEST_ASSERT(prods[0][0]->getAtomWithIdx(2)->getAtomicNum() == 7);
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == "CNC=O");

  delete rxn;
  smi = "[C:1](=[O:2])O.[N:3]>>[C:1](=[*:2])[*:3]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);

  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 4);
  TEST_ASSERT(prods[0][0]->getNumBonds() == 3);

  TEST_ASSERT(prods[0][0]->getAtomWithIdx(2)->getAtomicNum() == 7);
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == "CNC=O");

  delete rxn;
  smi = "[C:1](=[O:2])O.[N:3]>>[*:1](=[*:2])[*:3]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);

  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 4);
  TEST_ASSERT(prods[0][0]->getNumBonds() == 3);

  TEST_ASSERT(prods[0][0]->getAtomWithIdx(2)->getAtomicNum() == 7);
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == "CNC=O");

  delete rxn;
  smi =
      "[*:1]1:[*:2]:[*:3]:[*:4]:[*:5]:[*:6]:1>>[*:1]1:[*:2]:[*:3]:[*:4]:[*:5]:["
      "*:6]:1C";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);

  reacts.clear();
  smi = "c1ccccc1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 12);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 7);

  TEST_ASSERT(prods[0][0]->getAtomWithIdx(0)->getIsAromatic());
  TEST_ASSERT(prods[0][0]->getAtomWithIdx(1)->getIsAromatic());
  TEST_ASSERT(prods[0][0]->getBondBetweenAtoms(0, 1)->getIsAromatic());

  delete rxn;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test10ChiralityDaylight() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Chirality handling" << std::endl;

  {  // default behavior, make no changes w.r.t. stereochem
    std::string smi =
        "[F:1][C:2]([Cl:3])([Br:4])[I:5]>>[*:1][C:2]([*:3])([*:4])[*:5]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();

    MOL_SPTR_VECT reacts;
    reacts.clear();
    smi = "F[C@](Cl)(Br)I";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    std::vector<MOL_SPTR_VECT> prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
    TEST_ASSERT(prods[0][0]->getNumBonds() == 4);
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "F[C@](Cl)(Br)I");

    reacts.clear();
    smi = "F[C@@](Cl)(Br)I";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
    TEST_ASSERT(prods[0][0]->getNumBonds() == 4);
    BOOST_LOG(rdInfoLog) << MolToSmiles(*prods[0][0], true) << std::endl;
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "F[C@@](Cl)(Br)I");

    reacts.clear();
    smi = "FC(Cl)(Br)I";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
    TEST_ASSERT(prods[0][0]->getNumBonds() == 4);
    BOOST_LOG(rdInfoLog) << MolToSmiles(*prods[0][0], true) << std::endl;
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "FC(Cl)(Br)I");

    delete rxn;
  }

  {  // a reaction that retains stereochem
    std::string smi =
        "[F:1][C@:2]([Cl:3])([Br:4])[I:5]>>[*:1][C@:2]([*:3])([*:4])[*:5]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();

    MOL_SPTR_VECT reacts;
    reacts.clear();
    smi = "F[C@](Cl)(Br)I";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    std::vector<MOL_SPTR_VECT> prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
    TEST_ASSERT(prods[0][0]->getNumBonds() == 4);
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "F[C@](Cl)(Br)I");

    reacts.clear();
    smi = "F[C@@](Cl)(Br)I";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
    TEST_ASSERT(prods[0][0]->getNumBonds() == 4);
    BOOST_LOG(rdInfoLog) << MolToSmiles(*prods[0][0], true) << std::endl;
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "F[C@@](Cl)(Br)I");

    reacts.clear();
    smi = "FC(Cl)(Br)I";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
    TEST_ASSERT(prods[0][0]->getNumBonds() == 4);
    BOOST_LOG(rdInfoLog) << MolToSmiles(*prods[0][0], true) << std::endl;
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "FC(Cl)(Br)I");

    delete rxn;
  }
  {  // a reaction that inverts stereochem
    std::string smi =
        "[F:1][C@:2]([Cl:3])([Br:4])[I:5]>>[*:1][C@@:2]([*:3])([*:4])[*:5]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();

    MOL_SPTR_VECT reacts;
    reacts.clear();
    smi = "F[C@](Cl)(Br)I";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    std::vector<MOL_SPTR_VECT> prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
    TEST_ASSERT(prods[0][0]->getNumBonds() == 4);
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "F[C@@](Cl)(Br)I");

    reacts.clear();
    smi = "F[C@@](Cl)(Br)I";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
    TEST_ASSERT(prods[0][0]->getNumBonds() == 4);
    BOOST_LOG(rdInfoLog) << MolToSmiles(*prods[0][0], true) << std::endl;
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "F[C@](Cl)(Br)I");

    reacts.clear();
    smi = "FC(Cl)(Br)I";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
    TEST_ASSERT(prods[0][0]->getNumBonds() == 4);
    BOOST_LOG(rdInfoLog) << MolToSmiles(*prods[0][0], true) << std::endl;
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "FC(Cl)(Br)I");

    delete rxn;
  }

  {  // a reaction that induces/sets stereochem
    std::string smi =
        "[F:1][C:2]([Cl:3])([Br:4])[I:5]>>[*:1][C@@:2]([*:3])([*:4])[*:5]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();

    MOL_SPTR_VECT reacts;
    reacts.clear();
    smi = "F[C@](Cl)(Br)I";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    std::vector<MOL_SPTR_VECT> prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
    TEST_ASSERT(prods[0][0]->getNumBonds() == 4);
    BOOST_LOG(rdInfoLog) << MolToSmiles(*prods[0][0], true) << std::endl;
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "F[C@@](Cl)(Br)I");

    reacts.clear();
    smi = "F[C@@](Cl)(Br)I";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
    TEST_ASSERT(prods[0][0]->getNumBonds() == 4);
    BOOST_LOG(rdInfoLog) << MolToSmiles(*prods[0][0], true) << std::endl;
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "F[C@@](Cl)(Br)I");

    reacts.clear();
    smi = "FC(Cl)(Br)I";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
    TEST_ASSERT(prods[0][0]->getNumBonds() == 4);
    BOOST_LOG(rdInfoLog) << MolToSmiles(*prods[0][0], true) << std::endl;
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "F[C@@](Cl)(Br)I");

    delete rxn;
  }
  {  // a reaction that induces/sets stereochem
    std::string smi =
        "[F:1][C:2]([Cl:3])([Br:4])[I:5]>>[*:1][C@:2]([*:3])([*:4])[*:5]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();

    MOL_SPTR_VECT reacts;
    reacts.clear();
    smi = "F[C@](Cl)(Br)I";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    std::vector<MOL_SPTR_VECT> prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
    TEST_ASSERT(prods[0][0]->getNumBonds() == 4);
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "F[C@](Cl)(Br)I");

    reacts.clear();
    smi = "F[C@@](Cl)(Br)I";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
    TEST_ASSERT(prods[0][0]->getNumBonds() == 4);
    BOOST_LOG(rdInfoLog) << MolToSmiles(*prods[0][0], true) << std::endl;
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "F[C@](Cl)(Br)I");

    reacts.clear();
    smi = "FC(Cl)(Br)I";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
    TEST_ASSERT(prods[0][0]->getNumBonds() == 4);
    BOOST_LOG(rdInfoLog) << MolToSmiles(*prods[0][0], true) << std::endl;
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "F[C@](Cl)(Br)I");

    delete rxn;
  }
  {  // a reaction that removes stereochem
    std::string smi =
        "[F:1][C@:2]([Cl:3])([Br:4])[I:5]>>[*:1][C:2]([*:3])([*:4])[*:5]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();

    MOL_SPTR_VECT reacts;
    reacts.clear();
    smi = "F[C@](Cl)(Br)I";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    std::vector<MOL_SPTR_VECT> prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
    TEST_ASSERT(prods[0][0]->getNumBonds() == 4);
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "FC(Cl)(Br)I");

    reacts.clear();
    smi = "F[C@@](Cl)(Br)I";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
    TEST_ASSERT(prods[0][0]->getNumBonds() == 4);
    BOOST_LOG(rdInfoLog) << MolToSmiles(*prods[0][0], true) << std::endl;
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "FC(Cl)(Br)I");

    reacts.clear();
    smi = "FC(Cl)(Br)I";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
    TEST_ASSERT(prods[0][0]->getNumBonds() == 4);
    BOOST_LOG(rdInfoLog) << MolToSmiles(*prods[0][0], true) << std::endl;
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "FC(Cl)(Br)I");

    delete rxn;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test11ChiralityRxn() {
  ROMol *mol = nullptr;
  ChemicalReaction *rxn;
  MOL_SPTR_VECT reacts;
  std::vector<MOL_SPTR_VECT> prods;
  std::string smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Chirality handling (rxn)" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName;

  fName = rdbase + "/Code/GraphMol/ChemReactions/testData/SN2_FlagProd.rxn";

  rxn = RxnFileToChemicalReaction(fName);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
  TEST_ASSERT(rxn->getNumProductTemplates() == 2);

  reacts.clear();
  smi = "F[C@@](Cl)(Br)I";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));
  smi = "[OH-]";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 2);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
  TEST_ASSERT(prods[0][0]->getNumBonds() == 4);
  TEST_ASSERT(prods[0][1]->getNumAtoms() == 1);
  TEST_ASSERT(prods[0][1]->getNumBonds() == 0);
  TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "O[C@](F)(Cl)I");
  TEST_ASSERT(MolToSmiles(*prods[0][1], true) == "[Br-]");

  fName =
      rdbase + "/Code/GraphMol/ChemReactions/testData/SN2_FlagProd.retain.rxn";

  delete rxn;
  rxn = RxnFileToChemicalReaction(fName);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
  TEST_ASSERT(rxn->getNumProductTemplates() == 2);
  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 2);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
  TEST_ASSERT(prods[0][0]->getNumBonds() == 4);
  TEST_ASSERT(prods[0][1]->getNumAtoms() == 1);
  TEST_ASSERT(prods[0][1]->getNumBonds() == 0);
  TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "O[C@@](F)(Cl)I");
  TEST_ASSERT(MolToSmiles(*prods[0][1], true) == "[Br-]");

  // this has the chirality flag in the reactants, which is ignored, so the
  // chirality is preserved:
  fName = rdbase + "/Code/GraphMol/ChemReactions/testData/SN2_FlagReact.rxn";
  delete rxn;
  rxn = RxnFileToChemicalReaction(fName);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
  TEST_ASSERT(rxn->getNumProductTemplates() == 2);
  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 2);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
  TEST_ASSERT(prods[0][0]->getNumBonds() == 4);
  TEST_ASSERT(prods[0][1]->getNumAtoms() == 1);
  TEST_ASSERT(prods[0][1]->getNumBonds() == 0);
  TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "O[C@@](F)(Cl)I");
  TEST_ASSERT(MolToSmiles(*prods[0][1], true) == "[Br-]");

  delete rxn;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test12DoubleBondStereochem() {
  ROMol *mol = nullptr;
  ChemicalReaction *rxn;
  MOL_SPTR_VECT reacts;
  std::vector<MOL_SPTR_VECT> prods;
  std::string smi, stereo;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing handling of double bond stereochemistry"
                       << std::endl;

  smi = "[C:1](=[O:2])-[O;H0]>>[C:1](=[O:2])[X]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);

  reacts.clear();
  smi = "COC(=O)/C=C/Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));
  MolOps::assignStereochemistry(*mol);

  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 6);
  TEST_ASSERT(prods[0][0]->getNumBonds() == 5);
  prods[0][0]->updatePropertyCache();
  BOOST_LOG(rdInfoLog) << MolToSmiles(*prods[0][0], true) << std::endl;
  TEST_ASSERT(mol->getBondWithIdx(4)->getStereo() == Bond::STEREOE);

  Bond *stereoBond = prods[0][0]->getBondWithIdx(3);
  INT_VECT stereoAtomsRef{{0, 5}};
  TEST_ASSERT(stereoBond->getStereoAtoms() == stereoAtomsRef);
  TEST_ASSERT(stereoBond->getStereo() == Bond::STEREOTRANS);

  reacts.clear();

  // FIX: note that the handling of this situation, where we match double bonds
  // that have stereochem indicated, is currently not completely correct since
  // the stereochem querying stuff doesn't match C/C=C\Cl when C\C=C/Cl is
  // provided as a query.
  delete rxn;
  smi = "[Cl:3]\\[C:1]=[C:2]/[C:4]>>[Cl:3]\\[C:1]=[C:2]\\[C:4]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);

  reacts.clear();
  smi = "Cl\\C=C/C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));
  MolOps::assignStereochemistry(*mol);

  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 4);
  TEST_ASSERT(prods[0][0]->getNumBonds() == 3);
  prods[0][0]->updatePropertyCache();
  BOOST_LOG(rdInfoLog) << MolToSmiles(*prods[0][0], true) << std::endl;
  TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);

  stereoBond = prods[0][0]->getBondWithIdx(1);
  stereoAtomsRef = {0, 3};
  TEST_ASSERT(stereoBond->getStereoAtoms() == stereoAtomsRef);
  TEST_ASSERT(stereoBond->getStereo() == Bond::STEREOTRANS);

  reacts.clear();
  delete rxn;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test13Issue1748846() {
  ROMol *mol = nullptr;
  ChemicalReaction *rxn;
  MOL_SPTR_VECT reacts;
  std::vector<MOL_SPTR_VECT> prods;
  std::string smi, stereo;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing sf.net Issue 1748846: bad bond orders in reaction products"
      << std::endl;

  smi = "c1ccccc1[C:1].[*:2][At]>>c1ccccc1[C:1][*:2]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);

  reacts.clear();
  smi = "c1ccccc1C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "[At]OC";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() > 0);
  BOOST_LOG(rdInfoLog) << prods[0].size() << std::endl;
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 9);
  BOOST_LOG(rdInfoLog) << MolToSmiles(*prods[0][0], true) << std::endl;
  TEST_ASSERT(prods[0][0]->getBondBetweenAtoms(0, 1)->getBondType() ==
              Bond::AROMATIC);
  TEST_ASSERT(prods[0][0]->getBondBetweenAtoms(1, 2)->getBondType() ==
              Bond::AROMATIC);
  smi = MolToSmiles(*prods[0][0], true);
  TEST_ASSERT(smi == "COCc1ccccc1");

  delete rxn;
  smi = "[c:3][C:1].[*:2][At]>>[c:3][C:1][*:2]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);

  reacts.clear();
  smi = "c1ccccc1C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "[At]OC";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() > 0);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 9);
  smi = MolToSmiles(*prods[0][0], true);
  TEST_ASSERT(smi == "COCc1ccccc1");

  delete rxn;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test14Issue1804420() {
  ROMol *mol = nullptr;
  ChemicalReaction *rxn;
  MOL_SPTR_VECT reacts;
  std::vector<MOL_SPTR_VECT> prods;
  std::string smi, stereo;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing sf.net Issue 1804420: bad handling of query "
                          "atoms in products."
                       << std::endl;

  // NOTE that this bug was actually in the smarts parser, so this is really
  // just another reaction test... still, more tests are better

  smi = "[N;D3;R:1]-!@[*:2]>>[At][N:1].[*:2][At]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
  TEST_ASSERT(rxn->getNumProductTemplates() == 2);

  TEST_ASSERT((*rxn->beginReactantTemplates())
                  ->getAtomWithIdx(0)
                  ->hasProp(common_properties::molAtomMapNumber));

  reacts.clear();
  smi = "C1CCN1CCC";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 2);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
  smi = MolToSmiles(*prods[0][0], true);
  TEST_ASSERT(smi == "[At]N1CCC1");
  TEST_ASSERT(prods[0][1]->getNumAtoms() == 4);
  smi = MolToSmiles(*prods[0][1], true);
  TEST_ASSERT(smi == "CCC[At]");

  delete rxn;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test15Issue1882749() {
  ROMol *mol = nullptr;
  ChemicalReaction *rxn;
  MOL_SPTR_VECT reacts;
  unsigned int nWarn, nError;
  std::vector<MOL_SPTR_VECT> prods;
  std::string smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing sf.net Issue 1882749: property handling in products."
      << std::endl;

  smi = "[N:1]-!@*>>[N;+1,+0:1][#0]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);
  TEST_ASSERT(rxn->validate(nWarn, nError, false));
  TEST_ASSERT(nWarn == 1);
  TEST_ASSERT(nError == 0);

  delete rxn;
  smi = "[N:1]-!@*>>[N;-1:1]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);
  TEST_ASSERT(rxn->validate(nWarn, nError, false));
  TEST_ASSERT(nWarn == 0);
  TEST_ASSERT(nError == 0);

  reacts.clear();
  smi = "C1CCN1CCC";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 4);
  TEST_ASSERT(prods[0][0]->getAtomWithIdx(0)->getFormalCharge() == -1);

  delete rxn;
  smi = "[N;D3:1]-!@*>>[N;H1:1]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);
  TEST_ASSERT(rxn->validate(nWarn, nError, false));
  TEST_ASSERT(nWarn == 0);
  TEST_ASSERT(nError == 0);
  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 4);
  TEST_ASSERT(prods[0][0]->getAtomWithIdx(0)->getNumExplicitHs() == 1);

  delete rxn;
  smi = "[N;D3:1]-!@*>>[15N:1]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);
  TEST_ASSERT(rxn->validate(nWarn, nError, false));
  TEST_ASSERT(nWarn == 0);
  TEST_ASSERT(nError == 0);
  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 4);
  TEST_ASSERT(prods[0][0]->getAtomWithIdx(0)->getAtomicNum() == 7);
  TEST_ASSERT(prods[0][0]->getAtomWithIdx(0)->getIsotope() == 15);
  std::cerr << " mass: " << prods[0][0]->getAtomWithIdx(0)->getMass()
            << std::endl;
  TEST_ASSERT(feq(prods[0][0]->getAtomWithIdx(0)->getMass(), 15.0001));

  delete rxn;
  smi = "[N;D3:1]-!@*>>[15N;-1:1]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);
  TEST_ASSERT(rxn->validate(nWarn, nError, false));
  TEST_ASSERT(nWarn == 0);
  TEST_ASSERT(nError == 0);
  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 4);
  TEST_ASSERT(prods[0][0]->getAtomWithIdx(0)->getFormalCharge() == -1);
  TEST_ASSERT(prods[0][0]->getAtomWithIdx(0)->getIsotope() == 15);
  TEST_ASSERT(feq(prods[0][0]->getAtomWithIdx(0)->getMass(), 15.0001));

  reacts.clear();
  smi = "CS(=O)C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  delete rxn;
  smi = "[S:1]=[O:2]>>[S;+2:1]-[O;-:2]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);
  TEST_ASSERT(rxn->validate(nWarn, nError, false));
  TEST_ASSERT(nWarn == 0);
  TEST_ASSERT(nError == 0);
  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 4);
  TEST_ASSERT(prods[0][0]->getAtomWithIdx(0)->getFormalCharge() == +2);
  TEST_ASSERT(prods[0][0]->getAtomWithIdx(1)->getFormalCharge() == -1);

  reacts.clear();
  smi = "CO";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  delete rxn;
  smi = "[O:1]>>[O:1][13C]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);
  TEST_ASSERT(rxn->validate(nWarn, nError, false));
  TEST_ASSERT(nWarn == 0);
  TEST_ASSERT(nError == 0);
  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 3);
  TEST_ASSERT(prods[0][0]->getAtomWithIdx(1)->getIsotope() == 13);
  TEST_ASSERT(feq(prods[0][0]->getAtomWithIdx(1)->getMass(), 13.00335));

  reacts.clear();
  smi = "CO";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  delete rxn;
  smi = "[O:1]>>[O:1][3#0]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);
  TEST_ASSERT(rxn->validate(nWarn, nError, false));
  TEST_ASSERT(nWarn == 0);
  TEST_ASSERT(nError == 0);
  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 3);

  MolOps::sanitizeMol(*(static_cast<RWMol *>(prods[0][0].get())));
  TEST_ASSERT(prods[0][0]->getAtomWithIdx(1)->getIsotope() == 3);
  TEST_ASSERT(feq(prods[0][0]->getAtomWithIdx(1)->getMass(), 0.000));
  TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "[3*]OC");

  delete rxn;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test16Exceptions() {
  ChemicalReaction *rxn;
  std::string rxnB;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing parser exception handling" << std::endl;

  rxnB =
      "$RXN\n"
      "\n"
      "      ISIS     082120061354\n"
      "\n"
      "  2  1\n"
      "$MOL\n"
      "\n"
      "  -ISIS-  08210613542D\n"
      "\n"
      "  4  2  0  0  0  0  0  0  0  0999 V2000\n"
      "   -1.4340   -0.6042    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0\n"
      "   -0.8639   -0.9333    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -1.4340    0.0542    0.0000 O   0  0  0  0  0  0  0  0  0  1  0  0\n"
      "  1  2  1  0  0  0  0\n"
      "  1  3  2  0  0  0  0\n"
      "M  END\n"
      "$MOL\n"
      "\n"
      "  -ISIS-  08210613542D\n"
      "\n"
      "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
      "    2.2125   -0.7833    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0\n"
      "M  END\n"
      "$MOL\n"
      "\n"
      "  -ISIS-  08210613542D\n"
      "\n"
      "  3  2  0  0  0  0  0  0  0  0999 V2000\n"
      "    9.5282   -0.8083    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0\n"
      "    8.9579   -0.4792    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0\n"
      "    8.9579    0.1792    0.0000 O   0  0  0  0  0  0  0  0  0  1  0  0\n"
      "  1  2  1  0  0  0  0\n"
      "  2  3  2  0  0  0  0\n"
      "M  END\n";

  rxn = (ChemicalReaction *)0x1;
  try {
    rxn = RxnBlockToChemicalReaction(rxnB);
  } catch (ChemicalReactionParserException &) {
    rxn = (ChemicalReaction *)nullptr;
  }
  TEST_ASSERT(!rxn);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test17Issue1920627() {
  ROMol *mol = nullptr;
  ChemicalReaction *rxn;
  MOL_SPTR_VECT reacts;
  std::vector<MOL_SPTR_VECT> prods;
  ROMOL_SPTR prod;
  std::string smi, cip;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing sf.net Issue 1920627: chirality flip in reactions."
      << std::endl;

  smi = "[C:1](=[O:2])>>[C:1](=[S:2])";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);

#if 1
  reacts.clear();
  smi = "C[C@](Cl)(CO)CC(=O)NC";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  reacts.push_back(ROMOL_SPTR(mol));
  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  prod = prods[0][0];
  MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));
  MolOps::assignStereochemistry(*prod);
  TEST_ASSERT(prod->getNumAtoms() == 10);
  TEST_ASSERT(prod->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));
  prod->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  reacts.clear();
  smi = "C[C@H](CO)CC(=O)NC";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

  reacts.push_back(ROMOL_SPTR(mol));
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  prod = prods[0][0];
  MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));
  MolOps::assignStereochemistry(*prod);
  TEST_ASSERT(prod->getNumAtoms() == 9);
  TEST_ASSERT(prod->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));
  prod->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

  // make sure the above two tests work "backwards":
  reacts.clear();
  smi = "C[C@@](Cl)(CO)CC(=O)NC";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

  reacts.push_back(ROMOL_SPTR(mol));
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  prod = prods[0][0];
  MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));
  MolOps::assignStereochemistry(*prod);
  TEST_ASSERT(prod->getNumAtoms() == 10);
  TEST_ASSERT(prod->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));
  prod->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

  reacts.clear();
  smi = "C[C@@H](CO)CC(=O)NC";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  reacts.push_back(ROMOL_SPTR(mol));
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  prod = prods[0][0];
  MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));
  MolOps::assignStereochemistry(*prod);
  TEST_ASSERT(prod->getNumAtoms() == 9);
  TEST_ASSERT(prod->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));
  prod->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  // do some "restructuring" and make sure things still work:
  reacts.clear();
  smi = "[C@@](C)(Cl)(CO)CC(=O)NC";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

  reacts.push_back(ROMOL_SPTR(mol));
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  prod = prods[0][0];
  MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));
  MolOps::assignStereochemistry(*prod);
  TEST_ASSERT(prod->getNumAtoms() == 10);
  TEST_ASSERT(prod->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));
  prod->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

  reacts.clear();
  smi = "[C@H](C)(CO)CC(=O)NC";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  reacts.push_back(ROMOL_SPTR(mol));
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  prod = prods[0][0];
  MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));
  MolOps::assignStereochemistry(*prod);
  TEST_ASSERT(prod->getNumAtoms() == 9);
  TEST_ASSERT(prod->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));
  prod->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
#endif

  reacts.clear();
  smi = "C(=O)N[C@@H](CC)C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  reacts.push_back(ROMOL_SPTR(mol));
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  prod = prods[0][0];
  MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));
  MolOps::assignStereochemistry(*prod);
  TEST_ASSERT(prod->getNumAtoms() == 7);
  TEST_ASSERT(prod->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
  prod->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  reacts.clear();
  smi = "C(=O)N[C@@H](CC)C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  reacts.push_back(ROMOL_SPTR(mol));
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  prod = prods[0][0];
  MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));
  MolOps::assignStereochemistry(*prod);
  TEST_ASSERT(prod->getNumAtoms() == 7);
  TEST_ASSERT(prod->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
  prod->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  delete rxn;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test18PropertyTransfer() {
  ROMol *mol = nullptr;
  ChemicalReaction *rxn;
  MOL_SPTR_VECT reacts;
  std::vector<MOL_SPTR_VECT> prods;
  ROMOL_SPTR prod;
  std::string smi, cip;
  unsigned int nWarn, nError;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing property transfer in reactions."
                       << std::endl;

  // ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ -----
  // start with sf.net Issue 1934052: propagation of isotope information
  smi = "[C:1]>>[C:1]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);
  TEST_ASSERT(rxn->validate(nWarn, nError, false));
  TEST_ASSERT(nWarn == 0);
  TEST_ASSERT(nError == 0);

  reacts.clear();
  smi = "C";
  mol = SmilesToMol(smi);
  reacts.push_back(ROMOL_SPTR(mol));
  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  prod = prods[0][0];
  TEST_ASSERT(prod->getNumAtoms() == 1);
  TEST_ASSERT(feq(prod->getAtomWithIdx(0)->getMass(), 12, .1));

  reacts.clear();
  smi = "[13CH4]";
  mol = SmilesToMol(smi);
  reacts.push_back(ROMOL_SPTR(mol));
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  prod = prods[0][0];
  TEST_ASSERT(prod->getNumAtoms() == 1);
  TEST_ASSERT(feq(prod->getAtomWithIdx(0)->getMass(), 13, .1));

  delete rxn;
  smi = "[12C:1]>>[13C:1]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);
  TEST_ASSERT(rxn->validate(nWarn, nError, false));
  TEST_ASSERT(nWarn == 0);
  TEST_ASSERT(nError == 0);

  reacts.clear();
  smi = "[12CH4]";
  mol = SmilesToMol(smi);
  reacts.push_back(ROMOL_SPTR(mol));
  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  prod = prods[0][0];
  TEST_ASSERT(prod->getNumAtoms() == 1);
  TEST_ASSERT(prod->getAtomWithIdx(0)->getIsotope() == 13);
  TEST_ASSERT(feq(prod->getAtomWithIdx(0)->getMass(), 13.003, .001));

  reacts.clear();
  smi = "[13CH4]";
  mol = SmilesToMol(smi);
  reacts.push_back(ROMOL_SPTR(mol));
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 0);

  delete rxn;
  smi = "[13C:1]>>[12C:1]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);
  TEST_ASSERT(rxn->validate(nWarn, nError, false));
  TEST_ASSERT(nWarn == 0);
  TEST_ASSERT(nError == 0);

  reacts.clear();
  smi = "[13CH4]";
  mol = SmilesToMol(smi);
  reacts.push_back(ROMOL_SPTR(mol));
  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  prod = prods[0][0];
  TEST_ASSERT(prod->getNumAtoms() == 1);
  TEST_ASSERT(prod->getAtomWithIdx(0)->getIsotope() == 12);
  TEST_ASSERT(feq(prod->getAtomWithIdx(0)->getMass(), 12, .001));

  reacts.clear();
  smi = "C";
  mol = SmilesToMol(smi);
  reacts.push_back(ROMOL_SPTR(mol));
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 0);

  delete rxn;
  smi = "[C:1]>>[12C:1]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);
  TEST_ASSERT(rxn->validate(nWarn, nError, false));
  TEST_ASSERT(nWarn == 0);
  TEST_ASSERT(nError == 0);

  reacts.clear();
  smi = "C";
  mol = SmilesToMol(smi);
  reacts.push_back(ROMOL_SPTR(mol));
  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  prod = prods[0][0];
  TEST_ASSERT(prod->getNumAtoms() == 1);
  TEST_ASSERT(prod->getAtomWithIdx(0)->getIsotope() == 12);
  TEST_ASSERT(feq(prod->getAtomWithIdx(0)->getMass(), 12, .001));

  reacts.clear();
  smi = "[13CH4]";
  mol = SmilesToMol(smi);
  reacts.push_back(ROMOL_SPTR(mol));
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  prod = prods[0][0];
  TEST_ASSERT(prod->getNumAtoms() == 1);
  TEST_ASSERT(prod->getAtomWithIdx(0)->getIsotope() == 12);
  TEST_ASSERT(feq(prod->getAtomWithIdx(0)->getMass(), 12, .001));

  delete rxn;
  smi = "[C:1](=[O:2])>>[C:1](=[S:2])";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);
  TEST_ASSERT(rxn->validate(nWarn, nError, false));
  TEST_ASSERT(nWarn == 0);
  TEST_ASSERT(nError == 0);

  reacts.clear();
  smi = "C=O";
  mol = SmilesToMol(smi);
  reacts.push_back(ROMOL_SPTR(mol));
  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  prod = prods[0][0];
  TEST_ASSERT(prod->getNumAtoms() == 2);
  TEST_ASSERT(feq(prod->getAtomWithIdx(0)->getMass(), 12, .1));
  TEST_ASSERT(feq(prod->getAtomWithIdx(1)->getMass(), 32, .1));

  // ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ -----
  // now look at some other properties
  delete rxn;
  smi = "[c;H1:1][n:2]>>[c:1](O)[n:2]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);
  TEST_ASSERT(rxn->validate(nWarn, nError, false));
  TEST_ASSERT(nWarn == 0);
  TEST_ASSERT(nError == 0);

  reacts.clear();
  smi = "c1ccccn1";
  mol = SmilesToMol(smi);
  reacts.push_back(ROMOL_SPTR(mol));
  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 2);
  TEST_ASSERT(prods[0].size() == 1);
  prod = prods[0][0];
  TEST_ASSERT(prod->getNumAtoms() == 7);
  smi = MolToSmiles(*prod);
  std::cerr << smi << std::endl;
  TEST_ASSERT(smi == "Oc1ccccn1");

  reacts.clear();
  smi = "c1ccc[nH]1";
  mol = SmilesToMol(smi);
  reacts.push_back(ROMOL_SPTR(mol));
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 2);
  TEST_ASSERT(prods[0].size() == 1);
  prod = prods[0][0];
  TEST_ASSERT(prod->getNumAtoms() == 6);
  smi = MolToSmiles(*prod);
  std::cerr << smi << std::endl;
  TEST_ASSERT(smi == "Oc1ccc[nH]1");

  delete rxn;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test19Issue2050085() {
  ROMol *mol = nullptr;
  ChemicalReaction *rxn;
  MOL_SPTR_VECT reacts;
  std::vector<MOL_SPTR_VECT> prods;
  ROMOL_SPTR prod;
  std::string smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing sf.net Issue 2050085: recap failing on chiral molecules."
      << std::endl;

  smi = "[N;!D1;+0](-!@[*:1])-!@[*:2]>>*[*:1].[*:2]*";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
  TEST_ASSERT(rxn->getNumProductTemplates() == 2);

  reacts.clear();
  smi = "C1CC1N[C@H](C)F";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  MolOps::assignStereochemistry(*mol);

  reacts.push_back(ROMOL_SPTR(mol));
  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 2);
  TEST_ASSERT(prods[0].size() == 2);
  MolOps::sanitizeMol(*(static_cast<RWMol *>(prods[0][0].get())));
  MolOps::sanitizeMol(*(static_cast<RWMol *>(prods[0][1].get())));

  prod = prods[0][0];
  TEST_ASSERT(prod->getNumAtoms() == 4);
  prod = prods[0][1];
  TEST_ASSERT(prod->getNumAtoms() == 4);

  delete rxn;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test20BondQueriesInProduct() {
  ROMol *mol = nullptr;
  ChemicalReaction *rxn;
  MOL_SPTR_VECT reacts;
  std::vector<MOL_SPTR_VECT> prods;
  ROMOL_SPTR prod;
  std::string smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing bond queries in the product." << std::endl;

  smi = "[O:1]~[C:2]>>[O:1]~[C:2]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);

  reacts.clear();
  smi = "OCC";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);

  reacts.push_back(ROMOL_SPTR(mol));
  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  MolOps::sanitizeMol(*(static_cast<RWMol *>(prods[0][0].get())));

  prod = prods[0][0];
  TEST_ASSERT(prod->getNumAtoms() == 3);
  TEST_ASSERT(prod->getBondBetweenAtoms(0, 1)->getBondType() == Bond::SINGLE);

  reacts.clear();
  smi = "O=CC";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);

  reacts.push_back(ROMOL_SPTR(mol));
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  MolOps::sanitizeMol(*(static_cast<RWMol *>(prods[0][0].get())));

  prod = prods[0][0];
  TEST_ASSERT(prod->getNumAtoms() == 3);
  TEST_ASSERT(prod->getBondBetweenAtoms(0, 1)->getBondType() == Bond::DOUBLE);
  delete rxn;

  smi = "[O:1]~[C:2]>>[O:1][C:2]";
  rxn = RxnSmartsToChemicalReaction(smi);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);

  reacts.clear();
  smi = "OCC";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);

  reacts.push_back(ROMOL_SPTR(mol));
  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  MolOps::sanitizeMol(*(static_cast<RWMol *>(prods[0][0].get())));

  prod = prods[0][0];
  TEST_ASSERT(prod->getNumAtoms() == 3);
  TEST_ASSERT(prod->getBondBetweenAtoms(0, 1)->getBondType() == Bond::SINGLE);

  reacts.clear();
  smi = "O=CC";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);

  reacts.push_back(ROMOL_SPTR(mol));
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  MolOps::sanitizeMol(*(static_cast<RWMol *>(prods[0][0].get())));

  prod = prods[0][0];
  TEST_ASSERT(prod->getNumAtoms() == 3);
  TEST_ASSERT(prod->getBondBetweenAtoms(0, 1)->getBondType() == Bond::SINGLE);
  delete rxn;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test21Issue2540021() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing sf.net Issue 2540021: bad handling of atoms with explicit Hs."
      << std::endl;

  {
    std::string smi = "[#7;!H0:1]>>[#7:1]C";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);

    MOL_SPTR_VECT reacts;
    reacts.clear();
    smi = "C1=CNC=C1";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    rxn->initReactantMatchers();

    std::vector<MOL_SPTR_VECT> prods;
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);

    ROMOL_SPTR prod = prods[0][0];
    MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));
    TEST_ASSERT(prod->getNumAtoms() == 6);
    TEST_ASSERT(prod->getAtomWithIdx(0)->getAtomicNum() == 7);
    TEST_ASSERT(prod->getAtomWithIdx(0)->getImplicitValence() == 0);
    TEST_ASSERT(prod->getAtomWithIdx(0)->getExplicitValence() == 3);
    TEST_ASSERT(prod->getAtomWithIdx(0)->getNoImplicit() == false);

    delete rxn;
  }

  {
    std::string smi =
        "[c:1]1[c:2][n;H1:3][c:4][n:5]1>>[c:1]1[c:2][n:3][c:4](C)[n:5]1";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);

    MOL_SPTR_VECT reacts;
    reacts.clear();
    smi = "C1=CNC=N1";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    rxn->initReactantMatchers();

    std::vector<MOL_SPTR_VECT> prods;
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);

    ROMOL_SPTR prod = prods[0][0];
    MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));
    TEST_ASSERT(prod->getNumAtoms() == 6);
    TEST_ASSERT(prod->getAtomWithIdx(2)->getAtomicNum() == 7);
    TEST_ASSERT(prod->getAtomWithIdx(2)->getNumExplicitHs() == 1);
    TEST_ASSERT(prod->getAtomWithIdx(5)->getAtomicNum() == 7);
    TEST_ASSERT(prod->getAtomWithIdx(5)->getNumExplicitHs() == 0);

    delete rxn;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test22DotsToRemoveBonds() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing using dots in the products to remove bonds."
                       << std::endl;

  {
    /* 08/05/14
     * This test is changed due to a new behavior of the smarts reaction parser
     * which now
     * allows using parenthesis in products as well
     * original smiles:
     * std::string smi  = "[C:1]1[O:2][N:3]1>>[C:1]1[O:2].[N:3]1"; */
    std::string smi = "[C:1]1[O:2][N:3]1>>([C:1]1[O:2].[N:3]1)";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();

    MOL_SPTR_VECT reacts;
    reacts.clear();
    smi = "C1ON1";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));

    std::vector<MOL_SPTR_VECT> prods;
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);

    ROMOL_SPTR prod = prods[0][0];
    MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));
    TEST_ASSERT(prod->getNumAtoms() == 3);
    TEST_ASSERT(prod->getAtomWithIdx(0)->getAtomicNum() == 6);
    TEST_ASSERT(prod->getAtomWithIdx(1)->getAtomicNum() == 8);
    TEST_ASSERT(prod->getAtomWithIdx(2)->getAtomicNum() == 7);
    TEST_ASSERT(prod->getBondBetweenAtoms(0, 1));
    TEST_ASSERT(prod->getBondBetweenAtoms(0, 2));
    TEST_ASSERT(!prod->getBondBetweenAtoms(1, 2));

    delete rxn;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test23Pickling() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing pickling and depickling reactions."
                       << std::endl;

  {
    std::string smi;
    ChemicalReaction *rxn;
    MOL_SPTR_VECT reacts;
    std::vector<MOL_SPTR_VECT> prods;

    /* 08/05/14
     * This test is changed due to a new behavior of the smarts reaction parser
     * which now allows using parenthesis in products as well original smiles:
     * std::string smi  = "[C:1]1[O:2][N:3]1>>[C:1]1[O:2].[N:3]1"; */
    smi = "[C:1]1[O:2][N:3]1>>([C:1]1[O:2].[N:3]1)";
    rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();

    std::string pkl;
    ReactionPickler::pickleReaction(rxn, pkl);
    delete rxn;
    rxn = new ChemicalReaction();
    ReactionPickler::reactionFromPickle(pkl, rxn);
    rxn->initReactantMatchers();

    reacts.clear();
    smi = "C1ON1";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));

    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);

    ROMOL_SPTR prod = prods[0][0];
    MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));
    TEST_ASSERT(prod->getNumAtoms() == 3);
    TEST_ASSERT(prod->getAtomWithIdx(0)->getAtomicNum() == 6);
    TEST_ASSERT(prod->getAtomWithIdx(1)->getAtomicNum() == 8);
    TEST_ASSERT(prod->getAtomWithIdx(2)->getAtomicNum() == 7);
    TEST_ASSERT(prod->getBondBetweenAtoms(0, 1));
    TEST_ASSERT(prod->getBondBetweenAtoms(0, 2));
    TEST_ASSERT(!prod->getBondBetweenAtoms(1, 2));

    // test construction from a pickle:
    delete rxn;
    rxn = new ChemicalReaction(pkl);
    rxn->initReactantMatchers();
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    prod = prods[0][0];
    MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));
    TEST_ASSERT(prod->getNumAtoms() == 3);
    TEST_ASSERT(prod->getAtomWithIdx(0)->getAtomicNum() == 6);
    TEST_ASSERT(prod->getAtomWithIdx(1)->getAtomicNum() == 8);
    TEST_ASSERT(prod->getAtomWithIdx(2)->getAtomicNum() == 7);
    TEST_ASSERT(prod->getBondBetweenAtoms(0, 1));
    TEST_ASSERT(prod->getBondBetweenAtoms(0, 2));
    TEST_ASSERT(!prod->getBondBetweenAtoms(1, 2));

    delete rxn;
  }

  {
    std::string rdbase = getenv("RDBASE");
    std::string fName;
    std::string smi;
    ChemicalReaction *rxn;
    MOL_SPTR_VECT reacts;
    std::vector<MOL_SPTR_VECT> prods;

    fName = rdbase + "/Code/GraphMol/ChemReactions/testData/cyclization1.rxn";
    rxn = RxnFileToChemicalReaction(fName);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);

    std::string pkl;
    ReactionPickler::pickleReaction(rxn, pkl);
    delete rxn;
    rxn = new ChemicalReaction();
    ReactionPickler::reactionFromPickle(pkl, rxn);
    TEST_ASSERT(!rxn->isInitialized());
    rxn->initReactantMatchers();
    TEST_ASSERT(rxn->isInitialized());

    // quick test of github issue #249
    ReactionPickler::pickleReaction(rxn, pkl);
    delete rxn;
    rxn = new ChemicalReaction();
    ReactionPickler::reactionFromPickle(pkl, rxn);
    TEST_ASSERT(rxn->isInitialized());

    smi = "OC(=O)CN";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));

    smi = "OC(=O)CN";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));

    rxn->initReactantMatchers();
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 8);
    TEST_ASSERT(MolToSmiles(*prods[0][0]) == "O=C1CNC(=O)CN1");
    delete rxn;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test24AtomFlags() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing preservation of atom flags from rxn files."
                       << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName;
  unsigned int nWarn, nError;

  {
    fName = rdbase + "/Code/GraphMol/ChemReactions/testData/atomflags.rxn";
    ChemicalReaction *rxn = RxnFileToChemicalReaction(fName);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);
    TEST_ASSERT(
        rxn->beginReactantTemplates()->get()->getAtomWithIdx(0)->hasProp(
            common_properties::molAtomMapNumber));
    TEST_ASSERT((++(rxn->beginReactantTemplates()))
                    ->get()
                    ->getAtomWithIdx(0)
                    ->hasProp(common_properties::molAtomMapNumber));
    TEST_ASSERT((++(rxn->beginReactantTemplates()))
                    ->get()
                    ->getAtomWithIdx(1)
                    ->hasProp(common_properties::molAtomMapNumber));

    delete rxn;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test25Conformers() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing transfer of conformer data from reactants->products."
      << std::endl;

  unsigned int nWarn, nError;

  {
    std::string smi;
    smi = "[C:1]=[O:2].[N:3]>>[O:2]=[C:1][N:3]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);

    ROMol *mol;
    Conformer *conf;
    MOL_SPTR_VECT reacts;
    std::vector<MOL_SPTR_VECT> prods;

    reacts.clear();
    smi = "C(=O)C";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    conf = new Conformer(3);
    conf->setAtomPos(0, RDGeom::Point3D(1, 0, 0));
    conf->setAtomPos(1, RDGeom::Point3D(1, 1, 0));
    conf->setAtomPos(2, RDGeom::Point3D(1, 0, 1));
    mol->addConformer(conf, true);
    reacts.push_back(ROMOL_SPTR(mol));

    smi = "CN";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    conf = new Conformer(2);
    conf->setAtomPos(0, RDGeom::Point3D(-1, 0, 0));
    conf->setAtomPos(1, RDGeom::Point3D(-1, 1, 0));
    mol->addConformer(conf, true);
    reacts.push_back(ROMOL_SPTR(mol));

    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
    TEST_ASSERT(prods[0][0]->getNumConformers() == 1);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(0).x == 1.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(0).y == 1.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(0).z == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(1).x == 1.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(1).y == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(1).z == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(2).x == -1.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(2).y == 1.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(2).z == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(3).x == 1.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(3).y == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(3).z == 1.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(4).x == -1.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(4).y == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(4).z == 0.0);

    // test when only the first reactant has a conf:
    reacts.clear();
    smi = "C(=O)C";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    conf = new Conformer(3);
    conf->setAtomPos(0, RDGeom::Point3D(1, 0, 0));
    conf->setAtomPos(1, RDGeom::Point3D(1, 1, 0));
    conf->setAtomPos(2, RDGeom::Point3D(1, 0, 1));
    mol->addConformer(conf, true);
    reacts.push_back(ROMOL_SPTR(mol));

    smi = "CN";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));

    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
    TEST_ASSERT(prods[0][0]->getNumConformers() == 1);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(0).x == 1.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(0).y == 1.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(0).z == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(1).x == 1.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(1).y == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(1).z == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(2).x == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(2).y == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(2).z == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(3).x == 1.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(3).y == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(3).z == 1.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(4).x == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(4).y == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(4).z == 0.0);

    // test when only the second reactant has a conf:
    reacts.clear();
    smi = "C(=O)C";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));

    smi = "CN";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    conf = new Conformer(2);
    conf->setAtomPos(0, RDGeom::Point3D(-1, 0, 0));
    conf->setAtomPos(1, RDGeom::Point3D(-1, 1, 0));
    mol->addConformer(conf, true);
    reacts.push_back(ROMOL_SPTR(mol));

    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
    TEST_ASSERT(prods[0][0]->getNumConformers() == 1);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(0).x == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(0).y == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(0).z == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(1).x == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(1).y == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(1).z == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(2).x == -1.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(2).y == 1.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(2).z == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(3).x == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(3).y == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(3).z == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(4).x == -1.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(4).y == 0.0);
    TEST_ASSERT(prods[0][0]->getConformer().getAtomPos(4).z == 0.0);

    // and, of course, when neither has a conformer:
    reacts.clear();
    smi = "C(=O)C";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));

    smi = "CN";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));

    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 5);
    TEST_ASSERT(prods[0][0]->getNumConformers() == 0);

    delete rxn;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test26V3000MDLParser() {
  ROMol *mol = nullptr;
  ChemicalReaction *rxn;
  MOL_SPTR_VECT reacts;
  std::vector<MOL_SPTR_VECT> prods;
  std::string smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing v3000 MDL parser" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName;

  fName = rdbase + "/Code/GraphMol/ChemReactions/testData/v3k.AmideBond.rxn";

  rxn = RxnFileToChemicalReaction(fName);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);

  smi = "C(=O)O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 4);
  TEST_ASSERT(prods[0][0]->getNumBonds() == 3);

  delete rxn;
  reacts.clear();
  fName = rdbase + "/Code/GraphMol/ChemReactions/testData/v3k.cyclization1.rxn";
  rxn = RxnFileToChemicalReaction(fName);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);

  smi = "OC(=O)CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "OC(=O)CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 8);
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == "O=C1CNC(=O)CN1");

  delete rxn;
  reacts.clear();
  fName = rdbase + "/Code/GraphMol/ChemReactions/testData/v3k.cyclization1.rxn";
  rxn = RxnFileToChemicalReaction(fName);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);

  smi = "OC(=O)CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "OC(=O)CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 8);
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == "O=C1CNC(=O)CN1");

  delete rxn;
  reacts.clear();
  fName = rdbase + "/Code/GraphMol/ChemReactions/testData/v3k.cyclization2.rxn";
  rxn = RxnFileToChemicalReaction(fName);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);

  smi = "CC=C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "C=CC=C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 4);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 7);
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[1][0]));
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[2][0]));
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[3][0]));

  reacts.clear();
  smi = "CC=C[Cl]";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "[F]C=CC=C[Br]";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 4);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 10);
  for (unsigned int i = 0; i < prods.size(); i++) {
    unsigned int nDifferent = 0;
    for (unsigned int j = 0; j < prods.size(); j++) {
      if (MolToSmiles(*prods[i][0]) != MolToSmiles(*prods[j][0])) {
        nDifferent += 1;
      }
    }
    TEST_ASSERT(nDifferent == 2);
  }

  reacts.clear();
  smi = "C1C=CCCC1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "C=CC=C";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 4);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 10);
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[1][0]));
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[2][0]));
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == MolToSmiles(*prods[3][0]));

  delete rxn;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test27SmartsWriter() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing reaction SMARTS writer." << std::endl;

  unsigned int nWarn, nError;

  {
    std::string smi;
    smi = "[C:1]=[O:2].[N:3]>>[O:2]=[C:1]~[N:3]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);

    smi = ChemicalReactionToRxnSmarts(*rxn);
    TEST_ASSERT(smi == "[C:1]=[O:2].[N:3]>>[O:2]=[C:1]~[N:3]");

    delete rxn;
    rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);
    delete rxn;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test28RxnDepictor() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing reaction depictor." << std::endl;

  unsigned int nWarn, nError;

  {
    std::string smi;
    smi =
        "[C:10]1[C:11][C:12]1[C:1]=[O:2].[N:3][C:20]([C:21])[C:22]>>[O:2]=[C:1]"
        "([C:10]1[C:11][C:12]1)[N:3][C:20]([C:21])[C:22]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);

    for (RDKit::MOL_SPTR_VECT::const_iterator templIt =
             rxn->beginReactantTemplates();
         templIt != rxn->endReactantTemplates(); ++templIt) {
      TEST_ASSERT((*templIt)->getNumConformers() == 0);
    }
    for (RDKit::MOL_SPTR_VECT::const_iterator templIt =
             rxn->beginProductTemplates();
         templIt != rxn->endProductTemplates(); ++templIt) {
      TEST_ASSERT((*templIt)->getNumConformers() == 0);
    }

    RDDepict::compute2DCoordsForReaction(*rxn);
    for (RDKit::MOL_SPTR_VECT::const_iterator templIt =
             rxn->beginReactantTemplates();
         templIt != rxn->endReactantTemplates(); ++templIt) {
      TEST_ASSERT((*templIt)->getNumConformers() == 1);
      TEST_ASSERT(!(*templIt)->getConformer().is3D());
    }
    for (RDKit::MOL_SPTR_VECT::const_iterator templIt =
             rxn->beginProductTemplates();
         templIt != rxn->endProductTemplates(); ++templIt) {
      TEST_ASSERT((*templIt)->getNumConformers() == 1);
      TEST_ASSERT(!(*templIt)->getConformer().is3D());
    }
    delete rxn;
  }
  {  // make sure the depiction doesn't screw up the reaction itself
    std::string rdbase = getenv("RDBASE");
    std::string fName;
    std::string smi;
    ChemicalReaction *rxn;
    MOL_SPTR_VECT reacts;
    std::vector<MOL_SPTR_VECT> prods;

    fName = rdbase + "/Code/GraphMol/ChemReactions/testData/cyclization1.rxn";
    rxn = RxnFileToChemicalReaction(fName);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);

    RDDepict::compute2DCoordsForReaction(*rxn);

    smi = "OC(=O)CN";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));

    smi = "OC(=O)CN";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));

    rxn->initReactantMatchers();
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 8);
    TEST_ASSERT(MolToSmiles(*prods[0][0]) == "O=C1CNC(=O)CN1");
    delete rxn;
  }
  {
    std::string smi = "[#7;!H0:1]>>[#7:1]C";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);

    MOL_SPTR_VECT reacts;
    reacts.clear();
    smi = "C1=CNC=C1";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    RDDepict::compute2DCoordsForReaction(*rxn);

    rxn->initReactantMatchers();

    std::vector<MOL_SPTR_VECT> prods;
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);

    ROMOL_SPTR prod = prods[0][0];
    MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));
    TEST_ASSERT(prod->getNumAtoms() == 6);
    TEST_ASSERT(prod->getAtomWithIdx(0)->getAtomicNum() == 7);
    TEST_ASSERT(prod->getAtomWithIdx(0)->getImplicitValence() == 0);
    TEST_ASSERT(prod->getAtomWithIdx(0)->getExplicitValence() == 3);
    TEST_ASSERT(prod->getAtomWithIdx(0)->getNoImplicit() == false);

    delete rxn;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test29RxnWriter() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing RXN file writer." << std::endl;

  unsigned int nWarn, nError;

  {
    std::string rdbase = getenv("RDBASE");
    std::string fName;
    std::string smi;

    fName = rdbase + "/Code/GraphMol/ChemReactions/testData/AmideBond.rxn";
    ChemicalReaction *rxn = RxnFileToChemicalReaction(fName);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);

    smi = ChemicalReactionToRxnSmarts(*rxn);
    TEST_ASSERT(smi == "[#6:2](-[#8])=[#8:1].[#7:3]>>[#7:3]-[#6:2]=[#8:1]")

    std::string mb;
    mb = ChemicalReactionToRxnBlock(*rxn);

    delete rxn;
    rxn = RxnBlockToChemicalReaction(mb);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);
    delete rxn;
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test30ReactProdQueries() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing reactant and product queries." << std::endl;

  {
    ROMol *mol;
    unsigned int nWarn, nError, which;
    std::string smi;
    smi =
        "[c;H:1]:[c:2](:[c;H:3])Br.[C:4](=[O:5])Cl>>[c:1]:[c:2]([c:3])-[C:4]=["
        "O:5]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);

    smi = "c1ccccc1Br";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(isMoleculeReactantOfReaction(*rxn, *mol, which));
    TEST_ASSERT(which == 0);
    TEST_ASSERT(isMoleculeReactantOfReaction(*rxn, *mol));
    delete (mol);

    smi = "c1ccccc1Cl";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(!isMoleculeReactantOfReaction(*rxn, *mol, which));
    TEST_ASSERT(which == rxn->getNumReactantTemplates());
    TEST_ASSERT(!isMoleculeReactantOfReaction(*rxn, *mol));
    delete (mol);

    smi = "c1ccncc1Br";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(isMoleculeReactantOfReaction(*rxn, *mol, which));
    TEST_ASSERT(which == 0);
    delete (mol);

    smi = "c1cccnc1Br";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(!isMoleculeReactantOfReaction(*rxn, *mol, which));
    TEST_ASSERT(which == rxn->getNumReactantTemplates());
    delete (mol);

    smi = "c1cccc(C)c1Br";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(!isMoleculeReactantOfReaction(*rxn, *mol, which));
    TEST_ASSERT(which == rxn->getNumReactantTemplates());
    delete (mol);

    smi = "ClC(=O)C";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(isMoleculeReactantOfReaction(*rxn, *mol, which));
    TEST_ASSERT(which == 1);
    delete (mol);

    smi = "C(=O)C";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(!isMoleculeReactantOfReaction(*rxn, *mol, which));
    TEST_ASSERT(which == rxn->getNumReactantTemplates());
    delete (mol);

    smi = "c1ccccc1C(=O)C";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(isMoleculeProductOfReaction(*rxn, *mol, which));
    TEST_ASSERT(which == 0);
    TEST_ASSERT(isMoleculeProductOfReaction(*rxn, *mol));
    delete (mol);

    smi = "c1cccnc1C(=O)C";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(!isMoleculeProductOfReaction(*rxn, *mol, which));
    TEST_ASSERT(which == rxn->getNumProductTemplates());
    TEST_ASSERT(!isMoleculeProductOfReaction(*rxn, *mol));
    delete (mol);

    delete (rxn);
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test31Issue3140490() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Issue 3140490." << std::endl;

  {
    unsigned int nWarn, nError;
    std::string smi;
    smi = "[O:1]>>[N:1]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);

    MOL_SPTR_VECT reacts;
    reacts.clear();
    smi = "OC";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    std::vector<MOL_SPTR_VECT> prods;
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(!prods[0][0]->getAtomWithIdx(0)->hasProp(
        common_properties::molAtomMapNumber));
    TEST_ASSERT(!prods[0][0]->getAtomWithIdx(1)->hasProp(
        common_properties::molAtomMapNumber));
    delete (rxn);
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test32Replacements() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing string replacement in parsing." << std::endl;

  {
    std::map<std::string, std::string> repls;
    repls["{amine}"] = "$([N;!H0;$(N-[#6]);!$(N-[!#6;!#1]);!$(N-C=[O,N,S])])";
    unsigned int nWarn, nError;
    std::string smi;
    smi = "[{amine}:1]>>[*:1]-C";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi, &repls);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);

    MOL_SPTR_VECT reacts;
    reacts.clear();
    smi = "CCN";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    std::vector<MOL_SPTR_VECT> prods;
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 4);
    delete (rxn);
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test33ReactingAtoms1() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing getReactingAtoms() 1." << std::endl;

  {  // basics
    std::string smi;
    smi = "[O:1][C:2]>>[N:1][C:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn);
    TEST_ASSERT(ratoms.size() == 1);
    TEST_ASSERT(ratoms[0].size() == 1);
    delete rxn;
  }
  {  // no changes
    std::string smi;
    smi = "[O:1][C:2]>>[O:1][C:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn);
    TEST_ASSERT(ratoms.size() == 1);
    TEST_ASSERT(ratoms[0].size() == 0);
    delete rxn;
  }
  {  // make sure atomic number queries work:
    std::string smi;
    smi = "[#8:1][C:2]>>[O:1][C:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn);
    TEST_ASSERT(ratoms.size() == 1);
    TEST_ASSERT(ratoms[0].size() == 0);
    delete rxn;
  }
  {  // query in the reactants, dummy in product
    std::string smi;
    smi = "[O,N:1][C:2]>>[*:1][C:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn);
    TEST_ASSERT(ratoms.size() == 1);
    TEST_ASSERT(ratoms[0].size() == 0);
    delete rxn;
  }
  {  // recursive query in the reactants without an atomic number query
    std::string smi;
    smi = "[$(O):1][C:2]>>[O:1][C:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn);
    TEST_ASSERT(ratoms.size() == 1);
    TEST_ASSERT(ratoms[0].size() == 1);
    delete rxn;
  }
  {  // recursive query with atomic number query
    std::string smi;
    smi = "[O;$(O):1][C:2]>>[O:1][C:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn);
    TEST_ASSERT(ratoms.size() == 1);
    TEST_ASSERT(ratoms[0].size() == 0);
    delete rxn;
  }
  {  // recursive query with atomic number query, alternate ordering
    // FIX: this returns a changed atom (since we don't know the atomic
    // number of the atom) but probably shouldn't
    std::string smi;
    smi = "[$(O);O:1][C:2]>>[O:1][C:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn);
    TEST_ASSERT(ratoms.size() == 1);
    TEST_ASSERT(ratoms[0].size() == 1);
    delete rxn;
  }
  {  // recursive query in the reactants, dummy in products
    std::string smi;
    smi = "[$(O):1][C:2]>>[*:1][C:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn);
    TEST_ASSERT(ratoms.size() == 1);
    TEST_ASSERT(ratoms[0].size() == 0);
    delete rxn;
  }
  {  // query with degree/H info in the reactants:
    std::string smi;
    smi = "[O;H1:1][C:2]>>[O:1][C:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn);
    TEST_ASSERT(ratoms.size() == 1);
    TEST_ASSERT(ratoms[0].size() == 0);
    delete rxn;
  }
  {  // dummy in the reactants, dummy in products
    std::string smi;
    smi = "[*:1][C:2]>>[*:1][C:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn);
    TEST_ASSERT(ratoms.size() == 1);
    TEST_ASSERT(ratoms[0].size() == 0);
    delete rxn;
  }
  {  // two reactants, one changes
    std::string smi;
    smi = "[O:1][C:2].[N:3]>>[N:1][C:2].[N:3]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn);
    TEST_ASSERT(ratoms.size() == 2);
    TEST_ASSERT(ratoms[0].size() == 1);
    TEST_ASSERT(ratoms[1].size() == 0);
    delete rxn;
  }
  {  // two reactants, one changes, reordered
    std::string smi;
    smi = "[O:1][C:2].[N:3]>>[N:3].[N:1][C:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn);
    TEST_ASSERT(ratoms.size() == 2);
    TEST_ASSERT(ratoms[0].size() == 1);
    TEST_ASSERT(ratoms[1].size() == 0);
    delete rxn;
  }
  {  // dummies for duplicating atom properties
    std::string smi;
    smi = "[O:1][C:2].[N:3]>>[N:1][C:2].[*:3]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn);
    TEST_ASSERT(ratoms.size() == 2);
    TEST_ASSERT(ratoms[0].size() == 1);
    TEST_ASSERT(ratoms[1].size() == 0);
    delete rxn;
  }
  {  // query in the reactants, dummy in products
    std::string smi;
    smi = "[O,N:1][C:2]>>[*:1][C:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn);
    TEST_ASSERT(ratoms.size() == 1);
    TEST_ASSERT(ratoms[0].size() == 0);
    delete rxn;
  }
  {  // changed degree
    std::string smi;
    smi = "[O:1][C:2].[N:3]>>[N:1][C:2].[*:3]C";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn);
    TEST_ASSERT(ratoms.size() == 2);
    TEST_ASSERT(ratoms[0].size() == 1);
    TEST_ASSERT(ratoms[1].size() == 1);
    delete rxn;
  }
  {  // unmapped atoms in reactants:
    std::string smi;
    smi = "[O:1]C>>[O:1]C";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn);
    TEST_ASSERT(ratoms.size() == 1);
    TEST_ASSERT(ratoms[0].size() == 2);
    delete rxn;
  }
  {  // don't return info about unmapped atoms:
    std::string smi;
    smi = "[O:1]C>>[O:1]C";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn, true);
    TEST_ASSERT(ratoms.size() == 1);
    TEST_ASSERT(ratoms[0].size() == 1);
    delete rxn;
  }
  {  // changing atom order
    std::string smi;
    smi = "[C:1]-[C:2]>>[C:2]-[C:1]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn);
    TEST_ASSERT(ratoms.size() == 1);
    TEST_ASSERT(ratoms[0].size() == 0);
    delete rxn;
  }
  {  // changing bond order
    std::string smi;
    smi = "[C:1]-[C:2]>>[C:1]=[C:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn);
    TEST_ASSERT(ratoms.size() == 1);
    TEST_ASSERT(ratoms[0].size() == 2);
    delete rxn;
  }
  {  // changing bond orders
    std::string smi;
    smi = "[C:1]-[C:2]=[C:3]>>[C:1]=[C:2]-[C:3]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn);
    TEST_ASSERT(ratoms.size() == 1);
    TEST_ASSERT(ratoms[0].size() == 3);
    delete rxn;
  }
  {  // query bond order
    std::string smi;
    smi = "[C:1]-[C:2]>>[C:1]~[C:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn);
    TEST_ASSERT(ratoms.size() == 1);
    TEST_ASSERT(ratoms[0].size() == 0);
    delete rxn;
  }
  {  // changing connectivity
    std::string smi;
    smi = "[C:1]-[C:2]-[C:3]>>[C:1]~[C:3]~[C:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn);
    TEST_ASSERT(ratoms.size() == 1);
    TEST_ASSERT(ratoms[0].size() == 3);
    delete rxn;
  }
  {  // changing connectivity 2
    std::string smi;
    smi = "[C:1]1[C:2][C:3]1>>[C:1][C:2][C:3]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn);
    TEST_ASSERT(ratoms.size() == 1);
    TEST_ASSERT(ratoms[0].size() == 2);
    delete rxn;
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test34ReactingAtoms2() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing getReactingAtoms() 2" << std::endl;

  {
    std::string rdbase = getenv("RDBASE");
    std::string fName =
        rdbase + "/Code/GraphMol/ChemReactions/testData/AmideBond.rxn";
    ChemicalReaction *rxn = RxnFileToChemicalReaction(fName);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn);
    TEST_ASSERT(ratoms.size() == 2);
    TEST_ASSERT(ratoms[0].size() == 2);
    TEST_ASSERT(ratoms[1].size() == 1);
    delete rxn;
  }
  {
    std::string rdbase = getenv("RDBASE");
    std::string fName =
        rdbase + "/Code/GraphMol/ChemReactions/testData/cyclization1.rxn";
    ChemicalReaction *rxn = RxnFileToChemicalReaction(fName);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    VECT_INT_VECT ratoms = getReactingAtoms(*rxn);
    TEST_ASSERT(ratoms.size() == 2);
    TEST_ASSERT(ratoms[0].size() == 3);
    TEST_ASSERT(ratoms[1].size() == 3);
    delete rxn;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test35ParensInReactants1() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing grouping parens in reactants1" << std::endl;

  {
    std::string smi = "[C:1].[C:2]>>[C:1][C:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    delete rxn;
  }

  {
    std::string smi = "([C:1].[C:2])>>[C:1][C:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    delete rxn;
  }
  {
    std::string smi = "([C:1].C(=O)O).[C:2]>>[C:1][C:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    delete rxn;
  }

  {
    std::string smi = "[C:1](=O)O.[C:2]>>[C:1][C:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    delete rxn;
  }

  {
    std::string smi = "[C:1](=O)O.([C:2])>>[C:1][C:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    delete rxn;
  }
  {
    std::string smi = "[C:1](=O)O.([C:2].N)>>[C:1][C:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    delete rxn;
  }
  {
    std::string smi = "([C:1](=O)O.[C:2]>>[C:1][C:2]";
    ChemicalReaction *rxn = nullptr;
    try {
      rxn = RxnSmartsToChemicalReaction(smi);
    } catch (const ChemicalReactionParserException &) {
      rxn = nullptr;
    }
    TEST_ASSERT(!rxn);
    delete rxn;
  }
  {
    std::string smi = "[C:1](=O)O.([C:2]>>[C:1][C:2]";
    ChemicalReaction *rxn = nullptr;
    try {
      rxn = RxnSmartsToChemicalReaction(smi);
      TEST_ASSERT(!rxn);
    } catch (const ChemicalReactionParserException &) {
      ;
    }
    delete rxn;
  }
  {
    std::string smi = "[C:1](=O)O).[C:2]>>[C:1][C:2]";
    ChemicalReaction *rxn = nullptr;
    try {
      rxn = RxnSmartsToChemicalReaction(smi);
      TEST_ASSERT(!rxn);
    } catch (const ChemicalReactionParserException &) {
      ;
    }
    delete rxn;
  }
  {
    std::string smi = "[C:1](=O)O.[C:2])>>[C:1][C:2]";
    ChemicalReaction *rxn = nullptr;
    try {
      rxn = RxnSmartsToChemicalReaction(smi);
      TEST_ASSERT(!rxn);
    } catch (const ChemicalReactionParserException &) {
      ;
    }
    delete rxn;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test36ParensInReactants2() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing grouping parens in reactants2" << std::endl;

  {
    std::string smi = "([C:1].[O:2])>>[C:1][O:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();
    MOL_SPTR_VECT reacts;
    reacts.clear();
    smi = "CNO";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    std::vector<MOL_SPTR_VECT> prods;
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 3);
    TEST_ASSERT(prods[0][0]->getNumBonds() == 3);

    delete rxn;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test37ProtectOption() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing use of _protected property on atoms."
                       << std::endl;

  {
    unsigned int nWarn, nError;
    std::string smi = "[O:1]>>[N:1]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);

    MOL_SPTR_VECT reacts;
    reacts.clear();
    smi = "OCO";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    std::vector<MOL_SPTR_VECT> prods;
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 2);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[1].size() == 1);

    reacts[0]->getAtomWithIdx(0)->setProp(common_properties::_protected, 1);
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);

    delete (rxn);
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test38AddRecursiveQueriesToReaction() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing use of adding recursive queries to a reaction." << std::endl;

  {
    ROMol *mol = nullptr;
    ChemicalReaction rxn;
    MOL_SPTR_VECT reacts;
    std::string smi;

    smi = "[C:1](=[O:2])";
    mol = SmartsToMol(smi);
    mol->getAtomWithIdx(0)->setProp("replaceme", "foo");
    TEST_ASSERT(mol);
    rxn.addReactantTemplate(ROMOL_SPTR(mol));
    TEST_ASSERT(rxn.getNumReactantTemplates() == 1);

    smi = "[N:3][C:4]";
    mol = SmartsToMol(smi);
    TEST_ASSERT(mol);
    rxn.addReactantTemplate(ROMOL_SPTR(mol));
    TEST_ASSERT(rxn.getNumReactantTemplates() == 2);

    smi = "[C:1](=[O:2])[N:3][C:4]";
    mol = SmartsToMol(smi);
    TEST_ASSERT(mol);
    rxn.addProductTemplate(ROMOL_SPTR(mol));
    TEST_ASSERT(rxn.getNumProductTemplates() == 1);

    std::string smi2 = "CCl";
    ROMOL_SPTR q1(SmilesToMol(smi2));
    std::map<std::string, ROMOL_SPTR> mp;
    mp["foo"] = q1;

    bool ok = false;
    try {
      addRecursiveQueriesToReaction(rxn, mp, "replaceme");
    } catch (ChemicalReactionException &) {
      ok = true;
    }
    TEST_ASSERT(ok);

    rxn.initReactantMatchers();
    addRecursiveQueriesToReaction(rxn, mp, "replaceme");
    MatchVectType mv;
    std::string msmi = "C(=O)Cl";
    ROMol *mmol = SmilesToMol(msmi);
    TEST_ASSERT(SubstructMatch(*mmol, *(*(rxn.beginReactantTemplates())), mv));
    delete mmol;
    msmi = "C(=O)O";
    mmol = SmilesToMol(msmi);
    TEST_ASSERT(!SubstructMatch(*mmol, *(*(rxn.beginReactantTemplates())), mv));
    delete mmol;
  }

  {
    ROMol *mol = nullptr;
    ChemicalReaction rxn;
    MOL_SPTR_VECT reacts;
    std::string smi;

    smi = "[C:1](=[O:2])";
    mol = SmartsToMol(smi);
    mol->getAtomWithIdx(0)->setProp("replaceme", "foo");
    TEST_ASSERT(mol);
    rxn.addReactantTemplate(ROMOL_SPTR(mol));
    TEST_ASSERT(rxn.getNumReactantTemplates() == 1);

    smi = "[N:3][C:4]";
    mol = SmartsToMol(smi);
    TEST_ASSERT(mol);
    rxn.addReactantTemplate(ROMOL_SPTR(mol));
    TEST_ASSERT(rxn.getNumReactantTemplates() == 2);

    smi = "[C:1](=[O:2])[N:3][C:4]";
    mol = SmartsToMol(smi);
    TEST_ASSERT(mol);
    rxn.addProductTemplate(ROMOL_SPTR(mol));
    TEST_ASSERT(rxn.getNumProductTemplates() == 1);

    std::string smi2 = "CO";
    ROMOL_SPTR q1(SmilesToMol(smi2));
    std::map<std::string, ROMOL_SPTR> mp;
    mp["foo"] = q1;

    rxn.initReactantMatchers();
    std::vector<std::vector<std::pair<unsigned int, std::string>>> labels;
    addRecursiveQueriesToReaction(rxn, mp, "replaceme", &labels);
    TEST_ASSERT(labels.size() == 2);
    TEST_ASSERT(labels[0][0].second == "foo");
  }

  {
    ROMol *mol = nullptr;
    ChemicalReaction rxn;
    MOL_SPTR_VECT reacts;
    std::string smi;

    smi = "[N:3][C:4]";
    mol = SmartsToMol(smi);
    TEST_ASSERT(mol);
    rxn.addReactantTemplate(ROMOL_SPTR(mol));
    TEST_ASSERT(rxn.getNumReactantTemplates() == 1);

    smi = "[C:1](=[O:2])";
    mol = SmartsToMol(smi);
    mol->getAtomWithIdx(0)->setProp("replaceme", "foo");
    TEST_ASSERT(mol);
    rxn.addReactantTemplate(ROMOL_SPTR(mol));
    TEST_ASSERT(rxn.getNumReactantTemplates() == 2);

    smi = "[C:1](=[O:2])[N:3][C:4]";
    mol = SmartsToMol(smi);
    TEST_ASSERT(mol);
    rxn.addProductTemplate(ROMOL_SPTR(mol));
    TEST_ASSERT(rxn.getNumProductTemplates() == 1);

    std::string smi2 = "CO";
    ROMOL_SPTR q1(SmilesToMol(smi2));
    std::map<std::string, ROMOL_SPTR> mp;
    mp["foo"] = q1;

    rxn.initReactantMatchers();
    std::vector<std::vector<std::pair<unsigned int, std::string>>> labels;
    addRecursiveQueriesToReaction(rxn, mp, "replaceme", &labels);
    TEST_ASSERT(labels.size() == 2);
    TEST_ASSERT(labels[1][0].second == "foo");
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test39InnocentChiralityLoss() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing loss of chirality on innocent atoms."
                       << std::endl;

  {
    ChemicalReaction *rxn =
        RxnSmartsToChemicalReaction("[C:2][C:1]=O>>[C:2][C:1]=S");
    unsigned int nWarn, nError;
    rxn->initReactantMatchers();
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);

    {
      std::string smi = "Cl[C@H](F)C=O";
      ROMol *mol = SmilesToMol(smi);
      MOL_SPTR_VECT reacts;
      reacts.push_back(ROMOL_SPTR(mol));
      std::vector<MOL_SPTR_VECT> prods;
      prods = rxn->runReactants(reacts);
      TEST_ASSERT(prods.size() == 1);
      TEST_ASSERT(prods[0].size() == 1);
      smi = MolToSmiles(*prods[0][0], false);
      TEST_ASSERT(smi == "FC(Cl)C=S");
      smi = MolToSmiles(*prods[0][0], true);
      TEST_ASSERT(smi == "F[C@@H](Cl)C=S");
    }
    {
      std::string smi = "O=C[C@@H](F)Cl";
      ROMol *mol = SmilesToMol(smi);
      MOL_SPTR_VECT reacts;
      reacts.push_back(ROMOL_SPTR(mol));
      std::vector<MOL_SPTR_VECT> prods;
      prods = rxn->runReactants(reacts);
      TEST_ASSERT(prods.size() == 1);
      TEST_ASSERT(prods[0].size() == 1);
      smi = MolToSmiles(*prods[0][0], false);
      TEST_ASSERT(smi == "FC(Cl)C=S");
      smi = MolToSmiles(*prods[0][0], true);
      TEST_ASSERT(smi == "F[C@@H](Cl)C=S");
    }
    {
      std::string smi = "F[C@H](C=O)Cl";
      ROMol *mol = SmilesToMol(smi);
      MOL_SPTR_VECT reacts;
      reacts.push_back(ROMOL_SPTR(mol));
      std::vector<MOL_SPTR_VECT> prods;
      prods = rxn->runReactants(reacts);
      TEST_ASSERT(prods.size() == 1);
      TEST_ASSERT(prods[0].size() == 1);
      smi = MolToSmiles(*prods[0][0], false);
      TEST_ASSERT(smi == "FC(Cl)C=S");
      smi = MolToSmiles(*prods[0][0], true);
      TEST_ASSERT(smi == "F[C@@H](Cl)C=S");
    }
    delete rxn;
  }

  {
    // in this case the atom isn't actually innocent, but we can handle it
    // anyway because
    // only one bond changes
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction("[C:1]-O>>[C:1]-S");
    unsigned int nWarn, nError;
    rxn->initReactantMatchers();
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);

    {
      std::string smi = "Cl[C@H](F)O";
      ROMol *mol = SmilesToMol(smi);
      MOL_SPTR_VECT reacts;
      reacts.push_back(ROMOL_SPTR(mol));
      std::vector<MOL_SPTR_VECT> prods;
      prods = rxn->runReactants(reacts);
      TEST_ASSERT(prods.size() == 1);
      TEST_ASSERT(prods[0].size() == 1);
      smi = MolToSmiles(*prods[0][0], false);
      TEST_ASSERT(smi == "FC(S)Cl");
      smi = MolToSmiles(*prods[0][0], true);
      TEST_ASSERT(smi == "F[C@H](S)Cl");
    }
    delete rxn;
  }
  {
    // another non-innocent atom
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction("[C:1]-O>>[C:1]");
    unsigned int nWarn, nError;
    rxn->initReactantMatchers();
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);

    {
      std::string smi = "Cl[C@H](F)O";
      ROMol *mol = SmilesToMol(smi);
      MOL_SPTR_VECT reacts;
      reacts.push_back(ROMOL_SPTR(mol));
      std::vector<MOL_SPTR_VECT> prods;
      prods = rxn->runReactants(reacts);
      TEST_ASSERT(prods.size() == 1);
      TEST_ASSERT(prods[0].size() == 1);
      smi = MolToSmiles(*prods[0][0], false);
      TEST_ASSERT(smi == "FCCl");
      smi = MolToSmiles(*prods[0][0], true);
      TEST_ASSERT(smi == "FCCl");
    }
    delete rxn;
  }
  {
    // one of the original bug report cases:
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(
        "[C:1](=[O:2])[C:3][C:4]([OH:5])[#6:6]>>[C:1](=[O:2])[C:3][H].[C:4](=["
        "O:5])[#6:6]");
    unsigned int nWarn, nError;
    rxn->initReactantMatchers();
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);

    {
      std::string smi = "CC(O)[C@](N)(F)C(C)=O";
      ROMol *mol = SmilesToMol(smi);
      MOL_SPTR_VECT reacts;
      reacts.push_back(ROMOL_SPTR(mol));
      std::vector<MOL_SPTR_VECT> prods;
      prods = rxn->runReactants(reacts);
      TEST_ASSERT(prods.size() == 1);
      TEST_ASSERT(prods[0].size() == 2);
      smi = MolToSmiles(*prods[0][0], false);
      TEST_ASSERT(smi == "[H]C(N)(F)C(C)=O");
      smi = MolToSmiles(*prods[0][0], true);
      TEST_ASSERT(smi == "[H][C@](N)(F)C(C)=O");
    }
    delete rxn;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test40AgentsInSmarts() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing handling of agents in reaction smarts"
                       << std::endl;
  // This was github #222
  {
    ChemicalReaction *rxn;
    std::string smi;
    smi = "[C:1](=[O:2])O.[N:3][C:4]>[Pd]>[C:1](=[O:2])[N:3][C:4]";
    rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    delete rxn;
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test41Github233() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing github 233: chirality not preserved in nonmapped atoms"
      << std::endl;
  {
    ChemicalReaction *rxn;
    std::string smi;
    smi = "[F:1][C:2]([C:3])[I:4]>>[F:1][C:2]([C:3][C@H]([OH])Br)[Cl:4]";
    rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);

    unsigned int nWarn, nError;
    rxn->initReactantMatchers();
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);

    smi = "FC(C)I";
    ROMol *mol = SmilesToMol(smi);
    MOL_SPTR_VECT reacts;
    reacts.push_back(ROMOL_SPTR(mol));
    std::vector<MOL_SPTR_VECT> prods;
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);

    MolOps::sanitizeMol(*(static_cast<RWMol *>(prods[0][0].get())));
    smi = MolToSmiles(*prods[0][0], true);
    TEST_ASSERT(smi == "O[C@H](Br)CC(F)Cl");
    delete rxn;
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test42ReactionSmiles() {
  ROMol *mol = nullptr;
  ChemicalReaction *rxn;
  MOL_SPTR_VECT reacts;
  std::vector<MOL_SPTR_VECT> prods;
  std::string smi;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Daylight parser" << std::endl;

  smi = "[C:1](=[O:2])O.[N:3][C:4]>>[C:1](=[O:2])[N:3][C:4]";
  rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);

  smi = "C(=O)O";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 4);
  TEST_ASSERT(prods[0][0]->getNumBonds() == 3);

  delete rxn;
  reacts.clear();
  smi =
      "[N:1][C:2][C:3](=[O:4])[O:5].[N:6][C:7][C:8](=[O:9])[O:10]>>[N:1]1[C:2]["
      "C:3](=[O:4])[N:6][C:7][C:8]1=[O:9].[O:5][O:10]";
  rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
  TEST_ASSERT(rxn);

  TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
  TEST_ASSERT(rxn->getNumProductTemplates() == 2);

  smi = "OC(=O)CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  smi = "OC(=O)CN";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol);
  reacts.push_back(ROMOL_SPTR(mol));

  rxn->initReactantMatchers();
  prods = rxn->runReactants(reacts);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 2);
  TEST_ASSERT(prods[0][0]->getNumAtoms() == 8);
  TEST_ASSERT(prods[0][1]->getNumAtoms() == 2);
  TEST_ASSERT(MolToSmiles(*prods[0][0]) == "O=C1CNC(=O)CN1");
  TEST_ASSERT(MolToSmiles(*prods[0][1]) == "OO");

  delete rxn;
  reacts.clear();
  smi =
      "[N:1][C:2][C:3](=[O:4])[O:5].[N:6][C:7][C:8](=[O:9])[O:10].[N:1]1[C:2]["
      "C:3](=[O:4])[N:6][C:7][C:8]1=[O:9].[O:5][O:10]";
  try {
    rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
  } catch (ChemicalReactionParserException &) {
  }

  smi =
      "[N:1][C:2][C:3](=[O:4])[O:5].[N:6][C:7][C:8](=[O:9])[O:10]>>[N:1]1[C:2]["
      "C:3](=[O:4])[N:6][C:7][C:8]1=[O:9]>>[O:5][O:10]";
  try {
    rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
  } catch (ChemicalReactionParserException &) {
  }

  smi =
      "[Q:1][C:2][C:3](=[O:4])[O:5].[N:6][C:7][C:8](=[O:9])[O:10]>>[N:1]1[C:2]["
      "C:3](=[O:4])[N:6][C:7][C:8]1=[O:9].[O:5][O:10]";
  try {
    rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
  } catch (ChemicalReactionParserException &) {
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test43Github243() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing github 243: dummy labels copied into products" << std::endl;
  {
    std::string rdbase = getenv("RDBASE");
    std::string fName =
        rdbase + "/Code/GraphMol/ChemReactions/testData/github243.rxn";
    ChemicalReaction *rxn = RxnFileToChemicalReaction(fName);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();

    MOL_SPTR_VECT reacts;
    reacts.clear();
    std::string smi = "CCCN";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    std::vector<MOL_SPTR_VECT> prods;
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getAtomWithIdx(9)->getAtomicNum() == 6);
    TEST_ASSERT(!prods[0][0]->getAtomWithIdx(9)->hasProp(
        common_properties::dummyLabel));
    TEST_ASSERT(!prods[0][0]->getAtomWithIdx(9)->hasProp(
        common_properties::_MolFileRLabel));
    TEST_ASSERT(prods[0][0]->getAtomWithIdx(9)->getIsotope() == 0);

    TEST_ASSERT(prods[0][0]->getAtomWithIdx(10)->getAtomicNum() == 0);
    TEST_ASSERT(prods[0][0]->getAtomWithIdx(10)->hasProp(
        common_properties::dummyLabel));
    TEST_ASSERT(prods[0][0]->getAtomWithIdx(10)->hasProp(
        common_properties::_MolFileRLabel));

    delete (rxn);
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test44Github290() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing github 290: seg fault while parsing rxn"
                       << std::endl;
  {
    auto *rxn = new ChemicalReaction();
    delete rxn;
  }
  {
    std::string rdbase = getenv("RDBASE");
    std::string fName =
        rdbase + "/Code/GraphMol/ChemReactions/testData/bogus_github290.rxn";
    bool failed = false;
    try {
      RxnFileToChemicalReaction(fName);
    } catch (ChemicalReactionParserException &) {
      failed = true;
    }
    TEST_ASSERT(failed);
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test45SmilesWriter() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing reaction SMILES writer." << std::endl;
  unsigned int nWarn, nError;

  {
    std::string smi;
    smi = "[C:1]=[O:2].[N:3]>>[N:3]~[C:1]=[O:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);

    std::string res = "";
    for (MOL_SPTR_VECT::const_iterator iter = rxn->beginReactantTemplates();
         iter != rxn->endReactantTemplates(); ++iter) {
      if (iter != rxn->beginReactantTemplates()) {
        res += ".";
      }
      res += MolToSmiles(**iter, true);
    }
    res += ">>";
    for (MOL_SPTR_VECT::const_iterator iter = rxn->beginProductTemplates();
         iter != rxn->endProductTemplates(); ++iter) {
      if (iter != rxn->beginProductTemplates()) {
        res += ".";
      }
      res += MolToSmiles(**iter, true);
    }

    smi = ChemicalReactionToRxnSmiles(*rxn);
    TEST_ASSERT(smi == res)
    TEST_ASSERT(smi == "[C:1]=[O:2].[N:3]>>[C:1](=[O:2])~[N:3]");
    delete rxn;

    rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);
    delete rxn;
  }
  {
    std::string smi;
    smi = "C=O.N>>N~C=O";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);

    std::string res = "";
    for (MOL_SPTR_VECT::const_iterator iter = rxn->beginReactantTemplates();
         iter != rxn->endReactantTemplates(); ++iter) {
      if (iter != rxn->beginReactantTemplates()) {
        res += ".";
      }
      res += MolToSmiles(**iter, true);
    }
    res += ">>";
    for (MOL_SPTR_VECT::const_iterator iter = rxn->beginProductTemplates();
         iter != rxn->endProductTemplates(); ++iter) {
      if (iter != rxn->beginProductTemplates()) {
        res += ".";
      }
      res += MolToSmiles(**iter, true);
    }

    smi = ChemicalReactionToRxnSmiles(*rxn);
    TEST_ASSERT(smi == res)
    TEST_ASSERT(smi == "C=O.N>>N~C=O");

    delete rxn;
    rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    delete rxn;
  }
  {
    std::string smi;
    smi = "[S,C]=O.N>>N~[S,C]=O";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi, nullptr, false);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);

    std::string res = "";
    for (MOL_SPTR_VECT::const_iterator iter = rxn->beginReactantTemplates();
         iter != rxn->endReactantTemplates(); ++iter) {
      if (iter != rxn->beginReactantTemplates()) {
        res += ".";
      }
      res += MolToSmiles(**iter, true);
    }
    res += ">>";
    for (MOL_SPTR_VECT::const_iterator iter = rxn->beginProductTemplates();
         iter != rxn->endProductTemplates(); ++iter) {
      if (iter != rxn->beginProductTemplates()) {
        res += ".";
      }
      res += MolToSmiles(**iter, true);
    }

    smi = ChemicalReactionToRxnSmiles(*rxn, false);
    TEST_ASSERT(smi == res)
    TEST_ASSERT(smi != "C=O.N>>N~C=O");
    TEST_ASSERT(smi == "*=O.N>>N~*=O");

    delete rxn;
    rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    delete rxn;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test46Agents() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing handling of reaction agents." << std::endl;
  unsigned int nWarn, nError;

  {
    std::string smi;
    ROMol *mol = nullptr;

    smi = "[C:1]=[O:2].[N:3]>[OH2].[Na].[Cl]>[N:3]~[C:1]=[O:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 3);
    rxn->initReactantMatchers();
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);

    smi = "C(=O)O";
    mol = SmartsToMol(smi);
    TEST_ASSERT(mol);
    rxn->addAgentTemplate(ROMOL_SPTR(mol));
    TEST_ASSERT(rxn->getNumAgentTemplates() == 4);

    delete rxn;
  }
  {
    std::string smi;

    smi = ">[OH2].[Na].[Cl]>";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 0);
    TEST_ASSERT(rxn->getNumProductTemplates() == 0);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 3);

    delete rxn;
  }
  {
    std::string smi;
    ROMol *mol = nullptr;
    MOL_SPTR_VECT agents;

    smi = "[C:1]=[O:2].[N:3]>>[N:3]~[C:1]=[O:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 0);
    rxn->initReactantMatchers();
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);

    smi = "C(=O)O";
    mol = SmartsToMol(smi);
    TEST_ASSERT(mol);
    rxn->addAgentTemplate(ROMOL_SPTR(mol));
    TEST_ASSERT(rxn->getNumAgentTemplates() == 1);

    delete rxn;
  }
  {
    std::string smi;
    ROMol *mol = nullptr;
    unsigned int nWarn, nError, which;
    MOL_SPTR_VECT agents;

    smi = "[C:1]=[O:2].[N:3]>[OH2].[Na].[Cl]>[N:3]~[C:1]=[O:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 3);
    rxn->initReactantMatchers();
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);

    smi = "O";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(isMoleculeAgentOfReaction(*rxn, *mol, which));
    TEST_ASSERT(which == 0);
    TEST_ASSERT(isMoleculeAgentOfReaction(*rxn, *mol));
    delete (mol);

    smi = "C(=O)N";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(!isMoleculeAgentOfReaction(*rxn, *mol, which));
    TEST_ASSERT(which == rxn->getNumAgentTemplates());
    TEST_ASSERT(!isMoleculeAgentOfReaction(*rxn, *mol));
    delete (mol);

    smi = "C(=O)O";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(!isMoleculeAgentOfReaction(*rxn, *mol, which));
    TEST_ASSERT(which == rxn->getNumAgentTemplates());
    TEST_ASSERT(!isMoleculeAgentOfReaction(*rxn, *mol));
    delete (mol);

    smi = "C(=O)O";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    rxn->addAgentTemplate(ROMOL_SPTR(mol));
    TEST_ASSERT(rxn->getNumAgentTemplates() == 4);

    smi = "C(=O)O";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    TEST_ASSERT(isMoleculeAgentOfReaction(*rxn, *mol, which));
    TEST_ASSERT(which == 3);
    delete (mol);

    delete rxn;
  }
  {
    std::string smi1, smi2;

    smi1 = "[C:1]=[O:2].[N:3]>[OH2].[Na].[Cl]>[N:3]~[C:1]=[O:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi1, nullptr, true);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 3);
    rxn->initReactantMatchers();
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);

    smi2 = "[C:1]=[O:2].[N:3]>>[N:3]~[C:1]=[O:2]";
    ChemicalReaction *rxnq = RxnSmartsToChemicalReaction(smi2, nullptr, true);
    TEST_ASSERT(rxnq);
    TEST_ASSERT(rxnq->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxnq->getNumProductTemplates() == 1);
    TEST_ASSERT(rxnq->getNumAgentTemplates() == 0);
    rxnq->initReactantMatchers();
    TEST_ASSERT(rxnq->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);

    TEST_ASSERT(hasReactionSubstructMatch(*rxn, *rxnq, true));

    delete rxn;
    delete rxnq;
  }
  {
    std::string smi1;

    smi1 = "[C:1]=[O:2].[N:3].C(=O)O>[OH2].[Na].[Cl]>[N:3]~[C:1]=[O:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi1, nullptr, true);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 3);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 3);
    rxn->initReactantMatchers();

    rxn->removeUnmappedReactantTemplates();
    rxn->removeUnmappedProductTemplates();
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 4);
    delete rxn;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test47TestReactionMoleculeConversion() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing the conversion of a molecule with Rxn role "
                          "to a reaction and vice versa."
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/ChemReactions/testData/";

  {
    std::string fName;
    fName = rdbase + "rxn1.mol";
    ROMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 18);
    TEST_ASSERT(m->getNumBonds() == 16);

    ChemicalReaction *rxn = RxnMolToChemicalReaction(*m);

    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 0);

    ROMol *mol = ChemicalReactionToRxnMol(*rxn);

    TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp(common_properties::molRxnRole));
    TEST_ASSERT(mol->getAtomWithIdx(0)->getProp<int>(
                    common_properties::molRxnRole) == 1);
    TEST_ASSERT(
        mol->getAtomWithIdx(10)->hasProp(common_properties::molRxnRole));
    TEST_ASSERT(mol->getAtomWithIdx(10)->getProp<int>(
                    common_properties::molRxnRole) == 2);

    std::string smi1 = MolToSmiles(*m);
    std::string smi2 = MolToSmiles(*mol);

    TEST_ASSERT(smi1 == smi2)

    std::string smi3 = ChemicalReactionToRxnSmiles(*rxn);
    ChemicalReaction *rxn2 = RxnMolToChemicalReaction(*mol);
    std::string smi4 = ChemicalReactionToRxnSmiles(*rxn2);

    TEST_ASSERT(smi3 == smi4)

    delete rxn;
    delete rxn2;
    delete m;
    delete mol;
  }
  {
    // test molecule with bond between reactant and product, same test as having
    // different rxn roles in one fragment
    std::string fName;
    fName = rdbase + "rxn1_1.mol";
    ROMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 18);
    TEST_ASSERT(m->getNumBonds() == 17);

    ChemicalReaction *rxn = RxnMolToChemicalReaction(*m);

    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 0);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 0);

    delete rxn;
    delete m;
  }
  {
    // test for molecule with correct rxn role for only one reactant
    std::string fName;
    fName = rdbase + "rxn1_2.mol";
    ROMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 18);
    TEST_ASSERT(m->getNumBonds() == 16);

    ChemicalReaction *rxn = RxnMolToChemicalReaction(*m);

    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 0);

    delete rxn;
    delete m;
    ;
  }
  {
    // test for molecule with correct rxn role for only one reactant
    std::string smi = "[C:1][C:2]>>[C:1].[C:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);

    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 2);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 0);

    ROMol *mol = ChemicalReactionToRxnMol(*rxn);

    TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp(common_properties::molRxnRole));
    TEST_ASSERT(mol->getAtomWithIdx(0)->getProp<int>(
                    common_properties::molRxnRole) == 1);
    TEST_ASSERT(mol->getAtomWithIdx(2)->hasProp(common_properties::molRxnRole));
    TEST_ASSERT(mol->getAtomWithIdx(2)->getProp<int>(
                    common_properties::molRxnRole) == 2);
    TEST_ASSERT(mol->getAtomWithIdx(3)->hasProp(common_properties::molRxnRole));
    TEST_ASSERT(mol->getAtomWithIdx(3)->getProp<int>(
                    common_properties::molRxnRole) == 2);
    delete rxn;
    delete mol;

    smi = "[C:1][C:2]>[Na]>[C:1].[C:2]";
    rxn = RxnSmartsToChemicalReaction(smi);

    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 2);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 1);

    mol = ChemicalReactionToRxnMol(*rxn);

    TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp(common_properties::molRxnRole));
    TEST_ASSERT(mol->getAtomWithIdx(0)->getProp<int>(
                    common_properties::molRxnRole) == 1);
    TEST_ASSERT(mol->getAtomWithIdx(2)->hasProp(common_properties::molRxnRole));
    TEST_ASSERT(mol->getAtomWithIdx(2)->getProp<int>(
                    common_properties::molRxnRole) == 2);
    TEST_ASSERT(mol->getAtomWithIdx(4)->hasProp(common_properties::molRxnRole));
    TEST_ASSERT(mol->getAtomWithIdx(4)->getProp<int>(
                    common_properties::molRxnRole) == 3);
    delete rxn;
    delete mol;

    smi = "[C:1][C:2]>>";
    rxn = RxnSmartsToChemicalReaction(smi);

    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 0);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 0);

    mol = ChemicalReactionToRxnMol(*rxn);

    TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp(common_properties::molRxnRole));
    TEST_ASSERT(mol->getAtomWithIdx(0)->getProp<int>(
                    common_properties::molRxnRole) == 1);
    delete rxn;
    delete mol;

    smi = ">>[C:1].[C:2]";
    rxn = RxnSmartsToChemicalReaction(smi);

    TEST_ASSERT(rxn->getNumReactantTemplates() == 0);
    TEST_ASSERT(rxn->getNumProductTemplates() == 2);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 0);

    mol = ChemicalReactionToRxnMol(*rxn);

    TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp(common_properties::molRxnRole));
    TEST_ASSERT(mol->getAtomWithIdx(0)->getProp<int>(
                    common_properties::molRxnRole) == 2);
    TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::molRxnRole));
    TEST_ASSERT(mol->getAtomWithIdx(1)->getProp<int>(
                    common_properties::molRxnRole) == 2);

    delete rxn;
    delete mol;
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test48ParensInProducts1() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing grouping parens in products1" << std::endl;

  {
    std::string smi = "[C:1][C:2]>>[C:1].[C:2]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 2);
    delete rxn;
  }

  {
    std::string smi = "[C:1][C:2]>>([C:1].[C:2])";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    delete rxn;
  }
  {
    std::string smi = "([C:1].C(=O)O).[C:2]>>([C:1].[C:2]).C(=O)O";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 2);
    delete rxn;
  }

  {
    std::string smi = "[C:1](=O)O.[C:2]>>[C:1][C:2].C(=O)O";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 2);
    delete rxn;
  }

  {
    std::string smi = "[C:1](=O)O.([C:2])>>[C:1][C:2].(C(=O)O)";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 2);
    delete rxn;
  }
  {
    std::string smi = "[C:1](=O)O.([C:2].N)>>[C:1][C:2].(C(=O)O.N)";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 2);
    delete rxn;
  }
  {
    std::string smi = "[C:1](=O)O.[C:2]>>([C:1][C:2]";
    ChemicalReaction *rxn = nullptr;
    try {
      rxn = RxnSmartsToChemicalReaction(smi);
      TEST_ASSERT(!rxn);
    } catch (const ChemicalReactionParserException &) {
      ;
    }
    delete rxn;
  }
  {
    std::string smi = "[C:1](=O)O.[C:2]>>[C:1].([C:2]";
    ChemicalReaction *rxn = nullptr;
    try {
      rxn = RxnSmartsToChemicalReaction(smi);
      TEST_ASSERT(!rxn);
    } catch (const ChemicalReactionParserException &) {
      ;
    }
    delete rxn;
  }
  {
    std::string smi = "[C:1](=O)O.[C:2]>>[C:1][C:2])";
    ChemicalReaction *rxn = nullptr;
    try {
      rxn = RxnSmartsToChemicalReaction(smi);
      TEST_ASSERT(!rxn);
    } catch (const ChemicalReactionParserException &) {
      ;
    }
    delete rxn;
  }
  {
    std::string smi = "[C:1](=O)O.[C:2]>>[C:1]).[C:2]";
    ChemicalReaction *rxn = nullptr;
    try {
      rxn = RxnSmartsToChemicalReaction(smi);
      TEST_ASSERT(!rxn);
    } catch (const ChemicalReactionParserException &) {
      ;
    }
    delete rxn;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test49ParensInProducts2() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing grouping parens in products2" << std::endl;

  {
    std::string smi = "[C:1][O:2]>>([C:1].[O:2])";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();
    MOL_SPTR_VECT reacts;
    reacts.clear();
    smi = "C1NO1";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    std::vector<MOL_SPTR_VECT> prods;
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 3);
    TEST_ASSERT(prods[0][0]->getNumBonds() == 2);

    delete rxn;
  }
  {
    std::string smi = "([N:1].[O:2])>>([N:1]C.[O:2]CCC)";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();
    MOL_SPTR_VECT reacts;
    reacts.clear();
    smi = "Nn1ccc(O)c1";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    TEST_ASSERT(reacts.size() == 1);
    TEST_ASSERT(reacts[0]->getNumAtoms() == 7);
    TEST_ASSERT(reacts[0]->getNumBonds() == 7);
    std::vector<MOL_SPTR_VECT> prods;
    prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 11);
    TEST_ASSERT(prods[0][0]->getNumBonds() == 11);

    smi = "CCCOc1ccn(NC)c1";
    auto *p00 = static_cast<RWMol *>(prods[0][0].get());
    MolOps::sanitizeMol(*p00);
    TEST_ASSERT(MolToSmiles(*p00) == smi);

    delete rxn;
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test50RNXFileParserWithEmptyAgentColumn() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing RNX file parser with empty agent column"
                       << std::endl;

  {
    std::string rdbase = getenv("RDBASE");
    std::string fName;

    fName = rdbase + "/Code/GraphMol/ChemReactions/testData/rxn2.mol";

    ChemicalReaction *rxn = RxnFileToChemicalReaction(fName);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 0);

    delete rxn;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test51RNXSmilesFromPatentData() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing RNX from patent data" << std::endl;
  unsigned int nWarn, nError;
  {
    // product atom-mapping numbers found multiple times, validation should fail
    std::string smi =
        "[Na+].[Na+].[NH2:23][CH2:22][CH2:21][CH2:20][CH2:19][CH2:18][CH2:17]["
        "CH2:16][CH2:15][CH2:14][NH2:13].O=C([O-])[O-].[Cl:12][c:8]1[cH:9][cH:"
        "10][c:11]2[c:2](Cl)[n:3][cH:4][n:5][c:6]2[cH:7]1>CC(C)O>[Cl:12][c:8]1["
        "cH:9][cH:10][c:11]2[c:6]([cH:7]1)[n:5][cH:4][n:3][c:2]2[NH:23][CH2:22]"
        "[CH2:21][CH2:20][CH2:19][CH2:18][CH2:17][CH2:16][CH2:15][CH2:14][NH:"
        "13][c:2]1[n:3][cH:4][n:5][c:6]2[cH:7][c:8]([Cl:12])[cH:9][cH:10][c:11]"
        "21";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn != 0);
    TEST_ASSERT(nError == 0);

    delete rxn;
  }
  {
    // test removing reactant without heavyatoms to agents
    std::string smi =
        "[S].[H][H].[CH3:2][CH2:1][O:3][C:4](=[O:20])[c:5]1[c:6]([O:17][CH2:18]"
        "[CH3:19])[cH:7][c:8]([C:14](=[O:15])Cl)[cH:9][c:10]1[O:11][CH2:12]["
        "CH3:13].c1ccc2ncccc2c1>[Pd+2].[Ba+2].O=S(=O)([O-])[O-].O=S(=O)([O-])["
        "O-]>[CH3:2][CH2:1][O:3][C:4](=[O:20])[c:5]1[c:6]([O:17][CH2:18][CH3:"
        "19])[cH:7][c:8]([CH:14]=[O:15])[cH:9][c:10]1[O:11][CH2:12][CH3:13]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    rxn->initReactantMatchers();
    rxn->removeUnmappedReactantTemplates();
    rxn->removeUnmappedProductTemplates();
    delete rxn;
  }
  {
    // reactant atom-mapping numbers found multiple times, validation should
    // fail
    std::string smi =
        "[BrH:1].[CH3:24][CH2:23][CH2:22][CH2:21][CH2:20][C:19](=[O:25])O[C:19]"
        "(=[O:25])[CH2:20][CH2:21][CH2:22][CH2:23][CH3:24].[CH3:16][CH:15]1[N:"
        "14]([CH3:17])[CH:13]2[CH2:18][C:9]1([c:5]1[cH:6][cH:7][cH:8][c:3]([OH:"
        "2])[cH:4]1)[CH2:10][CH2:11][CH2:12]2>c1ccncc1>[BrH:1].[CH3:24][CH2:23]"
        "[CH2:22][CH2:21][CH2:20][C:19](=[O:25])[O:2][c:3]1[cH:8][cH:7][cH:6]["
        "c:5]([C:9]23[CH2:18][CH:13]([CH2:12][CH2:11][CH2:10]2)[N:14]([CH3:17])"
        "[CH:15]3[CH3:16])[cH:4]1";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(!rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError != 0);

    delete rxn;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test52RedundantProductMappingNumbersAndRunReactants() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing products redundant atom mapping numbers in run_reactants"
      << std::endl;
  unsigned int nWarn, nError;
  {
    std::string smi = "[C:1]-[OH:2]>>[C:1]-[O:2]-[C:1]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn != 0);
    TEST_ASSERT(nError == 0);
    rxn->initReactantMatchers();

    MOL_SPTR_VECT reacts;
    reacts.clear();
    smi = "N[13CH2]O";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    TEST_ASSERT(reacts.size() == 1);
    std::vector<MOL_SPTR_VECT> prods;
    prods = rxn->runReactants(reacts);

    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);

    smi = "N[13CH2]O[13CH2]N";
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == smi);

    delete rxn;
  }
  {
    std::string smi = "[C:1]-[OH:2]>>[C:1]-[O:2]-[C:1]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn != 0);
    TEST_ASSERT(nError == 0);
    rxn->initReactantMatchers();

    MOL_SPTR_VECT reacts;
    reacts.clear();
    smi = "[13CH3]O";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    TEST_ASSERT(reacts.size() == 1);
    std::vector<MOL_SPTR_VECT> prods;
    prods = rxn->runReactants(reacts);

    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getNumAtoms() == 3);
    TEST_ASSERT(prods[0][0]->getNumBonds() == 2);

    smi = "[13CH3]O[13CH3]";
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == smi);

    delete rxn;
  }
  {
    std::string smi = "[O:1]-[C:2]>>[O:1]-[C:2]-C-[C:2]-[O:1]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn != 0);
    TEST_ASSERT(nError == 0);
    rxn->initReactantMatchers();

    MOL_SPTR_VECT reacts;
    reacts.clear();
    smi = "[13C]1OCN1";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    TEST_ASSERT(reacts.size() == 1);
    std::vector<MOL_SPTR_VECT> prods;
    prods = rxn->runReactants(reacts);

    TEST_ASSERT(prods.size() == 2);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[1].size() == 1);

    smi = "C1N[13C](C[13C]2NCO2)O1";
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == smi);
    smi = "C(C1N[13C]O1)C1N[13C]O1";
    TEST_ASSERT(MolToSmiles(*prods[1][0], true) == smi);

    delete rxn;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test53ReactionSubstructureMatching() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing reaction substructure matching w/o agents"
                       << std::endl;
  {
    std::string smi = "c1ccccc1>>c1ccncc1";
    ChemicalReaction *query_rxn =
        RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(query_rxn);
    TEST_ASSERT(query_rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(query_rxn->getNumProductTemplates() == 1);
    TEST_ASSERT(query_rxn->getNumAgentTemplates() == 0);

    smi = "CCC.c1ccccc1>>CCC.c1ccncc1";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 2);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 0);
    TEST_ASSERT(hasReactionSubstructMatch(*rxn, *query_rxn, false));
    TEST_ASSERT(hasReactionSubstructMatch(*rxn, *query_rxn, true));
    delete rxn;

    smi = "c1ccccc1.CCC>>CCC.c1ccncc1";
    rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(hasReactionSubstructMatch(*rxn, *query_rxn, false));
    TEST_ASSERT(hasReactionSubstructMatch(*rxn, *query_rxn, true));
    delete rxn;

    smi = "c1ccccc1.CCC>>c1ccncc1.CCC";
    rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(hasReactionSubstructMatch(*rxn, *query_rxn, false));
    TEST_ASSERT(hasReactionSubstructMatch(*rxn, *query_rxn, true));
    delete rxn;

    smi = "CCC.c1ccccc1.C(=O)O>>CCC.c1ccncc1.C(=O)O";
    rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 3);
    TEST_ASSERT(rxn->getNumProductTemplates() == 3);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 0);
    TEST_ASSERT(hasReactionSubstructMatch(*rxn, *query_rxn, false));
    TEST_ASSERT(hasReactionSubstructMatch(*rxn, *query_rxn, true));
    delete rxn;

    smi = "CCC.C(=O)O.c1ccccc1>>CCC.c1ccncc1.C(=O)O";
    rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(hasReactionSubstructMatch(*rxn, *query_rxn, false));
    TEST_ASSERT(hasReactionSubstructMatch(*rxn, *query_rxn, true));
    delete rxn;

    smi = "CCC.C(=O)O.c1ccccc1>>CCC.C(=O)O.c1ccncc1";
    rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(hasReactionSubstructMatch(*rxn, *query_rxn, false));
    TEST_ASSERT(hasReactionSubstructMatch(*rxn, *query_rxn, true));
    delete rxn;

    smi = "CCC.c1ccccc1.C(=O)O>NC=O>CCC.c1ccncc1.C(=O)O";
    rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 3);
    TEST_ASSERT(rxn->getNumProductTemplates() == 3);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 1);
    TEST_ASSERT(hasReactionSubstructMatch(*rxn, *query_rxn, false));
    TEST_ASSERT(hasReactionSubstructMatch(*rxn, *query_rxn, true));

    delete query_rxn;
    smi = "c1ccccc1>NC=O>c1ccncc1";
    query_rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(query_rxn);
    TEST_ASSERT(query_rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(query_rxn->getNumProductTemplates() == 1);
    TEST_ASSERT(query_rxn->getNumAgentTemplates() == 1);
    TEST_ASSERT(hasReactionSubstructMatch(*rxn, *query_rxn, false));
    TEST_ASSERT(hasReactionSubstructMatch(*rxn, *query_rxn, true));
    delete rxn;

    smi = "CCC.c1ccccc1.C(=O)O>>CCC.c1ccncc1.C(=O)O";
    rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 3);
    TEST_ASSERT(rxn->getNumProductTemplates() == 3);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 0);
    TEST_ASSERT(hasReactionSubstructMatch(*rxn, *query_rxn, false));
    TEST_ASSERT(!hasReactionSubstructMatch(*rxn, *query_rxn, true));
    delete rxn;

    smi = "CCC.c1ccccc1.C(=O)O>>CCC.c1ncncc1.C(=O)O";
    rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 3);
    TEST_ASSERT(rxn->getNumProductTemplates() == 3);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 0);
    TEST_ASSERT(!hasReactionSubstructMatch(*rxn, *query_rxn, false));
    TEST_ASSERT(!hasReactionSubstructMatch(*rxn, *query_rxn, true));

    delete rxn;
    delete query_rxn;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test54RedundantProductMappingNumbersAndRSChirality() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing products with redundant atom mapping numbers and chirality"
      << std::endl;
  unsigned int nWarn, nError;
  std::string smi, cip;
  {
    // perserve the stereo chemistry of the reactant in the product
    smi = "[C:1][O:2]>>[C:1][O:2]N[O:2][C:1]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn != 0);
    TEST_ASSERT(nError == 0);
    rxn->initReactantMatchers();

    MOL_SPTR_VECT reacts;
    reacts.clear();
    smi = "[OH][C@](F)(Cl)Br";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    MolOps::assignStereochemistry(*mol);
    TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    reacts.push_back(ROMOL_SPTR(mol));
    TEST_ASSERT(reacts.size() == 1);
    std::vector<MOL_SPTR_VECT> prods;
    prods = rxn->runReactants(reacts);

    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);

    std::cout << MolToSmiles(*prods[0][0], true) << std::endl;
    smi = "F[C@@](Cl)(Br)ONO[C@](F)(Cl)Br";
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == smi);

    ROMOL_SPTR prod = prods[0][0];
    MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));
    MolOps::assignStereochemistry(*prod);
    TEST_ASSERT(prod->getAtomWithIdx(0)->getAtomicNum() == 6);
    TEST_ASSERT(prod->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
    prod->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    TEST_ASSERT(prod->getAtomWithIdx(4)->getAtomicNum() == 6);
    TEST_ASSERT(prod->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));
    prod->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    delete rxn;
  }
  {
    // invert the stereo chemistry of one carbon of the reactant in the product
    smi = "[C@:1][O:2]>>[C@:1][O:2]N[O:2][C@@:1]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn != 0);
    TEST_ASSERT(nError == 0);
    rxn->initReactantMatchers();

    MOL_SPTR_VECT reacts;
    reacts.clear();
    smi = "[OH][C@](F)(Cl)Br";
    ROMol *mol = SmilesToMol(smi);
    MolOps::assignStereochemistry(*mol);
    TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    TEST_ASSERT(reacts.size() == 1);
    std::vector<MOL_SPTR_VECT> prods;
    prods = rxn->runReactants(reacts);

    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);

    std::cout << MolToSmiles(*prods[0][0], true) << std::endl;
    smi = "F[C@](Cl)(Br)ONO[C@](F)(Cl)Br";
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == smi);

    ROMOL_SPTR prod = prods[0][0];
    MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));
    MolOps::assignStereochemistry(*prod);
    TEST_ASSERT(prod->getAtomWithIdx(0)->getAtomicNum() == 6);
    TEST_ASSERT(prod->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
    prod->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    TEST_ASSERT(prod->getAtomWithIdx(4)->getAtomicNum() == 6);
    TEST_ASSERT(prod->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));
    prod->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");
    delete rxn;
  }
  {
    // Both carbons in the product (mapped with #2), should be (S)
    smi =
        "[OH:5][C@:2]([F:3])([Cl:1])[Br:4]>>[F:3][C@@:2]([Cl:1])([Br:4])[O:5]N["
        "O:5][C@:2]([F:3])([Cl:1])[Br:4]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn != 0);
    TEST_ASSERT(nError == 0);

    MOL_SPTR_VECT reacts;
    rxn->initReactantMatchers();
    reacts.clear();

    smi = "[OH][C@](F)(Cl)Br";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    MolOps::assignStereochemistry(*mol);
    TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    reacts.push_back(ROMOL_SPTR(mol));
    TEST_ASSERT(reacts.size() == 1);
    std::vector<MOL_SPTR_VECT> prods;
    prods = rxn->runReactants(reacts);

    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);

    std::cout << MolToSmiles(*prods[0][0], true) << std::endl;
    smi = "F[C@@](Cl)(Br)ONO[C@](F)(Cl)Br";
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == smi);

    ROMOL_SPTR prod = prods[0][0];
    MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));
    MolOps::assignStereochemistry(*prod);
    TEST_ASSERT(prod->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    prod->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    TEST_ASSERT(prod->getAtomWithIdx(7)->hasProp(common_properties::_CIPCode));
    prod->getAtomWithIdx(7)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    delete rxn;
  }
  {
    // One carbon in the product (mapped with #2) should be (S), the other (R)
    smi =
        "[OH:5][C@:2]([F:3])([Cl:1])[Br:4]>>[F:3][C@@:2]([Cl:1])([Br:4])[O:5]N["
        "O:5][C@@:2]([F:3])([Cl:1])[Br:4]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn != 0);
    TEST_ASSERT(nError == 0);

    MOL_SPTR_VECT reacts;
    rxn->initReactantMatchers();
    reacts.clear();

    smi = "[OH][C@](F)(Cl)Br";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    MolOps::assignStereochemistry(*mol);
    TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    reacts.push_back(ROMOL_SPTR(mol));
    TEST_ASSERT(reacts.size() == 1);
    std::vector<MOL_SPTR_VECT> prods;
    prods = rxn->runReactants(reacts);

    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);

    std::cout << MolToSmiles(*prods[0][0], true) << std::endl;
    smi = "F[C@](Cl)(Br)ONO[C@](F)(Cl)Br";
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == smi);

    ROMOL_SPTR prod = prods[0][0];
    MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));
    MolOps::assignStereochemistry(*prod);
    TEST_ASSERT(prod->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    prod->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    TEST_ASSERT(prod->getAtomWithIdx(7)->hasProp(common_properties::_CIPCode));
    prod->getAtomWithIdx(7)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");

    delete rxn;
  }
  {
    // Both carbons in the product (mapped with #2), should be (S)
    smi =
        "[OH:5][C@:2]([F:3])([Cl:1])[Br:4]>>[F:3][C@@:2]([Cl:1])([Br:4])[O:5]["
        "C@:2]([F:3])([Cl:1])[Br:4]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn != 0);
    TEST_ASSERT(nError == 0);

    MOL_SPTR_VECT reacts;
    rxn->initReactantMatchers();
    reacts.clear();

    smi = "[OH][C@](F)(Cl)Br";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    MolOps::assignStereochemistry(*mol);
    TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    reacts.push_back(ROMOL_SPTR(mol));
    TEST_ASSERT(reacts.size() == 1);
    std::vector<MOL_SPTR_VECT> prods;
    prods = rxn->runReactants(reacts);

    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);

    std::cout << MolToSmiles(*prods[0][0], true) << std::endl;
    smi = "F[C@@](Cl)(Br)O[C@](F)(Cl)Br";
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == smi);

    ROMOL_SPTR prod = prods[0][0];
    MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));
    MolOps::assignStereochemistry(*prod);
    TEST_ASSERT(prod->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    prod->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    TEST_ASSERT(prod->getAtomWithIdx(5)->hasProp(common_properties::_CIPCode));
    prod->getAtomWithIdx(5)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    delete rxn;
  }
  {
    // One carbon in the product (mapped with #2) should be (S), the other (R)
    smi =
        "[OH:5][C@:2]([F:3])([Cl:1])[Br:4]>>[F:3][C@@:2]([Cl:1])([Br:4])[O:5]["
        "C@@:2]([F:3])([Cl:1])[Br:4]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn != 0);
    TEST_ASSERT(nError == 0);

    MOL_SPTR_VECT reacts;
    rxn->initReactantMatchers();
    reacts.clear();

    smi = "[OH][C@](F)(Cl)Br";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    MolOps::assignStereochemistry(*mol);
    TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    reacts.push_back(ROMOL_SPTR(mol));
    TEST_ASSERT(reacts.size() == 1);
    std::vector<MOL_SPTR_VECT> prods;
    prods = rxn->runReactants(reacts);

    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);

    std::cout << MolToSmiles(*prods[0][0], true) << std::endl;
    smi = "F[C@](Cl)(Br)O[C@](F)(Cl)Br";
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == smi);

    ROMOL_SPTR prod = prods[0][0];
    MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));
    MolOps::assignStereochemistry(*prod);
    TEST_ASSERT(prod->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    prod->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    TEST_ASSERT(prod->getAtomWithIdx(5)->hasProp(common_properties::_CIPCode));
    prod->getAtomWithIdx(5)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");

    delete rxn;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test55RedundantProductMappingNumbersAndEZStereochemistry() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing products with redundant atom mapping numbers and chirality"
      << std::endl;
  unsigned int nWarn, nError;
  std::string smi;
  {
    // both double bonds in the product are (E)
    smi =
        "[CH3:1]\\[CH:2]=[CH:3]\\[CH3:4]>>[CH3:1]\\[CH:2]=[CH:3]\\[CH2:4][CH2:"
        "4]\\[CH:3]=[CH:2]\\[CH3:1]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn != 0);
    TEST_ASSERT(nError == 0);

    MOL_SPTR_VECT reacts;
    rxn->initReactantMatchers();
    reacts.clear();

    smi = "[CH3]\\[CH]=[CH]\\[CH3]";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    MolOps::assignStereochemistry(*mol);
    TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
    reacts.push_back(ROMOL_SPTR(mol));
    TEST_ASSERT(reacts.size() == 1);
    std::vector<MOL_SPTR_VECT> prods;
    prods = rxn->runReactants(reacts);

    TEST_ASSERT(prods.size() == 2);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[1].size() == 1);

    smi = "C/C=C/CC/C=C/C";
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == smi);
    TEST_ASSERT(MolToSmiles(*prods[1][0], true) == smi);

    ROMOL_SPTR prod = prods[0][0];
    MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));

    Bond *stereoBond = prod->getBondWithIdx(1);
    INT_VECT stereoAtomsRef{{0, 3}};
    TEST_ASSERT(stereoBond->getStereoAtoms() == stereoAtomsRef);
    TEST_ASSERT(stereoBond->getStereo() == Bond::STEREOTRANS);

    stereoBond = prod->getBondWithIdx(5);
    stereoAtomsRef = {4, 7};
    TEST_ASSERT(stereoBond->getStereoAtoms() == stereoAtomsRef);
    TEST_ASSERT(stereoBond->getStereo() == Bond::STEREOTRANS);

    prod = prods[1][0];
    MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));

    stereoBond = prod->getBondWithIdx(1);
    stereoAtomsRef = {0, 3};
    TEST_ASSERT(stereoBond->getStereoAtoms() == stereoAtomsRef);
    TEST_ASSERT(stereoBond->getStereo() == Bond::STEREOTRANS);

    stereoBond = prod->getBondWithIdx(5);
    stereoAtomsRef = {4, 7};
    TEST_ASSERT(stereoBond->getStereoAtoms() == stereoAtomsRef);
    TEST_ASSERT(stereoBond->getStereo() == Bond::STEREOTRANS);

    delete rxn;
  }
  {
    // both double bonds in the product are (E)
    smi =
        "[CH3:1]\\[CH:2]=[CH:3]\\[CH3:4]>>[CH3:1]\\[CH:2]=[CH:3]\\[CH2:4][CH2:"
        "4]\\[CH:3]=[CH:2]/[CH3:1]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn != 0);
    TEST_ASSERT(nError == 0);

    MOL_SPTR_VECT reacts;
    rxn->initReactantMatchers();
    reacts.clear();

    smi = "[CH3]\\[CH]=[CH]\\[CH3]";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    MolOps::assignStereochemistry(*mol);
    TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
    reacts.push_back(ROMOL_SPTR(mol));
    TEST_ASSERT(reacts.size() == 1);
    std::vector<MOL_SPTR_VECT> prods;
    prods = rxn->runReactants(reacts);

    TEST_ASSERT(prods.size() == 2);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[1].size() == 1);

    smi = "C/C=C\\CC/C=C/C";
    // std::cerr<<MolToSmiles(*prods[0][0],true)<<std::endl;
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == smi);
    TEST_ASSERT(MolToSmiles(*prods[1][0], true) == smi);

    ROMOL_SPTR prod = prods[0][0];
    MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));

    Bond *stereoBond = prod->getBondWithIdx(1);
    INT_VECT stereoAtomsRef{{0, 3}};
    TEST_ASSERT(stereoBond->getStereoAtoms() == stereoAtomsRef);
    TEST_ASSERT(stereoBond->getStereo() == Bond::STEREOTRANS);

    stereoBond = prod->getBondWithIdx(5);
    stereoAtomsRef = {4, 7};
    TEST_ASSERT(stereoBond->getStereoAtoms() == stereoAtomsRef);
    TEST_ASSERT(stereoBond->getStereo() == Bond::STEREOCIS);

    prod = prods[1][0];
    MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));

    stereoBond = prod->getBondWithIdx(1);
    stereoAtomsRef = {0, 3};
    TEST_ASSERT(stereoBond->getStereoAtoms() == stereoAtomsRef);
    TEST_ASSERT(stereoBond->getStereo() == Bond::STEREOTRANS);

    stereoBond = prod->getBondWithIdx(5);
    stereoAtomsRef = {4, 7};
    TEST_ASSERT(stereoBond->getStereoAtoms() == stereoAtomsRef);
    TEST_ASSERT(stereoBond->getStereo() == Bond::STEREOCIS);

    delete rxn;
  }
  {
    // both double bonds in the product are (E)
    smi =
        "[CH3:1]\\[CH:2]=[CH:3]\\[CH3:4]>>[CH3:1]\\[CH:2]=[CH:3]\\[CH:3]=[CH:2]"
        "\\[CH3:1]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi, nullptr, true);
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn != 0);
    TEST_ASSERT(nError == 0);

    MOL_SPTR_VECT reacts;
    rxn->initReactantMatchers();
    reacts.clear();

    smi = "[CH3]\\[CH]=[CH]\\[CH3]";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    MolOps::assignStereochemistry(*mol);
    TEST_ASSERT(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
    reacts.push_back(ROMOL_SPTR(mol));
    TEST_ASSERT(reacts.size() == 1);
    std::vector<MOL_SPTR_VECT> prods;
    prods = rxn->runReactants(reacts);

    TEST_ASSERT(prods.size() == 2);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[1].size() == 1);

    smi = "C/C=C/C=C/C";
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == smi);
    TEST_ASSERT(MolToSmiles(*prods[1][0], true) == smi);

    ROMOL_SPTR prod = prods[0][0];
    MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));

    Bond *stereoBond = prod->getBondWithIdx(1);
    INT_VECT stereoAtomsRef{{0, 3}};
    TEST_ASSERT(stereoBond->getStereoAtoms() == stereoAtomsRef);
    TEST_ASSERT(stereoBond->getStereo() == Bond::STEREOTRANS);

    stereoBond = prod->getBondWithIdx(3);
    stereoAtomsRef = {2, 5};
    TEST_ASSERT(stereoBond->getStereoAtoms() == stereoAtomsRef);
    TEST_ASSERT(stereoBond->getStereo() == Bond::STEREOTRANS);

    prod = prods[1][0];
    MolOps::sanitizeMol(*(static_cast<RWMol *>(prod.get())));

    stereoBond = prod->getBondWithIdx(1);
    stereoAtomsRef = {0, 3};
    TEST_ASSERT(stereoBond->getStereoAtoms() == stereoAtomsRef);
    TEST_ASSERT(stereoBond->getStereo() == Bond::STEREOTRANS);

    stereoBond = prod->getBondWithIdx(3);
    stereoAtomsRef = {2, 5};
    TEST_ASSERT(stereoBond->getStereoAtoms() == stereoAtomsRef);
    TEST_ASSERT(stereoBond->getStereo() == Bond::STEREOTRANS);

    delete rxn;
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test56TestOldPickleVersion() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing old pickle file with new agent version"
                       << std::endl;
  unsigned int nWarn, nError;
  {
    std::string pklName = getenv("RDBASE");
    pklName += "/Code/GraphMol/ChemReactions/testData/testpickle.bin";
    std::ifstream inStream(pklName.c_str(), std::ios_base::binary);

    auto *rxn = new ChemicalReaction();
    ReactionPickler::reactionFromPickle(inStream, rxn);
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);
    rxn->initReactantMatchers();

    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);

    delete rxn;
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test57IntroductionOfNewChiralCenters() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing introduction of new atoms with chirality"
                       << std::endl;

  {  // a reaction that induces new atoms with stereochem
    std::string smi =
        "[F:1][C:2]([Cl:3])([Br:4])[I:5]>>[F:1][C@@:2]([Cl:3])([N:4][C@:6]([Cl:"
        "7])[Br:8])[I:5]";

    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();

    MOL_SPTR_VECT reacts;
    reacts.clear();
    smi = "FC(Cl)(Br)I";
    ROMol *mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    std::vector<MOL_SPTR_VECT> prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    BOOST_LOG(rdInfoLog) << MolToSmiles(*prods[0][0], true) << std::endl;
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "F[C@](Cl)(I)N[C@H](Cl)Br");

    delete rxn;
  }

  {  // a reaction that induces new atoms with stereochem
    std::string rdbase = getenv("RDBASE");
    std::string fName;

    fName =
        rdbase + "/Code/GraphMol/ChemReactions/testData/testRXNChirality.rxn";
    ChemicalReaction *rxn = RxnFileToChemicalReaction(fName);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    BOOST_LOG(rdInfoLog) << ChemicalReactionToRxnSmiles(*rxn) << std::endl;
    rxn->initReactantMatchers();

    MOL_SPTR_VECT reacts;
    reacts.clear();
    fName =
        rdbase + "/Code/GraphMol/ChemReactions/testData/testRXNChirality1.sdf";
    ROMol *mol = MolFileToMol(fName);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    std::vector<MOL_SPTR_VECT> prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    BOOST_LOG(rdInfoLog) << MolToSmiles(*prods[0][0], true) << std::endl;
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "C[C@H](F)CCC[C@@H](C)Cl");

    delete rxn;
  }

  {  // a reaction that induces new atoms with stereochem
    std::string rdbase = getenv("RDBASE");
    std::string fName;

    fName =
        rdbase + "/Code/GraphMol/ChemReactions/testData/testRXNChirality2.sdf";
    ROMol *mol = MolFileToMol(fName);
    auto *rxn = new ChemicalReaction();
    rxn->addReactantTemplate(ROMOL_SPTR(mol));
    fName =
        rdbase + "/Code/GraphMol/ChemReactions/testData/testRXNChirality3.sdf";
    mol = MolFileToMol(fName);
    rxn->addProductTemplate(ROMOL_SPTR(mol));
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    BOOST_LOG(rdInfoLog) << ChemicalReactionToRxnSmiles(*rxn) << std::endl;
    rxn->initReactantMatchers();
    updateProductsStereochem(rxn);

    MOL_SPTR_VECT reacts;
    reacts.clear();
    fName =
        rdbase + "/Code/GraphMol/ChemReactions/testData/testRXNChirality1.sdf";
    mol = MolFileToMol(fName);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));
    std::vector<MOL_SPTR_VECT> prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    BOOST_LOG(rdInfoLog) << MolToSmiles(*prods[0][0], true) << std::endl;
    TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "C[C@H](F)CCC[C@@H](C)Cl");

    delete rxn;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test58MolFileValueRoundTrip() {
  ChemicalReaction *rxn;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing round trip molFileValue" << std::endl;

  const char *rxnB =
      "$RXN\n"
      "\n"
      "      ISIS     090220091541\n"
      "\n"
      "  2  1\n"
      "$MOL\n"
      "\n"
      "  -ISIS-  09020915412D\n"
      "\n"
      "  3  2  0  0  0  0  0  0  0  0999 V2000\n"
      "   -2.9083   -0.4708    0.0000 R#  0  0  0  0  0  0  0  0  0  1  0  0\n"
      "   -2.3995   -0.1771    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0\n"
      "   -2.4042    0.4125    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "  1  2  1  0  0  0  0\n"
      "  2  3  2  0  0  0  0\n"
      "V    2 aldehyde\n"
      "M  RGP  1   1   1\n"
      "M  END\n"
      "$MOL\n"
      "\n"
      "  -ISIS-  09020915412D\n"
      "\n"
      "  2  1  0  0  0  0  0  0  0  0999 V2000\n"
      "    2.8375   -0.2500    0.0000 R#  0  0  0  0  0  0  0  0  0  3  0  0\n"
      "    3.3463    0.0438    0.0000 N   0  0  0  0  0  0  0  0  0  4  0  0\n"
      "  1  2  1  0  0  0  0\n"
      "V    2 aldehyde\n"
      "M  RGP  1   1   2\n"
      "M  END\n"
      "$MOL\n"
      "\n"
      "  -ISIS-  09020915412D\n"
      "\n"
      "  4  3  0  0  0  0  0  0  0  0999 V2000\n"
      "   13.3088    0.9436    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0\n"
      "   13.8206    1.2321    0.0000 R#  0  0  0  0  0  0  0  0  0  1  0  0\n"
      "   13.3028    0.3561    0.0000 N   0  0  0  0  0  0  0  0  0  4  0  0\n"
      "   12.7911    0.0676    0.0000 R#  0  0  0  0  0  0  0  0  0  3  0  0\n"
      "  1  3  1  0  0  0  0\n"
      "  1  2  1  0  0  0  0\n"
      "  3  4  1  0  0  0  0\n"
      "M  RGP  2   2   1   4   2\n"
      "M  END";

  rxn = RxnBlockToChemicalReaction(rxnB);
  // check the mol file values
  for (MOL_SPTR_VECT::const_iterator template_mol =
           rxn->beginReactantTemplates();
       template_mol != rxn->endReactantTemplates(); ++template_mol) {
    const Atom *at = (*template_mol)->getAtomWithIdx(1);
    TEST_ASSERT(at->hasProp(common_properties::molFileValue));
    TEST_ASSERT(at->getProp<std::string>(common_properties::molFileValue) ==
                "aldehyde");
  }

  ChemicalReaction *rxn2 =
      RxnBlockToChemicalReaction(ChemicalReactionToRxnBlock(*rxn));

  for (MOL_SPTR_VECT::const_iterator template_mol =
           rxn2->beginReactantTemplates();
       template_mol != rxn2->endReactantTemplates(); ++template_mol) {
    const Atom *at = (*template_mol)->getAtomWithIdx(1);
    TEST_ASSERT(at->hasProp(common_properties::molFileValue));
    TEST_ASSERT(at->getProp<std::string>(common_properties::molFileValue) ==
                "aldehyde");
  }

  delete rxn;
  delete rxn2;

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test59ReactionCanonicalization() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing canonicatization of reactions" << std::endl;

  {  // with mapping numbers
    std::string smi =
        "[K+].[NH4+].[C:2]([CH:4]([CH2:9][CH3:10])[CH2:5][C:6]([OH:8])=[O:7])#["
        "N:3]>Cl.[OH-]>[CH2:2]([CH:4]([CH2:9][CH3:10])[CH2:5][C:6]([OH:8])=[O:"
        "7])[NH2:3].[K+].[NH4+]";

    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 3);
    TEST_ASSERT(rxn->getNumProductTemplates() == 3);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 2);
    std::string rxnsmi = ChemicalReactionToRxnSmiles(*rxn);
    delete rxn;

    rxn = RxnSmartsToChemicalReaction(rxnsmi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 3);
    TEST_ASSERT(rxn->getNumProductTemplates() == 3);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 2);
    std::string rxnsmi2 = ChemicalReactionToRxnSmiles(*rxn);

    TEST_ASSERT(rxnsmi == rxnsmi2);

    delete rxn;
  }
  {  // without mapping numbers
    std::string smi =
        "[K+].[NH4+].CCC(CC(O)=O)C#N>Cl.[OH-]>CCC(CC(O)=O)CN.[K+].[NH4+]";

    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 3);
    TEST_ASSERT(rxn->getNumProductTemplates() == 3);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 2);
    std::string rxnsmi = ChemicalReactionToRxnSmiles(*rxn);
    delete rxn;

    rxn = RxnSmartsToChemicalReaction(rxnsmi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 3);
    TEST_ASSERT(rxn->getNumProductTemplates() == 3);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 2);
    std::string rxnsmi2 = ChemicalReactionToRxnSmiles(*rxn);

    TEST_ASSERT(rxnsmi == rxnsmi2);

    delete rxn;
  }
}

void test60RunSingleReactant() {
  const std::string smirks_thiourea =
      "[N;$(N-[#6]):3]=[C;$(C=S):1].[N;$(N[#6]);!$(N=*);!$([N-]);!$(N#*);!$(["
      "ND3]);!$([ND4]);!$(N[O,N]);!$(N[C,S]=[S,O,N]):2]>>[N:3]-[C:1]-[N+0:2]";
  ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smirks_thiourea);
  rxn->initReactantMatchers();
  const std::string smi1 = "C=CCN=C=S";
  const std::string smi2 = "NCc1ncc(Cl)cc1Br";

  ROMOL_SPTR reag1(SmilesToMol(smi1));
  ROMOL_SPTR reag2(SmilesToMol(smi2));

  std::vector<MOL_SPTR_VECT> prods;

  {
    ROMOL_SPTR expected_result(SmilesToMol("C=CCNC(N)=S"));
    ROMOL_SPTR expected_sidechain_result(SmilesToMol("[*:1]=S.[*:3]CC=C"));
    const bool isomericSmiles = true;
    std::string expected = MolToSmiles(*expected_result, isomericSmiles);
    std::string expected_sidechain =
        MolToSmiles(*expected_sidechain_result, isomericSmiles);

    prods = rxn->runReactant(reag1, 0);
    TEST_ASSERT(prods.size() > 0);

    for (auto &prod : prods) {
      for (size_t prodidx = 0; prodidx < prod.size(); prodidx++) {
        TEST_ASSERT(prodidx < 1);
        TEST_ASSERT(MolToSmiles(*prod[prodidx]) == expected);
        ROMol *sidechain = reduceProductToSideChains(prod[prodidx]);
        std::string smi = MolToSmiles(*sidechain, isomericSmiles);
        TEST_ASSERT(smi == expected_sidechain);
        delete sidechain;
      }
    }

    prods = rxn->runReactant(reag2, 0);
    TEST_ASSERT(prods.size() == 0)
  }

  {
    ROMOL_SPTR expected_result(SmilesToMol("NCNCc1ncc(Cl)cc1Br"));
    ROMOL_SPTR expected_sidechain_result(SmilesToMol("[*:2]Cc1ncc(Cl)cc1Br"));
    const bool isomericSmiles = true;
    std::string expected = MolToSmiles(*expected_result, isomericSmiles);
    std::string expected_sidechain =
        MolToSmiles(*expected_sidechain_result, isomericSmiles);

    prods = rxn->runReactant(reag2, 1);
    TEST_ASSERT(prods.size() > 0);

    for (auto &prod : prods) {
      for (size_t prodidx = 0; prodidx < prod.size(); prodidx++) {
        TEST_ASSERT(prodidx < 1);
        TEST_ASSERT(MolToSmiles(*prod[prodidx]) == expected);
        ROMol *sidechain = reduceProductToSideChains(prod[prodidx]);
        std::string smi = MolToSmiles(*sidechain, isomericSmiles);
        TEST_ASSERT(smi == expected_sidechain);
        delete sidechain;
      }
    }

    prods = rxn->runReactant(reag1, 1);
    TEST_ASSERT(prods.size() == 0)
  }

  delete rxn;
}

void test61Github685() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing github685: SMARTS reaction triggers "
                          "invariant violation on chiral compounds"
                       << std::endl;

  {
    std::string smi =
        "[C:1][C:2][C:3]([O:4])([N:5])[Cl:6]>>[Cl:6].[C:2][C:3]([O:4])([N:5])["
        "C:1]";

    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 2);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 0);
    rxn->initReactantMatchers();

    {
      MOL_SPTR_VECT reacts;
      reacts.clear();
      smi = "CCC(O)(N)Cl";
      RWMol *mol = SmilesToMol(smi);
      TEST_ASSERT(mol);
      reacts.push_back(ROMOL_SPTR(mol));
      std::vector<MOL_SPTR_VECT> prods = rxn->runReactants(reacts);
      TEST_ASSERT(prods.size() == 1);
      TEST_ASSERT(prods[0].size() == 2);
      TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "Cl");
      TEST_ASSERT(MolToSmiles(*prods[0][1], true) == "CC(C)(N)O");
    }
    {
      MOL_SPTR_VECT reacts;
      reacts.clear();
      smi = "CC[C@](O)(N)Cl";
      RWMol *mol = SmilesToMol(smi);
      TEST_ASSERT(mol);
      reacts.push_back(ROMOL_SPTR(mol));
      std::vector<MOL_SPTR_VECT> prods = rxn->runReactants(reacts);
      TEST_ASSERT(prods.size() == 1);
      TEST_ASSERT(prods[0].size() == 2);
      TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "Cl");
      TEST_ASSERT(MolToSmiles(*prods[0][1], true) ==
                  "CC(C)(N)O");  // we lose the chirality because there's an
                                 // unspecified bond in the products
    }
    delete rxn;
  }
  {  // the original example
    std::string smi =
        "[#1;D1R0:4][#8;H1D2R0:3][#6;H0D3R0:2](=[#8;H0D1R0:5])[#6;D4;H1,H0,H2:"
        "1]([#6:7])[#6,#1D1AR0,ClH0D1AR0,FH0D1AR0,BrH0D1AR0:6]>>[*:3]=[*:2]=[*:"
        "5].[*:4]-[*:1](-[*:6])-[*:7]";

    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 2);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 0);
    rxn->initReactantMatchers();

    {
      MOL_SPTR_VECT reacts;
      reacts.clear();
      smi = "CC(N)C(=O)O";
      RWMol *mol = SmilesToMol(smi);
      TEST_ASSERT(mol);
      MolOps::addHs(*mol);
      reacts.push_back(ROMOL_SPTR(mol));
      std::vector<MOL_SPTR_VECT> prods = rxn->runReactants(reacts);
      TEST_ASSERT(prods.size() == 1);
      TEST_ASSERT(prods[0].size() == 2);
      TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "O=C=O");
      TEST_ASSERT(MolToSmiles(*prods[0][1], true) ==
                  "[H]N([H])C([H])([H])C([H])([H])[H]");
    }
    {
      MOL_SPTR_VECT reacts;
      reacts.clear();
      smi = "C[C@@H](N)C(=O)O";
      RWMol *mol = SmilesToMol(smi);
      TEST_ASSERT(mol);
      MolOps::addHs(*mol);
      reacts.push_back(ROMOL_SPTR(mol));
      std::vector<MOL_SPTR_VECT> prods = rxn->runReactants(reacts);
      TEST_ASSERT(prods.size() == 1);
      TEST_ASSERT(prods[0].size() == 2);
      TEST_ASSERT(MolToSmiles(*prods[0][0], true) == "O=C=O");
      TEST_ASSERT(MolToSmiles(*prods[0][1], true) ==
                  "[H]N([H])C([H])([H])C([H])([H])[H]");
    }
    delete rxn;
  }

  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
}

void test62Github975() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing github975: Unnecessary warnings in rxn.validate()"
      << std::endl;

  {
    unsigned int nWarn, nError;
    std::string smi = "[N,O:1]>>[N+0,O+0:1]";

    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->validate(nWarn, nError, true));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);
    delete rxn;
    smi = "[N,O:1][C,N:2]>>[N+0,O+0:1][C+0,N+1:2]";
    rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->validate(nWarn, nError, true));
    TEST_ASSERT(nWarn == 1);
    TEST_ASSERT(nError == 0);
    delete rxn;
  }
  {
    unsigned int nWarn, nError;
    std::string smi = "[N,O:1]>>[NH,OH:1]";

    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->validate(nWarn, nError, true));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);
    delete rxn;
    smi = "[N,O:1][C,N:2]>>[NH,OH:1][CH2,NH:2]";
    rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->validate(nWarn, nError, true));
    TEST_ASSERT(nWarn == 1);
    TEST_ASSERT(nError == 0);
    delete rxn;
  }
  {
    unsigned int nWarn, nError;
    std::string smi = "[N,O:1]>>[14N,14O:1]";

    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->validate(nWarn, nError, true));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);
    delete rxn;
    smi = "[N,O:1][C,N:2]>>[14N,14O:1][12C,14N:2]";
    rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->validate(nWarn, nError, true));
    TEST_ASSERT(nWarn == 1);
    TEST_ASSERT(nError == 0);
    delete rxn;
  }

  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
}

void test63CopyConstructor() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing copy constructor" << std::endl;

  {
    std::string smi =
        "[C:1][C:2][C:3]([O:4])([N:5])[Cl:6]>>[Cl:6].[C:2][C:3]([O:4])([N:5])["
        "C:1]";

    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    std::string smi1 = ChemicalReactionToRxnSmiles(*rxn);
    TEST_ASSERT(rxn);
    auto *rxn_new = new ChemicalReaction(*rxn);
    removeMappingNumbersFromReactions(*rxn_new);
    std::string smi2 = ChemicalReactionToRxnSmiles(*rxn);
    std::string new_smi = ChemicalReactionToRxnSmiles(*rxn_new);
    std::cerr << "smi1 " << smi1 << std::endl;
    std::cerr << "smi2 " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);
    TEST_ASSERT(smi2 != new_smi);
    TEST_ASSERT(new_smi == "CCC(N)(O)Cl>>CC(C)(N)O.Cl");

    delete rxn;
    delete rxn_new;
  }

  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
}

void test64Github1266() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing github 1266: Reactions don't modify isotope "
                          "unless chemical element is specified for the product"
                       << std::endl;

  {
    std::string smi = "[13:1]>>[14:1]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();
    {
      MOL_SPTR_VECT reacts;
      smi = "C[13CH3]";
      RWMol *mol = SmilesToMol(smi);
      TEST_ASSERT(mol);
      reacts.push_back(ROMOL_SPTR(mol));
      std::vector<MOL_SPTR_VECT> prods = rxn->runReactants(reacts);
      TEST_ASSERT(prods.size() == 1);
      TEST_ASSERT(prods[0].size() == 1);
      TEST_ASSERT(prods[0][0]->getAtomWithIdx(0)->getIsotope() == 14);
    }

    delete rxn;
  }

  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
}

void test65SanitizeUnmappedHs() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Tests sanitize reaction (unmapped Hs) " << std::endl;

  const std::string unmappedHs =
      "$RXN\n"
      "\n"
      "  Marvin       031701170941\n"
      "\n"
      "  1  1\n"
      "$MOL\n"
      "\n"
      "  Mrv1583 03171709412D          \n"
      "\n"
      " 16 16  0  0  0  0            999 V2000\n"
      "   -2.5620    0.5265    0.0000 C   0  0  0  0  0  0  0  0  0  1  0  0\n"
      "   -2.4235   -0.2868    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0\n"
      "   -3.0587   -0.8133    0.0000 C   0  0  0  0  0  0  0  0  0  3  0  0\n"
      "   -3.8322   -0.5265    0.0000 C   0  0  0  0  0  0  0  0  0  4  0  0\n"
      "   -3.3355    0.8133    0.0000 C   0  0  0  0  0  0  0  0  0  6  0  0\n"
      "   -1.7581    0.7120    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -3.0942    1.6022    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -4.5124    0.9090    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -4.3663   -1.1553    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -3.3000   -1.6022    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -3.8242    1.4780    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -2.2307    1.2821    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -3.9706    0.2868    0.0000 C   0  0  0  0  0  0  0  0  0  5  0  0\n"
      "   -2.5700   -1.4780    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -1.6500   -0.5736    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "   -4.6278   -0.3082    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
      "  1  2  1  0  0  0  0\n"
      "  5  1  1  0  0  0  0\n"
      "  1  6  1  0  0  0  0\n"
      "  1 12  1  0  0  0  0\n"
      "  2  3  1  0  0  0  0\n"
      "  2 15  1  0  0  0  0\n"
      "  3  4  1  0  0  0  0\n"
      "  3 10  1  0  0  0  0\n"
      "  3 14  1  0  0  0  0\n"
      "  4  9  1  0  0  0  0\n"
      "  4 13  1  0  0  0  0\n"
      "  4 16  1  0  0  0  0\n"
      "  5  7  1  0  0  0  0\n"
      "  5 11  1  0  0  0  0\n"
      " 13  5  1  0  0  0  0\n"
      " 13  8  1  0  0  0  0\n"
      "M  END\n"
      "$MOL\n"
      "\n"
      "  Mrv1583 03171709412D          \n"
      "\n"
      "  6  6  0  0  0  0            999 V2000\n"
      "    3.8966    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  1  0  0\n"
      "    4.6111    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0\n"
      "    4.6111   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  3  0  0\n"
      "    3.8966   -0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  4  0  0\n"
      "    3.1821   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  5  0  0\n"
      "    3.1821    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  6  0  0\n"
      "  1  2  1  0  0  0  0\n"
      "  6  1  1  0  0  0  0\n"
      "  2  3  1  0  0  0  0\n"
      "  3  4  1  0  0  0  0\n"
      "  4  5  1  0  0  0  0\n"
      "  5  6  1  0  0  0  0\n"
      "M  END";

  ChemicalReaction *rxn = RxnBlockToChemicalReaction(unmappedHs);
  TEST_ASSERT(rxn);
  rxn->initReactantMatchers();
  MOL_SPTR_VECT reacts1, hreacts1, reacts2, hreacts2;
  std::vector<MOL_SPTR_VECT> prods;

  reacts1.push_back(ROMOL_SPTR(SmilesToMol("C1CCCCC1")));
  hreacts1.push_back(ROMOL_SPTR(MolOps::addHs(*reacts1[0].get())));

  reacts2.push_back(ROMOL_SPTR(SmilesToMol("C1CCCCC1Cl")));
  hreacts2.push_back(ROMOL_SPTR(MolOps::addHs(*reacts2[0].get())));

  // test with and without AddHs
  prods = rxn->runReactants(reacts1);
  TEST_ASSERT(prods.size() == 0);
  prods = rxn->runReactants(hreacts1);
  TEST_ASSERT(prods.size() == 768);

  prods = rxn->runReactants(reacts2);
  TEST_ASSERT(prods.size() == 0);

  prods = rxn->runReactants(hreacts2);
  TEST_ASSERT(prods.size() == 128);

  // Test after sanitization (way fewer matches than with AddHs..)
  RxnOps::sanitizeRxn(*rxn);
  prods = rxn->runReactants(reacts1);
  TEST_ASSERT(prods.size() == 12);

  prods = rxn->runReactants(reacts2);
  TEST_ASSERT(prods.size() == 4);

  delete rxn;

  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
}

void test66SanitizeMappedHs() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Tests sanitize reaction (mapped hs in react but not prod) "
      << std::endl;

  // H's are mapped in reactant but do not exist in product,
  //  they can be merged
  const std::string unmappedHs =
      "$RXN\n"
      "\n"
      "  Marvin       031701170941\n"
      "\n"
      "  1  1\n"
      "$MOL\n"
      "\n"
      "  Mrv1583 03171709412D          \n"
      "\n"
      " 16 16  0  0  0  0            999 V2000\n"
      "   -2.5620    0.5265    0.0000 C   0  0  0  0  0  0  0  0  0  1  0  0\n"
      "   -2.4235   -0.2868    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0\n"
      "   -3.0587   -0.8133    0.0000 C   0  0  0  0  0  0  0  0  0  3  0  0\n"
      "   -3.8322   -0.5265    0.0000 C   0  0  0  0  0  0  0  0  0  4  0  0\n"
      "   -3.3355    0.8133    0.0000 C   0  0  0  0  0  0  0  0  0  6  0  0\n"
      "   -1.7581    0.7120    0.0000 H   0  0  0  0  0  0  0  0  0  7  0  0\n"
      "   -3.0942    1.6022    0.0000 H   0  0  0  0  0  0  0  0  0  8  0  0\n"
      "   -4.5124    0.9090    0.0000 H   0  0  0  0  0  0  0  0  0  9  0  0\n"
      "   -4.3663   -1.1553    0.0000 H   0  0  0  0  0  0  0  0  0 10  0  0\n"
      "   -3.3000   -1.6022    0.0000 H   0  0  0  0  0  0  0  0  0 11  0  0\n"
      "   -3.8242    1.4780    0.0000 H   0  0  0  0  0  0  0  0  0 12  0  0\n"
      "   -2.2307    1.2821    0.0000 H   0  0  0  0  0  0  0  0  0 13  0  0\n"
      "   -3.9706    0.2868    0.0000 C   0  0  0  0  0  0  0  0  0  5  0  0\n"
      "   -2.5700   -1.4780    0.0000 H   0  0  0  0  0  0  0  0  0 14  0  0\n"
      "   -1.6500   -0.5736    0.0000 H   0  0  0  0  0  0  0  0  0 15  0  0\n"
      "   -4.6278   -0.3082    0.0000 H   0  0  0  0  0  0  0  0  0 16  0  0\n"
      "  1  2  1  0  0  0  0\n"
      "  5  1  1  0  0  0  0\n"
      "  1  6  1  0  0  0  0\n"
      "  1 12  1  0  0  0  0\n"
      "  2  3  1  0  0  0  0\n"
      "  2 15  1  0  0  0  0\n"
      "  3  4  1  0  0  0  0\n"
      "  3 10  1  0  0  0  0\n"
      "  3 14  1  0  0  0  0\n"
      "  4  9  1  0  0  0  0\n"
      "  4 13  1  0  0  0  0\n"
      "  4 16  1  0  0  0  0\n"
      "  5  7  1  0  0  0  0\n"
      "  5 11  1  0  0  0  0\n"
      " 13  5  1  0  0  0  0\n"
      " 13  8  1  0  0  0  0\n"
      "M  END\n"
      "$MOL\n"
      "\n"
      "  Mrv1583 03171709412D          \n"
      "\n"
      "  6  6  0  0  0  0            999 V2000\n"
      "    3.8966    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  1  0  0\n"
      "    4.6111    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0\n"
      "    4.6111   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  3  0  0\n"
      "    3.8966   -0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  4  0  0\n"
      "    3.1821   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  5  0  0\n"
      "    3.1821    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  6  0  0\n"
      "  1  2  1  0  0  0  0\n"
      "  6  1  1  0  0  0  0\n"
      "  2  3  1  0  0  0  0\n"
      "  3  4  1  0  0  0  0\n"
      "  4  5  1  0  0  0  0\n"
      "  5  6  1  0  0  0  0\n"
      "M  END";

  ChemicalReaction *rxn = RxnBlockToChemicalReaction(unmappedHs);
  TEST_ASSERT(rxn);
  rxn->initReactantMatchers();
  MOL_SPTR_VECT reacts1, hreacts1, reacts2, hreacts2;
  std::vector<MOL_SPTR_VECT> prods;

  reacts1.push_back(ROMOL_SPTR(SmilesToMol("C1CCCCC1")));
  hreacts1.push_back(ROMOL_SPTR(MolOps::addHs(*reacts1[0].get())));

  reacts2.push_back(ROMOL_SPTR(SmilesToMol("C1CCCCC1Cl")));
  hreacts2.push_back(ROMOL_SPTR(MolOps::addHs(*reacts2[0].get())));

  // test with and without AddHs
  prods = rxn->runReactants(reacts1);
  TEST_ASSERT(prods.size() == 0);
  prods = rxn->runReactants(hreacts1);
  TEST_ASSERT(prods.size() == 768);

  prods = rxn->runReactants(reacts2);
  TEST_ASSERT(prods.size() == 0);

  prods = rxn->runReactants(hreacts2);
  TEST_ASSERT(prods.size() == 128);

  // Test after sanitization (way fewer matches than with AddHs..)
  RxnOps::sanitizeRxn(*rxn);
  prods = rxn->runReactants(reacts1);
  TEST_ASSERT(prods.size() == 12);

  prods = rxn->runReactants(reacts2);
  TEST_ASSERT(prods.size() == 4);

  delete rxn;

  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
}

void test67SanitizeMappedHsInReactantAndProd() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Tests sanitize reaction (mapped hs in react and prod) " << std::endl;

  // H's are mapped in reactant and in prod
  //  they can be merged
  const std::string unmappedHs =
      "$RXN\n"
      "\n"
      "  Marvin       031701171002\n"
      "\n"
      "  1  1\n"
      "$MOL\n"
      "\n"
      "  Mrv1583 03171710022D          \n"
      "\n"
      " 16 16  0  0  0  0            999 V2000\n"
      "   -3.1881    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  1  0  0\n"
      "   -2.4736    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0\n"
      "   -2.4736   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  3  0  0\n"
      "   -3.1881   -0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  4  0  0\n"
      "   -3.9025    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  6  0  0\n"
      "   -2.8178    1.5622    0.0000 H   0  0  0  0  0  0  0  0  0  7  0  0\n"
      "   -4.3559    1.1018    0.0000 H   0  0  0  0  0  0  0  0  0  8  0  0\n"
      "   -4.6170   -0.8250    0.0000 H   0  0  0  0  0  0  0  0  0  9  0  0\n"
      "   -3.5583   -1.5622    0.0000 H   0  0  0  0  0  0  0  0  0 10  0  0\n"
      "   -2.0203   -1.1018    0.0000 H   0  0  0  0  0  0  0  0  0 11  0  0\n"
      "   -4.7262    0.4605    0.0000 H   0  0  0  0  0  0  0  0  0 12  0  0\n"
      "   -3.5583    1.5622    0.0000 H   0  0  0  0  0  0  0  0  0 13  0  0\n"
      "   -3.9025   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  5  0  0\n"
      "   -1.6500   -0.4605    0.0000 H   0  0  0  0  0  0  0  0  0 14  0  0\n"
      "   -1.7591    0.8250    0.0000 H   0  0  0  0  0  0  0  0  0 15  0  0\n"
      "   -2.8178   -1.5622    0.0000 H   0  0  0  0  0  0  0  0  0 16  0  0\n"
      "  1  2  1  0  0  0  0\n"
      "  5  1  1  0  0  0  0\n"
      "  1  6  1  0  0  0  0\n"
      "  1 12  1  0  0  0  0\n"
      "  2  3  1  0  0  0  0\n"
      "  2 15  1  0  0  0  0\n"
      "  3  4  1  0  0  0  0\n"
      "  3 10  1  0  0  0  0\n"
      "  3 14  1  0  0  0  0\n"
      "  4  9  1  0  0  0  0\n"
      "  4 13  1  0  0  0  0\n"
      "  4 16  1  0  0  0  0\n"
      "  5  7  1  0  0  0  0\n"
      "  5 11  1  0  0  0  0\n"
      " 13  5  1  0  0  0  0\n"
      " 13  8  1  0  0  0  0\n"
      "M  END\n"
      "$MOL\n"
      "\n"
      "  Mrv1583 03171710022D          \n"
      "\n"
      " 17 17  0  0  0  0            999 V2000\n"
      "    4.1309    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  1  0  0\n"
      "    4.8454    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0\n"
      "    4.8454   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  3  0  0\n"
      "    4.1309   -0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  4  0  0\n"
      "    3.4165   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  5  0  0\n"
      "    3.4165    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  6  0  0\n"
      "    4.5012    1.5622    0.0000 H   0  0  0  0  0  0  0  0  0 13  0  0\n"
      "    3.7607    1.5622    0.0000 H   0  0  0  0  0  0  0  0  0  7  0  0\n"
      "    5.6690    0.4605    0.0000 H   0  0  0  0  0  0  0  0  0 15  0  0\n"
      "    5.2987    1.1018    0.0000 H   0  0  0  0  0  0  0  0  0 17  0  0\n"
      "    5.2987   -1.1018    0.0000 H   0  0  0  0  0  0  0  0  0 14  0  0\n"
      "    3.7607   -1.5622    0.0000 H   0  0  0  0  0  0  0  0  0 16  0  0\n"
      "    2.9631   -1.1018    0.0000 H   0  0  0  0  0  0  0  0  0  9  0  0\n"
      "    2.5929   -0.4605    0.0000 H   0  0  0  0  0  0  0  0  0 18  0  0\n"
      "    2.7020    0.8250    0.0000 H   0  0  0  0  0  0  0  0  0  8  0  0\n"
      "    4.5012   -1.5622    0.0000 H   0  0  0  0  0  0  0  0  0 10  0  0\n"
      "    5.6690   -0.4605    0.0000 H   0  0  0  0  0  0  0  0  0 11  0  0\n"
      "  1  2  1  0  0  0  0\n"
      "  6  1  1  0  0  0  0\n"
      "  1  7  1  0  0  0  0\n"
      "  1  8  1  0  0  0  0\n"
      "  2  9  1  0  0  0  0\n"
      "  2 10  1  0  0  0  0\n"
      "  2  3  1  0  0  0  0\n"
      "  3  4  1  0  0  0  0\n"
      "  3 11  1  0  0  0  0\n"
      "  3 17  1  0  0  0  0\n"
      "  4  5  1  0  0  0  0\n"
      "  4 12  1  0  0  0  0\n"
      "  4 16  1  0  0  0  0\n"
      "  5  6  1  0  0  0  0\n"
      "  5 14  1  0  0  0  0\n"
      "  5 13  1  0  0  0  0\n"
      "  6 15  1  0  0  0  0\n"
      "M  END\n";

  ChemicalReaction *rxn = RxnBlockToChemicalReaction(unmappedHs);
  TEST_ASSERT(rxn);
  rxn->initReactantMatchers();
  MOL_SPTR_VECT reacts1, hreacts1, reacts2, hreacts2;
  std::vector<MOL_SPTR_VECT> prods;

  reacts1.push_back(ROMOL_SPTR(SmilesToMol("C1CCCCC1")));
  hreacts1.push_back(ROMOL_SPTR(MolOps::addHs(*reacts1[0].get())));

  reacts2.push_back(ROMOL_SPTR(SmilesToMol("C1CCCCC1Cl")));
  hreacts2.push_back(ROMOL_SPTR(MolOps::addHs(*reacts2[0].get())));

  // test with and without AddHs
  prods = rxn->runReactants(reacts1);
  TEST_ASSERT(prods.size() == 0);
  prods = rxn->runReactants(hreacts1);
  TEST_ASSERT(prods.size() == 768);

  prods = rxn->runReactants(reacts2);
  TEST_ASSERT(prods.size() == 0);

  prods = rxn->runReactants(hreacts2);
  TEST_ASSERT(prods.size() == 128);

  // Test after sanitization (way fewer matches than with AddHs..)
  RxnOps::sanitizeRxn(*rxn);
  prods = rxn->runReactants(reacts1);
  TEST_ASSERT(prods.size() == 12);

  prods = rxn->runReactants(reacts2);
  TEST_ASSERT(prods.size() == 4);

  delete rxn;

  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
}
void test68MappedHToHeavy() {
  const std::string rxnblock =
      "$RXN\n"
      "\n"
      "  Marvin       031701171005\n"
      "\n"
      "  1  1\n"
      "$MOL\n"
      "\n"
      "  Mrv1583 03171710052D          \n"
      "\n"
      "  3  2  0  0  0  0            999 V2000\n"
      "   -1.2721   -0.0116    0.0000 C   0  0  0  0  0  0  0  0  0  1  0  0\n"
      "   -1.9866   -0.4241    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0\n"
      "   -1.2721    0.8134    0.0000 H   0  0  0  0  0  0  0  0  0  3  0  0\n"
      "  2  1  1  0  0  0  0\n"
      "  1  3  1  0  0  0  0\n"
      "M  END\n"
      "$MOL\n"
      "\n"
      "  Mrv1583 03171710052D          \n"
      "\n"
      "  3  2  0  0  0  0            999 V2000\n"
      "    2.3886   -0.0563    0.0000 C   0  0  0  0  0  0  0  0  0  1  0  0\n"
      "    1.6741   -0.4688    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0\n"
      "    2.3886    0.7688    0.0000 C   0  0  0  0  0  0  0  0  0  3  0  0\n"
      "  2  1  1  0  0  0  0\n"
      "  1  3  1  0  0  0  0\n"
      "M  END\n";

  ChemicalReaction *rxn = RxnBlockToChemicalReaction(rxnblock);
  TEST_ASSERT(rxn);
  rxn->initReactantMatchers();
  MOL_SPTR_VECT reacts1, hreacts1, reacts2, hreacts2;
  std::vector<MOL_SPTR_VECT> prods;

  reacts1.push_back(ROMOL_SPTR(SmilesToMol("CC")));
  hreacts1.push_back(ROMOL_SPTR(MolOps::addHs(*reacts1[0].get())));

  // test with and without AddHs
  prods = rxn->runReactants(reacts1);
  TEST_ASSERT(prods.size() == 0);

  prods = rxn->runReactants(hreacts1);
  TEST_ASSERT(prods.size() == 6);

  std::stringstream sstrm;
  rdWarningLog->SetTee(sstrm);
  RxnOps::sanitizeRxn(*rxn);
  std::string s = sstrm.str();
  std::cerr << s << std::endl;
  TEST_ASSERT(s.find("Reaction has explicit hydrogens, reactants will need "
                     "explicit hydrogens (addHs)") != std::string::npos);
  rdWarningLog->ClearTee();

  delete rxn;
}

void test69Github1387() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Tests github1387: Bond from reactant not added to product"
      << std::endl;

  {
    const std::string smarts =
        "([*:1]-[O:2][CH3:3].[*:4]-[O:5][CH3:6])>>([*:1]-[O:2][C:3]C.[*:4]-[O:"
        "5][C:6]C)";

    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smarts);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();
    {  // this always worked
      MOL_SPTR_VECT reacts;
      std::vector<MOL_SPTR_VECT> prods;
      reacts.push_back(ROMOL_SPTR(SmilesToMol("COCCCOC")));
      prods = rxn->runReactants(reacts);
      TEST_ASSERT(prods.size() == 2);
      TEST_ASSERT(prods[0].size() == 1);
      std::vector<int> mapping;
      unsigned int nPieces = MolOps::getMolFrags(*prods[0][0], mapping);
      TEST_ASSERT(nPieces = 2);
      TEST_ASSERT(prods[0][0]->getBondBetweenAtoms(0, 4) == nullptr);
    }
    {  // the bug:
      MOL_SPTR_VECT reacts;
      std::vector<MOL_SPTR_VECT> prods;
      reacts.push_back(ROMOL_SPTR(SmilesToMol("COCCOC")));
      prods = rxn->runReactants(reacts);
      TEST_ASSERT(prods.size() == 2);
      TEST_ASSERT(prods[0].size() == 1);
      std::vector<int> mapping;
      unsigned int nPieces = MolOps::getMolFrags(*prods[0][0], mapping);
      TEST_ASSERT(nPieces = 2);
      TEST_ASSERT(prods[0][0]->getBondBetweenAtoms(0, 4) != nullptr);
    }
    {  // the bug, plus a ring closure:
      MOL_SPTR_VECT reacts;
      std::vector<MOL_SPTR_VECT> prods;
      reacts.push_back(ROMOL_SPTR(SmilesToMol("COC1CCC1OC")));
      prods = rxn->runReactants(reacts);
      TEST_ASSERT(prods.size() == 2);
      TEST_ASSERT(prods[0].size() == 1);
      std::vector<int> mapping;
      unsigned int nPieces = MolOps::getMolFrags(*prods[0][0], mapping);
      TEST_ASSERT(nPieces = 2);
      TEST_ASSERT(prods[0][0]->getBondBetweenAtoms(0, 4) != nullptr);
    }
    delete rxn;
  }
  {  // reorder the same reaction
    const std::string smarts =
        "([O:2](-[*:1])[CH3:3].[CH3:6][O:5]-[*:4])>>(C[C:3][O:2]-[*:1].[*:4]"
        "-[O:5][C:6]C)";

    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smarts);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();
    {  // this always worked
      MOL_SPTR_VECT reacts;
      std::vector<MOL_SPTR_VECT> prods;
      reacts.push_back(ROMOL_SPTR(SmilesToMol("COCCCOC")));
      prods = rxn->runReactants(reacts);
      TEST_ASSERT(prods.size() == 2);
      TEST_ASSERT(prods[0].size() == 1);
      std::vector<int> mapping;
      unsigned int nPieces = MolOps::getMolFrags(*prods[0][0], mapping);
      TEST_ASSERT(nPieces = 2);
      TEST_ASSERT(prods[0][0]->getBondBetweenAtoms(3, 4) == nullptr);
    }
    {  // the bug:
      MOL_SPTR_VECT reacts;
      std::vector<MOL_SPTR_VECT> prods;
      reacts.push_back(ROMOL_SPTR(SmilesToMol("COCCOC")));
      prods = rxn->runReactants(reacts);
      TEST_ASSERT(prods.size() == 2);
      TEST_ASSERT(prods[0].size() == 1);
      std::vector<int> mapping;
      unsigned int nPieces = MolOps::getMolFrags(*prods[0][0], mapping);
      TEST_ASSERT(nPieces = 2);
      TEST_ASSERT(prods[0][0]->getBondBetweenAtoms(3, 4) != nullptr);
    }
    {  // the bug, plus a ring closure:
      MOL_SPTR_VECT reacts;
      std::vector<MOL_SPTR_VECT> prods;
      reacts.push_back(ROMOL_SPTR(SmilesToMol("COC1CCC1OC")));
      prods = rxn->runReactants(reacts);
      TEST_ASSERT(prods.size() == 2);
      TEST_ASSERT(prods[0].size() == 1);
      std::vector<int> mapping;
      unsigned int nPieces = MolOps::getMolFrags(*prods[0][0], mapping);
      TEST_ASSERT(nPieces = 2);
      TEST_ASSERT(prods[0][0]->getBondBetweenAtoms(3, 4) != nullptr);
    }

    delete rxn;
  }
  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
}

void test70Github1544() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Tests github1544: ChemicalReaction code not calling "
                          "setNoImplicit() when H counts are set"
                       << std::endl;
  {
    const std::string smarts = "[CH3:1]>>[CH2:1]";

    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smarts);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();
    {
      MOL_SPTR_VECT reacts;
      std::vector<MOL_SPTR_VECT> prods;
      reacts.push_back(ROMOL_SPTR(SmilesToMol("Cc1ccccc1")));
      prods = rxn->runReactants(reacts);
      TEST_ASSERT(prods.size() == 1);
      TEST_ASSERT(prods[0].size() == 1);
      TEST_ASSERT(prods[0][0]->getAtomWithIdx(0)->getNoImplicit());
      TEST_ASSERT(MolToSmiles(*prods[0][0]) == "[CH2]c1ccccc1");
    }
    delete rxn;
  }
  {  // test that the change also works when the degree on the atom changes
    const std::string smarts = "[CH3:1][*:2]>>[CH2:1].[*:2]";

    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smarts);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();
    {
      MOL_SPTR_VECT reacts;
      std::vector<MOL_SPTR_VECT> prods;
      reacts.push_back(ROMOL_SPTR(SmilesToMol("Cc1ccccc1")));
      prods = rxn->runReactants(reacts);
      TEST_ASSERT(prods.size() == 1);
      TEST_ASSERT(prods[0].size() == 2);
      TEST_ASSERT(prods[0][0]->getAtomWithIdx(0)->getNoImplicit());
      TEST_ASSERT(MolToSmiles(*prods[0][0]) == "[CH2]");
      TEST_ASSERT(MolToSmiles(*prods[0][1]) == "c1ccccc1");
    }
    delete rxn;
  }
  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
}

void testSanitizeException() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing sanitizerxn exception" << std::endl;

  // This is a marvinjs conversion of the following smirks:
  // {indole}\t[*;Br,I;$(*c1ccccc1)]-[c:1]:[c:2]-[NH2:3].[CH1:5]#[C;$(C-[#6]):4]>>[c:1]1:[c:2]-[N:3]-[C:4]=[C:5]-1\tc1cc(I)c(N)cc1\tCC#C
  // From Hartenfeller et. al., J. Chem. Inf. Model., 2012, 52 (5), p 1167-1178
  std::string rxnB =
      "$RXN\n\n  Marvin       102501170854\n\n  2  1\n$MOL\n\n  Mrv1583 "
      "10251708542D          \n\n  4  3  0  0  0  0            999 V2000\n   "
      "-4.2429   -0.3572    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0\n   "
      "-5.0679   -0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  1  0  0\n   "
      "-5.4804    0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0\n   "
      "-6.3054    0.3572    0.0000 N   0  0  0  3  0  0  0  0  0  3  0  0\n  1 "
      " 2  1  0  0  0  0\n  2  3  4  0  0  0  0\n  3  4  1  0  0  0  0\nM  MRV "
      "SMA   1 [*;A]\nM  MRV SMA   4 [#7H2;A:3]\nM  END\n$MOL\n\n  Mrv1583 "
      "10251708542D          \n\n  2  1  0  0  0  0            999 V2000\n   "
      "-1.6500    0.0000    0.0000 C   0  0  0  2  0  0  0  0  0  5  0  0\n   "
      "-2.4750    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  4  0  0\n  1 "
      " 2  3  0  0  0  0\nM  MRV SMA   1 [#6H1:5]\nM  END\n$MOL\n\n  Mrv1583 "
      "10251708542D          \n\n  5  5  0  0  0  0            999 V2000\n    "
      "3.3782    0.6348    0.0000 C   0  0  0  0  0  0  0  0  0  1  0  0\n    "
      "4.0456    0.1498    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0\n    "
      "3.7907   -0.6348    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0\n    "
      "2.9657   -0.6348    0.0000 C   0  0  0  0  0  0  0  0  0  4  0  0\n    "
      "2.7107    0.1498    0.0000 C   0  0  0  0  0  0  0  0  0  5  0  0\n  1  "
      "2  4  0  0  0  0\n  1  5  1  0  0  0  0\n  2  3  1  0  0  0  0\n  3  4  "
      "1  0  0  0  0\n  4  5  2  0  0  0  0\nM  END\n";

  ChemicalReaction *rxn = RxnBlockToChemicalReaction(rxnB);
  RxnOps::sanitizeRxn(*rxn);
  delete rxn;
}

void testReactionProperties() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Tests reaction properties" << std::endl;
  {
    const std::string smarts = "[CH3:1]>>[CH2:1]";

    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smarts);
    TEST_ASSERT(rxn);
    TEST_ASSERT(!rxn->hasProp("fooprop"));
    rxn->setProp("fooprop", 3);
    TEST_ASSERT(rxn->hasProp("fooprop"));
    TEST_ASSERT(rxn->getProp<int>("fooprop") == 3);

    {
      std::string pkl;
      ReactionPickler::pickleReaction(rxn, pkl);
      auto *lrxn = new ChemicalReaction();
      ReactionPickler::reactionFromPickle(pkl, lrxn);
      TEST_ASSERT(!lrxn->hasProp("fooprop"));
      delete lrxn;
    }

    {
      std::string pkl;
      ReactionPickler::pickleReaction(rxn, pkl, PicklerOps::AllProps);
      auto *lrxn = new ChemicalReaction();
      ReactionPickler::reactionFromPickle(pkl, lrxn);
      TEST_ASSERT(lrxn->hasProp("fooprop"));
      TEST_ASSERT(lrxn->getProp<int>("fooprop") == 3);
      delete lrxn;
    }
    delete rxn;
  }
  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
}

void testGithub1950() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github #1950: Query features in products of "
                          "rxn files not properly handled"
                       << std::endl;

  std::string rdbase = getenv("RDBASE");

  {  // the original report
    auto fName =
        rdbase + "/Code/GraphMol/ChemReactions/testData/github1950_1.rxn";

    auto rxn = RxnFileToChemicalReaction(fName);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);

    MOL_SPTR_VECT reacts;
    std::string smi = "c1ccccc1Cl";
    auto mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));

    smi = "CCO";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));

    rxn->initReactantMatchers();
    auto prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 2);  // symmetric, so two products
    TEST_ASSERT(prods[0].size() == 1);

    smi = MolToSmiles(*prods[0][0]);
    TEST_ASSERT(smi == "CCOc1ccccc1");
    delete rxn;
  }
  {  // a modification using MRV SMA
    auto fName =
        rdbase + "/Code/GraphMol/ChemReactions/testData/github1950_2.rxn";
    std::cerr << "-------------" << std::endl;
    auto rxn = RxnFileToChemicalReaction(fName);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);

    MOL_SPTR_VECT reacts;
    std::string smi = "c1ccccc1Cl";
    auto mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));

    smi = "CCO";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));

    rxn->initReactantMatchers();
    auto prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 2);  // symmetric, so two products
    TEST_ASSERT(prods[0].size() == 1);

    smi = MolToSmiles(*prods[0][0]);
    std::cerr << smi << std::endl;
    TEST_ASSERT(smi == "CCOc1ccccc1");
    delete rxn;
  }

  {  // make sure we didn't break anything,
     // this is roughly the same reaction as SMARTS
    std::string smarts =
        "Cl[#6:2]=,:[#6,#7:3].[#8:5]-[#6:4]>>[#6:4]-[#8:5]-[#6:2]~[*:3]";
    auto rxn = RxnSmartsToChemicalReaction(smarts);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);

    MOL_SPTR_VECT reacts;
    std::string smi = "c1ccccc1Cl";
    auto mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));

    smi = "CCO";
    mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));

    rxn->initReactantMatchers();
    auto prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 2);  // symmetric, so two products
    TEST_ASSERT(prods[0].size() == 1);

    smi = MolToSmiles(*prods[0][0]);
    TEST_ASSERT(smi == "CCOc1ccccc1");
    delete rxn;
  }
  BOOST_LOG(rdInfoLog) << "Done" << std::endl;
}

void testGithub1869() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing github #1869 and #1955: parens around "
                          "reactant/product fragments."
                       << std::endl;
  {
    std::string smi;
    // nonsense test case
    smi = "([C:1].[Cl:2]).[N:3]-[C:4]>>[Cl:2]-[C:1].([N:3].[C:4])";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 2);
    rxn->initReactantMatchers();
    unsigned int nWarn, nError;
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);

    smi = ChemicalReactionToRxnSmarts(*rxn);
    TEST_ASSERT(smi ==
                "([C:1].[Cl:2]).[N:3]-[C:4]>>[Cl:2]-[C:1].([N:3].[C:4])");

    delete rxn;
  }
  {  // #1869 example
    std::string smi;
    smi = "[R:1]~;!@[*:2]>>([*:1].[*:2])";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();
    unsigned int nWarn, nError;
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);

    smi = ChemicalReactionToRxnSmarts(*rxn);
    TEST_ASSERT(smi == "[R:1]!@[*:2]>>([*:1].[*:2])");

    delete rxn;
  }
  {  // #1955 example
    std::string smi;
    smi = "([C:1].[C:2]).([C:3][C:4])>>[#6:1]-[#6:2]-[#6:3]-[#6:4]";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi);
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();
    unsigned int nWarn, nError;
    TEST_ASSERT(rxn->validate(nWarn, nError, false));
    TEST_ASSERT(nWarn == 0);
    TEST_ASSERT(nError == 0);

    smi = ChemicalReactionToRxnSmarts(*rxn);
    TEST_ASSERT(smi == "([C:1].[C:2]).[C:3][C:4]>>[#6:1]-[#6:2]-[#6:3]-[#6:4]");

    delete rxn;
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGithub1269() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing github #1269: preserve reactant atom idx in products"
      << std::endl;
  {
    // nonsense test case
    std::string sma = "[C:1]>>[O:1]Cl";
    ChemicalReaction *rxn = RxnSmartsToChemicalReaction(sma);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();

    MOL_SPTR_VECT reacts;
    std::string smi = "NC";
    auto mol = SmilesToMol(smi);
    TEST_ASSERT(mol);
    reacts.push_back(ROMOL_SPTR(mol));

    auto prods = rxn->runReactants(reacts);
    TEST_ASSERT(prods.size() == 1);
    TEST_ASSERT(prods[0].size() == 1);
    TEST_ASSERT(prods[0][0]->getAtomWithIdx(0)->getAtomicNum() == 8);
    TEST_ASSERT(prods[0][0]->getAtomWithIdx(0)->hasProp(
        common_properties::reactantAtomIdx));
    TEST_ASSERT(prods[0][0]->getAtomWithIdx(0)->getProp<unsigned int>(
                    "react_atom_idx") == 1);
    TEST_ASSERT(!prods[0][0]->getAtomWithIdx(1)->hasProp(
        common_properties::reactantAtomIdx));
    TEST_ASSERT(prods[0][0]->getAtomWithIdx(2)->hasProp(
        common_properties::reactantAtomIdx));
    TEST_ASSERT(prods[0][0]->getAtomWithIdx(2)->getProp<unsigned int>(
                    "react_atom_idx") == 0);

    delete rxn;
  }
}

void testGithub1868() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github #1868: Atom index out of range error"
                       << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName;

  fName = rdbase + "/Code/GraphMol/ChemReactions/testData/v3k.AmideBond.rxn";

  for (int i = 0; i < 1000; ++i) {
    std::unique_ptr<ChemicalReaction> rxn(RxnFileToChemicalReaction(fName));
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);

    for (auto v : rxn->getReactants()) {
      MolToSmiles(*v.get());
    }
    for (auto v : rxn->getProducts()) {
      MolToSmiles(*v.get());
    }
  }
}

bool check_bond_stereo(const ROMOL_SPTR &mol, unsigned bond_idx,
                       int start_anchor_idx, int end_anchor_idx,
                       Bond::BondStereo stereo) {
  auto *bond = mol->getBondWithIdx(bond_idx);
  std::vector<int> stereo_atoms_ref{{start_anchor_idx, end_anchor_idx}};
  return Bond::BondType::DOUBLE == bond->getBondType() &&
         stereo == bond->getStereo() &&
         stereo_atoms_ref == bond->getStereoAtoms();
}

ROMOL_SPTR run_simple_reaction(const std::string &reaction,
                               const ROMOL_SPTR &reactant) {
  std::unique_ptr<ChemicalReaction> rxn{RxnSmartsToChemicalReaction(reaction)};
  TEST_ASSERT(rxn);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);
  rxn->initReactantMatchers();

  const auto prods = rxn->runReactant(reactant, 0);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);

  return prods[0][0];
};

void run_gist_reaction_tests(const std::vector<RWMOL_SPTR> &mols) {
  TEST_ASSERT(check_bond_stereo(mols[0], 1, 0, 3, Bond::BondStereo::STEREOE));
  TEST_ASSERT(check_bond_stereo(mols[1], 1, 0, 3, Bond::BondStereo::STEREOE));
  TEST_ASSERT(check_bond_stereo(mols[2], 2, 2, 4, Bond::BondStereo::STEREOZ));
  TEST_ASSERT(check_bond_stereo(mols[3], 1, 0, 3, Bond::BondStereo::STEREOE));

  // StereoAtoms atoms get set such that TRANS == E and CIS == Z

  {  // Example 1: reaction that does not affect double bond or neighboring
     // single bonds
    const std::string reaction(R"([Br:1]>>[F:1])");

    {
      // 1a
      const auto product = run_simple_reaction(reaction, mols[0]);
      TEST_ASSERT(
          check_bond_stereo(product, 2, 4, 1, Bond::BondStereo::STEREOTRANS));
    }

    {
      // 1b
      const auto product = run_simple_reaction(reaction, mols[1]);
      TEST_ASSERT(
          check_bond_stereo(product, 1, 3, 0, Bond::BondStereo::STEREOTRANS));
    }
  }
  {  // Example 2: reaction that does affect neighboring single bonds
    const std::string reaction(R"([C:2][Br:1]>>[C:2][F:1])");

    {
      // 2a
      const auto product = run_simple_reaction(reaction, mols[1]);
      TEST_ASSERT(
          check_bond_stereo(product, 1, 3, 1, Bond::BondStereo::STEREOTRANS));
    }

    {
      // 2b
      const auto product = run_simple_reaction(reaction, mols[2]);
      TEST_ASSERT(
          check_bond_stereo(product, 2, 1, 4, Bond::BondStereo::STEREOCIS));
    }
  }
  {  // Example 3: reaction that affects the double bond itself
    const std::string reaction(R"([C:1]=[N:2]>>[C:1]=[C:2])");

    auto product = run_simple_reaction(reaction, mols[3]);
    TEST_ASSERT(
        check_bond_stereo(product, 0, 2, 3, Bond::BondStereo::STEREOTRANS));
  }
}

void strip_bond_directions(RWMOL_SPTR &mol) {
  for (auto &bond : mol->bonds()) {
    if (bond->getBondDir() == Bond::BondDir::ENDDOWNRIGHT ||
        bond->getBondDir() == Bond::BondDir::ENDUPRIGHT) {
      bond->setBondDir(Bond::BondDir::NONE);
    }
  }
}

void testBondStereoGistExamples() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing Bond Stereo propagation across reactions in gist examples"
      << std::endl;

  // Examples from
  // https://gist.github.com/greglandrum/1c7d933537fd97f14a5c1ec3419f8d7a

  // All these mols have Z stereo, except mol3, which is Z
  auto mol1 = RWMOL_SPTR(SmilesToMol(R"(C/C=C/C[Br])"));
  auto mol2 = RWMOL_SPTR(SmilesToMol(R"(C/C=C/[Br])"));
  auto mol3 = RWMOL_SPTR(SmilesToMol(R"(C/C(Br)=C/C)"));
  auto mol4 = RWMOL_SPTR(SmilesToMol(R"(C/C=N/F)"));

  std::vector<RWMOL_SPTR> mols{mol1, mol2, mol3, mol4};

  run_gist_reaction_tests(mols);

  // Now, strip bond directions and rerun the same tests to make sure tests
  // yield the same results without bond directions in the reactants

  for (auto &mol : mols) {
    strip_bond_directions(mol);
  }

  run_gist_reaction_tests(mols);
}

void testStereoBondIsomerization() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Stereo Bond isomerization" << std::endl;

  auto mol = RWMOL_SPTR(SmilesToMol(R"(C/C(Br)=C(Cl)/F)"));
  strip_bond_directions(mol);
  TEST_ASSERT(check_bond_stereo(mol, 2, 2, 4, Bond::BondStereo::STEREOE));

  // StereoAtoms will get set to the end of the directed bonds in the products

  {  // Try a null reaction
    const std::string reaction(
        R"([C:1]/[C:2](-[Br:3])=[C:4](\[Cl:5])-[F:6]>>[C:1]/[C:2](-[Br:3])=[C:4](\[Cl:5])-[F:6])");

    const auto product = run_simple_reaction(reaction, mol);
    TEST_ASSERT(
        check_bond_stereo(product, 2, 0, 4, Bond::BondStereo::STEREOCIS));
  }

  {  // Flip the bond stereo (second bond direction is flipped in the product)
    const std::string reaction(
        R"([C:1]/[C:2](-[Br:3])=[C:4](\[Cl:5])-[F:6]>>[C:1]/[C:2](-[Br:3])=[C:4](/[Cl:5])-[F:6])");

    const auto product = run_simple_reaction(reaction, mol);
    TEST_ASSERT(
        check_bond_stereo(product, 2, 0, 4, Bond::BondStereo::STEREOTRANS));
  }

  {  // This should fail, but doesn't:
    // we attempt a null reaction on the opposite isomer of the reactant (second
    // bond direction is flipped), but it sill works because RDKit cannot match
    // SMARTS bond directions in the reactant, but still can force them onto the
    // products of the reaction.
    const std::string reaction(
        R"([C:1]/[C:2](-[Br:3])=[C:4](/[Cl:5])-[F:6]>>[C:1]/[C:2](-[Br:3])=[C:4](/[Cl:5])-[F:6])");

    const auto product = run_simple_reaction(reaction, mol);
    TEST_ASSERT(
        check_bond_stereo(product, 2, 0, 4, Bond::BondStereo::STEREOTRANS));
  }

  {  // Isomerize swapping the directed bond to the other atom
    // Note that the stereo is the same, but we changed the stereoatom it is
    // referred to!
    const std::string reaction(
        R"([C:1]/[C:2](-[Br:3])=[C:4](\[Cl:5])-[F:6]>>[C:1]/[C:2](-[Br:3])=[C:4](-[Cl:5])\[F:6])");

    const auto product = run_simple_reaction(reaction, mol);
    TEST_ASSERT(
        check_bond_stereo(product, 2, 0, 5, Bond::BondStereo::STEREOCIS));
  }

  // From here on, StereoAtoms set to reactant's StereoAtoms (= highest CIPs)

  {  // One-side reaction isomerization (no double bond in reaction)
    const std::string reaction(R"([C:1](\[Cl:2])-[F:3]>>[C:1](/[Cl:2])-[F:3])");

    const auto product = run_simple_reaction(reaction, mol);
    TEST_ASSERT(
        check_bond_stereo(product, 2, 5, 1, Bond::BondStereo::STEREOTRANS));
  }

  {  // One-side reaction isomerization II (double bond in the reaction)
    const std::string reaction(
        R"([C:1]=[C:2](\[Cl:3])-[F:4]>>[C:1]=[C:2](/[Cl:3])-[F:4])");

    const auto product = run_simple_reaction(reaction, mol);
    TEST_ASSERT(
        check_bond_stereo(product, 0, 5, 2, Bond::BondStereo::STEREOTRANS));
  }
}

void testOtherBondStereo() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Other Bond Stereo reaction cases"
                       << std::endl;

  auto mol = RWMOL_SPTR(SmilesToMol(R"(C/C=C/[Br])"));
  strip_bond_directions(mol);
  TEST_ASSERT(check_bond_stereo(mol, 1, 0, 3, Bond::BondStereo::STEREOE));

  {  // Reaction changes order of the stereo bond
    const std::string reaction(R"([C:1]=[C:2]>>[C:1]-[C:2])");
    std::unique_ptr<ChemicalReaction> rxn{
        RxnSmartsToChemicalReaction(reaction)};
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();

    const auto product_sets = rxn->runReactant(mol, 0);
    TEST_ASSERT(!product_sets.empty());
    for (const auto &products : product_sets) {
      TEST_ASSERT(products.size() == 1);
      for (const Bond *bond : products[0]->bonds()) {
        TEST_ASSERT(Bond::BondType::SINGLE == bond->getBondType());
        TEST_ASSERT(Bond::BondStereo::STEREONONE == bond->getStereo());
      }
    }
    // and if we don't strip direction?
    auto lmol = RWMOL_SPTR(SmilesToMol(R"(C/C=C/[Br])"));
    auto product_sets2 = rxn->runReactant(lmol, 0);
    TEST_ASSERT(!product_sets2.empty());
    auto osmi = MolToSmiles(*product_sets2[0][0]);
    TEST_ASSERT(osmi == "CCCBr");
  }
  {  // Reaction explicitly destroys stereochemistry
     // (no directed bonds enclosing the double bond)
    const std::string reaction(
        R"([C:1]/[C:2]=[C:3]/[Br:4]>>[C:1][C:2]=[C:3]-[Br:4])");
    std::unique_ptr<ChemicalReaction> rxn{
        RxnSmartsToChemicalReaction(reaction)};
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();

    const auto product_sets = rxn->runReactant(mol, 0);
    TEST_ASSERT(!product_sets.empty());
    for (const auto &products : product_sets) {
      TEST_ASSERT(products.size() == 1);
      for (const Bond *bond : products[0]->bonds()) {
        TEST_ASSERT(Bond::BondStereo::STEREONONE == bond->getStereo());

        // Make sure the temporary mark set in the reaction has been removed.
        TEST_ASSERT(!bond->hasProp("_UnknownStereoRxnBond"));

        if (bond->getIdx() == 1) {
          TEST_ASSERT(Bond::BondType::DOUBLE == bond->getBondType());
        } else {
          TEST_ASSERT(Bond::BondType::SINGLE == bond->getBondType());
        }
      }
    }
    // and if we don't strip direction?
    auto lmol = RWMOL_SPTR(SmilesToMol(R"(C/C=C/[Br])"));
    auto product_sets2 = rxn->runReactant(lmol, 0);
    TEST_ASSERT(!product_sets.empty());
    auto osmi = MolToSmiles(*product_sets2[0][0]);
    TEST_ASSERT(osmi == "CC=CBr");
  }
  {  // Reaction with 2 product sets
    const std::string reaction(R"([C:1]=[C:2]>>[Si:1]=[C:2])");
    std::unique_ptr<ChemicalReaction> rxn{
        RxnSmartsToChemicalReaction(reaction)};
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    rxn->initReactantMatchers();

    const auto product_sets = rxn->runReactant(mol, 0);
    TEST_ASSERT(product_sets.size() == 2);
    for (const auto &products : product_sets) {
      TEST_ASSERT(products.size() == 1);
      TEST_ASSERT(check_bond_stereo(products[0], 0, 2, 3,
                                    Bond::BondStereo::STEREOTRANS));
    }
    // and if we don't strip direction?
    auto lmol = RWMOL_SPTR(SmilesToMol(R"(C/C=C/[Br])"));
    auto product_sets2 = rxn->runReactant(lmol, 0);
    TEST_ASSERT(!product_sets2.empty());
    auto osmi = MolToSmiles(*product_sets2[0][0]);
    TEST_ASSERT(osmi == "C/[SiH]=C/Br");
    osmi = MolToSmiles(*product_sets2[1][0]);
    TEST_ASSERT(osmi == "C/C=[SiH]/Br");
  }
  {  // Reactant stereo propagated by (stereoatom, anti-stereoatom) pair
    auto mol2 = RWMOL_SPTR(SmilesToMol(R"(Cl/C(C)=C(/Br)F)"));
    strip_bond_directions(mol2);
    TEST_ASSERT(check_bond_stereo(mol2, 2, 0, 4, Bond::BondStereo::STEREOE));

    const std::string reaction(R"([C:1]=[C:2][Br:3]>>[C:1]=[C:2].[Br:3])");
    std::unique_ptr<ChemicalReaction> rxn{
        RxnSmartsToChemicalReaction(reaction)};
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 2);
    rxn->initReactantMatchers();

    auto product_sets = rxn->runReactant(mol2, 0);
    TEST_ASSERT(product_sets.size() == 1);
    TEST_ASSERT(product_sets[0].size() == 2);
    // We are only interested in the product that keeps the double bond
    TEST_ASSERT(check_bond_stereo(product_sets[0][0], 0, 2, 4,
                                  Bond::BondStereo::STEREOCIS));

    // and if we don't strip direction?
    auto lmol = RWMOL_SPTR(SmilesToMol(R"(Cl/C(C)=C(/Br)F)"));
    auto product_sets2 = rxn->runReactant(lmol, 0);
    TEST_ASSERT(!product_sets2.empty());
    auto osmi = MolToSmiles(*product_sets2[0][0]);
    TEST_ASSERT(osmi == "C/C(Cl)=C/F");
  }
  {  // Reactant stereo propagated by anti-stereoatoms pair
    auto mol2 = RWMOL_SPTR(SmilesToMol(R"(Cl/C(C)=C(/Br)F)"));
    strip_bond_directions(mol2);
    TEST_ASSERT(check_bond_stereo(mol2, 2, 0, 4, Bond::BondStereo::STEREOE));

    const std::string reaction(
        R"([Cl:4][C:1]=[C:2][Br:3]>>[C:1]=[C:2].[Br:3].[Cl:4])");
    std::unique_ptr<ChemicalReaction> rxn{
        RxnSmartsToChemicalReaction(reaction)};
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 3);
    rxn->initReactantMatchers();

    auto product_sets = rxn->runReactant(mol2, 0);
    TEST_ASSERT(product_sets.size() == 1);
    TEST_ASSERT(product_sets[0].size() == 3);
    // We are only interested in the product that keeps the double bond
    TEST_ASSERT(check_bond_stereo(product_sets[0][0], 0, 2, 3,
                                  Bond::BondStereo::STEREOTRANS));

    // and if we don't strip direction?
    auto lmol = RWMOL_SPTR(SmilesToMol(R"(Cl/C(C)=C(/Br)F)"));
    auto product_sets2 = rxn->runReactant(lmol, 0);
    TEST_ASSERT(!product_sets2.empty());
    auto osmi = MolToSmiles(*product_sets2[0][0]);
    TEST_ASSERT(osmi == "C/C=C/F");
  }
}

void testGithub2547() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing Github #2547: Check kekulization issues in mdl rxn files"
      << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName;

  fName = rdbase + "/Code/GraphMol/ChemReactions/testData/kekule_reaction.rxn";
  BOOST_LOG(rdInfoLog) << "--- reading file" << std::endl;
  std::unique_ptr<ChemicalReaction> rxn(RxnFileToChemicalReaction(fName));
  BOOST_LOG(rdInfoLog) << "--- read file" << std::endl;
  ROMOL_SPTR mol("c1ccccc1"_smiles);
  rxn->initReactantMatchers();
  std::vector<ROMOL_SPTR> v{mol};
  auto prods = rxn->runReactants(v);
  TEST_ASSERT(prods.size() == 0);

  RxnOps::sanitizeRxn(*rxn);

  prods = rxn->runReactants(v);
  TEST_ASSERT(prods.size() > 1);
}

void testDblBondCrash() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing double bond stereo vanishing during a reaction" << std::endl;
  {
    const std::string reaction(R"([N;!H0:1]>>[N:1])");
    std::unique_ptr<ChemicalReaction> rxn(
        RxnSmartsToChemicalReaction(reaction));
    TEST_ASSERT(rxn);
    rxn->initReactantMatchers();
    {
      auto mol = RWMOL_SPTR(SmilesToMol("C/C=C(/C)CN"));
      std::vector<ROMOL_SPTR> v{mol};
      auto prods = rxn->runReactants(v);
      TEST_ASSERT(prods.size() == 1);
      TEST_ASSERT(prods[0].size() == 1);
      MolOps::sanitizeMol(*(RWMol *)prods[0][0].get());
      TEST_ASSERT(MolToSmiles(*prods[0][0]) == "C/C=C(/C)CN");
    }
    {
      // create an artificial situation in the molecule where bond stereo is
      // set, but neither stereoatoms nor CIP ranks are there (this can happen
      // with serialized molecules from old RDKit versions)
      auto mol = RWMOL_SPTR(SmilesToMol("C/C=C(/C)CN"));
      mol->getBondWithIdx(1)->getStereoAtoms().clear();
      for (auto atom : mol->atoms()) {
        if (atom->hasProp(common_properties::_CIPRank)) {
          atom->clearProp(common_properties::_CIPRank);
        }
      }
      std::vector<ROMOL_SPTR> v{mol};
      auto prods = rxn->runReactants(v);
      TEST_ASSERT(prods.size() == 1);
      TEST_ASSERT(prods[0].size() == 1);
      MolOps::sanitizeMol(*(RWMol *)prods[0][0].get());
      TEST_ASSERT(MolToSmiles(*prods[0][0]) == "CC=C(C)CN");
    }
    {
      SmilesParserParams ps;
      ps.sanitize = true;
      auto mol = RWMOL_SPTR(SmilesToMol("CC(=O)/N=C(\\C)c1nonc1N", ps));
      std::vector<ROMOL_SPTR> v{mol};
      auto prods = rxn->runReactants(v);
      TEST_ASSERT(prods.size() == 1);
      TEST_ASSERT(prods[0].size() == 1);
      MolOps::sanitizeMol(*(RWMol *)prods[0][0].get());
      TEST_ASSERT(MolToSmiles(*prods[0][0]) == "CC(=O)/N=C(\\C)c1nonc1N");
    }
    {
      SmilesParserParams ps;
      ps.sanitize = true;
      auto mol = RWMOL_SPTR(SmilesToMol("CC(=O)/N=C(\\C)c1nonc1N", ps));
      // make sure that sanitizing again doesn't screw things up.
      MolOps::sanitizeMol(*mol);
      std::vector<ROMOL_SPTR> v{mol};
      auto prods = rxn->runReactants(v);
      TEST_ASSERT(prods.size() == 1);
      TEST_ASSERT(prods[0].size() == 1);
      MolOps::sanitizeMol(*(RWMol *)prods[0][0].get());
      TEST_ASSERT(MolToSmiles(*prods[0][0]) == "CC(=O)/N=C(\\C)c1nonc1N");
    }
    {
      // make sure that sanitizing afterwards doesn't screw things up:
      SmilesParserParams ps;
      ps.sanitize = false;
      auto mol = RWMOL_SPTR(SmilesToMol("CC(=O)/N=C(\\C)c1nonc1N", ps));
      MolOps::sanitizeMol(*mol);
      std::vector<ROMOL_SPTR> v{mol};
      auto prods = rxn->runReactants(v);
      TEST_ASSERT(prods.size() == 1);
      TEST_ASSERT(prods[0].size() == 1);
      MolOps::sanitizeMol(*(RWMol *)prods[0][0].get());
      TEST_ASSERT(MolToSmiles(*prods[0][0]) == "CC(=O)/N=C(\\C)c1nonc1N");
    }
    {
      // another example
      const std::string reaction2(R"([N;!H0;$(N-c):1]>>[N:1]-c1cncnn1)");
      std::unique_ptr<ChemicalReaction> rxn2(
          RxnSmartsToChemicalReaction(reaction2));
      TEST_ASSERT(rxn2);
      rxn2->initReactantMatchers();
      auto mol = RWMOL_SPTR(SmilesToMol("CC(=O)/C=C(\\N)c1nonc1N"));
      std::vector<ROMOL_SPTR> v{mol};
      auto prods = rxn2->runReactants(v);
      TEST_ASSERT(prods.size() == 1);
      TEST_ASSERT(prods[0].size() == 1);
      MolOps::sanitizeMol(*(RWMol *)prods[0][0].get());
      TEST_ASSERT(MolToSmiles(*prods[0][0]) ==
                  "CC(=O)/C=C(\\N)c1nonc1Nc1cncnn1");
    }
  }
}

void testGithub3097() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github issue # 3097" << std::endl;

  std::string smi =
      "[c:6][n:7][c:8].[N:14]#[N:15]>>[C:6].[C:8][N:7]=[N+:14]=[N-:15]";
  std::unique_ptr<ChemicalReaction> rxn{RxnSmartsToChemicalReaction(smi)};
  TEST_ASSERT(rxn);
  rxn->initReactantMatchers();

  auto mol1 = RWMOL_SPTR(SmilesToMol("c1cc[nH]c1"));
  auto mol2 = RWMOL_SPTR(SmilesToMol("N#N"));
  std::vector<ROMOL_SPTR> v{mol1, mol2};
  auto prods = rxn->runReactants(v);
  // if products are not empty this is already a success
  TEST_ASSERT(prods.size() > 0);
}

void testRxnBlockRemoveHs() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing removing explicit Hs from RxnBlock"
                       << std::endl;
  std::string rxnB = R"RXN($RXN
Dummy 0
  Dummy        0123456789

  1  1
$MOL

  Dummy   01234567892D

 10 10  0  0  0  0            999 V2000
    7.0222  -11.1783    0.0000 C   0  0  0  0  0  0  0  0  0  1  0  0
    8.0615  -11.7783    0.0000 O   0  0  0  0  0  0  0  0  0  2  0  0
    7.0222   -9.6783    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0
    5.7231   -8.9283    0.0000 C   0  0  0  0  0  0  0  0  0  4  0  0
    5.7231   -7.7283    0.0000 A   0  0  0  0  0  0  0  0  0  5  0  0
    4.4242   -9.6783    0.0000 C   0  0  0  0  0  0  0  0  0  6  0  0
    4.4242  -11.1783    0.0000 C   0  0  0  0  0  0  0  0  0  7  0  0
    3.3849  -11.7783    0.0000 A   0  0  0  0  0  0  0  0  0  8  0  0
    5.7231  -11.9283    0.0000 N   0  0  0  0  0  0  0  0  0  9  0  0
    5.7231  -13.1094    0.0000 H   0  0
  1  2  2  0  0  0  8
  1  3  1  0  0  0  8
  3  4  2  0  0  0  8
  4  5  1  0  0  0  2
  4  6  1  0  0  0  8
  6  7  2  0  0  0  8
  7  8  1  0  0  0  2
  7  9  1  0  0  0  8
  9  1  1  0  0  0  8
  9 10  1  0
M  SUB  1   9   2
M  END
$MOL

  Dummy   01234567892D

  9  9  0  0  0  0            999 V2000
   17.0447  -11.1783    0.0000 C   0  0  0  0  0  0  0  0  0  1  0  0
   18.0840  -11.7783    0.0000 O   0  0  0  0  0  0  0  0  0  2  0  0
   17.0447   -9.6783    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0
   15.7457   -8.9283    0.0000 C   0  0  0  0  0  0  0  0  0  4  0  0
   15.7457   -7.7283    0.0000 A   0  0  0  0  0  0  0  0  0  5  0  0
   14.4467   -9.6783    0.0000 C   0  0  0  0  0  0  0  0  0  6  0  0
   14.4467  -11.1783    0.0000 C   0  0  0  0  0  0  0  0  0  7  0  0
   13.4074  -11.7783    0.0000 A   0  0  0  0  0  0  0  0  0  8  0  0
   15.7457  -11.9283    0.0000 N   0  0  0  0  0  0  0  0  0  9  0  0
  1  2  1  0  0  0  8
  1  3  1  0  0  0  8
  3  4  2  0  0  0  8
  4  5  1  0  0  0  2
  4  6  1  0  0  0  8
  6  7  2  0  0  0  8
  7  8  1  0  0  0  2
  7  9  1  0  0  0  8
  9  1  2  0  0  0  8
M  END
)RXN";

  ROMOL_SPTR mol(
      static_cast<ROMol *>(SmilesToMol("c1(=O)nc([Cl])cc([F])[nH]1")));
  std::vector<ROMOL_SPTR> v{mol};

  {
    std::unique_ptr<ChemicalReaction> rxn(RxnBlockToChemicalReaction(rxnB));
    TEST_ASSERT(rxn.get());
    rxn->initReactantMatchers();

    auto prods = rxn->runReactants(v);
    // if the explicit hydrogen is not removed and the reactant template
    // is not sanitized, the reactant template is not aromatic and our
    // aromatic reactant won't match
    TEST_ASSERT(prods.size() == 0);
  }
  {
    std::unique_ptr<ChemicalReaction> rxn(
        RxnBlockToChemicalReaction(rxnB, true, true));
    TEST_ASSERT(rxn.get());
    rxn->initReactantMatchers();

    auto prods = rxn->runReactants(v);
    TEST_ASSERT(prods.size() == 2);
  }
}

void testGithub3078() {
  BOOST_LOG(rdInfoLog) << "--------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Reaction-edge stereo bonds should default to STEREONONE" << std::endl;

  const std::string reaction(
      R"([C:1][C:2]1[C:3][N:4]1>>[C:1][C:2]1=[C:3][N:4]1)");
  std::unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction(reaction));
  TEST_ASSERT(rxn);
  rxn->initReactantMatchers();
  std::vector<ROMOL_SPTR> reactant{"CC1CN1"_smiles};

  auto prods = rxn->runReactants(reactant);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);

  const auto bond = prods[0][0]->getBondWithIdx(1);
  TEST_ASSERT(bond->getBondType() == Bond::DOUBLE);
  TEST_ASSERT(bond->getStereo() == Bond::STEREONONE);

  // Make sure the temporary mark set in the reaction has been removed.
  TEST_ASSERT(!bond->hasProp("_UnknownStereoRxnBond"));
}

void testGithub4162() {
  const std::string reaction(
      R"([C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3])");
  std::unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction(reaction));
  std::unique_ptr<ChemicalReaction> rxnCopy(new ChemicalReaction(*rxn));
  RxnOps::sanitizeRxn(*rxn);
  RxnOps::sanitizeRxn(*rxnCopy);
  std::string pkl;
  ReactionPickler::pickleReaction(*rxn, pkl);
  std::unique_ptr<ChemicalReaction> rxnFromPickle(new ChemicalReaction(pkl));
  RxnOps::sanitizeRxn(*rxnFromPickle);
  ReactionPickler::pickleReaction(*rxnFromPickle, pkl);
  rxnFromPickle.reset(new ChemicalReaction(pkl));
  RxnOps::sanitizeRxn(*rxnFromPickle);
}

void testGithub4114() {
  BOOST_LOG(rdInfoLog) << "--------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Reactions don't propagate bond properties"
                       << std::endl;

  const std::string DUMMY_PROP{"dummy_prop"};
  std::vector<ROMOL_SPTR> mol{"ClC(N=O)=C=C=C=C(Br)N=S"_smiles};

  // React the left side of the mol
  const std::string reaction(
      R"([Cl:1][C:2]([N:3]=[O:4])=[C:5]=[C:6]>>[F:1][Si:2]([P:3]~[S:4])(~[Si:5]=[Si:6]))");
  std::unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction(reaction));
  TEST_ASSERT(rxn);
  rxn->initReactantMatchers();

  // Set dummy properties on all bonds; mark the central double bonds
  // as EITHERDOUBLE
  for (unsigned int i = 0; i < mol[0]->getNumBonds(); ++i) {
    auto bnd = mol[0]->getBondWithIdx(i);
    bnd->setProp(DUMMY_PROP, i);

    if (i >= 3 && i <= 6) {
      TEST_ASSERT(bnd->getBondType() == Bond::DOUBLE);
      bnd->setBondDir(Bond::EITHERDOUBLE);
    }
  }

  auto prods = rxn->runReactants(mol);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);

  // Make sure the reaction did not alter the graph
  for (const auto &at : prods[0][0]->atoms()) {
    auto mapno = at->getProp<unsigned int>(common_properties::reactantAtomIdx);
    TEST_ASSERT(mapno == at->getIdx());
  }

  for (const auto &rBnd : mol[0]->bonds()) {
    auto pBnd = prods[0][0]->getBondBetweenAtoms(rBnd->getBeginAtomIdx(),
                                                 rBnd->getEndAtomIdx());

    // Bonds 0, 1 and 4 are overridden in the product template, and should
    // therefore override the properties.
    // Bonds 2 and 3 are "null bonds", and should get properties copied over.
    // Bonds > 4 are reactant bonds, and should keep their properties.
    auto rBndIdx = rBnd->getIdx();
    if (rBndIdx == 0 || rBndIdx == 1 || rBndIdx == 4) {
      TEST_ASSERT(!pBnd->hasProp(DUMMY_PROP));
    } else {
      unsigned int dummy_prop;
      TEST_ASSERT(pBnd->getPropIfPresent(DUMMY_PROP, dummy_prop));
      TEST_ASSERT(dummy_prop == rBndIdx);
    }

    // Bond 3 is a "null bond", and should preserve EITHERDOUBLE.
    // Bond 4 is specified in the reaction, and should not preserve
    // EITHERDOUBLE. Bonds 5 & 6 are not included in the reaction, so they
    // should keep EITHERDOUBLE dir.
    if (rBndIdx == 3 || rBndIdx == 5 || rBndIdx == 6) {
      TEST_ASSERT(pBnd->getBondType() == Bond::DOUBLE);
      TEST_ASSERT(pBnd->getBondDir() == Bond::EITHERDOUBLE);
    } else {
      TEST_ASSERT(pBnd->getBondDir() == Bond::NONE);
    }
  }
}

void testGithub4183() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github #4183: Reading a rxn file in v3000 "
                          "format that contains agents"
                       << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName;
  fName = rdbase + "/Code/GraphMol/ChemReactions/testData/v3k_with_agents.rxn";
  ChemicalReaction *rxn = RxnFileToChemicalReaction(fName);
  TEST_ASSERT(rxn);

  TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
  TEST_ASSERT(rxn->getNumProductTemplates() == 1);
  TEST_ASSERT(rxn->getNumAgentTemplates() == 3);

  delete rxn;
}

void testGithub4410() {
  BOOST_LOG(rdInfoLog) << "--------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Github #4410: wrong double bond stereochemistry"
                       << std::endl;

  const std::string reaction(
      R"([N:4][C:6](/[Cl:9])=[C:7](\[Cl:10])[C:11]>>[Br:4][CX3:6](\[Cl:9])=[CX3:7](/[Cl:10])[C:11])");
  std::unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction(reaction));
  TEST_ASSERT(rxn);
  rxn->initReactantMatchers();

  std::vector<ROMOL_SPTR> mol{R"(C\C(Cl)=C(\N)Cl)"_smiles};
  auto prods = rxn->runReactants(mol);
  TEST_ASSERT(prods.size() == 1);
  TEST_ASSERT(prods[0].size() == 1);

  auto dblBnd = prods[0][0]->getBondBetweenAtoms(1, 3);

  TEST_ASSERT(dblBnd->getBondType() == Bond::DOUBLE);
  TEST_ASSERT(dblBnd->getStereoAtoms() == std::vector<int>({2, 4}));
  TEST_ASSERT(dblBnd->getStereo() == Bond::STEREOTRANS);

  TEST_ASSERT(MolToSmiles(*prods[0][0]) == R"(C/C(Cl)=C(\Cl)Br)");
}

void testMultiTemplateRxnQueries() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing multi-template isMoleculeXOfReaction queries"
                       << std::endl;

  std::string rxn_smarts =
      "[S;v1&H0,v2&H1:1].[S;v2;H0,H1:2][S;v2;H0,H1:3]>>[S:3].[S:1][S:2]";
  ChemicalReaction *rxn = RxnSmartsToChemicalReaction(rxn_smarts);
  TEST_ASSERT(rxn->getNumReactantTemplates() == 2);
  TEST_ASSERT(rxn->getNumProductTemplates() == 2);
  rxn->initReactantMatchers();
  unsigned int nWarn, nError;
  TEST_ASSERT(rxn->validate(nWarn, nError, false));
  TEST_ASSERT(nWarn == 0 && nError == 0);

  ROMol *reactant = SmilesToMol("SC1=CC(CSSCC2=CC=CC=C2)=CC=C1");
  ROMol *product = reactant;
  ROMol *neither = SmilesToMol("c1ccccc1");

  std::vector<unsigned int> which;
  bool is_reactant = isMoleculeReactantOfReaction(*rxn, *reactant, which);
  TEST_ASSERT(is_reactant);
  TEST_ASSERT(which == std::vector<unsigned int>({0, 1}));
  unsigned int first_match;
  is_reactant = isMoleculeReactantOfReaction(*rxn, *reactant, first_match);
  TEST_ASSERT(is_reactant);
  TEST_ASSERT(first_match == 0);
  is_reactant = isMoleculeReactantOfReaction(*rxn, *reactant);
  TEST_ASSERT(is_reactant);

  is_reactant = isMoleculeReactantOfReaction(*rxn, *neither, which);
  TEST_ASSERT(!is_reactant);
  TEST_ASSERT(which.empty());
  is_reactant = isMoleculeReactantOfReaction(*rxn, *neither, first_match);
  TEST_ASSERT(!is_reactant);
  TEST_ASSERT(first_match == 2);
  is_reactant = isMoleculeReactantOfReaction(*rxn, *neither);
  TEST_ASSERT(!is_reactant);

  bool is_product = isMoleculeProductOfReaction(*rxn, *product, which);
  TEST_ASSERT(is_product);
  TEST_ASSERT(which == std::vector<unsigned int>({0, 1}));
  is_product = isMoleculeProductOfReaction(*rxn, *product, first_match);
  TEST_ASSERT(is_product);
  TEST_ASSERT(first_match == 0);
  is_product = isMoleculeProductOfReaction(*rxn, *product);
  TEST_ASSERT(is_product);

  is_product = isMoleculeProductOfReaction(*rxn, *neither, which);
  TEST_ASSERT(!is_product);
  TEST_ASSERT(which.empty());
  is_product = isMoleculeProductOfReaction(*rxn, *neither, first_match);
  TEST_ASSERT(!is_product);
  TEST_ASSERT(first_match == 2);
  is_product = isMoleculeProductOfReaction(*rxn, *neither);
  TEST_ASSERT(!is_product);

  delete reactant;
  delete neither;
  delete rxn;
}

void testChemicalReactionCopyAssignment() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing ChemicalReaction copy assignment operator"
                       << std::endl;

  std::string rxn_smarts1 =
      "[C;$(C=O):1][OH1].[N;$(N[#6]);!$(N=*);!$([N-]);!$(N#*);!$([ND3]);!$([ND4]);!$(N[O,N]);!$(N[C,S]=[S,O,N]):2]>>[C:1][N+0:2]";
  ChemicalReaction *rxn1 = RxnSmartsToChemicalReaction(rxn_smarts1);
  rxn1->setImplicitPropertiesFlag(true);
  rxn1->initReactantMatchers();
  unsigned int nWarn, nError;
  TEST_ASSERT(rxn1->validate(nWarn, nError, false));
  TEST_ASSERT(nWarn == 0 && nError == 0);

  std::string rxn_smarts2 = "[O:1]>>[N:1]";
  ChemicalReaction *rxn2 = RxnSmartsToChemicalReaction(rxn_smarts2);

  *rxn2 = *rxn1;

  // Check we copied the base class members
  TEST_ASSERT(rxn2->getPropList() == rxn1->getPropList());

  // Check we copied the flags
  TEST_ASSERT(rxn2->getImplicitPropertiesFlag());
  TEST_ASSERT(rxn2->isInitialized());

  // Check we copied the reactant/product templates
  TEST_ASSERT(rxn2->getNumReactantTemplates() == 2);
  TEST_ASSERT(rxn2->getNumProductTemplates() == 1);
  MOL_SPTR_VECT::const_iterator it1 = rxn1->beginReactantTemplates();
  MOL_SPTR_VECT::const_iterator it2 = rxn2->beginReactantTemplates();
  MOL_SPTR_VECT::const_iterator end_it1 = rxn1->endReactantTemplates();
  while (it1 != end_it1) {
    TEST_ASSERT(MolToSmiles(**it1) == MolToSmiles(**it2));
    ++it1;
    ++it2;
  }
  it1 = rxn1->beginProductTemplates();
  it2 = rxn2->beginProductTemplates();
  end_it1 = rxn1->endProductTemplates();
  while (it1 != end_it1) {
    TEST_ASSERT(MolToSmiles(**it1) == MolToSmiles(**it2));
    ++it1;
    ++it2;
  }

  // Check that the reactions don't share resources
  const RWMol &rxn1_reactant = *rxn1->getReactants().at(0);
  const_cast<RWMol &>(rxn1_reactant).clear();
  ROMOL_SPTR rxn2_reactant = rxn2->getReactants().at(0);
  TEST_ASSERT(rxn2_reactant->getNumAtoms() > 0);

  // Check the reaction works
  MOL_SPTR_VECT reactants;
  reactants.emplace_back(SmilesToMol("CC(=O)O"));
  reactants.emplace_back(SmilesToMol("CCN"));
  std::vector<MOL_SPTR_VECT> products = rxn2->runReactants(reactants);
  TEST_ASSERT(MolToSmiles(*products[0][0]) == "CCNC(C)=O");

  delete rxn1;
  delete rxn2;
}

void testGithub6138() {
  // Pickling reactions removed some of their properties set after reaction
  // initialization
  auto rxn_smarts = "[c:1]:[n&H1&+0&D2:3]:[n:2]>>[c:1]:[3n&H0&+0&D3:3]:[2n:2]";
  std::unique_ptr<ChemicalReaction> rxn(
      RxnSmartsToChemicalReaction(rxn_smarts));
  ROMOL_SPTR mol("c1cn[nH]c1"_smiles);
  rxn->initReactantMatchers();
  MOL_SPTR_VECT reacts;
  reacts.push_back(mol);
  auto prods = rxn->runReactants(reacts);
  std::string pkl;
  ReactionPickler::pickleReaction(*rxn, pkl);
  std::unique_ptr<ChemicalReaction> lrxn(new ChemicalReaction());
  ReactionPickler::reactionFromPickle(pkl, lrxn.get());
  auto prods2 = lrxn->runReactants(reacts);
  auto s1 = MolToSmiles(*prods[0][0]);
  auto s2 = MolToSmiles(*prods2[0][0]);
  TEST_ASSERT(s1 == s2);
}

void testReactionWithChiralAgent() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing introduction of new atoms with chirality"
                       << std::endl;

  {  // a reaction with a chiral agent - v3000
    std::string rdbase = getenv("RDBASE");
    std::string fName;

    fName =
        rdbase +
        "/Code/GraphMol/ChemReactions/testData/testRXNChiralityAgentV3000.rxn";
    ChemicalReaction *rxn =
        RxnFileToChemicalReaction(fName, false, false, false);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 1);
    auto outputRxn = ChemicalReactionToRxnSmiles(*rxn);
    TEST_ASSERT(
        outputRxn ==
        "[CH:1]([F:2])([CH3:3])[CH2:4][CH2:5][Br:6]>CC[C@H](C)Cl>[CH:1]([F:2])([CH3:3])[CH2:4][CH2:5][CH2:7][CH:8]([CH3:9])[Cl:10]");

    BOOST_LOG(rdInfoLog) << ChemicalReactionToRxnSmiles(*rxn) << std::endl;

    delete rxn;
  }

  {  // a reaction with a chiral agent - v2000
    std::string rdbase = getenv("RDBASE");
    std::string fName;

    fName =
        rdbase +
        "/Code/GraphMol/ChemReactions/testData/testRXNChiralityAgentV2000.rxn";
    ChemicalReaction *rxn =
        RxnFileToChemicalReaction(fName, false, false, false);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);
    TEST_ASSERT(rxn->getNumAgentTemplates() == 1);
    auto outputRxn = ChemicalReactionToRxnSmiles(*rxn);
    TEST_ASSERT(
        outputRxn ==
        "[CH:1]([F:2])([CH3:3])[CH2:4][CH2:5][Br:6]>CC[C@H](C)Cl>[CH:1]([F:2])([CH3:3])[CH2:4][CH2:5][CH2:7][CH:8]([CH3:9])[Cl:10]");

    BOOST_LOG(rdInfoLog) << ChemicalReactionToRxnSmiles(*rxn) << std::endl;

    delete rxn;
  }
}

void testGithub5890() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Github Issue 5890: Testing reaction with radicals"
                       << std::endl;

  {
    std::string rdbase = getenv("RDBASE");
    std::string fName;

    fName =
        rdbase +
        "/Code/GraphMol/ChemReactions/testData/v3k.radicals.rxn";
    ChemicalReaction *rxn =
        RxnFileToChemicalReaction(fName, false, false, false);
    TEST_ASSERT(rxn);
    TEST_ASSERT(rxn->getNumReactantTemplates() == 1);
    TEST_ASSERT(rxn->getNumProductTemplates() == 1);

    std::string pkl;
    ReactionPickler::pickleReaction(rxn, pkl);
    delete rxn;
    rxn = new ChemicalReaction();
    ReactionPickler::reactionFromPickle(pkl, rxn);

    auto outputRxn = ChemicalReactionToRxnSmiles(*rxn);
    BOOST_LOG(rdInfoLog) << outputRxn << std::endl;
    TEST_ASSERT(outputRxn == "[CH]1[CH][CH]1>>C");

    delete rxn;
  }
}

int main() {
  RDLog::InitLogs();

  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Testing Chemical Reactions \n";
  test1Basics();
  test2SimpleReactions();
  test3RingFormation();
  test4MultipleProducts();
  test5Salts();
  test6DaylightParser();
  test7MDLParser();
  test8Validation();
  test9ProductQueries();
  test10ChiralityDaylight();
  test11ChiralityRxn();
  test12DoubleBondStereochem();
  test13Issue1748846();
  test14Issue1804420();
  test15Issue1882749();
  test16Exceptions();
  test17Issue1920627();
  test18PropertyTransfer();
  test19Issue2050085();
  test20BondQueriesInProduct();
  test21Issue2540021();
  test22DotsToRemoveBonds();
  test23Pickling();
  test24AtomFlags();
  test25Conformers();
  test26V3000MDLParser();
  test27SmartsWriter();
  test28RxnDepictor();
  test29RxnWriter();
  test30ReactProdQueries();
  test31Issue3140490();
  test32Replacements();

  test33ReactingAtoms1();
  test34ReactingAtoms2();
  test35ParensInReactants1();
  test36ParensInReactants2();
  test37ProtectOption();

  test38AddRecursiveQueriesToReaction();

  test39InnocentChiralityLoss();
  test40AgentsInSmarts();
  test41Github233();
  test42ReactionSmiles();
  test44Github290();
  test45SmilesWriter();
  test46Agents();
  test47TestReactionMoleculeConversion();
  test48ParensInProducts1();
  test49ParensInProducts2();
  test50RNXFileParserWithEmptyAgentColumn();
  test51RNXSmilesFromPatentData();
  test52RedundantProductMappingNumbersAndRunReactants();
  test53ReactionSubstructureMatching();
  test54RedundantProductMappingNumbersAndRSChirality();
  test55RedundantProductMappingNumbersAndEZStereochemistry();
  test56TestOldPickleVersion();
  test57IntroductionOfNewChiralCenters();
  test58MolFileValueRoundTrip();
  test59ReactionCanonicalization();
  test60RunSingleReactant();
  test61Github685();
  test62Github975();
  test63CopyConstructor();
  test43Github243();
  test64Github1266();
  test65SanitizeUnmappedHs();
  test66SanitizeMappedHs();
  test67SanitizeMappedHsInReactantAndProd();
  test68MappedHToHeavy();
  test69Github1387();
  test70Github1544();
  testSanitizeException();
  testReactionProperties();
  testGithub1950();
  testGithub1869();
  testGithub1269();
  testGithub1868();
  testBondStereoGistExamples();
  testStereoBondIsomerization();
  testOtherBondStereo();
  testGithub2547();
  testGithub3097();
  testDblBondCrash();
  testRxnBlockRemoveHs();
  testGithub3078();
  testGithub4162();
  testGithub4114();
  testGithub4183();
  testGithub4410();
  testMultiTemplateRxnQueries();
  testChemicalReactionCopyAssignment();
  testGithub6138();
  testReactionWithChiralAgent();
  testGithub5890();

  BOOST_LOG(rdInfoLog)
      << "*******************************************************\n";
  return (0);
}
