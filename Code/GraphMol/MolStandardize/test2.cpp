//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// file used to test memory leakage

#include "MolStandardize.h"
#include "Validate.h"
#include "Metal.h"
#include "Normalize.h"
#include "Charge.h"
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/Invariant.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/ROMol.h>

#include <iostream>

using namespace RDKit;
using namespace MolStandardize;

void test1() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n test1" << std::endl;

  MolStandardize::CleanupParameters params;

  {
    //	std::shared_ptr<RWMol> m1( SmilesToMol(smi1) );
    RWMOL_SPTR m = "C1=CC=CC=C1"_smiles;
    RWMOL_SPTR res(cleanup(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "c1ccccc1");
  }

  //	std::string smi1 = "CCC(=O)O[Na]";
  ////	std::shared_ptr<RWMol> m1( SmilesToMol(smi1) );
  //	RWMol *m1 = SmilesToMol(smi1);
  //	cleanup(*m1, params);
  //	TEST_ASSERT(MolToSmiles(*m1) == "CCC(=O)[O-].[Na+]");
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testMetal() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n test metal" << std::endl;
  MetalDisconnector md;

  RWMOL_SPTR m = "CCC(=O)O[Na]"_smiles;
  TEST_ASSERT(m);
  md.disconnect(*m);
  TEST_ASSERT(MolToSmiles(*m) == "CCC(=O)[O-].[Na+]");
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testValidate() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n test validate"
                       << std::endl;
  // testing RDKitDefault
  {
    RDKitValidation vm;
    std::string smi = "CO(C)C";
    RWMOL_SPTR m(SmilesToMol(smi, 0, false));
    std::vector<ValidationErrorInfo> errout = vm.validate(*m, true);
    for (auto &query : errout) {
      std::string msg = query.what();
      TEST_ASSERT(
          msg ==
          "INFO: [ValenceValidation] Explicit valence for atom # 1 O, 3, "
          "is greater than permitted");
    }
  }
  //**************************
  {
    MolVSValidation vm;
    std::string smi = "O=C([O-])c1ccccc1";
    RWMOL_SPTR m(SmilesToMol(smi, 0, false));
    std::vector<ValidationErrorInfo> errout = vm.validate(*m, true);
    for (auto &query : errout) {
      std::string msg = query.what();
      TEST_ASSERT(
          msg ==
          "INFO: [NeutralValidation] Not an overall neutral system (-1)");
    }
  }
  // ************************
  {
    std::vector<unsigned int> atoms = {6, 7, 8};
    std::vector<std::shared_ptr<Atom>> atomList;

    for (auto &atom : atoms) {
      std::shared_ptr<Atom> a(new Atom(atom));
      atomList.push_back(a);
    }

    AllowedAtomsValidation vm(atomList);
    RWMOL_SPTR m = "CC(=O)CF"_smiles;
    std::vector<ValidationErrorInfo> errout = vm.validate(*m, true);
    for (auto &query : errout) {
      std::string msg = query.what();
      TEST_ASSERT(
          msg ==
          "INFO: [AllowedAtomsValidation] Atom F is not in allowedAtoms list");
    }
  }
  //********************************
  {
    std::vector<unsigned int> atoms = {9, 17, 35};
    std::vector<std::shared_ptr<Atom>> atomList;

    for (auto &atom : atoms) {
      std::shared_ptr<Atom> a(new Atom(atom));
      atomList.push_back(a);
    }

    DisallowedAtomsValidation vm(atomList);
    RWMOL_SPTR m = "CC(=O)CF"_smiles;
    std::vector<ValidationErrorInfo> errout = vm.validate(*m, true);
    for (auto &query : errout) {
      std::string msg = query.what();
      TEST_ASSERT(msg ==
                  "INFO: [DisallowedAtomsValidation] Atom F is in "
                  "disallowedAtoms list");
    }
  }
  //********************************
  // testing MolVSValidation fragmentValidation
  // FragmentValidation should identify 1,2-dichloroethane.
  {
    MolVSValidation vm;
    std::string smi = "ClCCCl.c1ccccc1O";
    RWMOL_SPTR m(SmilesToMol(smi, 0, false));
    std::vector<ValidationErrorInfo> errout = vm.validate(*m, true);
    for (auto &query : errout) {
      std::string msg = query.what();
      TEST_ASSERT(msg ==
                  "INFO: [FragmentValidation] 1,2-dichloroethane is present");
    }
    BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
  }
}

void testCharge() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n test charge" << std::endl;

  Reionizer reionizer;

  // Test table salt.
  {
    RWMOL_SPTR m = "[Na].[Cl]"_smiles;
    ROMOL_SPTR reionized(reionizer.reionize(*m));
    TEST_ASSERT(MolToSmiles(*reionized) == "[Cl-].[Na+]");
  }
  //*******************************
  MolStandardize::CleanupParameters params;
  // initialize CleanupParameters with preferOrganic=true
  MolStandardize::CleanupParameters params_preferorg;
  params_preferorg.preferOrganic = true;

  // Test neutralization of ionized acids and bases.
  {
    RWMOL_SPTR m = "C(C(=O)[O-])(Cc1n[n-]nn1)(C[NH3+])(C[N+](=O)[O-])"_smiles;
    RWMOL_SPTR res(MolStandardize::chargeParent(*m, params));
    TEST_ASSERT(MolToSmiles(*res) == "NCC(Cc1nn[nH]n1)(C[N+](=O)[O-])C(=O)O");
  }
  //**********************************
  // Testing MolStandardize::reionize
  {
    RWMOL_SPTR m = "C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O"_smiles;
    RWMOL_SPTR res(MolStandardize::reionize(m.get(), params));
    TEST_ASSERT(MolToSmiles(*res) == "O=S(O)c1ccc(S(=O)(=O)[O-])cc1");
    BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
  }
}

void testNormalize() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n test normalize"
                       << std::endl;

  Normalizer normalizer;

  // Test sulfoxide normalization.
  {
    RWMOL_SPTR m = "CS(C)=O"_smiles;
    ROMOL_SPTR normalized(normalizer.normalize(*m));
    TEST_ASSERT(MolToSmiles(*normalized) == "C[S+](C)[O-]");
  }

  // normalize sulfone.
  {
    RWMOL_SPTR m = "C[S+2]([O-])([O-])C"_smiles;
    ROMOL_SPTR normalized(normalizer.normalize(*m));
    TEST_ASSERT(MolToSmiles(*normalized) == "CS(C)(=O)=O");
    BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
  }
}

int main() {
  test1();  // cleanup test
  testMetal();
  testValidate();
  testCharge();  // TODO
  testNormalize();
  return 0;
}
