//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Validate.h"
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Chirality.h>

#include <iostream>

using namespace RDKit;
using namespace std;
using namespace MolStandardize;

void testRDKitValidation() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing RDKit validation"
                       << std::endl;

  string smi1, smi2, smi3, smi4;
  RDKitValidation vm;

  // testing RDKitDefault
  smi1 = "CO(C)C";
  unique_ptr<ROMol> m1(SmilesToMol(smi1, 0, false));
  vector<ValidationErrorInfo> errout1 = vm.validate(*m1, true);
  for (const auto &msg : errout1) {
    TEST_ASSERT(msg ==
                "INFO: [ValenceValidation] Explicit valence for atom # 1 O, 3, "
                "is greater than permitted");
  }

  // testing for molecule with no atoms
  smi2 = "";
  unique_ptr<ROMol> m2(SmilesToMol(smi2, 0, false));
  vector<ValidationErrorInfo> errout2 = vm.validate(*m2, true);
  for (const auto &msg : errout2) {
    TEST_ASSERT(msg == "ERROR: [NoAtomValidation] Molecule has no atoms");
  }

  // testing molecule with multiple valency errors
  smi3 = "CO(C)CCN(=O)=O";
  unique_ptr<ROMol> m3(SmilesToMol(smi3, 0, false));
  vector<ValidationErrorInfo> errout3 = vm.validate(*m3, true);
  std::vector<string> msgs1;
  std::vector<string> ans1 = {
      "INFO: [ValenceValidation] Explicit valence for atom # 1 O, 3, is "
      "greater than permitted",
      "INFO: [ValenceValidation] Explicit valence for atom # 5 N, 5, is "
      "greater than permitted"};
  for (const auto &msg : errout3) {
    msgs1.push_back(msg);
  }
  TEST_ASSERT(msgs1 == ans1);

  // testing molecule with multiple valency errors and only outputting
  // first error
  smi4 = "CO(C)CCN(=O)=O";
  unique_ptr<ROMol> m4(SmilesToMol(smi4, 0, false));
  vector<ValidationErrorInfo> errout4 = vm.validate(*m4, false);
  std::vector<string> msgs2;
  std::vector<string> ans2 = {
      "INFO: [ValenceValidation] Explicit valence for atom # 1 O, 3, is "
      "greater than permitted"};
  for (const auto &msg : errout4) {
    msgs2.push_back(msg);
  }
  TEST_ASSERT(msgs2 == ans2);
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testMolVSValidation() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing MolVS validation"
                       << std::endl;
  string smi1, smi2, smi3, smi4, smi5, smi6;
  MolVSValidation vm;

  // testing MolVSDefault
  // testing for molecule with no atoms
  smi1 = "";
  unique_ptr<ROMol> m1(SmilesToMol(smi1, 0, false));
  vector<ValidationErrorInfo> errout1 = vm.validate(*m1, true);
  for (const auto &msg : errout1) {
    TEST_ASSERT(msg == "ERROR: [NoAtomValidation] Molecule has no atoms");
  }

  smi2 = "O=C([O-])c1ccccc1";
  unique_ptr<ROMol> m2(SmilesToMol(smi2, 0, false));
  vector<ValidationErrorInfo> errout2 = vm.validate(*m2, true);
  for (const auto &msg : errout2) {
    TEST_ASSERT(msg ==
                "INFO: [NeutralValidation] Not an overall neutral system (-1)");
  }

  smi3 = "CN=[NH+]CN=N";
  unique_ptr<ROMol> m3(SmilesToMol(smi3, 0, false));
  vector<ValidationErrorInfo> errout3 = vm.validate(*m3, true);
  for (const auto &msg : errout3) {
    TEST_ASSERT(
        msg ==
        "INFO: [NeutralValidation] Not an overall neutral system (+1)");  // fix
                                                                          // to
                                                                          // show
                                                                          // +
                                                                          // sign
  }

  smi4 = "[13CH4]";
  unique_ptr<ROMol> m4(SmilesToMol(smi4, 0, false));
  vector<ValidationErrorInfo> errout4 = vm.validate(*m4, true);
  for (const auto &msg : errout4) {
    TEST_ASSERT(msg ==
                "INFO: [IsotopeValidation] Molecule contains isotope 13C");
  }

  smi5 = "[2H]C(Cl)(Cl)Cl";
  unique_ptr<ROMol> m5(SmilesToMol(smi5, 0, false));
  vector<ValidationErrorInfo> errout5 = vm.validate(*m5, true);
  for (const auto &msg : errout5) {
    TEST_ASSERT(msg ==
                "INFO: [IsotopeValidation] Molecule contains isotope 2H");
  }

  smi6 = "[2H]OC([2H])([2H])[2H]";
  unique_ptr<ROMol> m6(SmilesToMol(smi6, 0, false));
  vector<ValidationErrorInfo> errout6 = vm.validate(*m6, true);
  for (const auto &msg : errout6) {
    TEST_ASSERT(msg ==
                "INFO: [IsotopeValidation] Molecule contains isotope 2H");
  }

  std::string smi7 = "COc1cccc(C=N[N-]C(N)=O)c1[O-].O.O.O.O=[U+2]=O";
  unique_ptr<ROMol> m7(SmilesToMol(smi7, 0, false));
  vector<ValidationErrorInfo> errout7 = vm.validate(*m7, true);
  TEST_ASSERT(errout7.size() != 0);
  for (const auto &msg : errout7) {
    TEST_ASSERT(msg == "INFO: [FragmentValidation] water/hydroxide is present");
  }

  std::string smi8 = "CC(=O)O.NCC(=O)NCCCCCCCCCCNC(=O)CN";
  unique_ptr<ROMol> m8(SmilesToMol(smi8, 0, false));
  vector<ValidationErrorInfo> errout8 = vm.validate(*m8, true);
  TEST_ASSERT(errout8.size() != 0);
  for (const auto &msg : errout8) {
    TEST_ASSERT(msg ==
                "INFO: [FragmentValidation] acetate/acetic acid is present");
  }

  std::string smi9 = "N#CC(Br)(Br)C#N.[Br-].[K+]";
  unique_ptr<ROMol> m9(SmilesToMol(smi9, 0, false));
  vector<ValidationErrorInfo> errout9 = vm.validate(*m9, true);
  std::vector<std::string> ans = {
      "INFO: [FragmentValidation] bromine is present",
      "INFO: [FragmentValidation] potassium is present"};
  TEST_ASSERT(errout9.size() == ans.size());
  for (size_t i = 0; i < errout9.size(); ++i) {
    TEST_ASSERT(errout9[i] == ans[i]);
  }

  std::string smi10 = "C1COCCO1.O=C(NO)NO";
  unique_ptr<ROMol> m10(SmilesToMol(smi10, 0, false));
  vector<ValidationErrorInfo> errout10 = vm.validate(*m10, true);
  std::vector<std::string> ans10 = {
      "INFO: [FragmentValidation] 1,2-dimethoxyethane is present",
      "INFO: [FragmentValidation] 1,4-dioxane is present"};
  TEST_ASSERT(errout10.size() == ans10.size());
  for (size_t i = 0; i < errout10.size(); ++i) {
    TEST_ASSERT(errout10[i] == ans10[i]);
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testMolVSOptions() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing MolVS Options"
                       << std::endl;
  vector<std::shared_ptr<ValidationMethod>> validations = {
      std::make_shared<IsotopeValidation>()};
  MolVSValidation vm(validations);

  // testing MolVSDefault
  // testing for molecule with no atoms
  string smi1 = "";
  unique_ptr<ROMol> m1(SmilesToMol(smi1, 0, false));
  vector<ValidationErrorInfo> errout1 = vm.validate(*m1, true);
  // for (const auto &msg : errout1) {
  //   TEST_ASSERT(msg == "ERROR: [NoAtomValidation] Molecule has no atoms");
  // }
  TEST_ASSERT(errout1.empty());

  string smi2 = "O=C([O-])c1ccccc1";
  unique_ptr<ROMol> m2(SmilesToMol(smi2, 0, false));
  vector<ValidationErrorInfo> errout2 = vm.validate(*m2, true);
  // for (const auto &msg : errout2) {
  //   TEST_ASSERT(
  //     msg == "INFO: [NeutralValidation] Not an overall neutral system (-1)");
  // }
  TEST_ASSERT(errout2.empty());

  // test strict option of IsotopeValidation
  IsotopeValidation isotopeValidation(true);

  string smi3 = "[13CH3]";
  unique_ptr<ROMol> m3(SmilesToMol(smi3, 0, false));
  vector<ValidationErrorInfo> errout3 = isotopeValidation.validate(*m3, true);
  TEST_ASSERT(errout3.empty());

  string smi4 = "[3CH3]";
  unique_ptr<ROMol> m4(SmilesToMol(smi4, 0, false));
  vector<ValidationErrorInfo> errout4 = isotopeValidation.validate(*m4, true);
  TEST_ASSERT(errout4.size() == 1);
  TEST_ASSERT(
      errout4[0] ==
      "ERROR: [IsotopeValidation] The molecule contains an unknown isotope: 3C");

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testAllowedAtomsValidation() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing AllowedAtoms validation"
      << std::endl;

  //	std::vector<string> atoms = {"C", "N", "O"};
  std::vector<unsigned int> atoms = {6, 7, 8};
  std::vector<shared_ptr<Atom>> atomList;

  for (auto &atom : atoms) {
    shared_ptr<Atom> a(new Atom(atom));
    atomList.push_back(a);
  }

  AllowedAtomsValidation vm(atomList);
  std::string smi1;

  smi1 = "CC(=O)CF";
  unique_ptr<ROMol> m1(SmilesToMol(smi1));
  vector<ValidationErrorInfo> errout1 = vm.validate(*m1, true);
  for (const auto &msg : errout1) {
    TEST_ASSERT(
        msg ==
        "INFO: [AllowedAtomsValidation] Atom F is not in allowedAtoms list");
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testDisallowedAtomsValidation() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing DisallowedAtoms validation"
      << std::endl;

  //	std::vector<string> atoms = {"F", "Cl", "Br"};
  std::vector<unsigned int> atoms = {9, 17, 35};
  std::vector<shared_ptr<Atom>> atomList;

  for (auto &atom : atoms) {
    shared_ptr<Atom> a(new Atom(atom));
    atomList.push_back(a);
  }

  DisallowedAtomsValidation vm(atomList);
  std::string smi1;

  smi1 = "CC(=O)CF";
  unique_ptr<ROMol> m1(SmilesToMol(smi1));
  vector<ValidationErrorInfo> errout1 = vm.validate(*m1, true);
  for (const auto &msg : errout1) {
    TEST_ASSERT(
        msg ==
        "INFO: [DisallowedAtomsValidation] Atom F is in disallowedAtoms list");
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testFragment() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing fragment validation" << std::endl;

  string smi1, smi2, smi3, smi4, smi5, smi6;
  MolVSValidation vm;

  // testing MolVSValidation fragmentValidation
  // FragmentValidation should identify 1,2-dichloroethane.
  smi1 = "ClCCCl.c1ccccc1O";
  unique_ptr<ROMol> m1(SmilesToMol(smi1, 0, false));
  vector<ValidationErrorInfo> errout1 = vm.validate(*m1, true);
  for (const auto &msg : errout1) {
    TEST_ASSERT(msg ==
                "INFO: [FragmentValidation] 1,2-dichloroethane is present");
  }

  smi2 = "COCCOC.CCCBr";
  unique_ptr<ROMol> m2(SmilesToMol(smi2, 0, false));
  vector<ValidationErrorInfo> errout2 = vm.validate(*m2, true);
  for (const auto &msg : errout2) {
    TEST_ASSERT(msg ==
                "INFO: [FragmentValidation] 1,2-dimethoxyethane is present");
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testFeaturesValidation() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing FeaturesValidation"
                       << std::endl;

  unique_ptr<ROMol> mol;
  FeaturesValidation features;
  string mblock, errmsg;
  vector<ValidationErrorInfo> errout;

  mblock = R"(
                    2D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.8753 4.9367 0 0
M  V30 2 C -0.4583 4.1667 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  errout = features.validate(*mol, true);
  TEST_ASSERT(errout.empty());

  mblock = R"(
  Mrv2311 01162410492D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 A -20.708 9.9367 0 0
M  V30 2 C -22.0417 9.1667 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  errout = features.validate(*mol, true);
  TEST_ASSERT(errout.size() == 1);
  for (auto msg : errout) {
    cerr << msg << endl;
  }
  errmsg = errout[0];
  TEST_ASSERT(errmsg ==
              "ERROR: [FeaturesValidation] Query atom 0 is not allowed");
  {
    FeaturesValidation featuresCopy(features);
    featuresCopy.allowQueries = true;
    errout = featuresCopy.validate(*mol, true);
    TEST_ASSERT(errout.empty());
  }

  mblock = R"(
                    2D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.8753 4.9367 0 0
M  V30 2 C -0.4583 4.1667 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  mol->getAtomWithIdx(0)->setAtomicNum(0);
  errout = features.validate(*mol, true);
  TEST_ASSERT(errout.size() == 1);
  for (auto msg : errout) {
    cerr << msg << endl;
  }
  errmsg = errout[0];
  TEST_ASSERT(errmsg ==
              "ERROR: [FeaturesValidation] Dummy atom 0 is not allowed");
  {
    FeaturesValidation featuresCopy(features);
    featuresCopy.allowDummies = true;
    errout = featuresCopy.validate(*mol, true);
    TEST_ASSERT(errout.empty());
  }

  mblock = R"(
  Mrv2311 01162411522D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -20.708 9.9367 0 0
M  V30 2 C -22.0417 9.1667 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 8 2 1
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  errout = features.validate(*mol, true);
  TEST_ASSERT(errout.size() == 1);
  for (auto msg : errout) {
    cerr << msg << endl;
  }
  errmsg = errout[0];
  TEST_ASSERT(errmsg ==
              "ERROR: [FeaturesValidation] Query bond 0 is not allowed");
  {
    FeaturesValidation featuresCopy(features);
    featuresCopy.allowQueries = true;
    errout = featuresCopy.validate(*mol, true);
    TEST_ASSERT(errout.empty());
  }

  mblock = R"(
  Mrv2311 02272411562D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -10.3542 4.29 0 0
M  V30 2 C -11.6879 3.52 0 0
M  V30 3 C -11.6879 1.9798 0 0
M  V30 4 N -10.3542 1.21 0 0
M  V30 5 C -9.0204 1.9798 0 0
M  V30 6 C -9.0204 3.52 0 0
M  V30 7 C -10.3542 5.83 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 4 1 2
M  V30 2 4 1 6
M  V30 3 4 2 3
M  V30 4 4 5 6
M  V30 5 1 1 7
M  V30 6 4 3 4
M  V30 7 4 4 5
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  errout = features.validate(*mol, true);
  TEST_ASSERT(errout.size() == 6);
  for (auto msg : errout) {
    cerr << msg << endl;
  }
  errmsg = errout[0];
  TEST_ASSERT(
      errmsg ==
      "ERROR: [FeaturesValidation] Bond 0 of aromatic type is not allowed");
  {
    FeaturesValidation featuresCopy(features);
    featuresCopy.allowAromaticBondType = true;
    errout = featuresCopy.validate(*mol, true);
    TEST_ASSERT(errout.empty());
  }

  mblock = R"(
  Mrv2311 07222412542D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Pt -17.4792 5.75 0 0
M  V30 2 Cl -16.1042 6.8333 0 0
M  V30 3 Cl -16.1875 4.7917 0 0
M  V30 4 N -18.8958 6.8333 0 0
M  V30 5 N -18.8125 4.5833 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 3 1
M  V30 3 9 4 1
M  V30 4 9 5 1
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  errout = features.validate(*mol, true);
  TEST_ASSERT(errout.size() == 2);
  for (auto msg : errout) {
    cerr << msg << endl;
  }
  errmsg = errout[0];
  TEST_ASSERT(
      errmsg ==
      "ERROR: [FeaturesValidation] Bond 2 of dative type is not allowed");
  {
    FeaturesValidation featuresCopy(features);
    featuresCopy.allowDativeBondType = true;
    errout = featuresCopy.validate(*mol, true);
    TEST_ASSERT(errout.empty());
  }

  mblock = R"(
  MJ231601                      

  3  2  0  0  0  0  0  0  0  0999 V2000
   -8.7163    4.4303    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.4308    4.0178    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.7163    5.2553    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0  0  0  0
  1  3  1  0  0  0  0
A    3
CF3
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  errout = features.validate(*mol, true);
  TEST_ASSERT(errout.size() == 1);
  for (auto msg : errout) {
    cerr << msg << endl;
  }
  errmsg = errout[0];
  TEST_ASSERT(
      errmsg ==
      "ERROR: [FeaturesValidation] Atom 2 with alias 'CF3' is not allowed");
  {
    FeaturesValidation featuresCopy(features);
    featuresCopy.allowAtomAliases = true;
    errout = featuresCopy.validate(*mol, true);
    TEST_ASSERT(errout.empty());
  }

  mblock = R"(
  Mrv2311 01162411552D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -18.208 8.52 0 0 CFG=2
M  V30 2 F -19.5417 7.75 0 0
M  V30 3 C -16.8743 7.75 0 0
M  V30 4 Cl -18.208 10.06 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 3 CFG=1
M  V30 2 1 2 1
M  V30 3 1 1 4
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STERAC1 ATOMS=(1 1)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  errout = features.validate(*mol, true);
  TEST_ASSERT(errout.size() == 1);
  for (auto msg : errout) {
    cerr << msg << endl;
  }
  errmsg = errout[0];
  TEST_ASSERT(
      errmsg ==
      "ERROR: [FeaturesValidation] Enhanced stereochemistry features are not allowed");
  {
    FeaturesValidation featuresCopy(features);
    featuresCopy.allowEnhancedStereo = true;
    errout = featuresCopy.validate(*mol, true);
    TEST_ASSERT(errout.empty());
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testDisallowedRadicalValidation() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing DisallowedRadicalValidation"
      << std::endl;

  unique_ptr<ROMol> mol;
  DisallowedRadicalValidation radicalValidation;
  string mblock, errmsg;
  vector<ValidationErrorInfo> errout;

  mblock = R"(
  Mrv2311 02082417202D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -20.9372 7.145 0 0
M  V30 2 C -22.2708 6.375 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  errout = radicalValidation.validate(*mol, true);
  TEST_ASSERT(errout.empty());

  mblock = R"(
  Mrv2311 02082417212D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -20.9372 7.145 0 0 RAD=2
M  V30 2 C -22.2708 6.375 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  errout = radicalValidation.validate(*mol, true);
  TEST_ASSERT(errout.size() == 1);
  errmsg = errout[0];
  TEST_ASSERT(
      errmsg ==
      "ERROR: [DisallowedRadicalValidation] The radical at atom 0 is not allowed");

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testIs2DValidation() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Is2DValidation"
                       << std::endl;

  unique_ptr<ROMol> mol;
  Is2DValidation is2D;
  string mblock, errmsg;
  vector<ValidationErrorInfo> errout;

  mblock = R"(
                    2D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.8753 4.9367 0 0
M  V30 2 C -0.4583 4.1667 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  errout = is2D.validate(*mol, true);
  TEST_ASSERT(errout.empty());

  mblock = R"(
                    2D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.8753 4.9367 0 0
M  V30 2 C -0.4583 4.1667 0.2 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  errout = is2D.validate(*mol, true);
  TEST_ASSERT(errout.size() == 1);
  errmsg = errout[0];
  TEST_ASSERT(
      errmsg ==
      "ERROR: [Is2DValidation] The molecule includes non-null Z coordinates");

  mblock = R"(
                    2D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 2 C 0 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  errout = is2D.validate(*mol, true);
  TEST_ASSERT(errout.size() == 1);
  errmsg = errout[0];
  TEST_ASSERT(
      errmsg ==
      "ERROR: [Is2DValidation] All atoms have the same (x,y) coordinates");

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testLayout2DValidation() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Layout2DValidation"
                       << std::endl;

  unique_ptr<ROMol> mol;
  Layout2DValidation layout2D;
  string mblock, errmsg;
  vector<ValidationErrorInfo> errout;

  mblock = R"(
                    2D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.8753 4.9367 0 0
M  V30 2 C -0.4583 4.1667 0 0
M  V30 3 C -0.2691 5.9671 0 0
M  V30 4 C -1.7337 5.4912 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 1 3
M  V30 3 1 3 4
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  errout = layout2D.validate(*mol, true);
  TEST_ASSERT(errout.empty());

  mblock = R"(
                    2D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 5 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.6667 6.2067 0 0
M  V30 2 C -3.0004 5.4367 0 0
M  V30 3 C -3.0004 3.8965 0 0
M  V30 4 C -1.6667 3.1267 0 0
M  V30 5 C -0.3329 4.6000 0 0
M  V30 6 C -0.3329 4.7000 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 6
M  V30 3 1 2 3
M  V30 4 1 3 4
M  V30 5 1 4 5
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  errout = layout2D.validate(*mol, true);
  TEST_ASSERT(errout.size() == 1);
  errmsg = errout[0];
  TEST_ASSERT(errmsg ==
              "ERROR: [Layout2DValidation] Atom 4 is too close to atom 5");
  {
    Layout2DValidation layout2DCopy(layout2D);
    layout2DCopy.clashLimit = 1e-2;
    errout = layout2DCopy.validate(*mol, true);
    TEST_ASSERT(errout.empty());
  }

  mblock = R"(
                    2D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 5 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.3333 4.8733 0 0
M  V30 2 C -2.5837 4.7283 0 0
M  V30 3 C -2.7087 3.4798 0 0
M  V30 4 C -1.6667 3.1267 0 0
M  V30 5 C -0.3329 3.8965 0 0
M  V30 6 C -0.9913 3.495 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 6
M  V30 3 1 2 3
M  V30 4 1 3 4
M  V30 5 1 4 5
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  errout = layout2D.validate(*mol, true);
  TEST_ASSERT(errout.size() == 1);
  errmsg = errout[0];
  TEST_ASSERT(errmsg ==
              "ERROR: [Layout2DValidation] Atom 5 too close to bond 4");
  {
    Layout2DValidation layout2DCopy(layout2D);
    layout2DCopy.clashLimit = 1e-2;
    errout = layout2DCopy.validate(*mol, true);
    TEST_ASSERT(errout.empty());
  }

  mblock = R"(
          01112413352D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -28.1663 10.4367 0 0
M  V30 2 C -29.5 9.6667 0 0
M  V30 3 C -29.5 11.2067 0 0
M  V30 4 F 25.0 10.4367 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 2 3
M  V30 3 1 3 1
M  V30 4 1 1 4
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  errout = layout2D.validate(*mol, true);
  TEST_ASSERT(errout.size() == 1);
  errmsg = errout[0];
  TEST_ASSERT(
      errmsg ==
      "ERROR: [Layout2DValidation] The length of bond 3 between atoms 0 and 3 exceeds a configured limit");
  {
    Layout2DValidation layout2DCopy(layout2D);
    layout2DCopy.bondLengthLimit = 100.;
    errout = layout2DCopy.validate(*mol, true);
    TEST_ASSERT(errout.empty());
  }

  // Long bonds in rings
  mblock = R"(
  Mrv2311 02222409302D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 17 17 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -9.0205 2.1033 0 0
M  V30 2 C -10.3542 1.3333 0 0
M  V30 3 C -7.4805 2.1033 0 0
M  V30 4 C -5.9405 2.1033 0 0
M  V30 5 C -4.4005 2.1033 0 0
M  V30 6 C -2.8605 2.1033 0 0
M  V30 7 C -1.3205 2.1033 0 0
M  V30 8 C 0.2195 2.1033 0 0
M  V30 9 C 1.7595 2.1033 0 0
M  V30 10 C 3.2995 2.1033 0 0
M  V30 11 C 4.8395 2.1033 0 0
M  V30 12 C 6.3795 2.1033 0 0
M  V30 13 C 7.9195 2.1033 0 0
M  V30 14 C 9.4595 2.1033 0 0
M  V30 15 C 10.9995 2.1033 0 0
M  V30 16 C 12.5395 2.1033 0 0
M  V30 17 C 13.7854 1.1981 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 1 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 6 7
M  V30 7 1 7 8
M  V30 8 1 8 9
M  V30 9 1 9 10
M  V30 10 1 10 11
M  V30 11 1 11 12
M  V30 12 1 12 13
M  V30 13 1 13 14
M  V30 14 1 14 15
M  V30 15 1 15 16
M  V30 16 1 16 17
M  V30 17 1 17 2
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  errout = layout2D.validate(*mol, true);
  TEST_ASSERT(errout.size() == 0);

  Layout2DValidation customLayout2D(0.15, 10.0, false);
  errout = customLayout2D.validate(*mol, true);
  TEST_ASSERT(errout.size() == 1);
  errmsg = errout[0];
  TEST_ASSERT(
      errmsg ==
      "ERROR: [Layout2DValidation] The length of bond 16 between atoms 16 and 1 exceeds a configured limit");

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testValidateStereo() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing ValidateStereo"
                       << std::endl;

  unique_ptr<ROMol> mol;
  StereoValidation stereo;
  string mblock, errmsg;
  vector<ValidationErrorInfo> errout;

  // 4 ligands - no issues
  mblock = R"(
          10052309532D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Br 0.0003 7.27 0 0
M  V30 2 C -1.3333 6.5 0 0 CFG=1
M  V30 3 F -2.667 7.27 0 0
M  V30 4 O -1.3333 4.96 0 0
M  V30 5 C 0.0003 5.73 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 5 CFG=1
M  V30 2 1 2 3
M  V30 3 1 2 1
M  V30 4 1 2 4
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  errout = stereo.validate(*mol, true);
  TEST_ASSERT(errout.size() == 0);

  // 4 ligands - too many stereo bonds with the same wedge/dash direction
  mblock = R"(
          10052310002D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Br 0.0003 7.27 0 0
M  V30 2 C -1.3333 6.5 0 0 CFG=1
M  V30 3 F -2.667 7.27 0 0
M  V30 4 O -1.3333 4.96 0 0
M  V30 5 C 0.0003 5.73 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 5 CFG=1
M  V30 2 1 2 3 CFG=1
M  V30 3 1 2 1 CFG=1
M  V30 4 1 2 4
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  Chirality::reapplyMolBlockWedging(*mol);
  errout = stereo.validate(*mol, true);
  TEST_ASSERT(errout.size() == 2);
  errmsg = errout[0];
  TEST_ASSERT(
      errmsg ==
      "ERROR: [StereoValidation] Atom 1 has too many stereo bonds with like orientation");
  errmsg = errout[1];
  TEST_ASSERT(
      errmsg ==
      "ERROR: [StereoValidation] Atom 1 has adjacent stereo bonds with like orientation");

  // 4 ligands - mismatching opposed wedge/dash bonds
  mblock = R"(
          10052311582D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Br 0.0003 7.27 0 0
M  V30 2 C -1.3333 6.5 0 0
M  V30 3 F -2.667 7.27 0 0
M  V30 4 O -1.3333 4.96 0 0
M  V30 5 C 0.0003 5.73 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 5 CFG=1
M  V30 2 1 2 3 CFG=3
M  V30 3 1 2 1
M  V30 4 1 2 4
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  // 4 ligands - mismatching opposed wedge/dash bonds
  mol.reset(MolBlockToMol(mblock, false, false));
  Chirality::reapplyMolBlockWedging(*mol);
  errout = stereo.validate(*mol, true);
  TEST_ASSERT(errout.size() == 1);
  errmsg = errout[0];
  TEST_ASSERT(
      errmsg ==
      "ERROR: [StereoValidation] Atom 1 has opposing stereo bonds with different up/down orientation")

  // 4 ligands - potentially ambiguous umbrella configuration
  mblock = R"(
          10052313232D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Br 0.0003 7.27 0 0
M  V30 2 C -1.3333 6.5 0 0
M  V30 3 F -2.667 7.27 0 0
M  V30 4 O -1.3333 4.96 0 0
M  V30 5 C 0.0003 5.73 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 5
M  V30 2 1 2 3 CFG=3
M  V30 3 1 2 1
M  V30 4 1 2 4
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  Chirality::reapplyMolBlockWedging(*mol);
  errout = stereo.validate(*mol, true);
  TEST_ASSERT(errout.size() == 1);
  errmsg = errout[0];
  TEST_ASSERT(
      errmsg ==
      "ERROR: [StereoValidation] Atom 1 has a potentially ambiguous representation: all non-stereo bonds opposite to the only stereo bond")

  // 4 ligands - colinearity / triangle rule violation
  mblock = R"(
          10052313312D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 F -1.083 6.5617 0 0
M  V30 2 C -2.4167 5.7917 0 0
M  V30 3 O -3.7503 6.5617 0 0
M  V30 4 C -1.2083 5.0017 0 0
M  V30 5 Cl -2.4167 6.5617 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 4 CFG=1
M  V30 2 1 2 1
M  V30 3 1 2 5
M  V30 4 1 2 3
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  Chirality::reapplyMolBlockWedging(*mol);
  errout = stereo.validate(*mol, true);
  TEST_ASSERT(errout.size() == 1);
  errmsg = errout[0];
  TEST_ASSERT(
      errmsg ==
      "ERROR: [StereoValidation] Colinearity or triangle rule violation of non-stereo bonds at atom 1")

  // 4 ligands - wavy bond is allowed
  mblock = R"(
  Mrv2311 02022412452D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 12 13 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -14.6042 3.6233 0 0 CFG=3
M  V30 2 C -15.9379 2.8533 0 0
M  V30 3 C -15.9379 1.3131 0 0
M  V30 4 C -14.6042 0.5433 0 0
M  V30 5 C -13.2704 1.3131 0 0
M  V30 6 C -13.2704 2.8533 0 0 CFG=1
M  V30 7 C -11.9366 3.6235 0 0
M  V30 8 C -11.9369 5.1634 0 0
M  V30 9 C -13.2704 5.9335 0 0
M  V30 10 C -14.6042 5.1634 0 0
M  V30 11 C -14.6042 2.0833 0 0
M  V30 12 H -13.2704 4.3933 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 6
M  V30 3 1 2 3
M  V30 4 1 3 4
M  V30 5 1 4 5
M  V30 6 1 5 6
M  V30 7 1 7 8
M  V30 8 1 8 9
M  V30 9 1 9 10
M  V30 10 1 1 10
M  V30 11 1 6 7
M  V30 12 1 1 11 CFG=2
M  V30 13 1 6 12 CFG=1
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  Chirality::reapplyMolBlockWedging(*mol);
  errout = stereo.validate(*mol, true);
  TEST_ASSERT(errout.size() == 0);

  // 3 Ligands - No issues
  mblock = R"(
          10052313452D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Cl -0.9997 6.895 0 0
M  V30 2 C -2.3333 6.125 0 0 CFG=2
M  V30 3 F -3.667 6.895 0 0
M  V30 4 C -2.3333 4.585 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 4 CFG=1
M  V30 2 1 2 3
M  V30 3 1 2 1
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  Chirality::reapplyMolBlockWedging(*mol);
  errout = stereo.validate(*mol, true);
  TEST_ASSERT(errout.size() == 0);

  // 3 Ligands - multiple stereo bonds
  mblock = R"(
          10052313452D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Cl -0.9997 6.895 0 0
M  V30 2 C -2.3333 6.125 0 0 CFG=2
M  V30 3 F -3.667 6.895 0 0
M  V30 4 C -2.3333 4.585 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 4 CFG=1
M  V30 2 1 2 3 CFG=1
M  V30 3 1 2 1
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  Chirality::reapplyMolBlockWedging(*mol);
  errout = stereo.validate(*mol, true);
  TEST_ASSERT(errout.size() == 1);
  errmsg = errout[0];
  TEST_ASSERT(
      errmsg ==
      "ERROR: [StereoValidation] Atom 1 has 3 explicit substituents and multiple stereo bonds")

  // 3 Ligands - colinearity violation
  mblock = R"(
          10052313452D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Cl -0.9997 6.125 0 0
M  V30 2 C -2.3333 6.125 0 0 CFG=2
M  V30 3 F -3.667 6.125 0 0
M  V30 4 C -1.6667 4.71 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 4 CFG=1
M  V30 2 1 2 3
M  V30 3 1 2 1
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  Chirality::reapplyMolBlockWedging(*mol);
  errout = stereo.validate(*mol, true);
  TEST_ASSERT(errout.size() == 1);
  errmsg = errout[0];
  TEST_ASSERT(
      errmsg ==
      "ERROR: [StereoValidation] Colinearity of non-stereo bonds at atom 1");

  // 3 Ligands - either/unknown bond allowed
  mblock = R"(
  Mrv2311 02022410402D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -15.6875 8.3933 0 0
M  V30 2 N -16.9333 7.488 0 0
M  V30 3 C -16.4575 6.0234 0 0
M  V30 4 C -14.9175 6.0234 0 0 CFG=3
M  V30 5 C -14.4417 7.488 0 0
M  V30 6 C -14.0123 4.7775 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 5
M  V30 2 1 3 4
M  V30 3 1 4 5
M  V30 4 1 1 2
M  V30 5 1 2 3
M  V30 6 1 4 6 CFG=2
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  Chirality::reapplyMolBlockWedging(*mol);
  errout = stereo.validate(*mol, true);
  TEST_ASSERT(errout.size() == 0);

  // Double bond w/ 2 incident wavy bonds
  mblock = R"(
  Mrv2311 02022410492D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -15.6875 8.3933 0 0
M  V30 2 N -16.9333 7.488 0 0
M  V30 3 C -16.4575 6.0234 0 0
M  V30 4 C -14.9175 6.0234 0 0
M  V30 5 C -14.4417 7.488 0 0
M  V30 6 C -14.0123 4.7775 0 0
M  V30 7 C -14.6386 3.3706 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 5
M  V30 2 1 4 3 CFG=2
M  V30 3 1 4 5 CFG=2
M  V30 4 1 1 2
M  V30 5 1 2 3
M  V30 6 2 4 6
M  V30 7 1 6 7
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  Chirality::reapplyMolBlockWedging(*mol);
  errout = stereo.validate(*mol, true);
  for (auto msg : errout) {
    cerr << msg << endl;
  }
  TEST_ASSERT(errout.size() == 0);

  // Mixed wavy and wedged bonds
  mblock = R"(
  Mrv2311 02022412492D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -15.2292 8.6433 0 0
M  V30 2 N -16.475 7.738 0 0
M  V30 3 C -15.9992 6.2734 0 0
M  V30 4 C -14.4592 6.2734 0 0 CFG=3
M  V30 5 C -13.9834 7.738 0 0
M  V30 6 C -13.554 5.0275 0 0
M  V30 7 C -12.9381 6.5143 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 5
M  V30 2 1 3 4
M  V30 3 1 4 5
M  V30 4 1 1 2
M  V30 5 1 2 3
M  V30 6 1 4 6 CFG=1
M  V30 7 1 4 7 CFG=2
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  Chirality::reapplyMolBlockWedging(*mol);
  errout = stereo.validate(*mol, true);
  for (auto msg : errout) {
    cerr << msg << endl;
  }
  TEST_ASSERT(errout.size() == 1);
  errmsg = errout[0];
  TEST_ASSERT(
      errmsg ==
      "ERROR: [StereoValidation] Atom 3 has both unknown and wedged/dashed stereo bonds.");

  // Non-single bond with stereo bond orientation

  mblock = R"(
  Mrv2311 04232413302D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 S -11.583 11.3533 0 0
M  V30 2 C -12.9167 10.5833 0 0
M  V30 3 O -11.583 12.8933 0 0
M  V30 4 C -10.2493 10.5833 0 0
M  V30 5 C -10.2493 9.0433 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 4 5
M  V30 2 1 2 1
M  V30 3 1 1 4
M  V30 4 2 1 3 CFG=1
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  Chirality::reapplyMolBlockWedging(*mol);
  errout = stereo.validate(*mol, true);
  cerr << "here" << endl;
  for (auto msg : errout) {
    cerr << msg << endl;
  }
  cerr << "there" << endl;
  TEST_ASSERT(errout.size() == 1);
  errmsg = errout[0];
  TEST_ASSERT(
      errmsg ==
      "ERROR: [StereoValidation] Bond 3 has assigned stereo type, but unexpected bond order.");

  // Badly drawn perspective diagram
  mblock = R"(
  Mrv2311 02022413502D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 9 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -15.6872 7.1033 0 0
M  V30 2 C -16.7292 8.8333 0 0 CFG=2
M  V30 3 C -16.6861 7.8235 0 0
M  V30 4 C -15.1461 7.8235 0 0
M  V30 5 C -14.4805 7.27 0 0
M  V30 6 C -13.9002 8.7287 0 0 CFG=2
M  V30 7 C -13.5305 6.0579 0 0
M  V30 8 C -15.3955 9.6033 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1 CFG=1
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 1 5 CFG=1
M  V30 5 1 4 6
M  V30 6 1 6 5 CFG=1
M  V30 7 1 5 7
M  V30 8 1 2 8
M  V30 9 1 6 8
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  Chirality::reapplyMolBlockWedging(*mol);
  errout = stereo.validate(*mol, true);
  for (auto msg : errout) {
    cerr << msg << endl;
  }
  TEST_ASSERT(errout.size() == 1);
  errmsg = errout[0];
  TEST_ASSERT(
      errmsg ==
      "ERROR: [StereoValidation] Atom 0 has stereo bonds, but less than 3 explicit substituents.");

  // stereo bonds on allenes validate w/out errors
  mblock = R"(
  Mrv2311 02092408402D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -14.8958 5.5 0 0
M  V30 2 C -13.3558 5.5 0 0 CFG=2
M  V30 3 C -11.8158 5.5 0 0
M  V30 4 Cl -15.6658 6.8337 0 0
M  V30 5 Cl -11.0458 4.1663 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 2 2 3
M  V30 3 1 3 5 CFG=1
M  V30 4 1 1 4
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  Chirality::reapplyMolBlockWedging(*mol);
  errout = stereo.validate(*mol, true);
  TEST_ASSERT(errout.size() == 0);

  // stereo bonds in phenyl rings validate w/out errors
  mblock = R"(
  Mrv2311 02092408522D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -17.8543 8.7068 0 0
M  V30 2 C -19.1878 7.9368 0 0
M  V30 3 C -19.1878 6.3966 0 0
M  V30 4 C -17.8543 5.6266 0 0
M  V30 5 C -16.5205 6.3966 0 0
M  V30 6 C -16.5205 7.9368 0 0
M  V30 7 C -17.8543 4.0866 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5 CFG=1
M  V30 5 2 5 6
M  V30 6 1 6 1
M  V30 7 1 4 7
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  Chirality::reapplyMolBlockWedging(*mol);
  errout = stereo.validate(*mol, true);
  TEST_ASSERT(errout.size() == 0);

  // adjacent double bonds do not interfere with the validation of a
  // center's stereochemistry

  mblock = R"(
  Mrv2311 02142408162D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -6.9372 3.9783 0 0 CFG=2
M  V30 2 C -8.2708 3.2083 0 0
M  V30 3 F -6.9372 5.5183 0 0
M  V30 4 C -5.6035 3.2083 0 0
M  V30 5 H -6.9372 2.4383 0 0
M  V30 6 C -9.6045 3.9783 0 0
M  V30 7 C -10.9382 3.2083 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 4
M  V30 3 1 1 3
M  V30 4 1 1 5 CFG=1
M  V30 5 2 2 6
M  V30 6 1 6 7
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  Chirality::reapplyMolBlockWedging(*mol);
  errout = stereo.validate(*mol, true);

  for (auto msg : errout) {
    cerr << msg << endl;
  }

  TEST_ASSERT(errout.size() == 0);

  // stereo bonds between stereocenters do not interfere with the validation
  // of each individual center.
  // Note: according to IUPAC guidelines stereo bonds between stereocenters
  // should be avoided. The validation doesn't doesn't restrict the positioning
  // of stereo bonds on actual stereocenters and it would be complex to
  // determine if better options were available. The intention here is just to
  // try and make sure that in case any such stereo bonds are present, the
  // stereochemistry can be still interpreted correctly, even if their use may
  // be questionable.

  mblock = R"(
  Mrv2311 02142408012D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -6.9372 3.9783 0 0 CFG=2
M  V30 2 C -8.2708 3.2083 0 0
M  V30 3 F -6.9372 5.5183 0 0
M  V30 4 C -5.6035 3.2083 0 0
M  V30 5 H -6.9372 2.4383 0 0
M  V30 6 Br -9.6045 3.9783 0 0
M  V30 7 C -8.2708 1.6683 0 0
M  V30 8 F -8.2708 4.7483 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1 CFG=1
M  V30 2 1 1 4
M  V30 3 1 1 3
M  V30 4 1 2 7
M  V30 5 1 2 6
M  V30 6 1 1 5 CFG=1
M  V30 7 1 2 8
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  mol.reset(MolBlockToMol(mblock, false, false));
  Chirality::reapplyMolBlockWedging(*mol);
  errout = stereo.validate(*mol, true);

  for (auto msg : errout) {
    cerr << msg << endl;
  }

  TEST_ASSERT(errout.size() == 0);

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testValidateSmiles() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing ValidateSmiles"
                       << std::endl;

  // an invalid smiles should throw a ValueErrorException error
  try {
    vector<ValidationErrorInfo> errout1 = validateSmiles("3478q439g98h");
  } catch (const ValueErrorException &e) {
    std::string msg = e.what();
    TEST_ASSERT(msg ==
                "SMILES Parse Error: syntax error for input: 3478q439g98h")
  };

  vector<ValidationErrorInfo> errout2 = validateSmiles("");
  for (const auto &msg : errout2) {
    TEST_ASSERT(msg == "ERROR: [NoAtomValidation] Molecule has no atoms");
  }

  vector<ValidationErrorInfo> errout3 = validateSmiles("ClCCCl.c1ccccc1O");
  for (const auto &msg : errout3) {
    TEST_ASSERT(msg ==
                "INFO: [FragmentValidation] 1,2-dichloroethane is present");
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

int main() {
  testRDKitValidation();
  testMolVSValidation();
  testMolVSOptions();
  testAllowedAtomsValidation();
  testDisallowedAtomsValidation();
  testFragment();
  testFeaturesValidation();
  testDisallowedRadicalValidation();
  testIs2DValidation();
  testLayout2DValidation();
  testValidateStereo();
  testValidateSmiles();
  return 0;
}
