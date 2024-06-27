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
#include <iostream>
#include <fstream>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/StreamOps.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>

#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/Crippen.h>
#include <GraphMol/Descriptors/Property.h>

#include <GraphMol/PeriodicTable.h>
#include <GraphMol/atomic_data.h>
#include <GraphMol/Descriptors/USRDescriptor.h>

#include <DataStructs/BitVects.h>
#include <DataStructs/BitOps.h>
using namespace RDKit;
using namespace RDKit::Descriptors;

void test1() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Crippen parameter acquisition."
                        << std::endl;

  const CrippenParamCollection *params = CrippenParamCollection::getParams();
  TEST_ASSERT(params);

  CrippenParams p = *(params->begin());
  TEST_ASSERT(p.idx == 0);
  TEST_ASSERT(p.label == "C1");
  TEST_ASSERT(p.smarts == "[CH4]");

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void test2() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Crippen calculation." << std::endl;

  ROMol *mol;
  double logp, mr;

  mol = SmilesToMol("C");
  TEST_ASSERT(mol);
  calcCrippenDescriptors(*mol, logp, mr);
  TEST_ASSERT(feq(logp, 0.6361));
  TEST_ASSERT(feq(mr, 6.7310));
  // check singleton functions
  TEST_ASSERT(calcClogP(*mol) == logp);
  TEST_ASSERT(calcMR(*mol) == mr);
  // check that caching works:
  calcCrippenDescriptors(*mol, logp, mr);
  TEST_ASSERT(feq(logp, 0.6361));
  TEST_ASSERT(feq(mr, 6.7310));
  calcCrippenDescriptors(*mol, logp, mr, true, true);
  TEST_ASSERT(feq(logp, 0.6361));
  TEST_ASSERT(feq(mr, 6.7310));

  // check that things work when we don't add Hs:
  calcCrippenDescriptors(*mol, logp, mr, false, true);
  TEST_ASSERT(feq(logp, 0.1441));
  TEST_ASSERT(feq(mr, 2.503));
  delete mol;

  mol = SmilesToMol("C=C");
  TEST_ASSERT(mol);
  calcCrippenDescriptors(*mol, logp, mr, true);
  TEST_ASSERT(feq(logp, 0.8022));
  TEST_ASSERT(feq(mr, 11.2540));
  delete mol;

  mol = SmilesToMol("C#C");
  TEST_ASSERT(mol);
  calcCrippenDescriptors(*mol, logp, mr);
  TEST_ASSERT(feq(logp, 0.2494));
  TEST_ASSERT(feq(mr, 9.8900));
  delete mol;

  mol = SmilesToMol("CO");
  TEST_ASSERT(mol);
  calcCrippenDescriptors(*mol, logp, mr);
  TEST_ASSERT(feq(logp, -0.3915));
  TEST_ASSERT(feq(mr, 8.1428));
  delete mol;

  mol = SmilesToMol("C=O");
  TEST_ASSERT(mol);
  calcCrippenDescriptors(*mol, logp, mr);
  TEST_ASSERT(feq(logp, -0.1849));
  TEST_ASSERT(feq(mr, 7.121));
  delete mol;

  mol = SmilesToMol("C#[O+]");
  TEST_ASSERT(mol);
  calcCrippenDescriptors(*mol, logp, mr);
  TEST_ASSERT(feq(logp, 0.0059));
  TEST_ASSERT(feq(mr, 5.6315));
  delete mol;

  mol = SmilesToMol("C(C)(C)C");
  TEST_ASSERT(mol);
  calcCrippenDescriptors(*mol, logp, mr);
  TEST_ASSERT(feq(logp, 1.6623));
  TEST_ASSERT(feq(mr, 20.512));
  delete mol;

  mol = SmilesToMol("C(C)(C)(C)O");
  TEST_ASSERT(mol);
  calcCrippenDescriptors(*mol, logp, mr);
  TEST_ASSERT(feq(logp, 0.7772));
  TEST_ASSERT(feq(mr, 21.9718));
  delete mol;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testIssue262() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "    Test Issue262: problems with Crippen calculation from pickles."
      << std::endl;

  ROMol *mol, *mol2;
  RWMol *mol3;
  std::string pkl;
  double rlogp, rmr, logp, mr;

  mol = SmilesToMol("c1ncccc1");
  TEST_ASSERT(mol);
  calcCrippenDescriptors(*mol, rlogp, rmr);

  MolPickler::pickleMol(*mol, pkl);

  mol2 = new ROMol(pkl);
  TEST_ASSERT(mol2);
  calcCrippenDescriptors(*mol2, logp, mr);
  TEST_ASSERT(feq(logp, rlogp));
  TEST_ASSERT(feq(mr, rmr));

  mol3 = new RWMol();
  TEST_ASSERT(mol3);
  MolPickler::molFromPickle(pkl, mol3);

  calcCrippenDescriptors(*mol3, logp, mr);
  TEST_ASSERT(feq(logp, rlogp));
  TEST_ASSERT(feq(mr, rmr));

  delete mol;
  delete mol2;
  delete mol3;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void test3() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test AMW calculation." << std::endl;

  ROMol *mol, *mol2;
  double amw;

  mol = SmilesToMol("C");
  TEST_ASSERT(mol);
  amw = calcAMW(*mol);
  TEST_ASSERT(feq(amw, 16.043, .001));
  amw = calcAMW(*mol, true);
  TEST_ASSERT(feq(amw, 12.011, .001));
  mol2 = MolOps::addHs(*mol);
  amw = calcAMW(*mol2);
  TEST_ASSERT(feq(amw, 16.043, .001));
  amw = calcAMW(*mol2, true);
  TEST_ASSERT(feq(amw, 12.011, .001));
  delete mol;
  delete mol2;

  mol = SmilesToMol("[CH4]");
  TEST_ASSERT(mol);
  amw = calcAMW(*mol);
  TEST_ASSERT(feq(amw, 16.043, .001));
  amw = calcAMW(*mol, true);
  TEST_ASSERT(feq(amw, 12.011, .001));
  delete mol;

  mol = SmilesToMol("C[2H]");
  TEST_ASSERT(mol);
  amw = calcAMW(*mol);
  delete mol;

  TEST_ASSERT(feq(amw, 17.0, .1));

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void test3a() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Exact MW calculation." << std::endl;

  ROMol *mol, *mol2;
  double mw;

  mol = SmilesToMol("C");
  TEST_ASSERT(mol);
  mw = calcExactMW(*mol);
  TEST_ASSERT(feq(mw, 16.031, .001));
  mw = calcExactMW(*mol, true);
  TEST_ASSERT(feq(mw, 12.000, .001));
  mol2 = MolOps::addHs(*mol);
  mw = calcExactMW(*mol2);
  TEST_ASSERT(feq(mw, 16.031, .001));
  mw = calcExactMW(*mol2, true);
  TEST_ASSERT(feq(mw, 12.000, .001));
  delete mol;
  delete mol2;

  mol = SmilesToMol("[CH4]");
  TEST_ASSERT(mol);
  mw = calcExactMW(*mol);
  TEST_ASSERT(feq(mw, 16.031, .001));
  mw = calcExactMW(*mol, true);
  TEST_ASSERT(feq(mw, 12.000, .001));
  delete mol;

  mol = SmilesToMol("C[2H]");
  TEST_ASSERT(mol);
  mw = calcExactMW(*mol);
  TEST_ASSERT(feq(mw, 17.037, .001));
  mw = calcExactMW(*mol, true);
  TEST_ASSERT(feq(mw, 12.000, .001));
  delete mol;

  mol = SmilesToMol("Cl");
  TEST_ASSERT(mol);
  mw = calcAMW(*mol);
  TEST_ASSERT(feq(mw, 35.453 + 1.008, .001));
  mw = calcExactMW(*mol);
  TEST_ASSERT(feq(mw, 34.9688 + 1.0078, .001));
  delete mol;

  mol = SmilesToMol("[35ClH]");
  TEST_ASSERT(mol);
  mw = calcAMW(*mol);
  TEST_ASSERT(feq(mw, 34.9688 + 1.008, .001));
  mw = calcExactMW(*mol);
  TEST_ASSERT(feq(mw, 34.9688 + 1.0078, .001));
  delete mol;

  mol = SmilesToMol("[36ClH]");
  TEST_ASSERT(mol);
  mw = calcAMW(*mol);
  TEST_ASSERT(feq(mw, 35.9683 + 1.008, .001));
  mw = calcExactMW(*mol);
  TEST_ASSERT(feq(mw, 35.9683 + 1.0078, .001));
  delete mol;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testLabute() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Labute ASA descriptors." << std::endl;
  ROMol *mol;
  double asa;

  mol = SmilesToMol("CO");
  asa = calcLabuteASA(*mol);
  TEST_ASSERT(feq(asa, 13.5335, .0001));
  asa = calcLabuteASA(*mol);
  TEST_ASSERT(feq(asa, 13.5335, .0001));
  asa = calcLabuteASA(*mol, true, true);
  TEST_ASSERT(feq(asa, 13.5335, .0001));

  delete mol;
  mol = SmilesToMol("OC(=O)c1ccncc1C(=O)O");
  asa = calcLabuteASA(*mol);
  TEST_ASSERT(feq(asa, 67.2924, .0001));

  delete mol;
  mol = SmilesToMol("C1CCC(c2cccnc2)NC1");
  asa = calcLabuteASA(*mol);
  TEST_ASSERT(feq(asa, 73.0198, .0001));

  delete mol;
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testTPSA() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test TPSA descriptors." << std::endl;

  std::string fName = getenv("RDBASE");
  fName += "/Data/NCI/first_200.tpsa.csv";
  std::ifstream inf(fName.c_str());
  TEST_ASSERT(inf && !inf.bad());

  while (!inf.eof()) {
    std::string inl = getLine(inf);
    boost::trim(inl);
    if (inl.size() == 0 || inl[0] == '#') {
      continue;
    }
    std::vector<std::string> tokens;
    boost::split(tokens, inl, boost::is_any_of(","));
    if (tokens.size() != 2) {
      continue;
    }
    std::string smiles = tokens[0];
    auto oTPSA = boost::lexical_cast<double>(tokens[1]);
    ROMol *mol = SmilesToMol(smiles);
    TEST_ASSERT(mol);
    double nTPSA = calcTPSA(*mol);
    if (!feq(nTPSA, oTPSA, .0001)) {
      std::cerr << " TPSA ERR: " << smiles << " " << oTPSA << " " << nTPSA
                << std::endl;
      std::vector<double> contribs(mol->getNumAtoms());
      getTPSAAtomContribs(*mol, contribs);
      for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
        std::cerr << "\t" << i << "\t" << contribs[i] << std::endl;
      }
    }
    TEST_ASSERT(feq(nTPSA, oTPSA, .0001));

    // make sure that adding Hs doesn't affect the value
    // (this was issue 1969745)
    ROMol *mol2 = MolOps::addHs(*mol);
    double hTPSA = calcTPSA(*mol2);
    TEST_ASSERT(feq(nTPSA, hTPSA, .0001));

    delete mol2;
    delete mol;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

const std::string NonStrictRotProp = "NUM_ROTATABLEBONDS_O";
const std::string StrictRotProp = "NUM_ROTATABLEBONDS";

void testLipinski1() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Lipinski parameters." << std::endl;

  std::string fName = getenv("RDBASE");
  fName += "/Data/NCI/first_200.props.sdf";
  SDMolSupplier suppl(fName);
  int idx = -1;
  // Figure out which rotatable bond version we are using
  std::string rot_prop = StrictRotProp;
  {
    const bool sanitize = true;
    ROMol *test_mol =
        SmilesToMol("CC(C)(C)c1cc(O)c(cc1O)C(C)(C)C", 0, sanitize);
    if (calcNumRotatableBonds(*test_mol) == 2) {
      rot_prop = NonStrictRotProp;
    }
    delete test_mol;
  }
  while (!suppl.atEnd()) {
    ROMol *mol = nullptr;
    ++idx;
    try {
      mol = suppl.next();
    } catch (...) {
      continue;
    }
    if (!mol) {
      continue;
    }

    unsigned int oVal, nVal;
    std::string foo;

    mol->getProp("NUM_HACCEPTORS", foo);
    oVal = boost::lexical_cast<unsigned int>(foo);
    nVal = calcNumHBA(*mol);
    if (oVal != nVal) {
      std::cerr << "  failed: " << idx << " " << oVal << " " << nVal
                << std::endl;
    }
    TEST_ASSERT(oVal == nVal);

    mol->getProp("NUM_HDONORS", foo);
    oVal = boost::lexical_cast<unsigned int>(foo);
    nVal = calcNumHBD(*mol);
    if (oVal != nVal) {
      std::cerr << "  failed: " << idx << " " << oVal << " " << nVal
                << std::endl;
    }
    TEST_ASSERT(oVal == nVal);

    mol->getProp("NUM_LIPINSKIHDONORS", foo);
    oVal = boost::lexical_cast<unsigned int>(foo);
    nVal = calcLipinskiHBD(*mol);
    if (oVal != nVal) {
      std::cerr << "  failed: " << idx << " " << oVal << " " << nVal
                << std::endl;
    }
    TEST_ASSERT(oVal == nVal);

    mol->getProp("NUM_LIPINSKIHACCEPTORS", foo);
    oVal = boost::lexical_cast<unsigned int>(foo);
    nVal = calcLipinskiHBA(*mol);
    if (oVal != nVal) {
      std::cerr << "  failed: " << idx << " " << oVal << " " << nVal
                << std::endl;
    }
    TEST_ASSERT(oVal == nVal);

    mol->getProp("NUM_RINGS", foo);
    oVal = boost::lexical_cast<unsigned int>(foo);
    nVal = calcNumRings(*mol);
    if (oVal != nVal) {
      std::cerr << "  failed: " << idx << " " << oVal << " " << nVal
                << std::endl;
    }
    TEST_ASSERT(oVal == nVal);

    mol->getProp("NUM_HETEROATOMS", foo);
    oVal = boost::lexical_cast<unsigned int>(foo);
    nVal = calcNumHeteroatoms(*mol);
    if (oVal != nVal) {
      std::cerr << "  failed: " << idx << " " << oVal << " " << nVal
                << std::endl;
    }
    TEST_ASSERT(oVal == nVal);

    mol->getProp(rot_prop, foo);
    oVal = boost::lexical_cast<unsigned int>(foo);
    nVal = calcNumRotatableBonds(*mol);
    if (oVal != nVal) {
      std::cerr << "  failed: " << idx << " " << oVal << " " << nVal
                << " using stored sd prop " << rot_prop << std::endl;
    }
    TEST_ASSERT(oVal == nVal);

    mol->getProp(NonStrictRotProp, foo);
    oVal = boost::lexical_cast<unsigned int>(foo);
    nVal = calcNumRotatableBonds(*mol, NonStrict);
    if (oVal != nVal) {
      std::cerr << "  failed: " << idx << " " << oVal << " " << nVal
                << " using stored sd prop " << rot_prop << std::endl;
    }
    TEST_ASSERT(oVal == nVal);

    mol->getProp(StrictRotProp, foo);
    oVal = boost::lexical_cast<unsigned int>(foo);
    nVal = calcNumRotatableBonds(*mol, Strict);
    if (oVal != nVal) {
      std::cerr << "  failed: " << idx << " " << oVal << " " << nVal
                << " using stored sd prop " << rot_prop << std::endl;
    }
    TEST_ASSERT(oVal == nVal);

    delete mol;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testVSADescriptors() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test VSA descriptors." << std::endl;

  {
    ROMol *mol;
    std::vector<double> vals;

    mol = SmilesToMol("CO");
    vals = calcSlogP_VSA(*mol);
    TEST_ASSERT(vals.size() == 12);
    for (unsigned int i = 0; i < vals.size(); ++i) {
      switch (i) {
        case 1:
          TEST_ASSERT(feq(vals[i], 12.216, .001));
          break;
        default:
          TEST_ASSERT(feq(vals[i], 0, .001));
      }
    }
    delete mol;

    mol = SmilesToMol("CCO");
    vals = calcSlogP_VSA(*mol);
    TEST_ASSERT(vals.size() == 12);
    for (unsigned int i = 0; i < vals.size(); ++i) {
      switch (i) {
        case 1:
          TEST_ASSERT(feq(vals[i], 11.713, .001));
          break;
        case 4:
          TEST_ASSERT(feq(vals[i], 6.924, .001));
          break;
        default:
          TEST_ASSERT(feq(vals[i], 0, .001));
      }
    }
    delete mol;

    mol = SmilesToMol("Fc1ccccc1");
    vals = calcSlogP_VSA(*mol);
    TEST_ASSERT(vals.size() == 12);
    for (unsigned int i = 0; i < vals.size(); ++i) {
      switch (i) {
        case 3:
          TEST_ASSERT(feq(vals[i], 5.817, .001));
          break;
        case 5:
          TEST_ASSERT(feq(vals[i], 30.332, .001));
          break;
        case 9:
          TEST_ASSERT(feq(vals[i], 4.390, .001));
          break;
        default:
          TEST_ASSERT(feq(vals[i], 0, .001));
      }
    }
    delete mol;
  }

  {
    ROMol *mol;
    std::vector<double> vals;

    mol = SmilesToMol("CO");
    vals = calcSMR_VSA(*mol);
    TEST_ASSERT(vals.size() == 10);
    for (unsigned int i = 0; i < vals.size(); ++i) {
      switch (i) {
        case 0:
          TEST_ASSERT(feq(vals[i], 5.106, .001));
          break;
        case 5:
          TEST_ASSERT(feq(vals[i], 7.110, .001));
          break;
        default:
          TEST_ASSERT(feq(vals[i], 0, .001));
      }
    }
    delete mol;

    mol = SmilesToMol("CCO");
    vals = calcSMR_VSA(*mol);
    TEST_ASSERT(vals.size() == 10);
    for (unsigned int i = 0; i < vals.size(); ++i) {
      switch (i) {
        case 0:
          TEST_ASSERT(feq(vals[i], 5.106, .001));
          break;
        case 4:
          TEST_ASSERT(feq(vals[i], 6.924, .001));
          break;
        case 5:
          TEST_ASSERT(feq(vals[i], 6.607, .001));
          break;
        default:
          TEST_ASSERT(feq(vals[i], 0, .001));
      }
    }
    delete mol;
  }

  {
    ROMol *mol;
    std::vector<double> vals;

    mol = SmilesToMol("CO");
    vals = calcPEOE_VSA(*mol);
    TEST_ASSERT(vals.size() == 14);
    for (unsigned int i = 0; i < vals.size(); ++i) {
      switch (i) {
        case 0:
          TEST_ASSERT(feq(vals[i], 5.106, .001));
          break;
        case 7:
          TEST_ASSERT(feq(vals[i], 7.110, .001));
          break;
        default:
          TEST_ASSERT(feq(vals[i], 0, .001));
      }
    }
    delete mol;

    mol = SmilesToMol("CCO");
    vals = calcPEOE_VSA(*mol);
    TEST_ASSERT(vals.size() == 14);
    for (unsigned int i = 0; i < vals.size(); ++i) {
      switch (i) {
        case 0:
          TEST_ASSERT(feq(vals[i], 5.106, .001));
          break;
        case 6:
          TEST_ASSERT(feq(vals[i], 6.924, .001));
          break;
        case 7:
          TEST_ASSERT(feq(vals[i], 6.607, .001));
          break;
        default:
          TEST_ASSERT(feq(vals[i], 0, .001));
      }
    }
    delete mol;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMolFormula() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test molecular formula calculation."
                        << std::endl;

  ROMol *mol, *mol2;
  std::string formula;

  mol = SmilesToMol("C");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula == "CH4");
  mol2 = MolOps::addHs(*mol);
  formula = calcMolFormula(*mol2);
  TEST_ASSERT(formula == "CH4");
  delete mol;
  delete mol2;

  mol = SmilesToMol("[CH4]");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula == "CH4");
  delete mol;

  mol = SmilesToMol("CO");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula == "CH4O");
  mol2 = MolOps::addHs(*mol);
  formula = calcMolFormula(*mol2);
  TEST_ASSERT(formula == "CH4O");
  delete mol;
  delete mol2;

  mol = SmilesToMol("C(=O)N");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula == "CH3NO");
  mol2 = MolOps::addHs(*mol);
  formula = calcMolFormula(*mol2);
  TEST_ASSERT(formula == "CH3NO");
  delete mol;
  delete mol2;

  mol = SmilesToMol("C(=O)=O");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula == "CO2");
  delete mol;

  mol = SmilesToMol("C(=O)[O-]");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula == "CHO2-");
  delete mol;

  mol = SmilesToMol("C([O-])[O-]");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula == "CH2O2-2");
  delete mol;

  mol = SmilesToMol("C([NH3+])[O-]");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula == "CH5NO");
  delete mol;

  mol = SmilesToMol("C([NH3+])O");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula == "CH6NO+");
  delete mol;

  mol = SmilesToMol("C([NH3+])[NH3+]");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula == "CH8N2+2");
  delete mol;

  // H isotope tests
  mol = SmilesToMol("[2H]C([3H])O");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula == "CH4O");
  formula = calcMolFormula(*mol, true);
  TEST_ASSERT(formula == "CH2DTO");
  formula = calcMolFormula(*mol, true, false);
  TEST_ASSERT(formula == "CH2[2H][3H]O");
  delete mol;

  // isotope test
  mol = SmilesToMol("[13CH3]C([2H])O");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula == "C2H6O");
  formula = calcMolFormula(*mol, true);
  TEST_ASSERT(formula == "C[13C]H5DO");
  formula = calcMolFormula(*mol, true, false);
  TEST_ASSERT(formula == "C[13C]H5[2H]O");
  delete mol;

  // isotope test
  mol = SmilesToMol("[13CH3]C[13CH2]C");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula == "C4H10");
  formula = calcMolFormula(*mol, true);
  TEST_ASSERT(formula == "C2[13C]2H10");
  formula = calcMolFormula(*mol, true, false);
  TEST_ASSERT(formula == "C2[13C]2H10");
  delete mol;

  // order test
  mol = SmilesToMol("[13CH3]C[13CH2]CB(O)O[2H]");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula == "C4H11BO2");
  formula = calcMolFormula(*mol, true);
  TEST_ASSERT(formula == "C2[13C]2H10DBO2");
  formula = calcMolFormula(*mol, true, false);
  TEST_ASSERT(formula == "C2[13C]2H10[2H]BO2");
  delete mol;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testIssue3415534() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Issue 3415534." << std::endl;

  {
    ROMol *mol = SmilesToMol("CN");
    TEST_ASSERT(mol);
    int nHBD = calcLipinskiHBD(*mol);
    TEST_ASSERT(nHBD == 2);
    delete mol;
  }
  {
    ROMol *mol = SmilesToMol("CNC");
    TEST_ASSERT(mol);
    int nHBD = calcLipinskiHBD(*mol);
    TEST_ASSERT(nHBD == 1);
    delete mol;
  }
  {
    ROMol *mol = SmilesToMol("C[NH3+]");
    TEST_ASSERT(mol);
    int nHBD = calcLipinskiHBD(*mol);
    TEST_ASSERT(nHBD == 3);
    delete mol;
  }
  {
    ROMol *mol = SmilesToMol("CO");
    TEST_ASSERT(mol);
    int nHBD = calcLipinskiHBD(*mol);
    TEST_ASSERT(nHBD == 1);
    delete mol;
  }
  {
    ROMol *mol = SmilesToMol("C[OH2+]");
    TEST_ASSERT(mol);
    int nHBD = calcLipinskiHBD(*mol);
    TEST_ASSERT(nHBD == 2);
    delete mol;
  }
  {
    ROMol *mol = SmilesToMol("COC");
    TEST_ASSERT(mol);
    int nHBD = calcLipinskiHBD(*mol);
    TEST_ASSERT(nHBD == 0);
    delete mol;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testIssue3433771() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "    Test Issue3433771: Bad definition for Crippen atom type O11."
      << std::endl;

  ROMol *mol;
  double logp, mr;

  mol = SmilesToMol("O=C(NC)n1cccc1");
  TEST_ASSERT(mol);
  calcCrippenDescriptors(*mol, logp, mr);
  TEST_ASSERT(feq(logp, 0.6756, .001));

  delete mol;
  mol = SmilesToMol("O=C(n1cccc1)n1cccc1");
  TEST_ASSERT(mol);
  calcCrippenDescriptors(*mol, logp, mr);
  TEST_ASSERT(feq(logp, 1.806, .001));
  delete mol;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

#ifdef RDK_TEST_MULTITHREADED
namespace {
void runblock(const std::vector<ROMol *> &mols, unsigned int count,
              unsigned int idx) {
  for (unsigned int j = 0; j < 1000; j++) {
    for (unsigned int i = 0; i < mols.size(); ++i) {
      if (i % count != idx) {
        continue;
      }
      ROMol *mol = mols[i];
      int nHBD = calcNumHBD(*mol);
      int nHBA = calcNumHBA(*mol);

      int oVal;
      std::string foo;
      mol->getProp("NUM_HACCEPTORS", foo);
      oVal = boost::lexical_cast<int>(foo);
      TEST_ASSERT(oVal == nHBA);
      mol->getProp("NUM_HDONORS", foo);
      oVal = boost::lexical_cast<unsigned int>(foo);
      TEST_ASSERT(oVal == nHBD);

      unsigned int nAmide = calcNumAmideBonds(*mol);
      (void)nAmide;
      double logp, mr;
      calcCrippenDescriptors(*mol, logp, mr);
    }
  }
};
}  // namespace
#include <thread>
#include <future>
void testMultiThread() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test multithreading" << std::endl;

  std::string fName = getenv("RDBASE");
  fName += "/Data/NCI/first_200.props.sdf";
  SDMolSupplier suppl(fName);
  std::cerr << "reading molecules" << std::endl;
  std::vector<ROMol *> mols;
  while (!suppl.atEnd() && mols.size() < 100) {
    ROMol *mol = nullptr;
    try {
      mol = suppl.next();
    } catch (...) {
      continue;
    }
    if (!mol) {
      continue;
    }
    mols.push_back(mol);
  }
  std::vector<std::future<void>> tg;

  std::cerr << "processing" << std::endl;
  unsigned int count = 4;
  for (unsigned int i = 0; i < count; ++i) {
    std::cerr << " launch :" << i << std::endl;
    std::cerr.flush();
    tg.emplace_back(std::async(std::launch::async, runblock, mols, count, i));
  }
  for (auto &fut : tg) {
    fut.get();
  }
  for (auto &mol : mols) {
    delete mol;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}
#else
void testMultiThread() {}
#endif
void testCrippenContribs() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Crippen atom type calculations."
                        << std::endl;

  {
    ROMol *mol;
    mol = SmilesToMol("n1ccccc1CO");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());
    std::vector<unsigned int> ts(mol->getNumAtoms());
    std::vector<std::string> ls(mol->getNumAtoms());

    getCrippenAtomContribs(*mol, logp, mr, true, &ts, &ls);
    TEST_ASSERT(ts[0] == 59);
    TEST_ASSERT(ts[1] == 25);
    TEST_ASSERT(ts[2] == 25);
    TEST_ASSERT(ts[3] == 25);
    TEST_ASSERT(ts[4] == 25);
    TEST_ASSERT(ts[5] == 28);
    TEST_ASSERT(ts[6] == 17);
    TEST_ASSERT(ts[7] == 69);

    TEST_ASSERT(ls[0] == "N11");
    TEST_ASSERT(ls[7] == "O2");
    delete mol;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testIssue252() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "    Test Issue252: Bad definitions for Crippen atom types."
      << std::endl;

  {
    ROMol *mol;
    mol = SmilesToMol("O=[N+]([O-])C");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol, logp, mr, true);
    TEST_ASSERT(feq(logp[0], 0.0335, .001));
    TEST_ASSERT(feq(logp[1], -0.3396, .001));
    TEST_ASSERT(feq(logp[2], 0.0335, .001));
    TEST_ASSERT(feq(logp[3], -0.2035, .001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("CP");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol, logp, mr, true);
    TEST_ASSERT(feq(logp[0], -0.2035, .001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("C(C)P");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol, logp, mr, true);
    TEST_ASSERT(feq(logp[0], -0.2035, .001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("C(C)(C)P");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol, logp, mr, true);
    TEST_ASSERT(feq(logp[0], -0.2051, .001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("C(C)(C)(C)P");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol, logp, mr, true);
    TEST_ASSERT(feq(logp[0], -0.2051, .001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("C(=C)c1ccccc1");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol, logp, mr, true);
    TEST_ASSERT(feq(logp[0], 0.264, .001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("C(=C)(C)c1ccccc1");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol, logp, mr, true);
    TEST_ASSERT(feq(logp[0], 0.264, .001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("C=C");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol, logp, mr, true);
    TEST_ASSERT(feq(logp[0], 0.1551, .001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("O=S");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol, logp, mr, true);
    TEST_ASSERT(feq(logp[0], -0.3339, .001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("S=O");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol, logp, mr, true);
    TEST_ASSERT(feq(logp[0], -0.0024, .001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("C(Cl)C");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol, logp, mr, true);
    TEST_ASSERT(feq(logp[0], -0.2035, .001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("C(Cl)(Cl)C");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol, logp, mr, true);
    TEST_ASSERT(feq(logp[0], -0.2051, .001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("C(Cl)(Cl)(Cl)C");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol, logp, mr, true);
    TEST_ASSERT(feq(logp[0], -0.2051, .001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("C(Cl)c1ccccc1");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol, logp, mr, true);
    TEST_ASSERT(feq(logp[0], -0.0516, .001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("C(Cl)(Cl)c1ccccc1");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol, logp, mr, true);
    TEST_ASSERT(feq(logp[0], 0.1193, .001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("C(Cl)(Cl)(Cl)c1ccccc1");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol, logp, mr, true);
    TEST_ASSERT(feq(logp[0], -0.0967, .001));
    delete mol;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testChiVs() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test calculation of ChiVs." << std::endl;

  {
    std::string sdata[] = {"CCCCCC",     "CCC(C)CC", "CC(C)CCC", "CC(C)C(C)C",
                           "CC(C)(C)CC", "CCCCCO",   "CCC(O)CC", "CC(O)(C)CC",
                           "c1ccccc1O",  "CCCl",     "CCBr",     "CCI",
                           "EOS"};
    double ddata[] = {4.828, 4.992, 4.992, 5.155, 5.207, 4.276,
                      4.439, 4.654, 3.834, 2.841, 3.671, 4.242};
    unsigned int idx = 0;
    while (sdata[idx] != "EOS") {
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v = calcChi0v(*mol);
      TEST_ASSERT(feq(v, ddata[idx], 0.002));
      ++idx;
      delete mol;
    }
  }
  {
    std::string sdata[] = {"CCCCCC",     "CCC(C)CC", "CC(C)CCC", "CC(C)C(C)C",
                           "CC(C)(C)CC", "CCCCCO",   "CCC(O)CC", "CC(O)(C)CC",
                           "c1ccccc1O",  "EOS"};

    double ddata[] = {2.914, 2.808, 2.770, 2.643, 2.561,
                      2.523, 2.489, 2.284, 2.134};
    unsigned int idx = 0;
    while (sdata[idx] != "EOS") {
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v = calcChi1v(*mol);
      TEST_ASSERT(feq(v, ddata[idx], 0.002));
      ++idx;
      delete mol;
    }
  }
  {
    std::string sdata[] = {"CCCCCC",     "CCC(C)CC", "CC(C)CCC", "CC(C)C(C)C",
                           "CC(C)(C)CC", "CCCCCO",   "CCC(O)CC", "CC(O)(C)CC",
                           "c1ccccc1O",  "EOS"};

    double ddata[] = {1.707, 1.922, 2.183, 2.488, 2.914,
                      1.431, 1.470, 2.166, 1.336};

    unsigned int idx = 0;
    while (sdata[idx] != "EOS") {
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v = calcChi2v(*mol);
      TEST_ASSERT(feq(v, ddata[idx], 0.002));
      ++idx;
      delete mol;
    }
  }

  {
    std::string sdata[] = {"CCCCCC",     "CCC(C)CC", "CC(C)CCC", "CC(C)C(C)C",
                           "CC(C)(C)CC", "CCCCCO",   "CCC(O)CC", "CC(O)(C)CC",
                           "c1ccccc1O",  "EOS"};

    double ddata[] = {0.957, 1.394, 0.866, 1.333, 1.061,
                      0.762, 0.943, 0.865, 0.756};

    unsigned int idx = 0;
    while (sdata[idx] != "EOS") {
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v = calcChi3v(*mol);
      TEST_ASSERT(feq(v, ddata[idx], 0.002));
      ++idx;
      delete mol;
    }
  }

  {
    std::string sdata[] = {"CCCCCC",     "CCC(C)CC", "CC(C)CCC", "CC(C)C(C)C",
                           "CC(C)(C)CC", "CCCCCO",   "CCC(O)CC", "CC(O)(C)CC",
                           "c1ccccc1O",  "EOS"};

    double ddata[] = {0.500, 0.289, 0.577, 0.000, 0.000,
                      0.362, 0.289, 0.000, 0.428};

    unsigned int idx = 0;
    while (sdata[idx] != "EOS") {
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v = calcChi4v(*mol);
      TEST_ASSERT(feq(v, ddata[idx], 0.002));
      ++idx;
      delete mol;
    }
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testChiNs() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test calculation of ChiNs." << std::endl;

  {
    std::string sdata[] = {"CCCCCC",     "CCC(C)CC", "CC(C)CCC", "CC(C)C(C)C",
                           "CC(C)(C)CC", "CCCCCO",   "CCC(O)CC", "CC(O)(C)CC",
                           "c1ccccc1O",  "CCCl",     "CCBr",     "CCI",
                           "EOS"};
    double ddata[] = {4.828, 4.992, 4.992, 5.155, 5.207, 4.276,
                      4.439, 4.654, 3.834, 2.085, 2.085, 2.085};
    unsigned int idx = 0;
    while (sdata[idx] != "EOS") {
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v = calcChi0n(*mol);
      TEST_ASSERT(feq(v, ddata[idx], 0.002));
      ++idx;
      delete mol;
    }
  }
  {
    std::string sdata[] = {"CCCCCC",     "CCC(C)CC", "CC(C)CCC", "CC(C)C(C)C",
                           "CC(C)(C)CC", "CCCCCO",   "CCC(O)CC", "CC(O)(C)CC",
                           "c1ccccc1O",  "C=S",      "EOS"};

    double ddata[] = {2.914, 2.808, 2.770, 2.643, 2.561,
                      2.523, 2.489, 2.284, 2.134, 0.289};
    unsigned int idx = 0;
    while (sdata[idx] != "EOS") {
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v = calcChi1n(*mol);
      TEST_ASSERT(feq(v, ddata[idx], 0.002));
      ++idx;
      delete mol;
    }
  }
  {
    std::string sdata[] = {"CCCCCC",     "CCC(C)CC", "CC(C)CCC", "CC(C)C(C)C",
                           "CC(C)(C)CC", "CCCCCO",   "CCC(O)CC", "CC(O)(C)CC",
                           "c1ccccc1O",  "EOS"};

    double ddata[] = {1.707, 1.922, 2.183, 2.488, 2.914,
                      1.431, 1.470, 2.166, 1.336};

    unsigned int idx = 0;
    while (sdata[idx] != "EOS") {
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v = calcChi2n(*mol);
      TEST_ASSERT(feq(v, ddata[idx], 0.002));
      ++idx;
      delete mol;
    }
  }

  {
    std::string sdata[] = {"CCCCCC",     "CCC(C)CC", "CC(C)CCC", "CC(C)C(C)C",
                           "CC(C)(C)CC", "CCCCCO",   "CCC(O)CC", "CC(O)(C)CC",
                           "c1ccccc1O",  "EOS"};

    double ddata[] = {0.957, 1.394, 0.866, 1.333, 1.061,
                      0.762, 0.943, 0.865, 0.756};

    unsigned int idx = 0;
    while (sdata[idx] != "EOS") {
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v = calcChi3n(*mol);
      TEST_ASSERT(feq(v, ddata[idx], 0.002));
      ++idx;
      delete mol;
    }
  }

  {
    std::string sdata[] = {"CCCCCC",     "CCC(C)CC", "CC(C)CCC", "CC(C)C(C)C",
                           "CC(C)(C)CC", "CCCCCO",   "CCC(O)CC", "CC(O)(C)CC",
                           "c1ccccc1O",  "EOS"};

    double ddata[] = {0.500, 0.289, 0.577, 0.000, 0.000,
                      0.362, 0.289, 0.000, 0.428};

    unsigned int idx = 0;
    while (sdata[idx] != "EOS") {
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v = calcChi4n(*mol);
      TEST_ASSERT(feq(v, ddata[idx], 0.002));
      ++idx;
      delete mol;
    }
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testHallKierAlpha() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test calculation of HallKierAlpha."
                        << std::endl;

  {
    std::string sdata[] = {"C=O",
                           "CCC1(CC)C(=O)NC(=O)N(C)C1=O",
                           "OCC(O)C(O)C(O)C(O)CO",
                           "OCC1OC(O)C(O)C(O)C1O",
                           "Fc1c[nH]c(=O)[nH]c1=O",
                           "OC1CNC(C(=O)O)C1",
                           "CCCc1[nH]c(=S)[nH]c(=O)c1",
                           "CN(CCCl)CCCl",
                           "CBr",
                           "CI",
                           "EOS"};
    double ddata[] = {
        -0.3300, -1.3900, -0.2400, -0.2400, -1.3900,
        -0.6100, -0.9000, 0.5400,  0.480,   0.730,
    };
    unsigned int idx = 0;
    while (sdata[idx] != "EOS") {
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v = calcHallKierAlpha(*mol);
      TEST_ASSERT(feq(v, ddata[idx], 0.002));
      ++idx;
      delete mol;
    }
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testKappa1() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test calculation of Kappa1." << std::endl;

  {
    std::string sdata[] = {
        "C12CC2C3CC13",      "C1CCC12CC2",        "C1CCCCC1",          "CCCCCC",
        "CCC(C)C1CCC(C)CC1", "CC(C)CC1CCC(C)CC1", "CC(C)C1CCC(C)CCC1", "EOS"};
    double ddata[] = {2.344, 3.061, 4.167, 6.000, 9.091, 9.091, 9.091};
    unsigned int idx = 0;
    while (sdata[idx] != "EOS") {
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v = calcKappa1(*mol);
      TEST_ASSERT(feq(v, ddata[idx], 0.002));
      ++idx;
      delete mol;
    }
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testKappa2() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test calculation of Kappa2." << std::endl;

  {
    std::string sdata[] = {"[C+2](C)(C)(C)(C)(C)C",
                           "[C+](C)(C)(C)(C)(CC)",
                           "C(C)(C)(C)(CCC)",
                           "CC(C)CCCC",
                           "CCCCCCC",
                           "CCCCCC",
                           "CCCCCCC",
                           "C1CCCC1",
                           "C1CCCC1C",
                           "C1CCCCC1",
                           "C1CCCCCC1",
                           "CCCCC",
                           "CC=CCCC",
                           "C1=CN=CN1",
                           "c1ccccc1",
                           "c1cnccc1",
                           "n1ccncc1",
                           "CCCCF",
                           "CCCCCl",
                           "CCCCBr",
                           "CCC(C)C1CCC(C)CC1",
                           "CC(C)CC1CCC(C)CC1",
                           "CC(C)C1CCC(C)CCC1",
                           "EOS"};
    double ddata[] = {0.667000, 1.240000, 2.344400, 4.167000, 6.000000,
                      5.000000, 6.000000, 1.440000, 1.633000, 2.222000,
                      3.061000, 4.000000, 4.740000, 0.884000, 1.606000,
                      1.552000, 1.500000, 3.930000, 4.290000, 4.480000,
                      4.133000, 4.133000, 4.133000};
    unsigned int idx = 0;
    while (sdata[idx] != "EOS") {
      v2::SmilesParse::SmilesParserParams ps;
      ps.sanitize = false;
      auto mol = v2::SmilesParse::MolFromSmiles(sdata[idx], ps);
      TEST_ASSERT(mol);
      // we have some structures with unhappy valences, so be careful about the
      // sanitization:
      mol->updatePropertyCache(false);
      unsigned int opThatFailed;
      MolOps::sanitizeMol(*mol, opThatFailed,
                          MolOps::SANITIZE_ALL ^ MolOps::SANITIZE_PROPERTIES);
      double v = calcKappa2(*mol);
      TEST_ASSERT(feq(v, ddata[idx], 0.002));
      ++idx;
    }
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}
void testKappa3() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test calculation of Kappa3." << std::endl;

  {
    std::string sdata[] = {
        "C[C+](C)(C)(C)C(C)(C)C", "CCC(C)C(C)(C)(CC)", "CCC(C)CC(C)CC",
        "CC(C)CCC(C)CC",          "CC(C)CCCC(C)C",     "CCC(C)C1CCC(C)CC1",
        "CC(C)CC1CCC(C)CC1",      "CC(C)C1CCC(C)CCC1", "EOS"};
    double ddata[] = {2.000000, 2.380000, 4.500000, 5.878000,
                      8.000000, 2.500000, 3.265000, 2.844000};
    unsigned int idx = 0;
    while (sdata[idx] != "EOS") {
      v2::SmilesParse::SmilesParserParams ps;
      ps.sanitize = false;
      auto mol = v2::SmilesParse::MolFromSmiles(sdata[idx], ps);
      TEST_ASSERT(mol);
      // we have some structures with unhappy valences, so be careful about the
      // sanitization:
      mol->updatePropertyCache(false);
      unsigned int opThatFailed;
      MolOps::sanitizeMol(*mol, opThatFailed,
                          MolOps::SANITIZE_ALL ^ MolOps::SANITIZE_PROPERTIES);
      double v = calcKappa3(*mol);
      TEST_ASSERT(feq(v, ddata[idx], 0.002));
      ++idx;
    }
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testRingDescriptors() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test ring descriptors" << std::endl;

  std::string fName = getenv("RDBASE");
  fName += "/Code/GraphMol/Descriptors/test_data/aid466.trunc.sdf";
  SDMolSupplier suppl(fName);
  while (!suppl.atEnd()) {
    ROMol *mol = suppl.next();
    TEST_ASSERT(mol);
    unsigned int iv;
    mol->getProp("NumRings", iv);
    TEST_ASSERT(iv == calcNumRings(*mol));
    mol->getProp("NumAromaticRings", iv);
    TEST_ASSERT(iv == calcNumAromaticRings(*mol));
    mol->getProp("NumSaturatedRings", iv);
    TEST_ASSERT(iv == calcNumSaturatedRings(*mol));
    mol->getProp("NumAromaticHeterocycles", iv);
    TEST_ASSERT(iv == calcNumAromaticHeterocycles(*mol));
    mol->getProp("NumAromaticCarbocycles", iv);
    TEST_ASSERT(iv == calcNumAromaticCarbocycles(*mol));
    mol->getProp("NumSaturatedHeterocycles", iv);
    TEST_ASSERT(iv == calcNumSaturatedHeterocycles(*mol));
    mol->getProp("NumSaturatedCarbocycles", iv);
    TEST_ASSERT(iv == calcNumSaturatedCarbocycles(*mol));
    mol->getProp("NumAliphaticRings", iv);
    TEST_ASSERT(iv == calcNumAliphaticRings(*mol));
    mol->getProp("NumAliphaticHeterocycles", iv);
    TEST_ASSERT(iv == calcNumAliphaticHeterocycles(*mol));
    mol->getProp("NumAliphaticCarbocycles", iv);
    TEST_ASSERT(iv == calcNumAliphaticCarbocycles(*mol));
    mol->getProp("NumHeterocycles", iv);
    TEST_ASSERT(iv == calcNumHeterocycles(*mol));

    delete mol;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMiscCountDescriptors() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test other count descriptors." << std::endl;

  {
    ROMol *mol;
    mol = SmilesToMol("OCCO");
    TEST_ASSERT(feq(calcFractionCSP3(*mol), 1.0, 0.001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("OO");
    TEST_ASSERT(feq(calcFractionCSP3(*mol), 0.0, 0.001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("OC=CO");
    TEST_ASSERT(feq(calcFractionCSP3(*mol), 0.0, 0.001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("CCC=C");
    TEST_ASSERT(feq(calcFractionCSP3(*mol), 0.5, 0.001));
    delete mol;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMQNs() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test MQN" << std::endl;

  {
    unsigned int tgt[42] = {98, 0,  4,  0, 0,  1,  0,  3,  9, 5, 4, 124, 29, 3,
                            0,  66, 35, 0, 25, 30, 21, 2,  2, 0, 0, 6,   12, 6,
                            0,  70, 26, 0, 0,  0,  2,  16, 0, 0, 0, 0,   10, 5};
    // Figure out which rotatable bond version we are using
    //  and update the test accordingly
    {
      const bool sanitize = true;
      ROMol *test_mol =
          SmilesToMol("CC(C)(C)c1cc(O)c(cc1O)C(C)(C)C", 0, sanitize);
      if (calcNumRotatableBonds(*test_mol) == 2) {
        tgt[18] = 26;
      }
      delete test_mol;
    }
    std::vector<unsigned int> accum(42, 0);

    std::string fName = getenv("RDBASE");
    fName += "/Code/GraphMol/Descriptors/test_data/aid466.trunc.sdf";
    SDMolSupplier suppl(fName);
    while (!suppl.atEnd()) {
      ROMol *mol = suppl.next();
      TEST_ASSERT(mol);
      std::vector<unsigned int> v = calcMQNs(*mol);
      TEST_ASSERT(v.size() == 42);
      for (unsigned int i = 0; i < 42; ++i) {
        accum[i] += v[i];
      }
      delete mol;
    }
    for (unsigned int i = 0; i < 42; ++i) {
      if (accum[i] != tgt[i]) {
        std::cerr << " !! " << i << " " << accum[i] << "!=" << tgt[i]
                  << std::endl;
      }
      TEST_ASSERT(accum[i] == tgt[i]);
    }
  }
  {  // github #623
    ROMol *mol;
    mol = SmilesToMol("CC*");
    TEST_ASSERT(mol);
    std::vector<unsigned int> v = calcMQNs(*mol);
    TEST_ASSERT(v[11] == 2);
    delete mol;
  }
  {  // github #623
    ROMol *mol;
    mol = SmilesToMol("[2H][2H]");
    TEST_ASSERT(mol);
    std::vector<unsigned int> v = calcMQNs(*mol);
    TEST_ASSERT(v[11] == 0);
    delete mol;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testGitHubIssue56() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test GitHub Issue 56." << std::endl;

  {
    ROMol *mol;
    mol = SmilesToMol("[H+]");
    double mw;
    mw = calcAMW(*mol);
    TEST_ASSERT(feq(mw, 1.008, 0.001));
    mw = calcExactMW(*mol);
    TEST_ASSERT(feq(mw, 1.0078, 0.001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("[2H+]");
    double mw;
    mw = calcAMW(*mol);
    TEST_ASSERT(feq(mw, 2.014, 0.001));
    mw = calcExactMW(*mol);
    TEST_ASSERT(feq(mw, 2.014, 0.001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("[H]");
    double mw;
    mw = calcAMW(*mol);
    TEST_ASSERT(feq(mw, 1.008, 0.001));
    mw = calcExactMW(*mol);
    TEST_ASSERT(feq(mw, 1.0078, 0.001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("[2H]");
    double mw;
    mw = calcAMW(*mol);
    TEST_ASSERT(feq(mw, 2.014, 0.001));
    mw = calcExactMW(*mol);
    TEST_ASSERT(feq(mw, 2.014, 0.001));
    delete mol;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testGitHubIssue92() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "    Test Github92: Bad Crippen atom type for pyrrole H." << std::endl;

  {
    RWMol *mol;
    mol = SmilesToMol("c1cccn1[H]", 0, 0);
    TEST_ASSERT(mol);
    MolOps::sanitizeMol(*mol);
    TEST_ASSERT(mol->getNumAtoms() == 6);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());
    getCrippenAtomContribs(*mol, logp, mr, true);
    TEST_ASSERT(feq(logp[5], 0.2142, .001));
    delete mol;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testGitHubIssue463() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "    Test Github362: order dependence  in Kier-Hall descriptors."
      << std::endl;

  {  // start with the hall-kier delta values:
    RWMol *mol;
    mol = SmilesToMol("O=C(Nc1nccs1)NC(C1CC1)C");
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 14);
    unsigned int order[] = {0, 11, 8, 7, 2, 4, 5, 13, 10, 12, 9, 3, 1, 6};
    std::vector<unsigned int> nVect(
        order, order + sizeof(order) / sizeof(unsigned int));
    ROMol *nm = MolOps::renumberAtoms(*mol, nVect);
    TEST_ASSERT(nm);

    std::vector<double> hkds(mol->getNumAtoms());
    Descriptors::detail::hkDeltas(*mol, hkds, true);
    std::vector<double> nhkds(mol->getNumAtoms());
    Descriptors::detail::hkDeltas(*nm, nhkds, true);

    for (unsigned int j = 0; j < mol->getNumAtoms(); ++j) {
      TEST_ASSERT(feq(hkds[nVect[j]], nhkds[j]));
    }

    delete mol;
    delete nm;
  }

  {  // now chiNv values, where the problem was reported:

    RWMol *mol;
    mol = SmilesToMol("O=C(Nc1nccs1)NC(C1CC1)C");
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 14);
    unsigned int order[] = {0, 11, 8, 7, 2, 4, 5, 13, 10, 12, 9, 3, 1, 6};
    std::vector<unsigned int> nVect(
        order, order + sizeof(order) / sizeof(unsigned int));
    ROMol *nm = MolOps::renumberAtoms(*mol, nVect);
    TEST_ASSERT(nm);

    double cv = calcChi3v(*mol);
    double ncv = calcChi3v(*nm);

    TEST_ASSERT(feq(cv, ncv));

    delete mol;
    delete nm;
  }

  {  // now chiNn values, where the problem was reported:

    RWMol *mol;
    mol = SmilesToMol("O=C(Nc1nccs1)NC(C1CC1)C");
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 14);
    unsigned int order[] = {0, 11, 8, 7, 2, 4, 5, 13, 10, 12, 9, 3, 1, 6};
    std::vector<unsigned int> nVect(
        order, order + sizeof(order) / sizeof(unsigned int));
    ROMol *nm = MolOps::renumberAtoms(*mol, nVect);
    TEST_ASSERT(nm);

    double cv = calcChi3n(*mol);
    double ncv = calcChi3n(*nm);

    TEST_ASSERT(feq(cv, ncv));

    delete mol;
    delete nm;
  }

  {  // the root cause was handling of rings
    RWMol *mol;
    mol = SmilesToMol("C1CC1");
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 3);
    std::vector<double> hkds(mol->getNumAtoms());
    Descriptors::detail::hkDeltas(*mol, hkds, true);

    double cv = calcChi3v(*mol);
    TEST_ASSERT(feq(hkds[0], hkds[1]));
    TEST_ASSERT(feq(hkds[1], hkds[2]));
    TEST_ASSERT(feq(cv, hkds[0] * hkds[0] * hkds[0]));

    delete mol;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testSpiroAndBridgeheads() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "    Test calculation of spiro and bridgehead counts." << std::endl;

  {
    RWMol *mol;
    mol = SmilesToMol("C1CC2CCC1CC2");
    TEST_ASSERT(mol);

    unsigned int nSpiro = Descriptors::calcNumSpiroAtoms(*mol);
    TEST_ASSERT(nSpiro == 0);

    unsigned int nBridgehead = Descriptors::calcNumBridgeheadAtoms(*mol);
    TEST_ASSERT(nBridgehead == 2);

    delete mol;
  }

  {
    RWMol *mol;
    mol = SmilesToMol("C1CCC2(C1)CC1CCC2CC1");
    TEST_ASSERT(mol);

    unsigned int nSpiro = Descriptors::calcNumSpiroAtoms(*mol);
    TEST_ASSERT(nSpiro == 1);

    unsigned int nBridgehead = Descriptors::calcNumBridgeheadAtoms(*mol);
    TEST_ASSERT(nBridgehead == 2);

    delete mol;
  }

  {
    RWMol *mol;
    mol = SmilesToMol("CC1(C)CC2(C)CCC1(C)CC2");
    TEST_ASSERT(mol);

    unsigned int nSpiro = Descriptors::calcNumSpiroAtoms(*mol);
    TEST_ASSERT(nSpiro == 0);

    unsigned int nBridgehead = Descriptors::calcNumBridgeheadAtoms(*mol);
    TEST_ASSERT(nBridgehead == 2);

    delete mol;
  }

  {  // test the atoms parameter
    RWMol *mol;
    mol = SmilesToMol("C1CCC2(C1)CC1CCC2CC1");
    TEST_ASSERT(mol);

    std::vector<unsigned int> atoms;

    unsigned int nSpiro = Descriptors::calcNumSpiroAtoms(*mol, &atoms);
    TEST_ASSERT(nSpiro == 1);
    TEST_ASSERT(atoms.size() == nSpiro);
    TEST_ASSERT(std::find(atoms.begin(), atoms.end(), 3) != atoms.end());

    atoms.clear();
    unsigned int nBridgehead =
        Descriptors::calcNumBridgeheadAtoms(*mol, &atoms);
    TEST_ASSERT(nBridgehead == 2);
    TEST_ASSERT(atoms.size() == nBridgehead);
    TEST_ASSERT(std::find(atoms.begin(), atoms.end(), 9) != atoms.end());
    TEST_ASSERT(std::find(atoms.begin(), atoms.end(), 6) != atoms.end());

    delete mol;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testGitHubIssue694() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "    Test Github694: ExactMolWt ignoring the mass of the electron"
      << std::endl;

  {
    ROMol *mol = SmilesToMol("[35Cl]");
    TEST_ASSERT(mol);
    double mw = calcExactMW(*mol);
    TEST_ASSERT(
        feq(mw, PeriodicTable::getTable()->getMassForIsotope(17, 35), .000001));
    delete mol;
    mol = SmilesToMol("[35Cl-]");
    TEST_ASSERT(mol);
    mw = calcExactMW(*mol);
    TEST_ASSERT(feq(mw,
                    PeriodicTable::getTable()->getMassForIsotope(17, 35) +
                        constants::electronMass,
                    .000001));
    delete mol;
    mol = SmilesToMol("[35Cl+]");
    TEST_ASSERT(mol);
    mw = calcExactMW(*mol);
    TEST_ASSERT(feq(mw,
                    PeriodicTable::getTable()->getMassForIsotope(17, 35) -
                        constants::electronMass,
                    .000001));
    delete mol;
    mol = SmilesToMol("[35Cl+2]");
    TEST_ASSERT(mol);
    mw = calcExactMW(*mol);
    TEST_ASSERT(feq(mw,
                    PeriodicTable::getTable()->getMassForIsotope(17, 35) -
                        2 * constants::electronMass,
                    .000001));
    delete mol;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testProperties() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing Properties and KitchenSink"
                        << std::endl;

  {
    std::vector<std::string> all_names = Properties::getAvailableProperties();

    Properties sink;
    std::vector<std::string> names = sink.getPropertyNames();

    TEST_ASSERT(names == all_names);

    RWMol *mol;
    mol = SmilesToMol("C1CCC2(C1)CC1CCC2CC1");
    std::vector<double> props = sink.computeProperties(*mol);

    std::vector<double> res;
    for (const auto &prop : all_names) {
      std::vector<std::string> props;
      props.push_back(prop);
      Properties property(props);
      std::vector<double> computedProps = property.computeProperties(*mol);
      TEST_ASSERT(computedProps.size() == 1);
      res.push_back(computedProps[0]);
    }

    TEST_ASSERT(props == res);
    delete mol;
  }
  {
    std::vector<std::string> names;
    names.emplace_back("NumSpiroAtoms");
    names.emplace_back("NumBridgeheadAtoms");
    Properties sink(names);
    std::vector<std::string> sink_names = sink.getPropertyNames();
    TEST_ASSERT(names == sink_names);
    std::vector<double> res;
    res.push_back(1.);
    res.push_back(2.);

    RWMol *mol;
    mol = SmilesToMol("C1CCC2(C1)CC1CCC2CC1");
    TEST_ASSERT(mol);
    // Test annotation as well
    std::vector<double> props = sink.computeProperties(*mol, true);

    TEST_ASSERT(props == res);
    TEST_ASSERT(mol->getProp<double>("NumSpiroAtoms") == 1.);
    TEST_ASSERT(mol->getProp<double>("NumBridgeheadAtoms") == 2.);
    delete mol;

    mol = SmilesToMol("C1CCC2(C1)CC1CCC2CC1");
    TEST_ASSERT(mol);
    // Test annotation as well
    sink.annotateProperties(*mol);
    TEST_ASSERT(mol->getProp<double>("NumSpiroAtoms") == 1.);
    TEST_ASSERT(mol->getProp<double>("NumBridgeheadAtoms") == 2.);
    delete mol;
  }

  {
    try {
      std::vector<std::string> names;
      names.emplace_back("FakeName");
      Properties sink(names);
      TEST_ASSERT(0);  // should throw
    } catch (KeyErrorException &) {
      BOOST_LOG(rdErrorLog)
          << "---Caught keyerror (bad property name)---" << std::endl;
    }
  }
}

void testPropertyQueries() {
  RWMol *mol;
  mol = SmilesToMol("C1CCC2(C1)CC1CCC2CC1");
  {
    PROP_RANGE_QUERY *query = makePropertyRangeQuery("exactmw", 50, 300);
    TEST_ASSERT(query->Match(*mol));
    delete query;
  }
  {
    PROP_RANGE_QUERY *query = makePropertyRangeQuery("exactmw", 1000, 10300);
    TEST_ASSERT(!query->Match(*mol));
    delete query;
  }

  {
    auto pq = std::unique_ptr<PROP_EQUALS_QUERY>(
        makePropertyQuery<PROP_EQUALS_QUERY>("exactmw", calcExactMW(*mol)));
    TEST_ASSERT(pq->Match(*mol));

    pq = std::unique_ptr<PROP_EQUALS_QUERY>(
        makePropertyQuery<PROP_EQUALS_QUERY>("NumHBA", calcNumHBA(*mol)));
    TEST_ASSERT(pq->Match(*mol));

    pq =
        std::unique_ptr<PROP_EQUALS_QUERY>(makePropertyQuery<PROP_EQUALS_QUERY>(
            "lipinskiHBA", calcLipinskiHBA(*mol)));
    TEST_ASSERT(pq->Match(*mol));
  }
  delete mol;
}

void testStereoCounting() {
  const bool debugParse = false;
  const bool sanitize = false;
  ROMol *m = SmilesToMol("NC(C)(F)C(=O)O", debugParse, sanitize);
  TEST_ASSERT(!m->hasProp(common_properties::_StereochemDone));

  std::vector<std::string> names;
  names.emplace_back("NumAtomStereoCenters");
  names.emplace_back("NumUnspecifiedAtomStereoCenters");
  Properties prop(names);

  try {
    prop.computeProperties(*m);
    TEST_ASSERT(0);  // didn't catch exception
  } catch (ValueErrorException &) {
    BOOST_LOG(rdErrorLog) << "---Caught stereo value error---" << std::endl;
  }

  delete m;
  m = SmilesToMol("NC(C)(F)C(=O)O", sanitize);
  std::vector<double> res = prop.computeProperties(*m);
  TEST_ASSERT(res[0] == 1);
  TEST_ASSERT(res[1] == 1);
  delete m;
}

void testUSRDescriptor() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test USR Descriptor" << std::endl;
  std::vector<double> descriptor(12);

  // no conformers
  ROMol *mol = SmilesToMol("C1CCCCC1");
  bool ok = false;
  try {
    USR(*mol, descriptor);
  } catch (ConformerException &) {
    ok = true;
  }
  delete mol;
  TEST_ASSERT(ok);

  // number of atoms < 3
  mol = SmilesToMol("CC");
  ok = false;
  try {
    USR(*mol, descriptor);
  } catch (ValueErrorException &) {
    ok = true;
  }
  delete mol;
  TEST_ASSERT(ok);

  // DESCRIPTOR
  // comparing to results produced by Adrian Schreyer's code
  // http://hg.adrianschreyer.eu/usrcat/src/70e075d93cd25370e7ef93301d0e28d49a0851c2/usrcat/geometry.py?at=default
  double refValues[12] = {2.37938524,  0.62181927, -0.89089872, 2.63773456,
                          1.1577952,   -0.6937349, 3.38248245,  1.59816952,
                          -0.72933115, 3.38248245, 1.59816952,  -0.72933115};
  std::vector<double> refUSR(refValues,
                             refValues + sizeof(refValues) / sizeof(double));
  std::string rdbase = getenv("RDBASE");
  std::string fname1 =
      rdbase + "/Code/GraphMol/Descriptors/test_data/cyclohexane.mol";
  mol = MolFileToMol(fname1, true, false, true);
  std::vector<double> myUSR(12);
  USR(*mol, myUSR);
  for (unsigned int i = 0; i < myUSR.size(); ++i) {
    TEST_ASSERT(feq(myUSR[i], refUSR[i]));
  }
  delete mol;
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testUSRScore() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test USR Score" << std::endl;
  // SCORE
  // descriptors and reference score from JCC (2007), 28, 1711-1723.
  double m1[12] = {4.44, 2.98,  1.04,   4.55, 4.70,  0.23,
                   8.30, 16.69, -22.97, 7.37, 15.64, 0.51};
  double m2[12] = {4.39, 3.11,  1.36,   4.50, 4.44,  0.09,
                   8.34, 16.78, -23.20, 7.15, 16.52, 0.13};
  std::vector<double> d1(m1, m1 + sizeof(m1) / sizeof(double));
  std::vector<double> d2(m2, m2 + sizeof(m2) / sizeof(double));
  std::vector<double> weights(1, 1.0);
  TEST_ASSERT(feq(calcUSRScore(d1, d2, weights), 0.812, 0.001));

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testUSRCATDescriptor() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test USRCAT Descriptor" << std::endl;
  std::vector<double> descriptor(12);

  // no conformers
  ROMol *mol = SmilesToMol("C1CCCCC1");
  bool ok = false;
  try {
    USR(*mol, descriptor);
  } catch (ConformerException &) {
    ok = true;
  }
  TEST_ASSERT(ok);

  // number of atoms < 3
  mol = SmilesToMol("CC");
  ok = false;
  try {
    USR(*mol, descriptor);
  } catch (ValueErrorException &) {
    ok = true;
  }
  TEST_ASSERT(ok);

  // DESCRIPTOR
  // comparing to results produced by Adrian Schreyer's code
  // http://hg.adrianschreyer.eu/usrcat/src/70e075d93cd2?at=default
  double refValues[60] = {
      2.37938524, 0.62181927, -0.89089872, 2.63773456, 1.1577952,  -0.6937349,
      3.38248245, 1.59816952, -0.72933115, 3.38248245, 1.59816952, -0.72933115,
      1.50000000, 0.00000000, -0.77827171, 1.86602540, 1.00893468, -0.89059233,
      3.02358914, 1.02718486, -0.64081261, 3.02358914, 1.02718486, -0.64081261,
      0.00000000, 0.00000000, 0.00000000,  0.00000000, 0.00000000, 0.00000000,
      0.00000000, 0.00000000, 0.00000000,  0.00000000, 0.00000000, 0.00000000,
      0.00000000, 0.00000000, 0.00000000,  0.00000000, 0.00000000, 0.00000000,
      0.00000000, 0.00000000, 0.00000000,  0.00000000, 0.00000000, 0.00000000,
      0.00000000, 0.00000000, 0.00000000,  0.00000000, 0.00000000, 0.00000000,
      0.00000000, 0.00000000, 0.00000000,  0.00000000, 0.00000000, 0.00000000};

  std::vector<double> refUSR(refValues,
                             refValues + sizeof(refValues) / sizeof(double));
  std::string rdbase = getenv("RDBASE");
  std::string fname1 =
      rdbase + "/Code/GraphMol/Descriptors/test_data/cyclohexane.mol";
  mol = MolFileToMol(fname1, true, false, true);
  std::vector<std::vector<unsigned int>> atomIds;
  std::vector<double> myUSR(12);
  USRCAT(*mol, myUSR, atomIds);
  for (unsigned int i = 0; i < myUSR.size(); ++i) {
    TEST_ASSERT(feq(myUSR[i], refUSR[i]));
  }

  atomIds.resize(4);
  unsigned int h[6] = {0, 1, 2, 3, 4, 5};
  std::vector<unsigned int> hydrophobic(h,
                                        h + sizeof(h) / sizeof(unsigned int));
  atomIds[0] = hydrophobic;
  USRCAT(*mol, myUSR, atomIds);
  for (unsigned int i = 0; i < myUSR.size(); ++i) {
    TEST_ASSERT(feq(myUSR[i], refUSR[i]));
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testGithub1702() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Github #1702: AUTOCORR2D.h not installed "
                           "unless RDK_BUILD_DESCRIPTORS3D but is required"
                        << std::endl;
  std::vector<double> descriptor(12);

  ROMol *mol = SmilesToMol("C1CCCCC1");
  std::vector<double> dvals;
  AUTOCORR2D(*mol, dvals);
  TEST_ASSERT(dvals.size() == 192);
  delete mol;
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testGithub1973() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "    Test Github #1973: support P and S terms for TPSA calculation"
      << std::endl;

  {
    // Some examples from the original publication, including the non S and P
    // values calculated with RDKit (which match what's in table 3 of the
    // original publication) and then values with S and P contributions added by
    // hand.
    std::vector<std::string> smiles = {
        "Cc1ccccc1N1C(=O)c2cc(S(N)(=O)=O)c(Cl)cc2NC1C", "O=C(O)P(=O)(O)O",
        "O=C(O)c1cc(N=Nc2ccc(S(=O)(=O)Nc3ccccn3)cc2)ccc1O"};
    std::vector<double> orig_tpsa = {92.50, 94.83, 141.31};
    std::vector<double> new_tpsa = {100.88, 104.64, 149.69};
    for (unsigned int i = 0; i < smiles.size(); ++i) {
      std::unique_ptr<ROMol> mol(SmilesToMol(smiles[i]));
      TEST_ASSERT(mol);
      auto oTPSA = calcTPSA(*mol);
      TEST_ASSERT(feq(oTPSA, orig_tpsa[i], 0.01));
      auto nTPSA = calcTPSA(*mol, true, true);
      // std::cerr << smiles[i] << " " << new_tpsa[i] << " " << nTPSA <<
      // std::endl;
      TEST_ASSERT(feq(nTPSA, new_tpsa[i], 0.01));
    }
  }
  {
    // Some examples constructed by hand
    std::vector<std::string> smiles = {"c1ccccc1S", "c1cscc1",    "CC(=S)C",
                                       "CSC",       "CS(=O)C",    "CP(C)C",
                                       "CP=O",      "CP(C)(C)=O", "C[PH](C)=O"};
    std::vector<double> orig_tpsa = {0,   0,     0,     0,    17.07,
                                     0.0, 17.07, 17.07, 17.07};
    std::vector<double> new_tpsa = {38.8,  28.24, 32.09, 25.30, 36.28,
                                    13.59, 51.21, 26.88, 40.54};
    for (unsigned int i = 0; i < smiles.size(); ++i) {
      std::unique_ptr<ROMol> mol(SmilesToMol(smiles[i]));
      TEST_ASSERT(mol);
      auto oTPSA = calcTPSA(*mol);
      TEST_ASSERT(feq(oTPSA, orig_tpsa[i], 0.01));
      auto nTPSA = calcTPSA(*mol, true, true);
      // std::cerr << smiles[i] << " " << new_tpsa[i] << " " << nTPSA <<
      // std::endl;
      TEST_ASSERT(feq(nTPSA, new_tpsa[i], 0.01));
    }
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testGithub2948() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "    Test Github #2948: Empty molecule has non-zero LabuteASA"
      << std::endl;

  {
    ROMol m;
    auto asa = calcLabuteASA(m);
    TEST_ASSERT(feq(asa, 0, 0.0001))
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
int main() {
  RDLog::InitLogs();
#if 1
  test1();
  test2();
  testIssue262();
  test3();
  test3a();
  testLabute();
  testTPSA();
  testLipinski1();
  testVSADescriptors();
  testMolFormula();
  testIssue3415534();
  testIssue3433771();
  testMultiThread();
  testIssue252();
  testChiVs();
  testChiNs();
  testHallKierAlpha();
  testKappa1();
  testKappa2();
  testKappa3();
  testCrippenContribs();
  testRingDescriptors();
  testMiscCountDescriptors();
  testMQNs();
  testGitHubIssue56();
  testGitHubIssue92();
  testGitHubIssue463();
  testSpiroAndBridgeheads();
  testGitHubIssue694();
  testProperties();
  testPropertyQueries();
  testStereoCounting();
  testUSRDescriptor();
  testGithub1702();
#endif
  testGithub1973();
  testGithub2948();
}
