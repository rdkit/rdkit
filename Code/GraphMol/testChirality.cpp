//
//   Copyright (C) 2007-2017 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
// There are chirality test cases spread all over the place. Many of the
// tests here are repeats, but it's good to have everything in one place.
#include <RDGeneral/test.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
//#include <boost/log/functions.hpp>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Canon.h>
#include <GraphMol/new_canon.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/FileParsers/MolSupplier.h>

#include <iostream>

using namespace RDKit;
using namespace std;

void testMol1() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "CIP codes from a mol file (1)" << std::endl;
  std::string rdbase = getenv("RDBASE");
  RWMol *m;
  std::string fName, smi;
  std::string cip;

  // start with SMILES:
  BOOST_LOG(rdInfoLog) << " >>>>>>>>>>>>> smiles 1 <<<<<<<<<<<<<< "
                       << std::endl;
  smi = "O[C@@H](N)I";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  smi = "[C@H](O)(N)I";
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  MolOps::removeStereochemistry(*m);
  TEST_ASSERT(!m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));

  BOOST_LOG(rdInfoLog) << " >>>>>>>>>>>>> mol file <<<<<<<<<<<<<< "
                       << std::endl;
  delete m;
  fName =
      rdbase + "/Code/GraphMol/FileParsers/test_data/ChiralityAndBondDir1a.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  delete m;
  fName =
      rdbase + "/Code/GraphMol/FileParsers/test_data/ChiralityAndBondDir1b.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  MolOps::removeStereochemistry(*m);
  TEST_ASSERT(!m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));

  delete m;
  fName =
      rdbase + "/Code/GraphMol/FileParsers/test_data/ChiralityAndBondDir2a.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  MolOps::removeStereochemistry(*m);
  TEST_ASSERT(!m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));

  delete m;
  fName =
      rdbase + "/Code/GraphMol/FileParsers/test_data/ChiralityAndBondDir2b.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  MolOps::removeStereochemistry(*m);
  TEST_ASSERT(!m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  delete m;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
};

void testRoundTrip() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "CIP codes from a mol->smiles conversion (1)"
                       << std::endl;
  std::string rdbase = getenv("RDBASE");
  RWMol *m;
  std::string fName, smi, smi2;
  std::string cip;

  // start with SMILES:
  smi = "O[C@@H](N)I";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  smi = MolToSmiles(*m, true);
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  smi2 = MolToSmiles(*m, true);
  TEST_ASSERT(smi == smi2);

  smi = "[C@H](O)(N)I";
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
#if 1
  smi = MolToSmiles(*m, true);
  BOOST_LOG(rdInfoLog) << "smiout: " << smi << std::endl;
  TEST_ASSERT(smi == "N[C@H](O)I");
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  smi2 = MolToSmiles(*m, true);
  TEST_ASSERT(smi == smi2);
#endif

  BOOST_LOG(rdInfoLog) << " >>>>>>>>>>>>> mol file <<<<<<<<<<<<<< "
                       << std::endl;
  delete m;
  fName =
      rdbase + "/Code/GraphMol/FileParsers/test_data/ChiralityAndBondDir1a.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
#if 1
  smi = MolToSmiles(*m, true);
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  smi2 = MolToSmiles(*m, true);
  TEST_ASSERT(smi == smi2);
#endif

  delete m;
  fName =
      rdbase + "/Code/GraphMol/FileParsers/test_data/ChiralityAndBondDir1b.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
#if 1
  smi = MolToSmiles(*m, true);
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  smi2 = MolToSmiles(*m, true);
  TEST_ASSERT(smi == smi2);
#endif

  delete m;
  fName =
      rdbase + "/Code/GraphMol/FileParsers/test_data/ChiralityAndBondDir2a.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
#if 1
  smi = MolToSmiles(*m, true);
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  smi2 = MolToSmiles(*m, true);
  TEST_ASSERT(smi == smi2);
#endif

  delete m;
  fName =
      rdbase + "/Code/GraphMol/FileParsers/test_data/ChiralityAndBondDir2b.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
#if 1
  smi = MolToSmiles(*m, true);
  delete m;
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  smi2 = MolToSmiles(*m, true);
  TEST_ASSERT(smi == smi2);
#endif

  delete m;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
};

void testMol2() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "CIP codes from a mol file (2)" << std::endl;
  std::string rdbase = getenv("RDBASE");
  RWMol *m;
  std::string fName, smi;
  std::string cip;

  // start with SMILES:
  BOOST_LOG(rdInfoLog) << " >>>>>>>>>>>>> smiles 1 <<<<<<<<<<<<<< "
                       << std::endl;
  smi = "[C@]1(SC[C@@]([H])(F)[C@]1(Br)O)([I])[H]";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 9);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  TEST_ASSERT(m->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  TEST_ASSERT(m->getAtomWithIdx(5)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(5)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  // same molecule, H combined with the first atom (reproduces
  // exact situation in upcoming mol file)
  BOOST_LOG(rdInfoLog) << " >>>>>>>>>>>>> smiles 2 <<<<<<<<<<<<<< "
                       << std::endl;
  delete m;
  smi = "[C@@H]1(SC[C@@]([H])(F)[C@]1(Br)O)([I])";
  m = SmilesToMol(smi);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 9);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  TEST_ASSERT(m->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  TEST_ASSERT(m->getAtomWithIdx(5)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(5)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  delete m;
  BOOST_LOG(rdInfoLog) << " >>>>>>>>>>>>> mol file <<<<<<<<<<<<<< "
                       << std::endl;
  fName = rdbase + "/Code/GraphMol/FileParsers/test_data/Issue142b.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 9);
  MolOps::assignStereochemistry(*m);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");
  TEST_ASSERT(m->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  delete m;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
};

void testSmiles1() {
  ROMol *mol;
  std::string smi, cip;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "CIP codes from SMILES" << std::endl;

  smi = "F[C@](Cl)(Br)I";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

  delete mol;
  smi = "F[C@](Br)(I)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

  delete mol;
  smi = "F[C@](I)(Cl)Br";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

  delete mol;
  smi = "Cl[C@](Br)(F)I";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

  delete mol;
  smi = "Cl[C@](F)(I)Br";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

  delete mol;
  smi = "I[C@](F)(Br)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

  delete mol;
  smi = "I[C@](Br)(Cl)F";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

  delete mol;
  smi = "F[C@@](Br)(Cl)I";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

  delete mol;
  smi = "F[C@@](Cl)(I)Br";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

  delete mol;
  smi = "Cl[C@@](Br)(I)F";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

  delete mol;
  smi = "Cl[C@@](F)(Br)I";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

  delete mol;
  smi = "[C@@](Cl)(F)(Br)I";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

  delete mol;
  smi = "F[C@H](Cl)Br";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  delete mol;
  smi = "Br[C@H](F)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  delete mol;
  smi = "Br[C@]([H])(F)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  delete mol;
  smi = "Br[C@](F)(Cl)[H]";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  delete mol;
  smi = "Br[C@]1(F)(Cl).[H]1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  delete mol;
  smi = "Br[C@H]1Cl.F1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  delete mol;
  smi = "Br[C@]12Cl.F2.[H]1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  delete mol;
  smi = "Br[C@]21Cl.F1.[H]2";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  delete mol;
  smi = "Br[C@]12Cl.F1.[H]2";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

  delete mol;
  smi = "[C@@](C)(Br)(F)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  delete mol;
  smi = "[C@@]([H])(Br)(F)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  delete mol;
  smi = "[C@@H](Br)(F)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  delete mol;
  smi = "[H][C@@](Br)(F)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getAtomWithIdx(0)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol);
  TEST_ASSERT(mol->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  mol->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  delete mol;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testChiralityCleanup() {
  ROMol *mol, *mol2;
  Atom *chiral_center;
  std::string smi, cip;

  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "chirality cleanup" << std::endl;

  smi = "F[C@H+](Cl)(Br)I";
  mol = SmilesToMol(smi, false, false);
  mol2 = MolOps::removeHs(*mol, false, false);
  delete mol;
  mol = mol2;
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol, true);
  TEST_ASSERT(!mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  delete mol;

  smi = "F[C@+](C)(Cl)(Br)I";
  mol = SmilesToMol(smi, false, false);
  mol2 = MolOps::removeHs(*mol, false, false);
  delete mol;
  mol = mol2;
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  MolOps::assignStereochemistry(*mol, true);
  TEST_ASSERT(!mol->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  TEST_ASSERT(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  delete mol;

  // The remaining examples are for github #1614 :
  // AssignStereochemistry incorrectly removing CIS/TRANS bond stereo

  // cleanIt=true, force=true should NOT remove this manual cis/trans stereo
  // assignment
  smi = "CC=CC(Cl)C";
  mol = SmilesToMol(smi);
  mol->getAtomWithIdx(4)->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
  TEST_ASSERT(mol->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
  mol->getBondWithIdx(1)->setStereoAtoms(0, 3);
  mol->getBondWithIdx(1)->setStereo(Bond::STEREOTRANS);
  std::cerr << "--------------------------" << std::endl;
  MolOps::setDoubleBondNeighborDirections(*mol);
  std::cerr << "--------------------------" << std::endl;
  MolOps::assignStereochemistry(*mol, true, true);
  std::cerr << MolToSmiles(*mol, true) << std::endl;
  TEST_ASSERT(MolToSmiles(*mol, true) == "C/C=C/C(C)Cl");
  delete mol;

  // cleanIt=true, force=true should NOT remove this manual cis/trans stereo
  // assignment
  smi = "CC(F)=CC(Cl)C";
  mol = SmilesToMol(smi);
  mol->getAtomWithIdx(4)->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
  TEST_ASSERT(mol->getBondWithIdx(2)->getBondType() == Bond::DOUBLE);
  mol->getBondWithIdx(2)->setStereoAtoms(0, 4);
  mol->getBondWithIdx(2)->setStereo(Bond::STEREOTRANS);
  MolOps::setDoubleBondNeighborDirections(*mol);
  MolOps::assignStereochemistry(*mol, true, true);
  // std::cerr << MolToSmiles(*mol, true) << std::endl;
  TEST_ASSERT(MolToSmiles(*mol, true) == "C/C(F)=C/[C@H](C)Cl");
  delete mol;

  // cleanIt=true, force=true should NOT remove this manual cis/trans stereo
  // assignment
  smi = "CC(F)=CC(Cl)C";
  mol = SmilesToMol(smi);
  mol->getAtomWithIdx(4)->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
  TEST_ASSERT(mol->getBondWithIdx(2)->getBondType() == Bond::DOUBLE);
  mol->getBondWithIdx(2)->setStereoAtoms(2, 4);
  mol->getBondWithIdx(2)->setStereo(Bond::STEREOCIS);
  MolOps::setDoubleBondNeighborDirections(*mol);
  MolOps::assignStereochemistry(*mol, true, true);
  // std::cerr << MolToSmiles(*mol, true) << std::endl;
  TEST_ASSERT(MolToSmiles(*mol, true) == "C/C(F)=C/[C@H](C)Cl");
  delete mol;

  // cleanIt=true, force=true should clean up these manual assignments
  smi = "FC(F)=C(Cl)Cl";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getBondWithIdx(2)->getBondType() == Bond::DOUBLE);
  mol->getBondWithIdx(2)->setStereoAtoms(2, 4);
  mol->getBondWithIdx(2)->setStereo(Bond::STEREOCIS);
  MolOps::setDoubleBondNeighborDirections(*mol);
  MolOps::assignStereochemistry(*mol, true, true);
  TEST_ASSERT(MolToSmiles(*mol, true) == "FC(F)=C(Cl)Cl");
  delete mol;

  smi = "FC=C1CCCCC1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getBondWithIdx(1)->getBondType() == Bond::DOUBLE);
  mol->getBondWithIdx(1)->setStereoAtoms(0, 3);
  mol->getBondWithIdx(1)->setStereo(Bond::STEREOTRANS);
  MolOps::setDoubleBondNeighborDirections(*mol);
  MolOps::assignStereochemistry(*mol, true, true);
  TEST_ASSERT(MolToSmiles(*mol, true) == "FC=C1CCCCC1");
  delete mol;

  smi = "C1CCCCC1=CF";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getBondWithIdx(5)->getBondType() == Bond::DOUBLE);
  mol->getBondWithIdx(5)->setStereoAtoms(4, 7);
  mol->getBondWithIdx(5)->setStereo(Bond::STEREOCIS);
  MolOps::setDoubleBondNeighborDirections(*mol);
  MolOps::assignStereochemistry(*mol, true, true);
  TEST_ASSERT(MolToSmiles(*mol, true) == "FC=C1CCCCC1");
  delete mol;

///////////////////////////
// Pseudo-stereo test cases
// - bond stereo dependent on other bond stereo breaking symmetry - DOESN'T
// CURRENTLY WORK
// - bond stereo dependent on atom stereo breaking symmetry - DOESN'T CURRENTLY
// WORK
// - atom stereo dependent on other atom stereo breaking symmetry - WORKS
// - atom stereo dependent on bond stereo breaking symmetry - WORKS
///////////////////////////

// Everything ifdef'd USE_NEW_STEREOCHEMISTRY are test cases that
// should 'hopefully' work when switching to the new_canon.h ranking
// for stereo perception. However, making USE_NEW_STEREOCHEMISTRY
// copasetic with the current behavior is non-trivial, possibly
// impossible?

// bond stereo dependent on other bond stereo breaking symmetry
// cleanIt=true, force=true should NOT remove this trickier case of
// pseudo-stereo
#ifdef USE_NEW_STEREOCHEMISTRY
  smi = "CCC=CC(C=CCC)=C(CC)CO";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getBondWithIdx(2)->getBondType() == Bond::DOUBLE);
  mol->getBondWithIdx(2)->setStereoAtoms(1, 4);
  mol->getBondWithIdx(2)->setStereo(Bond::STEREOTRANS);

  TEST_ASSERT(mol->getBondWithIdx(5)->getBondType() == Bond::DOUBLE);
  mol->getBondWithIdx(5)->setStereoAtoms(4, 7);
  mol->getBondWithIdx(5)->setStereo(Bond::STEREOCIS);

  TEST_ASSERT(mol->getBondWithIdx(8)->getBondType() == Bond::DOUBLE);
  mol->getBondWithIdx(8)->setStereoAtoms(3, 10);
  mol->getBondWithIdx(8)->setStereo(Bond::STEREOCIS);

  MolOps::setDoubleBondNeighborDirections(*mol);
  MolOps::assignStereochemistry(*mol, true, true);
  BOOST_LOG(rdInfoLog) << MolToSmiles(*mol, true) << std::endl;
  TEST_ASSERT(MolToSmiles(*mol, true) == "CC/C=C\\C(\\C=C\\CC)=C(\\CC)CO");
  delete mol;

  // make sure there isn't a difference when the bond stereo is set from SMILES
  smi = "CC/C=C\\C(\\C=C\\CC)=C(CC)CO";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getBondWithIdx(8)->getBondType() == Bond::DOUBLE);
  TEST_ASSERT(mol->getBondWithIdx(8)->getStereo() == Bond::STEREONONE);
  mol->getBondWithIdx(8)->setStereoAtoms(3, 10);
  mol->getBondWithIdx(8)->setStereo(Bond::STEREOTRANS);

  MolOps::setDoubleBondNeighborDirections(*mol);
  MolOps::assignStereochemistry(*mol, true, true);
  BOOST_LOG(rdInfoLog) << MolToSmiles(*mol, true) << std::endl;
  TEST_ASSERT(MolToSmiles(*mol, true) == "CC/C=C\\C(\\C=C\\CC)=C(\\CC)CO");
  delete mol;
#endif

  // cleanIt=true, force=true should remove this trickier case of pseudo-stereo
  smi = "CCC=CC(=C(CC)CO)C=CCC";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getBondWithIdx(2)->getBondType() == Bond::DOUBLE);
  mol->getBondWithIdx(2)->setStereoAtoms(1, 4);
  mol->getBondWithIdx(2)->setStereo(Bond::STEREOCIS);

  TEST_ASSERT(mol->getBondWithIdx(4)->getBondType() == Bond::DOUBLE);
  mol->getBondWithIdx(4)->setStereoAtoms(10, 8);
  mol->getBondWithIdx(4)->setStereo(Bond::STEREOCIS);

  TEST_ASSERT(mol->getBondWithIdx(10)->getBondType() == Bond::DOUBLE);
  mol->getBondWithIdx(10)->setStereoAtoms(4, 12);
  mol->getBondWithIdx(10)->setStereo(Bond::STEREOCIS);
  std::cerr << "3>>>--------------------------" << std::endl;
  MolOps::setDoubleBondNeighborDirections(*mol);
  std::cerr << "<<<--------------------------" << std::endl;
  MolOps::assignStereochemistry(*mol, true, true);
  BOOST_LOG(rdInfoLog) << MolToSmiles(*mol, true) << std::endl;
  TEST_ASSERT(MolToSmiles(*mol, true) == "CC/C=C\\C(/C=C\\CC)=C(CC)CO");
  delete mol;

  // make sure there isn't a difference when the bond stereo is set from SMILES
  smi = "CC/C=C\\C(\\C=C/CC)=C(CC)CO";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getBondWithIdx(8)->getBondType() == Bond::DOUBLE);
  TEST_ASSERT(mol->getBondWithIdx(8)->getStereo() == Bond::STEREONONE);

  mol->getBondWithIdx(8)->setStereoAtoms(3, 10);
  mol->getBondWithIdx(8)->setStereo(Bond::STEREOCIS);
  MolOps::setDoubleBondNeighborDirections(*mol);
  MolOps::assignStereochemistry(*mol, true, true);
  // BOOST_LOG(rdInfoLog) << MolToSmiles(*mol, true) << std::endl;
  TEST_ASSERT(MolToSmiles(*mol, true) == "CC/C=C\\C(/C=C\\CC)=C(CC)CO");
  delete mol;

#ifdef USE_NEW_STEREOCHEMISTRY
  // bond stereo dependent on atom stereo breaking symmetry
  // cleanIt=true, force=true should NOT remove this trickier case of
  // pseudo-stereo
  smi = "CCC(CO)=C([C@H](C)F)[C@@H](C)F";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getBondWithIdx(4)->getBondType() == Bond::DOUBLE);
  mol->getBondWithIdx(4)->setStereoAtoms(1, 6);
  mol->getBondWithIdx(4)->setStereo(Bond::STEREOCIS);

  MolOps::setDoubleBondNeighborDirections(*mol);
  MolOps::assignStereochemistry(*mol, true, true);
  BOOST_LOG(rdInfoLog) << MolToSmiles(*mol, true) << std::endl;
  TEST_ASSERT(MolToSmiles(*mol, true) == "CC/C(CO)=C(/[C@@H](C)F)[C@H](C)F");
  delete mol;
#endif

  // cleanIt=true, force=true should remove this manual assignment of bond
  // stereochemistry since it's a pseudo-stereo center
  smi = "CCC(CO)=C([C@@H](C)F)[C@@H](C)F";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getBondWithIdx(4)->getBondType() == Bond::DOUBLE);
  mol->getBondWithIdx(4)->setStereoAtoms(1, 6);
  mol->getBondWithIdx(4)->setStereo(Bond::STEREOCIS);

  MolOps::setDoubleBondNeighborDirections(*mol);
  MolOps::assignStereochemistry(*mol, true, true);
  TEST_ASSERT(MolToSmiles(*mol, true) == "CCC(CO)=C([C@@H](C)F)[C@@H](C)F");
  delete mol;

  // atom stereo dependent on other atom stereo breaking symmetry
  // cleanIt=true, force=true should remove this manual assignment of atom
  // stereochemistry since it's a pseudo-stereo center
  smi = "C[C@H]1CC(N2CCc3ccc(N)cc3C2)C[C@H](C)C1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getRingInfo()->isAtomInRingOfSize(3, 6));
  chiral_center = mol->getAtomWithIdx(3);
  TEST_ASSERT(chiral_center->getAtomicNum() == 6);
  TEST_ASSERT(chiral_center->getDegree() == 3);
  TEST_ASSERT(!chiral_center->getIsAromatic());

  chiral_center->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
  MolOps::setDoubleBondNeighborDirections(*mol);
  MolOps::assignStereochemistry(*mol, true, true);
  TEST_ASSERT(chiral_center->getChiralTag() == Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(MolToSmiles(*mol, true) ==
              "C[C@H]1CC(N2CCc3ccc(N)cc3C2)C[C@H](C)C1");
  delete mol;

  // cleanIt=true, force=true should NOT remove this manual assignment of atom
  // stereochemistry since it's a pseudo-stereo center
  smi = "C[C@@H]1CC(N2CCc3ccc(N)cc3C2)C[C@H](C)C1";
  mol = SmilesToMol(smi);
  TEST_ASSERT(mol->getRingInfo()->isAtomInRingOfSize(3, 6));
  chiral_center = mol->getAtomWithIdx(3);
  TEST_ASSERT(chiral_center->getAtomicNum() == 6);
  TEST_ASSERT(chiral_center->getDegree() == 3);
  TEST_ASSERT(!chiral_center->getIsAromatic());

  chiral_center->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
  MolOps::setDoubleBondNeighborDirections(*mol);
  MolOps::assignStereochemistry(*mol, true, true);
  TEST_ASSERT(chiral_center->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);
  TEST_ASSERT(MolToSmiles(*mol, true) ==
              "C[C@H]1C[C@@H](C)C[C@H](N2CCc3ccc(N)cc3C2)C1");
  delete mol;

  // atom stereo dependent on bond stereo breaking symmetry
  // cleanIt=true, force=true should remove this manual assignment of atom
  // stereochemistry since it's a pseudo-stereo center
  smi = "C/C=C/C(/C=C/C)(CC)CO";
  mol = SmilesToMol(smi);
  chiral_center = mol->getAtomWithIdx(3);
  TEST_ASSERT(chiral_center->getAtomicNum() == 6);
  TEST_ASSERT(chiral_center->getDegree() == 4);

  chiral_center->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
  MolOps::setDoubleBondNeighborDirections(*mol);
  MolOps::assignStereochemistry(*mol, true, true);
  TEST_ASSERT(chiral_center->getChiralTag() == Atom::CHI_UNSPECIFIED);
  TEST_ASSERT(MolToSmiles(*mol, true) == "C/C=C/C(/C=C/C)(CC)CO");
  delete mol;

  // atom stereo dependent on bond stereo breaking symmetry
  // cleanIt=true, force=true should NOT remove this manual assignment of atom
  // stereochemistry since it's a pseudo-stereo center
  smi = "C/C=C\\C(/C=C/C)(CC)CO";
  mol = SmilesToMol(smi);
  chiral_center = mol->getAtomWithIdx(3);
  TEST_ASSERT(chiral_center->getAtomicNum() == 6);
  TEST_ASSERT(chiral_center->getDegree() == 4);

  chiral_center->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
  MolOps::setDoubleBondNeighborDirections(*mol);
  MolOps::assignStereochemistry(*mol, true, true);
  TEST_ASSERT(chiral_center->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);
  TEST_ASSERT(MolToSmiles(*mol, true) == "C/C=C\\[C@@](/C=C/C)(CC)CO");
  delete mol;

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testRingStereochemistry() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "test ring stereochemistry " << std::endl;
  // NOTE: this test is for correctness, not canonicality
  {
    std::string smi = "B[C@H]1CC[C@H](C)CC1";
    RWMol *m = SmilesToMol(smi);
    std::string smi1 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi << " " << smi1 << std::endl;
    TEST_ASSERT(smi1 == "B[C@H]1CC[C@H](C)CC1");

    delete m;
#if 0
    smi="B[C@@H]1CC[C@@H](C)CC1";
    m = SmilesToMol(smi);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi2<<" "<<smi1<<std::endl;
    TEST_ASSERT(smi2==smi);
    delete m;
#endif
  }

  {
    std::string smi = "C1[C@@H](B)CC[C@H](C)C1";
    RWMol *m = SmilesToMol(smi);
    std::string smi1 = MolToSmiles(*m, true);
    smi = "B[C@H]1CC[C@H](C)CC1";
    BOOST_LOG(rdInfoLog) << " : " << smi << " " << smi1 << std::endl;
    TEST_ASSERT(smi1 == smi);
    delete m;
#if 0
    smi="C1[C@H](B)CC[C@@H](C)C1";
    m = SmilesToMol(smi);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi2<<" "<<smi1<<std::endl;
    TEST_ASSERT(smi2==smi1);
    delete m;
#endif
  }

  {
    std::string smi = "C[C@H]1CC[C@H](F)CC1";
    RWMol *m = SmilesToMol(smi);
    std::string smi1 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi << " " << smi1 << std::endl;
    TEST_ASSERT(smi1 == "C[C@H]1CC[C@H](F)CC1");
    delete m;
#if 0
    smi="C[C@@H]1CC[C@@H](F)CC1";
    m = SmilesToMol(smi);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi2<<" "<<smi1<<std::endl;
    TEST_ASSERT(smi2==smi1);
    delete m;
#endif
  }

  {
    std::string smi = "F[C@H]1CC[C@H](C)CC1";
    RWMol *m = SmilesToMol(smi);
    std::string smi1 = MolToSmiles(*m, true);
    delete m;
#if 0
    smi="F[C@@H]1CC[C@@H](C)CC1";
    m = SmilesToMol(smi);
    std::string smi2=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi2<<" "<<smi1<<std::endl;
    TEST_ASSERT(smi2==smi1);
    delete m;
#endif
  }

  {
    std::string smi = "F[C@H]1CC[C@](C)(C)CC1";
    RWMol *m = SmilesToMol(smi);
    std::string smi1 = MolToSmiles(*m, true);
    delete m;
    smi = "FC1CCC(C)(C)CC1";
    m = SmilesToMol(smi);
    std::string smi2 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi2 << " " << smi1 << std::endl;
    TEST_ASSERT(smi2 == smi1);
    delete m;
  }

  {
    std::string smi = "C1C[C@H]2CC[C@@H]1CC2";
    RWMol *m = SmilesToMol(smi);
    std::string smi1 = MolToSmiles(*m, true);
    delete m;
    smi = "C1CC2CCC1CC2";
    m = SmilesToMol(smi);
    std::string smi2 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi2 << " " << smi1 << std::endl;
    TEST_ASSERT(smi2 == smi1);
    delete m;
  }

  {
    std::string smi = "C[C@]12CC[C@](C)(CC1)CC2";
    RWMol *m = SmilesToMol(smi);
    std::string smi1 = MolToSmiles(*m, true);
    delete m;
    smi = "CC12CCC(C)(CC1)CC2";
    m = SmilesToMol(smi);
    std::string smi2 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi2 << " " << smi1 << std::endl;
    TEST_ASSERT(smi2 != smi1);
    delete m;
  }

  {
    // make sure we aren't removing stereochem that should still be there
    std::string smi = "C[C@@]12CC[C@@](C)(NC1)OC2";
    RWMol *m = SmilesToMol(smi);
    std::string smi1 = MolToSmiles(*m, true);
    delete m;
    smi = "CC12CCC(C)(NC1)OC2";
    m = SmilesToMol(smi);
    std::string smi2 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi2 << " " << smi1 << std::endl;
    TEST_ASSERT(smi2 != smi1);
    delete m;
  }

#if 0
  // FIX : these tests do not pass
  {
    std::string smi = "C[C@H]1CC[C@H](C)CC1";
    RWMol *m = SmilesToMol(smi);
    std::string smi1=MolToSmiles(*m,true);
    BOOST_LOG(rdInfoLog)<<" : "<<smi<<" "<<smi1<<std::endl;
    TEST_ASSERT(smi1==smi);
    delete m;
  }

  {
    std::string smi = "C1[C@@H](C)CC[C@H](C)C1";
    RWMol *m = SmilesToMol(smi);
    m->debugMol(std::cerr);
    std::string smi1=MolToSmiles(*m,true);
    smi = "C[C@H]1CC[C@H](C)CC1";
    BOOST_LOG(rdInfoLog)<<" : "<<smi<<" "<<smi1<<std::endl;
    m->debugMol(std::cerr);
    TEST_ASSERT(smi1==smi);
    delete m;
  }
#endif

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testChiralityFrom3D() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "chirality perception from 3D coordinates: "
                       << std::endl;

  std::string rdbase = getenv("RDBASE");
  RWMol *m;
  std::string fName, smi;
  std::string cip;

  fName = rdbase + "/Code/GraphMol/test_data/chi3d_r1.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 5);

  MolOps::assignChiralTypesFrom3D(*m);
  MolOps::assignStereochemistry(*m, true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  delete m;
  fName = rdbase + "/Code/GraphMol/test_data/chi3d_s1.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 5);

  MolOps::assignChiralTypesFrom3D(*m);
  MolOps::assignStereochemistry(*m, true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

  delete m;
  fName = rdbase + "/Code/GraphMol/test_data/chi3d_r2.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);

  MolOps::assignChiralTypesFrom3D(*m);
  MolOps::assignStereochemistry(*m, true);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW);
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  delete m;
  fName = rdbase + "/Code/GraphMol/test_data/chi3d_s2.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 4);

  MolOps::assignChiralTypesFrom3D(*m);
  MolOps::assignStereochemistry(*m, true);
  TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
  TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);
  m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "S");

  delete m;
  fName = rdbase + "/Code/GraphMol/test_data/chi3d_r1_bad.mol";
  m = MolFileToMol(fName);
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 5);

  // this molecule starts out with incorrect stereochemistry (e.g. the bond
  // wedging does not match the 3D structure.
  // This is handled automatically by the mol file parser as of github #1679,
  // so we don't need to worry about it anymore
  MolOps::assignStereochemistry(*m, true);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
  TEST_ASSERT(cip == "R");

  delete m;
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIterativeChirality() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "iterative chirality (sf.net issue 1931470): "
                       << std::endl;

  std::string rdbase = getenv("RDBASE");

// unless otherwise noted, the R/S and Z/E assignments here
// match Marvin and ChemDraw.
#if 1
  {  // atom-chirality -> atom-chirality
    std::string cip;

    std::string fName = rdbase + "/Code/GraphMol/test_data/iChi1a.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 9);

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");

    TEST_ASSERT(m->getAtomWithIdx(5)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(5)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    std::string smi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi1 << " " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);

    delete m;
  }

  {  // atom-chirality -> atom-chirality
    std::string cip;

    std::string fName = rdbase + "/Code/GraphMol/test_data/iChi1b.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 9);

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");

    TEST_ASSERT(m->getAtomWithIdx(5)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(5)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");

    std::string smi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi1 << " " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);

    delete m;
  }

  {  // atom-chirality -> atom-chirality
    std::string cip;

    std::string fName = rdbase + "/Code/GraphMol/test_data/iChi1c.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 9);

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");

    TEST_ASSERT(m->getAtomWithIdx(5)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(5)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");

    TEST_ASSERT(!m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));

#if 1  // this fails due to sf.net bug 1896935
    std::string smi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi1 << " " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);
#endif

    delete m;
  }

  {  // atom-chirality -> atom-chirality
    std::string cip;

    std::string fName = rdbase + "/Code/GraphMol/test_data/iChi1d.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 9);

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");

    TEST_ASSERT(m->getAtomWithIdx(5)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(5)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");

    TEST_ASSERT(!m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));

#if 1  // this fails due to sf.net bug 1896935
    std::string smi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi1 << " " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);
#endif
    delete m;
  }

  {  // atom-chirality -> atom-chirality
    std::string cip;

    std::string fName = rdbase + "/Code/GraphMol/test_data/iChi1e.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 9);

    TEST_ASSERT(!m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
    TEST_ASSERT(!m->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
    TEST_ASSERT(!m->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));

#if 0  // this fails due to sf.net bug 1896935
    std::cerr<<"m pre -----"<<std::endl;
    m->debugMol(std::cerr);
    std::cerr<<"-----"<<std::endl;
    std::string smi1=MolToSmiles(*m,true);
    std::cerr<<"m post -----"<<std::endl;
    m->debugMol(std::cerr);
    std::cerr<<"-----"<<std::endl;
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::cerr<<"m2 pre -----"<<std::endl;
    m->debugMol(std::cerr);
    std::cerr<<"-----"<<std::endl;
    std::string smi2=MolToSmiles(*m,true);
    std::cerr<<"m post -----"<<std::endl;
    m->debugMol(std::cerr);
    std::cerr<<"-----"<<std::endl;
    BOOST_LOG(rdInfoLog)<<" : "<<smi1<<" "<<smi2<<std::endl;
    TEST_ASSERT(smi1==smi2);
#endif
    delete m;
  }

  {  // bond-stereochem -> atom-chirality
    std::string cip;

    std::string fName = rdbase + "/Code/GraphMol/test_data/iChi2a.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 8);

    TEST_ASSERT(m->getBondBetweenAtoms(2, 5)->getStereo() == Bond::STEREOE);
    TEST_ASSERT(m->getBondBetweenAtoms(0, 4)->getStereo() == Bond::STEREOZ);

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip ==
                "S");  // this value is from ChemDraw, Marvin doesn't tag it.

    std::string smi1 = MolToSmiles(*m, true);

    MolOps::removeStereochemistry(*m);
    TEST_ASSERT(!m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    TEST_ASSERT(m->getBondBetweenAtoms(0, 4)->getStereo() == Bond::STEREONONE);
    TEST_ASSERT(m->getBondBetweenAtoms(2, 5)->getStereo() == Bond::STEREONONE);

    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi1 << " " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);

    delete m;
  }

  {  // bond-stereochem -> atom-chirality
    std::string cip;

    std::string fName = rdbase + "/Code/GraphMol/test_data/iChi2b.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 8);

    TEST_ASSERT(m->getBondBetweenAtoms(2, 5)->getStereo() == Bond::STEREOE);
    TEST_ASSERT(m->getBondBetweenAtoms(0, 4)->getStereo() == Bond::STEREOZ);

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip ==
                "R");  // this value is from ChemDraw, Marvin doesn't tag it.

    std::string smi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi1 << " " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);

    delete m;
  }

  {  // bond-stereochem -> atom-chirality
    std::string cip;

    std::string fName = rdbase + "/Code/GraphMol/test_data/iChi2c.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 8);

    TEST_ASSERT(m->getBondBetweenAtoms(2, 5)->getStereo() == Bond::STEREOE);
    TEST_ASSERT(m->getBondBetweenAtoms(0, 4)->getStereo() == Bond::STEREOE);

    TEST_ASSERT(!m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));

    std::string smi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi1 << " " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);

    delete m;
  }

  {  // bond-stereochem -> atom-chirality
    std::string cip;

    std::string fName = rdbase + "/Code/GraphMol/test_data/iChi2d.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 8);

    TEST_ASSERT(m->getBondBetweenAtoms(2, 5)->getStereo() == Bond::STEREOE);
    TEST_ASSERT(m->getBondBetweenAtoms(0, 4)->getStereo() == Bond::STEREOANY);

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip ==
                "R");  // this value is from ChemDraw, Marvin doesn't tag it.

    std::string smi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi1 << " " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);

    delete m;
  }

  {  // bond-stereochem -> atom-chirality
    std::string cip;

    std::string fName = rdbase + "/Code/GraphMol/test_data/iChi2e.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 8);

    TEST_ASSERT(m->getBondBetweenAtoms(2, 5)->getStereo() == Bond::STEREOANY);
    TEST_ASSERT(m->getBondBetweenAtoms(0, 4)->getStereo() == Bond::STEREOZ);

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip ==
                "S");  // this value is from ChemDraw, Marvin doesn't tag it.

    std::string smi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi1 << " " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);

    delete m;
  }

  {  // atom chirality -> bond stereochemistry
    std::string cip;

    std::string fName = rdbase + "/Code/GraphMol/test_data/iChi3a.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 11);

    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");

    TEST_ASSERT(m->getAtomWithIdx(7)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(7)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    TEST_ASSERT(
        m->getBondBetweenAtoms(1, 2)->getStereo() ==
        Bond::STEREOZ);  // this value is from ChemDraw, Marvin doesn't tag it.

    std::string smi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi1 << " " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);

    delete m;
  }

  {  // atom chirality -> bond stereochemistry
    std::string cip;

    std::string fName = rdbase + "/Code/GraphMol/test_data/iChi3b.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 11);

    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    TEST_ASSERT(m->getAtomWithIdx(7)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(7)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");

    TEST_ASSERT(
        m->getBondBetweenAtoms(1, 2)->getStereo() ==
        Bond::STEREOE);  // this value is from ChemDraw, Marvin doesn't tag it.

    std::string smi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi1 << " " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);

    delete m;
  }

  {  // atom chirality -> bond stereochemistry
    std::string cip;

    std::string fName = rdbase + "/Code/GraphMol/test_data/iChi3c.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 11);

    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");

    TEST_ASSERT(m->getAtomWithIdx(7)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(7)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");

    TEST_ASSERT(m->getBondBetweenAtoms(1, 2)->getStereo() == Bond::STEREONONE);

    std::string smi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi1 << " " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);

    delete m;
  }
#endif

  {  // bond stereochemistry -> bond stereochemistry
    std::string cip;

    std::string fName = rdbase + "/Code/GraphMol/test_data/iChi4a.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 9);

    TEST_ASSERT(m->getBondBetweenAtoms(4, 5)->getStereo() == Bond::STEREOE);
    TEST_ASSERT(m->getBondBetweenAtoms(3, 7)->getStereo() == Bond::STEREOZ);

    TEST_ASSERT(m->getBondBetweenAtoms(0, 1)->getStereo() == Bond::STEREOE);

    std::string smi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi1 << " " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);

    delete m;
  }

  {  // bond stereochemistry -> bond stereochemistry
    std::string cip;

    std::string fName = rdbase + "/Code/GraphMol/test_data/iChi4b.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 9);

    TEST_ASSERT(m->getBondBetweenAtoms(4, 5)->getStereo() == Bond::STEREOZ);
    TEST_ASSERT(m->getBondBetweenAtoms(3, 7)->getStereo() == Bond::STEREOE);

    TEST_ASSERT(m->getBondBetweenAtoms(0, 1)->getStereo() == Bond::STEREOZ);

    std::string smi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi1 << " " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);

    delete m;
  }

  {  // bond stereochemistry -> bond stereochemistry
    std::string cip;

    std::string fName = rdbase + "/Code/GraphMol/test_data/iChi4c.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 9);

    TEST_ASSERT(m->getBondBetweenAtoms(4, 5)->getStereo() == Bond::STEREOE);
    TEST_ASSERT(m->getBondBetweenAtoms(3, 7)->getStereo() == Bond::STEREOE);

    TEST_ASSERT(m->getBondBetweenAtoms(0, 1)->getStereo() == Bond::STEREONONE);

    std::string smi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi1 << " " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);

    delete m;
  }

  {  // bond stereochemistry -> bond stereochemistry
    std::string cip;

    std::string fName = rdbase + "/Code/GraphMol/test_data/iChi4d.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 9);

    TEST_ASSERT(m->getBondBetweenAtoms(4, 5)->getStereo() == Bond::STEREOZ);
    TEST_ASSERT(m->getBondBetweenAtoms(3, 7)->getStereo() == Bond::STEREOZ);

    TEST_ASSERT(m->getBondBetweenAtoms(0, 1)->getStereo() == Bond::STEREONONE);

    std::string smi1 = MolToSmiles(*m, true);
    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    std::string smi2 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi1 << " " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);

    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testBondDirRemoval() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "testing that the removal of bond directions is correct: "
      << std::endl;

  std::string rdbase = getenv("RDBASE");

  {
    std::string cip;

    std::string fName = rdbase + "/Code/GraphMol/test_data/stereoOrder1.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 7);

    TEST_ASSERT(m->getBondBetweenAtoms(1, 2)->getStereo() == Bond::STEREOZ);
    TEST_ASSERT(m->getBondBetweenAtoms(4, 5)->getStereo() == Bond::STEREOE);

    // on input all the single bonds are in the same direction:
    TEST_ASSERT(m->getBondBetweenAtoms(0, 1)->getBondDir() ==
                m->getBondBetweenAtoms(1, 4)->getBondDir());
    TEST_ASSERT(m->getBondBetweenAtoms(2, 3)->getBondDir() ==
                m->getBondBetweenAtoms(1, 4)->getBondDir());
    TEST_ASSERT(m->getBondBetweenAtoms(5, 6)->getBondDir() ==
                m->getBondBetweenAtoms(1, 4)->getBondDir());

    std::string smi1 = MolToSmiles(*m, true);

    // check removal of redundant bond direction information:
    std::vector<unsigned int> oranks(m->getNumAtoms(), 0);
    Canon::rankMolAtoms(*m, oranks);
    std::vector<Canon::AtomColors> colors(m->getNumAtoms());
    Canon::MolStack stack;
    std::vector<unsigned int> ranks(oranks.size());
    for (unsigned int i = 0; i < ranks.size(); ++i) ranks[i] = oranks[i];
    Canon::canonicalizeFragment(*m, 0, colors, ranks, stack);

    TEST_ASSERT(m->getBondBetweenAtoms(0, 1)->getBondDir() == Bond::NONE);
    TEST_ASSERT(m->getBondBetweenAtoms(2, 3)->getBondDir() ==
                m->getBondBetweenAtoms(1, 4)->getBondDir());
    TEST_ASSERT(m->getBondBetweenAtoms(5, 6)->getBondDir() ==
                m->getBondBetweenAtoms(1, 4)->getBondDir());

    std::string smi2 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi1 << " " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);

    delete m;
    m = SmilesToMol(smi1);
    TEST_ASSERT(m);
    smi2 = MolToSmiles(*m, true);
    BOOST_LOG(rdInfoLog) << " : " << smi1 << " " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);

    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIssue2705543() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Issue 2705543: " << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName;
  RWMol *m;
  std::string cip;

  {
    fName = rdbase + "/Code/GraphMol/test_data/Issue2705543.1h.mol";
    m = MolFileToMol(fName, true, false);
    TEST_ASSERT(m);

    MolOps::assignChiralTypesFrom3D(*m);
    MolOps::assignStereochemistry(*m, true);

    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");

    TEST_ASSERT(m->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    TEST_ASSERT(m->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");

    delete m;
  }
  {
    fName = rdbase + "/Code/GraphMol/test_data/Issue2705543.1.mol";
    m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 13);

    MolOps::assignChiralTypesFrom3D(*m);
    MolOps::assignStereochemistry(*m, true);

    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");

    TEST_ASSERT(m->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    TEST_ASSERT(m->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");

    delete m;
  }
  {
    fName = rdbase + "/Code/GraphMol/test_data/Issue2705543.2h.mol";
    m = MolFileToMol(fName, true, false);
    TEST_ASSERT(m);

    MolOps::assignChiralTypesFrom3D(*m);
    MolOps::assignStereochemistry(*m, true);

    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
    TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CCW);
    m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CW);
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");

    TEST_ASSERT(m->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
    TEST_ASSERT(m->getAtomWithIdx(2)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CW);
    m->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
    TEST_ASSERT(m->getAtomWithIdx(3)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CCW);
    m->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    TEST_ASSERT(m->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));
    TEST_ASSERT(m->getAtomWithIdx(4)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CCW);
    m->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");

    delete m;
  }
  {
    fName = rdbase + "/Code/GraphMol/test_data/Issue2705543.2.mol";
    m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 13);

    MolOps::assignChiralTypesFrom3D(*m);
    MolOps::assignStereochemistry(*m, true);

#if 0
    for(unsigned int i=0;i<m->getNumAtoms();++i){
      if(m->getAtomWithIdx(i)->hasProp(common_properties::_CIPCode)){
        m->getAtomWithIdx(i)->getProp(common_properties::_CIPCode,cip);
        std::cerr<<"  >> "<<i<<" "<<cip<<std::endl;
      }
    }
#endif
    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
    TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CCW);
    m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CW);
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");

    TEST_ASSERT(m->getAtomWithIdx(2)->hasProp(common_properties::_CIPCode));
    TEST_ASSERT(m->getAtomWithIdx(2)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CW);
    m->getAtomWithIdx(2)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
    TEST_ASSERT(m->getAtomWithIdx(3)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CCW);
    m->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    TEST_ASSERT(m->getAtomWithIdx(4)->hasProp(common_properties::_CIPCode));
    TEST_ASSERT(m->getAtomWithIdx(4)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CCW);
    m->getAtomWithIdx(4)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");

    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIssue2762917() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Issue 2762917: chirality swap on addHs()"
                       << std::endl;

  std::string rdbase = getenv("RDBASE");

  {
    RWMol *m;
    std::string cip;
    std::string smiles = "[C@@H](C)(Cl)O";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);

    TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CCW);

    MolOps::assignStereochemistry(*m, true);

    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CCW);
    MolOps::assignStereochemistry(*m, true);

    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    delete m;
  }

  {
    RWMol *m;
    std::string cip;
    std::string smiles = "CCC.[C@@H](C)(Cl)O";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);

    TEST_ASSERT(m->getAtomWithIdx(3)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CCW);

    MolOps::assignStereochemistry(*m, true);

    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    MolOps::addHs(*m);

    TEST_ASSERT(m->getAtomWithIdx(3)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CCW);
    MolOps::assignStereochemistry(*m, true);

    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    delete m;
  }

  {
    RWMol *m;
    std::string cip;
    std::string smiles = "[C@@H]([C@H](C)O)(C)O";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);

    TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CCW);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CCW);

    MolOps::assignStereochemistry(*m, true);

    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    MolOps::addHs(*m);

    TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CCW);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CCW);

    MolOps::assignStereochemistry(*m, true);

    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(0)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    delete m;
  }

  {
    RWMol *m;
    std::string cip;
    std::string smiles = "C1CC.[C@@H]1(Cl)O";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);

    TEST_ASSERT(m->getAtomWithIdx(3)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CCW);

    MolOps::assignStereochemistry(*m, true);

    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    MolOps::addHs(*m);

    TEST_ASSERT(m->getAtomWithIdx(3)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CCW);
    MolOps::assignStereochemistry(*m, true);

    TEST_ASSERT(m->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(3)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIssue3009911() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Issue 3009911: bad atom priorities" << std::endl;

  {
    RWMol *m;
    std::string smiles = "F[C@](O)(c1ccccc1)C(=C)CO";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int *ranks;
    ranks = new int[m->getNumAtoms()];
    MolOps::assignStereochemistry(*m, true);
    for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
      unsigned int rank;
      TEST_ASSERT(m->getAtomWithIdx(i)->hasProp(common_properties::_CIPRank))
      m->getAtomWithIdx(i)->getProp(common_properties::_CIPRank, rank);
      ranks[i] = rank;
    }
    // basics:
    TEST_ASSERT(ranks[0] > ranks[1]);
    TEST_ASSERT(ranks[2] > ranks[1]);
    TEST_ASSERT(ranks[0] > ranks[2]);
    // now the key point:
    TEST_ASSERT(ranks[3] < ranks[9]);

    std::string cip;
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");
    delete m;
    delete[] ranks;
  }
  {
    RWMol *m;
    std::string smiles = "COC(C)(OC)[C@](O)(F)C(C)=O";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    int *ranks;
    ranks = new int[m->getNumAtoms()];
    MolOps::assignStereochemistry(*m, true);
    for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
      unsigned int rank;
      TEST_ASSERT(m->getAtomWithIdx(i)->hasProp(common_properties::_CIPRank))
      m->getAtomWithIdx(i)->getProp(common_properties::_CIPRank, rank);
      ranks[i] = rank;
    }
    // basics:
    TEST_ASSERT(ranks[8] > ranks[7]);
    TEST_ASSERT(ranks[7] > ranks[9]);
    TEST_ASSERT(ranks[7] > ranks[2]);
    // FIX: these are the key points, but at the moment they are not handled
    // correctly
    // due to a weakness in the CIP-ranking algorithm.
    // TEST_ASSERT(ranks[2]>ranks[9]);
    std::string cip;
    TEST_ASSERT(m->getAtomWithIdx(6)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(6)->getProp(common_properties::_CIPCode, cip);
    // TEST_ASSERT(cip=="R");

    delete m;
    delete[] ranks;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIssue3139534() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Issue 3139534: stereochemistry in larger rings"
                       << std::endl;

  // the smiles generation part of this is in SmilesParse/test.cpp

  // tests that the creation and assignment are correct:
  {
    RWMol *m;
    std::string smiles = "C1COCC/C=C\\CC1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(5)->getStereo() == Bond::STEREOZ);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C1COCC/C=C/CC1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(5)->getStereo() == Bond::STEREOE);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C/1=C/OCCC=CCC1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo() == Bond::STEREOZ);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C1=C/OCCC=CCC\\1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo() == Bond::STEREOZ);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C\\1=C/OCCC=CCC1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo() == Bond::STEREOE);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C1=C/OCCC=CCC/1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo() == Bond::STEREOE);
    delete m;
  }

  {
    RWMol *m;
    std::string smiles = "C/1=C/OCC/C=C\\CC1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo() == Bond::STEREOZ);
    TEST_ASSERT(m->getBondWithIdx(5)->getStereo() == Bond::STEREOZ);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C\\1=C/OCC/C=C\\CC1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo() == Bond::STEREOE);
    TEST_ASSERT(m->getBondWithIdx(5)->getStereo() == Bond::STEREOZ);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C1=C/OCC/C=C\\CC\\1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo() == Bond::STEREOZ);
    TEST_ASSERT(m->getBondWithIdx(5)->getStereo() == Bond::STEREOZ);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C1=C/OCC/C=C\\CC/1";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(0)->getStereo() == Bond::STEREOE);
    TEST_ASSERT(m->getBondWithIdx(5)->getStereo() == Bond::STEREOZ);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testFindChiralAtoms() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Test findChiralAtoms." << std::endl;

  {
    // by default the chirality possible flag is not assigned:
    RWMol *m;
    std::string smiles = "F[C@H](Cl)C(Cl)(Br)C(F)(F)F";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    MolOps::assignStereochemistry(*m, true, true);
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    TEST_ASSERT(!(m->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode)));
    TEST_ASSERT(!(m->getAtomWithIdx(6)->hasProp(common_properties::_CIPCode)));
    TEST_ASSERT(
        !m->getAtomWithIdx(1)->hasProp(common_properties::_ChiralityPossible));
    TEST_ASSERT(
        !m->getAtomWithIdx(3)->hasProp(common_properties::_ChiralityPossible));
    TEST_ASSERT(!(
        m->getAtomWithIdx(6)->hasProp(common_properties::_ChiralityPossible)));

    // but we can force it:
    MolOps::assignStereochemistry(*m, true, true, true);
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    TEST_ASSERT(!(m->getAtomWithIdx(3)->hasProp(common_properties::_CIPCode)));
    TEST_ASSERT(!(m->getAtomWithIdx(6)->hasProp(common_properties::_CIPCode)));
    TEST_ASSERT(
        m->getAtomWithIdx(1)->hasProp(common_properties::_ChiralityPossible));
    TEST_ASSERT(
        m->getAtomWithIdx(3)->hasProp(common_properties::_ChiralityPossible));
    TEST_ASSERT(!(
        m->getAtomWithIdx(6)->hasProp(common_properties::_ChiralityPossible)));

    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIssue3453172() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Issue 3453172: stereochemistry at three-coordinate S and Se"
      << std::endl;

  {
    RWMol *m;
    std::string smiles = "C=[S@](F)Br";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C[S@+](F)Br";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C=[Se@](F)Br";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C[Se@+](F)Br";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C=[S@](Br)Br";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C[S@+](Br)Br";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C=[Se@](Br)Br";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    delete m;
  }
  {
    RWMol *m;
    std::string smiles = "C[Se@+](Br)Br";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    delete m;
  }

  {
    // this was issue 254
    RWMol *m;
    std::string smiles = "O=[S@](c1ccccc1)C";
    m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getHybridization() == Atom::SP3);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub87() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing github issue 87: removal of bond wedging"
                       << std::endl;

  std::string rdbase = getenv("RDBASE");
  {
    std::string fName = rdbase + "/Code/GraphMol/test_data/github87.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 5);
    TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    WedgeMolBonds(*m, &m->getConformer());
    MolOps::assignStereochemistry(*m, true, true);
    TEST_ASSERT(m->getBondBetweenAtoms(0, 1)->getBondDir() == Bond::BEGINWEDGE);
    m->getAtomWithIdx(0)->setChiralTag(Atom::CHI_UNSPECIFIED);
    MolOps::assignStereochemistry(*m, true, true);
    TEST_ASSERT(m->getBondBetweenAtoms(0, 1)->getBondDir() == Bond::NONE);

    delete m;
  }
  {
    std::string fName = rdbase + "/Code/GraphMol/test_data/github87.2.mol";
    RWMol *m = MolFileToMol(fName);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 5);
    TEST_ASSERT(m->getAtomWithIdx(0)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    WedgeMolBonds(*m, &m->getConformer());
    MolOps::assignStereochemistry(*m, true, true);
    TEST_ASSERT(m->getBondBetweenAtoms(0, 1)->getBondDir() == Bond::BEGINDASH);
    m->getAtomWithIdx(0)->setChiralTag(Atom::CHI_UNSPECIFIED);
    MolOps::assignStereochemistry(*m, true, true);
    TEST_ASSERT(m->getBondBetweenAtoms(0, 1)->getBondDir() == Bond::NONE);

    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub90() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing github issue 90: isotopes and chirality"
                       << std::endl;

  {
    std::string smi = "C[C@@H](F)[13CH3]";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 4);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    MolOps::assignStereochemistry(*m, true, true);

    std::string cip;
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");
    delete m;
  }
  {
    std::string smi = "[13CH3][C@@H](F)C";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 4);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    MolOps::assignStereochemistry(*m, true, true);

    std::string cip;
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    delete m;
  }
  {
    std::string smi = "[CH3][C@@H](F)C";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 4);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    MolOps::assignStereochemistry(*m, true, true);
    TEST_ASSERT(!m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    delete m;
  }
  {
    std::string smi = "C\\C([13CH3])=C(/C)[13CH3]";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 6);
    TEST_ASSERT(m->getBondWithIdx(2)->getStereo() == Bond::STEREOZ);
    delete m;
  }
  {
    std::string smi = "C\\C([CH3])=C(/C)[13CH3]";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 6);
    TEST_ASSERT(m->getBondWithIdx(2)->getStereo() == Bond::STEREONONE);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub553() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing github issue 553: Chirality not affected by atom-map index"
      << std::endl;

  {
    std::string smi = "[*:1][C@H]([*:2])[*:3]";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 4);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    MolOps::assignStereochemistry(*m, true, true);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);

    std::string cip;
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");
    delete m;
  }

  {
    std::string smi = "*[C@H]([*:2])[*:3]";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 4);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    MolOps::assignStereochemistry(*m, true, true);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);

    std::string cip;
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");
    delete m;
  }

  {
    std::string smi = "[*:1][C@@H]([*:2])[*:3]";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 4);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    MolOps::assignStereochemistry(*m, true, true);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);

    std::string cip;
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub803() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing github issue 803: Support larger isotope "
                          "deltas in the chirality assignment"
                       << std::endl;

  {
    std::string smi = "*[C@H]([9*])[8*]";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 4);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    MolOps::assignStereochemistry(*m, true, true);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);

    std::string cip;
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    smi = MolToSmiles(*m, true);
    TEST_ASSERT(smi == "*[C@@H]([8*])[9*]");

    delete m;
  }
  {
    std::string smi = "*[C@H]([15*])[9*]";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 4);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    MolOps::assignStereochemistry(*m, true, true);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);

    std::string cip;
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "S");

    smi = MolToSmiles(*m, true);
    TEST_ASSERT(smi == "*[C@@H]([9*])[15*]");

    delete m;
  }
  {
    std::string smi = "[100U][C@H]([101U])[102U]";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 4);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    MolOps::assignStereochemistry(*m, true, true);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);

    std::string cip;
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    m->getAtomWithIdx(1)->getProp(common_properties::_CIPCode, cip);
    TEST_ASSERT(cip == "R");

    smi = MolToSmiles(*m, true);
    TEST_ASSERT(smi == "[100U][C@H]([101U])[102U]");

    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub1294() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing github issue 1294: ring stereochemistry "
                          "perception failing for spiro centers"
                       << std::endl;

  {  // the original example from the bug report
    std::string smi = "O[C@H]1CC[C@]11CC[C@@](Cl)(Br)CC1";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 12);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    TEST_ASSERT(m->getAtomWithIdx(4)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    TEST_ASSERT(m->getAtomWithIdx(7)->getChiralTag() != Atom::CHI_UNSPECIFIED);

    delete m;
  }

  {  // not spiro, but affected by same bug
    std::string smi = "C[C@H]1CC2CCCC3CCCC(C1)[C@@H]23";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 14);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    TEST_ASSERT(m->getAtomWithIdx(13)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    delete m;
  }

  {
    std::string smi = "C[C@H]1CC[C@@]2(CC[C@H](F)CC2)OC1";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 13);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    TEST_ASSERT(m->getAtomWithIdx(4)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    TEST_ASSERT(m->getAtomWithIdx(7)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub1423() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing github issue 1423: Generate a warning for "
                          "conflicting bond directions"
                       << std::endl;

  {  // this one is ok:
    std::stringstream warns;
    rdWarningLog->SetTee(warns);
    std::string smi = "C/C(/F)=C/C";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 5);
    TEST_ASSERT(m->getBondWithIdx(2)->getStereo() == Bond::STEREOZ);
    TEST_ASSERT(warns.str() == "");
    delete m;
    rdWarningLog->ClearTee();
  }

  {  // this one has a conflict:
    std::stringstream warns;
    rdWarningLog->SetTee(warns);
    std::string smi = "C/C(\\F)=C/C";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 5);
    TEST_ASSERT(m->getBondWithIdx(2)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(2)->getStereo() == Bond::STEREONONE);
    TEST_ASSERT(m->getBondWithIdx(0)->getBondType() == Bond::SINGLE);
    TEST_ASSERT(m->getBondWithIdx(0)->getBondDir() == Bond::NONE);
    TEST_ASSERT(m->getBondWithIdx(1)->getBondType() == Bond::SINGLE);
    TEST_ASSERT(m->getBondWithIdx(1)->getBondDir() == Bond::NONE);

    TEST_ASSERT(warns.str() != "");
    TEST_ASSERT(warns.str().find("BondStereo set to STEREONONE") !=
                std::string::npos);
    delete m;
    rdWarningLog->ClearTee();
  }
  {  // from the question that prompted this
    std::stringstream warns;
    rdWarningLog->SetTee(warns);
    std::string smi = "CCCO\\C(=C/c1ccccc1)/C(\\OCC)=C\\c1ccccc1";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(4)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(4)->getStereo() == Bond::STEREONONE);
    TEST_ASSERT(m->getBondWithIdx(15)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(15)->getStereo() == Bond::STEREONONE);
    TEST_ASSERT(warns.str() != "");
    TEST_ASSERT(warns.str().find("BondStereo set to STEREONONE") !=
                std::string::npos);
    delete m;
    rdWarningLog->ClearTee();
  }

  {  // a problem that came up during testing
    std::stringstream warns;
    rdWarningLog->SetTee(warns);
    std::string smi = "C/C(\\F)=C/[C@H](F)C=C(F)C";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondWithIdx(2)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondWithIdx(2)->getStereo() == Bond::STEREONONE);
    TEST_ASSERT(m->getAtomWithIdx(4)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    TEST_ASSERT(warns.str() != "");
    TEST_ASSERT(warns.str().find("BondStereo set to STEREONONE") !=
                std::string::npos);
    delete m;
    rdWarningLog->ClearTee();
  }

  {  // a problem that came up during testing
    std::stringstream warns;
    rdWarningLog->SetTee(warns);
    std::string smi = "C/C1=C/C=C=C=C2C(=C([Si](C)(C)C)\\C=C/1)C(=O)c1ccccc12";
    ROMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    TEST_ASSERT(warns.str() != "");
    TEST_ASSERT(warns.str().find("BondStereo set to STEREONONE") !=
                std::string::npos);
    delete m;
    rdWarningLog->ClearTee();
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

namespace {
void stereochemTester(RWMol *m, std::string expectedCIP,
                      Bond::BondStereo expectedStereo) {
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumAtoms() == 9)
  MolOps::sanitizeMol(*m);
  TEST_ASSERT(!m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  TEST_ASSERT(m->getBondWithIdx(3)->getStereo() == Bond::STEREONONE);
  // the mol file parser assigned bond dirs, get rid of them
  for (ROMol::BondIterator bIt = m->beginBonds(); bIt != m->endBonds(); ++bIt) {
    (*bIt)->setBondDir(Bond::NONE);
  }
  MolOps::assignStereochemistryFrom3D(*m);
  TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
  TEST_ASSERT(m->getAtomWithIdx(1)->getProp<std::string>(
                  common_properties::_CIPCode) == expectedCIP);
  TEST_ASSERT(m->getBondWithIdx(3)->getStereo() == expectedStereo);
}
}  // namespace
void testAssignStereochemistryFrom3D() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing assignStereochemistryFrom3D" << std::endl;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/test_data/";
  {
    SDMolSupplier suppl(pathName + "stereochem.sdf", false);  // don't sanitize
    {
      RWMol *m = (RWMol *)suppl.next();
      TEST_ASSERT(m->getProp<std::string>(common_properties::_Name) == "R-Z");
      stereochemTester(m, "R", Bond::STEREOZ);
      delete m;
    }
    {
      RWMol *m = (RWMol *)suppl.next();
      TEST_ASSERT(m->getProp<std::string>(common_properties::_Name) == "R-E");
      stereochemTester(m, "R", Bond::STEREOE);
      delete m;
    }
    {
      RWMol *m = (RWMol *)suppl.next();
      TEST_ASSERT(m->getProp<std::string>(common_properties::_Name) == "S-Z");
      stereochemTester(m, "S", Bond::STEREOZ);
      delete m;
    }
    {
      RWMol *m = (RWMol *)suppl.next();
      TEST_ASSERT(m->getProp<std::string>(common_properties::_Name) == "S-E");
      stereochemTester(m, "S", Bond::STEREOE);
      delete m;
    }
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testDoubleBondStereoInRings() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing double bond stereochemistry in rings"
                       << std::endl;
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/test_data/";
  {
    RWMol *m = MolFileToMol(pathName + "cyclohexene_3D.mol", true,
                            false);  // don't remove Hs
    MolOps::detectBondStereochemistry(*m);
    MolOps::assignChiralTypesFrom3D(*m);
    MolOps::assignStereochemistry(*m, true, true);
    const Bond *b = m->getBondBetweenAtoms(0, 1);
    TEST_ASSERT(b);
    TEST_ASSERT(b->getStereo() == Bond::STEREONONE);
    delete m;
  }
  {
    RWMol *m = MolFileToMol(pathName + "CHEMBL501674_3D.mol", true,
                            false);  // don't remove Hs
    MolOps::detectBondStereochemistry(*m);
    MolOps::assignChiralTypesFrom3D(*m);
    MolOps::assignStereochemistry(*m, true, true);
    const Bond *b = m->getBondBetweenAtoms(6, 7);
    TEST_ASSERT(b);
    TEST_ASSERT(b->getStereo() == Bond::STEREOE);
    b = m->getBondBetweenAtoms(11, 12);
    TEST_ASSERT(b);
    TEST_ASSERT(b->getStereo() == Bond::STEREOE);
    delete m;
  }
  {
    RWMol *m = MolFileToMol(pathName + "CHEMBL501674.mol");  // don't remove Hs
    const Bond *b = m->getBondBetweenAtoms(6, 7);
    TEST_ASSERT(b);
    TEST_ASSERT(b->getStereo() == Bond::STEREOE);
    b = m->getBondBetweenAtoms(11, 12);
    TEST_ASSERT(b);
    TEST_ASSERT(b->getStereo() == Bond::STEREOE);
    delete m;
  }
  {
    RWMol *m = MolFileToMol(pathName + "CHEMBL501674.mol");  // don't remove Hs
    std::string smi = MolToSmiles(*m, true);
    unsigned int nTrans = 0;
    size_t pos = 0;
    size_t len = smi.length();
    bool keepCounting = true;
    while (keepCounting) {
      pos = smi.find("/C=C/", pos);
      keepCounting = (pos != std::string::npos && ++nTrans && ++pos < len);
    }
    TEST_ASSERT(nTrans == 2);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testIssue1735() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing empty mols in rankMolAtoms" << std::endl;

  RWMol m;
  std::vector<unsigned int> oranks;
  Canon::rankMolAtoms(m, oranks);
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testStereoGroupUpdating() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------------------------------------------"
      << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Are stereo groups updated when atoms and bonds are deleted?"
      << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase + "/Code/GraphMol/FileParsers/test_data/two_centers_or.mol";
  std::unique_ptr<RWMol> m(MolFileToMol(fName));
  TEST_ASSERT(m.get());

  TEST_ASSERT(m->getStereoGroups().size() == 2);
  m->removeAtom(3);
  TEST_ASSERT(m->getStereoGroups().size() == 1);
  m->removeAtom(m->getAtomWithIdx(0));
  TEST_ASSERT(m->getStereoGroups().size() == 0u);
}

int main() {
  RDLog::InitLogs();
  // boost::logging::enable_logs("rdApp.debug");

#if 1
  testSmiles1();
  testMol1();
  testMol2();
  testRoundTrip();
  testChiralityCleanup();
  testChiralityFrom3D();
  testIterativeChirality();
  testBondDirRemoval();
  testIssue2762917();
  testIssue3009911();
  testIssue3139534();
  testFindChiralAtoms();
  testIssue3453172();
  testRingStereochemistry();
  testGithub87();
  testGithub90();
  testIssue2705543();
  testGithub553();
  testGithub803();
  testGithub1294();
#endif
  testGithub1423();
  testAssignStereochemistryFrom3D();
  testDoubleBondStereoInRings();
  testIssue1735();
  testStereoGroupUpdating();
  return 0;
}
