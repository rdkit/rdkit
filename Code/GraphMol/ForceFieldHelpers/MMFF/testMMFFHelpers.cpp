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
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>

#include <GraphMol/ForceFieldHelpers/FFConvenience.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <ForceField/ForceField.h>
#include <ForceField/MMFF/Params.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/Substruct/SubstructMatch.h>

using namespace RDKit;
void testMMFFTyper1() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test MMFF atom types." << std::endl;

  {
    std::uint8_t type;
    ROMol *mol = SmilesToMol("[SiH3]CC(=O)NC");
    TEST_ASSERT(mol);
    MMFF::MMFFMolProperties mmffMolProperties(*mol);
    TEST_ASSERT(mmffMolProperties.isValid());

    type = mmffMolProperties.getMMFFAtomType(0);
    TEST_ASSERT(type == 19);
    type = mmffMolProperties.getMMFFAtomType(1);
    TEST_ASSERT(type == 1);
    type = mmffMolProperties.getMMFFAtomType(2);
    TEST_ASSERT(type == 3);
    type = mmffMolProperties.getMMFFAtomType(3);
    TEST_ASSERT(type == 7);
    type = mmffMolProperties.getMMFFAtomType(4);
    TEST_ASSERT(type == 10);
    delete mol;
  }

  {
    std::uint8_t type;
    ROMol *mol = SmilesToMol("CC(=O)C");
    TEST_ASSERT(mol);
    MMFF::MMFFMolProperties mmffMolProperties(*mol);
    TEST_ASSERT(mmffMolProperties.isValid());

    type = mmffMolProperties.getMMFFAtomType(0);
    TEST_ASSERT(type == 1);
    type = mmffMolProperties.getMMFFAtomType(1);
    TEST_ASSERT(type == 3);
    type = mmffMolProperties.getMMFFAtomType(2);
    TEST_ASSERT(type == 7);
    type = mmffMolProperties.getMMFFAtomType(3);
    TEST_ASSERT(type == 1);
    delete mol;
  }

  {
    std::uint8_t type;
    ROMol *mol = SmilesToMol("C(=O)S");
    TEST_ASSERT(mol);
    MMFF::MMFFMolProperties mmffMolProperties(*mol);
    TEST_ASSERT(mmffMolProperties.isValid());

    type = mmffMolProperties.getMMFFAtomType(0);
    TEST_ASSERT(type == 3);
    type = mmffMolProperties.getMMFFAtomType(1);
    TEST_ASSERT(type == 7);
    type = mmffMolProperties.getMMFFAtomType(2);
    TEST_ASSERT(type == 15);
    delete mol;
  }
  {
    std::uint8_t type;
    ROMol *mol = SmilesToMol("SCS(=O)S(=O)(=O)O");
    TEST_ASSERT(mol);
    MMFF::MMFFMolProperties mmffMolProperties(*mol);
    TEST_ASSERT(mmffMolProperties.isValid());

    type = mmffMolProperties.getMMFFAtomType(0);
    TEST_ASSERT(type == 15);
    type = mmffMolProperties.getMMFFAtomType(1);
    TEST_ASSERT(type == 1);
    type = mmffMolProperties.getMMFFAtomType(2);
    TEST_ASSERT(type == 17);
    type = mmffMolProperties.getMMFFAtomType(4);
    TEST_ASSERT(type == 18);
    delete mol;
  }
  {
    std::uint8_t type;
    ROMol *mol = SmilesToMol("PCP(O)CP(=O)(=O)");
    TEST_ASSERT(mol);
    MMFF::MMFFMolProperties mmffMolProperties(*mol);
    TEST_ASSERT(mmffMolProperties.isValid());

    type = mmffMolProperties.getMMFFAtomType(0);
    TEST_ASSERT(type == 26);
    type = mmffMolProperties.getMMFFAtomType(1);
    TEST_ASSERT(type == 1);
    type = mmffMolProperties.getMMFFAtomType(2);
    TEST_ASSERT(type == 26);
    type = mmffMolProperties.getMMFFAtomType(5);
    TEST_ASSERT(type == 26);
    delete mol;
  }
  {
    std::uint8_t type;
    ROMol *mol = SmilesToMol("C(F)(Cl)(Br)I");
    TEST_ASSERT(mol);
    MMFF::MMFFMolProperties mmffMolProperties(*mol);
    TEST_ASSERT(mmffMolProperties.isValid());

    type = mmffMolProperties.getMMFFAtomType(0);
    TEST_ASSERT(type == 1);
    type = mmffMolProperties.getMMFFAtomType(1);
    TEST_ASSERT(type == 11);
    type = mmffMolProperties.getMMFFAtomType(2);
    TEST_ASSERT(type == 12);
    type = mmffMolProperties.getMMFFAtomType(3);
    TEST_ASSERT(type == 13);
    type = mmffMolProperties.getMMFFAtomType(4);
    TEST_ASSERT(type == 14);
    delete mol;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMMFFBuilder1() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing MMFF builder tools." << std::endl;

  ROMol *mol, *mol2;

  ForceFields::ForceField *field;
  boost::shared_array<std::uint8_t> nbrMat;

  mol = SmilesToMol("CC(O)C");
  auto *conf = new Conformer(mol->getNumAtoms());
  int cid = static_cast<int>(mol->addConformer(conf, true));
  TEST_ASSERT(mol);
  MMFF::MMFFMolProperties *mmffMolProperties =
      new MMFF::MMFFMolProperties(*mol);
  TEST_ASSERT(mmffMolProperties);
  TEST_ASSERT(mmffMolProperties->isValid());
  field = new ForceFields::ForceField();
  // add the atomic positions:
  for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
    field->positions().push_back(&(conf->getAtomPos(i)));
  }

  MMFF::Tools::addBonds(*mol, mmffMolProperties, field);

  TEST_ASSERT(field->contribs().size() == 3);

  nbrMat = MMFF::Tools::buildNeighborMatrix(*mol);
  // the neighbor matrix is an upper triangular matrix
  // position indices are as follows:
  //  0  1  2  3
  //     4  5  6
  //        7  8
  //           9
  TEST_ASSERT(MMFF::Tools::twoBitCellPos(mol->getNumAtoms(), 1, 1) == 4);
  TEST_ASSERT(MMFF::Tools::twoBitCellPos(mol->getNumAtoms(), 2, 1) == 5);
  TEST_ASSERT(MMFF::Tools::twoBitCellPos(mol->getNumAtoms(), 1, 2) == 5);
  TEST_ASSERT(MMFF::Tools::getTwoBitCell(nbrMat, 0) ==
              MMFF::Tools::RELATION_1_X);
  TEST_ASSERT(MMFF::Tools::getTwoBitCell(nbrMat, 1) ==
              MMFF::Tools::RELATION_1_2);
  TEST_ASSERT(MMFF::Tools::getTwoBitCell(nbrMat, 2) ==
              MMFF::Tools::RELATION_1_3);
  TEST_ASSERT(MMFF::Tools::getTwoBitCell(nbrMat, 3) ==
              MMFF::Tools::RELATION_1_3);

  MMFF::Tools::addAngles(*mol, mmffMolProperties, field);
  TEST_ASSERT(field->contribs().size() == 6);

  // there are no non-bonded terms here:
  MMFF::Tools::addVdW(*mol, cid, mmffMolProperties, field, nbrMat);
  TEST_ASSERT(field->contribs().size() == 6);

  // and no torsions either, until we add Hs:
  MMFF::Tools::addTorsions(*mol, mmffMolProperties, field);
  TEST_ASSERT(field->contribs().size() == 6);

  delete mol;
  delete field;
  delete mmffMolProperties;
  mol = SmilesToMol("CCOC");
  auto *conf2 = new Conformer(mol->getNumAtoms());
  cid = static_cast<int>(mol->addConformer(conf2, true));
  TEST_ASSERT(mol);
  mmffMolProperties = new MMFF::MMFFMolProperties(*mol);
  TEST_ASSERT(mmffMolProperties);
  TEST_ASSERT(mmffMolProperties->isValid());
  field = new ForceFields::ForceField();
  // add the atomic positions:
  for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
    field->positions().push_back(&(conf2->getAtomPos(i)));
  }

  MMFF::Tools::addBonds(*mol, mmffMolProperties, field);

  TEST_ASSERT(field->contribs().size() == 3);

  nbrMat = MMFF::Tools::buildNeighborMatrix(*mol);
  TEST_ASSERT(MMFF::Tools::getTwoBitCell(nbrMat, 0) ==
              MMFF::Tools::RELATION_1_X);
  TEST_ASSERT(MMFF::Tools::getTwoBitCell(nbrMat, 1) ==
              MMFF::Tools::RELATION_1_2);
  TEST_ASSERT(MMFF::Tools::getTwoBitCell(nbrMat, 2) ==
              MMFF::Tools::RELATION_1_3);
  TEST_ASSERT(MMFF::Tools::getTwoBitCell(nbrMat, 3) ==
              MMFF::Tools::RELATION_1_4);

  MMFF::Tools::addAngles(*mol, mmffMolProperties, field);
  TEST_ASSERT(field->contribs().size() == 5);
  MMFF::Tools::addVdW(*mol, cid, mmffMolProperties, field, nbrMat);
  TEST_ASSERT(field->contribs().size() == 6);
  MMFF::Tools::addTorsions(*mol, mmffMolProperties, field);
  TEST_ASSERT(field->contribs().size() == 7);

  delete mol;
  delete field;
  delete mmffMolProperties;
  mol = SmilesToMol("CO");
  auto *conf3 = new Conformer(mol->getNumAtoms());
  cid = static_cast<int>(mol->addConformer(conf3, true));
  TEST_ASSERT(mol);
  mmffMolProperties = new MMFF::MMFFMolProperties(*mol);
  TEST_ASSERT(mmffMolProperties);
  TEST_ASSERT(mmffMolProperties->isValid());

  field = new ForceFields::ForceField();
  // add the atomic positions:
  for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
    field->positions().push_back(&(conf3->getAtomPos(i)));
  }

  MMFF::Tools::addBonds(*mol, mmffMolProperties, field);

  TEST_ASSERT(field->contribs().size() == 1);

  nbrMat = MMFF::Tools::buildNeighborMatrix(*mol);
  TEST_ASSERT(MMFF::Tools::getTwoBitCell(nbrMat, 0) ==
              MMFF::Tools::RELATION_1_X);
  TEST_ASSERT(MMFF::Tools::getTwoBitCell(nbrMat, 1) ==
              MMFF::Tools::RELATION_1_2);

  MMFF::Tools::addAngles(*mol, mmffMolProperties, field);
  TEST_ASSERT(field->contribs().size() == 1);
  MMFF::Tools::addVdW(*mol, cid, mmffMolProperties, field, nbrMat);
  TEST_ASSERT(field->contribs().size() == 1);
  MMFF::Tools::addTorsions(*mol, mmffMolProperties, field);
  TEST_ASSERT(field->contribs().size() == 1);

  mol2 = MolOps::addHs(*mol);
  TEST_ASSERT(mol2->getNumAtoms() == 6);
  auto *conf4 = new Conformer(mol2->getNumAtoms());
  cid = static_cast<int>(mol2->addConformer(conf4, true));

  delete field;
  delete mmffMolProperties;

  mmffMolProperties = new MMFF::MMFFMolProperties(*mol2);
  TEST_ASSERT(mmffMolProperties);
  TEST_ASSERT(mmffMolProperties->isValid());

  field = new ForceFields::ForceField();
  // add the atomic positions:
  for (unsigned int i = 0; i < mol2->getNumAtoms(); ++i) {
    field->positions().push_back(&(conf4->getAtomPos(i)));
  }

  MMFF::Tools::addBonds(*mol2, mmffMolProperties, field);
  TEST_ASSERT(field->contribs().size() == 5);

  nbrMat = MMFF::Tools::buildNeighborMatrix(*mol2);
  MMFF::Tools::addAngles(*mol2, mmffMolProperties, field);
  TEST_ASSERT(field->contribs().size() == 12);
  MMFF::Tools::addVdW(*mol2, cid, mmffMolProperties, field, nbrMat);
  TEST_ASSERT(field->contribs().size() == 15);
  MMFF::Tools::addTorsions(*mol2, mmffMolProperties, field);
  TEST_ASSERT(field->contribs().size() == 18);
  delete mol2;

  delete mol;
  delete field;
  delete mmffMolProperties;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMMFFBuilder2() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing MMFF builder+minimization."
                        << std::endl;

  RWMol *mol;
  int needMore;

  ForceFields::ForceField *field;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/MMFF/test_data";
  mol = MolFileToMol(pathName + "/small1.mol", false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);

  field = MMFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  needMore = field->minimize();
  TEST_ASSERT(!needMore);
  // std::cout << MolToMolBlock(mol);
  delete mol;
  delete field;

  mol = MolFileToMol(pathName + "/small2.mol", false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);

  field = MMFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  needMore = field->minimize(100);
  TEST_ASSERT(!needMore);
  // std::cout << MolToMolBlock(mol);
  delete mol;
  delete field;

  mol = MolFileToMol(pathName + "/benzene.mol", false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);

  field = MMFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  needMore = field->minimize();
  TEST_ASSERT(!needMore);
  // std::cout << MolToMolBlock(mol);
  delete mol;
  delete field;

  mol = MolFileToMol(pathName + "/toluene.mol", false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);

  field = MMFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  needMore = field->minimize();
  TEST_ASSERT(!needMore);
  // std::cout << MolToMolBlock(mol);
  delete mol;
  delete field;

  mol = MolFileToMol(pathName + "/complex1.mol", false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);

  field = MMFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  needMore = field->minimize();
  TEST_ASSERT(!needMore);
  // std::cout << MolToMolBlock(mol);
  delete mol;
  delete field;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

#ifdef RDK_TEST_MULTITHREADED
// we do the equivalent tests below
void testMMFFBatch() {}
#else
void testMMFFBatch() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "    Testing bulk MMFF (regression to check that things run)."
      << std::endl;

  ROMol *mol;

  ForceFields::ForceField *field;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/MMFF/test_data";
  SDMolSupplier suppl(pathName + "/bulk.sdf");

  int count = 0;
  mol = suppl.next();
  while (mol && (!(suppl.atEnd()))) {
    ++count;
    std::string origMolBlock = MolToMolBlock(*mol);

    BOOST_LOG(rdErrorLog) << "Mol:" << count << std::endl;
    try {
      field = MMFF::constructForceField(*mol);
    } catch (...) {
      field = nullptr;
    }
    if (field) {
      field->initialize();
      int failed = field->minimize(500);
      delete field;
      if (failed) {
        BOOST_LOG(rdErrorLog)
            << " not converged (code = " << failed << ")" << std::endl;
        std::cout << origMolBlock << "$$$$" << std::endl;
        std::cout << MolToMolBlock(*mol) << "$$$$" << std::endl;
      }
    }
    delete mol;
    mol = suppl.next();
  }
  delete mol;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}
#endif
void testIssue239() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing Issue239." << std::endl;

  RWMol *mol;
  int needMore;
  (void)needMore;  // Add test later
  ForceFields::ForceField *field;
  double e1, e2;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/MMFF/test_data";
  mol = MolFileToMol(pathName + "/Issue239.mol", false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);

  field = MMFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  needMore = field->minimize(200, 1.0e-6, 1.0e-3);
  e1 = field->calcEnergy();
  needMore = field->minimize(200, 1.0e-6, 1.0e-3);
  e2 = field->calcEnergy();
  TEST_ASSERT(feq(e2, e1, 0.1));

  delete mol;
  delete field;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testCalcEnergyPassedCoords() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing calcEnergy with passed coords."
                        << std::endl;

  RWMol *mol;
  ForceFields::ForceField *field;
  double e1, e2, e3;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/MMFF/test_data";
  mol = MolFileToMol(pathName + "/Issue239.mol", false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);

  field = MMFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  size_t l = 3 * field->numPoints();
  auto *savedPos = new double[l];
  size_t i = 0;
  for (const auto pptr : field->positions()) {
    for (size_t j = 0; j < 3; ++j) {
      savedPos[i++] = (*pptr)[j];
    }
  }
  TEST_ASSERT(i == l);
  e1 = field->calcEnergy();
  field->minimize(10000, 1.0e-6, 1.0e-3);
  e2 = field->calcEnergy();
  TEST_ASSERT(e2 < e1);
  e3 = field->calcEnergy(savedPos);
  TEST_ASSERT(feq(e3, e1, 0.01));

  delete[] savedPos;
  delete mol;
  delete field;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testCalcGrad() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing calcGrad." << std::endl;

  RWMol *mol;
  ForceFields::ForceField *field;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/MMFF/test_data";
  mol = MolFileToMol(pathName + "/Issue239.mol", false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);

  field = MMFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  size_t l = 3 * field->numPoints();
  auto *savedPos = new double[l];
  auto *grad1 = new double[l];
  auto *grad2 = new double[l];
  size_t i = 0;
  for (const auto pptr : field->positions()) {
    for (size_t j = 0; j < 3; ++j) {
      savedPos[i++] = (*pptr)[j];
    }
  }
  TEST_ASSERT(i == l);

  std::memset(grad1, 0, l * sizeof(double));
  field->calcGrad(grad1);
  for (i = 0; i < l; ++i) {
    TEST_ASSERT(!feq(grad1[i], 0.0, 0.001));
  }

  field->minimize(10000, 1.0e-6, 1.0e-3);
  std::memset(grad2, 0, l * sizeof(double));
  field->calcGrad(grad2);
  for (i = 0; i < l; ++i) {
    TEST_ASSERT(feq(grad2[i], 0.0, 0.001));
  }

  field->initialize();
  std::memset(grad2, 0, l * sizeof(double));
  field->calcGrad(savedPos, grad2);
  for (i = 0; i < l; ++i) {
    TEST_ASSERT(feq(grad1[i], grad2[i], 0.001));
  }

  delete[] savedPos;
  delete[] grad1;
  delete[] grad2;
  delete mol;
  delete field;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testIssue242() {
  // FIX: Changes to the forcefield (connected to Issue 408) have
  // made it so that this particular problem no longer manifests
  // in this molecule/embedding. A new test case is needed.
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing Issue242." << std::endl;

  RWMol *mol, *mol2;
  int needMore;
  (void)needMore;  // add test later
  ForceFields::ForceField *field = nullptr, *field2 = nullptr;
  std::string mb1, mb2;
  double e1, e2;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/MMFF/test_data";

  mol = MolFileToMol(pathName + "/Issue242-2.mol");
  TEST_ASSERT(mol);

  mol2 = MolFileToMol(pathName + "/Issue242-2.mol");
  TEST_ASSERT(mol2);

  TEST_ASSERT(DGeomHelpers::EmbedMolecule(*mol2, 30, 2370) >= 0);
  mb1 = MolToMolBlock(*mol);
  mb2 = MolToMolBlock(*mol2);

  // std::cout << mb1 << "------\n";
  // std::cout << mb2 << "------\n";

  field = MMFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  field2 = MMFF::constructForceField(*mol2);
  TEST_ASSERT(field2);
  field2->initialize();
  e1 = field->calcEnergy();
  e2 = field2->calcEnergy();
  BOOST_LOG(rdInfoLog) << "E1: " << e1 << std::endl;
  BOOST_LOG(rdInfoLog) << "E2: " << e2 << std::endl;
  // TEST_ASSERT(feq(e2,e1,0.1));

  needMore = field->minimize(200, 1.0e-4);
  needMore = field2->minimize(200, 1.0e-4);
  e1 = field->calcEnergy();
  e2 = field2->calcEnergy();
  BOOST_LOG(rdInfoLog) << "E1: " << e1 << std::endl;
  BOOST_LOG(rdInfoLog) << "E2: " << e2 << std::endl;
  TEST_ASSERT(feq(e2, e1, 0.1));

  delete mol;
  delete mol2;
  delete field;
  delete field2;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testGithub308() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "    Testing github 308: crash during MMFF parameterization ."
      << std::endl;
  ROMol *mol = SmilesToMol("FF");
  TEST_ASSERT(DGeomHelpers::EmbedMolecule(*mol) >= 0);
  int needMore;
  ForceFields::ForceField *field = nullptr;
  TEST_ASSERT(mol);
  MMFF::MMFFMolProperties mmffMolProperties(*mol);
  TEST_ASSERT(mmffMolProperties.isValid());
  field = MMFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  needMore = field->minimize(200, 1.0e-6, 1.0e-3);
  TEST_ASSERT(!needMore);
  delete mol;
  delete field;
}

void testSFIssue1653802() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing SFIssue1653802." << std::endl;

  RWMol *mol;
  int needMore;
  ForceFields::ForceField *field;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/MMFF/test_data";

  mol = MolFileToMol(pathName + "/cyclobutadiene.mol", false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);
  MMFF::MMFFMolProperties *mmffMolProperties =
      new MMFF::MMFFMolProperties(*mol);
  TEST_ASSERT(mmffMolProperties);

  boost::shared_array<std::uint8_t> nbrMat;
  field = new ForceFields::ForceField();
  // add the atomic positions:
  for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
    field->positions().push_back(&((mol->getConformer().getAtomPos(i))));
  }

  MMFF::Tools::addBonds(*mol, mmffMolProperties, field);
  TEST_ASSERT(field->contribs().size() == 8);

  nbrMat = MMFF::Tools::buildNeighborMatrix(*mol);
  MMFF::Tools::addAngles(*mol, mmffMolProperties, field);
  TEST_ASSERT(field->contribs().size() == 20);
  MMFF::Tools::addTorsions(*mol, mmffMolProperties, field);
  // std::cout << field->contribs().size() << std::endl;
  TEST_ASSERT(field->contribs().size() == 36);
  MMFF::Tools::addVdW(*mol, 0, mmffMolProperties, field, nbrMat);
  delete field;

  field = MMFF::constructForceField(*mol);
  field->initialize();
  needMore = field->minimize(200, 1.0e-6, 1.0e-3);
  TEST_ASSERT(!needMore);

  delete mol;
  delete field;
  delete mmffMolProperties;
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testSFIssue2378119() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing SFIssue2378119." << std::endl;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/MMFF/test_data";
  {
    RWMol *mol = MolFileToMol(pathName + "/Issue2378119.mol");
    TEST_ASSERT(mol);
    ForceFields::ForceField *field = MMFF::constructForceField(*mol);
    TEST_ASSERT(field);
    field->initialize();
    double e1 = field->calcEnergy();
    TEST_ASSERT(e1 > 0.0 && e1 < 1.0e12);

    int needMore = field->minimize(200, 1.0e-6, 1.0e-3);
    TEST_ASSERT(!needMore);
    double e2 = field->calcEnergy();
    TEST_ASSERT(e2 < e1);

    delete mol;
    delete field;
  }
  {
    RWMol *mol = MolFileToMol(pathName + "/Issue2378119.2.mol");
    TEST_ASSERT(mol);
    ForceFields::ForceField *field = MMFF::constructForceField(*mol);
    TEST_ASSERT(field);
    field->initialize();
    double e1 = field->calcEnergy();
    TEST_ASSERT(e1 == 0.0);

    int needMore = field->minimize(200, 1.0e-6, 1.0e-3);
    TEST_ASSERT(!needMore);
    double e2 = field->calcEnergy();
    TEST_ASSERT(e2 == e1);

    delete mol;
    delete field;
  }
  {
    RWMol *mol = MolFileToMol(pathName + "/Issue2378119.2.mol");
    TEST_ASSERT(mol);
    ForceFields::ForceField *field =
        MMFF::constructForceField(*mol, 100.0, -1, false);
    TEST_ASSERT(field);
    field->initialize();
    double e1 = field->calcEnergy();
    TEST_ASSERT(e1 > 0.0 && e1 < 1.0e12);

    int needMore = field->minimize(200, 1.0e-6, 1.0e-3);
    TEST_ASSERT(!needMore);
    double e2 = field->calcEnergy();
    TEST_ASSERT(e2 < e1);

    delete mol;
    delete field;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testGithub162() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing github 162: Incorrect SMILES after "
                           "MMFF parameterization ."
                        << std::endl;
  {
    ROMol *mol = SmilesToMol("C1=CNC=C1");
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getAtomWithIdx(2)->getNumExplicitHs() == 1);
    MMFF::MMFFMolProperties *mmffMolProperties =
        new MMFF::MMFFMolProperties(*mol);
    TEST_ASSERT(mmffMolProperties);
    TEST_ASSERT(mmffMolProperties->isValid());
    TEST_ASSERT(mol->getAtomWithIdx(2)->getNumExplicitHs() == 1);

    delete mol;
    delete mmffMolProperties;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMMFFParamGetters() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test MMFF force-field parameter getters."
                        << std::endl;
  {
    ROMol *mol = SmilesToMol("c1ccccc1CCNN");
    TEST_ASSERT(mol);
    ROMol *molH = MolOps::addHs(*mol);
    delete mol;
    TEST_ASSERT(molH);
    MMFF::MMFFMolProperties *mmffMolProperties =
        new MMFF::MMFFMolProperties(*molH);
    TEST_ASSERT(mmffMolProperties);
    TEST_ASSERT(mmffMolProperties->isValid());
    unsigned int bondType;
    ForceFields::MMFF::MMFFBond mmffBondStretchParams[2];
    TEST_ASSERT(mmffMolProperties->getMMFFBondStretchParams(
        *molH, 6, 7, bondType, mmffBondStretchParams[0]));
    TEST_ASSERT((bondType == 0) &&
                ((int)std::round(mmffBondStretchParams[0].r0 * 1000) == 1508) &&
                ((int)std::round(mmffBondStretchParams[0].kb * 1000) == 4258));
    TEST_ASSERT(!(mmffMolProperties->getMMFFBondStretchParams(
        *molH, 0, 7, bondType, mmffBondStretchParams[0])));
    unsigned int angleType;
    ForceFields::MMFF::MMFFAngle mmffAngleBendParams;
    TEST_ASSERT(mmffMolProperties->getMMFFAngleBendParams(
        *molH, 6, 7, 8, angleType, mmffAngleBendParams));
    TEST_ASSERT(
        (angleType == 0) &&
        ((int)std::round(mmffAngleBendParams.theta0 * 1000) == 108290) &&
        ((int)std::round(mmffAngleBendParams.ka * 1000) == 777));
    TEST_ASSERT(!(mmffMolProperties->getMMFFAngleBendParams(
        *molH, 0, 7, 8, angleType, mmffAngleBendParams)));
    unsigned int stretchBendType;
    ForceFields::MMFF::MMFFStbn mmffStretchBendParams;
    TEST_ASSERT(mmffMolProperties->getMMFFStretchBendParams(
        *molH, 6, 7, 8, stretchBendType, mmffStretchBendParams,
        mmffBondStretchParams, mmffAngleBendParams));
    TEST_ASSERT(
        (stretchBendType == 0) &&
        ((int)std::round(mmffStretchBendParams.kbaIJK * 1000) == 136) &&
        ((int)std::round(mmffStretchBendParams.kbaKJI * 1000) == 282) &&
        ((int)std::round(mmffAngleBendParams.theta0 * 1000) == 108290) &&
        ((int)std::round(mmffAngleBendParams.ka * 1000) == 777) &&
        ((int)std::round(mmffBondStretchParams[0].r0 * 1000) == 1508) &&
        ((int)std::round(mmffBondStretchParams[0].kb * 1000) == 4258) &&
        ((int)std::round(mmffBondStretchParams[1].r0 * 1000) == 1451) &&
        ((int)std::round(mmffBondStretchParams[1].kb * 1000) == 5084));
    TEST_ASSERT(!(mmffMolProperties->getMMFFStretchBendParams(
        *molH, 0, 7, 8, stretchBendType, mmffStretchBendParams,
        mmffBondStretchParams, mmffAngleBendParams)));
    unsigned int torType;
    ForceFields::MMFF::MMFFTor mmffTorsionParams;
    TEST_ASSERT(mmffMolProperties->getMMFFTorsionParams(
        *molH, 6, 7, 8, 9, torType, mmffTorsionParams));
    TEST_ASSERT((torType == 0) &&
                ((int)std::round(mmffTorsionParams.V1 * 1000) == 0) &&
                ((int)std::round(mmffTorsionParams.V2 * 1000) == -300) &&
                ((int)std::round(mmffTorsionParams.V3 * 1000) == 500));
    TEST_ASSERT(!(mmffMolProperties->getMMFFTorsionParams(
        *molH, 0, 7, 8, 9, torType, mmffTorsionParams)));
    ForceFields::MMFF::MMFFOop mmffOopBendParams;
    TEST_ASSERT(mmffMolProperties->getMMFFOopBendParams(*molH, 6, 5, 4, 0,
                                                        mmffOopBendParams));
    TEST_ASSERT(((int)std::round(mmffOopBendParams.koop * 1000) == 40));
    TEST_ASSERT(!(mmffMolProperties->getMMFFOopBendParams(*molH, 6, 5, 4, 1,
                                                          mmffOopBendParams)));
    ForceFields::MMFF::MMFFVdWRijstarEps mmffVdWParams;
    RWMol *patt = SmartsToMol("NN[H]");
    MatchVectType matchVect;
    TEST_ASSERT(SubstructMatch(*molH, (ROMol &)*patt, matchVect));
    delete patt;
    unsigned int nIdx = matchVect[0].second;
    unsigned int hIdx = matchVect[2].second;
    TEST_ASSERT(mmffMolProperties->getMMFFVdWParams(nIdx, hIdx, mmffVdWParams));
    TEST_ASSERT(
        ((int)std::round(mmffVdWParams.R_ij_starUnscaled * 1000) == 3321) &&
        ((int)std::round(mmffVdWParams.epsilonUnscaled * 1000) == 34) &&
        ((int)std::round(mmffVdWParams.R_ij_star * 1000) == 2657) &&
        ((int)std::round(mmffVdWParams.epsilon * 1000) == 17));
    delete molH;
    delete mmffMolProperties;
  }
}
#ifdef RDK_TEST_MULTITHREADED
namespace {
void runblock_mmff(const std::vector<ROMol *> &mols,
                   const std::vector<double> &energies, unsigned int count,
                   unsigned int idx) {
  for (unsigned int rep = 0; rep < 100; ++rep) {
    for (unsigned int i = 0; i < mols.size(); ++i) {
      if (i % count != idx) {
        continue;
      }
      ROMol *mol = mols[i];
      ForceFields::ForceField *field = nullptr;
      if (!(rep % 20)) {
        BOOST_LOG(rdErrorLog) << "Rep: " << rep << " Mol:" << i << std::endl;
      }
      try {
        field = MMFF::constructForceField(*mol);
      } catch (...) {
        field = nullptr;
      }
      TEST_ASSERT(field);
      field->initialize();
      int failed = field->minimize(500);
      TEST_ASSERT(!failed);
      double eng = field->calcEnergy();
      TEST_ASSERT(feq(eng, energies[i]));
      delete field;
    }
  }
}
}  // namespace
#include <thread>
#include <future>
void testMMFFMultiThread() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test MMFF multithreading" << std::endl;

  // ForceFields::ForceField *field;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/MMFF/test_data";
  SDMolSupplier suppl(pathName + "/bulk.sdf");
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

  std::cerr << "generating reference data" << std::endl;
  std::vector<double> energies(mols.size(), 0.0);
  for (unsigned int i = 0; i < mols.size(); ++i) {
    ROMol mol(*mols[i]);
    ForceFields::ForceField *field = nullptr;
    try {
      field = MMFF::constructForceField(mol);
    } catch (...) {
      field = nullptr;
    }
    TEST_ASSERT(field);
    field->initialize();
    int failed = field->minimize(500);
    TEST_ASSERT(!failed);
    energies[i] = field->calcEnergy();
    delete field;
  }

  std::vector<std::future<void>> tg;

  std::cerr << "processing" << std::endl;
  unsigned int count = 4;
  for (unsigned int i = 0; i < count; ++i) {
    std::cerr << " launch :" << i << std::endl;
    std::cerr.flush();
    tg.emplace_back(std::async(std::launch::async, runblock_mmff, mols,
                               energies, count, i));
  }
  for (auto &fut : tg) {
    fut.get();
  }

  for (auto *mol : mols) {
    delete mol;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMMFFMultiThread2() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test MMFF multithreading2" << std::endl;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/UFF/test_data";
  SDMolSupplier suppl(pathName + "/bulk.sdf");
  ROMol *m = suppl[4];
  TEST_ASSERT(m);
  auto *om = new ROMol(*m);
  for (unsigned int i = 0; i < 200; ++i) {
    m->addConformer(new Conformer(m->getConformer()), true);
  }
  std::vector<std::pair<int, double>> res;

  MMFF::MMFFOptimizeMolecule(*om);
  MMFF::MMFFOptimizeMoleculeConfs(*m, res, 0);
  for (unsigned int i = 1; i < res.size(); ++i) {
    TEST_ASSERT(!res[i].first);
    TEST_ASSERT(feq(res[i].second, res[0].second, .00001));
  }
  for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
    RDGeom::Point3D p0 = om->getConformer().getAtomPos(i);
    RDGeom::Point3D np0 = m->getConformer().getAtomPos(i);
    TEST_ASSERT(feq(p0.x, np0.x));
    TEST_ASSERT(feq(p0.y, np0.y));
    TEST_ASSERT(feq(p0.z, np0.z));
    np0 =
        m->getConformer(11).getAtomPos(i);  // pick some random other conformer
    TEST_ASSERT(feq(p0.x, np0.x));
    TEST_ASSERT(feq(p0.y, np0.y));
    TEST_ASSERT(feq(p0.z, np0.z));
  }
  delete m;
  delete om;
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMMFFMultiThread3() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test MMFF multithreading3" << std::endl;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/UFF/test_data";
  SDMolSupplier suppl(pathName + "/bulk.sdf");
  ROMol *m = suppl[4];
  TEST_ASSERT(m);
  auto *om = new ROMol(*m);
  for (unsigned int i = 0; i < 200; ++i) {
    m->addConformer(new Conformer(m->getConformer()), true);
  }
  std::vector<std::pair<int, double>> res;

  ForceFields::ForceField *omField = MMFF::constructForceField(*om);
  TEST_ASSERT(omField);
  omField->initialize();
  ForceFields::ForceField *mField = MMFF::constructForceField(*m);
  TEST_ASSERT(mField);
  mField->initialize();

  ForceFieldsHelper::OptimizeMolecule(*omField);
  ForceFieldsHelper::OptimizeMoleculeConfs(*m, *mField, res, 0);
  for (unsigned int i = 1; i < res.size(); ++i) {
    TEST_ASSERT(!res[i].first);
    TEST_ASSERT(feq(res[i].second, res[0].second, .00001));
  }
  for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
    RDGeom::Point3D p0 = om->getConformer().getAtomPos(i);
    RDGeom::Point3D np0 = m->getConformer().getAtomPos(i);
    TEST_ASSERT(feq(p0.x, np0.x));
    TEST_ASSERT(feq(p0.y, np0.y));
    TEST_ASSERT(feq(p0.z, np0.z));
    np0 =
        m->getConformer(11).getAtomPos(i);  // pick some random other conformer
    TEST_ASSERT(feq(p0.x, np0.x));
    TEST_ASSERT(feq(p0.y, np0.y));
    TEST_ASSERT(feq(p0.z, np0.z));
  }
  delete m;
  delete om;
  delete mField;
  delete omField;
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}
#endif

void testGithub224() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "    Testing github 224: crash during MMFF parameterization ."
      << std::endl;
  {
    ROMol *mol = SmilesToMol("[1*]C");
    TEST_ASSERT(mol);
    MMFF::MMFFMolProperties *mmffMolProperties =
        new MMFF::MMFFMolProperties(*mol);
    TEST_ASSERT(mmffMolProperties);
    TEST_ASSERT(!mmffMolProperties->isValid());

    delete mol;
    delete mmffMolProperties;
  }
  {
    ROMol *mol = SmilesToMol("[1*]c1ccc(S(=O)(=O)Nc2ccc([2*])cc2)cc1");
    TEST_ASSERT(mol);
    MMFF::MMFFMolProperties *mmffMolProperties =
        new MMFF::MMFFMolProperties(*mol);
    TEST_ASSERT(mmffMolProperties);
    TEST_ASSERT(!mmffMolProperties->isValid());

    delete mol;
    delete mmffMolProperties;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testGithub6728() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "    Testing github6728: crash due to missing stretch-bend params."
      << std::endl;
  RWMol *mol = SmilesToMol("CC(C)(C)CO[H-]F");
  TEST_ASSERT(mol);
  MolOps::addHs(*mol);
  TEST_ASSERT(DGeomHelpers::EmbedMolecule(*mol) >= 0);
  auto field = MMFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  delete mol;
  delete field;
}

//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
int main() {
  RDLog::InitLogs();
  // we get a ton of warnings here about missing Hs... disable them
  boost::logging::disable_logs("rdApp.warning");
#if 1
  testMMFFTyper1();
  testMMFFBuilder1();
  testMMFFBuilder2();
  testIssue239();
  testCalcEnergyPassedCoords();
  testCalcGrad();
  testIssue242();
  testGithub308();
  testSFIssue1653802();
  testSFIssue2378119();
  testMMFFBatch();
  testMMFFParamGetters();
#ifdef RDK_TEST_MULTITHREADED
  testMMFFMultiThread();
  testMMFFMultiThread2();
  testMMFFMultiThread3();
#endif
  testGithub162();
#endif
  testGithub224();
  testGithub6728();
}
