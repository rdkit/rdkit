// $Id$
//
//  Copyright (C) 2004-2010 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <iostream>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>

#include <GraphMol/ForceFieldHelpers/UFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/UFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>
#include <ForceField/ForceField.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <boost/math/special_functions/round.hpp>

using namespace RDKit;
#if 1
void testUFFTyper1() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test UFF atom labels." << std::endl;

  ROMol *mol;
  std::string key;

  mol = SmilesToMol("[SiH3]CC(=O)NC");
  TEST_ASSERT(mol);

  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(0));
  TEST_ASSERT(key == "Si3");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(1));
  TEST_ASSERT(key == "C_3");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(2));
  TEST_ASSERT(key == "C_R");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(3));
  TEST_ASSERT(key == "O_R");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(4));
  TEST_ASSERT(key == "N_R");

  delete mol;
  mol = SmilesToMol("CC(=O)C");
  TEST_ASSERT(mol);

  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(0));
  TEST_ASSERT(key == "C_3");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(1));
  TEST_ASSERT(key == "C_2");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(2));
  TEST_ASSERT(key == "O_2");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(3));
  TEST_ASSERT(key == "C_3");

  delete mol;
  mol = SmilesToMol("C(=O)S");
  TEST_ASSERT(mol);

  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(0));
  TEST_ASSERT(key == "C_2");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(1));
  TEST_ASSERT(key == "O_2");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(2));
  TEST_ASSERT(key == "S_3+2");

  delete mol;
  mol = SmilesToMol("SCS(=O)S(=O)(=O)O");
  TEST_ASSERT(mol);

  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(0));
  TEST_ASSERT(key == "S_3+2");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(1));
  TEST_ASSERT(key == "C_3");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(2));
  TEST_ASSERT(key == "S_3+4");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(4));
  TEST_ASSERT(key == "S_3+6");

  delete mol;
  mol = SmilesToMol("PCP(O)CP(=O)(=O)");
  TEST_ASSERT(mol);

  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(0));
  TEST_ASSERT(key == "P_3+3");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(1));
  TEST_ASSERT(key == "C_3");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(2));
  TEST_ASSERT(key == "P_3+3");
  // FIX: I am not 100% convinced that this should be SP2
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(5));
  TEST_ASSERT(key == "P_3+5");

  delete mol;
  mol = SmilesToMol("C(F)(Cl)(Br)I");
  TEST_ASSERT(mol);

  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(0));
  TEST_ASSERT(key == "C_3");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(1));
  TEST_ASSERT(key == "F_");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(2));
  TEST_ASSERT(key == "Cl");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(3));
  TEST_ASSERT(key == "Br");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(4));
  TEST_ASSERT(key == "I_");

  delete mol;
  mol = SmilesToMol("[Li].[Na].[K].[Rb].[Cs]");
  TEST_ASSERT(mol);

  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(0));
  TEST_ASSERT(key == "Li");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(1));
  TEST_ASSERT(key == "Na");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(2));
  TEST_ASSERT(key == "K_");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(3));
  TEST_ASSERT(key == "Rb");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(4));
  TEST_ASSERT(key == "Cs");

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testUFFTyper2() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test UFF atom typer." << std::endl;

  ROMol *mol, *mol2;
  std::string key;

  mol = SmilesToMol("[SiH3]CC(=O)NC");
  TEST_ASSERT(mol);

  UFF::AtomicParamVect types;
  bool foundAll;
  boost::tie(types, foundAll) = UFF::getAtomTypes(*mol);
  TEST_ASSERT(foundAll);
  TEST_ASSERT(types.size() == mol->getNumAtoms());
  for (UFF::AtomicParamVect::const_iterator it = types.begin();
       it != types.end(); it++) {
    TEST_ASSERT((*it));
  }
  mol2 = MolOps::addHs(*mol);
  delete mol;
  boost::tie(types, foundAll) = UFF::getAtomTypes(*mol2);
  TEST_ASSERT(foundAll);
  TEST_ASSERT(types.size() == mol2->getNumAtoms());
  for (UFF::AtomicParamVect::const_iterator it = types.begin();
       it != types.end(); it++) {
    TEST_ASSERT((*it));
  }

  // connected with sf.net bug 2094445
  mol = SmilesToMol("[SiH2]=C");
  TEST_ASSERT(mol);
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(0));
  TEST_ASSERT(key == "Si3");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(1));
  TEST_ASSERT(key == "C_2");

  mol = SmilesToMol("[AlH]=C");
  TEST_ASSERT(mol);
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(0));
  TEST_ASSERT(key == "Al3");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(1));
  TEST_ASSERT(key == "C_2");

  mol = SmilesToMol("[Mg]=C");
  TEST_ASSERT(mol);
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(0));
  TEST_ASSERT(key == "Mg3+2");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(1));
  TEST_ASSERT(key == "C_2");

  mol = SmilesToMol("[SiH3][Si]([SiH3])=C");
  TEST_ASSERT(mol);
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(0));
  TEST_ASSERT(key == "Si3");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(1));
  TEST_ASSERT(key == "Si3");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(2));
  TEST_ASSERT(key == "Si3");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(3));
  TEST_ASSERT(key == "C_2");

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testUFFBuilder1() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing UFF builder tools." << std::endl;

  ROMol *mol, *mol2;
  std::string key;

  UFF::AtomicParamVect types;
  bool foundAll;
  ForceFields::ForceField *field;
  boost::shared_array<boost::uint8_t> nbrMat;

  mol = SmilesToMol("CC(O)C");
  Conformer *conf = new Conformer(mol->getNumAtoms());
  int cid = static_cast<int>(mol->addConformer(conf, true));
  TEST_ASSERT(mol);
  boost::tie(types, foundAll) = UFF::getAtomTypes(*mol);
  TEST_ASSERT(foundAll);
  TEST_ASSERT(types.size() == mol->getNumAtoms());
  field = new ForceFields::ForceField();
  // add the atomic positions:
  for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
    field->positions().push_back(&((conf->getAtomPos(i))));
  }
  
  UFF::Tools::addBonds(*mol, types, field);

  TEST_ASSERT(field->contribs().size() == 3);

  nbrMat = UFF::Tools::buildNeighborMatrix(*mol);
  // the neighbor matrix is an upper triangular matrix
  // position indices are as follows:
  //  0  1  2  3
  //     4  5  6
  //        7  8
  //           9
  TEST_ASSERT(UFF::Tools::twoBitCellPos(mol->getNumAtoms(), 1, 1) == 4);
  TEST_ASSERT(UFF::Tools::twoBitCellPos(mol->getNumAtoms(), 2, 1) == 5);
  TEST_ASSERT(UFF::Tools::twoBitCellPos(mol->getNumAtoms(), 1, 2) == 5);
  TEST_ASSERT(UFF::Tools::getTwoBitCell(nbrMat, 0) == UFF::Tools::RELATION_1_X);
  TEST_ASSERT(UFF::Tools::getTwoBitCell(nbrMat, 1) == UFF::Tools::RELATION_1_2);
  TEST_ASSERT(UFF::Tools::getTwoBitCell(nbrMat, 2) == UFF::Tools::RELATION_1_3);
  TEST_ASSERT(UFF::Tools::getTwoBitCell(nbrMat, 3) == UFF::Tools::RELATION_1_3);

  UFF::Tools::addAngles(*mol, types, field);
  TEST_ASSERT(field->contribs().size() == 6);

  // there are no non-bonded terms here:
  UFF::Tools::addNonbonded(*mol, cid, types, field, nbrMat);
  TEST_ASSERT(field->contribs().size() == 6);

  // and no torsions either, until we add Hs:
  UFF::Tools::addTorsions(*mol, types, field);
  TEST_ASSERT(field->contribs().size() == 6);

  delete mol;
  delete field;
  mol = SmilesToMol("CCOC");
  Conformer *conf2 = new Conformer(mol->getNumAtoms());
  cid = static_cast<int>(mol->addConformer(conf2, true));
  TEST_ASSERT(mol);
  boost::tie(types, foundAll) = UFF::getAtomTypes(*mol);
  TEST_ASSERT(foundAll);
  TEST_ASSERT(types.size() == mol->getNumAtoms());
  field = new ForceFields::ForceField();
  // add the atomic positions:
  for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
    field->positions().push_back(&(conf2->getAtomPos(i)));
  }
  
  UFF::Tools::addBonds(*mol, types, field);

  TEST_ASSERT(field->contribs().size() == 3);

  nbrMat = UFF::Tools::buildNeighborMatrix(*mol);
  TEST_ASSERT(UFF::Tools::getTwoBitCell(nbrMat, 0) == UFF::Tools::RELATION_1_X);
  TEST_ASSERT(UFF::Tools::getTwoBitCell(nbrMat, 1) == UFF::Tools::RELATION_1_2);
  TEST_ASSERT(UFF::Tools::getTwoBitCell(nbrMat, 2) == UFF::Tools::RELATION_1_3);
  TEST_ASSERT(UFF::Tools::getTwoBitCell(nbrMat, 3) == UFF::Tools::RELATION_1_X);

  UFF::Tools::addAngles(*mol, types, field);
  TEST_ASSERT(field->contribs().size() == 5);
  UFF::Tools::addNonbonded(*mol, cid, types, field, nbrMat);
  TEST_ASSERT(field->contribs().size() == 6);
  UFF::Tools::addTorsions(*mol, types, field);
  TEST_ASSERT(field->contribs().size() == 7);

  delete mol;
  delete field;
  mol = SmilesToMol("CO");
  Conformer *conf3 = new Conformer(mol->getNumAtoms());
  cid = static_cast<int>(mol->addConformer(conf3, true));
  TEST_ASSERT(mol);
  boost::tie(types, foundAll) = UFF::getAtomTypes(*mol);
  TEST_ASSERT(foundAll);
  TEST_ASSERT(types.size() == mol->getNumAtoms());

  field = new ForceFields::ForceField();
  // add the atomic positions:
  for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
    field->positions().push_back(&(conf3->getAtomPos(i)));
  }
  
  UFF::Tools::addBonds(*mol, types, field);

  TEST_ASSERT(field->contribs().size() == 1);

  nbrMat = UFF::Tools::buildNeighborMatrix(*mol);
  TEST_ASSERT(UFF::Tools::getTwoBitCell(nbrMat, 0) == UFF::Tools::RELATION_1_X);
  TEST_ASSERT(UFF::Tools::getTwoBitCell(nbrMat, 1) == UFF::Tools::RELATION_1_2);

  UFF::Tools::addAngles(*mol, types, field);
  TEST_ASSERT(field->contribs().size() == 1);
  UFF::Tools::addNonbonded(*mol, cid, types, field, nbrMat);
  TEST_ASSERT(field->contribs().size() == 1);
  UFF::Tools::addTorsions(*mol, types, field);
  TEST_ASSERT(field->contribs().size() == 1);

  mol2 = MolOps::addHs(*mol);
  TEST_ASSERT(mol2->getNumAtoms() == 6);
  delete field;

  boost::tie(types, foundAll) = UFF::getAtomTypes(*mol2);
  TEST_ASSERT(foundAll);
  TEST_ASSERT(types.size() == mol2->getNumAtoms());
  Conformer *conf4 = new Conformer(mol2->getNumAtoms());
  cid = static_cast<int>(mol2->addConformer(conf4, true));

  field = new ForceFields::ForceField();
  // add the atomic positions:
  for (unsigned int i = 0; i < mol2->getNumAtoms(); ++i) {
    field->positions().push_back(&(conf4->getAtomPos(i)));
  }
  
  UFF::Tools::addBonds(*mol2, types, field);
  TEST_ASSERT(field->contribs().size() == 5);

  nbrMat = UFF::Tools::buildNeighborMatrix(*mol2);
  UFF::Tools::addAngles(*mol2, types, field);
  TEST_ASSERT(field->contribs().size() == 12);
  UFF::Tools::addNonbonded(*mol2, cid, types, field, nbrMat);
  TEST_ASSERT(field->contribs().size() == 15);
  UFF::Tools::addTorsions(*mol2, types, field);
  TEST_ASSERT(field->contribs().size() == 18);
  delete mol2;

  delete mol;
  delete field;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testUFFBuilder2() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing UFF builder+minimization." << std::endl;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/UFF/test_data";
  {
    RWMol *mol = MolFileToMol(pathName + "/small1.mol", false);
    TEST_ASSERT(mol);
    MolOps::sanitizeMol(*mol);

    ForceFields::ForceField *field;
    field = UFF::constructForceField(*mol);
    TEST_ASSERT(field);
    field->initialize();
    int needMore = field->minimize();
    TEST_ASSERT(!needMore);
    // std::cout << MolToMolBlock(mol);
    delete mol;
    delete field;
  }

  {  // make sure the confId argument works
    RWMol *mol = MolFileToMol(pathName + "/small1.mol", false);
    TEST_ASSERT(mol);
    MolOps::sanitizeMol(*mol);
    Conformer *newConf = new Conformer(mol->getConformer());
    newConf->setId(111);
    mol->addConformer(newConf, false);
    RDGeom::Point3D p0 = mol->getConformer().getAtomPos(0);
    RDGeom::Point3D p1 = mol->getConformer().getAtomPos(1);

    ForceFields::ForceField *field;
    field = UFF::constructForceField(*mol, 100, 111);
    TEST_ASSERT(field);
    field->initialize();
    int needMore = field->minimize();
    TEST_ASSERT(!needMore);
    // std::cout << MolToMolBlock(mol);

    RDGeom::Point3D np0 = mol->getConformer().getAtomPos(0);
    RDGeom::Point3D np1 = mol->getConformer().getAtomPos(1);
    TEST_ASSERT(feq(p0.x, np0.x));
    TEST_ASSERT(feq(p0.y, np0.y));
    TEST_ASSERT(feq(p0.z, np0.z));
    TEST_ASSERT(feq(p1.x, np1.x));
    TEST_ASSERT(feq(p1.y, np1.y));
    TEST_ASSERT(feq(p1.z, np1.z));

    delete mol;
    delete field;
  }

  {
    RWMol *mol = MolFileToMol(pathName + "/small2.mol", false);
    TEST_ASSERT(mol);
    MolOps::sanitizeMol(*mol);

    ForceFields::ForceField *field;
    field = UFF::constructForceField(*mol);
    TEST_ASSERT(field);
    field->initialize();
    int needMore = field->minimize(150);
    TEST_ASSERT(!needMore);
    // std::cout << MolToMolBlock(mol);
    delete mol;
    delete field;
  }
  {
    RWMol *mol = MolFileToMol(pathName + "/benzene.mol", false);
    TEST_ASSERT(mol);
    MolOps::sanitizeMol(*mol);

    ForceFields::ForceField *field;
    field = UFF::constructForceField(*mol);
    TEST_ASSERT(field);
    field->initialize();
    int needMore = field->minimize();
    TEST_ASSERT(!needMore);
    // std::cout << MolToMolBlock(mol);
    delete mol;
    delete field;
  }
  {
    RWMol *mol = MolFileToMol(pathName + "/toluene.mol", false);
    TEST_ASSERT(mol);
    MolOps::sanitizeMol(*mol);

    ForceFields::ForceField *field;
    field = UFF::constructForceField(*mol);
    TEST_ASSERT(field);
    field->initialize();
    int needMore = field->minimize();
    TEST_ASSERT(!needMore);
    // std::cout << MolToMolBlock(mol);
    delete mol;
    delete field;
  }
  {
    RWMol *mol = MolFileToMol(pathName + "/complex1.mol", false);
    TEST_ASSERT(mol);
    MolOps::sanitizeMol(*mol);

    ForceFields::ForceField *field;
    field = UFF::constructForceField(*mol);
    TEST_ASSERT(field);
    field->initialize();
    int needMore = field->minimize();
    TEST_ASSERT(!needMore);
    // std::cout << MolToMolBlock(mol);
    delete mol;
    delete field;
  }

  {  // test the convenience function
    RWMol *mol = MolFileToMol(pathName + "/small1.mol", false);
    TEST_ASSERT(mol);
    MolOps::sanitizeMol(*mol);

    UFF::AtomicParamVect types;
    bool foundAll;
    boost::shared_array<boost::uint8_t> nbrMat;
    boost::tie(types, foundAll) = UFF::getAtomTypes(*mol);
    
    ForceFields::ForceField *field;
    field = new ForceFields::ForceField();

    // add the atomic positions:
    for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
      field->positions().push_back(&((mol->getConformer().getAtomPos(i))));
    }
    
    UFF::Tools::addBonds(*mol, types, field);
    nbrMat = UFF::Tools::buildNeighborMatrix(*mol);
    UFF::Tools::addAngles(*mol, types, field);
    UFF::Tools::addTorsions(*mol, types, field);
    // std::cout << field->contribs().size() << std::endl;
    UFF::Tools::addNonbonded(*mol, 0, types, field, nbrMat);
    delete field;

    field = UFF::constructForceField(*mol);
    field->initialize();
    field->minimize();
    delete field;
        
    RWMol *mol2 = new RWMol(*mol);

    field = UFF::constructForceField(*mol2);
    TEST_ASSERT(field);
    field->initialize();
    int needMore = field->minimize();
    TEST_ASSERT(!needMore);
    delete field;

    needMore = UFF::UFFOptimizeMolecule(*mol2).first;
    TEST_ASSERT(!needMore);
    for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
      const RDGeom::Point3D p1 = mol->getConformer().getAtomPos(i);
      const RDGeom::Point3D p2 = mol2->getConformer().getAtomPos(i);
      TEST_ASSERT(feq(p1.x, p2.x));
      TEST_ASSERT(feq(p1.y, p2.y));
      TEST_ASSERT(feq(p1.z, p2.z));
    }

    delete mol;
    delete mol2;
  }

  {  // test the convenience function for all confs
    RWMol *mol = MolFileToMol(pathName + "/complex1.mol", false);
    TEST_ASSERT(mol);
    MolOps::sanitizeMol(*mol);
    RWMol *mol2 = new RWMol(*mol);
    Conformer *newConf = new Conformer(mol->getConformer());
    newConf->setId(111);
    mol->addConformer(newConf, false);

    ForceFields::ForceField *field;
    field = UFF::constructForceField(*mol, 100, 111);
    TEST_ASSERT(field);
    field->initialize();
    int needMore = field->minimize();
    TEST_ASSERT(!needMore);
    delete field;

    // the first conf is the same as above,
    // but we add a second that's already minimized
    newConf = new Conformer(mol->getConformer(111));
    newConf->setId(112);
    mol2->addConformer(newConf, false);

    std::vector<std::pair<int, double> > res;
    UFF::UFFOptimizeMoleculeConfs(*mol2, res);
    TEST_ASSERT(res.size() == 2);
    TEST_ASSERT(!res[0].first);
    TEST_ASSERT(!res[1].first);
    // we expect the energy to go down at least a little bit.
    TEST_ASSERT(res[1].second < res[0].second);

    for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
      const RDGeom::Point3D p1 = mol->getConformer(111).getAtomPos(i);
      const RDGeom::Point3D p2 = mol2->getConformer(0).getAtomPos(i);
      TEST_ASSERT(feq(p1.x, p2.x));
      TEST_ASSERT(feq(p1.y, p2.y));
      TEST_ASSERT(feq(p1.z, p2.z));
    }

    delete mol;
    delete mol2;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testUFFBatch() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "    Testing bulk UFF (regression to check that things run)."
      << std::endl;

  ROMol *mol;
  std::string key;

  ForceFields::ForceField *field;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/UFF/test_data";
  SDMolSupplier suppl(pathName + "/bulk.sdf", false);

  int count = 0;
  mol = suppl.next();
  while (mol && !suppl.atEnd()) {
    count++;
    MolOps::sanitizeMol(*(RWMol *)mol);
    std::string origMolBlock = MolToMolBlock(*mol);

    BOOST_LOG(rdErrorLog) << "Mol:" << count << std::endl;
    try {
      field = UFF::constructForceField(*mol);
    } catch (...) {
      field = 0;
    }
    if (field) {
      field->initialize();
      int failed = field->minimize(500);
      if (failed) {
        BOOST_LOG(rdErrorLog) << " not converged (code=" << failed << ")"
                              << std::endl;
        std::cout << origMolBlock << "$$$$" << std::endl;
        std::cout << MolToMolBlock(*mol) << "$$$$" << std::endl;
      }
      delete field;
    }
    delete mol;
    mol = suppl.next();
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testUFFBuilderSpecialCases() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing UFF special cases." << std::endl;

  RWMol *mol;
  std::string key;
  int needMore;
  RDGeom::Point3D v1, v2;
  ForceFields::ForceField *field;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/UFF/test_data";
  // ----------
  //  Trigonal bipyramid
  // ----------
  mol = MolFileToMol(pathName + "/tbp.mol", false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);

  const Conformer &conf = mol->getConformer();
  field = UFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  needMore = field->minimize(200, 1e-8, 1e-4);
  TEST_ASSERT(!needMore);
  v1 = conf.getAtomPos(0).directionVector(conf.getAtomPos(1));
  v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(2));
  TEST_ASSERT(feq(v1.dotProduct(v2), -1.0, 1e-3));
  v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(3));
  TEST_ASSERT(feq(v1.dotProduct(v2), 0.0, 1e-3));
  v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(4));
  TEST_ASSERT(feq(v1.dotProduct(v2), 0.0, 1e-3));
  v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(5));
  TEST_ASSERT(feq(v1.dotProduct(v2), 0.0, 1e-3));

  v1 = conf.getAtomPos(0).directionVector(conf.getAtomPos(2));
  v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(3));
  TEST_ASSERT(feq(v1.dotProduct(v2), 0.0, 1e-3));
  v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(4));
  TEST_ASSERT(feq(v1.dotProduct(v2), 0.0, 1e-3));
  v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(5));
  TEST_ASSERT(feq(v1.dotProduct(v2), 0.0, 1e-3));

  v1 = conf.getAtomPos(0).directionVector(conf.getAtomPos(3));
  v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(4));
  TEST_ASSERT(feq(v1.dotProduct(v2), -0.5, 1e-3));
  v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(5));
  TEST_ASSERT(feq(v1.dotProduct(v2), -0.5, 1e-3));

  v1 = conf.getAtomPos(0).directionVector(conf.getAtomPos(4));
  v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(5));
  TEST_ASSERT(feq(v1.dotProduct(v2), -0.5, 1e-3));

  delete mol;
  delete field;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testIssue239() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing Issue239." << std::endl;

  RWMol *mol;
  int needMore;
  (void)needMore;  // add test later
  ForceFields::ForceField *field;
  double e1, e2;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/UFF/test_data";
  mol = MolFileToMol(pathName + "/Issue239.mol", false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);

  field = UFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  needMore = field->minimize(200, 1e-6, 1e-3);
  e1 = field->calcEnergy();
  needMore = field->minimize(200, 1e-6, 1e-3);
  e2 = field->calcEnergy();
  TEST_ASSERT(feq(e2, e1, 0.1));

  delete mol;
  delete field;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}
#endif

void testIssue242() {
#if 0
// FIX: Changes to the forcefield (connected to Issue 408) have
// made it so that this particular problem no longer manifests
// in this molecule/embedding. A new test case is needed.
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing Issue242." << std::endl;

  RWMol *mol,*mol2;
  int needMore;
  ForceFields::ForceField *field=0,*field2=0;
  std::string mb1,mb2;
  double e1,e2;

  std::string pathName=getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/UFF/test_data";

  mol = MolFileToMol(pathName+"/Issue242.mol");
  TEST_ASSERT(mol);

  mol2 = MolFileToMol(pathName+"/Issue242.mol");
  TEST_ASSERT(mol2);

  TEST_ASSERT(DGeomHelpers::EmbedMolecule(*mol2,30,2300)>=0);
  mb1 = MolToMolBlock(*mol);
  mb2 = MolToMolBlock(*mol2);

  //BOOST_LOG(rdInfoLog) << "\nMol1\n" << mb1 << std::endl;
  //BOOST_LOG(rdInfoLog) << "\nMol2\n" << mb2 << std::endl;

  
  field=UFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  field2=UFF::constructForceField(*mol2);
  TEST_ASSERT(field2);
  field2->initialize();
  e1 = field->calcEnergy();
  e2 = field2->calcEnergy();
  BOOST_LOG(rdInfoLog) << "E1: " << e1 << std::endl;
  BOOST_LOG(rdInfoLog) << "E2: " << e2 << std::endl;
  //TEST_ASSERT(feq(e2,e1,0.1));

  needMore = field->minimize(600,1e-4);
  TEST_ASSERT(!needMore)
  needMore = field2->minimize(600,1e-4);
  TEST_ASSERT(!needMore)
  e1 = field->calcEnergy();
  e2 = field2->calcEnergy();
  BOOST_LOG(rdInfoLog) << "E1: " << e1 << std::endl;
  BOOST_LOG(rdInfoLog) << "E2: " << e2 << std::endl;
  TEST_ASSERT(feq(e2,e1,1.0));
  
  needMore = field->minimize(600,1e-4);
  TEST_ASSERT(!needMore)
  needMore = field2->minimize(600,1e-4);
  TEST_ASSERT(!needMore)
  e1 = field->calcEnergy();
  e2 = field2->calcEnergy();
  BOOST_LOG(rdInfoLog) << "rE1: " << e1 << std::endl;
  BOOST_LOG(rdInfoLog) << "rE2: " << e2 << std::endl;
  TEST_ASSERT(feq(e2,e1,1.0));
  
  delete mol;
  delete mol2;
  delete field;
  delete field2;

  mol = MolFileToMol(pathName+"/Issue242-2.mol");
  TEST_ASSERT(mol);

  mol2 = MolFileToMol(pathName+"/Issue242-2.mol");
  TEST_ASSERT(mol2);

  TEST_ASSERT(DGeomHelpers::EmbedMolecule(*mol2,30,2370)>=0);
  mb1 = MolToMolBlock(*mol);
  mb2 = MolToMolBlock(*mol2);

  //std::cout << mb1 << "------\n";
  //std::cout << mb2 << "------\n";

  field=UFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  field2=UFF::constructForceField(*mol2);
  TEST_ASSERT(field2);
  field2->initialize();
  e1 = field->calcEnergy();
  e2 = field2->calcEnergy();
  BOOST_LOG(rdInfoLog) << "E1: " << e1 << std::endl;
  BOOST_LOG(rdInfoLog) << "E2: " << e2 << std::endl;
  //TEST_ASSERT(feq(e2,e1,0.1));

  needMore = field->minimize(200,1e-4);
  needMore = field2->minimize(200,1e-4);
  e1 = field->calcEnergy();
  e2 = field2->calcEnergy();
  BOOST_LOG(rdInfoLog) << "E1: " << e1 << std::endl;
  BOOST_LOG(rdInfoLog) << "E2: " << e2 << std::endl;
  TEST_ASSERT(feq(e2,e1,0.1));

  delete mol;
  delete mol2;
  delete field;
  delete field2;
  

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
#endif
}

void testSFIssue1653802() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing SFIssue1653802." << std::endl;

  RWMol *mol;
  int needMore;
  ForceFields::ForceField *field;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/UFF/test_data";

  mol = MolFileToMol(pathName + "/cyclobutadiene.mol", false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);

  UFF::AtomicParamVect types;
  bool foundAll;
  boost::shared_array<boost::uint8_t> nbrMat;
  boost::tie(types, foundAll) = UFF::getAtomTypes(*mol);
  TEST_ASSERT(foundAll);
  TEST_ASSERT(types.size() == mol->getNumAtoms());
  field = new ForceFields::ForceField();
  // add the atomic positions:
  for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
    field->positions().push_back(&((mol->getConformer().getAtomPos(i))));
  }
  
  UFF::Tools::addBonds(*mol, types, field);
  TEST_ASSERT(field->contribs().size() == 8);

  nbrMat = UFF::Tools::buildNeighborMatrix(*mol);
  UFF::Tools::addAngles(*mol, types, field);
  TEST_ASSERT(field->contribs().size() == 20);
  UFF::Tools::addTorsions(*mol, types, field);
  // std::cout << field->contribs().size() << std::endl;
  TEST_ASSERT(field->contribs().size() == 36);
  UFF::Tools::addNonbonded(*mol, 0, types, field, nbrMat);
  delete field;

  field = UFF::constructForceField(*mol);
  field->initialize();
  needMore = field->minimize(200, 1e-6, 1e-3);
  TEST_ASSERT(!needMore);

  delete mol;
  delete field;
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testSFIssue2378119() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing SFIssue2378119." << std::endl;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/UFF/test_data";
  {
    RWMol *mol = MolFileToMol(pathName + "/Issue2378119.mol");
    TEST_ASSERT(mol);
    ForceFields::ForceField *field = UFF::constructForceField(*mol);
    TEST_ASSERT(field);
    field->initialize();
    double e1 = field->calcEnergy();
    TEST_ASSERT(e1 > 0.0 && e1 < 1e12);

    int needMore = field->minimize(200, 1e-6, 1e-3);
    TEST_ASSERT(!needMore);
    double e2 = field->calcEnergy();
    TEST_ASSERT(e2 < e1);

    delete mol;
    delete field;
  }
  {
    RWMol *mol = MolFileToMol(pathName + "/Issue2378119.2.mol");
    TEST_ASSERT(mol);
    ForceFields::ForceField *field = UFF::constructForceField(*mol);
    TEST_ASSERT(field);
    field->initialize();
    double e1 = field->calcEnergy();
    TEST_ASSERT(e1 == 0.0);

    int needMore = field->minimize(200, 1e-6, 1e-3);
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
        UFF::constructForceField(*mol, 100.0, -1, false);
    TEST_ASSERT(field);
    field->initialize();
    double e1 = field->calcEnergy();
    TEST_ASSERT(e1 > 0.0 && e1 < 1e12);

    int needMore = field->minimize(200, 1e-6, 1e-3);
    TEST_ASSERT(!needMore);
    double e2 = field->calcEnergy();
    TEST_ASSERT(e2 < e1);

    delete mol;
    delete field;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMissingParams() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing handling missing parameters."
                        << std::endl;

  {
    UFF::AtomicParamVect types;
    bool foundAll;

    ROMol *mol = SmilesToMol("[Cu](C)(C)(C)(C)C");
    TEST_ASSERT(mol);

    ROMol *mol2 = MolOps::addHs(*mol);
    delete mol;

    TEST_ASSERT(DGeomHelpers::EmbedMolecule(*mol2) >= 0);

    boost::tie(types, foundAll) = UFF::getAtomTypes(*mol2);
    TEST_ASSERT(!foundAll);
    TEST_ASSERT(types.size() == mol2->getNumAtoms());
    TEST_ASSERT(!types[0]);

    // make sure we can optimize anyway:
    ForceFields::ForceField *field = UFF::constructForceField(*mol2, types);
    TEST_ASSERT(field);
    field->initialize();
    double e1 = field->calcEnergy();
    int needMore = field->minimize();
    TEST_ASSERT(needMore);
    double e2 = field->calcEnergy();
    TEST_ASSERT(e2 < e1);
    // std::cerr<<" DE: "<<e1<<" -> "<<e2<<std::endl;
    delete mol2;
    delete field;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}
void testSFIssue3009337() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing SFIssue3009337." << std::endl;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/UFF/test_data";
  {
    RWMol *mol = MolFileToMol(pathName + "/Issue3009337.mol", true, false);
    TEST_ASSERT(mol);
    ForceFields::ForceField *field = UFF::constructForceField(*mol);
    TEST_ASSERT(field);
    field->initialize();
    double e1 = field->calcEnergy();
    TEST_ASSERT(e1 > 0.0 && e1 < 1e12);

    int needMore = field->minimize(200, 1e-6, 1e-3);
    TEST_ASSERT(!needMore);
    double e2 = field->calcEnergy();
    TEST_ASSERT(e2 < e1);

    delete mol;
    delete field;
  }
  {
    RWMol *mol = MolFileToMol(pathName + "/Issue3009337.2.mol", true, false);
    TEST_ASSERT(mol);
    ForceFields::ForceField *field = UFF::constructForceField(*mol);
    TEST_ASSERT(field);
    field->initialize();
    double e1 = field->calcEnergy();
    TEST_ASSERT(e1 > 0.0 && e1 < 1e12);

    int needMore = field->minimize(200, 1e-6, 1e-3);
    TEST_ASSERT(!needMore);
    double e2 = field->calcEnergy();
    TEST_ASSERT(e2 < e1);

    delete mol;
    delete field;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}
void testGitHubIssue62() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing GitHubIssue62." << std::endl;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/UFF/test_data";
  {
    double energyValues[] = {
        38.687, 174.698, 337.986, 115.248, 2.482,   1.918,  10.165,  98.469,
        39.078, 267.236, 15.747,  202.121, 205.539, 20.044, 218.986, 79.627};
    SmilesMolSupplier smiSupplier(pathName + "/Issue62.smi");
    SDWriter *sdfWriter = new SDWriter(pathName + "/Issue62.sdf");
    for (unsigned int i = 0; i < smiSupplier.length(); ++i) {
      ROMol *mol = MolOps::addHs(*(smiSupplier[i]));
      TEST_ASSERT(mol);
      std::string molName = "";
      if (mol->hasProp(common_properties::_Name)) {
        mol->getProp(common_properties::_Name, molName);
      }
      DGeomHelpers::EmbedMolecule(*mol);
      ForceFields::ForceField *field = UFF::constructForceField(*mol);
      TEST_ASSERT(field);
      field->initialize();
      int needMore = field->minimize(200, 1.e-6, 1.e-3);
      TEST_ASSERT(!needMore);
      sdfWriter->write(*mol);
      double e = field->calcEnergy();
      BOOST_LOG(rdErrorLog) << molName << " " << e << std::endl;
      TEST_ASSERT(fabs(e - energyValues[i]) < 1.);
    }
    sdfWriter->close();
    BOOST_LOG(rdErrorLog) << "  done" << std::endl;
  }
}
void testUFFParamGetters() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test UFF force-field parameter getters."
                        << std::endl;
  {
    ROMol *mol = SmilesToMol("c1ccccc1CCNN");
    TEST_ASSERT(mol);
    ROMol *molH = MolOps::addHs(*mol);
    TEST_ASSERT(molH);
    ForceFields::UFF::UFFBond uffBondStretchParams;
    TEST_ASSERT(
        UFF::getUFFBondStretchParams(*molH, 6, 7, uffBondStretchParams));
    TEST_ASSERT(
        ((int)boost::math::round(uffBondStretchParams.kb * 1000) == 699592) &&
        ((int)boost::math::round(uffBondStretchParams.r0 * 1000) == 1514));
    TEST_ASSERT(
        !UFF::getUFFBondStretchParams(*molH, 0, 7, uffBondStretchParams));
    ForceFields::UFF::UFFAngle uffAngleBendParams;
    TEST_ASSERT(UFF::getUFFAngleBendParams(*molH, 6, 7, 8, uffAngleBendParams));
    TEST_ASSERT(
        ((int)boost::math::round(uffAngleBendParams.ka * 1000) == 303297) &&
        ((int)boost::math::round(uffAngleBendParams.theta0 * 1000) == 109470));
    TEST_ASSERT(
        !UFF::getUFFAngleBendParams(*molH, 0, 7, 8, uffAngleBendParams));
    ForceFields::UFF::UFFTor uffTorsionParams;
    TEST_ASSERT(UFF::getUFFTorsionParams(*molH, 6, 7, 8, 9, uffTorsionParams));
    TEST_ASSERT(((int)boost::math::round(uffTorsionParams.V * 1000) == 976));
    TEST_ASSERT(!UFF::getUFFTorsionParams(*molH, 0, 7, 8, 9, uffTorsionParams));
    ForceFields::UFF::UFFInv uffInversionParams;
    TEST_ASSERT(
        UFF::getUFFInversionParams(*molH, 6, 5, 4, 0, uffInversionParams));
    TEST_ASSERT(((int)boost::math::round(uffInversionParams.K * 1000) == 2000));
    TEST_ASSERT(
        !UFF::getUFFInversionParams(*molH, 6, 5, 4, 1, uffInversionParams));
    ForceFields::UFF::UFFVdW uffVdWParams;
    TEST_ASSERT(UFF::getUFFVdWParams(*molH, 0, 9, uffVdWParams));
    TEST_ASSERT(((int)boost::math::round(uffVdWParams.x_ij * 1000) == 3754) &&
                ((int)boost::math::round(uffVdWParams.D_ij * 1000) == 85));
  }
}

#ifdef RDK_TEST_MULTITHREADED
namespace {
void runblock_uff(const std::vector<ROMol *> &mols,
                  const std::vector<double> &energies, unsigned int count,
                  unsigned int idx) {
  for (unsigned int rep = 0; rep < 1000; ++rep) {
    for (unsigned int i = 0; i < mols.size(); ++i) {
      if (i % count != idx) continue;
      ROMol *mol = mols[i];
      ForceFields::ForceField *field = 0;
      if (!(rep % 100)) {
        BOOST_LOG(rdErrorLog) << "Rep: " << rep << " Mol:" << i << std::endl;
      }
      try {
        field = UFF::constructForceField(*mol);
      } catch (...) {
        field = 0;
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
}
#include <boost/thread.hpp>
void testUFFMultiThread() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test UFF multithreading" << std::endl;

  // ForceFields::ForceField *field;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/UFF/test_data";
  SDMolSupplier suppl(pathName + "/bulk.sdf");
  std::vector<ROMol *> mols;
  while (!suppl.atEnd() && mols.size() < 100) {
    ROMol *mol = 0;
    try {
      mol = suppl.next();
    } catch (...) {
      continue;
    }
    if (!mol) continue;
    mols.push_back(mol);
  }

  std::cerr << "generating reference data" << std::endl;
  std::vector<double> energies(mols.size(), 0.0);
  for (unsigned int i = 0; i < mols.size(); ++i) {
    ROMol mol(*mols[i]);
    ForceFields::ForceField *field = 0;
    try {
      field = UFF::constructForceField(mol);
    } catch (...) {
      field = 0;
    }
    TEST_ASSERT(field);
    field->initialize();
    int failed = field->minimize(500);
    TEST_ASSERT(!failed);
    energies[i] = field->calcEnergy();
    delete field;
  }

  boost::thread_group tg;

  std::cerr << "processing" << std::endl;
  unsigned int count = 4;
  for (unsigned int i = 0; i < count; ++i) {
    std::cerr << " launch :" << i << std::endl;
    std::cerr.flush();
    tg.add_thread(new boost::thread(runblock_uff, mols, energies, count, i));
  }
  tg.join_all();

  BOOST_FOREACH (ROMol *mol, mols) { delete mol; }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testUFFMultiThread2() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test UFF multithreading2" << std::endl;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/UFF/test_data";
  SDMolSupplier suppl(pathName + "/bulk.sdf");
  ROMol *m = suppl[4];
  TEST_ASSERT(m);
  ROMol *om = new ROMol(*m);
  for (unsigned int i = 0; i < 1000; ++i) {
    m->addConformer(new Conformer(m->getConformer()), true);
  }
  std::vector<std::pair<int, double> > res;

  UFF::UFFOptimizeMolecule(*om);
  UFF::UFFOptimizeMoleculeConfs(*m, res, 0);
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
#endif

void testGitHubIssue613() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Github Issue 613: UFF Atom type not "
                           "properly assigned to lanthanides." << std::endl;
  {
    ROMol *mol =
        SmilesToMol("[Eu+3]123456.[Cl]1.[Cl]2.[Cl]3.[Cl]4.[Cl]5.[Cl]6");
    TEST_ASSERT(mol);
    mol->getAtomWithIdx(0)->setHybridization(Atom::SP3D2);

    UFF::AtomicParamVect types;
    bool foundAll;
    boost::tie(types, foundAll) = UFF::getAtomTypes(*mol);
    TEST_ASSERT(foundAll);
    TEST_ASSERT(types.size() == mol->getNumAtoms());

    ForceFields::UFF::ParamCollection *params =
        ForceFields::UFF::ParamCollection::getParams();
    const ForceFields::UFF::AtomicParams *ap = (*params)("Eu6+3");
    TEST_ASSERT(ap);
    TEST_ASSERT(ap->r1 == types[0]->r1);
    TEST_ASSERT(ap->theta0 == types[0]->theta0);
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
int main() {
  RDLog::InitLogs();
#if 1
  testUFFTyper1();
  testUFFTyper2();
  testUFFBuilder1();
  testUFFBuilder2();
  testUFFBatch();
  testUFFBuilderSpecialCases();
  testIssue239();
  testIssue242();
  testSFIssue1653802();
  testSFIssue2378119();
  testUFFParamGetters();
  testMissingParams();
  testSFIssue3009337();
#ifdef RDK_TEST_MULTITHREADED
  testUFFMultiThread();
  testUFFMultiThread2();
#endif
  testGitHubIssue62();
#endif
  testGitHubIssue613();
}
