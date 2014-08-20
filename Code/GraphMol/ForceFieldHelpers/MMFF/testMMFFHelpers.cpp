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

#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <ForceField/ForceField.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>

using namespace RDKit;
void testMMFFTyper1()
{
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test MMFF atom types." << std::endl;

  {
    boost::uint8_t type;
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
    boost::uint8_t type;
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
    boost::uint8_t type;
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
    boost::uint8_t type;
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
    boost::uint8_t type;
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
    boost::uint8_t type;
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

void testMMFFBuilder1()
{
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing MMFF builder tools." << std::endl;

  ROMol *mol,*mol2;

  ForceFields::ForceField *field;
  boost::shared_array<boost::uint8_t> nbrMat;

  mol = SmilesToMol("CC(O)C");
  Conformer *conf = new Conformer(mol->getNumAtoms());
  int cid = static_cast<int>(mol->addConformer(conf, true));
  TEST_ASSERT(mol);
  MMFF::MMFFMolProperties *mmffMolProperties = new MMFF::MMFFMolProperties(*mol);
  TEST_ASSERT(mmffMolProperties);
  TEST_ASSERT(mmffMolProperties->isValid());
  field = new ForceFields::ForceField();
  MMFF::Tools::addBonds(*mol, mmffMolProperties, field);

  TEST_ASSERT(field->contribs().size() == 3);

  nbrMat = MMFF::Tools::buildNeighborMatrix(*mol);
  TEST_ASSERT(MMFF::Tools::getTwoBitCell(nbrMat, 0) == MMFF::Tools::RELATION_1_X);
  TEST_ASSERT(MMFF::Tools::getTwoBitCell(nbrMat, 1) == MMFF::Tools::RELATION_1_2);
  TEST_ASSERT(MMFF::Tools::getTwoBitCell(nbrMat, 2) == MMFF::Tools::RELATION_1_3);
  TEST_ASSERT(MMFF::Tools::getTwoBitCell(nbrMat, 3) == MMFF::Tools::RELATION_1_3);

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
  Conformer *conf2 = new Conformer(mol->getNumAtoms());
  cid = static_cast<int>(mol->addConformer(conf2, true));
  TEST_ASSERT(mol);
  mmffMolProperties = new MMFF::MMFFMolProperties(*mol);
  TEST_ASSERT(mmffMolProperties);
  TEST_ASSERT(mmffMolProperties->isValid());
  field = new ForceFields::ForceField();
  MMFF::Tools::addBonds(*mol, mmffMolProperties, field);

  TEST_ASSERT(field->contribs().size() == 3);

  nbrMat = MMFF::Tools::buildNeighborMatrix(*mol);
  TEST_ASSERT(MMFF::Tools::getTwoBitCell(nbrMat, 0) == MMFF::Tools::RELATION_1_X);
  TEST_ASSERT(MMFF::Tools::getTwoBitCell(nbrMat, 1) == MMFF::Tools::RELATION_1_2);
  TEST_ASSERT(MMFF::Tools::getTwoBitCell(nbrMat, 2) == MMFF::Tools::RELATION_1_3);
  TEST_ASSERT(MMFF::Tools::getTwoBitCell(nbrMat, 3) == MMFF::Tools::RELATION_1_4);

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
  Conformer *conf3 = new Conformer(mol->getNumAtoms());
  cid = static_cast<int>(mol->addConformer(conf3, true));
  TEST_ASSERT(mol);
  mmffMolProperties = new MMFF::MMFFMolProperties(*mol);
  TEST_ASSERT(mmffMolProperties);
  TEST_ASSERT(mmffMolProperties->isValid());

  field = new ForceFields::ForceField();
  MMFF::Tools::addBonds(*mol, mmffMolProperties, field);

  TEST_ASSERT(field->contribs().size() == 1);

  nbrMat = MMFF::Tools::buildNeighborMatrix(*mol);
  TEST_ASSERT(MMFF::Tools::getTwoBitCell(nbrMat, 0) == MMFF::Tools::RELATION_1_X);
  TEST_ASSERT(MMFF::Tools::getTwoBitCell(nbrMat, 1) == MMFF::Tools::RELATION_1_2);

  MMFF::Tools::addAngles(*mol,mmffMolProperties, field);
  TEST_ASSERT(field->contribs().size() == 1);
  MMFF::Tools::addVdW(*mol, cid, mmffMolProperties, field, nbrMat);
  TEST_ASSERT(field->contribs().size() == 1);
  MMFF::Tools::addTorsions(*mol, mmffMolProperties, field);
  TEST_ASSERT(field->contribs().size() == 1);

  
  mol2 = MolOps::addHs(*mol);
  TEST_ASSERT(mol2->getNumAtoms() == 6);
  delete field;
  delete mmffMolProperties;
  
  mmffMolProperties = new MMFF::MMFFMolProperties(*mol2);
  TEST_ASSERT(mmffMolProperties);
  TEST_ASSERT(mmffMolProperties->isValid());

  field = new ForceFields::ForceField();
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

void testMMFFBuilder2()
{
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing MMFF builder+minimization." << std::endl;

  RWMol *mol;
  int needMore;

  ForceFields::ForceField *field;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/MMFF/test_data";
  mol = MolFileToMol(pathName + "/small1.mol",false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);

  field = MMFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  needMore = field->minimize();
  TEST_ASSERT(!needMore);
  //std::cout << MolToMolBlock(mol);
  delete mol;
  delete field;

  mol = MolFileToMol(pathName+"/small2.mol",false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);

  field = MMFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  needMore = field->minimize(100);
  TEST_ASSERT(!needMore);
  //std::cout << MolToMolBlock(mol);
  delete mol;
  delete field;

  mol = MolFileToMol(pathName + "/benzene.mol",false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);

  field = MMFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  needMore = field->minimize();
  TEST_ASSERT(!needMore);
  //std::cout << MolToMolBlock(mol);
  delete mol;
  delete field;
  
  mol = MolFileToMol(pathName + "/toluene.mol",false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);

  field = MMFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  needMore = field->minimize();
  TEST_ASSERT(!needMore);
  //std::cout << MolToMolBlock(mol);
  delete mol;
  delete field;

  mol = MolFileToMol(pathName + "/complex1.mol",false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);

  field = MMFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  needMore = field->minimize();
  TEST_ASSERT(!needMore);
  //std::cout << MolToMolBlock(mol);
  delete mol;
  delete field;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}


void testMMFFBatch()
{
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing bulk MMFF (regression to check that things run)." << std::endl;

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
      field = 0;
    }
    if (field) {
      field->initialize();
      int failed = field->minimize(500);
      if (failed) {
        BOOST_LOG(rdErrorLog) << " not converged (code = " << failed << ")" << std::endl;
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


void testIssue239()
{
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing Issue239." << std::endl;

  RWMol *mol;
  int needMore;
  ForceFields::ForceField *field;
  double e1, e2;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/MMFF/test_data";
  mol = MolFileToMol(pathName + "/Issue239.mol",false);
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

void testIssue242()
{
// FIX: Changes to the forcefield (connected to Issue 408) have
// made it so that this particular problem no longer manifests
// in this molecule/embedding. A new test case is needed.
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing Issue242." << std::endl;

  RWMol *mol, *mol2;
  int needMore;
  ForceFields::ForceField *field = 0, *field2 = 0;
  std::string mb1,mb2;
  double e1,e2;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/MMFF/test_data";

  mol = MolFileToMol(pathName + "/Issue242-2.mol");
  TEST_ASSERT(mol);

  mol2 = MolFileToMol(pathName + "/Issue242-2.mol");
  TEST_ASSERT(mol2);

  TEST_ASSERT(DGeomHelpers::EmbedMolecule(*mol2, 30, 2370) >= 0);
  mb1 = MolToMolBlock(*mol);
  mb2 = MolToMolBlock(*mol2);

  //std::cout << mb1 << "------\n";
  //std::cout << mb2 << "------\n";

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
  //TEST_ASSERT(feq(e2,e1,0.1));

  needMore = field->minimize(200,1.0e-4);
  needMore = field2->minimize(200,1.0e-4);
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

void testGithub308()
{
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing github 308: crash during MMFF parameterization ." << std::endl;
  ROMol *mol = SmilesToMol("FF");
  TEST_ASSERT(DGeomHelpers::EmbedMolecule(*mol) >= 0);
  int needMore;
  ForceFields::ForceField *field = 0;
  TEST_ASSERT(mol);
  MMFF::MMFFMolProperties mmffMolProperties(*mol);
  TEST_ASSERT(mmffMolProperties.isValid());
  field = MMFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  needMore = field->minimize(200, 1.0e-6, 1.0e-3);
  TEST_ASSERT(!needMore);
}

void testSFIssue1653802()
{
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
  MMFF::MMFFMolProperties *mmffMolProperties = new MMFF::MMFFMolProperties(*mol);
  TEST_ASSERT(mmffMolProperties);

  boost::shared_array<boost::uint8_t> nbrMat;
  field = new ForceFields::ForceField();
  MMFF::Tools::addBonds(*mol, mmffMolProperties, field);
  TEST_ASSERT(field->contribs().size() == 8);

  nbrMat = MMFF::Tools::buildNeighborMatrix(*mol);
  MMFF::Tools::addAngles(*mol, mmffMolProperties, field);
  TEST_ASSERT(field->contribs().size() == 20);
  MMFF::Tools::addTorsions(*mol, mmffMolProperties, field);
  //std::cout << field->contribs().size() << std::endl;
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

void testSFIssue2378119()
{
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
    ForceFields::ForceField *field = MMFF::constructForceField(*mol, 100.0, -1, false);
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

void testGithub162()
{
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing github 162: Incorrect SMILES after MMFF parameterization ." << std::endl;
  {
    ROMol *mol = SmilesToMol("C1=CNC=C1");
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getAtomWithIdx(2)->getNumExplicitHs()==1);
    MMFF::MMFFMolProperties *mmffMolProperties = new MMFF::MMFFMolProperties(*mol);
    TEST_ASSERT(mmffMolProperties);
    TEST_ASSERT(mmffMolProperties->isValid());
    TEST_ASSERT(mol->getAtomWithIdx(2)->getNumExplicitHs()==1);
    
    delete mol;
    delete mmffMolProperties;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}



#ifdef RDK_TEST_MULTITHREADED
namespace {
  void runblock_mmff(const std::vector<ROMol *> &mols,const std::vector<double> &energies,
                unsigned int count,unsigned int idx){
    for(unsigned int rep=0;rep<500;++rep){
      for(unsigned int i=0;i<mols.size();++i){
        if(i%count != idx) continue;
        ROMol *mol = mols[i];
        ForceFields::ForceField *field = 0;
        if(!(rep%100)) BOOST_LOG(rdErrorLog) << "Rep: "<<rep<<" Mol:" << i << std::endl;
        try {
          field = MMFF::constructForceField(*mol);
        } catch (...) {
          field = 0;
        }
        TEST_ASSERT(field);
        field->initialize();
        int failed = field->minimize(500);
        TEST_ASSERT(!failed);
        double eng=field->calcEnergy();
        TEST_ASSERT(feq(eng,energies[i]));
        delete field;
      }
    }
  }
}
#include <boost/thread.hpp>  
void testMMFFMultiThread(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test MMFF multithreading" << std::endl;

  ForceFields::ForceField *field;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/MMFF/test_data";
  SDMolSupplier suppl(pathName + "/bulk.sdf");
  std::vector<ROMol *> mols;
  while(!suppl.atEnd()&&mols.size()<100){
    ROMol *mol=0;
    try{
      mol=suppl.next();
    } catch(...){
      continue;
    }
    if(!mol) continue;
    mols.push_back(mol);
  }

  std::cerr<<"generating reference data"<<std::endl;
  std::vector<double> energies(mols.size(),0.0);
  for(unsigned int i=0;i<mols.size();++i){
    ROMol mol(*mols[i]);
    ForceFields::ForceField *field = 0;
    try {
      field = MMFF::constructForceField(mol);
    } catch (...) {
      field = 0;
    }
    TEST_ASSERT(field);
    field->initialize();
    int failed = field->minimize(500);
    TEST_ASSERT(!failed);
    energies[i]=field->calcEnergy();
    delete field;
  }
  
  boost::thread_group tg;

  std::cerr<<"processing"<<std::endl;
  unsigned int count=4;
  for(unsigned int i=0;i<count;++i){
    std::cerr<<" launch :"<<i<<std::endl;std::cerr.flush();
    tg.add_thread(new boost::thread(runblock_mmff,mols,energies,count,i));
  }
  tg.join_all();
  
  BOOST_FOREACH(ROMol *mol,mols){
    delete mol;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}
#endif

void testGithub224()
{
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing github 224: crash during MMFF parameterization ." << std::endl;
  {
    ROMol *mol = SmilesToMol("[1*]C");
    TEST_ASSERT(mol);
    MMFF::MMFFMolProperties *mmffMolProperties = new MMFF::MMFFMolProperties(*mol);
    TEST_ASSERT(mmffMolProperties);
    TEST_ASSERT(!mmffMolProperties->isValid());
    
    delete mol;
    delete mmffMolProperties;
  }
  {
    ROMol *mol = SmilesToMol("[1*]c1ccc(S(=O)(=O)Nc2ccc([2*])cc2)cc1");
    TEST_ASSERT(mol);
    MMFF::MMFFMolProperties *mmffMolProperties = new MMFF::MMFFMolProperties(*mol);
    TEST_ASSERT(mmffMolProperties);
    TEST_ASSERT(!mmffMolProperties->isValid());
    
    delete mol;
    delete mmffMolProperties;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}




//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
int main()
{
  RDLog::InitLogs();
#if 1
  testMMFFTyper1();
  testMMFFBuilder1();
  testMMFFBuilder2();
  testIssue239();
  testIssue242();
  testGithub308();
  testSFIssue1653802();
  testSFIssue2378119();
  testMMFFBatch();
#ifdef RDK_TEST_MULTITHREADED
  testMMFFMultiThread();
#endif
  testGithub162();
#endif
  testGithub224();
}
