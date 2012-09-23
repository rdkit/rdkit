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

#include <GraphMol/ForceFieldHelpers/UFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/UFF/Builder.h>
#include <ForceField/ForceField.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>

using namespace RDKit;
#if 1
void testUFFTyper1(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test UFF atom labels." << std::endl;

  ROMol *mol;
  std::string key;

  mol = SmilesToMol("[SiH3]CC(=O)NC");
  TEST_ASSERT(mol);

  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(0));
  TEST_ASSERT(key=="Si3");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(1));
  TEST_ASSERT(key=="C_3");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(2));
  TEST_ASSERT(key=="C_R");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(3));
  TEST_ASSERT(key=="O_R");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(4));
  TEST_ASSERT(key=="N_R");

  delete mol;
  mol = SmilesToMol("CC(=O)C");
  TEST_ASSERT(mol);

  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(0));
  TEST_ASSERT(key=="C_3");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(1));
  TEST_ASSERT(key=="C_2");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(2));
  TEST_ASSERT(key=="O_2");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(3));
  TEST_ASSERT(key=="C_3");
  
 
  delete mol;
  mol = SmilesToMol("C(=O)S");
  TEST_ASSERT(mol);

  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(0));
  TEST_ASSERT(key=="C_2");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(1));
  TEST_ASSERT(key=="O_2");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(2));
  TEST_ASSERT(key=="S_3+2");

  delete mol;
  mol = SmilesToMol("SCS(=O)S(=O)(=O)O");
  TEST_ASSERT(mol);

  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(0));
  TEST_ASSERT(key=="S_3+2");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(1));
  TEST_ASSERT(key=="C_3");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(2));
  TEST_ASSERT(key=="S_3+4");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(4));
  TEST_ASSERT(key=="S_3+6");
  
  delete mol;
  mol = SmilesToMol("PCP(O)CP(=O)(=O)");
  TEST_ASSERT(mol);

  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(0));
  TEST_ASSERT(key=="P_3+3");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(1));
  TEST_ASSERT(key=="C_3");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(2));
  TEST_ASSERT(key=="P_3+3");
  // FIX: I am not 100% convinced that this should be SP2
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(5));
  TEST_ASSERT(key=="P_3+5");
  

  delete mol;
  mol = SmilesToMol("C(F)(Cl)(Br)I");
  TEST_ASSERT(mol);

  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(0));
  TEST_ASSERT(key=="C_3");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(1));
  TEST_ASSERT(key=="F_");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(2));
  TEST_ASSERT(key=="Cl");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(3));
  TEST_ASSERT(key=="Br");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(4));
  TEST_ASSERT(key=="I_");
  
  delete mol;
  mol = SmilesToMol("[Li].[Na].[K].[Rb].[Cs]");
  TEST_ASSERT(mol);

  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(0));
  TEST_ASSERT(key=="Li");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(1));
  TEST_ASSERT(key=="Na");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(2));
  TEST_ASSERT(key=="K_");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(3));
  TEST_ASSERT(key=="Rb");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(4));
  TEST_ASSERT(key=="Cs");
  

  
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testUFFTyper2(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test UFF atom typer." << std::endl;

  ROMol *mol,*mol2;
  std::string key;

  mol = SmilesToMol("[SiH3]CC(=O)NC");
  TEST_ASSERT(mol);

  UFF::AtomicParamVect types;
  bool foundAll;
  boost::tie(types,foundAll)=UFF::getAtomTypes(*mol);
  TEST_ASSERT(foundAll);
  TEST_ASSERT(types.size()==mol->getNumAtoms());
  for(UFF::AtomicParamVect::const_iterator it=types.begin();
      it!=types.end();it++){
    TEST_ASSERT((*it));
  }
  mol2 = MolOps::addHs(*mol);
  delete mol;
  boost::tie(types,foundAll)=UFF::getAtomTypes(*mol2);
  TEST_ASSERT(foundAll);
  TEST_ASSERT(types.size()==mol2->getNumAtoms());
  for(UFF::AtomicParamVect::const_iterator it=types.begin();
      it!=types.end();it++){
    TEST_ASSERT((*it));
  }
  
  // connected with sf.net bug 2094445 
  mol = SmilesToMol("[SiH2]=C");
  TEST_ASSERT(mol);
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(0));
  TEST_ASSERT(key=="Si3");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(1));
  TEST_ASSERT(key=="C_2");

  mol = SmilesToMol("[AlH]=C");
  TEST_ASSERT(mol);
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(0));
  TEST_ASSERT(key=="Al3");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(1));
  TEST_ASSERT(key=="C_2");

  mol = SmilesToMol("[Mg]=C");
  TEST_ASSERT(mol);
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(0));
  TEST_ASSERT(key=="Mg3+2");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(1));
  TEST_ASSERT(key=="C_2");

  mol = SmilesToMol("[SiH3][Si]([SiH3])=C");
  TEST_ASSERT(mol);
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(0));
  TEST_ASSERT(key=="Si3");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(1));
  TEST_ASSERT(key=="Si3");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(2));
  TEST_ASSERT(key=="Si3");
  key = UFF::Tools::getAtomLabel(mol->getAtomWithIdx(3));
  TEST_ASSERT(key=="C_2");


  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}


void testUFFBuilder1(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing UFF builder tools." << std::endl;

  ROMol *mol,*mol2;
  std::string key;

  UFF::AtomicParamVect types;
  bool foundAll;
  ForceFields::ForceField *field;
  boost::shared_array<int> nbrMat;

  mol = SmilesToMol("CC(O)C");
  Conformer *conf = new Conformer(mol->getNumAtoms());
  int cid = static_cast<int>(mol->addConformer(conf, true));
  TEST_ASSERT(mol);
  boost::tie(types,foundAll)=UFF::getAtomTypes(*mol);
  TEST_ASSERT(foundAll);
  TEST_ASSERT(types.size()==mol->getNumAtoms());
  field=new ForceFields::ForceField();
  UFF::Tools::addBonds(*mol,types,field);

  TEST_ASSERT(field->contribs().size()==3);

  nbrMat = UFF::Tools::buildNeighborMatrix(*mol);
  TEST_ASSERT(nbrMat[0]==-2);
  TEST_ASSERT(nbrMat[1]==-1);
  TEST_ASSERT(nbrMat[2]==1);
  TEST_ASSERT(nbrMat[3]==1);

  UFF::Tools::addAngles(*mol,types,field,nbrMat);
  TEST_ASSERT(field->contribs().size()==6);

  // there are no non-bonded terms here:
  UFF::Tools::addNonbonded(*mol,cid,types,field,nbrMat);
  TEST_ASSERT(field->contribs().size()==6);

  // and no torsions either, until we add Hs:
  UFF::Tools::addTorsions(*mol,types,field);
  TEST_ASSERT(field->contribs().size()==6);

  delete mol;
  delete field;
  mol = SmilesToMol("CCOC");
  Conformer *conf2 = new Conformer(mol->getNumAtoms());
  cid = static_cast<int>(mol->addConformer(conf2, true));
  TEST_ASSERT(mol);
  boost::tie(types,foundAll)=UFF::getAtomTypes(*mol);
  TEST_ASSERT(foundAll);
  TEST_ASSERT(types.size()==mol->getNumAtoms());
  field=new ForceFields::ForceField();
  UFF::Tools::addBonds(*mol,types,field);

  TEST_ASSERT(field->contribs().size()==3);

  nbrMat = UFF::Tools::buildNeighborMatrix(*mol);
  TEST_ASSERT(nbrMat[0]==-2);
  TEST_ASSERT(nbrMat[1]==-1);
  TEST_ASSERT(nbrMat[2]==1);
  TEST_ASSERT(nbrMat[3]==-2);

  UFF::Tools::addAngles(*mol,types,field,nbrMat);
  TEST_ASSERT(field->contribs().size()==5);
  UFF::Tools::addNonbonded(*mol,cid,types,field,nbrMat);
  TEST_ASSERT(field->contribs().size()==6);
  UFF::Tools::addTorsions(*mol,types,field);
  TEST_ASSERT(field->contribs().size()==7);



  delete mol;
  delete field;
  mol = SmilesToMol("CO");
  Conformer *conf3 = new Conformer(mol->getNumAtoms());
  cid = static_cast<int>(mol->addConformer(conf3, true));
  TEST_ASSERT(mol);
  boost::tie(types,foundAll)=UFF::getAtomTypes(*mol);
  TEST_ASSERT(foundAll);
  TEST_ASSERT(types.size()==mol->getNumAtoms());

  field=new ForceFields::ForceField();
  UFF::Tools::addBonds(*mol,types,field);

  TEST_ASSERT(field->contribs().size()==1);

  nbrMat = UFF::Tools::buildNeighborMatrix(*mol);
  TEST_ASSERT(nbrMat[0]==-2);
  TEST_ASSERT(nbrMat[1]==-1);

  UFF::Tools::addAngles(*mol,types,field,nbrMat);
  TEST_ASSERT(field->contribs().size()==1);
  UFF::Tools::addNonbonded(*mol,cid,types,field,nbrMat);
  TEST_ASSERT(field->contribs().size()==1);
  UFF::Tools::addTorsions(*mol,types,field);
  TEST_ASSERT(field->contribs().size()==1);

  
  mol2 = MolOps::addHs(*mol);
  TEST_ASSERT(mol2->getNumAtoms()==6);
  delete field;
  
  boost::tie(types,foundAll)=UFF::getAtomTypes(*mol2);
  TEST_ASSERT(foundAll);
  TEST_ASSERT(types.size()==mol2->getNumAtoms());

  field=new ForceFields::ForceField();
  UFF::Tools::addBonds(*mol2,types,field);
  TEST_ASSERT(field->contribs().size()==5);

  nbrMat = UFF::Tools::buildNeighborMatrix(*mol2);
  UFF::Tools::addAngles(*mol2,types,field,nbrMat);
  TEST_ASSERT(field->contribs().size()==12);
  UFF::Tools::addNonbonded(*mol2,cid,types,field,nbrMat);
  TEST_ASSERT(field->contribs().size()==15);
  UFF::Tools::addTorsions(*mol2,types,field);
  TEST_ASSERT(field->contribs().size()==18);
  delete mol2;

  delete mol;
  delete field;
  
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}


void testUFFBuilder2(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing UFF builder+minimization." << std::endl;

  RWMol *mol;
  std::string key;
  int needMore;

  ForceFields::ForceField *field;

  std::string pathName=getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/test_data";
  mol = MolFileToMol(pathName+"/small1.mol",false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);

  field=UFF::constructForceField(*mol);
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

  field=UFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  needMore = field->minimize(100);
  TEST_ASSERT(!needMore);
  //std::cout << MolToMolBlock(mol);
  delete mol;
  delete field;

  mol = MolFileToMol(pathName+"/benzene.mol",false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);

  field=UFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  needMore = field->minimize();
  TEST_ASSERT(!needMore);
  //std::cout << MolToMolBlock(mol);
  delete mol;
  delete field;
  
  mol = MolFileToMol(pathName+"/toluene.mol",false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);

  field=UFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  needMore = field->minimize();
  TEST_ASSERT(!needMore);
  //std::cout << MolToMolBlock(mol);
  delete mol;
  delete field;

  mol = MolFileToMol(pathName+"/complex1.mol",false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);

  field=UFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  needMore = field->minimize();
  TEST_ASSERT(!needMore);
  //std::cout << MolToMolBlock(mol);
  delete mol;
  delete field;

  

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}


void testUFFBatch(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing bulk UFF (regression to check that things run)." << std::endl;

  ROMol *mol;
  std::string key;

  ForceFields::ForceField *field;

  std::string pathName=getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/test_data";
  SDMolSupplier suppl(pathName+"/bulk.sdf",false);

  int count=0;
  mol = suppl.next();
  while(mol && !suppl.atEnd() ){
    count++;
    MolOps::sanitizeMol(*(RWMol *)mol);
    std::string origMolBlock = MolToMolBlock(*mol);

    BOOST_LOG(rdErrorLog) << "Mol:" << count << std::endl;
    try {
      field=UFF::constructForceField(*mol);
    } catch (...) {
      field = 0;
    }
    if(field){
      field->initialize();
      int failed=field->minimize(500);
      if(failed){
        BOOST_LOG(rdErrorLog) << " not converged (code="<<failed<<")" << std::endl;
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


void testUFFBuilderSpecialCases(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing UFF special cases." << std::endl;

  RWMol *mol;
  std::string key;
  int needMore;
  RDGeom::Point3D v1,v2;
  ForceFields::ForceField *field;

  std::string pathName=getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/test_data";
  // ----------
  //  Trigonal bipyramid
  // ----------
  mol = MolFileToMol(pathName+"/tbp.mol",false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);

  const Conformer &conf = mol->getConformer();
  field=UFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  needMore = field->minimize(200,1e-8,1e-4);
  TEST_ASSERT(!needMore);
  v1 = conf.getAtomPos(0).directionVector(conf.getAtomPos(1));
  v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(2));
  TEST_ASSERT(feq(v1.dotProduct(v2),-1.0,1e-3));
  v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(3));
  TEST_ASSERT(feq(v1.dotProduct(v2),0.0,1e-3));
  v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(4));
  TEST_ASSERT(feq(v1.dotProduct(v2),0.0,1e-3));
  v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(5));
  TEST_ASSERT(feq(v1.dotProduct(v2),0.0,1e-3));

  v1 = conf.getAtomPos(0).directionVector(conf.getAtomPos(2));
  v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(3));
  TEST_ASSERT(feq(v1.dotProduct(v2),0.0,1e-3));
  v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(4));
  TEST_ASSERT(feq(v1.dotProduct(v2),0.0,1e-3));
  v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(5));
  TEST_ASSERT(feq(v1.dotProduct(v2),0.0,1e-3));

  v1 = conf.getAtomPos(0).directionVector(conf.getAtomPos(3));
  v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(4));
  TEST_ASSERT(feq(v1.dotProduct(v2),-0.5,1e-3));
  v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(5));
  TEST_ASSERT(feq(v1.dotProduct(v2),-0.5,1e-3));

  v1 = conf.getAtomPos(0).directionVector(conf.getAtomPos(4));
  v2 = conf.getAtomPos(0).directionVector(conf.getAtomPos(5));
  TEST_ASSERT(feq(v1.dotProduct(v2),-0.5,1e-3));

  
  delete mol;
  delete field;
  

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}


void testIssue239(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing Issue239." << std::endl;

  RWMol *mol;
  int needMore;
  ForceFields::ForceField *field;
  double e1,e2;

  std::string pathName=getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/test_data";
  mol = MolFileToMol(pathName+"/Issue239.mol",false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);

  field=UFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  needMore = field->minimize(200,1e-6,1e-3);
  e1 = field->calcEnergy();
  needMore = field->minimize(200,1e-6,1e-3);
  e2 = field->calcEnergy();
  TEST_ASSERT(feq(e2,e1,0.1));
  
  delete mol;
  delete field;
  

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}
#endif

void testIssue242(){
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
  pathName += "/Code/GraphMol/ForceFieldHelpers/test_data";

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

void testSFIssue1653802(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing SFIssue1653802." << std::endl;

  RWMol *mol;
  int needMore;
  ForceFields::ForceField *field;

  std::string pathName=getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/test_data";

  mol = MolFileToMol(pathName+"/cyclobutadiene.mol",false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);


  UFF::AtomicParamVect types;
  bool foundAll;
  boost::shared_array<int> nbrMat;
  boost::tie(types,foundAll)=UFF::getAtomTypes(*mol);
  TEST_ASSERT(foundAll);
  TEST_ASSERT(types.size()==mol->getNumAtoms());
  field=new ForceFields::ForceField();
  UFF::Tools::addBonds(*mol,types,field);
  TEST_ASSERT(field->contribs().size()==8);

  nbrMat = UFF::Tools::buildNeighborMatrix(*mol);
  UFF::Tools::addAngles(*mol,types,field,nbrMat);
  TEST_ASSERT(field->contribs().size()==18);
  UFF::Tools::addTorsions(*mol,types,field);
  //std::cout << field->contribs().size() << std::endl;
  TEST_ASSERT(field->contribs().size()==34);
  UFF::Tools::addNonbonded(*mol,0,types,field,nbrMat);
  delete field;

  field = UFF::constructForceField(*mol);
  field->initialize();
  needMore = field->minimize(200,1e-6,1e-3);
  TEST_ASSERT(!needMore);
  
  delete mol;
  delete field;
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testSFIssue2378119(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing SFIssue2378119." << std::endl;

  std::string pathName=getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/test_data";
  {
    RWMol *mol = MolFileToMol(pathName+"/Issue2378119.mol");
    TEST_ASSERT(mol);
    ForceFields::ForceField *field=UFF::constructForceField(*mol);
    TEST_ASSERT(field);
    field->initialize();
    double e1 = field->calcEnergy();
    TEST_ASSERT(e1>0.0 && e1<1e12);

    int needMore = field->minimize(200,1e-6,1e-3);
    TEST_ASSERT(!needMore);
    double e2 = field->calcEnergy();
    TEST_ASSERT(e2<e1);

    delete mol;
    delete field;
  }
  {
    RWMol *mol = MolFileToMol(pathName+"/Issue2378119.2.mol");
    TEST_ASSERT(mol);
    ForceFields::ForceField *field=UFF::constructForceField(*mol);
    TEST_ASSERT(field);
    field->initialize();
    double e1 = field->calcEnergy();
    TEST_ASSERT(e1==0.0);

    int needMore = field->minimize(200,1e-6,1e-3);
    TEST_ASSERT(!needMore);
    double e2 = field->calcEnergy();
    TEST_ASSERT(e2==e1);

    delete mol;
    delete field;
  }
  {
    RWMol *mol = MolFileToMol(pathName+"/Issue2378119.2.mol");
    TEST_ASSERT(mol);
    ForceFields::ForceField *field=UFF::constructForceField(*mol,100.0,-1,false);
    TEST_ASSERT(field);
    field->initialize();
    double e1 = field->calcEnergy();
    TEST_ASSERT(e1>0.0 && e1<1e12);

    int needMore = field->minimize(200,1e-6,1e-3);
    TEST_ASSERT(!needMore);
    double e2 = field->calcEnergy();
    TEST_ASSERT(e2<e1);

    delete mol;
    delete field;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}



void testMissingParams(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing handling missing parameters." << std::endl;

  {
    UFF::AtomicParamVect types;
    bool foundAll;

    RWMol *mol = SmilesToMol("[Cu](C)(C)(C)(C)C");
    TEST_ASSERT(mol);

    ROMol *mol2 = MolOps::addHs(*mol);
    delete mol;

    TEST_ASSERT(DGeomHelpers::EmbedMolecule(*mol2)>=0);
    
    boost::tie(types,foundAll)=UFF::getAtomTypes(*mol2);
    TEST_ASSERT(!foundAll);
    TEST_ASSERT(types.size()==mol2->getNumAtoms());
    TEST_ASSERT(!types[0]);

    // make sure we can optimize anyway:
    ForceFields::ForceField *field=UFF::constructForceField(*mol2,types);
    TEST_ASSERT(field);
    field->initialize();
    double e1=field->calcEnergy();
    int needMore = field->minimize();
    TEST_ASSERT(needMore);
    double e2 = field->calcEnergy();
    TEST_ASSERT(e2<e1);
    //std::cerr<<" DE: "<<e1<<" -> "<<e2<<std::endl;
    delete mol2;
    delete field;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}
void testSFIssue3009337(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing SFIssue3009337." << std::endl;

  std::string pathName=getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/test_data";
  {
    RWMol *mol = MolFileToMol(pathName+"/Issue3009337.mol",true,false);
    TEST_ASSERT(mol);
    ForceFields::ForceField *field=UFF::constructForceField(*mol);
    TEST_ASSERT(field);
    field->initialize();
    double e1 = field->calcEnergy();
    TEST_ASSERT(e1>0.0 && e1<1e12);

    int needMore = field->minimize(200,1e-6,1e-3);
    TEST_ASSERT(!needMore);
    double e2 = field->calcEnergy();
    TEST_ASSERT(e2<e1);

    delete mol;
    delete field;
  }
  {
    RWMol *mol = MolFileToMol(pathName+"/Issue3009337.2.mol",true,false);
    TEST_ASSERT(mol);
    ForceFields::ForceField *field=UFF::constructForceField(*mol);
    TEST_ASSERT(field);
    field->initialize();
    double e1 = field->calcEnergy();
    TEST_ASSERT(e1>0.0 && e1<1e12);

    int needMore = field->minimize(200,1e-6,1e-3);
    TEST_ASSERT(!needMore);
    double e2 = field->calcEnergy();
    TEST_ASSERT(e2<e1);

    delete mol;
    delete field;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
int main(){
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
#endif
  testMissingParams();
  testSFIssue3009337();

}
