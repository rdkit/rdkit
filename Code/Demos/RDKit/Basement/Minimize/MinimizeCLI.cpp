// $Id$
//
//  Copyright (C) 2004 Rational Discovery LLC
//   All Rights Reserved
//
#include <iostream>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>

#include <GraphMol/ForceFieldHelpers/UFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/UFF/Builder.h>
#include <ForceField/ForceField.h>

using namespace RDKit;


void runMol(ROMol *mol,int checkEvery=10,bool verbose=true){
  ForceFields::ForceField *field;

  std::cout << MolToMolBlock(*mol) << "$$$$" << std::endl;

  try{
    field=UFF::constructForceField(*mol,2.5);
  } catch (...) {
    field=0;
  }
  if(field){
    field->initialize();
    int needMore=1;
    int nPasses=0;
    while(needMore){
#if 1
      needMore = field->minimize(checkEvery);
      if(verbose) std::cerr << "\t" << ++nPasses << std::endl;
#else
      needMore = field->minimize(1);
      std::cout << MolToMolBlock(mol) << "$$$$" << std::endl;
#endif
    }
    std::cout << MolToMolBlock(*mol) << "$$$$" << std::endl;
    delete field;
  } else {
    std::cerr << "failed";
  }
  
}

void runMolFile(std::string fileName,int checkEvery=10){
  RWMol *mol=MolFileToMol(fileName,false);
  TEST_ASSERT(mol);
  MolOps::sanitizeMol(*mol);

  ROMol *mol2=MolOps::addHs(*mol,false,true);

  runMol(mol2,checkEvery);

  delete mol;
  delete mol2;
}

void runSDFile(std::string fileName,int checkEvery=10){
  SDMolSupplier suppl(fileName,false);
  
  RWMol *mol;

  mol = (RWMol *)suppl.next();
  while(mol){
    std::string name;
    mol->getProp("_Name",name);
    std::cerr << "Mol: " << name << std::endl;
    try{
      MolOps::sanitizeMol(*mol);
    } catch (...) {
      std::cerr << " sanitization failed" << std::endl;
      delete mol;
      mol = 0;
    }
    if(mol){
      ROMol *mol2=MolOps::addHs(*mol,false,true);
      delete mol;
      runMol(mol2,checkEvery,false);
      delete mol2;
    }
    mol = (RWMol *)suppl.next();

  }

}

//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
int main(int argc,char *argv[]){
  PRECONDITION(argc>1,"bad arguments");
  std::string fileName=argv[1];
  int checkEvery=10;
  std::cerr << ">" << fileName<< " " << fileName.find(".sdf") << std::endl;
  if(fileName.find(".sdf")==std::string::npos){
    runMolFile(fileName,checkEvery);
  } else {
    runSDFile(fileName,checkEvery);
  }

  std::cerr << "done" << std::endl;
}
