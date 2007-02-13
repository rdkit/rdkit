// $Id$
//
//  Copyright (C) 2006 Greg Landrum
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include "MolCatalog.h"
#include "MolCatalogEntry.h"
#include "MolCatalogParams.h"

#include <stdlib.h>
#include <iostream>
#include <fstream>

#include <GraphMol/FileParsers/MolSupplier.h>
#include <RDGeneral/types.h>
#include <RDGeneral/FileParseException.h>
#include <DataStructs/BitVects.h>

using namespace RDKit;

void test1(){
  BOOST_LOG(rdInfoLog) << ">>------------- Test 1" << std::endl;

  // MolCatalogParams are currently unused, so testing is easy:
  MolCatalogParams *mparams = new MolCatalogParams();
  std::string pkl = mparams->Serialize();
  TEST_ASSERT(pkl=="");

  MolCatalog *mcat = new MolCatalog(mparams);
  TEST_ASSERT(mcat->getNumEntries()==0);
  TEST_ASSERT(mcat->getFPLength()==0);

  MolCatalogEntry *entry;
  std::string smi;
  RWMol *mol;

  mol = SmilesToMol("c1ccc(O)cc1C(=O)O");
  entry = new MolCatalogEntry(mol);
  TEST_ASSERT(entry->getMol()==mol);
  entry->setDescription("top");
  TEST_ASSERT(entry->getDescription()=="top");
  entry->setOrder(0);
  TEST_ASSERT(entry->getOrder()==0);


  mcat->addEntry(entry);
  TEST_ASSERT(mcat->getNumEntries()==1);
  TEST_ASSERT(mcat->getFPLength()==1);
  TEST_ASSERT(mcat->getEntryWithBitId(0)==entry);
  
  mol = SmilesToMol("O");
  entry = new MolCatalogEntry(mol);
  entry->setDescription("child1");
  entry->setOrder(1);
  mcat->addEntry(entry);
  mol = SmilesToMol("C(=O)O");
  entry = new MolCatalogEntry(mol);
  entry->setDescription("child2");
  entry->setOrder(1);
  mcat->addEntry(entry);

  TEST_ASSERT(mcat->getNumEntries()==3);
  TEST_ASSERT(mcat->getFPLength()==3);

  TEST_ASSERT(mcat->getEntryWithBitId(0)->getMol()->getNumAtoms()==10);
  TEST_ASSERT(mcat->getEntryWithBitId(1)->getMol()->getNumAtoms()==1);
  TEST_ASSERT(mcat->getEntryWithBitId(2)->getMol()->getNumAtoms()==3);
  
  TEST_ASSERT(mcat->getDownEntryList(0).size()==0);
  TEST_ASSERT(mcat->getDownEntryList(1).size()==0);
  TEST_ASSERT(mcat->getDownEntryList(2).size()==0);

  TEST_ASSERT(mcat->getEntriesOfOrder(0).size()==1);
  TEST_ASSERT(mcat->getEntriesOfOrder(1).size()==2);
  TEST_ASSERT(mcat->getEntriesOfOrder(2).size()==0);
 
  mcat->addEdge(0,1);
  mcat->addEdge(0,2);

  TEST_ASSERT(mcat->getDownEntryList(0).size()==2);
  TEST_ASSERT(mcat->getDownEntryList(1).size()==0);
  TEST_ASSERT(mcat->getDownEntryList(2).size()==0);
  
  TEST_ASSERT(mcat->getEntriesOfOrder(0).size()==1);
  TEST_ASSERT(mcat->getEntriesOfOrder(1).size()==2);
  TEST_ASSERT(mcat->getEntriesOfOrder(2).size()==0);

  pkl = mcat->Serialize();
  delete mcat;
  mcat = new MolCatalog(pkl);
  TEST_ASSERT(mcat->getNumEntries()==3);
  TEST_ASSERT(mcat->getFPLength()==3);
  TEST_ASSERT(mcat->getDownEntryList(0).size()==2);
  TEST_ASSERT(mcat->getDownEntryList(1).size()==0);
  TEST_ASSERT(mcat->getDownEntryList(2).size()==0);
  
  TEST_ASSERT(mcat->getEntriesOfOrder(0).size()==1);
  TEST_ASSERT(mcat->getEntriesOfOrder(1).size()==2);
  TEST_ASSERT(mcat->getEntriesOfOrder(2).size()==0);

  TEST_ASSERT(mcat->getEntryWithBitId(0)->getMol()->getNumAtoms()==10);
  TEST_ASSERT(mcat->getEntryWithBitId(1)->getMol()->getNumAtoms()==1);
  TEST_ASSERT(mcat->getEntryWithBitId(2)->getMol()->getNumAtoms()==3);
  
  
  BOOST_LOG(rdInfoLog) << "<<-------------- Done" << std::endl;
}

int main() {
  RDLog::InitLogs();
  test1();
  return 0;

}


    
