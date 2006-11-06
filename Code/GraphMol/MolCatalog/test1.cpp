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

  MolCatalogParams *mparams = new MolCatalogParams();
  MolCatalog mcat(mparams);


  BOOST_LOG(rdInfoLog) << "<<-------------- Done" << std::endl;
}

int main() {
  RDLog::InitLogs();
  test1();
  return 0;

}


    
