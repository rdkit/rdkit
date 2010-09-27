//  $Id$
// 
//   Copyright (C) 2002-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <RDGeneral/types.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SanitException.h>

using namespace RDKit;

int testArom(RWMol m) {
  //BOOST_LOG(rdErrorLog) << "********************************************\n";
  //BOOST_LOG(rdErrorLog) << "Testing aromaticity :\n";

  int res = MolOps::setAromaticity(m);

  //BOOST_LOG(rdErrorLog) << "\nDone testing aromaticity \n";
  //BOOST_LOG(rdErrorLog) << "********************************************\n";
  return res;
}

void testKekule(RWMol &m) {
  MolOps::Kekulize(m);
}

void getSmiles(std::string line, std::string &smi) {
  int x = line.find("\t");
  if (x == -1) {
    x = line.find(" ");
  }

  smi = line.substr(0,x);
  std::string rem = line.substr(x+1, line.length() - x);
  x = rem.find(" ");
  std::string name;
  name = rem.substr(0, x);
  
  BOOST_LOG(rdInfoLog) << smi << " " << name << " " ;
}

int main(int argc, char *argv[]) {
  RDLog::InitLogs();
  std::string fname;
  if (argc > 1) {
    fname = argv[1];
  }
  else {
    BOOST_LOG(rdErrorLog) << "Pass in the list of smiles\n";
  }
  
  std::ifstream inStream(fname.c_str());
  const int MAX_LINE_LEN = 256;
  char inLine[MAX_LINE_LEN];
  std::string tmpstr;
  std::string smi;
  inStream.getline(inLine, MAX_LINE_LEN,'\n');
  tmpstr = inLine;
  //MolOps molop;
  while (tmpstr.size() > 0) {
    getSmiles(tmpstr, smi);
    RWMol *m = SmilesToMol(smi, 0, 0);
    
    try {
      //testKekule(*m);
      
      MolOps::sanitizeMol(*m);

      int nar;
      m->getProp("numArom", nar);
      BOOST_LOG(rdInfoLog)<< nar << "\n";

      //MolOps::setHybridization(*m);
      delete m;
    }
    catch (MolSanitizeException){
      BOOST_LOG(rdErrorLog) << smi << "\n";
      delete m;
    }
    inStream.getline(inLine, MAX_LINE_LEN,'\n');
    tmpstr = inLine;
  }
  return 1;
}


