// $Id$
//
//  Copyright (C) 2001-2006 Randal Henne, Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <iostream>
#include <RDGeneral/RDLog.h>
//#include <boost/log/functions.hpp>

#include "SmilesParse.h"
#include "SmilesWrite.h"

using namespace RDKit;
using namespace std;
typedef ROMol Mol;
int
main(int argc, char *argv[])
{
  RDLog::InitLogs();
  boost::logging::enable_logs("rdApp.debug");

  int i=0;
  Mol *mol;
	
  if(argc < 2){
	string smis[]={
		       "C",
		       "CC",
		       "C-C",
		       "C=C",
		       "[CH2+]C[CH+2]",
		       "C1CC=1",
		       "C=1CC1",
		       "CCC",
		       "C=C-O",
		       "C1CC1",
		       "C1NC1",
		       "C1=CC1",
		       "C1CCC1",
		       "CC(C)CC",
		       "CC(=O)O",
		       "C1C(=O)C1",
		       "C1C(N)C1",
		       "CC(O)C",
		       //"(OC)CCC",
		       "CC=(CO)C",
		       //"(OC)=CCC",
		       "OC=CCC",
		       "CC([O-])O",
		       "C1CC2C1CC2",
		       "Cl/C=C/Cl",
		       "Cl/C=C\\Cl",
		       "CCC=",
		       "Cl/C=C/Cl",
		       "Cl/C=C\\Cl",
 			   "EOS"};
	while( smis[i] != "EOS" ){

	  std::cout << "Doing: " << smis[i] << std::endl;

	  mol = SmilesToMol(smis[i]);
	  if (mol) {
	    std::cout << "\t" << mol->getNumAtoms() << " atoms" << std::endl;
	    //mol->debugMol(cout);
	    //std::cout << std::endl;

	    std::cout << "\tSMILES: ";
	    std::cout << MolToSmiles(*mol) << std::endl;
	  
	    delete mol;
	  }
	  i++;
	}
  } else {
	string smi;
	std::cout << "Doing: " << smi << std::endl;
	bool debugParse = false,doSmarts=false;
	int startP=1;
	while(startP<argc){
	  string arg(argv[startP]);
	  if(arg=="-d"){
	    debugParse = true;
	    startP++;
	  }else if(arg=="-sma"){
	    doSmarts = true;
	    startP++;
	  } else {
	    break;
	  }
	}
	while(startP<argc){
	  smi = argv[startP++];
	  if(doSmarts){
	    std::cout << "In SMARTS: " << smi << std::endl;
	    mol = SmartsToMol(smi,debugParse);
	  } else {
	    std::cout << "In SMILES: " << smi << std::endl;
	    mol = SmilesToMol(smi,debugParse);
	    std::cout << "Out SMILES: " << std::endl;
	    smi = MolToSmiles(*mol,true);
	    std::cout << "  " << smi << std::endl;
	  }
	}
  }	

	return 1;
}
