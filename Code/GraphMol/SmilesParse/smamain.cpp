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
#include "SmilesParse.h"
#include "SmilesWrite.h"
#include <RDGeneral/RDLog.h>

using namespace RDKit;
using namespace std;
typedef ROMol Mol;
int
main(int argc, char *argv[])
{
  RDLog::InitLogs();
  int i=0;
  Mol *mol;
	
  if(argc < 2){
    return 1;
  } else {

    string sma;
    bool debugParse = false;
    int startP=1;
    if(argc>2){
      string arg(argv[1]);
      if(arg=="-d"){
	debugParse = true;
	startP = 2;
      }
    }

    while(startP<argc){
      sma = argv[startP++];
      std::cout << "In SMARTS: " << sma << std::endl;
      mol = SmartsToMol(sma,debugParse);
      if(!mol){
	BOOST_LOG(rdErrorLog) << "FAILED PARSE: " << sma << std::endl;
      }
    }
  }	
  
  return 1;
}
