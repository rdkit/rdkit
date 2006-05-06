// $Id$
//
//  Copyright (C) 2001-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
//

// std bits
#include <iostream>


// vflib bits
#include <argedit.h>
#include <argraph.h>
#include <vf_sub_state.h>
#include <match.h>

// RD bits
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <SmilesParse/SmilesParse.h>
#include <CDXMLParse/ParseCDXML.h>

// local
#include "SubstructMatch.h"

using namespace RDKit;

int main(int argc,char *argv[])
{
  std::string m1Name(argv[1]);
  std::string m2Smi(argv[2]);

  XML4CInit();
  Mol::MOL_VECT *molList=CDXML2Mol(m1Name);
  XML4CTerminate();
  Mol *m1 = &((*molList)[0]);
  
  Mol *m2 = SmilesToMol(m2Smi);

  MatchVectType matchV;
  if(SubstructMatch(*m2,*m1,matchV)){
    std::cout << "Got a match: " << std::endl;
    MatchVectType::iterator i;
    for(i=matchV.begin();i!=matchV.end();i++){
      std::cout << "\t" << i->first << " -> " << i->second << std::endl;
    }
    
  } else {
    std::cout << "No match" << std::endl;
  }
  return 0;
}

  
