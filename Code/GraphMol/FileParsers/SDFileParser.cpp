// $Id$
//
//  Copyright (C) 2002-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#include "FileParsers.h"

#include <RDGeneral/StreamOps.h>
#include <fstream>

namespace RDKit{
    
  //*************************************
  //
  // Every effort has been made to adhere to MDL's standard
  // for SD files
  //  
  //*************************************
  
  RWMOL_SPTR_VECT SDDataStreamToMols(std::istream *inStream,bool sanitize,bool removeHs){
    unsigned int line = 0;
    std::string tempStr;

    RWMOL_SPTR_VECT res;
    res.clear();

    while(1){
      // read a molecule
      RWMOL_SPTR nextMol = RWMOL_SPTR(MolDataStreamToMol(inStream, line, sanitize, removeHs));
      if(!nextMol) break;

      // now read the data
      tempStr = getLine(inStream);
      
      // FIX: report files missing the $$$$ marker
      while(!(inStream->eof()) && tempStr.find("$$$$") != 0){
	// parse the data itself
	tempStr = getLine(inStream);
	
      }      

      res.push_back(nextMol);
    }
    return res;
  };

  RWMOL_SPTR_VECT SDDataStreamToMols(std::istream &inStream,bool sanitize,bool removeHs){
    return SDDataStreamToMols(&inStream,sanitize,removeHs);
  };

  
  RWMOL_SPTR_VECT SDFileToMols(std::string fName,bool sanitize,bool removeHs){
    std::ifstream inStream(fName.c_str());
    
    // FIX: we find a better way of reporting errors than this.  The
    // rest of the code already deals with the stream not existing, so
    // it's just a matter of reporting.
    RWMOL_SPTR_VECT res;
    res.clear();
    if(inStream){
      res = SDDataStreamToMols(inStream,sanitize,removeHs);
    }

    return res;
  };
  
    
}









