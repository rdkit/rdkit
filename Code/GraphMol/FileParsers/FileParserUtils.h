//
//  Copyright (C) 2010 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_FILEPARSERUTILS_H
#define _RD_FILEPARSERUTILS_H

#include <string>
#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

namespace RDKit{
  class RWMol;
  class Conformer;

  namespace FileParserUtils {
    template <typename T>
    T stripSpacesAndCast(const std::string &input,bool acceptSpaces=false){
      std::string trimmed=boost::trim_copy(input);
      if(acceptSpaces && trimmed==""){
	return 0;
      } else {
	return boost::lexical_cast<T>(trimmed);
      }
    }
    int toInt(const std::string &input,bool acceptSpaces=false);
    double toDouble(const std::string &input,bool acceptSpaces=true);


    // reads a line from an MDL v3K CTAB
    std::string getV3000Line(std::istream *inStream,unsigned int &line);

    // nAtoms and nBonds are ignored on input, set on output
    bool ParseV3000CTAB(std::istream *inStream,unsigned int &line,
			RWMol *mol, Conformer *&conf,
			bool &chiralityPossible,
			unsigned int &nAtoms,unsigned int &nBonds,
                        bool strictParsing=true,bool expectMEND=true);

    // nAtoms and nBonds are used
    bool ParseV2000CTAB(std::istream *inStream,unsigned int &line,
			RWMol *mol, Conformer *&conf,
			bool &chiralityPossible,
			unsigned int &nAtoms,unsigned int &nBonds,
                        bool strictParsing=true);

    Atom *replaceAtomWithQueryAtom(RWMol *mol,Atom *atom);

  }
}

#endif


















