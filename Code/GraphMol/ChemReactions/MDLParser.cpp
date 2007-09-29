// $Id$
//
//  Copyright (c) 2007, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met: 
//
//     * Redistributions of source code must retain the above copyright 
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following 
//       disclaimer in the documentation and/or other materials provided 
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
//       nor the names of its contributors may be used to endorse or promote 
//       products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <sstream>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/FileParseException.h>

#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>


namespace RDKit {
  // it's kind of stinky that we have to do this, but as of g++3.2 and
  // boost 1.30, on linux calls to lexical_cast<int>(std::string)
  // crash if the string starts with spaces.
  template <typename T>
  T stripSpacesAndCast(const std::string &input,bool acceptSpaces=false){
    T res;
    std::string trimmed=boost::trim_copy(input);
    if(acceptSpaces && trimmed==""){
      return 0;
    } else {
      return boost::lexical_cast<T>(trimmed);
    }
    // if this fails we throw a boost::bad_lexical_cast exception:
    res = boost::lexical_cast<T>(trimmed);
    return res;
  }


  //! Parse a text stream in MDL rxn format into a ChemicalReaction 
  ChemicalReaction * RxnDataStreamToChemicalReaction(std::istream &inStream,unsigned int &line) {
    std::string tempStr;

    // header line
    line++;
    tempStr = getLine(inStream);
    if(inStream.eof()){
      throw ChemicalReactionParserException("premature EOF hit.");
    }
    if(tempStr.substr(0,4)!="$RXN"){
      throw ChemicalReactionParserException("$RXN header not found");
    }
    
    // FIX: parse name and comment fields
    line++;
    tempStr = getLine(inStream);
    line++;
    tempStr = getLine(inStream);
    line++;
    tempStr = getLine(inStream);
    line++;
    tempStr = getLine(inStream);
    if(inStream.eof()){
      throw ChemicalReactionParserException("premature EOF hit.");
    }
  
    ChemicalReaction *res=new ChemicalReaction();
  
    unsigned int nReacts=0,nProds=0;
    unsigned int spos = 0;
    if(tempStr.size()<6){
      throw ChemicalReactionParserException("rxn counts line is too short");
    }
    try {
      // it *sucks* that the lexical_cast stuff above doesn't work on linux        
      nReacts = stripSpacesAndCast<unsigned int>(tempStr.substr(0,3));
      spos = 3;
      nProds = stripSpacesAndCast<unsigned int>(tempStr.substr(spos,3));
      spos = 6;
    } catch (boost::bad_lexical_cast &) {
      if (res) {
        delete res;
      }
      std::ostringstream errout;
      errout << "Cannot convert " << tempStr.substr(spos,3) << " to int";
      throw ChemicalReactionParserException(errout.str()) ;
    }
    
    for(unsigned int i=0;i<nReacts;++i){
      line++;
      tempStr = getLine(inStream);
      if(inStream.eof()){
        throw ChemicalReactionParserException("premature EOF hit.");
      }
      if(tempStr.substr(0,4)!="$MOL"){
        throw ChemicalReactionParserException("$MOL header not found");
      }
      ROMol *react;
      try {
        react=MolDataStreamToMol(inStream,line,false);
      } catch (FileParseException &e){
        std::ostringstream errout;
        errout << "Cannot parse reactant " << i << ". The error was:\n\t" << e.message();
        throw ChemicalReactionParserException(errout.str()) ;
      }
      if(!react){
        throw ChemicalReactionParserException("Null reactant in reaction file.");
      }
      res->addReactantTemplate(ROMOL_SPTR(react));
    }
    for(unsigned int i=0;i<nProds;++i){
      line++;
      tempStr = getLine(inStream);
      if(inStream.eof()){
        throw ChemicalReactionParserException("premature EOF hit.");
      }
      if(tempStr.substr(0,4)!="$MOL"){
        throw ChemicalReactionParserException("$MOL header not found");
      }
      ROMol *prod;
      try{
        prod=MolDataStreamToMol(inStream,line,false);
      } catch (FileParseException &e){
        std::ostringstream errout;
        errout << "Cannot parse product " << i << ". The error was:\n\t" << e.message();
        throw ChemicalReactionParserException(errout.str()) ;
      }
      if(!prod){
        throw ChemicalReactionParserException("Null product in reaction file.");
      }
      res->addProductTemplate(ROMOL_SPTR(prod));
    }
    return res;  
  };

  ChemicalReaction * RxnBlockToChemicalReaction(const std::string &rxnBlock) {
    std::istringstream inStream(rxnBlock);
    unsigned int line = 0;
    return RxnDataStreamToChemicalReaction(inStream, line);
  };

  ChemicalReaction * RxnFileToChemicalReaction(const std::string &fName) {
    std::ifstream inStream(fName.c_str());
    if(!inStream){
      return NULL;
    }
    ChemicalReaction *res=NULL;
    if(!inStream.eof()){
      unsigned int line = 0;
      res=RxnDataStreamToChemicalReaction(inStream, line);
    }
    return res;
  
  };


}
