// $Id$
//
//  Copyright (c) 2007-2014, Novartis Institutes for BioMedical Research Inc.
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
#include <GraphMol/FileParsers/FileParserUtils.h>
#include <sstream>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/FileParseException.h>

#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/tokenizer.hpp>


namespace RDKit {

  namespace {
    void ParseV2000RxnBlock(std::istream &inStream,unsigned int &line,
			    ChemicalReaction *&rxn)
    {
      std::string tempStr;
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
      rxn=new ChemicalReaction();
  
      unsigned int nReacts=0,nProds=0,nAgents=0;
      unsigned int spos = 0;
      if(tempStr.size()<6){
	throw ChemicalReactionParserException("rxn counts line is too short");
      }
      try {
	nReacts = FileParserUtils::stripSpacesAndCast<unsigned int>(tempStr.substr(0,3));
	spos = 3;
	nProds = FileParserUtils::stripSpacesAndCast<unsigned int>(tempStr.substr(spos,3));
	spos = 6;;
	if(tempStr.size()>6){
	  std::string trimmed=boost::trim_copy(tempStr.substr(spos,3));
	  if(trimmed.size() > 0){
	    nAgents = FileParserUtils::stripSpacesAndCast<unsigned int>(tempStr.substr(spos,3));
	    spos = 9;
	  }
	}
      } catch (boost::bad_lexical_cast &) {
        delete rxn;
        rxn=0;
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
	rxn->addReactantTemplate(ROMOL_SPTR(react));
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
	rxn->addProductTemplate(ROMOL_SPTR(prod));
      }

      for(unsigned int i=0;i<nAgents;++i){
	line++;
	tempStr = getLine(inStream);
	if(inStream.eof()){
	  throw ChemicalReactionParserException("premature EOF hit.");
	}
	if(tempStr.substr(0,4)!="$MOL"){
	  throw ChemicalReactionParserException("$MOL header not found");
	}
	ROMol *agent;
	try{
		agent=MolDataStreamToMol(inStream,line,false);
	} catch (FileParseException &e){
	  std::ostringstream errout;
	  errout << "Cannot parse agent " << i << ". The error was:\n\t" << e.message();
	  throw ChemicalReactionParserException(errout.str()) ;
	}
	rxn->addAgentTemplate(ROMOL_SPTR(agent));
      }

    }


    void ParseV3000RxnBlock(std::istream &inStream,unsigned int &line,
			    ChemicalReaction *&rxn)
    {
      std::string tempStr;
  
      // skip the header block:
      line++;
      tempStr = getLine(inStream);
      line++;
      tempStr = getLine(inStream);
      line++;
      tempStr = getLine(inStream);
      line++;

      rxn=new ChemicalReaction();
  
      tempStr = FileParserUtils::getV3000Line(&inStream,line);
      boost::to_upper(tempStr);
      tempStr = boost::trim_copy(tempStr);
      std::vector<std::string> tokens;
      boost::split(tokens,tempStr,boost::is_any_of(" \t"),
		   boost::token_compress_on);
      if(tokens.size()<3 || tokens[0]!="COUNTS" ){
	throw ChemicalReactionParserException("bad counts line");
      }
      unsigned int nReacts=FileParserUtils::stripSpacesAndCast<unsigned int>(tokens[1]);
      unsigned int nProds=FileParserUtils::stripSpacesAndCast<unsigned int>(tokens[2]);
      unsigned int nAgents=0;
      if(tokens.size()>3){
        nAgents=FileParserUtils::stripSpacesAndCast<unsigned int>(tokens[3]);
      }
      
      tempStr = FileParserUtils::getV3000Line(&inStream,line);
      boost::to_upper(tempStr);
      if(tempStr.length()<14 || tempStr.substr(0,14) != "BEGIN REACTANT"){
	throw FileParseException("BEGIN REACTANT line not found") ;
      }
      for(unsigned int i=0;i<nReacts;++i){
	RWMol *react;
	unsigned int natoms,nbonds;
	bool chiralityPossible;
	Conformer *conf=0;
	react= new RWMol();
	try {
	  FileParserUtils::ParseV3000CTAB(&inStream,line,react,conf,
					  chiralityPossible,natoms,nbonds,true,false);
	} catch (FileParseException &e){
	  delete react;
	  react=0;
	  std::ostringstream errout;
	  errout << "Cannot parse reactant " << i << ". The error was:\n\t" << e.message();
	  throw ChemicalReactionParserException(errout.str()) ;
	}
	if(!react){
	  throw ChemicalReactionParserException("Null reactant in reaction file.");
	}
	rxn->addReactantTemplate(ROMOL_SPTR(dynamic_cast<ROMol *>(react)));
      }
      tempStr = FileParserUtils::getV3000Line(&inStream,line);
      boost::to_upper(tempStr);
      if(tempStr.length()<12 || tempStr.substr(0,12) != "END REACTANT"){
	throw FileParseException("END REACTANT line not found") ;
      }
      tempStr = FileParserUtils::getV3000Line(&inStream,line);
      boost::to_upper(tempStr);
      if(tempStr.length()<13 || tempStr.substr(0,13) != "BEGIN PRODUCT"){
	throw FileParseException("BEGIN PRODUCT line not found") ;
      }
      for(unsigned int i=0;i<nProds;++i){
	RWMol *prod;
	unsigned int natoms,nbonds;
	bool chiralityPossible;
	Conformer *conf=0;
	prod= new RWMol();
	try {
	  FileParserUtils::ParseV3000CTAB(&inStream,line,prod,conf,
					  chiralityPossible,natoms,nbonds,true,false);
	} catch (FileParseException &e){
	  delete prod;
	  prod=0;
	  std::ostringstream errout;
	  errout << "Cannot parse product " << i << ". The error was:\n\t" << e.message();
	  throw ChemicalReactionParserException(errout.str()) ;
	}
	if(!prod){
	  throw ChemicalReactionParserException("Null product in reaction file.");
	}
	rxn->addProductTemplate(ROMOL_SPTR(dynamic_cast<ROMol *>(prod)));
      }
      tempStr = FileParserUtils::getV3000Line(&inStream,line);
      boost::to_upper(tempStr);
      if(tempStr.length()<11 || tempStr.substr(0,11) != "END PRODUCT"){
	throw FileParseException("END PRODUCT line not found") ;
      }

      if(nAgents){
        tempStr = FileParserUtils::getV3000Line(&inStream,line);
        boost::to_upper(tempStr);
        if(tempStr.length()<14 || tempStr.substr(0,14) != "BEGIN AGENT"){
          throw FileParseException("BEGIN AGENT line not found") ;
        }
      }
      for(unsigned int i=0;i<nAgents;++i){
        RWMol *agent;
        unsigned int natoms,nbonds;
        bool chiralityPossible;
        Conformer *conf=0;
        agent= new RWMol();
        try {
          FileParserUtils::ParseV3000CTAB(&inStream,line,agent,conf,
        		  chiralityPossible,natoms,nbonds,true,false);
        } catch (FileParseException &e){
          delete agent;
          agent=0;
          std::ostringstream errout;
          errout << "Cannot parse agent " << i << ". The error was:\n\t" << e.message();
          throw ChemicalReactionParserException(errout.str()) ;
        }
        rxn->addAgentTemplate(ROMOL_SPTR(dynamic_cast<ROMol *>(agent)));
      }
      if(nAgents){
        tempStr = FileParserUtils::getV3000Line(&inStream,line);
        boost::to_upper(tempStr);
        if(tempStr.length()<12 || tempStr.substr(0,12) != "END AGENT"){
    	    throw FileParseException("END AGENT line not found") ;
        }
      }
    }
  } // end of local namespace

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
    int version=2000;
    if(tempStr.size()>=10&&tempStr.substr(5,5)=="V3000") version=3000;

    ChemicalReaction *res=0;
    try{
      if(version==2000){
	ParseV2000RxnBlock(inStream,line,res);
      } else {
	ParseV3000RxnBlock(inStream,line,res);
      }
    }
    catch (ChemicalReactionParserException &e) { 
      // catch our exceptions and throw them back after cleanup
      delete res;
      res=0;
      throw e;
    }
    // convert atoms to queries:
    for(MOL_SPTR_VECT::const_iterator iter=res->beginReactantTemplates();
	iter != res->endReactantTemplates();++iter){
      // to write the mol block, we need ring information:
      for(ROMol::AtomIterator atomIt=(*iter)->beginAtoms();
          atomIt!=(*iter)->endAtoms();++atomIt){
        FileParserUtils::replaceAtomWithQueryAtom((RWMol *)iter->get(),(*atomIt));
      }
    }
    for(MOL_SPTR_VECT::const_iterator iter=res->beginProductTemplates();
	iter != res->endProductTemplates();++iter){
      // to write the mol block, we need ring information:
      for(ROMol::AtomIterator atomIt=(*iter)->beginAtoms();
          atomIt!=(*iter)->endAtoms();++atomIt){
        FileParserUtils::replaceAtomWithQueryAtom((RWMol *)iter->get(),(*atomIt));
      }
    }

    
    // RXN-based reactions do not have implicit properties
    res->setImplicitPropertiesFlag(false);
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
