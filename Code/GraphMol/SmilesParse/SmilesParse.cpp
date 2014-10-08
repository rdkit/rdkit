// $Id$
//
//  Copyright (C) 2001-2014 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// ----------------------------------------------------------------------------------
//  Despite the name of this file, both SMILES and SMARTS parsers are exposed here
//
//  General comments about the parsers:
//   - Atom numbering will be preserved, so input order of atoms==internal order
//
//   - Bond ordering is not, in general, preserved.  Specifically, ring closure
//     bonds will occur at the end of the bond list in general.  Basically ring
//     closure bonds are not constructed until fragments are closed.  This forces
//     some form of reordering.
//     
//
//
#include <GraphMol/RDKitBase.h>
#include "SmilesParse.h"
#include "SmilesParseOps.h"
#include <RDGeneral/RDLog.h>
#include <RDGeneral/Invariant.h>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <list>

int yysmiles_parse (const char *,std::vector<RDKit::RWMol *>*,std::list<unsigned int> *,void *);
int yysmiles_lex_init (void **);
int yysmiles_lex_destroy (void *);
void setup_smiles_string(const std::string &text,void *);
extern int yysmiles_debug; 

int yysmarts_parse (const char *,std::vector<RDKit::RWMol *>*,void *);
int yysmarts_lex_init (void **);
int yysmarts_lex_destroy (void *);
void setup_smarts_string(const std::string &text,void *);
extern int yysmarts_debug; 
namespace RDKit{
  namespace {
    int smiles_parse(const std::string &inp,
		     std::vector<RDKit::RWMol *> &molVect){
      std::list<unsigned int> branchPoints;
      void *scanner;
      int res;

      TEST_ASSERT(!yysmiles_lex_init(&scanner));
      try {
        setup_smiles_string(inp,scanner);
        res = yysmiles_parse(inp.c_str(),&molVect,&branchPoints,scanner);
      } catch(...) {
        yysmiles_lex_destroy(scanner);
        throw;
      }
      yysmiles_lex_destroy(scanner);
      if(!branchPoints.empty()){
        throw SmilesParseException("extra open parentheses");
      }
      return res;
    }
    int smarts_parse(const std::string &inp,
		     std::vector<RDKit::RWMol *> &molVect){
      void *scanner;
      int res;
      TEST_ASSERT(!yysmarts_lex_init(&scanner));
      try {
        setup_smarts_string(inp,scanner);
        res=yysmarts_parse(inp.c_str(),&molVect,scanner);
      } catch(...) {
        yysmarts_lex_destroy(scanner);
        throw;
      }
      yysmarts_lex_destroy(scanner);
      return res;
    }

    typedef enum {
      BASE=0,
      BRANCH,
      RECURSE
    } SmaState;

    std::string labelRecursivePatterns(std::string sma){
#ifndef NO_AUTOMATIC_SMARTS_RELABELLING
      std::list<SmaState> state;
      std::list<unsigned int> startRecurse;
      std::map<std::string, std::string> patterns;
      std::string res="";

      state.push_back(BASE);

      unsigned int pos=0;
      while(pos<sma.size()){
	res += sma[pos];
	if(sma[pos]=='$' && pos+1<sma.size() && sma[pos+1]=='('){
	  state.push_back(RECURSE);
	  startRecurse.push_back(pos);
	  ++pos;
	  res += sma[pos];
	} else if(sma[pos]=='('){
	  state.push_back(BRANCH);
	} else if(sma[pos]==')'){
	  SmaState currState=state.back();
	  state.pop_back();
	  if(currState==RECURSE){
	    unsigned int dollarPos=startRecurse.back();
	    startRecurse.pop_back();
	    if(pos+1>=sma.size() || sma[pos+1] !='_'){
	      std::string recurs = sma.substr(dollarPos,pos-dollarPos+1);
	      std::string label;
	      if(patterns.find(recurs)!=patterns.end()){
		// seen this one before, add the label
		label=patterns[recurs];
	      } else {
		label=boost::lexical_cast<std::string>(patterns.size()+100);
		patterns[recurs]=label;
	      }
	      res += "_" + label;
	    }
	  } else if(currState==BRANCH) {
	    // no need to do anything here.
	  }
	}
	++pos;
      }
      //std::cerr<< " >"<<sma<<"->"<<res<<std::endl;
      return res;
#else
      return sma;
#endif
    }
  } // end of local namespace

  RWMol *toMol(std::string inp,int func(const std::string &,
					std::vector<RDKit::RWMol *> &),
               std::string origInp){
    RWMol *res = 0;
    std::vector<RDKit::RWMol *> molVect;
    try {
      func(inp,molVect);
      if(molVect.size()>0){
	res = molVect[0];
	SmilesParseOps::CloseMolRings(res,false);
	SmilesParseOps::AdjustAtomChiralityFlags(res);
	// No sense leaving this bookmark intact:
	if(res->hasAtomBookmark(ci_RIGHTMOST_ATOM)){
	  res->clearAtomBookmark(ci_RIGHTMOST_ATOM);
	}
        molVect[0]=0; // NOTE: to avoid leaks on failures, this should occur last in this if.
      }
    } catch (SmilesParseException &e) {
      std::string nm="SMILES";
      if(func==smarts_parse){
        nm="SMARTS";
        
      }
      BOOST_LOG(rdErrorLog) << nm<<" Parse Error: "<< e.message() << " for input: "<< origInp << std::endl;
      res = 0;
    }
    BOOST_FOREACH(RDKit::RWMol *molPtr,molVect){
      if (molPtr) {
        // Clean-up the bond bookmarks when not calling CloseMolRings
        SmilesParseOps::CleanupAfterParseError(molPtr);
        delete molPtr;
      }
    }

    return res;
  }

  RWMol *SmilesToMol(std::string smi,int debugParse,bool sanitize,
                     std::map<std::string,std::string> *replacements){
    yysmiles_debug = debugParse;
    // strip any leading/trailing whitespace:
    boost::trim_if(smi,boost::is_any_of(" \t\r\n"));

    if(replacements){
      bool loopAgain=true;
      while(loopAgain){
        loopAgain=false;
        for(std::map<std::string, std::string>::const_iterator replIt=replacements->begin();
            replIt!=replacements->end();++replIt){
          if(boost::find_first(smi,replIt->first)){
            loopAgain=true;
            boost::replace_all(smi,replIt->first,replIt->second);
          }
        }
      }
    }

    RWMol *res = toMol(smi,smiles_parse,smi);
    if(sanitize && res){
      // we're going to remove explicit Hs from the graph,
      // this triggers a sanitization, so we do not need to
      // worry about doing one here:
      try {
        MolOps::removeHs(*res,false,false);
        // figure out stereochemistry:
        MolOps::assignStereochemistry(*res,true,true,true);
      } catch (...) {
        delete res;
        throw;
      }
    }

    return res;
  };
  RWMol *SmartsToMol(std::string sma,int debugParse,bool mergeHs,
                     std::map<std::string, std::string> *replacements){
    yysmarts_debug = debugParse;
    boost::trim_if(sma,boost::is_any_of(" \t\r\n"));
    if(replacements){
      bool loopAgain=true;
      while(loopAgain){
        loopAgain=false;
        for(std::map<std::string, std::string>::const_iterator replIt=replacements->begin();
            replIt!=replacements->end();++replIt){
          if(boost::find_first(sma,replIt->first)){
            loopAgain=true;
            boost::replace_all(sma,replIt->first,replIt->second);
          }
        }
      }
    }
    std::string oInput=sma;
    sma=labelRecursivePatterns(sma);

    RWMol *res = toMol(sma,smarts_parse,oInput);
    if(res && mergeHs){
      try {
        MolOps::mergeQueryHs(*res);
      } catch(...) {
        delete res;
        throw;
      }
    }
    return res;
  };
}
