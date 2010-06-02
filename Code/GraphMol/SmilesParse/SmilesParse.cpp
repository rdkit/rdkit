// $Id$
//
//  Copyright (C) 2001-2010 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
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
int yysmiles_parse (std::vector<RDKit::RWMol *>*,void *);
int yysmiles_lex_init (void **);
int yysmiles_lex_destroy (void *);
void setup_smiles_string(const std::string &text,void *);
extern int yysmiles_debug; 

int yysmarts_parse (std::vector<RDKit::RWMol *>*,void *);
int yysmarts_lex_init (void **);
int yysmarts_lex_destroy (void *);
void setup_smarts_string(const std::string &text,void *);
extern int yysmarts_debug; 

int smiles_parse(const std::string &inp,
    std::vector<RDKit::RWMol *> &molVect){
  void *scanner;
  TEST_ASSERT(!yysmiles_lex_init(&scanner));
  setup_smiles_string(inp,scanner);
  int res=yysmiles_parse(&molVect,scanner);
  yysmiles_lex_destroy(scanner);
  return res;
}

int smarts_parse(const std::string &inp,
    std::vector<RDKit::RWMol *> &molVect){
  void *scanner;
  TEST_ASSERT(!yysmarts_lex_init(&scanner));
  setup_smarts_string(inp,scanner);
  int res=yysmarts_parse(&molVect,scanner);
  yysmarts_lex_destroy(scanner);
  return res;
}

namespace RDKit{
RWMol *toMol(std::string inp,int func(const std::string &,
    std::vector<RDKit::RWMol *> &)){
  RWMol *res;
  std::vector<RDKit::RWMol *> molVect;
  try {
    func(inp,molVect);
    if(molVect.size()<=0){
      res = 0;
    } else {
      res = molVect[0];
      molVect[0]=0;
      SmilesParseOps::CloseMolRings(res,false);
      SmilesParseOps::AdjustAtomChiralityFlags(res);
      // No sense leaving this bookmark intact:
      if(res->hasAtomBookmark(ci_RIGHTMOST_ATOM)){
        res->clearAtomBookmark(ci_RIGHTMOST_ATOM);
      }
    }
  } catch (SmilesParseException &e) {
    BOOST_LOG(rdErrorLog) << e.message() << std::endl;
    res = 0;
  }
  BOOST_FOREACH(RDKit::RWMol *molPtr,molVect){
    if(molPtr) delete molPtr;
  }

  return res;
}

RWMol *SmilesToMol(std::string smi,int debugParse,bool sanitize){
  yysmiles_debug = debugParse;
  // strip any leading/trailing whitespace:
  boost::trim_if(smi,boost::is_any_of(" \t\r\n"));
  RWMol *res = toMol(smi,smiles_parse);
  if(sanitize && res){
    // we're going to remove explicit Hs from the graph,
    // this triggers a sanitization, so we do not need to
    // worry about doing one here:
    ROMol *tmp = MolOps::removeHs(*res,false,false);
    delete res;
    res = static_cast<RWMol *>(tmp);
    // figure out stereochemistry:
    MolOps::assignStereochemistry(*res,true);
  }

  return res;
};
RWMol *SmartsToMol(std::string sma,int debugParse,bool mergeHs){
  yysmarts_debug = debugParse;
  boost::trim_if(sma,boost::is_any_of(" \t\r\n"));
  RWMol *res = toMol(sma,smarts_parse);
  if(res && mergeHs){
    ROMol *tmp = MolOps::mergeQueryHs(*res);
    delete res;
    res = static_cast<RWMol *>(tmp);
  }
  return res;
};
}
