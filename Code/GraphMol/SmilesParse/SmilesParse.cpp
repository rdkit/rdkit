// $Id$
//
//  Copyright (C) 2001-2006 Greg Landrum and Rational Discovery LLC
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
#include "InputFiller.h"
#include <RDGeneral/RDLog.h>
#include <RDGeneral/Invariant.h>
#include <boost/algorithm/string.hpp>
int yysmiles_parse (void);
extern int yysmiles_debug; 

int yysmarts_parse (void);
extern int yysmarts_debug; 


namespace RDKit{
  namespace SmilesParse{
    std::vector<RDKit::RWMol *> molList_g;
  }

  RWMol *toMol(std::string inp,int func(void)){
    RWMol *res; 
    if(SmilesParse::molList_g.size()){
      for(std::vector<RDKit::RWMol *>::iterator iter=SmilesParse::molList_g.begin();
          iter!=SmilesParse::molList_g.end();++iter){
        delete *iter;
      }
      SmilesParse::molList_g.clear();
      SmilesParse::molList_g.resize(0);
    }
    setInputCharPtr( (char *)inp.c_str());
    try {
      func();
      if(SmilesParse::molList_g.size()<=0){
        res = 0;
      } else {
        res = SmilesParse::molList_g[0];
        SmilesParse::molList_g.resize(0);
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
    charPtrCleanup();

    return res;
  }
  
RWMol *SmilesToMol(std::string smi,int debugParse,bool sanitize){
  yysmiles_debug = debugParse;
  // strip any leading/trailing whitespace:
  boost::trim_if(smi,boost::is_any_of(" \t\r\n"));
  RWMol *res = toMol(smi,yysmiles_parse);
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
  RWMol *res = toMol(sma,yysmarts_parse);
  if(res && mergeHs){
    ROMol *tmp = MolOps::mergeQueryHs(*res);
    delete res;
    res = static_cast<RWMol *>(tmp);
  }
  return res;
};
}
