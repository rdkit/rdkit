//
//  Copyright (C) 2001-2016 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// ----------------------------------------------------------------------------------
//  Despite the name of this file, both SMILES and SMARTS parsers are exposed
//  here
//
//  General comments about the parsers:
//   - Atom numbering will be preserved, so input order of atoms==internal order
//
//   - Bond ordering is not, in general, preserved.  Specifically, ring closure
//     bonds will occur at the end of the bond list in general.  Basically ring
//     closure bonds are not constructed until fragments are closed.  This
//     forces
//     some form of reordering.
//
//
//
#include <RDGeneral/BoostStartInclude.h>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <GraphMol/RDKitBase.h>
#include "SmilesParse.h"
#include "SmilesParseOps.h"
#include <RDGeneral/RDLog.h>
#include <RDGeneral/Invariant.h>
#include <list>

int yysmiles_parse(const char *, std::vector<RDKit::RWMol *> *,
                   std::list<unsigned int> *, void *);
int yysmiles_lex_init(void **);
int yysmiles_lex_destroy(void *);
size_t setup_smiles_string(const std::string &text, void *);
extern int yysmiles_debug;

int yysmarts_parse(const char *, std::vector<RDKit::RWMol *> *, void *);
int yysmarts_lex_init(void **);
int yysmarts_lex_destroy(void *);
size_t setup_smarts_string(const std::string &text, void *);
extern int yysmarts_debug;
namespace RDKit {
namespace {
int smiles_parse(const std::string &inp, std::vector<RDKit::RWMol *> &molVect) {
  std::list<unsigned int> branchPoints;
  void *scanner;
  int res;

  TEST_ASSERT(!yysmiles_lex_init(&scanner));
  try {
    size_t ltrim = setup_smiles_string(inp, scanner);
    res = yysmiles_parse(inp.c_str() + ltrim, &molVect, &branchPoints, scanner);
  } catch (...) {
    yysmiles_lex_destroy(scanner);
    throw;
  }
  yysmiles_lex_destroy(scanner);
  if (!branchPoints.empty()) {
    throw SmilesParseException("extra open parentheses");
  }
  return res;
}
int smarts_parse(const std::string &inp, std::vector<RDKit::RWMol *> &molVect) {
  void *scanner;
  int res;
  TEST_ASSERT(!yysmarts_lex_init(&scanner));
  try {
    size_t ltrim = setup_smarts_string(inp, scanner);
    res = yysmarts_parse(inp.c_str() + ltrim, &molVect, scanner);
  } catch (...) {
    yysmarts_lex_destroy(scanner);
    throw;
  }
  yysmarts_lex_destroy(scanner);
  return res;
}

typedef enum { BASE = 0, BRANCH, RECURSE } SmaState;

std::string labelRecursivePatterns(const std::string &sma) {
#ifndef NO_AUTOMATIC_SMARTS_RELABELLING
  std::list<SmaState> state;
  std::list<unsigned int> startRecurse;
  std::map<std::string, std::string> patterns;
  std::string res = "";

  state.push_back(BASE);

  unsigned int pos = 0;
  while (pos < sma.size()) {
    res += sma[pos];
    if (sma[pos] == '$' && pos + 1 < sma.size() && sma[pos + 1] == '(') {
      state.push_back(RECURSE);
      startRecurse.push_back(pos);
      ++pos;
      res += sma[pos];
    } else if (sma[pos] == '(') {
      state.push_back(BRANCH);
    } else if (sma[pos] == ')') {
      SmaState currState = state.back();
      state.pop_back();
      if (currState == RECURSE) {
        unsigned int dollarPos = startRecurse.back();
        startRecurse.pop_back();
        if (pos + 1 >= sma.size() || sma[pos + 1] != '_') {
          std::string recurs = sma.substr(dollarPos, pos - dollarPos + 1);
          std::string label;
          if (patterns.find(recurs) != patterns.end()) {
            // seen this one before, add the label
            label = patterns[recurs];
          } else {
            label = boost::lexical_cast<std::string>(patterns.size() + 100);
            patterns[recurs] = label;
          }
          res += "_" + label;
        }
      } else if (currState == BRANCH) {
        // no need to do anything here.
      }
    }
    ++pos;
  }
  // std::cerr<< " >"<<sma<<"->"<<res<<std::endl;
  return res;
#else
  return sma;
#endif
}
}  // end of local namespace

RWMol *toMol(const std::string &inp,
             int func(const std::string &, std::vector<RDKit::RWMol *> &),
             const std::string &origInp) {
  // empty strings produce empty molecules:
  if (inp == "") return new RWMol();
  RWMol *res = 0;
  std::vector<RDKit::RWMol *> molVect;
  try {
    func(inp, molVect);
    if (molVect.size() > 0) {
      res = molVect[0];
      SmilesParseOps::CloseMolRings(res, false);
      SmilesParseOps::AdjustAtomChiralityFlags(res);
      // No sense leaving this bookmark intact:
      if (res->hasAtomBookmark(ci_RIGHTMOST_ATOM)) {
        res->clearAtomBookmark(ci_RIGHTMOST_ATOM);
      }
      SmilesParseOps::CleanupAfterParsing(res);
      molVect[0] = 0;  // NOTE: to avoid leaks on failures, this should occur
                       // last in this if.
    }
  } catch (SmilesParseException &e) {
    std::string nm = "SMILES";
    if (func == smarts_parse) {
      nm = "SMARTS";
    }
    BOOST_LOG(rdErrorLog) << nm << " Parse Error: " << e.message()
                          << " for input: '" << origInp << "'" << std::endl;
    res = 0;
  }
  BOOST_FOREACH (RDKit::RWMol *molPtr, molVect) {
    if (molPtr) {
      // Clean-up the bond bookmarks when not calling CloseMolRings
      SmilesParseOps::CleanupAfterParseError(molPtr);
      delete molPtr;
    }
  }

  return res;
}

RWMol *SmilesToMol(const std::string &smiles, const SmilesParserParams &params) {
  yysmiles_debug = params.debugParse;

  std::string lsmiles = "", name = "", cxPart = "";
  if(params.parseName && !params.allowCXSMILES){
    std::vector<std::string> tokens;
    boost::split(tokens, smiles, boost::is_any_of(" \t"),
                 boost::token_compress_on);
    lsmiles = tokens[0];
    if(tokens.size()>1) name = tokens[1];
  } else if (params.allowCXSMILES) {
    size_t sidx = smiles.find_first_of(" \t");
    if (sidx != std::string::npos && sidx != 0)  {
      lsmiles = smiles.substr(0, sidx);
      cxPart = boost::trim_copy(smiles.substr(sidx, smiles.size() - sidx));
    }
  }

  if(lsmiles=="") {
    lsmiles = smiles;
  }
  // strip any leading/trailing whitespace:
  // boost::trim_if(smi,boost::is_any_of(" \t\r\n"));
  RWMol *res=NULL;
  if (params.replacements) {
    std::string smi = lsmiles;
    bool loopAgain = true;
    while (loopAgain) {
      loopAgain = false;
      for (std::map<std::string, std::string>::const_iterator replIt =
               params.replacements->begin();
           replIt != params.replacements->end(); ++replIt) {
        if (boost::find_first(smi, replIt->first)) {
          loopAgain = true;
          boost::replace_all(smi, replIt->first, replIt->second);
        }
      }
    }
    res = toMol(smi, smiles_parse, smi);
  } else {
    res = toMol(lsmiles, smiles_parse, lsmiles);
  }
  if ( res && (params.sanitize || params.removeHs)) {
    try {
      if(params.removeHs) {
        bool implicitOnly=false,updateExplicitCount=true;
        MolOps::removeHs(*res, implicitOnly, updateExplicitCount,
          params.sanitize);
      } else if(params.sanitize) {
        MolOps::sanitizeMol(*res);
      }
    } catch (...) {
      delete res;
      throw;
    }
    // figure out stereochemistry:
    bool cleanIt = true, force = true, flagPossible = true;
    MolOps::assignStereochemistry(*res, cleanIt, force, flagPossible);
  }
  if (res && params.allowCXSMILES && cxPart != "") {
    // it's goofy that we have to do this, but c++0x doesn't seem to have
    // a way to get a const_iterator from a non-const std::string. In C++11
    // we could just use .cend()
    const std::string &cxcopy = cxPart;
    std::string::const_iterator pos;
    SmilesParseOps::parseCXExtensions(*res, cxcopy, pos);
    if (params.parseName && pos != cxcopy.end()) {
      std::string nmpart(pos, cxcopy.end());
      name = boost::trim_copy(nmpart);
    }
  }
  if (res && name != "") res->setProp(common_properties::_Name, name);
  return res;
};
RWMol *SmartsToMol(const std::string &smarts, int debugParse, bool mergeHs,
                   std::map<std::string, std::string> *replacements) {
  yysmarts_debug = debugParse;
  // boost::trim_if(sma,boost::is_any_of(" \t\r\n"));
  std::string sma;
  RWMol *res;

  if (replacements) {
    sma = smarts;

    bool loopAgain = true;
    while (loopAgain) {
      loopAgain = false;
      for (std::map<std::string, std::string>::const_iterator replIt =
               replacements->begin();
           replIt != replacements->end(); ++replIt) {
        if (boost::find_first(sma, replIt->first)) {
          loopAgain = true;
          boost::replace_all(sma, replIt->first, replIt->second);
        }
      }
    }
    std::string oInput = sma;
    res = toMol(labelRecursivePatterns(sma), smarts_parse, oInput);
  } else {
    res = toMol(labelRecursivePatterns(smarts), smarts_parse, smarts);
  }
  if (res && mergeHs) {
    try {
      MolOps::mergeQueryHs(*res);
    } catch (...) {
      delete res;
      throw;
    }
  }
  return res;
};
}
