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
#include "SmilesParse.h"
#include <RDGeneral/BoostStartInclude.h>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <GraphMol/RDKitBase.h>
#include "SmilesParseOps.h"
#include <RDGeneral/RDLog.h>
#include <RDGeneral/Invariant.h>
#include "smiles.tab.hpp"
// NOTE: this is a bit fragile since a lot of the #defines in smiles.tab.hpp
// could prevent the same #defines in smarts.tab.hpp from being read.
// Fortunately if there are actually any problems here, they will inevitably
// show up very quickly in the tests.
#include "smarts.tab.hpp"
#include <list>

int yysmiles_lex_init(void **);
int yysmiles_lex_destroy(void *);
size_t setup_smiles_string(const std::string &text, void *);
extern int yysmiles_debug;

int yysmarts_lex_init(void **);
int yysmarts_lex_destroy(void *);
size_t setup_smarts_string(const std::string &text, void *);
extern int yysmarts_debug;
namespace RDKit {
namespace {
int smarts_bond_parse(const std::string &inp, Bond *&bond) {
  void *scanner;
  int res;

  TEST_ASSERT(!yysmarts_lex_init(&scanner));
  try {
    size_t ltrim = setup_smarts_string(inp, scanner);
    int start_tok = static_cast<int>(START_BOND);
    std::vector<RWMol *> molVect;
    Atom *lastAtom = nullptr;
    res = yysmarts_parse(inp.c_str() + ltrim, &molVect, lastAtom, bond, scanner,
                         start_tok);
  } catch (...) {
    yysmarts_lex_destroy(scanner);
    throw;
  }
  yysmarts_lex_destroy(scanner);
  return res;
}
int smarts_atom_parse(const std::string &inp, Atom *&atom) {
  void *scanner;
  int res;

  TEST_ASSERT(!yysmarts_lex_init(&scanner));
  try {
    size_t ltrim = setup_smarts_string(inp, scanner);
    int start_tok = static_cast<int>(START_ATOM);
    std::vector<RWMol *> molVect;
    Bond *lastBond = nullptr;
    res = yysmarts_parse(inp.c_str() + ltrim, &molVect, atom, lastBond, scanner,
                         start_tok);
  } catch (...) {
    yysmarts_lex_destroy(scanner);
    throw;
  }
  yysmarts_lex_destroy(scanner);
  return res;
}

int smiles_bond_parse(const std::string &inp, Bond *&bond) {
  std::list<unsigned int> branchPoints;
  void *scanner;
  int res;

  TEST_ASSERT(!yysmiles_lex_init(&scanner));
  try {
    size_t ltrim = setup_smiles_string(inp, scanner);
    int start_tok = static_cast<int>(START_BOND);
    std::vector<RWMol *> molVect;
    Atom *lastAtom = nullptr;
    res = yysmiles_parse(inp.c_str() + ltrim, &molVect, lastAtom, bond,
                         &branchPoints, scanner, start_tok);
  } catch (...) {
    yysmiles_lex_destroy(scanner);
    throw;
  }
  yysmiles_lex_destroy(scanner);
  return res;
}
int smiles_atom_parse(const std::string &inp, Atom *&atom) {
  std::list<unsigned int> branchPoints;
  void *scanner;
  int res;

  TEST_ASSERT(!yysmiles_lex_init(&scanner));
  try {
    size_t ltrim = setup_smiles_string(inp, scanner);
    int start_tok = static_cast<int>(START_ATOM);
    std::vector<RWMol *> molVect;
    Bond *lastBond = nullptr;
    res = yysmiles_parse(inp.c_str() + ltrim, &molVect, atom, lastBond,
                         &branchPoints, scanner, start_tok);
  } catch (...) {
    yysmiles_lex_destroy(scanner);
    throw;
  }
  yysmiles_lex_destroy(scanner);
  return res;
}
int smiles_parse(const std::string &inp, std::vector<RDKit::RWMol *> &molVect) {
  std::list<unsigned int> branchPoints;
  void *scanner;
  int res;

  TEST_ASSERT(!yysmiles_lex_init(&scanner));
  try {
    size_t ltrim = setup_smiles_string(inp, scanner);
    int start_tok = static_cast<int>(START_MOL);
    Atom *lastAtom = nullptr;
    Bond *lastBond = nullptr;
    res = yysmiles_parse(inp.c_str() + ltrim, &molVect, lastAtom, lastBond,
                         &branchPoints, scanner, start_tok);
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
    int start_tok = static_cast<int>(START_MOL);
    Atom *lastAtom = nullptr;
    Bond *lastBond = nullptr;
    res = yysmarts_parse(inp.c_str() + ltrim, &molVect, lastAtom, lastBond,
                         scanner, start_tok);
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
  std::string res;

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
            label = std::to_string(patterns.size() + 100);
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
}  // namespace

RWMol *toMol(const std::string &inp,
             int func(const std::string &, std::vector<RDKit::RWMol *> &),
             const std::string &origInp) {
  // empty strings produce empty molecules:
  if (inp.empty()) return new RWMol();
  RWMol *res = nullptr;
  std::vector<RDKit::RWMol *> molVect;
  try {
    func(inp, molVect);
    if (!molVect.empty()) {
      res = molVect[0];
      SmilesParseOps::CloseMolRings(res, false);
      SmilesParseOps::AdjustAtomChiralityFlags(res);
      // No sense leaving this bookmark intact:
      if (res->hasAtomBookmark(ci_RIGHTMOST_ATOM)) {
        res->clearAtomBookmark(ci_RIGHTMOST_ATOM);
      }
      SmilesParseOps::CleanupAfterParsing(res);
      molVect[0] = nullptr;  // NOTE: to avoid leaks on failures, this should
                             // occur last in this if.
    }
  } catch (SmilesParseException &e) {
    std::string nm = "SMILES";
    if (func == smarts_parse) {
      nm = "SMARTS";
    }
    BOOST_LOG(rdErrorLog) << nm << " Parse Error: " << e.message()
                          << " for input: '" << origInp << "'" << std::endl;
    res = nullptr;
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

Atom *toAtom(const std::string &inp, int func(const std::string &, Atom *&)) {
  // empty strings produce empty molecules:
  if (inp.empty()) return nullptr;
  Atom *res = nullptr;
  try {
    func(inp, res);
  } catch (SmilesParseException &e) {
    std::string nm = "SMILES";
    if (func != smiles_atom_parse) {
      nm = "SMARTS";
    }
    BOOST_LOG(rdErrorLog) << nm << " Parse Error: " << e.message()
                          << " for input: '" << inp << "'" << std::endl;
    res = nullptr;
  }
  return res;
}

Bond *toBond(const std::string &inp, int func(const std::string &, Bond *&)) {
  // empty strings produce empty molecules:
  if (inp.empty()) return nullptr;
  Bond *res = nullptr;
  try {
    func(inp, res);
  } catch (SmilesParseException &e) {
    std::string nm = "SMILES";
    if (func != smiles_bond_parse) {
      nm = "SMARTS";
    }
    BOOST_LOG(rdErrorLog) << nm << " Parse Error: " << e.message()
                          << " for input: '" << inp << "'" << std::endl;
    res = nullptr;
  }
  return res;
}

namespace {
void preprocessSmiles(const std::string &smiles,
                      const SmilesParserParams &params, std::string &lsmiles,
                      std::string &name, std::string &cxPart) {
  if (params.parseName && !params.allowCXSMILES) {
    std::vector<std::string> tokens;
    boost::split(tokens, smiles, boost::is_any_of(" \t"),
                 boost::token_compress_on);
    lsmiles = tokens[0];
    if (tokens.size() > 1) name = tokens[1];
  } else if (params.allowCXSMILES) {
    size_t sidx = smiles.find_first_of(" \t");
    if (sidx != std::string::npos && sidx != 0) {
      lsmiles = smiles.substr(0, sidx);
      cxPart = boost::trim_copy(smiles.substr(sidx, smiles.size() - sidx));
    }
  }

  if (lsmiles.empty()) {
    lsmiles = smiles;
  }

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
    lsmiles = smi;
  }
}
}  // namespace

Atom *SmilesToAtom(const std::string &smiles) {
  yysmiles_debug = false;

  Atom *res = nullptr;
  res = toAtom(smiles, smiles_atom_parse);
  return res;
};

Bond *SmilesToBond(const std::string &smiles) {
  yysmiles_debug = false;

  Bond *res = nullptr;
  res = toBond(smiles, smiles_bond_parse);
  return res;
};

RWMol *SmilesToMol(const std::string &smiles,
                   const SmilesParserParams &params) {
  yysmiles_debug = params.debugParse;

  std::string lsmiles, name, cxPart;
  preprocessSmiles(smiles, params, lsmiles, name, cxPart);
  // strip any leading/trailing whitespace:
  // boost::trim_if(smi,boost::is_any_of(" \t\r\n"));
  RWMol *res = nullptr;
  res = toMol(lsmiles, smiles_parse, lsmiles);

  if (res && params.allowCXSMILES && !cxPart.empty()) {
    std::string::const_iterator pos = cxPart.cbegin();
    SmilesParseOps::parseCXExtensions(*res, cxPart, pos);
    if (params.parseName && pos != cxPart.cend()) {
      std::string nmpart(pos, cxPart.cend());
      name = boost::trim_copy(nmpart);
    }
  }
  if (res && (params.sanitize || params.removeHs)) {
    try {
      if (params.removeHs) {
        bool implicitOnly = false, updateExplicitCount = true;
        MolOps::removeHs(*res, implicitOnly, updateExplicitCount,
                         params.sanitize);
      } else if (params.sanitize) {
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
  if (res && !name.empty()) res->setProp(common_properties::_Name, name);
  return res;
};

Atom *SmartsToAtom(const std::string &smiles) {
  yysmarts_debug = false;

  Atom *res = nullptr;
  res = toAtom(smiles, smarts_atom_parse);
  return res;
};

Bond *SmartsToBond(const std::string &smiles) {
  yysmarts_debug = false;

  Bond *res = nullptr;
  res = toBond(smiles, smarts_bond_parse);
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
}  // namespace RDKit
