//
//  Copyright (C) 2001-2021 Greg Landrum and other RDKit contributors
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
#include <boost/lexical_cast.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/Chirality.h>

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

int smarts_parse_helper(const std::string &inp,
                        std::vector<RDKit::RWMol *> &molVect, Atom *&atom,
                        Bond *&bond, int start_tok) {
  std::list<unsigned int> branchPoints;
  void *scanner;
  int res = 1;  // initialize with fail code

  TEST_ASSERT(!yysmarts_lex_init(&scanner));
  try {
    size_t ltrim = setup_smarts_string(inp, scanner);
    unsigned numAtomsParsed = 0;
    unsigned numBondsParsed = 0;
    res = yysmarts_parse(inp.c_str() + ltrim, &molVect, atom, bond,
                         numAtomsParsed, numBondsParsed, &branchPoints, scanner,
                         start_tok);
  } catch (...) {
    yysmarts_lex_destroy(scanner);
    throw;
  }
  yysmarts_lex_destroy(scanner);

  if (res == 1) {
    std::stringstream errout;
    errout << "Failed parsing SMARTS '" << inp << "'";
    throw SmilesParseException(errout.str());
  }
  if (!branchPoints.empty()) {
    throw SmilesParseException("extra open parentheses");
  }

  return res;
}
int smarts_bond_parse(const std::string &inp, Bond *&bond) {
  auto start_tok = static_cast<int>(START_BOND);
  std::vector<RWMol *> molVect;
  Atom *atom = nullptr;
  return smarts_parse_helper(inp, molVect, atom, bond, start_tok);
}

int smarts_atom_parse(const std::string &inp, Atom *&atom) {
  auto start_tok = static_cast<int>(START_ATOM);
  std::vector<RWMol *> molVect;
  Bond *bond = nullptr;
  return smarts_parse_helper(inp, molVect, atom, bond, start_tok);
}

int smarts_parse(const std::string &inp, std::vector<RDKit::RWMol *> &molVect) {
  auto start_tok = static_cast<int>(START_MOL);
  Atom *atom = nullptr;
  Bond *bond = nullptr;
  return smarts_parse_helper(inp, molVect, atom, bond, start_tok);
}

int smiles_parse_helper(const std::string &inp,
                        std::vector<RDKit::RWMol *> &molVect, Atom *&atom,
                        Bond *&bond, int start_tok) {
  std::list<unsigned int> branchPoints;
  void *scanner;
  int res = 1;  // initialize with fail code
  unsigned numAtomsParsed = 0;
  unsigned numBondsParsed = 0;
  TEST_ASSERT(!yysmiles_lex_init(&scanner));
  try {
    size_t ltrim = setup_smiles_string(inp, scanner);
    res = yysmiles_parse(inp.c_str() + ltrim, &molVect, atom, bond,
                         numAtomsParsed, numBondsParsed, &branchPoints, scanner,
                         start_tok);
  } catch (...) {
    yysmiles_lex_destroy(scanner);
    throw;
  }
  yysmiles_lex_destroy(scanner);

  if (res == 1) {
    std::stringstream errout;
    errout << "Failed parsing SMILES '" << inp << "'";
    throw SmilesParseException(errout.str());
  }

  if (!branchPoints.empty()) {
    throw SmilesParseException("extra open parentheses");
  }
  return res;
}

int smiles_bond_parse(const std::string &inp, Bond *&bond) {
  auto start_tok = static_cast<int>(START_BOND);
  std::vector<RWMol *> molVect;
  Atom *atom = nullptr;
  return smiles_parse_helper(inp, molVect, atom, bond, start_tok);
}
int smiles_atom_parse(const std::string &inp, Atom *&atom) {
  auto start_tok = static_cast<int>(START_ATOM);
  std::vector<RWMol *> molVect;
  Bond *bond = nullptr;
  return smiles_parse_helper(inp, molVect, atom, bond, start_tok);
}

int smiles_parse(const std::string &inp, std::vector<RDKit::RWMol *> &molVect) {
  auto start_tok = static_cast<int>(START_MOL);
  Atom *atom = nullptr;
  Bond *bond = nullptr;
  return smiles_parse_helper(inp, molVect, atom, bond, start_tok);
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
      if (state.empty() || state.back() == BASE) {
        // seriously bogus input. Just return the input
        // and let the SMARTS parser itself report the error
        return sma;
      }
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
  if (inp.empty()) {
    return new RWMol();
  }
  RWMol *res = nullptr;
  std::vector<RDKit::RWMol *> molVect;
  try {
    func(inp, molVect);
    if (!molVect.empty()) {
      res = molVect[0];
      SmilesParseOps::CloseMolRings(res, false);
      SmilesParseOps::CheckChiralitySpecifications(res, true);
      SmilesParseOps::SetUnspecifiedBondTypes(res);
      SmilesParseOps::AdjustAtomChiralityFlags(res);
      // No sense leaving this bookmark intact:
      if (res->hasAtomBookmark(ci_RIGHTMOST_ATOM)) {
        res->clearAtomBookmark(ci_RIGHTMOST_ATOM);
      }
      molVect[0] = nullptr;  // NOTE: to avoid leaks on failures, this should
                             // occur last in this if.
    }
  } catch (SmilesParseException &e) {
    std::string nm = "SMILES";
    if (func == smarts_parse) {
      nm = "SMARTS";
    }
    BOOST_LOG(rdErrorLog) << nm << " Parse Error: " << e.what()
                          << " for input: '" << origInp << "'" << std::endl;
    res = nullptr;
  }
  for (auto *molPtr : molVect) {
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
  if (inp.empty()) {
    return nullptr;
  }
  Atom *res = nullptr;
  try {
    func(inp, res);
  } catch (SmilesParseException &e) {
    std::string nm = "SMILES";
    if (func != smiles_atom_parse) {
      nm = "SMARTS";
    }
    BOOST_LOG(rdErrorLog) << nm << " Parse Error: " << e.what()
                          << " for input: '" << inp << "'" << std::endl;
    res = nullptr;
  }
  return res;
}

Bond *toBond(const std::string &inp, int func(const std::string &, Bond *&)) {
  // empty strings produce empty molecules:
  if (inp.empty()) {
    return nullptr;
  }
  Bond *res = nullptr;
  try {
    func(inp, res);
  } catch (SmilesParseException &e) {
    std::string nm = "SMILES";
    if (func != smiles_bond_parse) {
      nm = "SMARTS";
    }
    BOOST_LOG(rdErrorLog) << nm << " Parse Error: " << e.what()
                          << " for input: '" << inp << "'" << std::endl;
    res = nullptr;
  }
  return res;
}

namespace {
// despite the name: works for both SMILES and SMARTS
template <typename T>
void preprocessSmiles(const std::string &smiles, const T &params,
                      std::string &lsmiles, std::string &name,
                      std::string &cxPart) {
  cxPart = "";
  name = "";
  if (params.parseName && !params.allowCXSMILES) {
    size_t sidx = smiles.find_first_of(" \t");
    if (sidx != std::string::npos && sidx != 0) {
      lsmiles = smiles.substr(0, sidx);
      name = boost::trim_copy(smiles.substr(sidx, smiles.size() - sidx));
    }
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
      for (const auto &pr : *(params.replacements)) {
        if (smi.find(pr.first) != std::string::npos) {
          loopAgain = true;
          boost::replace_all(smi, pr.first, pr.second);
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

namespace {
template <typename T>
void handleCXPartAndName(RWMol *res, const T &params, const std::string &cxPart,
                         std::string &name) {
  if (!res || cxPart.empty()) {
    return;
  }
  std::string::const_iterator pos = cxPart.cbegin();
  bool cxfailed = false;
  if (params.allowCXSMILES) {
    if (*pos == '|') {
      try {
        SmilesParseOps::parseCXExtensions(*res, cxPart, pos);
      } catch (...) {
        cxfailed = true;
        if (params.strictCXSMILES) {
          delete res;
          throw;
        }
      }
      res->setProp("_CXSMILES_Data", std::string(cxPart.cbegin(), pos));
    } else if (params.strictCXSMILES && !params.parseName &&
               pos != cxPart.cend()) {
      delete res;
      throw RDKit::SmilesParseException(
          "CXSMILES extension does not start with | and parseName=false");
    }
  }
  if (!cxfailed && params.parseName && pos != cxPart.end()) {
    std::string nmpart(pos, cxPart.cend());
    name = boost::trim_copy(nmpart);
  }
}
}  // namespace

RWMol *SmilesToMol(const std::string &smiles,
                   const SmilesParserParams &params) {
  // Calling SmilesToMol in a multithreaded context is generally safe *unless*
  // the value of debugParse is different for different threads. The if
  // statement below avoids a TSAN warning in the case where multiple threads
  // all use the same value for debugParse.
  if (yysmiles_debug != params.debugParse) {
    yysmiles_debug = params.debugParse;
  }

  std::string lsmiles, name, cxPart;
  preprocessSmiles(smiles, params, lsmiles, name, cxPart);
  // strip any leading/trailing whitespace:
  // boost::trim_if(smi,boost::is_any_of(" \t\r\n"));
  RWMol *res = nullptr;
  res = toMol(lsmiles, smiles_parse, lsmiles);
  handleCXPartAndName(res, params, cxPart, name);
  if (res && (params.sanitize || params.removeHs)) {
    if (res->hasProp(SmilesParseOps::detail::_needsDetectAtomStereo)) {
      // we encountered a wedged bond in the CXSMILES,
      // these need to be handled the same way they were in mol files
      res->clearProp(SmilesParseOps::detail::_needsDetectAtomStereo);
      if (res->getNumConformers()) {
        const auto &conf = res->getConformer();
        if (!conf.is3D()) {
          MolOps::assignChiralTypesFromBondDirs(*res, conf.getId());
        }
      }
    }
    // if we read a 3D conformer, set the stereo:
    if (res->getNumConformers() && res->getConformer().is3D()) {
      res->updatePropertyCache(false);
      MolOps::assignChiralTypesFrom3D(*res, res->getConformer().getId(), true);
    }
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
    if (res->hasProp(SmilesParseOps::detail::_needsDetectBondStereo)) {
      // we encountered either wiggly bond in the CXSMILES,
      // these need to be handled the same way they were in mol files
      res->clearProp(SmilesParseOps::detail::_needsDetectBondStereo);
      MolOps::clearSingleBondDirFlags(*res);
      MolOps::detectBondStereochemistry(*res);
    }
    // figure out stereochemistry:
    bool oVal = Chirality::getUseLegacyStereoPerception();
    if (!params.useLegacyStereo) {
      Chirality::setUseLegacyStereoPerception(false);
    }
    bool cleanIt = true, force = true, flagPossible = true;
    MolOps::assignStereochemistry(*res, cleanIt, force, flagPossible);
    Chirality::setUseLegacyStereoPerception(oVal);
  }
  if (res && res->hasProp(common_properties::_NeedsQueryScan)) {
    res->clearProp(common_properties::_NeedsQueryScan);
    if (!params.sanitize) {
      // we know that this can be the ring bond query, do ring perception if we
      // need to:
      MolOps::fastFindRings(*res);
    }
    QueryOps::completeMolQueries(res, 0xDEADBEEF);
  }

  if (res) {
    if (!params.skipCleanup) {
      SmilesParseOps::CleanupAfterParsing(res);
    }
    if (!name.empty()) {
      res->setProp(common_properties::_Name, name);
    }
  }
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

RWMol *SmartsToMol(const std::string &smarts,
                   const SmartsParserParams &params) {
  // Calling SmartsToMol in a multithreaded context is generally safe *unless*
  // the value of debugParse is different for different threads. The if
  // statement below avoids a TSAN warning in the case where multiple threads
  // all use the same value for debugParse.
  if (yysmarts_debug != params.debugParse) {
    yysmarts_debug = params.debugParse;
  }

  std::string lsmarts, name, cxPart;
  preprocessSmiles(smarts, params, lsmarts, name, cxPart);

  RWMol *res = nullptr;
  res = toMol(labelRecursivePatterns(lsmarts), smarts_parse, lsmarts);
  handleCXPartAndName(res, params, cxPart, name);
  if (res) {
    if (params.mergeHs) {
      try {
        MolOps::mergeQueryHs(*res);
      } catch (...) {
        delete res;
        throw;
      }
    }
    MolOps::setBondStereoFromDirections(*res);
    if (!params.skipCleanup) {
      SmilesParseOps::CleanupAfterParsing(res);
    }
    if (!name.empty()) {
      res->setProp(common_properties::_Name, name);
    }
  }
  return res;
};
}  // namespace RDKit
