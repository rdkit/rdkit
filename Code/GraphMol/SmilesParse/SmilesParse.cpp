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

class AbstractParser {
 public:
  virtual ~AbstractParser() {
    for (auto molPtr : d_molVect) {
      if (molPtr) {
        // Clean-up the bond bookmarks when not calling CloseMolRings
        SmilesParseOps::CleanupAfterParseError(molPtr);
        delete molPtr;
      }
    }
  }
  virtual int yyParse(const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, std::list<unsigned int> *branchPoints, int& start_token) = 0;
  virtual size_t setupString(const std::string &text) = 0;

  RWMol *toMol(const std::string &inp, const std::string &origInp) {
    RWMol *res = nullptr;
    // empty strings produce empty molecules:
    if (inp.empty()) {
      return new RWMol();
    }
    auto start_tok = static_cast<int>(START_MOL);
    Atom *atom = nullptr;
    Bond *bond = nullptr;
    auto errorMsg = parse(inp, atom, bond, start_tok);
    if (errorMsg.empty() && !d_molVect.empty() && d_molVect[0]) {
      res = d_molVect[0];
      errorMsg = SmilesParseOps::CloseMolRings(*res, false);
      if (errorMsg.empty()) {
        errorMsg = SmilesParseOps::CheckChiralitySpecifications(*res, true);
      }
      if (errorMsg.empty()) {
        SmilesParseOps::SetUnspecifiedBondTypes(*res);
        SmilesParseOps::AdjustAtomChiralityFlags(*res);
        // No sense leaving this bookmark intact:
        if (res->hasAtomBookmark(ci_RIGHTMOST_ATOM)) {
          res->clearAtomBookmark(ci_RIGHTMOST_ATOM);
        }
      }
    }
    if (errorMsg.empty()) {
      d_molVect[0] = nullptr;  // NOTE: to avoid leaks on failures, this should
                               // occur last in this if.
    } else {
      logErrorMsg(errorMsg, origInp);
      res = nullptr;
    }
    return res;
  }

  Atom *toAtom(const std::string &inp) {
    // empty strings produce nullptr:
    if (inp.empty()) {
      return nullptr;
    }
    Atom *res = nullptr;
    auto start_tok = static_cast<int>(START_ATOM);
    Bond *bond = nullptr;
    auto errorMsg = parse(inp, res, bond, start_tok);
    if (!errorMsg.empty()) {
      logErrorMsg(errorMsg, inp);
      res = nullptr;
    }
    return res;
  }

  Bond *toBond(const std::string &inp) {
    // empty strings produce nullptr:
    if (inp.empty()) {
      return nullptr;
    }
    Bond *res = nullptr;
    auto start_tok = static_cast<int>(START_BOND);
    Atom *atom = nullptr;
    auto errorMsg = parse(inp, atom, res, start_tok);
    if (!errorMsg.empty()) {
      logErrorMsg(errorMsg, inp);
      res = nullptr;
    }
    return res;
  }

 protected:
  void *d_scanner;
  std::vector<RDKit::RWMol *> d_molVect;

 private:
  virtual const char *parserName() const = 0;
  void logErrorMsg(const std::string &errorMsg, const std::string &origInp) {
    BOOST_LOG(rdErrorLog) << parserName() << " Parse Error: " << errorMsg
                          << " for input: '" << origInp << "'" << std::endl;
  }
  std::string parse(const std::string &inp, Atom *&atom,
                        Bond *&bond, int start_tok) {
    std::list<unsigned int> branchPoints;
    std::stringstream errorMsg;
    int res = 1;  // initialize with fail code

    try {
      size_t ltrim = setupString(inp);
      unsigned int numAtomsParsed = 0;
      unsigned int numBondsParsed = 0;
      res = yyParse(inp.c_str() + ltrim, &d_molVect, atom, bond, numAtomsParsed, numBondsParsed, &branchPoints, start_tok);
    } catch (...) {
    }

    if (res == 1) {
      errorMsg << "Failed parsing " << parserName() << " '" << inp << "'";
    }

    if (!branchPoints.empty()) {
      errorMsg << (res == 1 ? "\n" : "") << "extra open parentheses";
    }
    return errorMsg.str();
  }
};

class SmilesParser : public AbstractParser {
 public:
  SmilesParser() {
    TEST_ASSERT(!yysmiles_lex_init(&d_scanner));
  }
  ~SmilesParser() {
    if (d_scanner) {
      yysmiles_lex_destroy(d_scanner);
    }
  }
 private:
  inline int yyParse(const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, std::list<unsigned int> *branchPoints, int& start_token) {
    return yysmiles_parse(input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, d_scanner, start_token);
  }
  inline size_t setupString(const std::string &text) {
    return setup_smiles_string(text, d_scanner);
  }
  const char *parserName() const {
    static const char *name = "SMILES";
    return name;
  };
};

class SmartsParser : public AbstractParser {
 public:
  SmartsParser() {
    TEST_ASSERT(!yysmarts_lex_init(&d_scanner));
  }
  ~SmartsParser() {
    if (d_scanner) {
      yysmarts_lex_destroy(d_scanner);
    }
  }
 private:
  inline int yyParse(const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, std::list<unsigned int> *branchPoints, int& start_token) {
    return yysmarts_parse(input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, d_scanner, start_token);
  }
  inline size_t setupString(const std::string &text) {
    return setup_smarts_string(text, d_scanner);
  }
  const char *parserName() const {
    static const char *name = "SMARTS";
    return name;
  };
};

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

  SmilesParser smilesParser;
  return smilesParser.toAtom(smiles);
};

Bond *SmilesToBond(const std::string &smiles) {
  yysmiles_debug = false;

  SmilesParser smilesParser;
  return smilesParser.toBond(smiles);
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
  SmilesParser smilesParser;
  res = smilesParser.toMol(lsmiles, lsmiles);
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
    bool cleanIt = true, force = true, flagPossible = true;
    MolOps::assignStereochemistry(*res, cleanIt, force, flagPossible);
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

Atom *SmartsToAtom(const std::string &smarts) {
  yysmarts_debug = false;

  SmartsParser smartsParser;
  return smartsParser.toAtom(smarts);
};

Bond *SmartsToBond(const std::string &smarts) {
  yysmarts_debug = false;

  SmartsParser smartsParser;
  return smartsParser.toBond(smarts);
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
  SmartsParser smartsParser;
  res = smartsParser.toMol(labelRecursivePatterns(lsmarts), lsmarts);
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
