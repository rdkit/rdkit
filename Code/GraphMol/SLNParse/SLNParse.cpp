//
//  Copyright (c) 2008, Novartis Institutes for BioMedical Research Inc.
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
//       products derived from this software without specific prior
//       written permission.
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
//  Created by Greg Landrum, Sept. 2006
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SLNParse/SLNParse.h>
#include <GraphMol/SLNParse/SLNAttribs.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/Invariant.h>
#include <boost/algorithm/string.hpp>
#include <regex>

int yysln_parse(const char *, std::vector<RDKit::RWMol *> *, bool, void *);
int yysln_lex_init(void **);
void yysln_set_extra(void *, void *);
int yysln_lex_destroy(void *);
void setup_sln_string(const std::string &text, void *);
extern int yysln_debug;

int sln_parse(const std::string &inp, bool doQueries,
              std::vector<RDKit::RWMol *> &molVect) {
  void *scanner;
  TEST_ASSERT(!yysln_lex_init(&scanner));
  setup_sln_string(inp, scanner);
  yysln_set_extra((void *)doQueries, scanner);
  int res = 1;  // initialize with failed state
  try {
    res = yysln_parse(inp.c_str(), &molVect, doQueries, scanner);
  } catch (...) {
    yysln_lex_destroy(scanner);
    throw;
  }
  yysln_lex_destroy(scanner);

  if (res == 1) {
    std::stringstream errout;
    errout << "Failed parsing SLN '" << inp << "'";
    throw RDKit::SLNParseException(errout.str());
  }

  return res;
}

namespace RDKit {
namespace SLNParse {
std::vector<RDKit::RWMol *> molList_g;

void finalizeQueryMol(ROMol *mol, bool mergeHs) {
  PRECONDITION(mol, "bad query molecule");

  // do we need to remove the Hs from the molecule?
  if (mergeHs) {
    for (ROMol::AtomIterator atomIt = mol->beginAtoms();
         atomIt != mol->endAtoms(); ++atomIt) {
      // set a query for the H count:
      if ((*atomIt)->getNumExplicitHs()) {
        (*atomIt)->expandQuery(
            makeAtomHCountQuery((*atomIt)->getNumExplicitHs()));
      }
    }
  }

  // we don't want to sanitize, but we do need to get
  // some ring info:
  VECT_INT_VECT sssr;
  MolOps::symmetrizeSSSR(*mol, sssr);
  int rootIdx = -1;
  for (ROMol::AtomIterator atomIt = mol->beginAtoms();
       atomIt != mol->endAtoms(); ++atomIt) {
    SLNParse::parseFinalAtomAttribs(*atomIt, true);
    if ((*atomIt)->hasProp(common_properties::_starred)) {
      if (rootIdx > -1) {
        BOOST_LOG(rdErrorLog) << "SLN Error: multiple starred atoms in a "
                                 "recursive query. Extra stars ignored"
                              << std::endl;
      } else {
        rootIdx = (*atomIt)->getIdx();
      }
    }
  }
  if (rootIdx > -1) {
    mol->setProp(common_properties::_queryRootAtom, rootIdx);
  }
}

std::string replaceSLNMacroAtoms(std::string inp) {
  const std::regex defn("\\{(.+?):(.+?)\\}");
  const char *empty = "";

  std::string res;
  // remove any macro definitions:
  res = std::regex_replace(
      inp, defn, empty,
      std::regex_constants::match_default | std::regex_constants::format_sed);

  if (res != inp) {
    // there are macro definitions, we're going to replace
    // the macro atoms in the input:
    std::string::const_iterator start, end;
    start = inp.begin();
    end = inp.end();
    std::match_results<std::string::const_iterator> what;
    std::regex_constants::match_flag_type flags =
        std::regex_constants::match_default;
    while (std::regex_search(start, end, what, defn, flags)) {
      std::string macroNm(what[1].first, what[1].second);
      std::string macroVal(what[2].first, what[2].second);
      res = std::regex_replace(res, std::regex(macroNm), macroVal.c_str(),
                               std::regex_constants::match_default |
                                   std::regex_constants::format_sed);
      // update search position:
      start = what[0].second;
      // update flags:
      flags |= std::regex_constants::match_prev_avail;
    }
  }
  return res;
}

RWMol *toMol(std::string inp, bool doQueries, int debugParse) {
  boost::trim_if(inp, boost::is_any_of(" \t\r\n"));
  inp = replaceSLNMacroAtoms(inp);
  if (debugParse) {
    std::cerr << "****** PARSING SLN: ->" << inp << "<-" << std::endl;
  }

  RWMol *res = nullptr;
  std::vector<RDKit::RWMol *> molVect;
  try {
    sln_parse(inp, doQueries, molVect);
    if (molVect.size() > 0) {
      res = molVect[0];

      for (ROMol::BOND_BOOKMARK_MAP::const_iterator bmIt =
               res->getBondBookmarks()->begin();
           bmIt != res->getBondBookmarks()->end(); ++bmIt) {
        if (bmIt->first > 0 &&
            bmIt->first < static_cast<int>(res->getNumAtoms())) {
          std::stringstream err;
          err << "SLN Parser error: Ring closure " << bmIt->first
              << " does not have a corresponding opener.";
          throw SLNParseException(err.str());
        }
      }

      molVect[0] = nullptr;
    }

  } catch (SLNParseException &e) {
    BOOST_LOG(rdErrorLog) << e.what() << std::endl;
    res = nullptr;
  }
  if (res) {
    // cleanup:
    res->clearAllAtomBookmarks();
    res->clearAllBondBookmarks();

    // set up the chirality flags as soon as the molecule is finished
    // since we'll be removing Hs later and that will break things:
    adjustAtomChiralities(res);
  }

  for (auto &iter : molVect) {
    if (iter) {
      CleanupAfterParse(iter);
      delete iter;
    }
  }

  return res;
};
}  // namespace SLNParse

RWMol *SLNToMol(const std::string &sln, bool sanitize, int debugParse) {
  // FIX: figure out how to reset lexer state
  yysln_debug = debugParse;
  // strip any leading/trailing whitespace:
  // boost::trim_if(sln,boost::is_any_of(" \t\r\n"));

  RWMol *res = SLNParse::toMol(sln, false, debugParse);
  if (res) {
    for (ROMol::AtomIterator atomIt = res->beginAtoms();
         atomIt != res->endAtoms(); ++atomIt) {
      SLNParse::parseFinalAtomAttribs(*atomIt, false);
    }
    if (sanitize) {
      // we're going to remove explicit Hs from the graph,
      // this triggers a sanitization, so we do not need to
      // worry about doing one here:
      try {
        MolOps::removeHs(*res);
      } catch (...) {
        delete res;
        throw;
      }
    }
  }
  return res;
};

RWMol *SLNQueryToMol(const std::string &sln, bool mergeHs, int debugParse) {
  yysln_debug = debugParse;
  // strip any leading/trailing whitespace:
  // boost::trim_if(sln,boost::is_any_of(" \t\r\n"));
  RWMol *res = SLNParse::toMol(sln, true, debugParse);
  if (res) {
    SLNParse::finalizeQueryMol(res, mergeHs);
  }
  return res;
};

void SLNParse::CleanupAfterParse(RWMol *mol) {
  PRECONDITION(mol, "no molecule");
  // blow out any partial bonds:
  RWMol::BOND_BOOKMARK_MAP *marks = mol->getBondBookmarks();
  auto markI = marks->begin();
  while (markI != marks->end()) {
    RWMol::BOND_PTR_LIST &bonds = markI->second;
    for (auto &bond : bonds) {
      delete bond;
    }
    ++markI;
  }
}

}  // namespace RDKit
