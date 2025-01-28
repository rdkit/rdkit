//
//  Copyright (c) 2007-2021, Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
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
//       products derived from this software without specific prior written
//       permission.
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
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesParseOps.h>
#include <boost/range/iterator_range.hpp>

#include <boost/algorithm/string.hpp>
#include <vector>
#include <string>
#include "ReactionUtils.h"

namespace RDKit {
namespace v2 {
namespace ReactionParser {

namespace DaylightParserUtils {
std::vector<std::string> splitSmartsIntoComponents(
    const std::string &reactText) {
  std::vector<std::string> res;
  unsigned int pos = 0;
  unsigned int blockStart = 0;
  unsigned int level = 0;
  unsigned int inBlock = 0;
  while (pos < reactText.size()) {
    if (reactText[pos] == '(') {
      if (pos == blockStart) {
        inBlock = 1;
      }
      ++level;
    } else if (reactText[pos] == ')') {
      if (level == 1 && inBlock) {
        // this closes a block
        inBlock = 2;
      }
      --level;
    } else if (level == 0 && reactText[pos] == '.') {
      if (inBlock == 2) {
        std::string element =
            reactText.substr(blockStart + 1, pos - blockStart - 2);
        res.push_back(element);
      } else {
        std::string element = reactText.substr(blockStart, pos - blockStart);
        res.push_back(element);
      }
      blockStart = pos + 1;
      inBlock = 0;
    }
    ++pos;
  }
  if (blockStart < pos) {
    if (inBlock == 2) {
      std::string element =
          reactText.substr(blockStart + 1, pos - blockStart - 2);
      res.push_back(element);
    } else {
      std::string element = reactText.substr(blockStart, pos - blockStart);
      res.push_back(element);
    }
  }
  return res;
}

std::unique_ptr<RWMol> constructMolFromString(
    const std::string &txt, const ReactionSmartsParserParams &params,
    bool useSmiles) {
  if (!useSmiles) {
    SmilesParse::SmartsParserParams ps;
    ps.replacements = params.replacements;
    ps.allowCXSMILES = false;
    ps.parseName = false;
    ps.mergeHs = false;
    ps.skipCleanup = true;
    return SmilesParse::MolFromSmarts(txt, ps);
  } else {
    SmilesParse::SmilesParserParams ps;
    ps.replacements = params.replacements;
    ps.allowCXSMILES = false;
    ps.parseName = false;
    ps.sanitize = params.sanitize;
    ps.removeHs = false;
    ps.skipCleanup = true;
    return SmilesParse::MolFromSmiles(txt, ps);
  }
}

}  // end of namespace DaylightParserUtils

namespace {
void removeSpacesAround(std::string &text, size_t pos) {
  auto nextp = pos + 1;
  while (nextp < text.size() && (text[nextp] == ' ' || text[nextp] == '\t')) {
    text.erase(nextp, 1);
  }
  if (pos > 0) {
    nextp = pos - 1;
    while (text[nextp] == ' ' || text[nextp] == '\t') {
      text.erase(nextp, 1);
      if (nextp > 0) {
        --nextp;
      } else {
        break;
      }
    }
  }
}

std::unique_ptr<ChemicalReaction> parseReaction(
    const std::string &origText, const ReactionSmartsParserParams &params,
    bool useSmiles) {
  std::string text = origText;
  std::string cxPart;
  if (params.allowCXSMILES) {
    auto sidx = origText.find_first_of("|");
    if (sidx != std::string::npos && sidx != 0) {
      text = origText.substr(0, sidx);
      cxPart = boost::trim_copy(origText.substr(sidx, origText.size() - sidx));
    }
  }
  // remove any spaces at the beginning, end, or before the '>'s
  boost::trim(text);
  std::vector<std::size_t> pos;
  for (std::size_t i = 0; i < text.length(); ++i) {
    if (text[i] == '>' && (i == 0 || text[i - 1] != '-')) {
      pos.push_back(i);
    }
  }
  if (pos.size() < 2) {
    throw ChemicalReactionParserException(
        "a reaction requires at least two > characters");
  }

  // remove spaces around ">" symbols
  for (auto p : boost::make_iterator_range(pos.rbegin(), pos.rend())) {
    removeSpacesAround(text, p);
  }

  // remove spaces around "." symbols
  pos.clear();
  for (std::size_t i = 0; i < text.length(); ++i) {
    if (text[i] == '.') {
      pos.push_back(i);
    }
  }
  for (auto p : boost::make_iterator_range(pos.rbegin(), pos.rend())) {
    removeSpacesAround(text, p);
  }

  // we shouldn't have whitespace left in the reaction string, so go ahead and
  // split and strip:
  auto sidx = text.find_first_of(" \t");
  if (sidx != std::string::npos && sidx != 0) {
    text = text.substr(0, sidx);
  }

  // re-find the '>' characters so that we can split on them
  pos.clear();
  for (std::size_t i = 0; i < text.length(); ++i) {
    if (text[i] == '>' && (i == 0 || text[i - 1] != '-')) {
      pos.push_back(i);
    }
  }

  // there's always the chance that one or more of the ">" was in the name
  // part, so verify that we have exactly two:
  if (pos.size() < 2) {
    throw ChemicalReactionParserException(
        "a reaction requires at least two > characters");
  }
  if (pos.size() > 2) {
    throw ChemicalReactionParserException("multi-step reactions not supported");
  }

  auto pos1 = pos[0];
  auto pos2 = pos[1];

  auto reactText = text.substr(0, pos1);
  std::string agentText;
  if (pos2 != pos1 + 1) {
    agentText = text.substr(pos1 + 1, (pos2 - pos1) - 1);
  }
  auto productText = text.substr(pos2 + 1);

  // recognize changes within the same molecules, e.g., intra molecular bond
  // formation therefore we need to correctly interpret parenthesis and dots
  // in the reaction smarts
  auto reactSmarts = DaylightParserUtils::splitSmartsIntoComponents(reactText);
  auto productSmarts =
      DaylightParserUtils::splitSmartsIntoComponents(productText);

  auto rxn = std::make_unique<ChemicalReaction>();

  for (const auto &txt : reactSmarts) {
    auto mol =
        DaylightParserUtils::constructMolFromString(txt, params, useSmiles);
    if (!mol) {
      std::string errMsg = "Problems constructing reactant from SMARTS: ";
      errMsg += txt;
      throw ChemicalReactionParserException(errMsg);
    }
    rxn->addReactantTemplate(ROMOL_SPTR(mol.release()));
  }

  for (const auto &txt : productSmarts) {
    auto mol =
        DaylightParserUtils::constructMolFromString(txt, params, useSmiles);
    if (!mol) {
      std::string errMsg = "Problems constructing product from SMARTS: ";
      errMsg += txt;
      throw ChemicalReactionParserException(errMsg);
    }
    rxn->addProductTemplate(ROMOL_SPTR(mol.release()));
  }
  updateProductsStereochem(rxn.get());

  // allow a reaction template to have no agent specified
  if (agentText.size() != 0) {
    auto agentMol = DaylightParserUtils::constructMolFromString(
        agentText, params, useSmiles);
    if (!agentMol) {
      std::string errMsg = "Problems constructing agent from SMARTS: ";
      errMsg += agentText;
      throw ChemicalReactionParserException(errMsg);
    }
    std::vector<ROMOL_SPTR> agents = MolOps::getMolFrags(*agentMol, false);
    for (auto &agent : agents) {
      rxn->addAgentTemplate(agent);
    }
  }

  if (params.allowCXSMILES && !cxPart.empty()) {
    unsigned int startAtomIdx = 0;
    unsigned int startBondIdx = 0;
    for (auto &mol : boost::make_iterator_range(rxn->beginReactantTemplates(),
                                                rxn->endReactantTemplates())) {
      SmilesParseOps::parseCXExtensions(*static_cast<RWMol *>(mol.get()),
                                        cxPart, startAtomIdx, startBondIdx);
      startAtomIdx += mol->getNumAtoms();
      startBondIdx += mol->getNumBonds();
    }
    for (auto &mol : boost::make_iterator_range(rxn->beginAgentTemplates(),
                                                rxn->endAgentTemplates())) {
      SmilesParseOps::parseCXExtensions(*static_cast<RWMol *>(mol.get()),
                                        cxPart, startAtomIdx, startBondIdx);
      startAtomIdx += mol->getNumAtoms();
      startBondIdx += mol->getNumBonds();
    }
    for (auto &mol : boost::make_iterator_range(rxn->beginProductTemplates(),
                                                rxn->endProductTemplates())) {
      SmilesParseOps::parseCXExtensions(*static_cast<RWMol *>(mol.get()),
                                        cxPart, startAtomIdx, startBondIdx);
      startAtomIdx += mol->getNumAtoms();
      startBondIdx += mol->getNumBonds();
    }
  }

  // final cleanups:
  for (auto &mol : boost::make_iterator_range(rxn->beginReactantTemplates(),
                                              rxn->endReactantTemplates())) {
    SmilesParseOps::CleanupAfterParsing(static_cast<RWMol *>(mol.get()));
  }
  for (auto &mol : boost::make_iterator_range(rxn->beginAgentTemplates(),
                                              rxn->endAgentTemplates())) {
    SmilesParseOps::CleanupAfterParsing(static_cast<RWMol *>(mol.get()));
  }
  for (auto &mol : boost::make_iterator_range(rxn->beginProductTemplates(),
                                              rxn->endProductTemplates())) {
    SmilesParseOps::CleanupAfterParsing(static_cast<RWMol *>(mol.get()));
  }

  // "SMARTS"-based reactions have implicit properties
  rxn->setImplicitPropertiesFlag(true);

  return rxn;
}
}  // namespace

std::unique_ptr<ChemicalReaction> ReactionFromSmarts(
    const std::string &origText, const ReactionSmartsParserParams &options) {
  return parseReaction(origText, options, false);
}
std::unique_ptr<ChemicalReaction> ReactionFromSmiles(
    const std::string &origText, const ReactionSmartsParserParams &options) {
  return parseReaction(origText, options, true);
}
}  // namespace ReactionParser
}  // namespace v2
}  // namespace RDKit
