//
//  Copyright (c) 2022 Brian P Kelley
//  All rights reserved.
//
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/QueryOps.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/FileParserUtils.h>
#include <sstream>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/FileParseException.h>
#include "SanitizeRxn.h"

#include <fstream>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/tokenizer.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include "ReactionUtils.h"

namespace RDKit {

namespace {
void make_query_atoms(RWMol &mol) {
  for (auto &atom : mol.atoms()) {
    QueryOps::replaceAtomWithQueryAtom(&mol, atom);
  }
}

void add_template(const std::string &prop, std::map<int, ROMOL_SPTR> &templates,
                  std::unique_ptr<RWMol> &mol) {
  auto reactant_idx = mol->getProp<int>(prop);
  if (templates.find(reactant_idx) != templates.end()) {
    templates[reactant_idx] =
        ROMOL_SPTR(combineMols(*templates[reactant_idx], *mol));
  } else {
    templates[reactant_idx] = ROMOL_SPTR(std::move(mol));
  }
}
}  // namespace

//! Parse a text stream with CDXML data into a ChemicalReaction
std::vector<std::unique_ptr<ChemicalReaction>>
CDXMLDataStreamToChemicalReactions(std::istream &inStream, bool sanitize,
                                   bool removeHs) {
  auto mols = CDXMLDataStreamToMols(inStream, sanitize, removeHs);
  std::vector<std::unique_ptr<ChemicalReaction>> result;

  std::map<std::pair<unsigned int, unsigned int>, std::vector<unsigned int>>
      schemes;
  std::set<unsigned int> used;
  std::map<int, ROMOL_SPTR> reactant_templates;
  std::map<int, ROMOL_SPTR> product_templates;
  std::map<int, ROMOL_SPTR> agent_templates;

  for (size_t i = 0; i < mols.size(); ++i) {
    unsigned int step = 0;
    unsigned int scheme = 0;
    if (mols[i]->getPropIfPresent("CDX_SCHEME_ID", scheme) &&
        mols[i]->getPropIfPresent("CDX_STEP_ID", step)) {
      auto schemestep = std::pair<unsigned int, unsigned int>(scheme, step);
      schemes[schemestep].push_back(i);
    }
  }
  if (schemes.empty()) {
    return result;
  }
  for (const auto &scheme : schemes) {
    // convert atoms to queries:
    ChemicalReaction *res = new ChemicalReaction;
    result.push_back(std::unique_ptr<ChemicalReaction>(res));
    for (auto idx : scheme.second) {
      CHECK_INVARIANT(
          used.find(idx) == used.end(),
          "Fragment used in twice in one or more reactions, this shouldn't happen");
      if (mols[idx]->hasProp("CDX_REAGENT_ID")) {
        used.insert(idx);
        make_query_atoms(*mols[idx]);
        add_template("CDX_REAGENT_ID", reactant_templates, mols[idx]);
      } else if (mols[idx]->hasProp("CDX_AGENT_ID")) {
        used.insert(idx);
        make_query_atoms(*mols[idx]);
        add_template("CDX_AGENT_ID", agent_templates, mols[idx]);
      } else if (mols[idx]->hasProp("CDX_PRODUCT_ID")) {
        used.insert(idx);
        make_query_atoms(*mols[idx]);
        add_template("CDX_PRODUCT_ID", product_templates, mols[idx]);
      }
    }
    for (auto reactant : reactant_templates) {
      res->addReactantTemplate(reactant.second);
    }
    for (auto reactant : agent_templates) {
      res->addAgentTemplate(reactant.second);
    }
    for (auto reactant : product_templates) {
      res->addProductTemplate(reactant.second);
    }
    updateProductsStereochem(res);
    // CDXML-based reactions do not have implicit properties
    res->setImplicitPropertiesFlag(false);

    if (!sanitize) {  // we still need to fix the reaction for smarts style
                      // matching
      unsigned int failed;
      RxnOps::sanitizeRxn(
          *res, failed,
          RxnOps::SANITIZE_ADJUST_REACTANTS | RxnOps::SANITIZE_ADJUST_PRODUCTS,
          RxnOps::MatchOnlyAtRgroupsAdjustParams());
    }
  }
  return result;
}

std::vector<std::unique_ptr<ChemicalReaction>> CDXMLToChemicalReactions(
    const std::string &rxnBlock, bool sanitize, bool removeHs) {
  std::istringstream inStream(rxnBlock);
  return CDXMLDataStreamToChemicalReactions(inStream, sanitize, removeHs);
}

std::vector<std::unique_ptr<ChemicalReaction>> CDXMLFileToChemicalReactions(
    const std::string &fName, bool sanitize, bool removeHs) {
  std::ifstream inStream(fName.c_str());
  std::vector<std::unique_ptr<ChemicalReaction>> res;
  ;

  if (!inStream || inStream.bad()) {
    return res;
  }
  if (!inStream.eof()) {
    return CDXMLDataStreamToChemicalReactions(inStream, sanitize, removeHs);
  }
  return res;
}
}  // namespace RDKit
