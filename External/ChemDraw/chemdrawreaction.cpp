//
//  Copyright (c) 2024, Glysade Inc
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

#include "chemdraw.h"
#include "chemdrawreaction.h"
#include "reaction.h"
#include "utils.h"

#include <GraphMol/QueryOps.h>
#include <GraphMol/ChemReactions/SanitizeRxn.h>
#include <GraphMol/ChemReactions/ReactionUtils.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>

namespace RDKit {
using namespace RDKit::v2;
using namespace RDKit::ChemDraw;

// ChemDraw reaction API
// Convert reaction information to RDKIT reactions
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

namespace v2 {
//! Parse a text stream with ChemDraw data into a ChemicalReaction
std::vector<std::unique_ptr<ChemicalReaction>>
ChemDrawDataStreamToChemicalReactions(std::istream &inStream, bool sanitize,
                                      bool removeHs) {
  ChemDrawParserParams params;
  params.sanitize = sanitize;
  params.removeHs = removeHs;
  auto mols = MolsFromChemDrawDataStream(inStream, params);
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
    if (mols[i]->getPropIfPresent(CDX_SCHEME_ID, scheme) &&
        mols[i]->getPropIfPresent(CDX_STEP_ID, step)) {
      auto schemestep = std::pair<unsigned int, unsigned int>(scheme, step);
      schemes[schemestep].push_back(i);
    }
  }
  if (schemes.empty()) {
    return result;
  }
  for (const auto &scheme : schemes) {
    // convert atoms to queries:
    auto *res = new ChemicalReaction;
    result.push_back(std::unique_ptr<ChemicalReaction>(res));
    for (auto idx : scheme.second) {
      CHECK_INVARIANT(
          used.find(idx) == used.end(),
          "Fragment used in twice in one or more reactions, this shouldn't happen");
      if (mols[idx]->hasProp(CDX_REAGENT_ID)) {
        used.insert(idx);
        make_query_atoms(*mols[idx]);
        add_template(CDX_REAGENT_ID, reactant_templates, mols[idx]);
      } else if (mols[idx]->hasProp(CDX_AGENT_ID)) {
        used.insert(idx);
        make_query_atoms(*mols[idx]);
        add_template(CDX_AGENT_ID, agent_templates, mols[idx]);
      } else if (mols[idx]->hasProp(CDX_PRODUCT_ID)) {
        used.insert(idx);
        make_query_atoms(*mols[idx]);
        add_template(CDX_PRODUCT_ID, product_templates, mols[idx]);
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
    // ChemDraw-based reactions do not have implicit properties
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

std::vector<std::unique_ptr<ChemicalReaction>> ChemDrawToChemicalReactions(
    const std::string &rxnBlock, bool sanitize, bool removeHs) {
  std::istringstream inStream(rxnBlock);
  return ChemDrawDataStreamToChemicalReactions(inStream, sanitize, removeHs);
}

std::vector<std::unique_ptr<ChemicalReaction>> ChemDrawFileToChemicalReactions(
    const std::string &fName, bool sanitize, bool removeHs) {
  std::ifstream inStream(fName.c_str());
  std::vector<std::unique_ptr<ChemicalReaction>> res;
  ;

  if (!inStream || inStream.bad()) {
    return res;
  }
  if (!inStream.eof()) {
    return ChemDrawDataStreamToChemicalReactions(inStream, sanitize, removeHs);
  }
  return res;
}

}  // namespace v2
}  // namespace RDKit
