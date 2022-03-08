//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "MolEnumerator.h"
#include "LinkNode.h"
#include <RDGeneral/Exceptions.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/format.hpp>
#include <algorithm>

typedef boost::tokenizer<boost::char_separator<char>> tokenizer;

namespace RDKit {
namespace MolEnumerator {

void LinkNodeOp::initFromMol(const ROMol &mol) {
  dp_mol.reset(new ROMol(mol));
  initFromMol();
}
void LinkNodeOp::initFromMol() {
  if (!dp_mol) {
    return;
  }
  d_variations.clear();
  d_pointRanges.clear();
  d_isotopeMap.clear();
  if (!dp_mol->hasProp(common_properties::molFileLinkNodes)) {
    return;
  }

  if (!dp_mol->hasProp(detail::idxPropName)) {
    detail::preserveOrigIndices(*dp_mol);
  }

  d_atomMap.clear();
  for (auto atom : dp_mol->atoms()) {
    unsigned int oidx;
    if (atom->getPropIfPresent(detail::idxPropName, oidx)) {
      d_atomMap[oidx] = atom;
    }
  }

  std::vector<int> mapping;
  if (MolOps::getMolFrags(*dp_mol, mapping) > 1) {
    throw ValueErrorException(
        "LINKNODE enumeration only supported for molecules with a single "
        "fragment.");
  }
  dp_frame.reset(new RWMol(*dp_mol));
  auto nodes = utils::getMolLinkNodes(*dp_frame, true, &d_atomMap);
  std::string attachSmarts = "";
  std::vector<std::string> linkEnums;
  std::vector<std::string> molEnums;
  for (auto node : nodes) {
    if (node.nBonds != 2) {
      UNDER_CONSTRUCTION(
          "only link nodes with 2 bonds are currently supported");
    }
    std::string productSmarts = "";
    d_countAtEachPoint.push_back(node.maxRep - node.minRep + 1);
    d_pointRanges.push_back(std::make_pair(node.minRep, node.maxRep));
    auto varAtom = dp_frame->getAtomWithIdx(node.bondAtoms[0].first);
    auto attach1 = dp_frame->getAtomWithIdx(node.bondAtoms[0].second);
    auto attach2 = dp_frame->getAtomWithIdx(node.bondAtoms[1].second);
    // save the isotope values:
    if (d_isotopeMap.find(1000 * (node.bondAtoms[0].first + 1)) ==
        d_isotopeMap.end()) {
      d_isotopeMap[1000 * (node.bondAtoms[0].first + 1)] =
          varAtom->getIsotope();
    }
    if (d_isotopeMap.find(1000 * (node.bondAtoms[0].second + 1)) ==
        d_isotopeMap.end()) {
      d_isotopeMap[1000 * (node.bondAtoms[0].second + 1)] =
          attach1->getIsotope();
    }
    if (d_isotopeMap.find(1000 * (node.bondAtoms[1].second + 1)) ==
        d_isotopeMap.end()) {
      d_isotopeMap[1000 * (node.bondAtoms[1].second + 1)] =
          attach2->getIsotope();
    }
    varAtom->setIsotope(1000 * (node.bondAtoms[0].first + 1));
    attach1->setIsotope(1000 * (node.bondAtoms[0].second + 1));
    attach2->setIsotope(1000 * (node.bondAtoms[1].second + 1));
    d_variations.push_back(
        std::make_tuple(1000 * (node.bondAtoms[0].first + 1),
                        1000 * (node.bondAtoms[0].second + 1),
                        1000 * (node.bondAtoms[1].second + 1)));
  }
}

std::vector<size_t> LinkNodeOp::getVariationCounts() const {
  return d_countAtEachPoint;
}

std::unique_ptr<ROMol> LinkNodeOp::operator()(
    const std::vector<size_t> &which) const {
  PRECONDITION(dp_mol, "no molecule");
  PRECONDITION(dp_frame, "not initialized");
  if (which.size() != d_countAtEachPoint.size()) {
    throw ValueErrorException("bad element choice in enumeration");
  }
  // quick error checking before we do any work:
  for (size_t i = 0; i < which.size(); ++i) {
    if (which[i] >= d_countAtEachPoint[i]) {
      throw ValueErrorException("bad element value in enumeration");
    }
  }
  // we do the enumeration of each of the variation points independantly
  ROMOL_SPTR res(new ROMol(*dp_frame));
  for (size_t i = 0; i < which.size(); ++i) {
    auto variationIdx = i + 1;
    auto variationCount = d_pointRanges[i].first + which[i];
    auto reactFormat =
        boost::format("[%d*:%d]-[%d*:%d]-[%d*:%d]") %
        (std::get<1>(d_variations[i])) % (variationIdx * 100 + 1) %
        (std::get<0>(d_variations[i])) % (variationIdx * 100) %
        (std::get<2>(d_variations[i])) % (variationIdx * 100 + 2);
    auto reacts = reactFormat.str();

    auto prodFormat = boost::format("[*:%d]") % (variationIdx * 100 + 1);
    auto prods = prodFormat.str();
    for (size_t j = 0; j < variationCount - 1; ++j) {
      prods += (boost::format("-[*:%d]") % (variationIdx * 100)).str();
    }
    prods +=
        (boost::format("-[%d*:%d]-[*:%d]") % (std::get<0>(d_variations[i])) %
         (variationIdx * 100) % (variationIdx * 100 + 2))
            .str();
    auto reactSmarts = reacts + ">>" + prods;
    // std::cerr << "variation index: " << variationIdx
    //           << " count: " << variationCount << " reaction " << reactSmarts
    //           << std::endl;
    std::unique_ptr<ChemicalReaction> rxn(
        RxnSmartsToChemicalReaction(reactSmarts));
    ASSERT_INVARIANT(rxn, "reaction could not be constructed");
    // we expect warnings for these alchemical reactions. :-)
    bool silent = true;
    rxn->initReactantMatchers(silent);
    ROMOL_SPTR reactant(new ROMol(*res));
    std::vector<MOL_SPTR_VECT> ps;
    {
      RDLog::LogStateSetter blocker;
      ps = rxn->runReactant(reactant, 0);
    }
    ASSERT_INVARIANT(!ps.empty(), "no products from reaction");
    ASSERT_INVARIANT(ps[0].size() == 1, "too many products from reaction");
    res = ps[0][0];
    // std::cerr << "   APPLICATION: " << MolToSmiles(*res) << std::endl;
  }
  // reset the original isotopes
  for (auto atom : res->atoms()) {
    if (atom->getIsotope() >= 1000) {
      auto iter = d_isotopeMap.find(atom->getIsotope());
      if (iter != d_isotopeMap.end()) {
        atom->setIsotope(iter->second);
      }
    }
  }
  return std::unique_ptr<ROMol>(new ROMol(*res));
}

}  // namespace MolEnumerator

}  // namespace RDKit
