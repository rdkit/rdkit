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
  std::string pval;
  if (!dp_mol->getPropIfPresent(common_properties::molFileLinkNodes, pval)) {
    return;
  }
  std::vector<int> mapping;
  if (MolOps::getMolFrags(*dp_mol, mapping) > 1) {
    throw ValueErrorException(
        "LINKNODE enumeration only supported for molecules with a single "
        "fragment.");
  }

  dp_frame.reset(new RWMol(*dp_mol));
  boost::char_separator<char> pipesep("|");
  boost::char_separator<char> spacesep(" ");
  std::string attachSmarts = "";
  std::vector<std::string> linkEnums;
  std::vector<std::string> molEnums;
  d_variations.clear();
  d_pointRanges.clear();
  d_isotopeMap.clear();
  for (auto linknode : tokenizer(pval, pipesep)) {
    std::string productSmarts = "";
    tokenizer tokens(linknode, spacesep);
    std::vector<unsigned int> data;
    try {
      std::transform(tokens.begin(), tokens.end(), std::back_inserter(data),
                     [](const std::string &token) -> unsigned int {
                       return boost::lexical_cast<unsigned int>(token);
                     });
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert values in LINKNODE '" << linknode
             << "' to unsigned ints";
      throw ValueErrorException(errout.str());
    }
    // the second test here is for the atom-pairs defining the bonds
    // data[2] contains the number of bonds
    if (data.size() < 5 || data.size() < 3 + 2 * data[2]) {
      std::ostringstream errout;
      errout << "not enough values in LINKNODE '" << linknode << "'";
      throw ValueErrorException(errout.str());
    }
    unsigned int minRep = data[0];
    unsigned int maxRep = data[1];
    if (minRep == 0 || maxRep < minRep) {
      std::ostringstream errout;
      errout << "bad counts in LINKNODE '" << linknode << "'";
      throw ValueErrorException(errout.str());
    }
    d_countAtEachPoint.push_back(maxRep - minRep + 1);
    d_pointRanges.push_back(std::make_pair(minRep, maxRep));
    unsigned int nBonds = data[2];
    if (nBonds != 2) {
      UNDER_CONSTRUCTION(
          "only link nodes with 2 bonds are currently supported");
    }
    // both bonds must start from the same atom:
    if (data[3] != data[5]) {
      std::ostringstream errout;
      errout << "bonds don't start at the same atom for LINKNODE '" << linknode
             << "'";
      throw ValueErrorException(errout.str());
    }
    auto varAtom = dp_frame->getAtomWithIdx(data[3] - 1);
    auto attach1 = dp_frame->getAtomWithIdx(data[4] - 1);
    auto attach2 = dp_frame->getAtomWithIdx(data[6] - 1);
    if (!dp_frame->getBondBetweenAtoms(data[4] - 1, data[3] - 1) ||
        !dp_frame->getBondBetweenAtoms(data[4] - 1, data[3] - 1)) {
      std::ostringstream errout;
      errout << "bond not found between atoms in LINKNODE '" << linknode << "'";
      throw ValueErrorException(errout.str());
    }
    // save the isotope values:
    if (d_isotopeMap.find(1000 * data[3]) == d_isotopeMap.end()) {
      d_isotopeMap[1000 * data[3]] = varAtom->getIsotope();
    }
    if (d_isotopeMap.find(1000 * data[4]) == d_isotopeMap.end()) {
      d_isotopeMap[1000 * data[4]] = attach1->getIsotope();
    }
    if (d_isotopeMap.find(1000 * data[6]) == d_isotopeMap.end()) {
      d_isotopeMap[1000 * data[6]] = attach2->getIsotope();
    }
    varAtom->setIsotope(1000 * data[3]);
    attach1->setIsotope(1000 * data[4]);
    attach2->setIsotope(1000 * data[6]);
    d_variations.push_back(
        std::make_tuple(1000 * data[3], 1000 * data[4], 1000 * data[6]));
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
      RDLog::BlockLogs blocker;
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
