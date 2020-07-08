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
    unsigned int variationIdx = d_countAtEachPoint.size();
    unsigned int nBonds = data[2];
    if (nBonds != 2) {
      UNDER_CONSTRUCTION(
          "only link nodes with 2 bonds are currently supported");
    }
    std::vector<std::string> bondSmarts;
    for (unsigned int idx = 0; idx < nBonds; ++idx) {
      auto begAtIdx = data[3 + 2 * idx] - 1;
      auto endAtIdx = data[3 + 2 * idx + 1] - 1;
      auto bnd = dp_frame->getBondBetweenAtoms(begAtIdx, endAtIdx);
      if (!bnd) {
        std::ostringstream errout;
        errout << "bond not found between atoms for LINKNODE '" << linknode
               << "'";
        throw ValueErrorException(errout.str());
      }
      bondSmarts.push_back(SmartsWrite::GetBondSmarts((QueryBond *)(bnd)));
      auto bondLabel = variationIdx * 100 + 10 * idx;
      Atom dummy(0);
      dummy.setIsotope(bondLabel);
      auto d1idx = dp_frame->addAtom(&dummy);
      dp_frame->addBond(begAtIdx, d1idx, Bond::SINGLE);
      dummy.setIsotope(bondLabel + 1);
      auto d2idx = dp_frame->addAtom(&dummy);
      dp_frame->addBond(endAtIdx, d2idx, Bond::SINGLE);
      dp_frame->removeBond(begAtIdx, endAtIdx);
    }
    // both bonds must start from the same atom:
    if (data[3] != data[5]) {
      std::ostringstream errout;
      errout << "bonds don't start at the same atom for LINKNODE '" << linknode
             << "'";
      throw ValueErrorException(errout.str());
    }

#if 0    
    auto fmt1 = boost::format("[%d#0]-[*:%d]-[%d#0]") % (variationIdx * 100) %
                (variationIdx) % (variationIdx * 100 + 10);
    std::string reactSmarts = fmt1.str();
    // Unfortunately the enumeration reactions for the first pass, where the
    // two attachment points are on the atom, and the second pass, where the two
    // attachment points are on different atoms, have to be different:
    auto fmt1 =
        boost::format("[%d#0]-[*:%d]-[%d#0]>>[%d#0]-[*:%d]-[*:%d]-[%d#0]") %
        (variationIdx * 100) % (variationIdx) % (variationIdx * 100 + 10) %
        (variationIdx * 100) % (variationIdx) % (variationIdx) %
        (variationIdx * 100 + 10);
    auto fmt2 =
        boost::format(
            "([%d#0]-[*:%d].[*:%d]-[%d#0])>>[%d#0]-[*:%d]-[*:%d]-[%d#0]") %
        (variationIdx * 100) % (variationIdx) % (variationIdx * 100 + 10) %
        (variationIdx * 100) % (variationIdx) % (variationIdx) %
        (variationIdx * 100 + 10);
    auto rxn1 = std::shared_ptr<ChemicalReaction>(
        RxnSmartsToChemicalReaction(fmt1.str()));
    ASSERT_INVARIANT(rxn1, "reaction is somehow bad");
    rxn1->initReactantMatchers();
    enumReactions1.push_back(rxn1);
    auto rxn1 = std::shared_ptr<ChemicalReaction>(
        RxnSmartsToChemicalReaction(fmt1.str()));
    ASSERT_INVARIANT(rxn1, "reaction is somehow bad");
    rxn1->initReactantMatchers();
    enumReactions1.push_back(rxn1);

    std::cerr << "  react: " << reactSmarts << std::endl;
    std::cerr << "  react enumerator: " << fmt2 << std::endl;
#endif
    // std::cerr << "  product: " << productSmarts << std::endl;
  }
  std::cerr << ">>>>  " << MolToSmiles(*dp_frame) << std::endl;
}

std::vector<size_t> LinkNodeOp::getVariationCounts() const {
  // std::vector<size_t> res(d_variationPoints.size());
  // std::transform(
  //     d_variationPoints.begin(), d_variationPoints.end(), res.begin(),
  //     [](std::pair<unsigned int, std::vector<unsigned int>> pr) -> size_t {
  //       return pr.second.size();
  //     });
  // return res;
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
  // we do the enumeration and attach one side in a single reaction.
  // we can't attach both sides because it's entirely possible that a
  // single atom has multiple attachment points and that would lead to
  // a mess.
  std::string reacts = "";
  std::string prods = "";
  for (size_t i = 0; i < which.size(); ++i) {
    auto variationIdx = i + 1;
    auto variationCount = d_pointRanges[i].first + which[i];
    auto reactFormat = boost::format("[%d#0]-[*:%d]-[%d#0].[%d#0][*:%d]") %
                       (variationIdx * 100) % (10 * variationIdx) %
                       (variationIdx * 100 + 10) % (variationIdx * 100 + 1) %
                       (10 * variationIdx + 1);
    auto react = reactFormat.str();
    if (!reacts.empty()) {
      reacts += ".";
    }
    reacts += reactFormat.str();
    auto prodFormat = boost::format("[*:%d]-[*:%d]") % (10 * variationIdx + 1) %
                      (10 * variationIdx);
    auto prod = prodFormat.str();
    for (size_t j = 1; j < variationCount; ++j) {
      prod += (boost::format("-[*:%d]") % variationIdx).str();
    }
    prod += (boost::format("-[%d#0]") % (variationIdx * 100 + 10)).str();
    if (!prods.empty()) {
      prods += ".";
    }
    prods += prod;
  }
  reacts = "(" + reacts + ")";
  prods = "(" + prods + ")";
  std::cerr << reacts << ">>" << prods << std::endl;
  RWMol *res = new RWMol(*dp_mol);
  return std::unique_ptr<ROMol>(static_cast<ROMol *>(res));
}

}  // namespace MolEnumerator

}  // namespace RDKit
