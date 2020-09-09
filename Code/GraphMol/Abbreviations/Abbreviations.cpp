//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Abbreviations.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/types.h>
#include <RDGeneral/Invariant.h>

#include <boost/dynamic_bitset.hpp>
#include <iostream>

namespace RDKit {

namespace Abbreviations {

void applyMatches(RWMol& mol, const std::vector<AbbreviationMatch>& matches) {
  boost::dynamic_bitset<> toRemove(mol.getNumAtoms());
  for (const auto& amatch : matches) {
    // throughout this remember that atom 0 in the match is the dummy

    // convert atom 1 to be the abbreviation so that we don't have to
    // worry about messing up chirality, etc.
    auto connectIdx = amatch.match[1].second;
    auto connectingAtom = mol.getAtomWithIdx(connectIdx);
    connectingAtom->setProp(RDKit::common_properties::atomLabel,
                            amatch.abbrev.llabel);
    connectingAtom->setFormalCharge(0);
    connectingAtom->setAtomicNum(0);
    connectingAtom->setIsotope(0);
    // set the hybridization so these are drawn linearly
    connectingAtom->setHybridization(Atom::HybridizationType::SP);

    for (unsigned int i = 2; i < amatch.match.size(); ++i) {
      const auto& pr = amatch.match[i];
      CHECK_INVARIANT(!toRemove[pr.second], "overlapping matches");
      toRemove.set(pr.second);
      if (mol.getAtomWithIdx(pr.second)->getDegree() >
          amatch.abbrev.mol->getAtomWithIdx(pr.first)->getDegree()) {
        for (const auto& nbri : boost::make_iterator_range(
                 mol.getAtomNeighbors(mol.getAtomWithIdx(pr.second)))) {
          const auto& nbr = mol[nbri];
          auto nbrIdx = nbr->getIdx();
          // if this neighbor isn't in the match:
          if (!std::any_of(amatch.match.begin(), amatch.match.end(),
                           [&](const std::pair<int, int>& tpr) {
                             return tpr.second == rdcast<int>(nbrIdx);
                           })) {
            mol.addBond(nbrIdx, connectIdx, Bond::BondType::SINGLE);
          }
        }
      }
    }
  }
  for (unsigned int i = toRemove.size(); i > 0; --i) {
    if (toRemove[i - 1]) {
      mol.removeAtom(i - 1);
    }
  }
}

std::vector<AbbreviationMatch> findApplicableAbbreviationMatches(
    const ROMol& mol, const std::vector<AbbreviationDefinition>& abbrevs,
    double maxCoverage) {
  std::vector<AbbreviationMatch> res;
  auto nAtoms = mol.getNumAtoms();
  if (!nAtoms || abbrevs.empty()) {
    return res;
  }

  MolOps::fastFindRings(mol);

  std::vector<AbbreviationMatch> tres;
  boost::dynamic_bitset<> dummies(mol.getNumAtoms());
  boost::dynamic_bitset<> firstAts(mol.getNumAtoms());
  boost::dynamic_bitset<> covered(mol.getNumAtoms());

  for (const auto& abbrev : abbrevs) {
    if (maxCoverage > 0) {
      unsigned int nDummies;
      abbrev.mol->getProp(common_properties::numDummies, nDummies);
      if (double(abbrev.mol->getNumAtoms() - nDummies) / nAtoms >=
          maxCoverage) {
        continue;
      }
    }
    auto matches = SubstructMatch(mol, *abbrev.mol);
    for (const auto& match : matches) {
      CHECK_INVARIANT(match.size() > 1, "bad match size");
      // if we've already covered the first non-dummy atom or used it as a first
      // atom skip this.
      if (firstAts[match[1].second] || covered[match[1].second]) {
        continue;
      }
      bool keepIt = true;
      for (unsigned int i = 2; i < match.size(); ++i) {
        const auto& pr = match[i];
        if (covered[pr.second]) {
          keepIt = false;
          break;
        }
      }
      if (!keepIt) {
        continue;
      }
      for (unsigned int i = 1; i < match.size(); ++i) {
        const auto& pr = match[i];
        covered.set(pr.second);
      }
      dummies.set(match[0].second);
      firstAts.set(match[1].second);
      if (!firstAts[match[0].second]) {
        tres.emplace_back(match, abbrev);
      }
    }
  }
  for (const auto& itm : tres) {
    // if the dummy in this wasn't a first atom anywhere
    if (!firstAts[itm.match[0].second]) {
      res.push_back(std::move(itm));
    }
  }

  return res;
}

void condenseMolAbbreviations(
    RWMol& mol, const std::vector<AbbreviationDefinition>& abbrevs,
    double maxCoverage, bool sanitize) {
  auto applicable =
      findApplicableAbbreviationMatches(mol, abbrevs, maxCoverage);
  applyMatches(mol, applicable);
  if (sanitize) {
    MolOps::symmetrizeSSSR(mol);
  }
};

}  // namespace Abbreviations
}  // namespace RDKit
