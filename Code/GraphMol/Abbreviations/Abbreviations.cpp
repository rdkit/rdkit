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
#include <boost/dynamic_bitset.hpp>
#include <iostream>

namespace RDKit {

namespace Abbreviations {

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
}  // namespace Abbreviations
}  // namespace RDKit
