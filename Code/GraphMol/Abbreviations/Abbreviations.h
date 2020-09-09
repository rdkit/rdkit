//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_ABBREVIATIONS_H
#define RD_ABBREVIATIONS_H
#include <vector>
#include <string>
#include <memory>

namespace RDKit {
class ROMol;
class RWMol;

namespace Abbreviations {
RDKIT_ABBREVIATIONS_EXPORT struct AbbreviationDefinition {
  std::string llabel;
  std::string rlabel;
  std::string smarts;
  std::shared_ptr<ROMol> mol;
};
RDKIT_ABBREVIATIONS_EXPORT struct AbbreviationMatch {
  std::vector<std::pair<int, int>> match;
  AbbreviationDefinition abbrev;
  AbbreviationMatch(const std::vector<std::pair<int, int>>& matchArg,
                    const AbbreviationDefinition& abbrevArg)
      : match(matchArg), abbrev(abbrevArg){};
};
namespace common_properties {
RDKIT_ABBREVIATIONS_EXPORT extern const std::string numDummies;
}
namespace Utils {
RDKIT_ABBREVIATIONS_EXPORT std::vector<AbbreviationDefinition>
getDefaultAbbreviations();
RDKIT_ABBREVIATIONS_EXPORT std::vector<AbbreviationDefinition>
getDefaultLinkers();
RDKIT_ABBREVIATIONS_EXPORT std::vector<AbbreviationDefinition>
parseAbbreviations(const std::string& text, bool removeExtraDummies = false,
                   bool allowConnectionToDummies = false);
}  // namespace Utils

std::vector<AbbreviationMatch> findApplicableAbbreviationMatches(
    const ROMol& mol, const std::vector<AbbreviationDefinition>& abbrevs,
    double maxCoverage = 0.4);

}  // namespace Abbreviations
}  // namespace RDKit
#endif
