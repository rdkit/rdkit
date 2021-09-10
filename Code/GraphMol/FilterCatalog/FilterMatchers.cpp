//  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
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

#include "FilterMatchers.h"
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <utility>

namespace RDKit {
const char *DEFAULT_FILTERMATCHERBASE_NAME = "Unnamed FilterMatcherBase";
const char *SMARTS_MATCH_NAME_DEFAULT = "Unnamed SmartsMatcher";
namespace {
const int debugParse = 0;
const bool mergeHs = true;
}  // namespace
SmartsMatcher::SmartsMatcher(const ROMol &pattern, unsigned int minCount,
                             unsigned int maxCount)
    : FilterMatcherBase(SMARTS_MATCH_NAME_DEFAULT),
      d_pattern(new ROMol(pattern)),
      d_min_count(minCount),
      d_max_count(maxCount) {}

SmartsMatcher::SmartsMatcher(const std::string &name, const ROMol &pattern,
                             unsigned int minCount, unsigned int maxCount)
    : FilterMatcherBase(name),
      d_pattern(new ROMol(pattern)),
      d_min_count(minCount),
      d_max_count(maxCount) {}

SmartsMatcher::SmartsMatcher(const std::string &name, const std::string &smarts,
                             unsigned int minCount, unsigned int maxCount)
    : FilterMatcherBase(name),
      d_pattern(SmartsToMol(smarts, debugParse, mergeHs)),
      d_min_count(minCount),
      d_max_count(maxCount) {}

SmartsMatcher::SmartsMatcher(const std::string &name, ROMOL_SPTR pattern,
                             unsigned int minCount, unsigned int maxCount)
    : FilterMatcherBase(name),
      d_pattern(std::move(pattern)),
      d_min_count(minCount),
      d_max_count(maxCount) {}

void SmartsMatcher::setPattern(const std::string &smarts) {
  d_pattern.reset(SmartsToMol(smarts, debugParse, mergeHs));
}

void SmartsMatcher::setPattern(const ROMol &mol) {
  d_pattern.reset(new ROMol(mol));
}

SmartsMatcher::SmartsMatcher(const SmartsMatcher &rhs)
    : FilterMatcherBase(rhs),
      d_pattern(rhs.d_pattern),
      d_min_count(rhs.d_min_count),
      d_max_count(rhs.d_max_count) {}

bool SmartsMatcher::getMatches(const ROMol &mol,
                               std::vector<FilterMatch> &matchVect) const {
  PRECONDITION(d_pattern.get(), "bad on pattern");

  bool onPatExists = false;
  std::vector<RDKit::MatchVectType> matches;

  if (d_min_count == 1 && d_max_count == UINT_MAX) {
    RDKit::MatchVectType match;
    onPatExists = RDKit::SubstructMatch(mol, *d_pattern.get(), match);
    if (onPatExists) {
      matchVect.emplace_back(copy(), match);
    }
  } else {  // need to count
    const bool uniquify = true;
    unsigned int count =
        RDKit::SubstructMatch(mol, *d_pattern.get(), matches, uniquify);
    onPatExists = (count >= d_min_count &&
                   (d_max_count == UINT_MAX || count <= d_max_count));
    if (onPatExists) {
      boost::shared_ptr<FilterMatcherBase> clone = copy();
      for (auto &match : matches) {
        matchVect.emplace_back(clone, match);
      }
    }
  }
  return onPatExists;
}

bool SmartsMatcher::hasMatch(const ROMol &mol) const {
  PRECONDITION(d_pattern.get(), "bad on pattern");

  if (d_min_count == 1 && d_max_count == UINT_MAX) {
    RDKit::MatchVectType matches;
    return SubstructMatch(mol, *d_pattern.get(), matches);
  } else {  // need to count
    const bool uniquify = true;
    std::vector<RDKit::MatchVectType> matches;
    unsigned int count =
        RDKit::SubstructMatch(mol, *d_pattern.get(), matches, uniquify);
    return (count >= d_min_count &&
            (d_max_count == UINT_MAX || count <= d_max_count));
  }
}

bool FilterHierarchyMatcher::getMatches(const ROMol &mol,
                                        std::vector<FilterMatch> &m) const {
  // a null matcher is root, just goes to the children

  std::vector<FilterMatch> temp;
  bool result = d_matcher->getMatches(mol, temp);

  if (result) {
    std::vector<FilterMatch> children;

    for (auto matcher : d_children) {
      matcher->getMatches(mol, children);
    }

    if (children.size()) {
      m.insert(m.end(), children.begin(), children.end());
    } else {
      m.insert(m.end(), temp.begin(), temp.end());
    }
  }

  return result;
}
}  // namespace RDKit
