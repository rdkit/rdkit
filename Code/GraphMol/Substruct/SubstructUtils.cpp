//
//  Copyright (C) 2003-2021 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "SubstructUtils.h"
#include <set>
#include <RDGeneral/utils.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/Substruct/SubstructUtils.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/dynamic_bitset.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {

namespace detail {
// Helper class used by the sortMatchesByDegreeOfCoreSubstitution
// and getMostSubstitutedCoreMatch functions. A penalty of 1.0 is assigned
// to matches for each terminal dummy atom matching hydrogen.
// To make the sort stable in case of ties, a fraction of 1.0
// is added to each score based on match indices.
class ScoreMatchesByDegreeOfCoreSubstitution {
 public:
  typedef std::pair<unsigned int, double> IdxScorePair;
  ScoreMatchesByDegreeOfCoreSubstitution(
      const RDKit::ROMol &mol, const RDKit::ROMol &query,
      const std::vector<RDKit::MatchVectType> &matches)
      : d_mol(mol),
        d_query(query),
        d_matches(matches),
        d_sumIndices(0.0),
        d_minIdx(-1),
        d_isSorted(false) {
    PRECONDITION(!matches.empty(), "matches must not be empty");
    auto na = d_mol.getNumAtoms();
    d_sumIndices = static_cast<double>(na * (na + 1) / 2);
    unsigned int i = 0;
    d_matchIdxVsScore.reserve(d_matches.size());
    for (const auto &match : d_matches) {
      d_matchIdxVsScore.emplace_back(i++, computeScore(match));
    }
  }
  const RDKit::MatchVectType &getMostSubstitutedCoreMatch() {
    if (d_minIdx == -1) {
      d_minIdx = std::min_element(d_matchIdxVsScore.begin(),
                                  d_matchIdxVsScore.end(), compare)
                     ->first;
    }
    return d_matches.at(d_minIdx);
  }
  std::vector<MatchVectType> sortMatchesByDegreeOfCoreSubstitution() {
    if (!d_isSorted) {
      std::sort(d_matchIdxVsScore.begin(), d_matchIdxVsScore.end(), compare);
      d_isSorted = true;
      d_minIdx = d_matchIdxVsScore.front().first;
    }
    std::vector<MatchVectType> res(d_matches.size());
    std::transform(
        d_matchIdxVsScore.begin(), d_matchIdxVsScore.end(), res.begin(),
        [this](const IdxScorePair &pair) { return d_matches.at(pair.first); });
    return res;
  }

 private:
  static bool compare(const IdxScorePair &aPair, const IdxScorePair &bPair) {
    return (aPair.second < bPair.second);
  }
  bool doesRGroupMatchHydrogen(const std::pair<int, int> &pair) const {
    const auto queryAtom = d_query.getAtomWithIdx(pair.first);
    const auto molAtom = d_mol.getAtomWithIdx(pair.second);
    return (molAtom->getAtomicNum() == 1 &&
            isAtomTerminalRGroupOrQueryHydrogen(queryAtom));
  }
  double computeScore(const RDKit::MatchVectType &match) const {
    double penalty = 0.0;
    double i = 0.0;
    for (const auto &pair : match) {
      i += static_cast<double>(pair.second);
      if (doesRGroupMatchHydrogen(pair)) {
        penalty += 1.0;
      }
    }
    penalty += i / d_sumIndices;
    return penalty;
  }
  const RDKit::ROMol &d_mol;
  const RDKit::ROMol &d_query;
  const std::vector<RDKit::MatchVectType> &d_matches;
  std::vector<IdxScorePair> d_matchIdxVsScore;
  double d_sumIndices;
  int d_minIdx;
  bool d_isSorted;
};
}  // namespace detail

bool atomCompat(const Atom *a1, const Atom *a2,
                const SubstructMatchParameters &ps) {
  PRECONDITION(a1, "bad atom");
  PRECONDITION(a2, "bad atom");
  // std::cerr << "\t\tatomCompat: "<< a1 << " " << a1->getIdx() << "-" << a2 <<
  // " " << a2->getIdx() << std::endl;
  bool res;
  if (ps.useQueryQueryMatches && a1->hasQuery() && a2->hasQuery()) {
    res = static_cast<const QueryAtom *>(a1)->QueryMatch(
        static_cast<const QueryAtom *>(a2));
  } else {
    res = a1->Match(a2);
  }
  return res;
}

bool chiralAtomCompat(const Atom *&a1, const Atom *&a2) {
  PRECONDITION(a1, "bad atom");
  PRECONDITION(a2, "bad atom");
  bool res = a1->Match(a2);
  if (res) {
    std::string s1, s2;
    bool hascode1 = a1->getPropIfPresent(common_properties::_CIPCode, s1);
    bool hascode2 = a2->getPropIfPresent(common_properties::_CIPCode, s2);
    if (hascode1 || hascode2) {
      res = hascode1 && hascode2 && s1 == s2;
    }
  }
  std::cerr << "\t\tchiralAtomCompat: " << a1 << " " << a1->getIdx() << "-"
            << a2 << " " << a2->getIdx() << std::endl;
  std::cerr << "\t\t    " << res << std::endl;
  return res;
}

bool bondCompat(const Bond *b1, const Bond *b2,
                const SubstructMatchParameters &ps) {
  PRECONDITION(b1, "bad bond");
  PRECONDITION(b2, "bad bond");
  bool res;
  if (ps.useQueryQueryMatches && b1->hasQuery() && b2->hasQuery()) {
    res = static_cast<const QueryBond *>(b1)->QueryMatch(
        static_cast<const QueryBond *>(b2));
  } else if (ps.aromaticMatchesConjugated && !b1->hasQuery() &&
             !b2->hasQuery() &&
             ((b1->getBondType() == Bond::AROMATIC &&
               b2->getBondType() == Bond::AROMATIC) ||
              (b1->getBondType() == Bond::AROMATIC && b2->getIsConjugated()) ||
              (b2->getBondType() == Bond::AROMATIC && b1->getIsConjugated()))) {
    res = true;
  } else {
    res = b1->Match(b2);
  }
  if (res && b1->getBondType() == Bond::DATIVE &&
      b2->getBondType() == Bond::DATIVE) {
    // for dative bonds we need to make sure that the direction also matches:
    if (!b1->getBeginAtom()->Match(b2->getBeginAtom()) ||
        !b1->getEndAtom()->Match(b2->getEndAtom())) {
      res = false;
    }
  }
  // std::cerr << "\t\tbondCompat: " << b1->getIdx() << "-" << b2->getIdx() <<
  // ":"
  //           << res << std::endl;
  return res;
}

void removeDuplicates(std::vector<MatchVectType> &matches,
                      unsigned int nAtoms) {
  //
  //  This works by tracking the indices of the atoms in each match vector.
  //  This can lead to unexpected behavior when looking at rings and queries
  //  that don't specify bond orders.  For example querying this molecule:
  //    C1CCC=1
  //  with the pattern constructed from SMARTS C~C~C~C will return a
  //  single match, despite the fact that there are 4 different paths
  //  when valence is considered.  The defense of this behavior is
  //  that the 4 paths are equivalent in the semantics of the query.
  //  Also, OELib returns the same results
  //
  std::set<boost::dynamic_bitset<>> seen;
  std::vector<MatchVectType> res;
  res.reserve(matches.size());
  for (auto &&match : matches) {
    boost::dynamic_bitset<> val(nAtoms);
    for (const auto &ci : match) {
      val.set(ci.second);
    }
    auto pos = seen.lower_bound(val);
    if (pos == seen.end() || *pos != val) {
      res.push_back(std::move(match));
      seen.insert(pos, std::move(val));
    }
  }
  res.shrink_to_fit();
  matches = std::move(res);
}

const MatchVectType &getMostSubstitutedCoreMatch(
    const ROMol &mol, const ROMol &core,
    const std::vector<MatchVectType> &matches) {
  detail::ScoreMatchesByDegreeOfCoreSubstitution matchScorer(mol, core,
                                                             matches);
  return matchScorer.getMostSubstitutedCoreMatch();
}

std::vector<MatchVectType> sortMatchesByDegreeOfCoreSubstitution(
    const ROMol &mol, const ROMol &core,
    const std::vector<MatchVectType> &matches) {
  detail::ScoreMatchesByDegreeOfCoreSubstitution matchScorer(mol, core,
                                                             matches);
  return matchScorer.sortMatchesByDegreeOfCoreSubstitution();
}

bool isAtomTerminalRGroupOrQueryHydrogen(const Atom *atom) {
  return atom->getDegree() == 1 &&
         (atom->getAtomicNum() == 0 ||
          (atom->hasQuery() &&
           describeQuery(atom).find("AtomAtomicNum 1 = val") !=
               std::string::npos));
}

#define PT_OPT_GET(opt) params.opt = pt.get(#opt, params.opt)
#define PT_OPT_PUT(opt) pt.put(#opt, params.opt);

void updateSubstructMatchParamsFromJSON(SubstructMatchParameters &params,
                                        const std::string &json) {
  if (json.empty()) {
    return;
  }
  std::istringstream ss;
  ss.str(json);
  boost::property_tree::ptree pt;
  boost::property_tree::read_json(ss, pt);
  PT_OPT_GET(useChirality);
  PT_OPT_GET(useEnhancedStereo);
  PT_OPT_GET(aromaticMatchesConjugated);
  PT_OPT_GET(useQueryQueryMatches);
  PT_OPT_GET(recursionPossible);
  PT_OPT_GET(uniquify);
  PT_OPT_GET(maxMatches);
  PT_OPT_GET(numThreads);
}

std::string substructMatchParamsToJSON(const SubstructMatchParameters &params) {
  boost::property_tree::ptree pt;

  PT_OPT_PUT(useChirality);
  PT_OPT_PUT(useEnhancedStereo);
  PT_OPT_PUT(aromaticMatchesConjugated);
  PT_OPT_PUT(useQueryQueryMatches);
  PT_OPT_PUT(recursionPossible);
  PT_OPT_PUT(uniquify);
  PT_OPT_PUT(maxMatches);
  PT_OPT_PUT(numThreads);

  std::stringstream ss;
  boost::property_tree::json_parser::write_json(ss, pt);
  return ss.str();
}

}  // namespace RDKit
