//
//  Copyright (C) 2003-2019 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "SubstructUtils.h"
#include <RDGeneral/utils.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <boost/dynamic_bitset.hpp>

namespace RDKit {

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
  std::cerr << "\t\tatomCompat: " << a1 << " " << a1->getIdx() << "-" << a2
            << " " << a2->getIdx() << std::endl;
  std::cerr << "\t\t    " << res << std::endl;
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
    if (!b1->getBeginAtom()->Match(b1->getBeginAtom()) ||
        !b1->getEndAtom()->Match(b2->getEndAtom())) {
      res = false;
    }
  }
  // std::cerr << "\t\tbondCompat: " << b1->getIdx() << "-" << b2->getIdx() <<
  // ":"
  //           << res << std::endl;
  return res;
}

void removeDuplicates(std::vector<MatchVectType> &v, unsigned int nAtoms) {
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
  std::vector<boost::dynamic_bitset<>> seen;
  std::vector<MatchVectType> res;
  for (std::vector<MatchVectType>::const_iterator i = v.begin(); i != v.end();
       ++i) {
    boost::dynamic_bitset<> val(nAtoms);
    for (const auto &ci : *i) {
      val.set(ci.second);
    }
    if (std::find(seen.begin(), seen.end(), val) == seen.end()) {
      // it's something new
      res.push_back(*i);
      seen.push_back(val);
    }
  }
  v = res;
}
}  // namespace RDKit
