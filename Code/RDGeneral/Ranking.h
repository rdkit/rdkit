//
//  Copyright (C) 2004-2015 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

//! \file Ranking.h
/*!
    \brief Utility functionality used to rank sequences

    Much of this used to be in GraphMol/RankAtoms.h
*/
#include <RDGeneral/export.h>
#ifndef RD_RANKING_H
#define RD_RANKING_H

#include <vector>
#include <functional>
#include <algorithm>
#include <cstdint>

namespace Rankers {
inline auto pairGreater = [](const auto &v1, const auto &v2) {
  return v1.first > v2.first;
};

//! function for implementing < on two std::pairs.  The first entries are
/// compared.
inline auto pairLess = [](const auto &v1, const auto &v2) {
  return v1.first < v2.first;
};

//! ranks the entries in a vector
/*!
  \param vect the vector to rank
  \param res  is used to return the ranks of each entry
*/
template <typename T1, typename T2>
void rankVect(const std::vector<T1> &vect, T2 &res) {
  PRECONDITION(res.size() >= vect.size(), "vector size mismatch");
  unsigned int nEntries = rdcast<unsigned int>(vect.size());

  std::vector<unsigned int> indices(nEntries);
  for (unsigned int i = 0; i < nEntries; ++i) {
    indices[i] = i;
  }
  std::sort(indices.begin(), indices.end(),
            [&](auto i1, auto i2) { return vect[i1] < vect[i2]; });

  int currRank = 0;
  unsigned int lastIdx = indices[0];
  for (auto idx : indices) {
    if (vect[idx] == vect[lastIdx]) {
      res[idx] = currRank;
    } else {
      res[idx] = ++currRank;
      lastIdx = idx;
    }
  }
}
}  // namespace Rankers
#endif
