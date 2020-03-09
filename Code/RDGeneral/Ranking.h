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
#include <boost/foreach.hpp>
#include <cstdint>

namespace Rankers {
//! functor for implementing > on two std::pairs.  The first entries are
// compared.
template <typename T1, typename T2>
struct pairGreater
    : public std::binary_function<std::pair<T1, T2>, std::pair<T1, T2>, bool> {
  bool operator()(const std::pair<T1, T2> &v1,
                  const std::pair<T1, T2> &v2) const {
    return v1.first > v2.first;
  }
};

//! function for implementing < on two std::pairs.  The first entries are
// compared.
template <typename T1, typename T2>
struct pairLess
    : public std::binary_function<std::pair<T1, T2>, std::pair<T1, T2>, bool> {
  bool operator()(const std::pair<T1, T2> &v1,
                  const std::pair<T1, T2> &v2) const {
    return v1.first < v2.first;
  }
};

template <typename T>
class argless : public std::binary_function<T, T, bool> {
 public:
  argless(const T &c) : std::binary_function<T, T, bool>(), container(c){};
  bool operator()(unsigned int v1, unsigned int v2) const {
    return container[v1] < container[v2];
  }
  const T &container;
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
  for (unsigned int i = 0; i < nEntries; ++i) indices[i] = i;
  std::sort(indices.begin(), indices.end(), argless<std::vector<T1>>(vect));

  int currRank = 0;
  T1 lastV = vect[indices[0]];
  BOOST_FOREACH (unsigned int idx, indices) {
    T1 v = vect[idx];
    if (v == lastV) {
      res[idx] = currRank;
    } else {
      res[idx] = ++currRank;
      lastV = v;
    }
  }
}
}  // namespace Rankers
#endif
