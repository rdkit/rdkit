//
//  Copyright (C) 2003-2007 Greg Landrum and Rational Discovery LLC
//  Copyright (C) 2017-2019 Greg Landrum and NextMove Software
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_LEADERPICKER_H
#define RD_LEADERPICKER_H

#include <RDGeneral/types.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/Exceptions.h>
#include <cstdlib>
#include "DistPicker.h"

namespace RDPickers {

/*! \brief Implements the Leader algorithm for picking a subset of item from a
 *pool
 *
 *  This class inherits from the DistPicker and implements a specific picking
 *strategy
 *  aimed at diversity. See documentation for "pick()" member function for the
 *algorithm details
 */
class RDKIT_SIMDIVPICKERS_EXPORT LeaderPicker : public DistPicker {
 public:
  double default_threshold;
  int default_nthreads;

  /*! \brief Default Constructor
   *
   */
  LeaderPicker() : default_threshold(0.0), default_nthreads(1) {}
  LeaderPicker(double threshold)
      : default_threshold(threshold), default_nthreads(1) {}
  LeaderPicker(double threshold, int nthreads)
      : default_threshold(threshold), default_nthreads(nthreads) {}

  /*! \brief Contains the implementation for a lazy Leader diversity picker
   *
   * See the documentation for the pick() method for details about the algorithm
   *
   *   \param func - a function (or functor) taking two unsigned ints as
   *arguments
   *              and returning the distance (as a double) between those two
   *elements.
   *   \param poolSize - the size of the pool to pick the items from. It is
   *assumed that the
   *              distance matrix above contains the right number of elements;
   *i.e.
   *              poolSize*(poolSize-1)
   *   \param pickSize - the number items to pick from pool (<= poolSize)
   *   \param firstPicks - (optional)the first items in the pick list
   *   \param seed - (optional) seed for the random number generator.
   *                 If this is <0 the generator will be seeded with a
   *                 random number.
   */
  template <typename T>
  RDKit::INT_VECT lazyPick(T &func, unsigned int poolSize,
                           unsigned int pickSize) const;

  template <typename T>
  RDKit::INT_VECT lazyPick(T &func, unsigned int poolSize,
                           unsigned int pickSize, double threshold) const;

  template <typename T>
  RDKit::INT_VECT lazyPick(T &func, unsigned int poolSize,
                           unsigned int pickSize,
                           const RDKit::INT_VECT &firstPicks,
                           double threshold) const;

  template <typename T>
  RDKit::INT_VECT lazyPick(T &func, unsigned int poolSize,
                           unsigned int pickSize,
                           const RDKit::INT_VECT &firstPicks, double threshold,
                           int nthreads) const;

  /*! \brief Contains the implementation for the Leader diversity picker
   *   \param distMat - distance matrix - a vector of double. It is assumed that
   *only the
   *              lower triangle element of the matrix are supplied in a 1D
   *array\n
   *   \param poolSize - the size of the pool to pick the items from. It is
   *assumed that the
   *              distance matrix above contains the right number of elements;
   *i.e.
   *              poolSize*(poolSize-1) \n
   *   \param pickSize - maximum number items to pick from pool (<= poolSize)
   *   \param firstPicks - indices of the items used to seed the pick set.
   */
  RDKit::INT_VECT pick(const double *distMat, unsigned int poolSize,
                       unsigned int pickSize, const RDKit::INT_VECT &firstPicks,
                       double threshold, int nthreads) const {
    CHECK_INVARIANT(distMat, "Invalid Distance Matrix");
    if (!poolSize) throw ValueErrorException("empty pool to pick from");
    if (poolSize < pickSize)
      throw ValueErrorException("pickSize cannot be larger than the poolSize");
    distmatFunctor functor(distMat);
    return this->lazyPick(functor, poolSize, pickSize, firstPicks, threshold,
                          nthreads);
  }

  /*! \overload */
  RDKit::INT_VECT pick(const double *distMat, unsigned int poolSize,
                       unsigned int pickSize) const {
    RDKit::INT_VECT iv;
    return pick(distMat, poolSize, pickSize, iv, default_threshold,
                default_nthreads);
  }
};

template <typename T>
struct LeaderPickerState {
  std::vector<int> v;
  unsigned int left;
  double threshold;
  int query;
  T *func;

  LeaderPickerState(unsigned int count, int) {
    v.resize(count);
    for (unsigned int i = 0; i < count; i++) v[i] = i;
    left = count;
  }

  bool empty() { return left == 0; }

  unsigned int compact(int *dst, int *src, unsigned int len) {
    unsigned int count = 0;
    for (unsigned int i = 0; i < len; i++) {
      double ld = (*func)(query, src[i]);
      std::cerr << query << "-" << src[i] << " " << ld << std::endl;
      if (ld > threshold) dst[count++] = src[i];
    }
    return count;
  }

  void compact(int pick) {
    query = pick;
    left = compact(&v[0], &v[0], left);
  }

  int compact_next() {
    query = v[0];
    left = compact(&v[0], &v[1], left - 1);
    return query;
  }
};

// we implement this here in order to allow arbitrary functors without link
// errors
template <typename T>
RDKit::INT_VECT LeaderPicker::lazyPick(T &func, unsigned int poolSize,
                                       unsigned int pickSize,
                                       const RDKit::INT_VECT &firstPicks,
                                       double threshold, int nthreads) const {
  if (!poolSize) throw ValueErrorException("empty pool to pick from");

  if (poolSize < pickSize)
    throw ValueErrorException("pickSize cannot be larger than the poolSize");

  RDKit::INT_VECT picks;

  LeaderPickerState<T> stat(poolSize, nthreads);
  stat.threshold = threshold;
  stat.func = &func;

  unsigned int picked = 0;  // picks.size()
  unsigned int pick = 0;

  if (!firstPicks.empty()) {
    for (RDKit::INT_VECT::const_iterator pIdx = firstPicks.begin();
         pIdx != firstPicks.end(); ++pIdx) {
      pick = static_cast<unsigned int>(*pIdx);
      if (pick >= poolSize) {
        throw ValueErrorException("pick index was larger than the poolSize");
      }
      picks.push_back(pick);
      stat.compact(pick);
      picked++;
    }
  }

  while (picked < pickSize && !stat.empty()) {
    pick = stat.compact_next();
    picks.push_back(pick);
    picked++;
    printf("\t %u %u %u\n", picked, pick, stat.left);
  }
  return picks;
}

template <typename T>
RDKit::INT_VECT LeaderPicker::lazyPick(T &func, unsigned int poolSize,
                                       unsigned int pickSize) const {
  RDKit::INT_VECT firstPicks;
  return LeaderPicker::lazyPick(func, poolSize, pickSize, firstPicks,
                                default_threshold, default_nthreads);
}

template <typename T>
RDKit::INT_VECT LeaderPicker::lazyPick(T &func, unsigned int poolSize,
                                       unsigned int pickSize,
                                       double threshold) const {
  RDKit::INT_VECT firstPicks;
  return LeaderPicker::lazyPick(func, poolSize, pickSize, firstPicks, threshold,
                                default_nthreads);
}
template <typename T>
RDKit::INT_VECT LeaderPicker::lazyPick(T &func, unsigned int poolSize,
                                       unsigned int pickSize,
                                       const RDKit::INT_VECT &firstPicks,
                                       double threshold) const {
  return LeaderPicker::lazyPick(func, poolSize, pickSize, firstPicks, threshold,
                                default_nthreads);
}

};  // namespace RDPickers

#endif
