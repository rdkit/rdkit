//
//  Copyright (C) 2003-2007 Greg Landrum and Rational Discovery LLC
//  Copyright (C) 2017 Greg Landrum and NextMove Software
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_MAXMINPICKER_H
#define RD_MAXMINPICKER_H

#include <RDGeneral/types.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/Exceptions.h>
#include <cstdlib>
#include "DistPicker.h"
#include <boost/random.hpp>
#include <random>

namespace RDPickers {

/*! \brief Implements the MaxMin algorithm for picking a subset of item from a
 *pool
 *
 *  This class inherits from the DistPicker and implements a specific picking
 *strategy
 *  aimed at diversity. See documentation for "pick()" member function for the
 *algorithm details
 */
class RDKIT_SIMDIVPICKERS_EXPORT MaxMinPicker : public DistPicker {
 public:
  /*! \brief Default Constructor
   *
   */
  MaxMinPicker(){};

  /*! \brief Contains the implementation for a lazy MaxMin diversity picker
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
                           unsigned int pickSize,
                           const RDKit::INT_VECT &firstPicks,
                           int seed = -1) const;

  template <typename T>
  RDKit::INT_VECT lazyPick(T &func, unsigned int poolSize,
                           unsigned int pickSize,
                           const RDKit::INT_VECT &firstPicks, int seed,
                           double &threshold) const;

  /*! \brief Contains the implementation for the MaxMin diversity picker
   *
   * Here is how the picking algorithm works, refer to
   *   Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604
   * for more detail:
   *
   * A subset of k items is to be selected from a pool containing N molecules.
   * Then the MaxMin method is as follows:
   *  -# Initialise Subset with some appropriately chosen seed
   *     compound and set x = 1.
   *  -# For each of the N-x remaining compounds in Dataset
   *     calculate its dissimilarity with each of the x compounds in
   *     Subset and retain the smallest of these x dissimilarities for
   *     each compound in Dataset.
   *  -# Select the molecule from Dataset with the largest value
   *     for the smallest dissimilarity calculated in Step 2 and
   *     transfer it to Subset.
   *  -# Set x = x + 1 and return to Step 2 if x < k.
   *
   *
   *
   *   \param distMat - distance matrix - a vector of double. It is assumed that
   *only the
   *              lower triangle element of the matrix are supplied in a 1D
   *array\n
   *   \param poolSize - the size of the pool to pick the items from. It is
   *assumed that the
   *              distance matrix above contains the right number of elements;
   *i.e.
   *              poolSize*(poolSize-1) \n
   *   \param pickSize - the number items to pick from pool (<= poolSize)
   *   \param firstPicks - indices of the items used to seed the pick set.
   *   \param seed - (optional) seed for the random number generator
   *                 If this is <0 the generator will be seeded with a
   *                 random number.
   */
  RDKit::INT_VECT pick(const double *distMat, unsigned int poolSize,
                       unsigned int pickSize, RDKit::INT_VECT firstPicks,
                       int seed = -1) const {
    CHECK_INVARIANT(distMat, "Invalid Distance Matrix");
    if (!poolSize) throw ValueErrorException("empty pool to pick from");
    if (poolSize < pickSize)
      throw ValueErrorException("pickSize cannot be larger than the poolSize");
    distmatFunctor functor(distMat);
    return this->lazyPick(functor, poolSize, pickSize, firstPicks, seed);
  }

  /*! \overload */
  RDKit::INT_VECT pick(const double *distMat, unsigned int poolSize,
                       unsigned int pickSize) const {
    RDKit::INT_VECT iv;
    return pick(distMat, poolSize, pickSize, iv);
  }
};

struct MaxMinPickInfo {
  double dist_bound;   // distance to closest reference
  unsigned int picks;  // number of references considered
  unsigned int next;   // singly linked list of candidates
};

// we implement this here in order to allow arbitrary functors without link
// errors
template <typename T>
RDKit::INT_VECT MaxMinPicker::lazyPick(T &func, unsigned int poolSize,
                                       unsigned int pickSize,
                                       const RDKit::INT_VECT &firstPicks,
                                       int seed, double &threshold) const {
  if (!poolSize) throw ValueErrorException("empty pool to pick from");

  if (poolSize < pickSize)
    throw ValueErrorException("pickSize cannot be larger than the poolSize");

  RDKit::INT_VECT picks;

  unsigned int memsize = (unsigned int)(poolSize * sizeof(MaxMinPickInfo));
  MaxMinPickInfo *pinfo = new MaxMinPickInfo[memsize];
  memset(pinfo, 0, memsize);

  picks.reserve(pickSize);
  unsigned int picked = 0;  // picks.size()
  unsigned int pick = 0;

  // pick the first entry
  if (firstPicks.empty()) {
    // get a seeded random number generator:
    typedef boost::mt19937 rng_type;
    typedef boost::uniform_int<> distrib_type;
    typedef boost::variate_generator<rng_type &, distrib_type> source_type;
    rng_type generator;
    distrib_type dist(0, poolSize - 1);
    if (seed >= 0) {
      generator.seed(static_cast<rng_type::result_type>(seed));
    } else {
      generator.seed(std::random_device()());
    }
    source_type randomSource(generator, dist);
    pick = randomSource();
    // add the pick to the picks
    picks.push_back(pick);
    // and remove it from the pool
    pinfo[pick].picks = 1;
    picked = 1;

  } else {
    for (RDKit::INT_VECT::const_iterator pIdx = firstPicks.begin();
         pIdx != firstPicks.end(); ++pIdx) {
      pick = static_cast<unsigned int>(*pIdx);
      if (pick >= poolSize) {
        delete[] pinfo;
        throw ValueErrorException("pick index was larger than the poolSize");
      }
      picks.push_back(pick);
      pinfo[pick].picks = 1;
      picked++;
    }
  }

  if (picked >= pickSize) {
    threshold = -1.0;
    delete[] pinfo;
    return picks;
  }

  unsigned int pool_list = 0;
  unsigned int *prev = &pool_list;
  // enter the pool into a list so that we can pick out of it easily
  for (unsigned int i = 0; i < poolSize; i++)
    if (pinfo[i].picks == 0) {
      *prev = i;
      prev = &pinfo[i].next;
    }
  *prev = 0;

  unsigned int poolIdx;
  unsigned int pickIdx;

  // First pass initialize dist_bound
  prev = &pool_list;
  pickIdx = picks[0];
  do {
    poolIdx = *prev;
    pinfo[poolIdx].dist_bound = func(poolIdx, pickIdx);
    pinfo[poolIdx].picks = 1;
    prev = &pinfo[poolIdx].next;
  } while (*prev != 0);

  // now pick 1 compound at a time
  double maxOFmin = -1.0;
  double tmpThreshold = -1.0;
  while (picked < pickSize) {
    unsigned int *pick_prev = nullptr;
    maxOFmin = -1.0;
    prev = &pool_list;
    do {
      poolIdx = *prev;
      double minTOi = pinfo[poolIdx].dist_bound;
      if (minTOi > maxOFmin) {
        unsigned int pi = pinfo[poolIdx].picks;
        while (pi < picked) {
          unsigned int picki = picks[pi];
          CHECK_INVARIANT(poolIdx != picki, "pool index != pick index");
          double dist = func(poolIdx, picki);
          pi++;
          if (dist <= minTOi) {
            minTOi = dist;
            if (minTOi <= maxOFmin) break;
          }
        }
        pinfo[poolIdx].dist_bound = minTOi;
        pinfo[poolIdx].picks = pi;
        if (minTOi > maxOFmin) {
          maxOFmin = minTOi;
          pick_prev = prev;
          pick = poolIdx;
        }
      }
      prev = &pinfo[poolIdx].next;
    } while (*prev != 0);

    // if the current distance is closer then threshold, we're done
    if (maxOFmin <= threshold && threshold >= 0.0) break;
    tmpThreshold = maxOFmin;
    // now add the new pick to picks and remove it from the pool
    *pick_prev = pinfo[pick].next;
    picks.push_back(pick);
    picked++;
  }

  threshold = tmpThreshold;
  delete[] pinfo;
  return picks;
}

template <typename T>
RDKit::INT_VECT MaxMinPicker::lazyPick(T &func, unsigned int poolSize,
                                       unsigned int pickSize,
                                       const RDKit::INT_VECT &firstPicks,
                                       int seed) const {
  double threshold = -1.0;
  return MaxMinPicker::lazyPick(func, poolSize, pickSize, firstPicks, seed,
                                threshold);
}

template <typename T>
RDKit::INT_VECT MaxMinPicker::lazyPick(T &func, unsigned int poolSize,
                                       unsigned int pickSize) const {
  RDKit::INT_VECT firstPicks;
  double threshold = -1.0;
  int seed = -1;
  return MaxMinPicker::lazyPick(func, poolSize, pickSize, firstPicks, seed,
                                threshold);
}
};  // namespace RDPickers

#endif
