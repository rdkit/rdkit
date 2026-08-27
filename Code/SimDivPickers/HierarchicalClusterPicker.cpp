//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "HierarchicalClusterPicker.h"
#include <RDGeneral/Invariant.h>
#include <RDGeneral/types.h>

typedef double real;
extern "C" void distdriver_(long int *n, long int *len, real *dists,
                            long int *toggle, long int *ia, long int *ib,
                            real *crit);

namespace RDPickers {

RDKit::VECT_INT_VECT HierarchicalClusterPicker::cluster(
    const double *distMat, unsigned int poolSize, unsigned int pickSize) const {
  PRECONDITION(distMat, "Invalid Distance Matrix");
  PRECONDITION((poolSize >= pickSize),
               "pickSize cannot be larger than the poolSize");

  // Do the clustering
  auto method = (long int)d_method;
  long int len = poolSize * (poolSize - 1);
  std::vector<long int> ia(poolSize, 0);
  std::vector<long int> ib(poolSize, 0);
  std::vector<real> crit(poolSize, 0.0);
  auto poolSize2 = static_cast<long int>(poolSize);

  distdriver_(&poolSize2,       // number of items in the pool
              &len,             // number of entries in the distance matrix
              (real *)distMat,  // distance matrix
              &method,          // the clustering method (ward, slink etc.)
              ia.data(),        // int vector with clustering history
              ib.data(),        // one more clustering history matrix
              crit.data()       // I believe this is a vector the difference in
                                // heights of two clusters
  );

  // we have the clusters now merge then until the number of clusters is same
  // as the number of picks we need
  // before we do that a bit of explanation on the vectors "ia" and "ib"
  //  - We with each item in the pool as an individual cluster
  //  - then we use the vectors ia and ib to merge them.
  //     ia and ib provides the ids of the clusters that need to be merged
  //     it is assumed that when a cluster ia[j] is merged with ib[j]
  //     ia[j] is replaced by the new cluster in the cluster list
  //
  RDKit::VECT_INT_VECT clusters;
  clusters.reserve(poolSize);
  for (unsigned int i = 0; i < poolSize; i++) {
    RDKit::INT_VECT cls = {static_cast<int>(i)};
    clusters.push_back(cls);
  }

  // do the merging, each round of of this loop eliminates one cluster
  RDKit::INT_VECT removed;
  removed.reserve(poolSize - pickSize);
  for (unsigned int i = 0; i < (poolSize - pickSize); i++) {
    int cx1 = ia[i] - 1;
    int cx2 = ib[i] - 1;

    // add the items from cluster cx2 to cx1
    // REVIEW: merge function???
    for (const auto &cx2i : clusters[cx2]) {
      clusters[cx1].push_back(cx2i);
    }

    // mark the second cluster as removed
    removed.push_back(cx2);
  }
  // sort removed so that looping will be easier later
  std::sort(removed.begin(), removed.end());

  // some error checking here, uniqueify removed and the vector should not
  // changed
  // REVIEW can we put this inside a #ifdef DEBUG?
  auto nEnd = std::unique(removed.begin(), removed.end());
  CHECK_INVARIANT(
      nEnd == removed.end(),
      "Somehow there are duplicates in the list of removed clusters");

  RDKit::VECT_INT_VECT res;
  unsigned int j = 0;
  for (unsigned int i = 0; i < poolSize; i++) {
    if (j < removed.size() && static_cast<int>(i) == removed[j]) {
      j++;
      continue;
    }
    res.push_back(clusters[i]);
  }
  return res;
}

RDKit::INT_VECT HierarchicalClusterPicker::pick(const double *distMat,
                                                unsigned int poolSize,
                                                unsigned int pickSize) const {
  PRECONDITION(distMat, "bad distance matrix");
  auto clusters = this->cluster(distMat, poolSize, pickSize);
  CHECK_INVARIANT(clusters.size() == pickSize, "");

  // the last step: find a representative element from each of the
  // remaining clusters
  RDKit::INT_VECT picks;
  for (unsigned int i = 0; i < pickSize; i++) {
    int pick;
    double minSumD2 = RDKit::MAX_DOUBLE;
    for (auto cxi1 = clusters[i].begin(); cxi1 != clusters[i].end(); ++cxi1) {
      int curPick = (*cxi1);
      double d2sum = 0.0;
      for (auto cxi2 = clusters[i].begin(); cxi2 != clusters[i].end(); ++cxi2) {
        if (cxi1 == cxi2) {
          continue;
        }
        double d = getDistFromLTM(distMat, curPick, (*cxi2));
        d2sum += (d * d);
      }
      if (d2sum < minSumD2) {
        pick = curPick;
        minSumD2 = d2sum;
      }
    }
    picks.push_back(pick);
  }
  return picks;
}
}  // namespace RDPickers
