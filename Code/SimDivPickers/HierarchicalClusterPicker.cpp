// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include "HierarchicalClusterPicker.h"
#include <RDGeneral/Invariant.h>
#include <RDGeneral/types.h>


typedef double real;
extern "C" void distdriver_(int *n,int *len,
                 real *dists,
                 int *toggle,
                 int *ia,int *ib,real *crit);

namespace RDPickers {

  RDKit::VECT_INT_VECT HierarchicalClusterPicker::cluster(const double *distMat,
							  unsigned int poolSize,
							  unsigned int pickSize) {
    PRECONDITION(distMat, "Invalid Distance Matrix");
    PRECONDITION((poolSize >= pickSize),
		 "pickSize cannot be larger than the poolSize");

    // Do the clustering 
    int *ia, *ib;
    real *crit;
    int method = (int)d_method;
    int len = poolSize*(poolSize-1);
    ia = (int *)calloc(poolSize, sizeof(int));
    ib = (int *)calloc(poolSize, sizeof(int));
    crit = (real *)calloc(poolSize,sizeof(real));
    CHECK_INVARIANT(ia,"failed to allocate memory");
    CHECK_INVARIANT(ib,"failed to allocate memory");
    CHECK_INVARIANT(crit,"failed to allocate memory");
    
    distdriver_(reinterpret_cast<int *>(&poolSize), // number of items in the pool
                &len, // number of entries in the distance matrix
                (real *)distMat, // distance matrix
                &method, // the clustering method (ward, slink etc.)
                ia, // int vector with clustering history
                ib, // one more clustering history matrix
                crit // I believe this is a vector the difference in heights of two clusters
                );

    // we have the clusters now merge then untill the number of clusters is same
    // as the number of picks we need
    // before we do that a bit of explanation on the vectors "ia" and "ib"
    //  - We with each item in the pool as an individual cluster
    //  - then we use the vectors ia and ib to merge them.
    //     ia and ib provides the ids of the clusters that need to be merged
    //     it is assumed that when a cluster ia[j] is merged with ib[j] 
    //     ia[j] is replaced by the new cluster in the cluster list
    // 
    RDKit::VECT_INT_VECT clusters;
    unsigned int i;
    for (i = 0; i < poolSize; i++) {
      RDKit::INT_VECT cls;
      cls.push_back(i);
      clusters.push_back(cls);
    }
    
    int cx1, cx2;
    RDKit::INT_VECT removed;

    // do the merging, each round of of this loop eleminates one cluster
    for (i = 0; i < (poolSize - pickSize); i++) {
      cx1 = ia[i] - 1;
      cx2 = ib[i] - 1;

      // add the items from cluster cx2 to cx1
      RDKit::INT_VECT_CI cx2i;
      // REVIEW: merge function???
      for (cx2i = clusters[cx2].begin(); cx2i != clusters[cx2].end(); cx2i++) {
        clusters[cx1].push_back(*cx2i);
      }
      
      // mark the second cluster as removed
      removed.push_back(cx2);
    }
    
    free(ia);
    free(ib);
    free(crit);

    // sort removed so that looping will ve easier later
    std::sort(removed.begin(), removed.end());

    //some error checking here, uniqueify removed and the vector should not changed
    // REVIEW can we put this inside a #ifdef DEBUG?
    RDKit::INT_VECT_CI nEnd = std::unique(removed.begin(), removed.end());
    CHECK_INVARIANT(nEnd == removed.end(), "Somehow there are duplicates in the list of removed clusters");

    RDKit::VECT_INT_VECT res;
    unsigned int j = 0;
    for (i = 0; i < poolSize; i++) {
      if (i == removed[j]) {
        j++;
        continue;
      }
      res.push_back(clusters[i]);
    }
    return res;
  }
    
  RDKit::INT_VECT HierarchicalClusterPicker::pick(const double *distMat, 
                                                  unsigned int poolSize,
						  unsigned int pickSize) {
    PRECONDITION(distMat,"bad distance matrix");
    RDKit::VECT_INT_VECT clusters = this->cluster(distMat, poolSize, pickSize);

    CHECK_INVARIANT(clusters.size() == pickSize, "");
    // the last step find the representative element each of the remaining clusters
    RDKit::INT_VECT picks;
    unsigned int i;
    for (i = 0; i < pickSize; i++) {
      int pick, curPick;
      double minSumD2 = RDKit::MAX_DOUBLE;
      RDKit::INT_VECT_CI cxi1, cxi2; 
      for (cxi1 = clusters[i].begin(); cxi1 != clusters[i].end(); cxi1++) {
        curPick = (*cxi1);
        double d, d2sum = 0.0;
        for (cxi2 = clusters[i].begin(); cxi2 != clusters[i].end(); cxi2++) {
          if (cxi1 == cxi2) {
            continue;
          }
          d = getDistFromLTM(distMat, curPick, (*cxi2));
          d2sum += (d*d);
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
}
    

      
    
