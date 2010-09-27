//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _HIERARCHCLUSTERPICKER_H
#define _HIERARCHCLUSTERPICKER_H

#include <RDGeneral/types.h>
#include "DistPicker.h"

namespace RDPickers {
  
  /*! \brief Diversity picker based on hierarchical clustering
   *  
   *  This class inherits from DistPicker since it uses the distance matrix
   *  for diversity picking. The clustering itself is done using the Murtagh 
   *  code in $RDBASE/Code/ML/Cluster/Mutagh/
   */
  class HierarchicalClusterPicker : public DistPicker {
  public:

    /*! \brief The type of hierarchical clustering algorithm to use
     */
    typedef enum {
      WARD=1,
      SLINK=2,
      CLINK=3,
      UPGMA=4,
      MCQUITTY=5,
      GOWER=6,
      CENTROID=7 } ClusterMethod;

    /*! \brief Constructor - takes a ClusterMethod as an argument
     *
     * Sets the hierarch clustering method
     */
    explicit HierarchicalClusterPicker(ClusterMethod clusterMethod) : d_method(clusterMethod) {;};

    /*! \brief This is the function that does the picking
     *
     * Here is how the algorithm works \n
     *  FIX: Supply reference
     *
     * - The entire pool is clustered using the distance matrix using one of the 
     *   hierachical clustering method (specified via the constructor). \n
     * - Starting with the individaul items in the pool, clusters are merged based 
     *   on the output from clustering method. \n
     * - The merging is stopped when the number of clusters is same as 
     *   the number of picks.
     * - For each item in a cluster the sum of square of the distances to the rest of
     *   of the items (in the cluster) is computed. The item with the smallest of values is
     *   picked as a representative of the cluster. Basically trying to pick the item closest
     *   to the centroid of the cluster. 
     *
     *
     *    \param distMat - distance matrix - a vector of double. It is assumed that only the 
     *              lower triangle element of the matrix are supplied in a 1D array\n
     *              NOTE: this matrix WILL BE ALTERED during the picking\n
     *    \param poolSize - the size of the pool to pick the items from. It is assumed that the
     *              distance matrix above contains the right number of elements; i.e.
     *              poolSize*(poolSize-1) \n
     *    \param pickSize - the number items to pick from pool (<= poolSize)
     */
    RDKit::INT_VECT pick(const double *distMat, unsigned int poolSize, unsigned int pickSize) const ;

    /*! \brief This is the function that does the clustering of the items - used by the picker
     *
     * ARGUMENTS:
     *
     *   \param distMat - distance matrix - a vector of double. It is assumed that only the 
     *              lower triangle element of the matrix are supplied in a 1D array\n
     *              NOTE: this matrix WILL BE ALTERED during the picking\n
     *   \param poolSize - the size of the pool to pick the items from. It is assumed that the
     *              distance matrix above contains the right number of elements; i.e.
     *              poolSize*(poolSize-1) \n
     *   \param pickSize - the number clusters to divide the pool into (<= poolSize)
     */
    RDKit::VECT_INT_VECT cluster(const double *distMat, unsigned int poolSize, unsigned int pickSize) const;

  private:
    ClusterMethod d_method;
  };
};

#endif
