//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef _RD_MAXMINPICKER_H_
#define _RD_MAXMINPICKER_H_

#include <RDGeneral/types.h>
#include "DistPicker.h"

namespace RDPickers {

  /*! \brief Implements the MaxMin algorithm for picking a subset of item from a pool
   *
   *  This class inherits from the DistPicker and implements a specific picking strategy
   *  aimed at diversity. See documentation for "pick()" member function for the algorithm details
   */
  class MaxMinPicker : public DistPicker {
  public:
    /*! \brief Default Constructor
     *
     */
    MaxMinPicker() {};

    /*! \brief Contains the implementation for the MaxMin diversity picker
     *
     * Here is how the picking algorithm works, refer to
     *   Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604
     * for more detail:
     *
     * A subset of k items is to be selected from a pool containing N molecules. 
     * Then the MaxMin method is as follows:
     *  1. Initialise Subset with a randomly chosen seed
     *     compound and set x = 1.
     *  2. For each of the N-x remaining compounds in Dataset
     *     calculate its dissimilarity with each of the x compounds in
     *     Subset and retain the smallest of these x dissimilarities for
     *     each compound in Dataset.
     *  3. Select the molecule from Dataset with the largest value
     *     for the smallest dissimilarity calculated in Step 2 and
     *     transfer it to Subset.
     *  4. Set x = x + 1 and return to Step 2 if x < k.
     *
     *   \param distMat - distance matrix - a vector of double. It is assumed that only the 
     *              lower triangle element of the matrix are supplied in a 1D array\n
     *   \param poolSize - the size of teh pool to pick the items from. It is assumed that the
     *              distance matrix above contains the right number of elements; i.e.
     *              poolSize*(poolSize-1) \n
     *   \param pickSize - the number items to pick from pool (<= poolSize)
     *
     *   \return a vector with indices of the picked items.
    */
    RDKit::INT_VECT pick(const double *distMat, unsigned int poolSize, unsigned int pickSize) const;
  };
};

#endif
