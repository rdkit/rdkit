//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_DISTPICKER_H
#define _RD_DISTPICKER_H

#include <RDGeneral/types.h>

namespace RDPickers {

  /*! \brief function to lookup distance from 1D lower triangular distance matrix
   *
   *
   *    \param distMat - a pointer to a 1D lower triangular distance matrix \n
   *    \param i - row index \n
   *    \param j - column index \n
   *
   *  RETURNS:
   *
   *    if (i == j) : 0.0
   *    if (i > j) : distMat[i*(i-1)/2 + j]
   *    if (j < i) : distMat[j*(j-1)/2 + i]
   */
  double getDistFromLTM(const double *distMat, unsigned int i, unsigned int j);
      
  /*! \brief Abstract base class to do perform item picking (typically molecules) using a 
   *         distance matrix
   *
   *  This class should never be instantiated by itself. One of the child classes need to be
   *  used. The picking algorithm itself is missing here and only the child calsses implement that
   *  This class contains a pointer to a distance matrix, but it is not responsible for cleaning it up
   */
  class DistPicker {
    
  public:
    /*! \brief Default constructor
     *
     */
    DistPicker(){};
    virtual ~DistPicker() {};
    
    /*! \brief this is a virtual function specific to the type of algorihtm used
     *
     *  The child classes need to implement this function
     *
     *  ARGUMENTS:
     *
     *    \param distMat - distance matrix - a vector of double. It is assumed that only the 
     *              lower triangle elements of the matrix are supplied in a 1D array
     *    \param poolSize - the size of teh pool to pick the items from. It is assumed that the
     *              distance matrix above contains the right number of elements; i.e.
     *              poolSize*(poolSize-1)
     *    \param pickSize - the number items to pick from pool (<= poolSize)
     * 
     *    \return a vector with indices of the picked items.
     */
    virtual RDKit::INT_VECT pick(const double *distMat, unsigned int poolSize,
                                 unsigned int pickSize) const = 0;
  };
};

#endif
