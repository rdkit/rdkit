//
//  Copyright (C) 2003-2007 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_MAXMINPICKER_H__
#define __RD_MAXMINPICKER_H__

#include <RDGeneral/types.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDBoost/Exceptions.h>
#include <cstdlib>
#include "DistPicker.h"
#include <boost/random.hpp>

namespace RDPickers {

  namespace {
    class distmatFunctor{
    public:
      distmatFunctor(const double *distMat) : dp_distMat(distMat) {};
      double operator()(unsigned int i,unsigned int j) {
        return getDistFromLTM(this->dp_distMat,i,j);
      }
    private:
      const double *dp_distMat;
    };
  }

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

    /*! \brief Contains the implementation for a lazy MaxMin diversity picker
     *
     * See the documentation for the pick() method for details about the algorithm
     *
     *   \param func - a function (or functor) taking two unsigned ints as arguments
     *              and returning the distance (as a double) between those two elements.   
     *   \param poolSize - the size of the pool to pick the items from. It is assumed that the
     *              distance matrix above contains the right number of elements; i.e.
     *              poolSize*(poolSize-1) 
     *   \param pickSize - the number items to pick from pool (<= poolSize)
     *   \param firstPicks - (optional)the first items in the pick list
     *   \param seed - (optional) seed for the random number generator
     */
    template <typename T>
    RDKit::INT_VECT lazyPick(T &func, 
                             unsigned int poolSize, unsigned int pickSize,
                             RDKit::INT_VECT firstPicks=RDKit::INT_VECT(),
                             int seed=-1) const;

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
     *   \param distMat - distance matrix - a vector of double. It is assumed that only the 
     *              lower triangle element of the matrix are supplied in a 1D array\n
     *   \param poolSize - the size of the pool to pick the items from. It is assumed that the
     *              distance matrix above contains the right number of elements; i.e.
     *              poolSize*(poolSize-1) \n
     *   \param pickSize - the number items to pick from pool (<= poolSize)
     *   \param firstPicks - indices of the items used to seed the pick set.
     *   \param seed - (optional) seed for the random number generator
    */
    RDKit::INT_VECT pick(const double *distMat, 
                         unsigned int poolSize, unsigned int pickSize,
                         RDKit::INT_VECT firstPicks,
                         int seed=-1) const {
      CHECK_INVARIANT(distMat, "Invalid Distance Matrix");
      if(poolSize<pickSize)
	throw ValueErrorException("pickSize cannot be larger than the poolSize");
      distmatFunctor functor(distMat);
      return this->lazyPick(functor,poolSize,pickSize,firstPicks,seed);
    }

    /*! \overload */
    RDKit::INT_VECT pick(const double *distMat, 
                         unsigned int poolSize, unsigned int pickSize) const {
      RDKit::INT_VECT iv;
      return pick(distMat,poolSize,pickSize,iv);
    }



  };
  // we implement this here in order to allow arbitrary functors without link errors
  template <typename T>
  RDKit::INT_VECT MaxMinPicker::lazyPick(T &func,
                                         unsigned int poolSize, unsigned int pickSize,
                                         RDKit::INT_VECT firstPicks,
                                         int seed) const {
    if(poolSize<pickSize)
      throw ValueErrorException("pickSize cannot be larger than the poolSize");

    RDKit::INT_LIST pool;

    RDKit::INT_VECT picks;
    picks.reserve(pickSize);
    unsigned int pick=0;

    // enter the pool into a list so that we can pick out of it easily
    for (unsigned int i = 0; i < poolSize; i++) {
      pool.push_back(i);
    }


    // get a seeded random number generator:
    typedef boost::mt19937 rng_type;
    typedef boost::uniform_int<> distrib_type;
    typedef boost::variate_generator<rng_type &,distrib_type> source_type;
    rng_type generator(42u);
    distrib_type dist(0,poolSize);
    source_type randomSource(generator,dist);
    if(seed>0) generator.seed(static_cast<rng_type::result_type>(seed));

    // pick the first entry
    if(!firstPicks.size()){
      pick = randomSource();
      // add the pick to the picks
      picks.push_back(pick);
      // and remove it from the pool
      pool.remove(pick);
    } else{
      for(RDKit::INT_VECT::const_iterator pIdx=firstPicks.begin();
          pIdx!=firstPicks.end();++pIdx){
        pick = static_cast<unsigned int>(*pIdx);
	if(pick>=poolSize)
	  throw ValueErrorException("pick index was larger than the poolSize");
        picks.push_back(pick);
        pool.remove(pick);
      }
    }
    // now pick 1 compound at a time
    while (picks.size() < pickSize) {
      double maxOFmin = -1.0;
      RDKit::INT_LIST_I plri=pool.end();
      for(RDKit::INT_LIST_I pli=pool.begin();
          pli!=pool.end(); ++pli){
        unsigned int poolIdx = (*pli);
        double minTOi = RDKit::MAX_DOUBLE;
        for (RDKit::INT_VECT_CI pi = picks.begin(); 
             pi != picks.end(); ++pi) {
          unsigned int pickIdx = (*pi);
          CHECK_INVARIANT(poolIdx!=pickIdx,"");
          double dist = func(poolIdx,pickIdx);
          if (dist <= minTOi) {
            minTOi = dist;
          }
        }
        if (minTOi > maxOFmin || (RDKit::feq(minTOi,maxOFmin) && poolIdx<pick) ) {
          maxOFmin = minTOi;
          pick = poolIdx;
          plri = pli;
        }
      }
      
      // now add the new pick to  picks and remove it from the pool
      picks.push_back(pick);
      CHECK_INVARIANT(plri!=pool.end(),"");
      pool.erase(plri);
    }
    return picks;
  }


};

#endif
