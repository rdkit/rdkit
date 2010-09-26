//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_BOUNDS_MATRIX_H__
#define __RD_BOUNDS_MATRIX_H__

#include <RDGeneral/Invariant.h>
#include <boost/smart_ptr.hpp>
#include <iostream>
#include <iomanip>
#include <Numerics/SquareMatrix.h>

namespace DistGeom {
  //! Class to store the distance bound
  /*! 
    Basically a N by N matrix 
    with lower distance bounds on the lower traingle and upper bounds in the upper 
    triangle
  */
  class BoundsMatrix : public RDNumeric::SquareMatrix<double> {
  public:
    typedef boost::shared_array<double> DATA_SPTR;

    explicit BoundsMatrix(unsigned int N) : RDNumeric::SquareMatrix<double>(N,0.0) {}; 
    BoundsMatrix(unsigned int N, DATA_SPTR data) : 
      RDNumeric::SquareMatrix<double>(N,data) {}; 

    //! Get the upper bound between points i and j
    inline double getUpperBound(unsigned int i, unsigned int j) const {
      RANGE_CHECK(0, i, d_nRows-1);
      RANGE_CHECK(0, j, d_nCols-1);

      if (i < j) {
        return getVal(i,j);
      } else {
        return getVal(j,i);
      }
    }
    
    //! Set the lower bound between points i and j
    inline void setUpperBound(unsigned int i, unsigned int j, double val) {
      RANGE_CHECK(0, i, d_nRows-1);
      RANGE_CHECK(0, j, d_nCols-1);
      CHECK_INVARIANT(val >= 0.0, "Negative upper bound");
      if (i < j) {
        setVal(i,j,val);
      } else {
        setVal(j,i,val);
      }
    }
    
    //! Set the upper bound between points i and j only if it is better than 
    //! previously existing value (i.e. the new value is smaller)
    inline void setUpperBoundIfBetter(unsigned int i, unsigned int j, double val) {
      if ((val < getUpperBound(i, j)) && (val > getLowerBound(i, j)) ) {
        setUpperBound(i, j, val);
      }
    }

    //! Set the lower bound between points i and j 
    inline void setLowerBound(unsigned int i, unsigned int j, double val) {
      RANGE_CHECK(0, i, d_nRows-1);
      RANGE_CHECK(0, j, d_nCols-1);
      CHECK_INVARIANT(val >= 0.0, "Negative lower bound");
      if (i < j) {
        setVal(j,i,val);
      } else {
        setVal(i,j,val);
      }
    }
    
    //! Set the lower bound between points i and j only if it is better than 
    //! previously existing value (i.e. the new value is larger)
    inline void setLowerBoundIfBetter(unsigned int i, unsigned int j, double val) {
      if ((val > getLowerBound(i,j)) && (val < getUpperBound(i,j))) {
        setLowerBound(i,j, val);
      } 
    } 
    
    //! Get the lower bound between points i and j
    inline double getLowerBound(unsigned int i, unsigned int j) const {
      RANGE_CHECK(0, i, d_nRows-1);
      RANGE_CHECK(0, j, d_nCols-1);

      if (i < j) {
        return getVal(j,i);
      } else {
        return getVal(i,j);
      }
    }

    //! Do a simple check of the current bounds - i.e. all lower bounds are
    //! smaller than the existing upper bounds
    inline bool checkValid() const {
      unsigned int i, j;
      for (i = 1; i < d_nRows; i++) {
        for (j = 0; j < i; j++) {
          if (getUpperBound(i,j) < getLowerBound(i,j)) {
            return false;
          }
        }
      }
      return true;
    }
  }; 

  typedef boost::shared_ptr<BoundsMatrix> BoundsMatPtr;
}

#endif
          
