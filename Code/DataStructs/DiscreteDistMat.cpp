// $Id$
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "DiscreteDistMat.h"
#include "DiscreteValueVect.h"
#include <iostream>
#include "DatastructsException.h"

namespace RDKit {
  void _fillDistMat(unsigned int dmat[], unsigned int nBits) {
    unsigned int i,j, a, b, ta, tb, dist;
    int temp;
    unsigned int mask = ((1<<nBits) -1);
    for (i = 0; i < 256; ++i) {
      for (j = 0; j < 256; ++j) {
        dist = 0;
        a = i;
        b = j;
        while (a || b) {
          ta = a&mask;
          tb = b&mask;
          temp = ta-tb;
          if (temp > 0) {
            dist += temp;
          } else {
            dist -= temp;
          }
          a >>= nBits;
          b >>= nBits;
        }
        dmat[i*256 + j] = dist;
      }
    }
  }

  DiscreteDistMat::DiscreteDistMat() {
    // fill in the distance matrix table

    // one bit per value table
    _fillDistMat(d_oneBitTab, 1);

    // two bits per value table
    _fillDistMat(d_twoBitTab, 2);

    // four bits per value table
    _fillDistMat(d_fourBitTab, 4);
  }

  unsigned int DiscreteDistMat::getDist(unsigned char v1, 
                                        unsigned char v2, 
                                        DiscreteValueVect::DiscreteValueType type) {
    unsigned int res=0;
    int temp;
    unsigned int id = static_cast<unsigned int>(v1)*256 + static_cast<unsigned int>(v2);
    switch(type) {
    case DiscreteValueVect::ONEBITVALUE :
      res = d_oneBitTab[id];
      break;
    case DiscreteValueVect::TWOBITVALUE :
      res = d_twoBitTab[id];
      break;
    case DiscreteValueVect::FOURBITVALUE :
      res = d_fourBitTab[id];
      break;
    case DiscreteValueVect::EIGHTBITVALUE :
      temp = static_cast<unsigned int>(v1) - static_cast<unsigned int>(v2);
      if (temp < 0) {
        res -= temp;
      } else {
        res += temp;
      }
      break;
    default:
      // ummm.. we shouldn't have come here
      throw DatastructsException("We shouldn't be here");
    }
    return res;
  }

  static DiscreteDistMat discreteDMat;
  DiscreteDistMat *getDiscreteDistMat() {
    return &discreteDMat;
  }

}  
