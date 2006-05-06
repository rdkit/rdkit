//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
#include "DiscreteValueVect.h"
#include <RDGeneral/Invariant.h>
#include "DatastructsException.h"
#include "DiscreteDistMat.h"
#include <RDBoost/Exceptions.h>

DiscreteValueVect::DiscreteValueVect(const DiscreteValueVect &other) {
  d_type = other.getValueType();
  d_bitsPerVal = other.getNumBitsPerVal();
  d_numInts = other.getNumInts();
  d_length = other.getLength();
  d_valsPerInt = BITS_PER_INT/d_bitsPerVal;
  const unsigned int *odata = other.getData();
  unsigned int *data = new unsigned int[d_numInts];
  memcpy(static_cast<void *>(data), static_cast<const void *>(odata),
         d_numInts*sizeof(unsigned int));
  d_data.reset(data);
}

unsigned int DiscreteValueVect::getVal(unsigned int i) const {
  RANGE_CHECK(0, i, d_length-1);
  unsigned int shift = d_bitsPerVal*(i%d_valsPerInt);
  unsigned int intId = i/d_valsPerInt;
  return ( (d_data[intId] >> shift) & d_mask);
}

void DiscreteValueVect::setVal(unsigned int i, unsigned int val) {
  RANGE_CHECK(0, i, d_length-1);
  if ((val & d_mask) != val) {
    throw ValueErrorException("Value out of range");
  }
  unsigned int shift = d_bitsPerVal*(i%d_valsPerInt);
  unsigned int intId = i/d_valsPerInt;
  unsigned int mask = ((1<<d_bitsPerVal) -1) << shift;
  mask = ~mask;
  d_data[intId] = (d_data[intId]&mask)|(val << shift);
}

unsigned int DiscreteValueVect::getTotalVal() const {
  unsigned int i, j, res = 0;
  
  for (i = 0; i < d_numInts; ++i) {
    for (j = 0; j < d_valsPerInt; ++j) {
      res += ((d_data[i] >> (j*d_bitsPerVal)) & d_mask);
    }
  }
  return res;
}

unsigned int DiscreteValueVect::getLength() const {
  return d_length;
}

const unsigned int *DiscreteValueVect::getData() const {
  return d_data.get();
}
  
unsigned int computeL1Norm(const DiscreteValueVect &v1, const DiscreteValueVect &v2) {
  if (v1.getLength() != v2.getLength()) {
    throw ValueErrorException("Comparing vectors of different lengths");
  }

  DiscreteValueVect::DiscreteValueType valType = v1.getValueType();

  if (valType != v2.getValueType()) {
    throw ValueErrorException("Comparing vector of different value types");
  }

  const unsigned int* data1 = v1.getData();
  const unsigned int* data2 = v2.getData();

  unsigned int res = 0;
  if (valType <= DiscreteValueVect::EIGHTBITVALUE) {
    DiscreteDistMat *dmat = getDiscreteDistMat();
    
    unsigned char *cd1 = (unsigned char *)(data1);
    unsigned char *cd2 = (unsigned char *)(data2);
    const unsigned char *cend = cd1 + (v1.getNumInts()*4);
    while (cd1 != cend) {
      if (*cd1 == *cd2) {
        cd1++;
        cd2++;
        continue;
      }
      res += dmat->getDist(*cd1, *cd2, valType);
      cd1++;
      cd2++;
    }
  } else {
    // we have a sixteen bits per value type
    // REVIEW: we are making an assumption here that a short 
    // is 16 bit - may fail on a different compiler
    const unsigned short int *sd1 = (unsigned short int *)(data1);
    const unsigned short int *sd2 = (unsigned short int *)(data2);
    
    const unsigned short int *send = sd1 + (v1.getNumInts()*2);
    while (sd1 != send) {
      if (*sd1 == *sd2) {
        sd1++;
        sd2++;
        continue;
      }
      res += abs((*sd1) - (*sd2));
      sd1++;
      sd2++;
    }
  }
  return res;
}


