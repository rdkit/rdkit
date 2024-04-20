// $Id$
//
//  Copyright (C) 2004-2012 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "DiscreteValueVect.h"
#include <RDGeneral/Invariant.h>
#include <RDGeneral/StreamOps.h>
#include "DatastructsException.h"
#include "DiscreteDistMat.h"
#include <RDGeneral/Exceptions.h>
#include <cstdint>
#include <algorithm>

namespace RDKit {
const int ci_DISCRETEVALUEVECTPICKLE_VERSION = 0x1;

DiscreteValueVect::DiscreteValueVect(const DiscreteValueVect &other) {
  d_type = other.getValueType();
  d_bitsPerVal = other.getNumBitsPerVal();
  d_numInts = other.getNumInts();
  d_length = other.getLength();
  d_valsPerInt = other.d_valsPerInt;
  d_mask = other.d_mask;
  const std::uint32_t *odata = other.getData();
  auto *data = new std::uint32_t[d_numInts];
  memcpy(static_cast<void *>(data), static_cast<const void *>(odata),
         d_numInts * sizeof(std::uint32_t));
  d_data.reset(data);
}

DiscreteValueVect &DiscreteValueVect::operator=(
    const DiscreteValueVect &other) {
  if (this == &other) {
    return *this;
  }

  d_type = other.getValueType();
  d_bitsPerVal = other.getNumBitsPerVal();
  d_numInts = other.getNumInts();
  d_length = other.getLength();
  d_valsPerInt = other.d_valsPerInt;
  d_mask = other.d_mask;
  const std::uint32_t *odata = other.getData();
  auto *data = new std::uint32_t[d_numInts];
  memcpy(static_cast<void *>(data), static_cast<const void *>(odata),
         d_numInts * sizeof(std::uint32_t));
  d_data.reset(data);

  return *this;
}

unsigned int DiscreteValueVect::getVal(unsigned int i) const {
  if (i >= d_length) {
    throw IndexErrorException(i);
  }
  unsigned int shift = d_bitsPerVal * (i % d_valsPerInt);
  unsigned int intId = i / d_valsPerInt;
  return ((d_data[intId] >> shift) & d_mask);
}

void DiscreteValueVect::setVal(unsigned int i, unsigned int val) {
  if (i >= d_length) {
    throw IndexErrorException(i);
  }
  if ((val & d_mask) != val) {
    throw ValueErrorException("Value out of range");
  }
  unsigned int shift = d_bitsPerVal * (i % d_valsPerInt);
  unsigned int intId = i / d_valsPerInt;
  unsigned int mask = ((1 << d_bitsPerVal) - 1) << shift;
  mask = ~mask;
  d_data[intId] = (d_data[intId] & mask) | (val << shift);
}

unsigned int DiscreteValueVect::getTotalVal() const {
  unsigned int i, j, res = 0;

  for (i = 0; i < d_numInts; ++i) {
    for (j = 0; j < d_valsPerInt; ++j) {
      res += ((d_data[i] >> (j * d_bitsPerVal)) & d_mask);
    }
  }
  return res;
}

unsigned int DiscreteValueVect::getLength() const { return d_length; }

const std::uint32_t *DiscreteValueVect::getData() const { return d_data.get(); }

unsigned int computeL1Norm(const DiscreteValueVect &v1,
                           const DiscreteValueVect &v2) {
  if (v1.getLength() != v2.getLength()) {
    throw ValueErrorException("Comparing vectors of different lengths");
  }

  DiscreteValueVect::DiscreteValueType valType = v1.getValueType();

  if (valType != v2.getValueType()) {
    throw ValueErrorException("Comparing vector of different value types");
  }

  const std::uint32_t *data1 = v1.getData();
  const std::uint32_t *data2 = v2.getData();

  unsigned int res = 0;
  if (valType <= DiscreteValueVect::EIGHTBITVALUE) {
    DiscreteDistMat *dmat = getDiscreteDistMat();

    auto *cd1 = (unsigned char *)(data1);
    auto *cd2 = (unsigned char *)(data2);
    const unsigned char *cend = cd1 + (v1.getNumInts() * 4);
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

    const unsigned short int *send = sd1 + (v1.getNumInts() * 2);
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

std::string DiscreteValueVect::toString() const {
  std::stringstream ss(std::ios_base::binary | std::ios_base::out |
                       std::ios_base::in);

  std::int32_t tVers = ci_DISCRETEVALUEVECTPICKLE_VERSION * -1;
  streamWrite(ss, tVers);
  std::uint32_t tInt;
  tInt = d_type;
  streamWrite(ss, tInt);
  tInt = d_bitsPerVal;
  streamWrite(ss, tInt);
  tInt = d_mask;
  streamWrite(ss, tInt);
  tInt = d_length;
  streamWrite(ss, tInt);
  tInt = d_numInts;
  streamWrite(ss, tInt);

#if defined(BOOST_BIG_ENDIAN)
  std::uint32_t *td = new std::uint32_t[d_numInts];
  for (unsigned int i = 0; i < d_numInts; ++i)
    td[i] = EndianSwapBytes<HOST_ENDIAN_ORDER, LITTLE_ENDIAN_ORDER>(
        d_data.get()[i]);
  ss.write((const char *)td, d_numInts * sizeof(tInt));
  delete[] td;
#else
  ss.write((const char *)d_data.get(), d_numInts * sizeof(tInt));
#endif
  std::string res(ss.str());
  return res;
};

void DiscreteValueVect::initFromText(const char *pkl, const unsigned int len) {
  std::stringstream ss(std::ios_base::binary | std::ios_base::in |
                       std::ios_base::out);
  ss.write(pkl, len);
  std::int32_t tVers;
  streamRead(ss, tVers);
  tVers *= -1;
  if (tVers == 0x1) {
  } else {
    throw ValueErrorException("bad version in DiscreteValueVect pickle");
  }
  std::uint32_t tInt;
  streamRead(ss, tInt);
  d_type = static_cast<DiscreteValueType>(tInt);

  streamRead(ss, tInt);
  d_bitsPerVal = tInt;
  d_valsPerInt = BITS_PER_INT / d_bitsPerVal;
  streamRead(ss, tInt);
  d_mask = tInt;
  streamRead(ss, tInt);
  d_length = tInt;
  streamRead(ss, tInt);
  d_numInts = tInt;
  auto *data = new std::uint32_t[d_numInts];
  ss.read((char *)data, d_numInts * sizeof(std::uint32_t));

#if defined(BOOST_BIG_ENDIAN)
  std::uint32_t *td = new std::uint32_t[d_numInts];
  for (unsigned int i = 0; i < d_numInts; ++i)
    td[i] = EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(data[i]);
  d_data.reset(td);
  delete[] data;
#else
  d_data.reset(data);
#endif
};

DiscreteValueVect DiscreteValueVect::operator&(
    const DiscreteValueVect &other) const {
  PRECONDITION(other.d_length == d_length, "length mismatch");
  DiscreteValueType typ = d_type;
  if (other.d_type < typ) {
    typ = other.d_type;
  }
  DiscreteValueVect ans(typ, d_length);
  for (unsigned int i = 0; i < d_length; ++i) {
    unsigned int v1 = getVal(i);
    unsigned int v2 = other.getVal(i);
    ans.setVal(i, std::min(v2, v1));
  }
  return (ans);
};

DiscreteValueVect DiscreteValueVect::operator|(
    const DiscreteValueVect &other) const {
  PRECONDITION(other.d_length == d_length, "length mismatch");
  DiscreteValueType typ = d_type;
  if (other.d_type > typ) {
    typ = other.d_type;
  }
  DiscreteValueVect ans(typ, d_length);
  for (unsigned int i = 0; i < d_length; ++i) {
    unsigned int v1 = getVal(i);
    unsigned int v2 = other.getVal(i);
    ans.setVal(i, std::max(v2, v1));
  }
  return (ans);
};

DiscreteValueVect &DiscreteValueVect::operator+=(
    const DiscreteValueVect &other) {
  PRECONDITION(other.d_length == d_length, "length mismatch");
  unsigned int maxVal = (1 << d_bitsPerVal) - 1;

  for (unsigned int i = 0; i < d_length; i++) {
    unsigned int v = getVal(i) + other.getVal(i);
    setVal(i, std::min(v, maxVal));
  }
  return *this;
}
DiscreteValueVect &DiscreteValueVect::operator-=(
    const DiscreteValueVect &other) {
  PRECONDITION(other.d_length == d_length, "length mismatch");

  for (unsigned int i = 0; i < d_length; i++) {
    unsigned int v1 = getVal(i);
    unsigned int v2 = other.getVal(i);
    setVal(i, v1 > v2 ? (v1 - v2) : 0);
  }
  return *this;
}

#if 0
  DiscreteValueVect DiscreteValueVect::operator~() const {
    DiscreteValueVect ans(d_type,d_length);
    unsigned int maxVal = (1<<d_bitsPerVal) - 1;
    for(unsigned int i=0;i<d_length;++i){
      unsigned int v1=getVal(i);
      ans.setVal(i,maxVal-v1);
    }
    return(ans);
  };
#endif

DiscreteValueVect operator+(const DiscreteValueVect &p1,
                            const DiscreteValueVect &p2) {
  DiscreteValueVect res(p1);
  res += p2;
  return res;
};
DiscreteValueVect operator-(const DiscreteValueVect &p1,
                            const DiscreteValueVect &p2) {
  DiscreteValueVect res(p1);
  res -= p2;
  return res;
};

}  // end of namespace RDKit
