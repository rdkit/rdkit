//
// Copyright (c) 2001-2024 greg Landrum and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "SparseBitVect.h"
#include <RDGeneral/Exceptions.h>

#include "base64.h"
#include <RDGeneral/StreamOps.h>
#include <sstream>

#ifdef WIN32
#include <ios>
#endif
#include <cstdint>

// """ -------------------------------------------------------
//
// Construct a SparseBitVect from a binary string.
//   The format of the string should be
//   the format produced by SparseBitVect::toString
//
// """ -------------------------------------------------------
SparseBitVect::SparseBitVect(const std::string &s) {
  d_size = 0;
  dp_bits = nullptr;
  initFromText(s.c_str(), s.length());
}

// """ -------------------------------------------------------
//
// Construct a SparseBitVect from a binary string stored as a char *.
//   The format of the string should be
//   the format produced by SparseBitVect::toString
//
// """ -------------------------------------------------------
SparseBitVect::SparseBitVect(const char *data, const unsigned int dataLen) {
  d_size = 0;
  dp_bits = nullptr;
  initFromText(data, dataLen);
}

// """ -------------------------------------------------------
//
//  Operator[]
//   Returns the state of the ith bit.
//
//   **C++ NOTE** The state of BitVects cannot be changed using
//    operator&.  So this will not work:
//     SBV[i] = 0;
//    In Python this type of assignment **is** valid.
//
// """ -------------------------------------------------------
bool SparseBitVect::operator[](const unsigned int which) const {
  if (which >= d_size) {
    throw IndexErrorException(which);
  }
  return dp_bits->count(which) > 0u;
}

// """ -------------------------------------------------------
//
//  Assignment operator
//    The bits of the other SBV are copied.
//
// """ -------------------------------------------------------
SparseBitVect &SparseBitVect::operator=(const SparseBitVect &other) {
  if (this == &other) {
    return *this;
  }
  IntSet *bv = other.dp_bits;
  delete dp_bits;
  d_size = other.getNumBits();
  dp_bits = new IntSet;
  std::copy(bv->begin(), bv->end(), std::inserter(*dp_bits, dp_bits->end()));

  return *this;
}

// -------------------------------------------------------
//
//  Operator|
//   allows SBV3 = SBV1|SBV2;
//
// -------------------------------------------------------
SparseBitVect SparseBitVect::operator|(const SparseBitVect &other) const {
  SparseBitVect ans(d_size);
  std::set_union(dp_bits->begin(), dp_bits->end(), other.dp_bits->begin(),
                 other.dp_bits->end(),
                 std::inserter(*(ans.dp_bits), ans.dp_bits->end()));
  return ans;
}

// """ -------------------------------------------------------
//
//  Operator&
//   allows SBV3 = SBV1&SBV2;
//
// """ -------------------------------------------------------
SparseBitVect SparseBitVect::operator&(const SparseBitVect &other) const {
  SparseBitVect ans(d_size);
  std::set_intersection(dp_bits->begin(), dp_bits->end(),
                        other.dp_bits->begin(), other.dp_bits->end(),
                        std::inserter(*(ans.dp_bits), ans.dp_bits->end()));
  return ans;
}

// """ -------------------------------------------------------
//
//  Operator^
//   allows SBV3 = SBV1^SBV2;
//
// """ -------------------------------------------------------
SparseBitVect SparseBitVect::operator^(const SparseBitVect &other) const {
  SparseBitVect ans(d_size);
  std::set_symmetric_difference(
      dp_bits->begin(), dp_bits->end(), other.dp_bits->begin(),
      other.dp_bits->end(), std::inserter(*(ans.dp_bits), ans.dp_bits->end()));
  return (ans);
}

// """ -------------------------------------------------------
//
//  Operator~  (the negation operator)
//   allows SBV2 = ~SBV1;
//
// """ -------------------------------------------------------
SparseBitVect SparseBitVect::operator~() const {
  SparseBitVect ans(d_size);
  for (unsigned int i = 0; i < d_size; i++) {
    if (!getBit(i)) {
      ans.setBit(i);
    }
  }

  return (ans);
}

// """ -------------------------------------------------------
//
// getBit(int which)
//   Returns the state of bit which
//
// """ -------------------------------------------------------
bool SparseBitVect::getBit(const unsigned int which) const {
  if (!checkIndex(which)) {
    throw IndexErrorException(which);
  }
  return dp_bits->count(which) > 0u;
}

// """ -------------------------------------------------------
//
// getBit(const IntVectIter which)  (C++ SPECIFIC)
//   Returns the state of bit which
//
// """ -------------------------------------------------------
bool SparseBitVect::getBit(const IntVectIter which) const {
  if (!checkIndex(which)) {
    throw IndexErrorException(*which);
  }
  return dp_bits->count(*which) > 0u;
}

// """ -------------------------------------------------------
//
// getBit(const IntSetIter which)  (C++ SPECIFIC)
//   Returns the state of bit which
//
// """ -------------------------------------------------------
bool SparseBitVect::getBit(const IntSetIter which) const {
  if (!checkIndex(which)) {
    throw IndexErrorException(*which);
  }
  return dp_bits->count(*which) > 0u;
}

// """ -------------------------------------------------------
//
//  Sets bit which to be on.
//
//  Returns the original state of the bit
//
// """ -------------------------------------------------------
bool SparseBitVect::setBit(const unsigned int which) {
  if (!dp_bits) {
    throw ValueErrorException("BitVect not properly initialized.");
  }
  if (!checkIndex(which)) {
    throw IndexErrorException(which);
  }
  auto res = dp_bits->insert(which);
  return !(res.second);
}

// """ -------------------------------------------------------
//
// setBit(const IntSetIter which)  (C++ SPECIFIC)
//  Sets bit which to be on.
//
//  Returns the original state of the bit
//
// """ -------------------------------------------------------
bool SparseBitVect::setBit(const IntSetIter which) {
  if (!dp_bits) {
    throw ValueErrorException("BitVect not properly initialized.");
  }
  if (!checkIndex(which)) {
    throw IndexErrorException(*which);
  }
  auto res = dp_bits->insert(*which);
  return !(res.second);
}

// """ -------------------------------------------------------
//
//  Sets bit which to be off.
//
//  Returns the original state of the bit
//
// """ -------------------------------------------------------
bool SparseBitVect::unsetBit(const unsigned int which) {
  if (!dp_bits) {
    throw ValueErrorException("BitVect not properly initialized.");
  }
  if (!checkIndex(which)) {
    throw IndexErrorException(which);
  }

  if (dp_bits->count(which)) {
    dp_bits->erase(dp_bits->find(which));
    return true;
  } else {
    return false;
  }
}

// """ -------------------------------------------------------
//
// getOnBits(IntVect &which)
//  C++: Passes the set of on bits out in the IntVect passed in.
//   The contents of IntVect are destroyed.
//
//  Python: Returns the tuple of on bits
//
// """ -------------------------------------------------------
void SparseBitVect::getOnBits(IntVect &v) const {
  if (!dp_bits) {
    throw ValueErrorException("BitVect not properly initialized.");
  }
  unsigned int nOn = getNumOnBits();
  if (!v.empty()) {
    IntVect().swap(v);
  }
  v.reserve(nOn);
  v.resize(nOn);
  std::copy(dp_bits->begin(), dp_bits->end(), v.begin());
};

// """ -------------------------------------------------------
//
// toString()
//  Returns a binary string with the contents of the BitVector
//
// """ -------------------------------------------------------
std::string SparseBitVect::toString() const {
  // This Function replaces the older version (version 16) of writing the onbits
  // to
  // a string
  // the old version does not perform any run length encoding, it only checks to
  // see if
  // the length of the bitvect can be short ints and writes the on bits as
  // shorts
  // other wise the onbits are all written as ints

  // here we do run length encoding and the version number has been bumped to 32
  // as well.
  // only the reader needs to take care of readinf all legacy versions
  // also in this scheme each bit number written to the string is checked to see
  // how many
  // bytes it needs
  std::stringstream ss(std::ios_base::binary | std::ios_base::out |
                       std::ios_base::in);

  std::int32_t tInt = ci_BITVECT_VERSION * -1;
  RDKit::streamWrite(ss, tInt);
  tInt = d_size;
  RDKit::streamWrite(ss, tInt);
  tInt = getNumOnBits();
  RDKit::streamWrite(ss, tInt);

  int prev = -1;
  unsigned int zeroes;
  for (int dp_bit : *dp_bits) {
    zeroes = dp_bit - prev - 1;
    RDKit::appendPackedIntToStream(ss, zeroes);
    prev = dp_bit;
  }
  zeroes = d_size - prev - 1;
  RDKit::appendPackedIntToStream(ss, zeroes);

  std::string res(ss.str());
  return res;
}

void SparseBitVect::_initForSize(unsigned int size) {
  d_size = size;
  delete dp_bits;
  dp_bits = new IntSet;
};
