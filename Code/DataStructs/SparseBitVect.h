//
//  Copyright (C) 2007-2024 Greg Landrum and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef __RD_SPARSEBITVECTS_H__
#define __RD_SPARSEBITVECTS_H__

#include "BitVect.h"

#include <set>
using std::set;
#include <iterator>
#include <algorithm>
#include <limits>

typedef set<int> IntSet;
typedef IntSet::iterator IntSetIter;
typedef IntSet::const_iterator IntSetConstIter;

//! a class for bit vectors that are sparsely occupied.
/*!
    SparseBitVect objects store only their on bits, in an
    std::set.

    They are, as you might expect, quite memory efficient for sparsely populated
    vectors but become rather a nightmare if they need to be negated.

 */
class RDKIT_DATASTRUCTS_EXPORT SparseBitVect : public BitVect {
 public:
  SparseBitVect() {}
  //! initialize with a particular size;
  explicit SparseBitVect(unsigned int size) : dp_bits(nullptr), d_size(0) {
    _initForSize(size);
  }

  //! copy constructor
  SparseBitVect(const SparseBitVect &other) : BitVect(other) {
    d_size = 0;
    dp_bits = nullptr;
    _initForSize(other.getNumBits());
    IntSet *bv = other.dp_bits;
    std::copy(bv->begin(), bv->end(), std::inserter(*dp_bits, dp_bits->end()));
  }
  //! construct from a string pickle
  SparseBitVect(const std::string &pkl);
  //! construct from a text pickle
  SparseBitVect(const char *data, const unsigned int dataLen);

  SparseBitVect &operator=(const SparseBitVect &);
  ~SparseBitVect() override { delete dp_bits; }

  bool operator[](const unsigned int which) const override;
  SparseBitVect operator|(const SparseBitVect &) const;
  SparseBitVect operator&(const SparseBitVect &) const;
  SparseBitVect operator^(const SparseBitVect &) const;
  SparseBitVect operator~() const;

  //! returns a (const) pointer to our raw storage
  const IntSet *getBitSet() const { return dp_bits; }

  unsigned int getNumBits() const override { return d_size; }
  bool setBit(const unsigned int which) override;
  bool setBit(const IntSetIter which);
  bool unsetBit(const unsigned int which) override;
  bool getBit(const unsigned int which) const override;
  bool getBit(const IntVectIter which) const;
  bool getBit(const IntSetIter which) const;

  unsigned int getNumOnBits() const override {
    return static_cast<unsigned int>(dp_bits->size());
  }
  unsigned int getNumOffBits() const override {
    return d_size - static_cast<unsigned int>(dp_bits->size());
  }

  std::string toString() const override;

  void getOnBits(IntVect &v) const override;
  void clearBits() override { dp_bits->clear(); }
  IntSet *dp_bits{
      nullptr};  //!< our raw data, exposed for the sake of efficiency

  bool operator==(const SparseBitVect &o) const {
    return *dp_bits == *o.dp_bits;
  }
  bool operator!=(const SparseBitVect &o) const {
    return *dp_bits != *o.dp_bits;
  }

 private:
  unsigned int d_size{0};
  void _initForSize(const unsigned int size) override;
  bool checkIndex(const unsigned int idx) const {
    return idx < d_size || (idx == d_size &&
                            d_size == std::numeric_limits<unsigned int>::max());
  }
  template <typename T>
  bool checkIndex(const T which) const {
    return *which >= 0 && static_cast<unsigned int>(*which) < d_size;
  }
};

#endif
