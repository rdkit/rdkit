//
// Copyright (c) 2003-2008 greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_SPARSEBITVECTS_H__
#define __RD_SPARSEBITVECTS_H__

#include "BitVect.h"

#include <set>
using std::set;
#include <iterator>
#include <algorithm>



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
class SparseBitVect : public BitVect{
public:
  SparseBitVect() : dp_bits(0), d_size(0) {};
  //! initialize with a particular size;
  explicit SparseBitVect(unsigned int size): dp_bits(0), d_size(0) {_initForSize(size); };

  //! copy constructor
  SparseBitVect(const SparseBitVect& other){
    d_size=0;dp_bits = 0;
    _initForSize(other.getNumBits());
    IntSet *bv=other.dp_bits;
    std::copy(bv->begin(),bv->end(),std::inserter(*dp_bits,dp_bits->end()));
  }
  //! construct from a string pickle
  SparseBitVect(const std::string &);
  //! construct from a text pickle
  SparseBitVect(const char *data,const unsigned int dataLen);

  SparseBitVect& operator=(const SparseBitVect&);
  ~SparseBitVect(){ delete dp_bits; };

  bool operator[](const unsigned int which) const;
  SparseBitVect operator| (const SparseBitVect&) const;
  SparseBitVect operator& (const SparseBitVect&) const;
  SparseBitVect operator^ (const SparseBitVect&) const;
  SparseBitVect operator~ () const;

  //! returns a (const) pointer to our raw storage
  const IntSet *getBitSet() const { return dp_bits;}

  unsigned int getNumBits() const { return d_size; };
  bool setBit(const unsigned int which);
  bool setBit(const IntSetIter which);
  bool unsetBit(const unsigned int which);
  bool getBit (const unsigned int which) const;
  bool getBit(const IntVectIter which) const;
  bool getBit(const IntSetIter which) const;

  unsigned int getNumOnBits() const { return dp_bits->size(); };
  unsigned int getNumOffBits() const { return d_size - dp_bits->size(); };

  std::string toString() const;

  void getOnBits (IntVect& v) const;
  void clearBits() { dp_bits->clear(); };
  IntSet *dp_bits; //!< our raw data, exposed for the sake of efficiency

  bool operator==(const SparseBitVect &o) const {
    return *dp_bits==*o.dp_bits;
  }
  bool operator!=(const SparseBitVect &o) const {
    return *dp_bits!=*o.dp_bits;
  }


private:
  unsigned int d_size;
  void _initForSize(const unsigned int size);
};

#endif
