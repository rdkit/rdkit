//
// Copyright (c) 2003-2006 greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
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
  explicit SparseBitVect(unsigned int size): dp_bits(0), d_size(0) {_InitForSize(size); };

  //! copy constructor
  SparseBitVect(const SparseBitVect& other){
    d_size=0;dp_bits = 0;
    _InitForSize(other.GetNumBits());
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
  const IntSet *GetBitSet() const { return dp_bits;}

  const unsigned int GetNumBits() const { return d_size; };
  bool SetBit(const unsigned int which);
  bool SetBit(const IntSetIter which);
  bool UnSetBit(const unsigned int which);
  bool GetBit (const unsigned int which) const;
  bool GetBit(const IntVectIter which) const;
  bool GetBit(const IntSetIter which) const;

  const unsigned int GetNumOnBits() const { return dp_bits->size(); };
  const unsigned int GetNumOffBits() const { return d_size - dp_bits->size(); };

  std::string ToString() const;

  void GetOnBits (IntVect& v) const;
  void ClearBits() { dp_bits->clear(); };
  IntSet *dp_bits; //!< our raw data, exposed for the sake of efficiency
private:
  unsigned int d_size;
  void _InitForSize(const unsigned int size);
};

#endif
