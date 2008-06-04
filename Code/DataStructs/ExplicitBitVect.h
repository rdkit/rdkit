//
// Copyright (c) 2003-2006 greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
#ifndef __RD_EXPLICITBITVECTS_H__
#define __RD_EXPLICITBITVECTS_H__

#include <boost/dynamic_bitset.hpp>
#include "BitVect.h"

//! a class for bit vectors that are densely occupied
/*!
    ExplicitBitVect objects store all of their bits using
    a boost::dynamic_bitset

    These are very fast, but can require large amounts of memory for large,
    sparsely occupied vectors.

 */
class ExplicitBitVect : public BitVect {
public:
  ExplicitBitVect() : dp_bits(0), d_size(0), d_numOnBits(0) {};
  //! initialize with a particular size;
  explicit ExplicitBitVect(unsigned int size) : dp_bits(0), d_size(0), d_numOnBits(0) {_InitForSize(size);};
  ExplicitBitVect(const ExplicitBitVect& other);
  //! construct from a string pickle
  ExplicitBitVect(const std::string &);
  //! construct from a text pickle
  ExplicitBitVect(const char *,const unsigned int);
  
  ~ExplicitBitVect();

  ExplicitBitVect& operator=(const ExplicitBitVect& other);
  bool operator[] (const unsigned int which) const;
  bool SetBit(const unsigned int which);
  bool UnSetBit(const unsigned int which);
  bool GetBit(const unsigned int which) const;

  ExplicitBitVect operator^ (const ExplicitBitVect &other) const;
  ExplicitBitVect operator& (const ExplicitBitVect &other) const;
  ExplicitBitVect operator| (const ExplicitBitVect &other) const;
  ExplicitBitVect operator~ () const;
  const unsigned int GetNumBits() const;
  const unsigned int GetNumOnBits() const;
  const unsigned int GetNumOffBits() const;

  void GetOnBits (IntVect& v) const;

  // FIX: complete these
  void ClearBits() { dp_bits->reset(); };
  std::string ToString() const;
  
  boost::dynamic_bitset<> *dp_bits; //!< our raw storage
private:
  unsigned int d_size;
  unsigned int d_numOnBits;
  void _InitForSize(const unsigned int size);
};


#endif
