//
// Copyright (c) 2003-208 greg Landrum and Rational Discovery LLC
//  Copyright (c) 2014, Novartis Institutes for BioMedical Research Inc.
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
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
  explicit ExplicitBitVect(unsigned int size) : dp_bits(0), d_size(0), d_numOnBits(0) {_initForSize(size);};
  //! initialize with a particular size and all bits set
  ExplicitBitVect(unsigned int size, bool bitsSet);
  ExplicitBitVect(const ExplicitBitVect& other);
  //! construct from a string pickle
  ExplicitBitVect(const std::string &);
  //! construct from a text pickle
  ExplicitBitVect(const char *,const unsigned int);
  
  ~ExplicitBitVect();

  ExplicitBitVect& operator=(const ExplicitBitVect& other);
  bool operator[] (const unsigned int which) const;
  bool setBit(const unsigned int which);
  bool unsetBit(const unsigned int which);
  bool getBit(const unsigned int which) const;

  ExplicitBitVect operator^ (const ExplicitBitVect &other) const;
  ExplicitBitVect operator& (const ExplicitBitVect &other) const;
  ExplicitBitVect operator| (const ExplicitBitVect &other) const;
  ExplicitBitVect operator~ () const;
  /* concatenate two ExplicitBitVects */
  ExplicitBitVect operator+ (const ExplicitBitVect &other) const;

  ExplicitBitVect& operator^= (const ExplicitBitVect &other);
  ExplicitBitVect& operator&= (const ExplicitBitVect &other);
  ExplicitBitVect& operator|= (const ExplicitBitVect &other);
  /* concatenate two ExplicitBitVects */
  ExplicitBitVect& operator+= (const ExplicitBitVect &other);

  unsigned int getNumBits() const;
  unsigned int getNumOnBits() const;
  unsigned int getNumOffBits() const;

  void getOnBits (IntVect& v) const;

  void clearBits() { dp_bits->reset(); };
  std::string toString() const;
  
  boost::dynamic_bitset<> *dp_bits; //!< our raw storage

  bool operator==(const ExplicitBitVect &o) const {
    return *dp_bits==*o.dp_bits;
  }
  bool operator!=(const ExplicitBitVect &o) const {
    return *dp_bits!=*o.dp_bits;
  }

private:
  unsigned int d_size;
  unsigned int d_numOnBits;
  void _initForSize(const unsigned int size);
};


#endif
