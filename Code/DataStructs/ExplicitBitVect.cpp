// $Id$
//
// Copyright (c) 2001-2008 greg Landrum and Rational Discovery LLC
//  Copyright (c) 2014, Novartis Institutes for BioMedical Research Inc.
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <iostream>
#include <RDBoost/Exceptions.h>
#include "ExplicitBitVect.h"
#include <RDGeneral/StreamOps.h>
#include "base64.h"
#include <sstream>
#include <limits>
#ifdef WIN32
#include <ios>
#endif
#include <boost/cstdint.hpp>

ExplicitBitVect::ExplicitBitVect(unsigned int size, bool bitsSet)
{
  d_size=0;dp_bits = 0;d_numOnBits=0;
  _initForSize(size);
  if (bitsSet) {
    dp_bits->set(); // set all bits to 1
    d_numOnBits = size;
  }
}
ExplicitBitVect::ExplicitBitVect(const std::string &s)
{
  d_size=0;dp_bits = 0;d_numOnBits=0;
  initFromText(s.c_str(),s.length());
}
ExplicitBitVect::ExplicitBitVect(const char *data,const unsigned int dataLen)
{
  d_size=0;dp_bits = 0;d_numOnBits=0;
  initFromText(data,dataLen);
}

  ExplicitBitVect::ExplicitBitVect(const ExplicitBitVect& other){
    d_size = other.d_size;
    dp_bits = new boost::dynamic_bitset<>(*(other.dp_bits));
    d_numOnBits=other.d_numOnBits;
  };

  ExplicitBitVect& ExplicitBitVect::operator=(const ExplicitBitVect& other){
    d_size = other.d_size;
    delete dp_bits;
    dp_bits = new boost::dynamic_bitset<>(*(other.dp_bits));
    d_numOnBits=other.d_numOnBits;
    return *this;
  };
  bool ExplicitBitVect::operator[] (const unsigned int which) const {
    if(which >= d_size){
      throw IndexErrorException(which);
    }
    return (bool)(*dp_bits)[which];
  };
  bool ExplicitBitVect::setBit(const unsigned int which){
    if(which >= d_size){
      throw IndexErrorException(which);
    }
    if((bool)(*dp_bits)[which]){
      return true;
    } else {
      (*dp_bits)[which] = 1;
      ++d_numOnBits;
      return false;
    }
  };
  bool ExplicitBitVect::unsetBit(const unsigned int which){
    if(which >= d_size){
      throw IndexErrorException(which);
    }
    if((bool)(*dp_bits)[which]){
      (*dp_bits)[which] = 0;
      --d_numOnBits;
      return true;
    } else {
      return false;
    }
  };
  bool ExplicitBitVect::getBit(const unsigned int which) const {
    if(which >= d_size){
      throw IndexErrorException(which);
    }
    return((bool)(*dp_bits)[which]);
  };

  ExplicitBitVect ExplicitBitVect::operator^ (const ExplicitBitVect &other) const {
    ExplicitBitVect ans(d_size);
    *(ans.dp_bits) = (*dp_bits) ^ *(other.dp_bits);
    ans.d_numOnBits=ans.dp_bits->count();
    return(ans);
  };

  ExplicitBitVect ExplicitBitVect::operator& (const ExplicitBitVect &other) const {
    ExplicitBitVect ans(d_size);
    *(ans.dp_bits) = (*dp_bits) & *(other.dp_bits);
    ans.d_numOnBits=ans.dp_bits->count();
    return(ans);
  };

  ExplicitBitVect ExplicitBitVect::operator| (const ExplicitBitVect &other) const {
    ExplicitBitVect ans(d_size);
    *(ans.dp_bits) = (*dp_bits) | *(other.dp_bits);
    ans.d_numOnBits=ans.dp_bits->count();
    return(ans);
  };
  
  ExplicitBitVect& ExplicitBitVect::operator^= (const ExplicitBitVect &other) {
    *(dp_bits) ^= *(other.dp_bits);
    d_numOnBits=dp_bits->count();
    return *this;
  };

  ExplicitBitVect& ExplicitBitVect::operator&= (const ExplicitBitVect &other) {
    *(dp_bits) &= *(other.dp_bits);
    d_numOnBits=dp_bits->count();
    return *this;
  };

  ExplicitBitVect& ExplicitBitVect::operator|= (const ExplicitBitVect &other) {
    *(dp_bits) |= *(other.dp_bits);
    d_numOnBits=dp_bits->count();
    return *this;
  };

  ExplicitBitVect ExplicitBitVect::operator~ () const {
    ExplicitBitVect ans(d_size);
    *(ans.dp_bits) = ~(*dp_bits);
    ans.d_numOnBits=ans.dp_bits->count();
    return(ans);
  };

  ExplicitBitVect& ExplicitBitVect::operator+= (const ExplicitBitVect &other) {
    dp_bits->resize(d_size+other.d_size);
    unsigned int original_size = d_size;
    d_size = dp_bits->size();
    for(unsigned i=0;i<other.d_size;i++){
      if(other[i]){
        setBit(i+original_size);
      }
    }
    d_numOnBits=dp_bits->count();
    return *this;
  };

  ExplicitBitVect ExplicitBitVect::operator+ (const ExplicitBitVect &other) const {
    ExplicitBitVect ans(*this);
    return ans+=other;
  };

  unsigned int ExplicitBitVect::getNumBits() const {
    return d_size;
  };
  unsigned int ExplicitBitVect::getNumOnBits() const {
    return d_numOnBits;
  };
  unsigned int ExplicitBitVect::getNumOffBits() const {
    return d_size - d_numOnBits;
  };

  // the contents of v are blown out
  void ExplicitBitVect::getOnBits (IntVect& v) const {
    unsigned int nOn = getNumOnBits();
    if(!v.empty()) IntVect().swap(v);
    v.reserve(nOn);
    for(unsigned int i=0;i<d_size;i++){
      if((bool)(*dp_bits)[i]) v.push_back(i);
    }
  };

  void ExplicitBitVect::_initForSize(unsigned int size) {
    d_size = size;
    delete dp_bits;
    dp_bits = new boost::dynamic_bitset<>(size);
    d_numOnBits=0;
  };

  ExplicitBitVect::~ExplicitBitVect() {
    delete dp_bits;
    dp_bits=NULL;
  };

std::string
ExplicitBitVect::toString() const
{
  // This Function replaces the older version (version 16) of writing the onbits to
  // a string
  // the old version does not perform any run length encoding, it only checks to see if 
  // the length of the bitvect can be short ints and writes the on bits as shorts 
  // other wise the onbits are all written as ints

  // here we do run length encoding and the version number has been bumped to 32 as well. 
  // only the reader needs to take care of readinf all legacy versions
  // also in this scheme each bit number written to the string is checked to see how many 
  // bytes it needs
  std::stringstream ss(std::ios_base::binary|std::ios_base::out|std::ios_base::in);

  boost::int32_t tInt = ci_BITVECT_VERSION*-1;
  RDKit::streamWrite(ss,tInt);
  tInt=d_size;
  RDKit::streamWrite(ss,tInt);
  tInt=getNumOnBits();
  RDKit::streamWrite(ss,tInt);

  int prev = -1;
  unsigned int zeroes;
  for(unsigned int i=0;i<d_size;i++){
    if( (bool)(*dp_bits)[i] ){
      zeroes = i - prev -1;
      RDKit::appendPackedIntToStream(ss, zeroes);
      prev = i;
    }
  }
  zeroes = d_size - prev -1;
  RDKit::appendPackedIntToStream(ss, zeroes);
  std::string res(ss.str());
  return res;
}


