// $Id: ExplicitBitVect.cpp 5062 2006-03-08 01:43:55Z glandrum $
//
// Copyright (c) 2001-2006 greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
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


ExplicitBitVect::ExplicitBitVect(const std::string s)
{
  d_size=0;dp_bits = 0;
  InitFromText(s.c_str(),s.length());
}
ExplicitBitVect::ExplicitBitVect(const char *data,const unsigned int dataLen)
{
  d_size=0;dp_bits = 0;
  InitFromText(data,dataLen);
}

  ExplicitBitVect::ExplicitBitVect(const ExplicitBitVect& other){
    d_size = other.d_size;
    dp_bits = new boost::dynamic_bitset<>(*(other.dp_bits));
  };

  ExplicitBitVect& ExplicitBitVect::operator=(const ExplicitBitVect& other){
    d_size = other.d_size;
    dp_bits = new boost::dynamic_bitset<>(*(other.dp_bits));
    return *this;
  };
  bool ExplicitBitVect::operator[] (const unsigned int which) const {
    if(which < 0 || which >= d_size){
      throw IndexErrorException(which);
    }
    return (bool)(*dp_bits)[which];
  };
  bool ExplicitBitVect::SetBit(const unsigned int which){
    if(which < 0 || which >= d_size){
      throw IndexErrorException(which);
    }
    if((bool)(*dp_bits)[which]){
      return true;
    } else {
      (*dp_bits)[which] = 1;
      return false;
    }
  };
  bool ExplicitBitVect::UnSetBit(const unsigned int which){
    if(which < 0 || which >= d_size){
      throw IndexErrorException(which);
    }
    if((bool)(*dp_bits)[which]){
      (*dp_bits)[which] = 0;
      return true;
    } else {
      return false;
    }
  };
  bool ExplicitBitVect::GetBit(const unsigned int which) const {
    if(which < 0 || which >= d_size){
      throw IndexErrorException(which);
    }
    return((bool)(*dp_bits)[which]);
  };

  ExplicitBitVect ExplicitBitVect::operator^ (const ExplicitBitVect &other) const {
    ExplicitBitVect ans(d_size);
    *(ans.dp_bits) = (*dp_bits) ^ *(other.dp_bits);
    return(ans);
  };

  ExplicitBitVect ExplicitBitVect::operator& (const ExplicitBitVect &other) const {
    ExplicitBitVect ans(d_size);
    *(ans.dp_bits) = (*dp_bits) & *(other.dp_bits);
    return(ans);
  };

  ExplicitBitVect ExplicitBitVect::operator| (const ExplicitBitVect &other) const {
    ExplicitBitVect ans(d_size);
    *(ans.dp_bits) = (*dp_bits) | *(other.dp_bits);
    return(ans);
  };
  
  ExplicitBitVect ExplicitBitVect::operator~ () const {
    ExplicitBitVect ans(d_size);
    *(ans.dp_bits) = ~(*dp_bits);
    return(ans);
  };

  const unsigned int ExplicitBitVect::GetNumBits() const {
    return d_size;
  };
  const unsigned int ExplicitBitVect::GetNumOnBits() const {
    return dp_bits->count();
  };
  const unsigned int ExplicitBitVect::GetNumOffBits() const {
    return d_size - dp_bits->count();
  };

  // the contents of v are blown out
  void ExplicitBitVect::GetOnBits (IntVect& v) const {
    unsigned int nOn = GetNumOnBits();
    if(!v.empty()) IntVect().swap(v);
    v.reserve(nOn);
    for(unsigned int i=0;i<d_size;i++){
      if((bool)(*dp_bits)[i]) v.push_back(i);
    }
  };

  void ExplicitBitVect::_InitForSize(unsigned int size) {
    d_size = size;
    if(dp_bits) delete dp_bits;
    dp_bits = new boost::dynamic_bitset<>(size);
  };


  ExplicitBitVect::~ExplicitBitVect() {
    if(dp_bits) delete dp_bits;
  };

std::string
ExplicitBitVect::ToString() const
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

  unsigned int nOnBits=GetNumOnBits();
  int tVers = ci_BITVECT_VERSION*-1;
  ss.write((const char *)&(tVers),sizeof(tVers));
  ss.write((const char *)&d_size,sizeof(d_size));
  ss.write((const char *)&nOnBits,sizeof(nOnBits));
  
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


