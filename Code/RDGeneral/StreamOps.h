//
//  Copyright (C) 2002-2008 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
#ifndef _RD_STREAMOPS_H
#define _RD_STREAMOPS_H

#include "types.h"
#include <string>
#include <sstream>
#include <iostream>
#include <boost/cstdint.hpp>
#include <boost/detail/endian.hpp>

namespace RDKit{
  // this code block for handling endian problems is from :
  // http://stackoverflow.com/questions/105252/how-do-i-convert-between-big-endian-and-little-endian-values-in-c
  enum EEndian
    {
      LITTLE_ENDIAN_ORDER,
      BIG_ENDIAN_ORDER,
#if defined(BOOST_LITTLE_ENDIAN)
      HOST_ENDIAN_ORDER = LITTLE_ENDIAN_ORDER
#elif defined(BOOST_BIG_ENDIAN)
      HOST_ENDIAN_ORDER = BIG_ENDIAN_ORDER
#else
#error "Failed to determine the system endian value"
#endif
    };

  // this function swap the bytes of values given it's size as a template
  // parameter (could sizeof be used?).
  template <class T, unsigned int size>
  inline T SwapBytes(T value)
  {
    union
    {
      T value;
      char bytes[size];
    } in, out;

    in.value = value;

    for (unsigned int i = 0; i < size / 2; ++i)
      {
        out.bytes[i] = in.bytes[size - 1 - i];
        out.bytes[size - 1 - i] = in.bytes[i];
      }

    return out.value;
  }
  
  // Here is the function you will use. Again there is two compile-time assertion
  // that use the boost librarie. You could probably comment them out, but if you
  // do be cautious not to use this function for anything else than integers
  // types. This function need to be calles like this :
  //
  //     int x = someValue;
  //     int i = EndianSwapBytes<HOST_ENDIAN_ORDER, BIG_ENDIAN_ORDER>(x);
  //
  template<EEndian from, EEndian to, class T>
  inline T EndianSwapBytes(T value)
  {
    // A : La donnée à swapper à une taille de 2, 4 ou 8 octets
    BOOST_STATIC_ASSERT(sizeof(T)==1 || sizeof(T) == 2 || sizeof(T) == 4 || sizeof(T) == 8);
    if(sizeof(T)==1) return value;

    // A : La donnée à swapper est d'un type arithmetic
    //BOOST_STATIC_ASSERT(boost::is_arithmetic<T>::value);

    // Si from et to sont du même type on ne swap pas.
    if (from == to)
      return value;

    return SwapBytes<T, sizeof(T)>(value);
  }
  template<EEndian from, EEndian to>
  inline char EndianSwapBytes(char value)
  {
    return value;
  }
  template<EEndian from, EEndian to>
  inline unsigned char EndianSwapBytes(unsigned char value)
  {
    return value;
  }
  template<EEndian from, EEndian to>
  inline signed char EndianSwapBytes(signed char value)
  {
    return value;
  }
  // --------------------------------------
  
  //! Packs an integer and outputs it to a stream
  inline void appendPackedIntToStream(std::stringstream &ss, boost::uint32_t num) {
    int nbytes, bix;
    unsigned int val, res;
    char tc;
    
    CHECK_INVARIANT(num >= 0, "");
    res=num;
    while (1) {
      if (res < (1<<7)) {
        val = (res<<1);
        nbytes = 1;
        break;
      }
      res -= (1<<7);
      if (res < (1<<14)) {
        val = ((res<<2) | 1);
        nbytes = 2;
        break;
      }
      res -= (1<<14);
      if (res < (1<<21)) {
        val = ((res<<3) | 3);
        nbytes = 3;
        break;
      }
      res -= (1<<21);
      if ( res < (1<<29)) {
        val = ((res<<3) | 7);
        nbytes = 4;
        break;
      }
      else {
        CHECK_INVARIANT(0, "ERROR: Integer too big to pack\n");
      }
    }
    //val = EndianSwapBytes<HOST_ENDIAN_ORDER,LITTLE_ENDIAN_ORDER>(val);
    
    for (bix = 0; bix < nbytes; bix++) {
      tc = (char) (val & 255);
      ss.write(&tc, 1);
      val >>= 8;
    }
  }
  
  //! Reads an integer from a stream in packed format and returns the result.
  inline boost::uint32_t readPackedIntFromStream(std::stringstream &ss) {
    boost::uint32_t val, num;
    int shift, offset;
    char tmp;
    ss.read(&tmp, sizeof(tmp));
    val = UCHAR(tmp);
    offset = 0;
    if ((val&1) == 0) {
      shift = 1;
    }
    else if ((val&3) == 1) {
      ss.read((char *)&tmp, sizeof(tmp));
      val |= (UCHAR(tmp) << 8);
      shift = 2;
      offset = (1<<7);
    }
    else if ((val&7) == 3) {
      ss.read((char *)&tmp, sizeof(tmp));
      val |= (UCHAR(tmp) << 8);
      ss.read((char *)&tmp, sizeof(tmp));
      val |= (UCHAR(tmp) << 16);
      shift = 3;
      offset = (1<<7) + (1<<14);
    }
    else {
      ss.read((char *)&tmp, sizeof(tmp));
      val |= (UCHAR(tmp) << 8);
      ss.read((char *)&tmp, sizeof(tmp));
      val |= (UCHAR(tmp) << 16);
      ss.read((char *)&tmp, sizeof(tmp));
      val |= (UCHAR(tmp) << 24);
      shift = 3;
      offset = (1<<7) + (1<<14) + (1<<21);
    }
    num = (val >> shift) + offset;
    //num = EndianSwapBytes<LITTLE_ENDIAN_ORDER,HOST_ENDIAN_ORDER>(num);
    return num;
  }

  //! Reads an integer from a char * in packed format and returns the result.
  //!  The argument is advanced
  inline boost::uint32_t pullPackedIntFromString(const char *&text) {
    boost::uint32_t val, num;
    int shift, offset;
    char tmp;
    tmp = *text;
    text++;
    val = UCHAR(tmp);
    offset = 0;
    if ((val&1) == 0) {
      shift = 1;
    }
    else if ((val&3) == 1) {
      tmp = *text;
      text++;
      val |= (UCHAR(tmp) << 8);
      shift = 2;
      offset = (1<<7);
    }
    else if ((val&7) == 3) {
      tmp = *text;
      text++;
      val |= (UCHAR(tmp) << 8);
      tmp = *text;
      text++;
      val |= (UCHAR(tmp) << 16);
      shift = 3;
      offset = (1<<7) + (1<<14);
    }
    else {
      tmp = *text;
      text++;
      val |= (UCHAR(tmp) << 8);
      tmp = *text;
      text++;
      val |= (UCHAR(tmp) << 16);
      tmp = *text;
      text++;
      val |= (UCHAR(tmp) << 24);
      shift = 3;
      offset = (1<<7) + (1<<14) + (1<<21);
    }
    num = (val >> shift) + offset;
    //num = EndianSwapBytes<LITTLE_ENDIAN_ORDER,HOST_ENDIAN_ORDER>(num);
    return num;
  }
  
  //! does a binary write of an object to a stream
  template <typename T>
    void streamWrite(std::ostream &ss,const T &val){
    T tval=EndianSwapBytes<HOST_ENDIAN_ORDER,LITTLE_ENDIAN_ORDER>(val);
    ss.write((const char *)&tval,sizeof(T));
  }
  //! does a binary read of an object from a stream
  template <typename T>
    void streamRead(std::istream &ss,T &loc){
    T tloc;
    ss.read((char *)&tloc,sizeof(T));
    loc = EndianSwapBytes<LITTLE_ENDIAN_ORDER,HOST_ENDIAN_ORDER>(tloc);
  }
 
  //! grabs the next line from an instream and returns it.
  inline std::string getLine(std::istream *inStream) {
    std::string res;
    std::getline(*inStream,res);
    if ((res.length() > 0) && (res[res.length()-1]=='\r')){
      res.erase(res.length()-1);
    }
    return res;
  }
  //! grabs the next line from an instream and returns it.
  inline std::string getLine(std::istream &inStream) {
    return getLine(&inStream);
  }
}



#endif
