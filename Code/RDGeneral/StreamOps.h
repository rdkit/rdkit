//
//  Copyright (C) 2002-2006 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
//
#ifndef _RD_STREAMOPS_H
#define _RD_STREAMOPS_H

#include "types.h"
#include <string>
#include <sstream>
#include <iostream>


namespace RDKit{
    
  //! Packs an integer and outputs it to a stream
  inline void appendPackedIntToStream(std::stringstream &ss, unsigned int num) {
    int nbytes, bix;
    unsigned int val, res;
    char tc;
    
    CHECK_INVARIANT(num >= 0, "");
    res = num;
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
        CHECK_INVARIANT(0, "ERROR: Integer to big to pack\n");
      }
    }
    
    for (bix = 0; bix < nbytes; bix++) {
      tc = (char) (val & 255);
      ss.write(&tc, 1);
      val >>= 8;
    }
  }
  
  //! Reads an integer from a stream in packed format and returns the result.
  inline unsigned int readPackedIntFromStream(std::stringstream &ss) {
    unsigned int val, num;
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
    return num;
  }

  //! Reads an integer from a char * in packed format and returns the result.
  //!  The argument is advanced
  inline unsigned int pullPackedIntFromString(const char *&text) {
    unsigned int val, num;
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
    return num;
  }
  
  //! does a binary write of an object to a stream
  template <typename T>
    void streamWrite(std::ostream &ss,const T &val){
    ss.write((const char *)&val,sizeof(T));
  }
  //! does a binary read of an object from a stream
  template <typename T>
    void streamRead(std::istream &ss,T &loc){
    ss.read((char *)&loc,sizeof(T));
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
