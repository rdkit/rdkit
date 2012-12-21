//
// Copyright (c) 2003-2008 greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_BITVECT_H__
#define __RD_BITVECT_H__

#include <vector>
#include <string>

typedef std::vector<int> IntVect;
typedef IntVect::iterator IntVectIter;
typedef std::vector<double> DoubleVect;
typedef DoubleVect::iterator DoubleVectIter;
const int ci_BITVECT_VERSION=0x0020; //!< version number to use in pickles

//! Abstract base class for storing BitVectors
class BitVect{
public:
  virtual ~BitVect() = 0;
  //! sets a particular bit and returns its original value
  virtual bool setBit(const unsigned int which) = 0;
  //! unsets a particular bit and returns its original value
  virtual bool unsetBit(const unsigned int which) = 0;
  //! returns the value of a particular bit
  virtual bool getBit(const unsigned int which) const = 0;
  //! returns the number of bits (the length of the BitVect)
  virtual unsigned int getNumBits() const = 0;
  //! returns the number of on bits
  virtual unsigned int getNumOnBits() const = 0;
  //! returns the number of off bits
  virtual unsigned int getNumOffBits() const =0;
  //! replaces the contents of \c v with indices of our on bits
  virtual void getOnBits (IntVect& v) const = 0;
  //! clears (sets to off) all of our bits
  virtual void clearBits() = 0;

  //! initializes this BitVect from a pickle
  /*!
    \param data     the raw pickle data
    \param dataLen  the length of \c data 
    \param isBase64 (optional) if this is set, \c data is assumed to
         be base64 encoded.
    \param allowOldFormat (optional) allows a very old form of the BitVect
         representation to be recognized. This argument disables a large
         amount of error checking and it is strongly suggested that it not
         be used in client code.
   */
  void initFromText(const char *data,const unsigned int dataLen,
                    bool isBase64=false,bool allowOldFormat=false);
                    
  //! returns a serialized (pickled) version of this BitVect
  virtual std::string toString() const = 0;

  virtual bool operator[] (const unsigned int which) const = 0;
  unsigned int size() const { return getNumBits(); }

private:
  virtual void _initForSize(const unsigned int size) = 0;
};


#endif
