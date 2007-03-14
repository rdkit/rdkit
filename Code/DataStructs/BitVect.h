//
// Copyright (c) 2003-2006 greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
#ifndef __RD_BITVECT_H__
#define __RD_BITVECT_H__

#include <vector>
using std::vector;
typedef vector<int> IntVect;
typedef IntVect::iterator IntVectIter;
typedef vector<double> DoubleVect;
typedef DoubleVect::iterator DoubleVectIter;
const int ci_BITVECT_VERSION=0x0020; //!< version number to use in pickles

//! Abstract base class for storing BitVectors
class BitVect{
public:
  virtual ~BitVect() = 0;
  //! sets a particular bit and returns its original value
  virtual bool SetBit(const unsigned int which) = 0;
  //! unsets a particular bit and returns its original value
  virtual bool UnSetBit(const unsigned int which) = 0;
  //! returns the value of a particular bit
  virtual bool GetBit(const unsigned int which) const = 0;
  //! returns the number of bits (the length of the BitVect)
  virtual const unsigned int GetNumBits() const = 0;
  //! returns the number of on bits
  virtual const unsigned int GetNumOnBits() const = 0;
  //! returns the number of off bits
  virtual const unsigned int GetNumOffBits() const =0;
  //! replaces the contents of \c v with indices of our on bits
  virtual void GetOnBits (IntVect& v) const = 0;
  //! clears (sets to off) all of our bits
  virtual void ClearBits() = 0;

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
  void InitFromText(const char *data,const unsigned int dataLen,
                    bool isBase64=false,bool allowOldFormat=false);
                    
  //! returns a serialized (pickled) version of this BitVect
  virtual std::string ToString() const = 0;

  virtual bool operator[] (const unsigned int which) const = 0;

private:
  virtual void _InitForSize(const unsigned int size) = 0;
};


#endif
