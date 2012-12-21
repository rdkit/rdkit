//
//  Copyright (C) 2004-2008 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_DISCRETE_VALUE_VECT_20050124__
#define __RD_DISCRETE_VALUE_VECT_20050124__

#include <boost/smart_ptr.hpp>
#include <string>
#include <cstring>
#include <boost/cstdint.hpp>

namespace RDKit{
  // we require 32bit unsigneds using the boost::uint32_t type:
  const unsigned int BITS_PER_INT=32;

  //! a class for efficiently storing vectors of discrete values
  class DiscreteValueVect {
  public:
    typedef boost::shared_array<boost::uint32_t> DATA_SPTR;
  
    //! used to define the possible range of the values
    typedef enum {
      ONEBITVALUE=0,
      TWOBITVALUE,
      FOURBITVALUE,
      EIGHTBITVALUE,
      SIXTEENBITVALUE,
    } DiscreteValueType;

    //! initialize with a particular type and size
    DiscreteValueVect(DiscreteValueType valType, unsigned int length) : d_type(valType), d_length(length) {
      d_bitsPerVal = (1 << static_cast<unsigned int>(valType));
      d_valsPerInt = BITS_PER_INT/d_bitsPerVal;
      d_numInts = (length + d_valsPerInt -1)/d_valsPerInt;
      d_mask = ((1<<d_bitsPerVal) -1);
      boost::uint32_t *data = new boost::uint32_t[d_numInts];
      memset(static_cast<void *>(data),0,d_numInts*sizeof(boost::uint32_t));
      d_data.reset(data);
    }

    //! Copy constructor
    DiscreteValueVect(const DiscreteValueVect& other);

    //! constructor from a pickle
    DiscreteValueVect(const std::string pkl){
      initFromText(pkl.c_str(),pkl.size());
    };
    //! constructor from a pickle
    DiscreteValueVect(const char *pkl,const unsigned int len){
      initFromText(pkl,len);
    };

    ~DiscreteValueVect() {}

    //! return the value at an index
    unsigned int getVal(unsigned int i) const;

    //! support indexing using []
    int operator[] (unsigned int idx) const { return getVal(idx); };

    //! set the value at an index
    /*!
      NOTE: it is an error to have val > the max value this
      DiscreteValueVect can accomodate 
    */
    void setVal(unsigned int i, unsigned int val);

    //! returns the sum of all the elements in the vect
    unsigned int getTotalVal() const;

    //! returns the length
    unsigned int getLength() const;
    //! returns the length
    unsigned int size() const { return getLength(); };

    //! return a pointer to our raw data storage
    const boost::uint32_t *getData() const;

    //! return the number of bits used to store each value
    unsigned int getNumBitsPerVal() const {
      return d_bitsPerVal;
    }

    //! return the type of value being stored
    DiscreteValueType getValueType() const {
      return d_type;
    }

    //! returns the size of our storage
    unsigned int getNumInts() const {
      return d_numInts;
    }

    //! support dvv3 = dvv1&dvv2
    /*!

       operator& returns the minimum value for each element.
       e.g.:
         [0,1,2,0] & [0,1,1,1] -> [0,1,1,0]

    */
    DiscreteValueVect operator& (const DiscreteValueVect &other) const;
    //! support dvv3 = dvv1|dvv2
    /*!

       operator& returns the maximum value for each element.
       e.g.:
         [0,1,2,0] | [0,1,1,1] -> [0,1,2,1]

    */
    DiscreteValueVect operator| (const DiscreteValueVect &other) const;
    //DiscreteValueVect operator^ (const DiscreteValueVect &other) const;
    //DiscreteValueVect operator~ () const;


    DiscreteValueVect& operator+=(const DiscreteValueVect &other);
    DiscreteValueVect& operator-=(const DiscreteValueVect &other);

    //! returns a binary string representation (pickle)
    std::string toString() const;

  private:
    DiscreteValueType d_type;
    unsigned int d_bitsPerVal;
    unsigned int d_valsPerInt;
    unsigned int d_numInts;
    unsigned int d_length;
    unsigned int d_mask;
    DATA_SPTR d_data;

    void initFromText(const char *pkl,const unsigned int len);
  };

  unsigned int computeL1Norm(const DiscreteValueVect &v1, const DiscreteValueVect &v2);

  DiscreteValueVect operator+ (const DiscreteValueVect& p1,
			       const DiscreteValueVect& p2);
  DiscreteValueVect operator- (const DiscreteValueVect& p1,
			       const DiscreteValueVect& p2);

} 



#endif
