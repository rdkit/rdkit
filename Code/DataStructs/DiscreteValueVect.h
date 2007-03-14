//
//  Copyright (C) 2004-2007 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
#ifndef __RD_DISCRETE_VALUE_VECT_20050124__
#define __RD_DISCRETE_VALUE_VECT_20050124__

#include <boost/smart_ptr.hpp>
#include <string>

namespace RDKit{
  // we are making an assumption here that and unsigned int is 32 bits long
  const unsigned int BITS_PER_INT=32;

  class DiscreteValueVect {
  public:
    typedef boost::shared_array<unsigned int> DATA_SPTR;
  
    typedef enum {
      ONEBITVALUE=0,
      TWOBITVALUE,
      FOURBITVALUE,
      EIGHTBITVALUE,
      SIXTEENBITVALUE,
    } DiscreteValueType;

    DiscreteValueVect(DiscreteValueType valType, unsigned int length) : d_type(valType), d_length(length) {
      d_bitsPerVal = (1 << static_cast<unsigned int>(valType));
      d_valsPerInt = BITS_PER_INT/d_bitsPerVal;
      d_numInts = (length + d_valsPerInt -1)/d_valsPerInt;
      d_mask = ((1<<d_bitsPerVal) -1);
      unsigned int *data = new unsigned int[d_numInts];
      memset(static_cast<void *>(data),0,d_numInts*sizeof(unsigned int));
      d_data.reset(data);
    }

    //! Copy constructor
    DiscreteValueVect(const DiscreteValueVect& other);

    //! constructors from pickles
    DiscreteValueVect(const std::string pkl){
      initFromText(pkl.c_str(),pkl.size());
    };
    DiscreteValueVect(const char *pkl,const unsigned int len){
      initFromText(pkl,len);
    };

    ~DiscreteValueVect() {}

    unsigned int getVal(unsigned int i) const;
    void setVal(unsigned int i, unsigned int val);
    unsigned int getTotalVal() const;

    unsigned int getLength() const;

    const unsigned int *getData() const;
    unsigned int getNumBitsPerVal() const {
      return d_bitsPerVal;
    }

    DiscreteValueType getValueType() const {
      return d_type;
    }

    unsigned int getNumInts() const {
      return d_numInts;
    }

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
}
#endif
