// $Id$
//
//  Copyright (C) 2007 Greg Landrum
//
//  @@ All Rights Reserved @@
//
#ifndef __RD_SPARSE_INT_VECT_20070921__
#define __RD_SPARSE_INT_VECT_20070921__

#include <map>
#include <string>
#include <RDGeneral/Invariant.h>
#include <sstream>
#include <RDBoost/Exceptions.h>

const int ci_SPARSEINTVECT_VERSION=0x0001; //!< version number to use in pickles
namespace RDKit{
  //! a class for efficiently storing sparse vectors of ints
  template <typename IndexType>
  class SparseIntVect {
  public:
    typedef std::map<IndexType,int> StorageType;
  
    SparseIntVect() : d_length(0) {};

    //! initialize with a particular length
    SparseIntVect(IndexType length) : d_length(length) {};

    //! Copy constructor
    SparseIntVect(const SparseIntVect<IndexType> &other){
      d_length=other.d_length;
      d_data.insert(other.d_data.begin(),other.d_data.end());
    }

    //! constructor from a pickle
    SparseIntVect(const std::string pkl){
      initFromText(pkl.c_str(),pkl.size());
    };
    //! constructor from a pickle
    SparseIntVect(const char *pkl,const unsigned int len){
      initFromText(pkl,len);
    };

    //! destructor (doesn't need to do anything)
    ~SparseIntVect() {}

    //! return the value at an index
    int getVal(IndexType idx) const {
      if(idx<0||idx>=d_length){
	throw IndexErrorException(idx);
      }
      int res=0;
      typename StorageType::const_iterator iter=d_data.find(idx);
      if(iter!=d_data.end()){
	res=iter->second;
      }
      return res;
    };
    int operator[] (IndexType idx) const { return getVal(idx); };

    //! set the value at an index
    void setVal(IndexType idx, int val){
      if(idx<0||idx>=d_length){
	throw IndexErrorException(idx);
      }
      if(val!=0){
	d_data[idx]=val;
      } else {
	d_data.erase(idx);
      }
    };
    //! returns the length
    IndexType getLength() const { return d_length; };

    //! returns the sum of all the elements in the vect
    int getTotalVal() const {
      int res=0;
      typename StorageType::const_iterator iter;
      for(iter=d_data.begin();iter!=d_data.end();++iter){
	res+=iter->second;
      }
      return res;
    };


    //! returns our nonzero elements as a map(IndexType->int)
    const StorageType &getNonzeroElements() const {
      return d_data;
    }


    //! this is a "fuzzy" intesection, the final value
    //! of each element is equal to the minimum from
    //! the two vects.
    SparseIntVect<IndexType> &
    operator&= (const SparseIntVect<IndexType> &other) {
      if(other.d_length!=d_length){
	throw ValueErrorException("SparseIntVect size mismatch");
      }

      typename StorageType::iterator iter=d_data.begin();
      typename StorageType::const_iterator oIter=other.d_data.begin();
      while(iter!=d_data.end()){
	// we're relying on the fact that the maps are sorted:
	while(oIter!=other.d_data.end() && oIter->first < iter->first){
	  ++oIter;
	}
	if(oIter!=d_data.end() && oIter->first==iter->first){
	  // found it:
	  if(oIter->second<iter->second){
	    iter->second=oIter->second;
	  }
	  ++oIter;
          ++iter;
	} else {
	  // not there; our value is zero, which means
	  // we should remove this value:
	  typename StorageType::iterator tmpIter=iter;
	  ++tmpIter;
	  d_data.erase(iter);
	  iter=tmpIter;
	}
      }
      return *this;
    };
    const SparseIntVect<IndexType> 
    operator& (const SparseIntVect<IndexType> &other) const {
      SparseIntVect<IndexType> res(*this);
      return res&=other;
    }

    //! this is a "fuzzy" union, the final value
    //! of each element is equal to the maximum from
    //! the two vects.
    SparseIntVect<IndexType> &
    operator|= (const SparseIntVect<IndexType> &other) {
      if(other.d_length!=d_length){
	throw ValueErrorException("SparseIntVect size mismatch");
      }

      typename StorageType::iterator iter=d_data.begin();
      typename StorageType::const_iterator oIter=other.d_data.begin();
      while(iter!=d_data.end()){
	// we're relying on the fact that the maps are sorted:
	while(oIter!=other.d_data.end() &&
	      oIter->first < iter->first){
	  d_data[oIter->first]=oIter->second;
	  ++oIter;
	}
	if(oIter!=other.d_data.end() && oIter->first==iter->first){
	  // found it:
	  if(oIter->second>iter->second){
	    iter->second=oIter->second;
	  }
	  ++oIter;
	}
	++iter;
      }
      // finish up the other vect:
      while(oIter!=other.d_data.end()){
	d_data[oIter->first]=oIter->second;
	++oIter;
      }
      return *this;
    };
    const SparseIntVect<IndexType> 
    operator| (const SparseIntVect<IndexType> &other) const {
      SparseIntVect<IndexType> res(*this);
      return res|=other;
    }

    SparseIntVect<IndexType> &
    operator+= (const SparseIntVect<IndexType> &other) {
      if(other.d_length!=d_length){
	throw ValueErrorException("SparseIntVect size mismatch");
      }
      typename StorageType::iterator iter=d_data.begin();
      typename StorageType::const_iterator oIter=other.d_data.begin();
      while(oIter!=other.d_data.end()){
	while(iter!=d_data.end() &&
	      iter->first < oIter->first){
	  ++iter;
	}
	if(iter!=d_data.end() && oIter->first==iter->first){
	  // found it:
	  iter->second+=oIter->second;
	  if(!iter->second){
	    typename StorageType::iterator tIter=iter;
	    ++tIter;
	    d_data.erase(iter);
	    iter=tIter;
	  } else {
	    ++iter;
	  }
	} else {
	  d_data[oIter->first]=oIter->second;
	}
	++oIter;
      }
      return *this;
    };
    const SparseIntVect<IndexType> 
    operator+ (const SparseIntVect<IndexType> &other) const {
      SparseIntVect<IndexType> res(*this);
      return res+=other;
    }

    SparseIntVect<IndexType> &
    operator-= (const SparseIntVect<IndexType> &other) {
      if(other.d_length!=d_length){
	throw ValueErrorException("SparseIntVect size mismatch");
      }
      typename StorageType::iterator iter=d_data.begin();
      typename StorageType::const_iterator oIter=other.d_data.begin();
      while(oIter!=other.d_data.end()){
	while(iter!=d_data.end() &&
	      iter->first < oIter->first){
	  ++iter;
	}
	if(iter!=d_data.end() && oIter->first==iter->first){
	  // found it:
	  iter->second-=oIter->second;
	  if(!iter->second){
	    typename StorageType::iterator tIter=iter;
	    ++tIter;
	    d_data.erase(iter);
	    iter=tIter;
	  } else {
	    ++iter;
	  }
	} else {
	  d_data[oIter->first]=-oIter->second;
	}
	++oIter;
      }
      return *this;
    };
    const SparseIntVect<IndexType> 
    operator- (const SparseIntVect<IndexType> &other) const {
      SparseIntVect<IndexType> res(*this);
      return res-=other;
    }

    bool operator==(const SparseIntVect<IndexType> &v2){
      if(d_length!=v2.d_length){
	return false;
      }
      return d_data==v2.d_data;
    }
    bool operator!=(const SparseIntVect<IndexType> &v2){
      return !(*this==v2);
    }

    //! returns a binary string representation (pickle)
    const std::string toString() const {
      std::stringstream ss(std::ios_base::binary|std::ios_base::out|std::ios_base::in);
      ss.write((const char *)&(ci_SPARSEINTVECT_VERSION),sizeof(ci_SPARSEINTVECT_VERSION));
      unsigned int pieceSize=sizeof(IndexType);
      ss.write((const char *)&pieceSize,sizeof(pieceSize));
      ss.write((const char *)&d_length,sizeof(d_length));
      IndexType nEntries=d_data.size();
      ss.write((const char *)&nEntries,sizeof(nEntries));

      typename StorageType::const_iterator iter=d_data.begin();
      while(iter!=d_data.end()){
	ss.write((const char *)&iter->first,sizeof(iter->first));
	ss.write((const char *)&iter->second,sizeof(iter->second));
	++iter;
      }
      return ss.str();
    };

    void fromString(std::string &txt) {
      initFromText(txt.c_str(),txt.length());
    }

  private:
    IndexType d_length;
    StorageType d_data;
    
    void initFromText(const char *pkl,const unsigned int len) {
      d_data.clear();
      std::stringstream ss(std::ios_base::binary|std::ios_base::out|std::ios_base::in);
      ss.write(pkl,len);
      
      int vers;
      ss.read((char *)&vers,sizeof(vers));
      if(vers==0x0001){
	unsigned int idxSize;
	ss.read((char *)&idxSize,sizeof(idxSize));
	if(idxSize>sizeof(IndexType)){
	  throw ValueErrorException("IndexType cannot accomodate index size in SparseIntVect pickle");
	}
	switch(idxSize){
	case 1:
	  readVals<unsigned char>(ss);break;
	case 4:
	  readVals<unsigned int>(ss);break;
	case 8:
	  readVals<unsigned long long>(ss);break;
	default:
	  throw ValueErrorException("unreadable format");
	}
      } else {
	throw ValueErrorException("bad version in SparseIntVect pickle");
      }
    };
    template <typename T>
    void readVals(std::stringstream &ss){
      PRECONDITION(sizeof(T)<=sizeof(IndexType),"invalid size");
      T tVal;
      ss.read((char *)&tVal,sizeof(T));
      d_length=tVal;
      T nEntries;
      ss.read((char *)&nEntries,sizeof(T));
      for(T i=0;i<nEntries;++i){
	ss.read((char *)&tVal,sizeof(tVal));
	int val;
	ss.read((char *)&val,sizeof(val));
	d_data[tVal]=val;
      }
    }
  };

  template <typename IndexType, typename SequenceType>
  void updateFromSequence(SparseIntVect<IndexType> &vect,
			  const SequenceType &seq){
    typename SequenceType::const_iterator seqIt;
    for(seqIt=seq.begin();seqIt!=seq.end();++seqIt){
      // EFF: probably not the most efficient approach
      IndexType idx=*seqIt;
      vect.setVal(idx,vect.getVal(idx)+1);
    }
  }

  template <typename IndexType>
  double DiceSimilarity(const SparseIntVect<IndexType> &v1,
			const SparseIntVect<IndexType> &v2,
			double bounds=0.0){
    double v1Sum=v1.getTotalVal();
    double v2Sum=v2.getTotalVal();
    double denom=v1Sum+v2Sum;
    if(fabs(denom)<1e-6){
      return 0.0;
    }
    if(bounds>0.0){
      double minV=v1Sum<v2Sum?v1Sum:v2Sum;
      if(2.*minV/denom<bounds){
	return 0.0;
      }
    }
    double numer=(v1&v2).getTotalVal();
    return 2.*numer/denom;
  }
} 



#endif
