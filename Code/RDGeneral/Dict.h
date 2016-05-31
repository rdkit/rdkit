//
// Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file Dict.h

  \brief Defines the Dict class

*/
#ifndef __RD_DICT_H__
#define __RD_DICT_H__

#include <map>
#include <string>
#include <vector>
#include "RDValue.h"
#include "Exceptions.h"
#include "tags.h"

#include "BoostStartInclude.h"
#include <boost/lexical_cast.hpp>
#include "BoostEndInclude.h"

namespace RDKit {
typedef std::vector<std::string> STR_VECT;

//! \brief The \c Dict class can be used to store objects of arbitrary
//!        type keyed by \c strings.
//!
//!  The actual storage is done using \c RDValue objects.
//!
class Dict {
public:
#ifdef RDVALUE_HAS_KEY
  // The RDValue has a .getKey implementation
  //  so use it to save some space.
  struct Pair {
    RDValue val;

   Pair() : val() {}
   Pair(int k, const RDValue &v) : val(v) {
     val.setKey(k);
     TEST_ASSERT(val.getKey() == k);
   }
   Pair(const std::string &s, RDValue_cast_t v) :
    val(v) { val.setKey(tagmap.get(s));}

    inline void setKey(int k) { val.setKey(k); }
    inline int  getKey() const { return val.getKey(); }
  };  
#else
  // The RDValue has no getKey implementation
  //  so implement our own
  struct Pair {
    int key;
    RDValue val;

   Pair() : key(), val() {}
   Pair(int k, const RDValue &v) : key(k), val(v) {}
   Pair(const std::string &s, RDValue_cast_t v) :
    key(tagmap.get(s)), val(v) {}

    inline void setKey(int k) { key = k; }
    inline int  getKey() const { return key; }
  };
#endif  
  typedef std::vector<Pair> DataType;
  static RDTags tagmap;

  Dict() : _data(), _hasNonPodData(false) {  };

  Dict(const Dict &other) : _data(other._data) {
    _hasNonPodData = other._hasNonPodData;
    if (_hasNonPodData) {
      std::vector<Pair> data(other._data.size());
      _data.swap(data);
      for (size_t i=0; i< _data.size(); ++i) {
        _data[i].setKey(other._data[i].getKey());
        copy_rdvalue(_data[i].val, other._data[i].val);
      }
    }   
  }
  
  ~Dict() {
    reset(); // to clear pointers if necessary
  }
  
  Dict &operator=(const Dict &other) {
    _hasNonPodData = other._hasNonPodData;
    if (_hasNonPodData) {
      std::vector<Pair> data(other._data.size());
      _data.swap(data);
      for (size_t i=0; i< _data.size(); ++i) {
        _data[i].setKey(other._data[i].getKey());
        copy_rdvalue(_data[i].val, other._data[i].val);
      }
    } else {
      _data = other._data;      
    }    
    return *this;
  };

  //----------------------------------------------------------
  //! \brief Returns whether or not the dictionary contains a particular
  //!        key.
  bool hasVal(int what) const {
    for(size_t i=0 ; i< _data.size(); ++i) {
      if (_data[i].getKey() == what ) return true;
    }
    return false;
  };

  bool hasVal(const std::string & what) const {
    return hasVal(tagmap.get(what));
  }

  //----------------------------------------------------------
  //! Returns the set of keys in the dictionary
  /*!
     \return  a \c STR_VECT
  */
  STR_VECT keys() const {
    STR_VECT res;
    DataType::const_iterator item;
    for (item = _data.begin(); item != _data.end(); item++) {
      res.push_back(tagmap.get(item->getKey()));
    }
    return res;
  }

  //----------------------------------------------------------
  //! \brief Gets the value associated with a particular key
  /*!
     \param what  the key to lookup
     \param res   a reference used to return the result

     <B>Notes:</b>
      - If \c res is a \c std::string, every effort will be made
        to convert the specified element to a string using the
        \c boost::lexical_cast machinery.
      - If the dictionary does not contain the key \c what,
        a KeyErrorException will be thrown.
  */
  template <typename T>
  void getVal(const std::string &what, T &res) const {
    res = getVal<T>(what);
  };
  
  template <typename T>
  void getVal(int tag, T &res) const {
    res = getVal<T>(tag);
  };

  //! \overload
  template <typename T>
  T getVal(int tag) const {
    for(size_t i=0; i< _data.size(); ++i) {
      if (_data[i].getKey() == tag) {
        return from_rdvalue<T>(_data[i].val);
      }
    }
    throw KeyErrorException(common_properties::getPropName(tag));
  }

  template <typename T>
  T getVal(const std::string &what) const {
    return getVal<T>(tagmap.get(what));
  }

  //! \overload
  void getVal(const std::string &what, std::string &res) const {
    return getVal(tagmap.get(what), res);
  }

  void getVal(int tag, std::string &res) const;

  //----------------------------------------------------------
  //! \brief Potentially gets the value associated with a particular key
  //!        returns true on success/false on failure.
  /*!
     \param what  the key to lookup
     \param res   a reference used to return the result

     <B>Notes:</b>
      - If \c res is a \c std::string, every effort will be made
        to convert the specified element to a string using the
        \c boost::lexical_cast machinery.
      - If the dictionary does not contain the key \c what,
        a KeyErrorException will be thrown.
  */

  template <typename T>
  bool getValIfPresent(const std::string &what, T &res) const {
    return getValIfPresent(tagmap.get(what), res);
  }
  
  template <typename T>
  bool getValIfPresent(int tag, T &res) const {
    for(size_t i=0; i< _data.size(); ++i) {
      if (_data[i].getKey() == tag) {
        res = from_rdvalue<T>(_data[i].val);
        return true;
      }
    }
    return false;
  };


  //! \overload
  bool getValIfPresent(int tag, std::string &res) const;
  
  bool getValIfPresent(const std::string &what, std::string &res) const {
    return getValIfPresent(tagmap.get(what), res);
  }

  //----------------------------------------------------------
  //! \brief Sets the value associated with a key
  /*!

     \param what the key to set
     \param val  the value to store

     <b>Notes:</b>
        - If \c val is a <tt>const char *</tt>, it will be converted
           to a \c std::string for storage.
        - If the dictionary already contains the key \c what,
          the value will be replaced.
  */
  template <typename T>
  void setVal(int what, T &val) {
    _hasNonPodData = true;
    for(size_t i=0; i< _data.size(); ++i) {
      if (_data[i].getKey() == what) {
        _data[i].val = val;
#ifdef RDVALUE_HAS_KEY
        _data[i].setKey(what);
#endif        
        return;
      }
    }
    _data.push_back(Pair(what, val));
  };

  template <typename T>
  void setVal(const std::string &what, T &val) {
    setVal(tagmap.get(what), val);
  }

  template <typename T>
  void setPODVal(int what, T val) {
    // don't change the hasNonPodData status
    for(size_t i=0; i< _data.size(); ++i) {
      if (_data[i].getKey() == what) {
        _data[i].val = val;
#ifdef RDVALUE_HAS_KEY
        _data[i].setKey(what);
#endif        
        return;
      }
    }
    _data.push_back(Pair(what, val));
  };

  template <typename T>
  void setPODVal(const std::string &what, T &val) {
    setPODVal(tagmap.get(what), val);
  }
  
  void setVal(const std::string &what, bool val) {
    setPODVal(what, val);
  }

  void setVal(int what, bool val) {
    setPODVal(what, val);
  }

  void setVal(const std::string &what, double val) {
    setPODVal(what, val);
  }

  void setVal(int what, double val) {
    setPODVal(what, val);
  }
  
  void setVal(const std::string &what, float val) {
    setPODVal(what, val);
  }

  void setVal(int what, float val) {
    setPODVal(what, val);
  }

  void setVal(const std::string &what, int val) {
    setPODVal(what, val);
  }

  void setVal(int what, int val) {
    setPODVal(what, val);
  }
  
  void setVal(const std::string &what, unsigned int val) {
    setPODVal(what, val);
  }

  void setVal(int what, unsigned int val) {
    setPODVal(what, val);
  }
  
  //! \overload
  void setVal(const std::string &what, const char *val) {
    std::string h(val);
    setVal(what, h);
  }

  //! \overload
  void setVal(int what, const char *val) {
    std::string h(val);
    setVal(what, h);
  }
  
  //----------------------------------------------------------
  //! \brief Clears the value associated with a particular key,
  //!     removing the key from the dictionary.
  /*!

     \param what the key to clear

   <b>Notes:</b>
      - If the dictionary does not contain the key \c what,
        a KeyErrorException will be thrown.
  */
  void clearVal(const std::string &what) {
    clearVal(tagmap.get(what));
  }
  
  void clearVal(int tag) {
    for(DataType::iterator it = _data.begin(); it < _data.end() ; ++it) {
      if (it->getKey() == tag) {
        _data.erase(it);
        return;
      }
    }
    throw KeyErrorException(common_properties::getPropName(tag));
  };
  
  //----------------------------------------------------------
  //! \brief Clears all keys (and values) from the dictionary.
  //!
  void reset() {
    if (_hasNonPodData) {
      for (size_t i=0; i< _data.size(); ++i) {
        RDValue::cleanup_rdvalue(_data[i].val);
      }
    }
    DataType data;
    _data.swap(data);
  };

 private:
  DataType _data;  //!< the actual dictionary
  bool     _hasNonPodData; // if true, need a deep copy
                           //  (copy_rdvalue)

 public:
  // Iteration API (over Pair)
  typedef DataType::iterator iterator;
  typedef DataType::const_iterator const_iterator;

  iterator begin() { return _data.begin(); }
  iterator end() { return _data.end(); }
  const_iterator begin() const { return _data.begin(); }
  const_iterator end()   const { return _data.end(); }
  
};
}
#endif
