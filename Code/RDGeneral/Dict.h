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
#include <cstring>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <boost/any.hpp>
#include <RDBoost/Exceptions.h>
#include <boost/lexical_cast.hpp>

namespace RDKit{
  typedef std::vector<std::string> STR_VECT;
  namespace detail {
    struct strcmpFunctor {
      bool operator()(const char *a, const char *b) const {
        return std::strcmp(a, b) < 0;
      }
    };
  }

  
  //! \brief The \c Dict class can be used to store objects of arbitrary
  //!        type keyed by \c strings.
  //!
  //!  The actual storage is done using \c boost::any objects.
  //!
  class Dict {
  public:
    typedef std::map<const char *, boost::any, detail::strcmpFunctor> DataType;
    Dict(){
      _data.clear();
    };

    Dict(const Dict &other) {
      _data = other._data;
    };

    Dict &operator=(const Dict &other) {
      _data = other._data;
      return *this;
    };

    //----------------------------------------------------------
    //! \brief Returns whether or not the dictionary contains a particular
    //!        key.
    bool hasVal(const char *what) const{
      return _data.find(what)!=_data.end();
    };
    bool hasVal(const std::string &what) const {
      return hasVal(what.c_str());
    };

    //----------------------------------------------------------
    //! Returns the set of keys in the dictionary
    /*!
       \return  a \c STR_VECT
    */
    STR_VECT keys() const {
      STR_VECT res;
      DataType::const_iterator item;
      for (item = _data.begin(); item != _data.end(); item++) {
        res.push_back(item->first);
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
    T getVal(const char *what) const {
      DataType::const_iterator pos=_data.find(what);
      if(pos==_data.end())
	throw KeyErrorException(what);
      const boost::any &val = pos->second;
      return fromany<T>(val);
    };
    //! \overload
    template <typename T>
    void getVal(const char *what,T &res) const {
      res = getVal<T>(what);
    };

    //! \overload
    template <typename T>
    void getVal(const std::string &what,T &res) const {
      res = getVal<T>(what.c_str());
    };
    //! \overload
    template <typename T>
    T getVal(const std::string &what) const {
      return getVal<T>(what.c_str());
    }

    //! \overload
    void getVal(const char *what, std::string &res) const;
    //! \overload
    void getVal(const std::string &what, std::string &res) const { getVal(what.c_str(),res); };

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
    void setVal(const char *what, T &val){
      _data[strdup(what)] = toany(val);
    };
    //! \overload
    template <typename T>
    void setVal(const std::string &what, T &val){
      setVal(what.c_str(),val);
    };
    //! \overload
    void setVal(const char *what, const char *val){
      std::string v(val);
      setVal(what,v);
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
    void clearVal(const char *what) {
      if(! this->hasVal(what) ) throw KeyErrorException(what);
      _data.erase(what);
    };
    //! \overload
    void clearVal(const std::string &what) {
      clearVal(what.c_str());
    };


    //----------------------------------------------------------
    //! \brief Clears all keys (and values) from the dictionary.
    //!
    void reset(){
      _data.clear();
    };

    //----------------------------------------------------------
    //! Converts a \c boost::any to type \c T
    /*!
       \param arg a \c boost::any reference

       \returns the converted object of type \c T
    */
    template <typename T>
      T fromany(const boost::any &arg) const;
    

    //----------------------------------------------------------
    //! Converts an instance of type \c T to \c boost::any
    /*!
       \param arg the object to be converted

       \returns a \c boost::any instance
    */
    template <typename T>
      boost::any toany(T arg) const;

  private:
    DataType _data; //!< the actual dictionary
  };
}
#endif
