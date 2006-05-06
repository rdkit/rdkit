//
// Copyright (C) 2003-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
#ifndef __RD_DICT_H__
#define __RD_DICT_H__


#include <map>
#include <string>
#include <vector>
#include <boost/any.hpp>
#include <RDBoost/Exceptions.h>

namespace RDKit{

  typedef std::vector<std::string>                  STR_VECT;

  //! \brief The \c Dict class can be used to store objects of arbitrary
  //!        type keyed by \c strings.
  //!
  //!  The actual storage is done using \c boost::any objects.
  //!
  class Dict {
    //! \brief this function is used solely to force the instantiation of particular
    //!    types of the \c Dict.toAny() and \c Dict.fromAny() methods in order to
    //!    avoid link errors.
    friend void force_types();
  public:
    typedef std::map<const std::string, boost::any> DataType;
    Dict();
    Dict(const Dict &other);
    Dict &operator=(const Dict &other);

    //----------------------------------------------------------
    //! \brief Returns whether or not the dictionary contains a particular
    //!        key.
    bool hasVal(const char *what) const;
    bool hasVal(const std::string &what) const;

    //----------------------------------------------------------
    //! Returns the set of keys in the dictionary
    /*!
       \return  a \c STR_VECT
    */
    STR_VECT keys() const;

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
    void getVal(const std::string &what,T &res) const {
      DataType::const_iterator pos=_data.find(what);
      if(pos==_data.end())
	throw KeyErrorException(what);
      const boost::any &val = pos->second;
      res = fromany<T>(val);
    };

    //! \overload
    template <typename T>
    void getVal(const char *what,T &res) const {
      std::string key(what);
      getVal(key, res);
    };

    //! \overload
    void getVal(const std::string &what, std::string &res) const;


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
    void setVal(const std::string &what, T &val){
      std::string key = what;
      _data[key] = toany(val);
    };
    //! \overload
    template <typename T>
    void setVal(const char *what, T &val){
      std::string key = what;
      setVal(key,val);
    };
    //! \overload
    void setVal(const std::string &what, const char *val){
      std::string h(val);
      setVal(what,h);
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
    void clearVal(const std::string &what);

    //! \overload
    void clearVal(const char *what) {
      std::string key=what;
      clearVal(key);
    };

    //----------------------------------------------------------
    //! \brief Clears all keys (and values) from the dictionary.
    //!
    void reset();

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
