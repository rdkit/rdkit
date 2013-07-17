// $Id$
//
// Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Dict.h"


#include <boost/shared_array.hpp>
#include <boost/cstdint.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/utility.hpp>
#include <vector>
#include <list>
#include <iostream>
#include <sstream>

namespace RDKit{
  namespace {
    template <class T>
    std::string vectToString(const boost::any &val){
      const std::vector<T> &tv=boost::any_cast<std::vector<T> >(val);
      std::ostringstream sstr;
      sstr<<"[";
      std::copy(tv.begin(),tv.end(),std::ostream_iterator<T>(sstr,","));
      sstr<<"]";
      return sstr.str();
    }
  }
  
  
  void Dict::getVal(const std::string &what, std::string &res) const {
    //
    //  We're going to try and be somewhat crafty about this getVal stuff to make these
    //  containers a little bit more generic.  The normal behavior here is that the
    //  value being queried must be directly castable to type T.  We'll robustify a
    //  little bit by trying that and, if the cast fails, attempting a couple of 
    //  other casts, which will then be lexically cast to type T.
    //
    if(!hasVal(what) ) throw KeyErrorException(what);
    const boost::any &val = _data.find(what)->second;
    try{
      res = boost::any_cast<std::string>(val);
    } catch (const boost::bad_any_cast &) {
      if(val.type()==typeid(int)){
        res = boost::lexical_cast<std::string>(boost::any_cast<int>(val));
      } else if(val.type()==typeid(unsigned int)){
        res = boost::lexical_cast<std::string>(boost::any_cast<unsigned int>(val));
      } else if(val.type()==typeid(long)){
        res = boost::lexical_cast<std::string>(boost::any_cast<long>(val));
      } else if(val.type()==typeid(unsigned long)){
        res = boost::lexical_cast<std::string>(boost::any_cast<unsigned long>(val));
      } else if(val.type()==typeid(float)){
        res = boost::lexical_cast<std::string>(boost::any_cast<float>(val));
      } else if(val.type()==typeid(double)){
        res = boost::lexical_cast<std::string>(boost::any_cast<double>(val));
      } else if(val.type()==typeid(const char *)){
        res = std::string(boost::any_cast<const char *>(val));
      } else if(val.type()==typeid(std::vector<unsigned int>)){
        res = vectToString<unsigned int>(val);
      } else if(val.type()==typeid(std::vector<int>)){
        res = vectToString<int>(val);
      } else {
        throw;
      }
    }
  };

  namespace {
    template <class T>
    typename boost::enable_if<boost::is_arithmetic<T>, T>::type 
    _fromany(const boost::any &arg) {
      T res;
      if(arg.type()==typeid(std::string) || arg.type()==typeid(const char *)){
        try {
          res = boost::any_cast<T>(arg);
        } catch (const boost::bad_any_cast &exc) {
          try{
            res = boost::lexical_cast<T>(boost::any_cast<std::string>(arg));
          } catch (...){
            throw exc;
          }
        }
      } else {
        res = boost::any_cast<T>(arg);
      }
      return res;
    }
    template <class T>
    typename boost::disable_if<boost::is_arithmetic<T>, T>::type 
    _fromany(const boost::any &arg) {
      return boost::any_cast<T>(arg);
    }
  }
  template <typename T>
  T Dict::fromany(const boost::any &arg) const {
    return _fromany<T>(arg);
  };
  template <typename T>
  boost::any Dict::toany(T arg) const {
    return boost::any(arg);
  };

#define ANY_FORCE(T) template T Dict::fromany<T>(const boost::any& arg) const; \
                     template boost::any Dict::toany<T>(T arg) const;
  
  ANY_FORCE(bool);
  ANY_FORCE(boost::shared_array<double>);
  ANY_FORCE(boost::shared_array<int>);
  ANY_FORCE(double);
  ANY_FORCE(int);
  ANY_FORCE(std::list<int>);
  ANY_FORCE(std::string);
  ANY_FORCE(std::vector<boost::shared_array<double> >);
  ANY_FORCE(std::vector<boost::shared_array<int> >);
  ANY_FORCE(std::vector<double>);
  ANY_FORCE(std::vector<int>);
  ANY_FORCE(std::vector<std::list<int> >);
  ANY_FORCE(std::vector<std::string>);
  ANY_FORCE(std::vector<std::vector<double> >);
  ANY_FORCE(std::vector<std::vector<int> >);
  ANY_FORCE(std::vector<unsigned int>);
  ANY_FORCE(std::vector<unsigned long long>);
  ANY_FORCE(unsigned int);

  template const std::string & Dict::fromany<const std::string &>(const boost::any &arg) const;
  
  typedef boost::tuples::tuple<boost::uint32_t, boost::uint32_t, boost::uint32_t> uint32_t_tuple ;
  typedef boost::tuples::tuple<double, double, double> double_tuple;

  template uint32_t_tuple Dict::fromany<uint32_t_tuple>(const boost::any& arg) const;
  template boost::any Dict::toany<uint32_t_tuple>(uint32_t_tuple arg) const;
  template double_tuple Dict::fromany<double_tuple>(const boost::any& arg) const;
  template boost::any Dict::toany<double_tuple>(double_tuple arg) const;
}
