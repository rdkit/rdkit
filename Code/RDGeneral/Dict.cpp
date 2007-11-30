// $Id$
//
// Copyright (C) 2003-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
#include "Dict.h"


#include <boost/lexical_cast.hpp>
#include <boost/shared_array.hpp>
#include <boost/tuple/tuple.hpp>
#include <vector>
#include <list>

namespace RDKit{
  Dict::Dict() {
    _data.clear();
  };
  Dict::Dict(const Dict &other) {
    _data = other._data;
  };
  Dict &Dict::operator=(const Dict &other) {
    _data = other._data;
    return *this;
  };
    
  bool Dict::hasVal(const char *what) const {
    std::string key(what);
    return hasVal(key);
  };
  bool Dict::hasVal(const std::string &what) const {
    return _data.count(what)>0;
  };
    

  STR_VECT Dict::keys() const {
    STR_VECT res;
    DataType::const_iterator item;
    for (item = _data.begin(); item != _data.end(); item++) {
      res.push_back(item->first);
    }
    return res;
  }

  //
  //  We're going to try and be somewhat crafty about this getVal stuff to make these
  //  containers a little bit more generic.  The normal behavior here is that the
  //  value being queried must be directly castable to type T.  We'll robustify a
  //  little bit by trying that and, if the cast fails, attempting a couple of 
  //  other casts, which will then be lexically cast to type T.
  //
  void Dict::getVal(const std::string &what, std::string &res) const {
    if(! this->hasVal(what) ) throw KeyErrorException(what);
    const boost::any &val = _data.find(what)->second;
    try{
      res = boost::any_cast<std::string>(val);
    } catch (const boost::bad_any_cast &) {
      if(val.type()==typeid(int)){
	res = boost::lexical_cast<std::string>(boost::any_cast<int>(val));
      } else if(val.type()==typeid(long)){
	res = boost::lexical_cast<std::string>(boost::any_cast<long>(val));
      } else if(val.type()==typeid(float)){
	res = boost::lexical_cast<std::string>(boost::any_cast<float>(val));
      } else if(val.type()==typeid(double)){
	res = boost::lexical_cast<std::string>(boost::any_cast<double>(val));
      } else if(val.type()==typeid(const char *)){
	res = std::string(boost::any_cast<const char *>(val));
      } else {
	throw;
      }
    }
  };

  void Dict::clearVal(const std::string &what) {
    if(! this->hasVal(what) ) throw KeyErrorException(what);
    _data.erase(what);
  };

  void Dict::reset() {
    _data.clear();
  };
  
  template <typename T>
  T Dict::fromany(const boost::any &arg) const {
    return boost::any_cast<T>(arg);
  };

  template <typename T>
  boost::any Dict::toany(T arg) const {
    return boost::any(arg);
  };

  
#define ANY_FORCE(T) {tD.fromany< T >(boost::any(1));tD.toany< T >( T() );}

  void force_types(){
    Dict tD;
    bool fooBool = tD.fromany<bool>(boost::any(1));
    tD.toany<bool>(false);

    int fooInt = tD.fromany<int>(boost::any(1));
    fooInt += 1;
    tD.toany<int>(1);

    unsigned int fooUnsigned = tD.fromany<unsigned int>(boost::any(1));
    fooUnsigned += 1;
    tD.toany<unsigned int>(1);

    double fooDouble = tD.fromany<double>(boost::any(1));
    tD.toany<double>(1.0);

    std::string fooString = tD.fromany<std::string>(boost::any(std::string("1")));
    tD.toany<std::string>(std::string("1"));



    ANY_FORCE(std::vector<int>);
    ANY_FORCE(std::vector<unsigned int>);
    ANY_FORCE(std::vector<unsigned long long>);
    ANY_FORCE(std::vector<double>);
    ANY_FORCE(std::vector<std::string>);

    ANY_FORCE(std::vector< std::vector<int> >);
    ANY_FORCE(std::vector< std::vector<double> >);

    ANY_FORCE(boost::shared_array<double>);
    ANY_FORCE(boost::shared_array<int>);

    ANY_FORCE(std::list<int>);

    
    // FIX: it's UGLY that we have to include things like this:
    //ANY_FORCE( boost::tuples::tuple<double,double,double> );
    tD.fromany< boost::tuples::tuple<double,double,double> >(boost::any(1));
    tD.toany< boost::tuples::tuple<double,double,double> >(boost::tuples::tuple<double,double,double>());

  }





}
