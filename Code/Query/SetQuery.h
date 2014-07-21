//
// Copyright (c) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_SETQUERY_H__
#define __RD_SETQUERY_H__
#include <set>
#include "Query.h"
#include <sstream>
#include <algorithm>
#include <iterator>

namespace Queries{
  //! \brief a Query implementing a set: arguments must 
  //!  one of a set of values
  //!
  template <class MatchFuncArgType, class DataFuncArgType=MatchFuncArgType,
    bool needsConversion=false>
  class SetQuery :
    public Query<MatchFuncArgType, DataFuncArgType,needsConversion> {

  public:
    typedef std::set<MatchFuncArgType> CONTAINER_TYPE;

    SetQuery() : Query<MatchFuncArgType,DataFuncArgType,needsConversion>() {};

    //! insert an entry into our \c set
    void insert(const MatchFuncArgType what){
      if(d_set.find(what) == this->d_set.end()) this->d_set.insert(what);
    }

    //! clears our \c set
    void clear(){
      this->d_set.clear();
    }

    bool Match(const DataFuncArgType what) const {
      MatchFuncArgType mfArg = this->TypeConvert(what,Int2Type<needsConversion>());
      return ( this->d_set.find(mfArg) != this->d_set.end() ) ^ this->getNegation();
    };

    Query<MatchFuncArgType,DataFuncArgType,needsConversion> *
    copy( ) const {
      SetQuery<MatchFuncArgType,DataFuncArgType,needsConversion> *res =
	new SetQuery<MatchFuncArgType,DataFuncArgType,needsConversion>();
      res->setDataFunc(this->d_dataFunc);
      typename std::set<MatchFuncArgType>::const_iterator i;
      for(i=this->d_set.begin();
	  i!=this->d_set.end();
	  ++i){
	res->insert(*i);
      }
      res->setNegation(this->getNegation());
      res->d_description = this->d_description;
      return res;
    };

    typename CONTAINER_TYPE::const_iterator beginSet() const {
      return d_set.begin();
    };
    typename CONTAINER_TYPE::const_iterator endSet() const {
      return d_set.end();
    };
    unsigned int size() const {
      return d_set.size();
    };

    std::string getFullDescription() const {
      std::ostringstream res;
      res<<this->getDescription()<<" val";
      if(this->getNegation()) res<<" not in ";
      else res<<" in (";
      std::copy(d_set.begin(),d_set.end(),std::ostream_iterator<MatchFuncArgType>(res,", "));
      res<<")";
      return res.str();
    }

    
  protected:
    CONTAINER_TYPE d_set;
  };

}
#endif
