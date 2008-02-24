//
// Copyright (c) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef __RD_SETQUERY_H__
#define __RD_SETQUERY_H__
#include <boost/serialization/set.hpp>

#include "Query.h"

namespace Queries{
  //! \brief a Query implementing a set: arguments must 
  //!  one of a set of values
  //!
  //!  There is also an optional tolerance to be used in comparisons
  template <class MatchFuncArgType, class DataFuncArgType=MatchFuncArgType,
    bool needsConversion=false>
  class SetQuery :
    public Query<MatchFuncArgType, DataFuncArgType,needsConversion> {

  public:
    SetQuery() : Query<MatchFuncArgType,DataFuncArgType,needsConversion>() {};

    //! insert an entry into our \c set
    void insert(const MatchFuncArgType what){
      //std::cout << "SET QUERY INSERT: " << what << std::endl;
      if(d_set.find(what) == this->d_set.end()) this->d_set.insert(what);
    }

    //! clears our \c set
    void clear(){
      //std::cout << "SET QUERY CLEAR " << std::endl;
      this->d_set.clear();
    }

    bool Match(const DataFuncArgType what) const {
      MatchFuncArgType mfArg = TypeConvert(what,Int2Type<needsConversion>());
      //std::cerr << "SET QUERY SEARCH: " << mfArg << ": "  << (d_set.find(mfArg)==d_set.end()) << std::endl;
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
  protected:
    std::set<MatchFuncArgType> d_set;
  private:    
    friend class boost::serialization::access;
    //! required by boost::serialization, handles both serialization and
    //! deserialization
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & boost::serialization::base_object<Query<MatchFuncArgType,DataFuncArgType,needsConversion> >(*this);
      ar & d_set;
    }

  };

}
#endif
