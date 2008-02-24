//
// Copyright (c) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef __RD_XORQUERY_H__
#define __RD_XORQUERY_H__

#include "Query.h"

namespace Queries{

  //! a Query implementing XOR: requires exactly one child to be \c true
  template <class MatchFuncArgType, class DataFuncArgType=MatchFuncArgType,
    bool needsConversion=false>
  class XOrQuery : public Query<MatchFuncArgType, DataFuncArgType,needsConversion> {

  public:
    typedef Query<MatchFuncArgType, DataFuncArgType,needsConversion> BASE;
    XOrQuery() {
      this->df_negate = false;
    };

    bool Match(const DataFuncArgType what) const {
      bool res = false;
      typename BASE::CHILD_VECT_CI it1;
      for(it1=this->beginChildren();
	  it1!=this->endChildren();
	  ++it1){
	bool tmp = (*it1)->Match(what);
	if( tmp ){
	  if( res ){
	    res = false;
	    break;
	  } else{
	    res = true;
	  }
	}
      }
      if( this->getNegation() ) res = !res;
      return res;
    };

    Query<MatchFuncArgType,DataFuncArgType,needsConversion> *
    copy( ) const {
      XOrQuery<MatchFuncArgType,DataFuncArgType,needsConversion> *res =
	new XOrQuery<MatchFuncArgType,DataFuncArgType,needsConversion>();

      typename BASE::CHILD_VECT_CI i;
      for(i=this->beginChildren();
	  i!=this->endChildren();
	  ++i){
	res->addChild(*i);
      }
      res->setNegation(this->getNegation());
      res->d_description = this->d_description;
      return res;
    };
  private:    
    friend class boost::serialization::access;
    //! required by boost::serialization, handles both serialization and
    //! deserialization
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & boost::serialization::base_object<Query<MatchFuncArgType,DataFuncArgType,needsConversion> >(*this);
    }

  };

}
#endif
