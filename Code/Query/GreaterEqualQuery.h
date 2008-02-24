//
// Copyright (c) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef __RD_GREATEREQUALQUERY_H__
#define __RD_GREATEREQUALQUERY_H__
#include "Query.h"
#include "EqualityQuery.h"

namespace Queries {
  //! \brief a Query implementing >= using a particular
  //!  value (and an optional tolerance)
  template <class MatchFuncArgType, class DataFuncArgType=MatchFuncArgType,
    bool needsConversion=false>
  class GreaterEqualQuery :
    public EqualityQuery<MatchFuncArgType, DataFuncArgType,needsConversion> {

  public:
    GreaterEqualQuery() { this->d_tol=0;};
    //! constructs with our target value
    explicit GreaterEqualQuery(DataFuncArgType what) {
      this->d_val = what;
      this->d_tol = 0;
      this->df_negate = false;
    };
    //! constructs with our target value and a tolerance
    GreaterEqualQuery(DataFuncArgType v,DataFuncArgType t) {
      this->d_val = v;
      this->d_tol = t;
      this->df_negate = false;
    };


    bool Match(const DataFuncArgType what) const {
      MatchFuncArgType mfArg = TypeConvert(what,Int2Type<needsConversion>());
      if( queryCmp(this->d_val,mfArg,this->d_tol) >= 0 ){
	if( this->getNegation() ) return false;
	else return true;
      } else {
	if( this->getNegation() ) return true;
	else return false;
      }
    };
    Query<MatchFuncArgType,DataFuncArgType,needsConversion> *
    copy( ) const {
      GreaterEqualQuery<MatchFuncArgType,DataFuncArgType,needsConversion> *res =
	new GreaterEqualQuery<MatchFuncArgType,DataFuncArgType,needsConversion>();
      res->setVal(this->d_val);
      res->setTol(this->d_tol);
      res->setNegation(this->getNegation());
      res->setDataFunc(this->d_dataFunc);
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
      ar & boost::serialization::base_object<EqualityQuery<MatchFuncArgType,DataFuncArgType,needsConversion> >(*this);
    }

  };

}
#endif
