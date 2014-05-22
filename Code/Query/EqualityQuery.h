//
// Copyright (c) 2003-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_EQUALITYQUERY_H__
#define __RD_EQUALITYQUERY_H__
#include "Query.h"
#include <sstream>

namespace Queries {


  //! \brief a Query implementing ==: arguments must match a particular
  //!  value (within an optional tolerance)
  template <typename MatchFuncArgType, typename DataFuncArgType=MatchFuncArgType,
    bool needsConversion=false>
  class EqualityQuery : public Query<MatchFuncArgType, DataFuncArgType,needsConversion> {

  public:
    EqualityQuery() : d_tol(0) {
      this->df_negate=false;
    };

    //! constructs with our target value
    explicit EqualityQuery(MatchFuncArgType v) : d_val(v), d_tol(0) {
      this->df_negate=false;
    };


    //! constructs with our target value and a tolerance
    EqualityQuery(MatchFuncArgType v,MatchFuncArgType t) : d_val(v), d_tol(t) {
      this->df_negate=false;
    };


    //! sets our target value
    void setVal(MatchFuncArgType what) { this->d_val = what; };
    //! returns our target value
    const MatchFuncArgType getVal() const { return this->d_val; };

    //! sets our tolerance
    void setTol(MatchFuncArgType what) { this->d_tol = what; };
    //! returns out tolerance
    const MatchFuncArgType getTol() const { return this->d_tol; };
  
    virtual bool Match(const DataFuncArgType what) const {
      MatchFuncArgType mfArg = this->TypeConvert(what,Int2Type<needsConversion>());
      if( queryCmp(this->d_val,mfArg,this->d_tol) == 0 ){
	if( this->getNegation() ){
	  return false;
	}
	else{
	  return true;
	}
      } else {
	if( this->getNegation() ){
	  return true;
	}
	else{
	  return false;
	}
      }
    };

    virtual Query<MatchFuncArgType,DataFuncArgType,needsConversion> *
    copy( ) const {
      EqualityQuery<MatchFuncArgType,DataFuncArgType,needsConversion> *res =
	new EqualityQuery<MatchFuncArgType,DataFuncArgType,needsConversion>();
      res->setNegation(this->getNegation());
      res->setVal(this->d_val);
      res->setTol(this->d_tol);
      res->setDataFunc(this->d_dataFunc);
      res->d_description = this->d_description;
      return res;
    };

    std::string getFullDescription() const {
      std::ostringstream res;
      res<<this->getDescription();
      res<<" "<<this->d_val;
      if(this->getNegation()) res<<" != ";
      else res<<" = ";
      res<<"val";
      return res.str();
    }

  protected:
    MatchFuncArgType d_val;
    MatchFuncArgType d_tol;
  };

}
#endif
