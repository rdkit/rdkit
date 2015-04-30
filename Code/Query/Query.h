//
// Copyright (c) 2003-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_QUERY_H__
#define __RD_QUERY_H__

#ifdef _MSC_VER
#pragma warning (disable: 4800) // warning: converting things to bool
#endif

#include <vector>
#include <string>
#include <boost/smart_ptr.hpp>
#include <RDGeneral/Invariant.h>

namespace Queries {

  //! class to allow integer values to pick templates
  template <int v> 
  class Int2Type
  {
    enum { value = v };
  };

  //! Base class for all queries
  /*!
    Query objects have one or two functions associated with them:
      - <tt>bool matchFunc(MatchFuncArgType other)</tt> returns true or false
        to indicate whether this query matches \c other.
	This is mandatory.

      - <tt>MatchFuncArgType dataFunc(DataFuncArgType other)</tt> converts
        the argument \c other from \c DataFuncArgType to \c MatchFuncArgType.
	This is optional if \c DataFuncArgType is the same as (or implicitly
	convertible to) \c MatchFuncArgType.
	
  */
  template <class MatchFuncArgType, class DataFuncArgType=MatchFuncArgType,
    bool needsConversion=false>
  class Query {
  public:
    typedef boost::shared_ptr< Query<MatchFuncArgType, DataFuncArgType, needsConversion> > CHILD_TYPE;
    typedef std::vector< CHILD_TYPE > CHILD_VECT;
    typedef typename CHILD_VECT::iterator CHILD_VECT_I;
    typedef typename CHILD_VECT::const_iterator CHILD_VECT_CI;
  
    Query() : d_description(""),df_negate(false),d_matchFunc(NULL),d_dataFunc(NULL){};
    virtual ~Query() { this->d_children.clear(); };

  
    //! sets whether or not we are negated
    void setNegation(bool what) { this->df_negate = what; };
    //! returns whether or not we are negated
    bool getNegation() const { return this->df_negate; };

    //! sets our text description
    void setDescription(const std::string &descr) { this->d_description = descr; };
    //! \overload
    void setDescription(const char *descr) { this->d_description = std::string(descr); };
    //! returns our text description
    const std::string &getDescription() const { return this->d_description; };
    //! returns a fuller text description
    virtual std::string getFullDescription() const {
      if(!getNegation())
        return getDescription();
      else
        return "not "+getDescription();
    }
    
    //! sets our match function
    void setMatchFunc(bool (*what)(MatchFuncArgType)) { this->d_matchFunc = what; };
    //! returns our match function:
    bool (*getMatchFunc() const)(MatchFuncArgType) { return this->d_matchFunc; };
    //! sets our data function
    void setDataFunc(MatchFuncArgType (*what)(DataFuncArgType)) { this->d_dataFunc = what; };
    //! returns our data function:
    MatchFuncArgType (*getDataFunc() const)(DataFuncArgType) { return this->d_dataFunc; };

    //! adds a child to our list of children
    void addChild(CHILD_TYPE child) { this->d_children.push_back(child); };
    //! returns an iterator for the beginning of our child vector
    CHILD_VECT_CI beginChildren() const { return this->d_children.begin(); }
    //! returns an iterator for the end of our child vector
    CHILD_VECT_CI endChildren() const { return this->d_children.end(); }
  
    //! returns whether or not we match the argument
    virtual bool Match(const DataFuncArgType arg) const{
      MatchFuncArgType mfArg = TypeConvert(arg,Int2Type<needsConversion>());
      bool tRes;
      if(this->d_matchFunc) tRes = this->d_matchFunc(mfArg);
      else tRes = static_cast<bool>(mfArg);

      if( this->getNegation() ) return !tRes;
      else return tRes;
    };

    //! returns a copy of this Query
    /*!
       <b>Notes:</b>
         - the caller is responsible for <tt>delete</tt>ing the result
     */
    virtual Query<MatchFuncArgType,DataFuncArgType,needsConversion> *
    copy( ) const {
      Query<MatchFuncArgType,DataFuncArgType,needsConversion> *res =
	new Query<MatchFuncArgType,DataFuncArgType,needsConversion>();
      typename Query<MatchFuncArgType,DataFuncArgType,needsConversion>::CHILD_VECT_CI iter;
      for(iter=this->beginChildren();
	  iter!=this->endChildren();
	  ++iter){
        res->addChild(CHILD_TYPE(iter->get()->copy()));
      }
      res->df_negate = this->df_negate;
      res->d_matchFunc = this->d_matchFunc;
      res->d_dataFunc = this->d_dataFunc;
      res->d_description = this->d_description;
      return res;
    };

  protected :
    std::string d_description;
    CHILD_VECT d_children;
    bool df_negate;
    bool (*d_matchFunc) (MatchFuncArgType);
    MatchFuncArgType (*d_dataFunc)(DataFuncArgType);

    //! \brief calls our \c dataFunc (if it's set) on \c what and returns
    //! the result, otherwise returns \c what
    MatchFuncArgType TypeConvert(MatchFuncArgType what,Int2Type<false> /*d*/) const{
      MatchFuncArgType mfArg;
      if( this->d_dataFunc != NULL ){
	mfArg = this->d_dataFunc(what);
      } else {
	mfArg = what;
      }
      return mfArg;
    }
    //! calls our \c dataFunc (which must be set) on \c what and returns the result
    MatchFuncArgType TypeConvert(DataFuncArgType what,Int2Type<true> /*d*/) const{
      PRECONDITION(this->d_dataFunc,"no data function");
      MatchFuncArgType mfArg;
      mfArg = this->d_dataFunc(what);
      return mfArg;
    }
    

  };


  //----------------------------
  //
  // Used within query functions to compare values
  //
  //----------------------------
  template <class T1, class T2>
  int queryCmp(const T1 v1, const T2 v2, const T1 tol) {
    T1 diff = v1 - v2;
    if( diff <= tol ){
      if( diff >= -tol ){
	return 0;
      } else {
	return -1;
      }
    } else {
      return 1;
    }
  };


}
#endif
