//
//  Copyright (C) 2001-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_QUERYATOM_H_
#define _RD_QUERYATOM_H_

#include "Atom.h"
#include <Query/QueryObjects.h>
#include <GraphMol/QueryOps.h>

namespace RDKit{

  //! Class for storing atomic queries
  /*!
    QueryAtom objects are derived from Atom objects, so they can be
    added to molecules and the like, but they have much fancier
    querying capabilities.
    
   */
  class QueryAtom : public Atom {
  public:
    typedef Queries::Query<int,Atom const *,true> QUERYATOM_QUERY;

    QueryAtom() : Atom(), dp_query(NULL) {};
    explicit QueryAtom(int num) : Atom(num), dp_query(makeAtomNumQuery(num)) {};
    explicit QueryAtom(const Atom &other) : Atom(other), dp_query(makeAtomNumQuery(other.getAtomicNum())) {};
    QueryAtom( const QueryAtom & other) : Atom(other){
      dp_query  = other.dp_query->copy();
    };
    ~QueryAtom();


    //! returns a copy of this query, owned by the caller
    Atom *copy() const;

    // This method can be used to distinguish query atoms from standard atoms:
    bool hasQuery() const { return dp_query!=0; };
    
    //! replaces our current query with the value passed in
    void setQuery(QUERYATOM_QUERY *what) { dp_query = what; }
    //! returns our current query
    QUERYATOM_QUERY *getQuery() const { return dp_query; };

    //! expands our current query
    /*!
      \param what          the Queries::Query to be added
      \param how           the operator to be used in the expansion
      \param maintainOrder (optional) flags whether the relative order of
                           the queries needs to be maintained, if this is
			   false, the order is reversed
      <b>Notes:</b>
        - \c what should probably be constructed using one of the functions
	   defined in QueryOps.h
	- the \c maintainOrder option can be useful because the combination
	  operators short circuit when possible.
      
    */
    void expandQuery(QUERYATOM_QUERY *what,
		     Queries::CompositeQueryType how=Queries::COMPOSITE_AND,
		     bool maintainOrder=true);

    //! returns true if we match Atom \c what
    bool Match(const Atom::ATOM_SPTR &what) const;
    //! \overload
    bool Match(Atom const *what) const;
    
    //! returns true if our query details match those of QueryAtom \c what
    bool QueryMatch(QueryAtom const *what) const;

  private:
    QUERYATOM_QUERY *dp_query;
  
  }; // end o' class

}; // end o' namespace


#endif
