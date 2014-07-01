//
//  Copyright (C) 2001-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_QUERYBOND_H
#define _RD_QUERYBOND_H

#include <Query/QueryObjects.h>
#include "Bond.h"
#include "QueryOps.h"


namespace RDKit{

  //! Class for storing Bond queries
  /*!
    QueryBond objects are derived from Bond objects, so they can be
    added to molecules and the like, but they have much fancier
    querying capabilities.
    
   */

  class QueryBond : public Bond {
  public:
    typedef Queries::Query<int,Bond const *,true> QUERYBOND_QUERY;

    QueryBond() : Bond(), dp_query(NULL) {};
    //! initialize with a particular bond order
    explicit QueryBond(BondType bT);
    //! initialize from a bond
    explicit QueryBond(const Bond &other) : Bond(other), dp_query(makeBondOrderEqualsQuery(other.getBondType())) {};
    QueryBond(const QueryBond &other) : Bond(other), dp_query(other.dp_query->copy()) {};

    ~QueryBond();

    
    //! returns a copy of this query, owned by the caller
    virtual Bond *copy() const;

    QueryBond &operator=(const QueryBond &other);

    //! sets the BondType of this query:
    void setBondType(BondType bT);
    //! sets the BondDir of this query:
    void setBondDir(BondDir bD);
  

    //! returns true if we match Bond \c what
    bool Match(const Bond::BOND_SPTR what) const;
    //! \overload
    bool Match(Bond const *what) const;

    //! returns true if our query details match those of QueryBond \c what
    bool QueryMatch(QueryBond const *what) const;


    // This method can be used to distinguish query bonds from standard bonds
    bool hasQuery() const { return dp_query!=0; };
    
    //! returns our current query
    QUERYBOND_QUERY *getQuery() const { return dp_query; };
    //! replaces our current query with the value passed in
    void setQuery(QUERYBOND_QUERY *what) {
      // free up any existing query (Issue255):
      if(dp_query)delete dp_query;
      dp_query = what;
    };

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
    void expandQuery(QUERYBOND_QUERY *what,
		     Queries::CompositeQueryType how=Queries::COMPOSITE_AND,
		     bool maintainOrder=true);

  protected:
    QUERYBOND_QUERY *dp_query;
  };

};


#endif
