// $Id$ 
//
//  Copyright (C) 2001-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryOps.h>

namespace RDKit{

QueryAtom::~QueryAtom(){
  if( dp_query ){
    delete dp_query;
    dp_query=0;
  }
};

Atom *QueryAtom::copy() const {
  QueryAtom *res = new QueryAtom(*this);
  return static_cast<Atom *>(res);
}
  
void QueryAtom::expandQuery(QUERYATOM_QUERY *what,
			    Queries::CompositeQueryType how,
			    bool maintainOrder) {
  QUERYATOM_QUERY *origQ = dp_query;
  std::string descrip;
  switch(how) {
  case Queries::COMPOSITE_AND:
    dp_query = new ATOM_AND_QUERY;
    descrip = "AtomAnd";
    break;
  case Queries::COMPOSITE_OR:
    dp_query = new ATOM_OR_QUERY;
    descrip = "AtomOr";
    break;
  case Queries::COMPOSITE_XOR:
    dp_query = new ATOM_XOR_QUERY;
    descrip = "AtomXor";
    break;
  default:
    UNDER_CONSTRUCTION("unrecognized combination query");
  }
  dp_query->setDescription(descrip);
  if(maintainOrder){
    dp_query->addChild(QUERYATOM_QUERY::CHILD_TYPE(origQ));
    dp_query->addChild(QUERYATOM_QUERY::CHILD_TYPE(what));
  } else {
    dp_query->addChild(QUERYATOM_QUERY::CHILD_TYPE(what));
    dp_query->addChild(QUERYATOM_QUERY::CHILD_TYPE(origQ));
  }
}

// FIX: the const crap here is all mucked up.
bool QueryAtom::Match(const Atom::ATOM_SPTR what) const {
  return Match(what.get());
}
  
bool QueryAtom::Match(Atom const *what) const {
  PRECONDITION(what,"bad query atom");
  PRECONDITION(dp_query,"no query set");
  return dp_query->Match(what);
}



}; // end o' namespace
