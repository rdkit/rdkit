// $Id: QueryBond.cpp 4961 2006-02-18 00:14:47Z glandrum $
//
//  Copyright (C) 2001-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <GraphMol/QueryBond.h>

namespace RDKit{

QueryBond::QueryBond(BondType bT) : Bond(bT) {
  if( bT != Bond::UNSPECIFIED ) dp_query = makeBondOrderEqualsQuery(bT);
  else dp_query = makeBondNullQuery();
};


QueryBond::~QueryBond(){
  if( dp_query ) {
    delete dp_query;
    dp_query=0;
  }
};

QueryBond &QueryBond::operator=(const QueryBond &other){
  // FIX: should we copy molecule ownership?  I don't think so.
  // FIX: how to deal with atom indices?
  dp_mol = 0;
  d_bondType = other.d_bondType;
  dp_query = other.dp_query->copy();
  if(dp_props) delete dp_props;
  if(other.dp_props)
    dp_props = new Dict(*other.dp_props);
  return *this;
}

Bond *QueryBond::copy() const {
  QueryBond *res = new QueryBond(*this);
  return res;
}


void QueryBond::setBondType(BondType bT) {
  // NOTE: calling this blows out any existing query
  d_bondType = bT;
  if(dp_query){
    delete dp_query;
    dp_query = 0;
  }
  dp_query = makeBondOrderEqualsQuery(bT);
}

void QueryBond::setBondDir(BondDir bD) {
  // NOTE: calling this blows out any existing query
  //
  //   Ignoring bond orders (which this implicitly does by blowing out
  //   any bond order query) is ok for organic molecules, where the
  //   only bonds assigned directions are single.  It'll fail in other
  //   situations, whatever those may be.
  //
  d_dirTag = bD;
  if(dp_query){
    delete dp_query;
    dp_query = 0;
  }
  dp_query = makeBondDirEqualsQuery(bD);
}


bool QueryBond::Match(const Bond::BOND_SPTR what) const {
  return Match(what.get());
}
  
bool QueryBond::Match(Bond const *what) const {
  PRECONDITION(dp_query,"no query set");
  return dp_query->Match(what);
}

void QueryBond::expandQuery(QUERYBOND_QUERY *what,
		 Queries::CompositeQueryType how,
		 bool maintainOrder){
  QUERYBOND_QUERY *origQ = dp_query;
  std::string descrip;
  switch(how) {
  case Queries::COMPOSITE_AND:
    dp_query = new BOND_AND_QUERY;
    descrip = "BondAnd";
    break;
  case Queries::COMPOSITE_OR:
    dp_query = new BOND_OR_QUERY;
    descrip = "BondOr";
    break;
  case Queries::COMPOSITE_XOR:
    dp_query = new BOND_XOR_QUERY;
    descrip = "BondXor";
    break;
  default:
    UNDER_CONSTRUCTION("unrecognized combination query");
  }
  dp_query->setDescription(descrip);
  if(maintainOrder){
    dp_query->addChild(QUERYBOND_QUERY::CHILD_TYPE(origQ));
    dp_query->addChild(QUERYBOND_QUERY::CHILD_TYPE(what));
  } else {
    dp_query->addChild(QUERYBOND_QUERY::CHILD_TYPE(what));
    dp_query->addChild(QUERYBOND_QUERY::CHILD_TYPE(origQ));
  }
}

} // end o' namespace
