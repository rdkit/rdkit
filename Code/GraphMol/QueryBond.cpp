// $Id$
//
//  Copyright (C) 2001-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/QueryBond.h>

namespace RDKit{

QueryBond::QueryBond(BondType bT) : Bond(bT) {
  if( bT != Bond::UNSPECIFIED ) dp_query = makeBondOrderEqualsQuery(bT);
  else dp_query = makeBondNullQuery();
};


QueryBond::~QueryBond(){
  delete dp_query;
  dp_query=NULL;
};

QueryBond &QueryBond::operator=(const QueryBond &other){
  // FIX: should we copy molecule ownership?  I don't think so.
  // FIX: how to deal with atom indices?
  dp_mol = 0;
  d_bondType = other.d_bondType;
  dp_query = other.dp_query->copy();
  delete dp_props;
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
  delete dp_query;
  dp_query = NULL;

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
#if 0  
  delete dp_query;
  dp_query = NULL;
  dp_query = makeBondDirEqualsQuery(bD);
#endif
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

bool QueryBond::Match(const Bond::BOND_SPTR what) const {
  return Match(what.get());
}
  

  namespace {
    bool localMatch(BOND_EQUALS_QUERY const *q1,BOND_EQUALS_QUERY const *q2){
      if(q1->getNegation()==q2->getNegation()){
        return q1->getVal()==q2->getVal();
      } else {
        return q1->getVal()!=q2->getVal();
      }
    }

    bool queriesMatch(QueryBond::QUERYBOND_QUERY const *q1,
                      QueryBond::QUERYBOND_QUERY const *q2){
      PRECONDITION(q1,"no q1");
      PRECONDITION(q2,"no q2");

      static const unsigned int nQueries=6;
      static std::string equalityQueries[nQueries]={"BondRingSize","BondMinRingSize","BondOrder","BondDir",\
                                                    "BondInRing","BondInNRings"};
                                                  
      bool res=false;
      std::string d1=q1->getDescription();
      std::string d2=q2->getDescription();
      if(d1=="BondNull"||d2=="BondNull"){
        res = true;
      } else if(d1=="BondOr"){
        // FIX: handle negation on BondOr and BondAnd
        for(QueryBond::QUERYBOND_QUERY::CHILD_VECT_CI iter1=q1->beginChildren();iter1!=q1->endChildren();++iter1){
          if(d2=="BondOr"){
            for(QueryBond::QUERYBOND_QUERY::CHILD_VECT_CI iter2=q2->beginChildren();iter2!=q2->endChildren();++iter2){
              if(queriesMatch(iter1->get(),iter2->get())){
                res=true;
                break;
              }
            }
          } else {
            if(queriesMatch(iter1->get(),q2)) {
              res=true;
            }
          }
          if(res) break;
        }
      } else if(d1=="BondAnd"){
          res=true;
          for(QueryBond::QUERYBOND_QUERY::CHILD_VECT_CI iter1=q1->beginChildren();iter1!=q1->endChildren();++iter1){
            bool matched=false;
            if(d2=="BondAnd"){
              for(QueryBond::QUERYBOND_QUERY::CHILD_VECT_CI iter2=q2->beginChildren();iter2!=q2->endChildren();++iter2){
                if(queriesMatch(iter1->get(),iter2->get())){
                  matched=true;
                  break;
                }
              }
            } else {
              matched = queriesMatch(iter1->get(),q2);
            }
            if(!matched){
              res=false;
              break;
            }
          }
      // FIX : handle BondXOr
      } else if(d2=="BondOr"){
        // FIX: handle negation on BondOr and BondAnd
        for(QueryBond::QUERYBOND_QUERY::CHILD_VECT_CI iter2=q2->beginChildren();iter2!=q2->endChildren();++iter2){
          if(queriesMatch(q1,iter2->get())) {
            res=true;
            break;
          }
        }
      } else if(d2=="BondAnd"){
        res=true;
        for(QueryBond::QUERYBOND_QUERY::CHILD_VECT_CI iter2=q2->beginChildren();iter2!=q2->endChildren();++iter2){
          if(queriesMatch(q1,iter2->get())){
            res=false;
            break;
          }
        }
      } else if(std::find(&equalityQueries[0],&equalityQueries[nQueries],d1)!= &equalityQueries[nQueries]){
          res=localMatch(static_cast<BOND_EQUALS_QUERY const *>(q1),
                         static_cast<BOND_EQUALS_QUERY const *>(q2));
      } 
      return res;
    }
  } //end of local namespace

bool QueryBond::Match(Bond const *what) const {
  PRECONDITION(what,"bad query bond");
  PRECONDITION(dp_query,"no query set");
  return dp_query->Match(what);
}
bool QueryBond::QueryMatch(QueryBond const *what) const {
  PRECONDITION(what,"bad query bond");
  PRECONDITION(dp_query,"no query set");
  if(!what->hasQuery()){
    return dp_query->Match(what);
  } else {
    return queriesMatch(dp_query,what->getQuery());
  }
}




} // end o' namespace
