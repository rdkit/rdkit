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
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryOps.h>

namespace RDKit{

QueryAtom::~QueryAtom(){
  delete dp_query;
  dp_query=NULL;
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

bool QueryAtom::Match(const Atom::ATOM_SPTR &what) const {
  return Match(what.get());
}
  
  namespace {
    bool localMatch(ATOM_EQUALS_QUERY const *q1,ATOM_EQUALS_QUERY const *q2){
      if(q1->getNegation()==q2->getNegation()){
        return q1->getVal()==q2->getVal();
      } else {
        return q1->getVal()!=q2->getVal();
      }
    }
    bool queriesMatch(QueryAtom::QUERYATOM_QUERY const *q1,
                      QueryAtom::QUERYATOM_QUERY const *q2){
      PRECONDITION(q1,"no q1");
      PRECONDITION(q2,"no q2");

      static const unsigned int nQueries=18;
      static std::string equalityQueries[nQueries]={"AtomRingBondCount","AtomRingSize","AtomMinRingSize","AtomImplicitValence",\
                                                    "AtomExplicitValence","AtomTotalValence","AtomAtomicNum","AtomExplicitDegree",\
                                                    "AtomTotalDegree","AtomHCount","AtomIsAromatic","AtomIsAliphatic","AtomUnsaturated",\
                                                    "AtomMass","AtomFormalCharge","AtomHybridization","AtomInRing","AtomInNRings"};
                                                  
      bool res=false;
      std::string d1=q1->getDescription();
      std::string d2=q2->getDescription();
      if(d1=="AtomNull"||d2=="AtomNull"){
        res = true;
      } else if(d1=="AtomOr"){
        // FIX: handle negation on AtomOr and AtomAnd
        for(QueryAtom::QUERYATOM_QUERY::CHILD_VECT_CI iter1=q1->beginChildren();iter1!=q1->endChildren();++iter1){
          if(d2=="AtomOr"){
            for(QueryAtom::QUERYATOM_QUERY::CHILD_VECT_CI iter2=q2->beginChildren();iter2!=q2->endChildren();++iter2){
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
      } else if(d1=="AtomAnd"){
          res=true;
          for(QueryAtom::QUERYATOM_QUERY::CHILD_VECT_CI iter1=q1->beginChildren();iter1!=q1->endChildren();++iter1){
            bool matched=false;
            if(d2=="AtomAnd"){
              for(QueryAtom::QUERYATOM_QUERY::CHILD_VECT_CI iter2=q2->beginChildren();iter2!=q2->endChildren();++iter2){
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
      // FIX : handle AtomXOr
      } else if(d2=="AtomOr"){
        // FIX: handle negation on AtomOr and AtomAnd
        for(QueryAtom::QUERYATOM_QUERY::CHILD_VECT_CI iter2=q2->beginChildren();iter2!=q2->endChildren();++iter2){
          if(queriesMatch(q1,iter2->get())) {
            res=true;
            break;
          }
        }
      } else if(d2=="AtomAnd"){
        res=true;
        for(QueryAtom::QUERYATOM_QUERY::CHILD_VECT_CI iter2=q2->beginChildren();iter2!=q2->endChildren();++iter2){
          if(!queriesMatch(q1,iter2->get())){
            res=false;
            break;
          }
        }
      } else if(d1==d2){
          if(std::find(&equalityQueries[0],&equalityQueries[nQueries],d1)!= &equalityQueries[nQueries]){
            res=localMatch(static_cast<ATOM_EQUALS_QUERY const *>(q1),
                           static_cast<ATOM_EQUALS_QUERY const *>(q2));
          }
      } else {
      }
      return res;
    }
  } //end of local namespace

bool QueryAtom::Match(Atom const *what) const {
  PRECONDITION(what,"bad query atom");
  PRECONDITION(dp_query,"no query set");
  return dp_query->Match(what);
}
bool QueryAtom::QueryMatch(QueryAtom const *what) const {
  PRECONDITION(what,"bad query atom");
  PRECONDITION(dp_query,"no query set");
  if(!what->hasQuery()){
    return dp_query->Match(what);
  } else {
    return queriesMatch(dp_query,what->getQuery());
  }
}

}; // end o' namespace
