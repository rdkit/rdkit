// $Id$
//
//  Copyright (C) 2001-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include "SubstructMatch.h"
#include "SubstructUtils.h"
#include <boost/smart_ptr.hpp>

#include <argedit.h>
#include <vf2_mono_state.h>
#include <match.h>

namespace RDKit{
  typedef ARGraph<const Atom, const Bond> AR_MOLGRAPH;


  AR_MOLGRAPH *getMolGraph(const ROMol &mol){
    AR_MOLGRAPH *molG;
    ARGEdit mEd;
    MolToVFGraph(mol,&mEd);
    molG = new AR_MOLGRAPH(&mEd);
    return molG;
  }

  void MatchSubqueries(AR_MOLGRAPH *molG,QueryAtom::QUERYATOM_QUERY *q);

  // ----------------------------------------------
  //
  // find one match
  //
  bool SubstructMatch(const ROMol &mol,const ROMol &query,MatchVectType &matchVect,
		      bool recursionPossible)
  {
    AR_MOLGRAPH *molG = getMolGraph(mol);
    bool res = SubstructMatch(molG,query,matchVect,recursionPossible);
    delete molG;
    return res;
  }

  bool SubstructMatch(AR_MOLGRAPH *molG,const ROMol &query,MatchVectType &matchVect,
		      bool recursionPossible){
    PRECONDITION(molG,"bad molecule");

    if(recursionPossible){
      ROMol::ConstAtomIterator atIt;
      for(atIt=query.beginAtoms();atIt!=query.endAtoms();atIt++){
	if((*atIt)->getQuery()){
	  MatchSubqueries(molG,(*atIt)->getQuery());
	}
      }
    }
    AR_MOLGRAPH *queryG = getMolGraph(query);
    queryG->SetNodeCompat(atomCompat);
    queryG->SetEdgeCompat(bondCompat);

    int nQueryAtoms;
    nQueryAtoms = query.getNumAtoms();
    node_id *ni1 = new node_id[nQueryAtoms];
    node_id *ni2 = new node_id[nQueryAtoms];
    MatcherState s0(queryG,molG);

    bool res;
    int n;
    matchVect.clear();
    matchVect.resize(0);
    if(match(&s0,&n,ni1,ni2)){
      matchVect.resize(nQueryAtoms);
      for(int i=0;i<nQueryAtoms;i++) matchVect[i]=std::pair<int,int>(ni1[i],ni2[i]);
      res = true;
    } else {
      res = false;
    }
    delete [] ni1;
    delete [] ni2;
    delete queryG;

    return res;
  }


  // ----------------------------------------------
  //
  // find all matches
  //
  //  NOTE: this blows out the contents of matches
  //
  unsigned int SubstructMatch(const ROMol &mol,const ROMol &query,
		     std::vector< MatchVectType > &matches,
		     bool uniquify,bool recursionPossible) {
    AR_MOLGRAPH *molG = getMolGraph(mol);
    unsigned int res = SubstructMatch(molG,query,matches,uniquify,recursionPossible);
    delete molG;
    return res;
  }
  unsigned int SubstructMatch(AR_MOLGRAPH *molG,const ROMol &query,
		     std::vector< MatchVectType > &matches,
		     bool uniquify,bool recursionPossible)
  {
    PRECONDITION(molG,"bad molecule pointer");

    if(recursionPossible){
      ROMol::ConstAtomIterator atIt;
      for(atIt=query.beginAtoms();atIt!=query.endAtoms();atIt++){
    	if((*atIt)->getQuery()){
    	  MatchSubqueries(molG,(*atIt)->getQuery());
    	}
      }
    }

    AR_MOLGRAPH *queryG = getMolGraph(query);
    queryG->SetNodeCompat(atomCompat);
    queryG->SetEdgeCompat(bondCompat);

    MatcherState s0(queryG,molG);

    matches.clear();
    matches.resize(0);
    unsigned int res;
    res = match(&s0,substructVisitor,(void *)&matches);

    if(res){
      if(uniquify){
	   removeDuplicates(matches);
	   res = matches.size();
      }
    } else {
      matches.clear();
      matches.resize(0);
    }

    delete queryG;

    return res;
  }


  // ----------------------------------------------
  //
  // Intended for internal use 
  //
  unsigned int RecursiveMatcher(AR_MOLGRAPH *molG,const ROMol &query,
		       std::vector< int > &matches)
  {
    PRECONDITION(molG,"bad molecule");
    //std::cout << " >>> RecursiveMatcher: " << int(query) << std::endl;
    ARGEdit qEd;

    ROMol::ConstAtomIterator atIt;
    for(atIt=query.beginAtoms();atIt!=query.endAtoms();atIt++){
      if((*atIt)->getQuery()){
    	MatchSubqueries(molG,(*atIt)->getQuery());
      }
    }
 
    AR_MOLGRAPH *queryG = getMolGraph(query);
    queryG->SetNodeCompat(atomCompat);
    queryG->SetEdgeCompat(bondCompat);

    MatcherState s0(queryG,molG);

    matches.clear();
    matches.resize(0);
    unsigned int res;
    res = match(&s0,substructHeadVisitor,(void *)&matches);

    if(!res){
      matches.clear();
      matches.resize(0);
    }
    //std::cout << " <<< RecursiveMatcher: " << int(query) << std::endl;
    delete queryG;
  
    return res;
  }


  void MatchSubqueries(AR_MOLGRAPH  *molG,QueryAtom::QUERYATOM_QUERY *query){
    PRECONDITION(molG,"bad molecule");
    PRECONDITION(query,"bad query");
    //std::cout << "*-*-* MS: " << (int)query << std::endl;
    //std::cout << "\t\t" << typeid(*query).name() << std::endl;
    if(query->getDescription()=="RecursiveStructure"){
      RecursiveStructureQuery *rsq=(RecursiveStructureQuery *)query;
      ROMol const *queryMol = rsq->getQueryMol();
      // in case we are reusing this query, clear its contents now.
      rsq->clear();
      if(queryMol){
    	std::vector< int > matchStarts;
    	unsigned int res = RecursiveMatcher(molG,*queryMol,matchStarts);
    	if(res){
    	  for(std::vector<int>::iterator i=matchStarts.begin();
    	      i!=matchStarts.end();
    	      i++){
    	    rsq->insert(*i);
    	  }
    	}
      }
    } else {
      //std::cout << "\tmsq1: "; 
    }
  
    // now recurse over our children (these things can be nested)
    Queries::Query<int,Atom const*,true>::CHILD_VECT_CI childIt;
    //std::cout << query << " " << query->endChildren()-query->beginChildren() <<  std::endl;
    for(childIt=query->beginChildren();childIt!=query->endChildren();childIt++){
      MatchSubqueries(molG,childIt->get());
    }
    //std::cout << "<<- back " << (int)query << std::endl;
  }

};
