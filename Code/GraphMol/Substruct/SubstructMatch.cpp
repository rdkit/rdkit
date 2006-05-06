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


#ifdef CACHE_ARMOLGRAPHS
#include <RDGeneral/Dict.h>
#endif

namespace RDKit{
  typedef ARGraph<const Atom, const Bond> AR_MOLGRAPH;

#ifdef CACHE_ARMOLGRAPHS
  // FIX: this stuff does not work
  void force_SubStructtypes(){
    boost::shared_ptr<AR_MOLGRAPH> molG_sptr;
    Dict tD;
    tD.fromany< boost::shared_ptr<AR_MOLGRAPH> >(boost::any(1));
    tD.toany< boost::shared_ptr<AR_MOLGRAPH> >(molG_sptr);
  }
#endif  

  AR_MOLGRAPH *getMolGraph(const ROMol &mol){
    AR_MOLGRAPH *molG;
#ifdef CACHE_ARMOLGRAPHS
    // FIX: this stuff does not work
    if(!mol.hasProp("SubstructGraphPtr")){
      ARGEdit mEd;
      MolToVFGraph(mol,&mEd);
      molG = new AR_MOLGRAPH(&mEd);
      boost::shared_ptr<AR_MOLGRAPH> molG_sptr(molG);
      boost::shared_ptr<void> void_sptr;
      void_sptr = reinterpret_cast< boost::shared_ptr<void> >(molG_sptr);
      mol.setProp("SubstructGraphPtr",void_sptr,true);
    } else {
      boost::shared_ptr<void> void_sptr;
      mol.getProp("SubstructGraphPtr",void_sptr);
      boost::shared_ptr<AR_MOLGRAPH> molG_sptr;
      molG_sptr = reinterpret_cast< boost::shared_ptr<AR_MOLGRAPH> >(void_sptr);
      molG = molG_sptr.get(); 
    }
#else
    ARGEdit mEd;
    MolToVFGraph(mol,&mEd);
    molG = new AR_MOLGRAPH(&mEd);
#endif
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
#ifndef CACHE_ARMOLGRAPHS
    delete queryG;
#endif

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

#ifndef CACHE_ARMOLGRAPHS
    delete queryG;
#endif

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
#ifndef CACHE_ARMOLGRAPHS
    delete queryG;
#endif
  
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
  
  ROMol *deleteSubstructs(const ROMol &mol, const ROMol &query,
			  bool onlyFrags) {
    RWMol *res = static_cast<RWMol*>(new ROMol(mol,false));
    std::vector<MatchVectType> fgpMatches;
    std::vector<MatchVectType>::const_iterator mati;
    std::pair<int, int> amat;
    VECT_INT_VECT matches; // all matches onto the molecule - list of list of atom ids
    MatchVectType::const_iterator mi;
    // do the substructure matching and get the atoms that match the query
    SubstructMatch(*res, query, fgpMatches);

    // if didn't find any matches nothing to be done here
    // simply return a copy of the molecule
    if (fgpMatches.size() == 0) {
      return res;
    }

    for (mati = fgpMatches.begin(); mati != fgpMatches.end(); mati++) {
      INT_VECT match; // each match onto the molecule - list of atoms ids
      for (mi = mati->begin(); mi != mati->end(); mi++) {
        match.push_back(mi->second);
      }
      matches.push_back(match);
    }

    // now loop over the list of matches and check if we can delete any of them
    INT_VECT delList;
    
    VECT_INT_VECT_I mxi, fi;
    if (onlyFrags) {
      VECT_INT_VECT frags;
      
      unsigned int nfrags = MolOps::getMolFrags(*res, frags);
      for (fi = frags.begin(); fi != frags.end(); fi++) {
        std::sort(fi->begin(), fi->end());
        for (mxi = matches.begin(); mxi != matches.end(); mxi++) {
          std::sort(mxi->begin(), mxi->end());
          if ((*fi) == (*mxi) ) {
            INT_VECT tmp; 
            Union((*mxi), delList, tmp);
            delList = tmp;
            break;
          } // end of if we found a matching fragment
        } // endof loop over matches
      } // end of loop over fragments
    } // end of if onlyFrags
    else {
      // in this case we want to delete any matches we find
      // simply loop over the matches and collect the atoms that need to 
      // be removes
      for (mxi = matches.begin(); mxi != matches.end(); mxi++) {
        INT_VECT tmp; 
        Union((*mxi), delList, tmp);
        delList = tmp;
      }
    }

    // now loop over the union list and delete the atoms
    // Will do this in the decreasing order of the atomIds
    // this is so that the AtomIds ids in the "delList" are
    // not invalidated by a previous removal (removing atom number i changes 
    // the atom indices only atoms with indices >i )
    std::sort(delList.begin(), delList.end());

    INT_VECT_RI dri;
    for (dri = delList.rbegin(); dri != delList.rend(); dri++) {
      res->removeAtom(*dri);
    }
    // if we removed any atoms, clear the computed properties:
    if(delList.size()){
      res->clearComputedProps(true);
      // update our properties, but allow unhappiness:
      res->updatePropertyCache(false);
    }
    return res;
  }

};
