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

#ifdef USE_VFLIB
#include <argedit.h>
#include <vf2_mono_state.h>
#include <match.h>

namespace RDKit{
  typedef ARGraph<const Atom, const Bond> AR_MOLGRAPH;

#ifdef CACHE_ARMOLGRAPHS
  // FIX: given that it's using a global and creates a core leak this is no 
  //      permanent solution. 
namespace SubstructLocal {
  std::vector<AR_MOLGRAPH *> molGraphCache;
  void clearMolGraphCache(){
    for(std::vector<AR_MOLGRAPH *>::iterator it=molGraphCache.begin();
        it!=molGraphCache.end();++it){
      delete *it;
    }
    molGraphCache.clear(); 
  }
}
#endif  

  AR_MOLGRAPH *getMolGraph(const ROMol &mol,bool registerIt){
    AR_MOLGRAPH *molG=0;
#ifdef CACHE_ARMOLGRAPHS
    if(mol.hasProp("_SubstructGraphPtr")){
      unsigned int idx;
      mol.getProp("_SubstructGraphPtr",idx);
      if(idx<SubstructLocal::molGraphCache.size()){
        molG = SubstructLocal::molGraphCache[idx];
      }
    } 
    if(!molG){
      ARGEdit mEd;
      MolToVFGraph(mol,&mEd);
      molG = new AR_MOLGRAPH(&mEd);
      if(registerIt){
        unsigned int idx = SubstructLocal::molGraphCache.size();
        mol.setProp("_SubstructGraphPtr",idx,true);
        SubstructLocal::molGraphCache.push_back(molG);
      }
    }
#else
    ARGEdit mEd;
    MolToVFGraph(mol,&mEd);
    molG = new AR_MOLGRAPH(&mEd);
#endif
    return molG;
  }

  void MatchSubqueries(AR_MOLGRAPH *molG,QueryAtom::QUERYATOM_QUERY *q,bool useChirality,bool registerQuery);

  // ----------------------------------------------
  //
  // find one match
  //
  bool SubstructMatch(const ROMol &mol,const ROMol &query,MatchVectType &matchVect,
                      bool recursionPossible,bool useChirality,bool registerQuery)
  {
    AR_MOLGRAPH *molG = getMolGraph(mol);
    bool res = SubstructMatch(molG,query,matchVect,recursionPossible,useChirality,registerQuery);
#ifndef CACHE_ARMOLGRAPHS
    delete molG;
#else
    if(std::find(SubstructLocal::molGraphCache.begin(),
                 SubstructLocal::molGraphCache.end(),molG)==
       SubstructLocal::molGraphCache.end())
      delete molG;
#endif
    return res;
  }

  bool SubstructMatch(AR_MOLGRAPH *molG,const ROMol &query,MatchVectType &matchVect,
                      bool recursionPossible,bool useChirality,bool registerQuery){
    PRECONDITION(molG,"bad molecule");

    if(recursionPossible){
      ROMol::ConstAtomIterator atIt;
      for(atIt=query.beginAtoms();atIt!=query.endAtoms();atIt++){
        if((*atIt)->getQuery()){
          MatchSubqueries(molG,(*atIt)->getQuery(),useChirality,registerQuery);
        }
      }
    }
    AR_MOLGRAPH *queryG = getMolGraph(query,registerQuery);
    if(!useChirality){
      queryG->SetNodeCompat(atomCompat);
    } else {
      queryG->SetNodeCompat(chiralAtomCompat);
    }
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
#else
    if(std::find(SubstructLocal::molGraphCache.begin(),
                 SubstructLocal::molGraphCache.end(),queryG)==
       SubstructLocal::molGraphCache.end())
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
			      bool uniquify,bool recursionPossible,
			      bool useChirality,bool registerQuery) {
    AR_MOLGRAPH *molG = getMolGraph(mol);
    unsigned int res = SubstructMatch(molG,query,matches,uniquify,recursionPossible,
                                      useChirality,registerQuery);
#ifndef CACHE_ARMOLGRAPHS
    delete molG;
#else
    if(std::find(SubstructLocal::molGraphCache.begin(),
                 SubstructLocal::molGraphCache.end(),molG)==
       SubstructLocal::molGraphCache.end())
      delete molG;
#endif
    return res;
  }
  unsigned int SubstructMatch(AR_MOLGRAPH *molG,const ROMol &query,
			      std::vector< MatchVectType > &matches,
			      bool uniquify,bool recursionPossible,
			      bool useChirality,bool registerQuery)
  {
    PRECONDITION(molG,"bad molecule pointer");

    if(recursionPossible){
      ROMol::ConstAtomIterator atIt;
      for(atIt=query.beginAtoms();atIt!=query.endAtoms();atIt++){
        if((*atIt)->getQuery()){
          MatchSubqueries(molG,(*atIt)->getQuery(),useChirality,registerQuery);
        }
      }
    }

    AR_MOLGRAPH *queryG = getMolGraph(query,registerQuery);
    if(!useChirality){
      queryG->SetNodeCompat(atomCompat);
    } else {
      queryG->SetNodeCompat(chiralAtomCompat);
    }
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
#else
    if(std::find(SubstructLocal::molGraphCache.begin(),
                 SubstructLocal::molGraphCache.end(),queryG)==
       SubstructLocal::molGraphCache.end())
      delete queryG;
#endif

    return res;
  }


  // ----------------------------------------------
  //
  // Intended for internal use 
  //
  unsigned int RecursiveMatcher(AR_MOLGRAPH *molG,const ROMol &query,
				std::vector< int > &matches,bool useChirality,
                                bool registerQuery)
  {
    PRECONDITION(molG,"bad molecule");
    //std::cout << " >>> RecursiveMatcher: " << int(query) << std::endl;
    ARGEdit qEd;

    ROMol::ConstAtomIterator atIt;
    for(atIt=query.beginAtoms();atIt!=query.endAtoms();atIt++){
      if((*atIt)->getQuery()){
        MatchSubqueries(molG,(*atIt)->getQuery(),useChirality,registerQuery);
      }
    }
 
    AR_MOLGRAPH *queryG = getMolGraph(query,registerQuery);
    if(!useChirality){
      queryG->SetNodeCompat(atomCompat);
    } else {
      queryG->SetNodeCompat(chiralAtomCompat);
    }
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
#else
    if(std::find(SubstructLocal::molGraphCache.begin(),
                 SubstructLocal::molGraphCache.end(),queryG)==
       SubstructLocal::molGraphCache.end())
      delete queryG;
#endif
  
    return res;
  }


  void MatchSubqueries(AR_MOLGRAPH  *molG,QueryAtom::QUERYATOM_QUERY *query,bool useChirality,
                       bool registerQuery){
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
        unsigned int res = RecursiveMatcher(molG,*queryMol,matchStarts,useChirality,registerQuery);
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
      MatchSubqueries(molG,childIt->get(),useChirality,registerQuery);
    }
    //std::cout << "<<- back " << (int)query << std::endl;
  }

};

#else
#include "ullmann.hpp"

namespace RDKit{
  void MatchSubqueries(const ROMol &mol,QueryAtom::QUERYATOM_QUERY *q,bool useChirality,bool registerQuery);
  namespace detail {
    typedef std::list<std::pair<MolGraph::vertex_descriptor,MolGraph::vertex_descriptor> > ssPairType;

    class AtomLabelFunctor{
    public:
      AtomLabelFunctor(const ROMol &query,const ROMol &mol, bool useChirality) :
        d_query(query), d_mol(mol), df_useChirality(useChirality) {};
      bool operator()(unsigned int i,unsigned int j){
        bool res=false;
        if(!df_useChirality){
          res=atomCompat(d_query[i],d_mol[j]);
        } else {
          res=chiralAtomCompat(d_query[i],d_mol[j]);
        }
        //std::cerr<<" alf: "<<i<<" - "<<j<<"? "<<res<<std::endl;
        return res;
      }
    private:
      const ROMol &d_query;
      const ROMol &d_mol;
      bool df_useChirality;
    };
    class BondLabelFunctor{
    public:
      BondLabelFunctor(const ROMol &query,const ROMol &mol) :
        d_query(query), d_mol(mol) {};
      bool operator()(MolGraph::edge_descriptor i,MolGraph::edge_descriptor j){
        bool res=bondCompat(d_query[i],d_mol[j]);
        return res;
      }
    private:
      const ROMol &d_query;
      const ROMol &d_mol;
    };
  }    
  
  // ----------------------------------------------
  //
  // find one match
  //
  bool SubstructMatch(const ROMol &mol,const ROMol &query,MatchVectType &matchVect,
                      bool recursionPossible,bool useChirality,bool registerQuery)
  {

    //std::cerr<<"begin match"<<std::endl;
    if(recursionPossible){
      ROMol::ConstAtomIterator atIt;
      for(atIt=query.beginAtoms();atIt!=query.endAtoms();atIt++){
        if((*atIt)->getQuery()){
          MatchSubqueries(mol,(*atIt)->getQuery(),useChirality,registerQuery);
        }
      }
    }
    //std::cerr<<"main matching"<<std::endl;

    matchVect.clear();
    matchVect.resize(0);

    
    detail::AtomLabelFunctor atomLabeler(query,mol,useChirality);
    detail::BondLabelFunctor bondLabeler(query,mol);

    detail::ssPairType match;
    bool res=boost::ullmann(query.getTopology(),mol.getTopology(),
                            atomLabeler,bondLabeler,match);
    if(res){
     matchVect.resize(query.getNumAtoms());
     for(detail::ssPairType::const_iterator iter=match.begin();iter!=match.end();++iter){
       matchVect[iter->first]=std::pair<int,int>(iter->first,iter->second);
     }
    }    

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
			      bool uniquify,bool recursionPossible,
			      bool useChirality,bool registerQuery) {

    if(recursionPossible){
      ROMol::ConstAtomIterator atIt;
      for(atIt=query.beginAtoms();atIt!=query.endAtoms();atIt++){
        if((*atIt)->getQuery()){
          //std::cerr<<"recurse from atom "<<(*atIt)->getIdx()<<std::endl;
          MatchSubqueries(mol,(*atIt)->getQuery(),useChirality,registerQuery);
        }
      }
    }

    matches.clear();
    matches.resize(0);

    detail::AtomLabelFunctor atomLabeler(query,mol,useChirality);
    detail::BondLabelFunctor bondLabeler(query,mol);
    
    std::list<detail::ssPairType> pms;
    bool found=boost::ullmann_all(query.getTopology(),mol.getTopology(),
                                  atomLabeler,bondLabeler,pms);
    unsigned int res=0;
    if(found){
      unsigned int nQueryAtoms=query.getNumAtoms();
      matches.reserve(pms.size());
      for(std::list<detail::ssPairType>::const_iterator iter1=pms.begin();
          iter1!=pms.end();++iter1){
        MatchVectType matchVect;
        matchVect.resize(nQueryAtoms);
        for(detail::ssPairType::const_iterator iter2=iter1->begin();
            iter2!=iter1->end();++iter2){
          matchVect[iter2->first]=std::pair<int,int>(iter2->first,iter2->second);
        }
        matches.push_back(matchVect);
      }
      if(uniquify){
        removeDuplicates(matches);
      }
      res = matches.size();
    } 
    return res;
  }

  // ----------------------------------------------
  //
  // Intended for internal use 
  //
  unsigned int RecursiveMatcher(const ROMol &mol,const ROMol &query,
				std::vector< int > &matches,bool useChirality,
                                bool registerQuery)
  {
    ROMol::ConstAtomIterator atIt;
    for(atIt=query.beginAtoms();atIt!=query.endAtoms();atIt++){
      if((*atIt)->getQuery()){
        MatchSubqueries(mol,(*atIt)->getQuery(),useChirality,registerQuery);
      }
    }
 
    detail::AtomLabelFunctor atomLabeler(query,mol,useChirality);
    detail::BondLabelFunctor bondLabeler(query,mol);

    matches.clear();
    matches.resize(0);
    std::list<detail::ssPairType> pms;
    bool found=boost::ullmann_all(query.getTopology(),mol.getTopology(),
                            atomLabeler,bondLabeler,pms);

    unsigned int res=0;
    if(found){
      matches.reserve(pms.size());
      for(std::list<detail::ssPairType>::const_iterator iter1=pms.begin();
          iter1!=pms.end();++iter1){
        matches.push_back(iter1->begin()->second);
      }
      res = matches.size();
    } 
    //std::cout << " <<< RecursiveMatcher: " << int(query) << std::endl;
    return res;
  }

  void MatchSubqueries(const ROMol &mol,QueryAtom::QUERYATOM_QUERY *query,bool useChirality,
                       bool registerQuery){
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
        unsigned int res = RecursiveMatcher(mol,*queryMol,matchStarts,useChirality,registerQuery);
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
      MatchSubqueries(mol,childIt->get(),useChirality,registerQuery);
    }
    //std::cout << "<<- back " << (int)query << std::endl;
  }

};
#endif
