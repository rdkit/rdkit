// $Id$
//
//  Copyright (C) 2001-2010 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include "SubstructMatch.h"
#include "SubstructUtils.h"
#include <boost/smart_ptr.hpp>
#include <boost/thread/mutex.hpp>
#include <map>

#include "ullmann.hpp"
#include "vf2.hpp"

namespace RDKit{
  namespace detail {
    typedef std::map<unsigned int,QueryAtom::QUERYATOM_QUERY *> SUBQUERY_MAP;
    
    void MatchSubqueries(const ROMol &mol,QueryAtom::QUERYATOM_QUERY *q,bool useChirality,
			 SUBQUERY_MAP &subqueryMap );
    void ClearSubqueryLocks(QueryAtom::QUERYATOM_QUERY *q);

    typedef std::list<std::pair<MolGraph::vertex_descriptor,MolGraph::vertex_descriptor> > ssPairType;

    class AtomLabelFunctor{
    public:
      AtomLabelFunctor(const ROMol &query,const ROMol &mol, bool useChirality) :
        d_query(query), d_mol(mol), df_useChirality(useChirality) {};
      bool operator()(unsigned int i,unsigned int j) const{
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
      bool operator()(MolGraph::edge_descriptor i,MolGraph::edge_descriptor j) const{
        bool res=bondCompat(d_query[i],d_mol[j]);
        //std::cerr<<" blf: "<<i<<" - "<<j<<"? "<<res<<std::endl;
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
                      bool recursionPossible,bool useChirality)
  {

    //std::cerr<<"begin match"<<std::endl;
    if(recursionPossible){
      ROMol::ConstAtomIterator atIt;
      detail::SUBQUERY_MAP subqueryMap;
      for(atIt=query.beginAtoms();atIt!=query.endAtoms();atIt++){
        if((*atIt)->getQuery()){
	  detail::MatchSubqueries(mol,(*atIt)->getQuery(),useChirality,
				  subqueryMap);
        }
      }
    }
    //std::cerr<<"main matching"<<std::endl;

    matchVect.clear();
    matchVect.resize(0);

    detail::AtomLabelFunctor atomLabeler(query,mol,useChirality);
    detail::BondLabelFunctor bondLabeler(query,mol);

    detail::ssPairType match;
#if 0
    bool res=boost::ullmann(query.getTopology(),mol.getTopology(),
                            atomLabeler,bondLabeler,match);
#else
    bool res=boost::vf2(query.getTopology(),mol.getTopology(),
                        atomLabeler,bondLabeler,match);
#endif
    if(res){
     matchVect.resize(query.getNumAtoms());
     for(detail::ssPairType::const_iterator iter=match.begin();iter!=match.end();++iter){
       matchVect[iter->first]=std::pair<int,int>(iter->first,iter->second);
     }
    }    

    if(recursionPossible){
      ROMol::ConstAtomIterator atIt;
      for(atIt=query.beginAtoms();atIt!=query.endAtoms();atIt++){
        if((*atIt)->getQuery()){
	  detail::ClearSubqueryLocks((*atIt)->getQuery());
        }
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
			      bool useChirality){

    if(recursionPossible){
      detail::SUBQUERY_MAP subqueryMap;
      ROMol::ConstAtomIterator atIt;
      for(atIt=query.beginAtoms();atIt!=query.endAtoms();atIt++){
        if((*atIt)->getQuery()){
          //std::cerr<<"recurse from atom "<<(*atIt)->getIdx()<<std::endl;
	  detail::MatchSubqueries(mol,(*atIt)->getQuery(),useChirality,subqueryMap);
        }
      }
    }

    matches.clear();
    matches.resize(0);

    detail::AtomLabelFunctor atomLabeler(query,mol,useChirality);
    detail::BondLabelFunctor bondLabeler(query,mol);
    
    std::list<detail::ssPairType> pms;
#if 0
    bool found=boost::ullmann_all(query.getTopology(),mol.getTopology(),
                                  atomLabeler,bondLabeler,pms);
#else
    bool found=boost::vf2_all(query.getTopology(),mol.getTopology(),
                                  atomLabeler,bondLabeler,pms);
#endif
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
        removeDuplicates(matches,mol.getNumAtoms());
      }
      res = matches.size();
    } 

    if(recursionPossible){
      ROMol::ConstAtomIterator atIt;
      for(atIt=query.beginAtoms();atIt!=query.endAtoms();atIt++){
        if((*atIt)->getQuery()){
	  detail::ClearSubqueryLocks((*atIt)->getQuery());
        }
      }
    }

    return res;
  }

  namespace detail {
    unsigned int RecursiveMatcher(const ROMol &mol,const ROMol &query,
				  std::vector< int > &matches,bool useChirality,
				  SUBQUERY_MAP &subqueryMap)
    {
      ROMol::ConstAtomIterator atIt;
      for(atIt=query.beginAtoms();atIt!=query.endAtoms();atIt++){
	if((*atIt)->getQuery()){
	  MatchSubqueries(mol,(*atIt)->getQuery(),useChirality,subqueryMap);
	}
      }
 
      detail::AtomLabelFunctor atomLabeler(query,mol,useChirality);
      detail::BondLabelFunctor bondLabeler(query,mol);

      matches.clear();
      matches.resize(0);
      std::list<detail::ssPairType> pms;
#if 0
      bool found=boost::ullmann_all(query.getTopology(),mol.getTopology(),
				    atomLabeler,bondLabeler,pms);
#else
      bool found=boost::vf2_all(query.getTopology(),mol.getTopology(),
				atomLabeler,bondLabeler,pms);
#endif
      unsigned int res=0;
      if(found){
	matches.reserve(pms.size());
	for(std::list<detail::ssPairType>::const_iterator iter1=pms.begin();
	    iter1!=pms.end();++iter1){
	  if(!query.hasProp("_queryRootAtom")){
	    matches.push_back(iter1->begin()->second);
	  } else {
	    int rootIdx;
	    query.getProp("_queryRootAtom",rootIdx);
	    bool found=false;
	    for(detail::ssPairType::const_iterator pairIter=iter1->begin();
		pairIter!=iter1->end();++pairIter){
	      if(pairIter->first==rootIdx){
		matches.push_back(pairIter->second);
		found=true;
		break;
	      }
	    }
	    if(!found){
	      BOOST_LOG(rdErrorLog)<<"no match found for queryRootAtom"<<std::endl;
	    }
	  }
	}
	res = matches.size();
      } 
      //std::cout << " <<< RecursiveMatcher: " << int(query) << std::endl;
      return res;
    }

    void MatchSubqueries(const ROMol &mol,QueryAtom::QUERYATOM_QUERY *query,bool useChirality,
			 SUBQUERY_MAP &subqueryMap){
      PRECONDITION(query,"bad query");
      //std::cout << "*-*-* MS: " << (int)query << std::endl;
      //std::cout << "\t\t" << typeid(*query).name() << std::endl;
      if(query->getDescription()=="RecursiveStructure"){
	RecursiveStructureQuery *rsq=(RecursiveStructureQuery *)query;
        //std::cerr<<"lock: "<<rsq<<std::endl;
        rsq->d_mutex.lock();
	rsq->clear();
	bool matchDone=false;
	if(rsq->getSerialNumber() &&
	   subqueryMap.find(rsq->getSerialNumber()) != subqueryMap.end()){
	  // we've matched an equivalent serial number before, just
	  // copy in the matches:
	  matchDone=true;
	  const RecursiveStructureQuery *orsq=
	    (const RecursiveStructureQuery *)subqueryMap[rsq->getSerialNumber()];
	  for(RecursiveStructureQuery::CONTAINER_TYPE::const_iterator setIter=orsq->beginSet();
	      setIter!=orsq->endSet();++setIter){
	    rsq->insert(*setIter);
	  }
	  //std::cerr<<" copying results for query serial number: "<<rsq->getSerialNumber()<<std::endl;
	}
	
	if(!matchDone){
	  ROMol const *queryMol = rsq->getQueryMol();
	  // in case we are reusing this query, clear its contents now.
	  if(queryMol){
	    std::vector< int > matchStarts;
	    unsigned int res = RecursiveMatcher(mol,*queryMol,matchStarts,useChirality,
						subqueryMap);
	    if(res){
	      for(std::vector<int>::iterator i=matchStarts.begin();
		  i!=matchStarts.end();
		  i++){
		rsq->insert(*i);
	      }
	    }
	  }
	  if(rsq->getSerialNumber()){
	    subqueryMap[rsq->getSerialNumber()]=query;
	    //std::cerr<<" storing results for query serial number: "<<rsq->getSerialNumber()<<std::endl;
	  }

	}
      } else {
	//std::cout << "\tmsq1: "; 
      }
  
      // now recurse over our children (these things can be nested)
      Queries::Query<int,Atom const*,true>::CHILD_VECT_CI childIt;
      //std::cout << query << " " << query->endChildren()-query->beginChildren() <<  std::endl;
      for(childIt=query->beginChildren();childIt!=query->endChildren();childIt++){
	MatchSubqueries(mol,childIt->get(),useChirality,subqueryMap);
      }
      //std::cout << "<<- back " << (int)query << std::endl;
    }

    void ClearSubqueryLocks(QueryAtom::QUERYATOM_QUERY *query){
      PRECONDITION(query,"bad query");
      if(query->getDescription()=="RecursiveStructure"){
	RecursiveStructureQuery *rsq=(RecursiveStructureQuery *)query;
        //std::cerr<<"unlock: "<<rsq<<std::endl;
        rsq->d_mutex.unlock();

        for(ROMol::ConstAtomIterator atIt=rsq->getQueryMol()->beginAtoms();
            atIt!=rsq->getQueryMol()->endAtoms();atIt++){
          if((*atIt)->getQuery()){
            ClearSubqueryLocks((*atIt)->getQuery());
          }
        }
      }

      // now recurse over our children (these things can be nested)
      for(Queries::Query<int,Atom const*,true>::CHILD_VECT_CI childIt=query->beginChildren();
          childIt!=query->endChildren();childIt++){
	ClearSubqueryLocks(childIt->get());
      }
    }

  } // end of namespace detail
}

