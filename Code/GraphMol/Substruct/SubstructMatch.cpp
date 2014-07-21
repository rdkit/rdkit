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
#include <map>
#ifdef RDK_THREADSAFE_SSS
#include <boost/thread/mutex.hpp>
#endif

#include "ullmann.hpp"
#include "vf2.hpp"

namespace RDKit{
  namespace detail {
    typedef std::map<unsigned int,QueryAtom::QUERYATOM_QUERY *> SUBQUERY_MAP;
    
    void MatchSubqueries(const ROMol &mol,QueryAtom::QUERYATOM_QUERY *q,bool useChirality,
			 SUBQUERY_MAP &subqueryMap,bool useQueryQueryMatches);
#ifdef RDK_THREADSAFE_SSS
    void ClearSubqueryLocks(QueryAtom::QUERYATOM_QUERY *q);
#endif
    typedef std::list<std::pair<MolGraph::vertex_descriptor,MolGraph::vertex_descriptor> > ssPairType;

    class MolMatchFinalCheckFunctor {
    public:
      MolMatchFinalCheckFunctor(const ROMol &query,const ROMol &mol, bool useChirality) :
        d_query(query), d_mol(mol), df_useChirality(useChirality) {};
      bool operator()(const boost::detail::node_id c1[], const boost::detail::node_id c2[]) const {
        //std::cerr<<"  check! "<<df_useChirality<<std::endl;
        if(!df_useChirality) return true;
        //for(unsigned int i=0;i<d_query.getNumAtoms();++i){
        //  std::cerr<<"    "<<c1[i]<<" "<<c2[i]<<std::endl;
        //}

        // check chiral atoms:
        for(unsigned int i=0;i<d_query.getNumAtoms();++i){
          const Atom *qAt=d_query.getAtomWithIdx(c1[i]);
          if(qAt->getChiralTag()!=Atom::CHI_TETRAHEDRAL_CW &&
             qAt->getChiralTag()!=Atom::CHI_TETRAHEDRAL_CCW) continue;
          const Atom *mAt=d_mol.getAtomWithIdx(c2[i]);
          if(mAt->getChiralTag()!=Atom::CHI_TETRAHEDRAL_CW &&
             mAt->getChiralTag()!=Atom::CHI_TETRAHEDRAL_CCW) return false;
          if(qAt->getDegree()!=mAt->getDegree()) return false;
          INT_LIST qOrder;
          for(unsigned int j=0;j<d_query.getNumAtoms();++j){
            const Bond *qB = d_query.getBondBetweenAtoms(c1[i],c1[j]);
            if(qB){
              qOrder.push_back(qB->getIdx());
              if(qOrder.size()==qAt->getDegree()) break;
            }
          }
          int qPermCount=qAt->getPerturbationOrder(qOrder);

          INT_LIST mOrder;
          for(unsigned int j=0;j<d_query.getNumAtoms();++j){
            const Bond *mB = d_mol.getBondBetweenAtoms(c2[i],c2[j]);
            if(mB){
              mOrder.push_back(mB->getIdx());
              if(mOrder.size()==mAt->getDegree()) break;
            }
          }
          int mPermCount=mAt->getPerturbationOrder(mOrder);

          if((qPermCount%2 == mPermCount%2 &&
              qAt->getChiralTag()!=mAt->getChiralTag()) ||
             (qPermCount%2 != mPermCount%2 &&
              qAt->getChiralTag()==mAt->getChiralTag())) return false;
          
        }

        // now check double bonds
        for(unsigned int i=0;i<d_query.getNumBonds();++i){
          const Bond *qBnd=d_query.getBondWithIdx(i);
          if(qBnd->getBondType()!=Bond::DOUBLE ||
             (qBnd->getStereo()!=Bond::STEREOZ &&
              qBnd->getStereo()!=Bond::STEREOE)) continue;

          // don't think this can actually happen, but check to be sure:
          if(qBnd->getStereoAtoms().size()!=2) continue;

          std::map<unsigned int,unsigned int> qMap;
          for(unsigned int j=0;j<d_query.getNumAtoms();++j){
            qMap[c1[j]]=j;
          }
          const Bond *mBnd=d_mol.getBondBetweenAtoms(c2[qMap[qBnd->getBeginAtomIdx()]],
                                                     c2[qMap[qBnd->getEndAtomIdx()]]);
          CHECK_INVARIANT(mBnd,"Matching bond not found");
          if(mBnd->getBondType()!=Bond::DOUBLE ||
             (mBnd->getStereo()!=Bond::STEREOZ &&
              mBnd->getStereo()!=Bond::STEREOE)) continue;
          // don't think this can actually happen, but check to be sure:
          if(mBnd->getStereoAtoms().size()!=2) continue;

          unsigned int end1Matches=0;
          unsigned int end2Matches=0;
          if(c2[qMap[qBnd->getBeginAtomIdx()]]==mBnd->getBeginAtomIdx()){
            // query Begin == mol Begin
            if(c2[qMap[qBnd->getStereoAtoms()[0]]]==mBnd->getStereoAtoms()[0]) end1Matches=1;
            if(c2[qMap[qBnd->getStereoAtoms()[1]]]==mBnd->getStereoAtoms()[1]) end2Matches=1;
          } else {
            // query End == mol Begin
            if(c2[qMap[qBnd->getStereoAtoms()[0]]]==mBnd->getStereoAtoms()[1]) end1Matches=1;
            if(c2[qMap[qBnd->getStereoAtoms()[1]]]==mBnd->getStereoAtoms()[0]) end2Matches=1;
          }
          //std::cerr<<"  bnd: "<<qBnd->getIdx()<<":"<<qBnd->getStereo()<<" - "<<mBnd->getIdx()<<":"<<mBnd->getStereo()<<"  --  "<<end1Matches<<" "<<end2Matches<<std::endl;
          if(mBnd->getStereo()==qBnd->getStereo() && (end1Matches+end2Matches)==1) return false;
          if(mBnd->getStereo()!=qBnd->getStereo() && (end1Matches+end2Matches)!=1) return false;
        }
        
          

        return true;
      }
    private:
      const ROMol &d_query;
      const ROMol &d_mol;
      bool df_useChirality;
    };

    class AtomLabelFunctor{
    public:
      AtomLabelFunctor(const ROMol &query,const ROMol &mol, bool useChirality,
                       bool useQueryQueryMatches) :
        d_query(query), d_mol(mol), df_useChirality(useChirality),
        df_useQueryQueryMatches(useQueryQueryMatches) {};
      bool operator()(unsigned int i,unsigned int j) const{
        bool res=false;
        if(df_useChirality){
          const Atom *qAt=d_query.getAtomWithIdx(i);
          if(qAt->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW ||
             qAt->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW) {
            const Atom *mAt=d_mol.getAtomWithIdx(j);
            if(mAt->getChiralTag()!=Atom::CHI_TETRAHEDRAL_CW &&
               mAt->getChiralTag()!=Atom::CHI_TETRAHEDRAL_CCW) return false;
          }
        }
        res=atomCompat(d_query[i],d_mol[j],df_useQueryQueryMatches);
        return res;
      }
    private:
      const ROMol &d_query;
      const ROMol &d_mol;
      bool df_useChirality;
      bool df_useQueryQueryMatches;
    };
    class BondLabelFunctor{
    public:
      BondLabelFunctor(const ROMol &query,const ROMol &mol,bool useChirality,
                       bool useQueryQueryMatches) :
        d_query(query), d_mol(mol),df_useChirality(useChirality),
        df_useQueryQueryMatches(useQueryQueryMatches) {};
      bool operator()(MolGraph::edge_descriptor i,MolGraph::edge_descriptor j) const{
        if(df_useChirality){
          const BOND_SPTR qBnd=d_query[i];
          if(qBnd->getBondType()==Bond::DOUBLE &&
             (qBnd->getStereo()==Bond::STEREOZ ||
              qBnd->getStereo()==Bond::STEREOE)){
            const BOND_SPTR mBnd=d_mol[j];
            if(mBnd->getBondType()==Bond::DOUBLE &&
               !(mBnd->getStereo()==Bond::STEREOZ ||
                 mBnd->getStereo()==Bond::STEREOE))
              return false;
          }
        }
        bool res=bondCompat(d_query[i],d_mol[j]);
        return res;
      }
    private:
      const ROMol &d_query;
      const ROMol &d_mol;
      bool df_useChirality;
      bool df_useQueryQueryMatches;
    };
  }    
  
  // ----------------------------------------------
  //
  // find one match
  //
  bool SubstructMatch(const ROMol &mol,const ROMol &query,MatchVectType &matchVect,
                      bool recursionPossible,bool useChirality,bool useQueryQueryMatches)
  {

    //std::cerr<<"begin match"<<std::endl;
    if(recursionPossible){
      ROMol::ConstAtomIterator atIt;
      detail::SUBQUERY_MAP subqueryMap;
      for(atIt=query.beginAtoms();atIt!=query.endAtoms();atIt++){
        if((*atIt)->getQuery()){
	  detail::MatchSubqueries(mol,(*atIt)->getQuery(),useChirality,subqueryMap,
                                  useQueryQueryMatches);
        }
      }
    }
    //std::cerr<<"main matching"<<std::endl;

    matchVect.clear();
    matchVect.resize(0);

    detail::MolMatchFinalCheckFunctor matchChecker(query,mol,useChirality);
    detail::AtomLabelFunctor atomLabeler(query,mol,useChirality,useQueryQueryMatches);
    detail::BondLabelFunctor bondLabeler(query,mol,useChirality,useQueryQueryMatches);

    detail::ssPairType match;
#if 0
    bool res=boost::ullmann(query.getTopology(),mol.getTopology(),
                            atomLabeler,bondLabeler,match);
#else
    bool res=boost::vf2(query.getTopology(),mol.getTopology(),
                        atomLabeler,bondLabeler,matchChecker,match);
#endif
    if(res){
     matchVect.resize(query.getNumAtoms());
     for(detail::ssPairType::const_iterator iter=match.begin();iter!=match.end();++iter){
       matchVect[iter->first]=std::pair<int,int>(iter->first,iter->second);
     }
    }    

#ifdef RDK_THREADSAFE_SSS
    if(recursionPossible){
      ROMol::ConstAtomIterator atIt;
      for(atIt=query.beginAtoms();atIt!=query.endAtoms();atIt++){
        if((*atIt)->getQuery()){
	  detail::ClearSubqueryLocks((*atIt)->getQuery());
        }
      }
    }
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
			      bool useChirality,bool useQueryQueryMatches){

    if(recursionPossible){
      detail::SUBQUERY_MAP subqueryMap;
      ROMol::ConstAtomIterator atIt;
      for(atIt=query.beginAtoms();atIt!=query.endAtoms();atIt++){
        if((*atIt)->getQuery()){
          //std::cerr<<"recurse from atom "<<(*atIt)->getIdx()<<std::endl;
	  detail::MatchSubqueries(mol,(*atIt)->getQuery(),useChirality,subqueryMap,
                                  useQueryQueryMatches);
        }
      }
    }

    matches.clear();
    matches.resize(0);

    detail::AtomLabelFunctor atomLabeler(query,mol,useChirality,useQueryQueryMatches);
    detail::BondLabelFunctor bondLabeler(query,mol,useChirality,useQueryQueryMatches);
    detail::MolMatchFinalCheckFunctor matchChecker(query,mol,useChirality);
    
    std::list<detail::ssPairType> pms;
#if 0
    bool found=boost::ullmann_all(query.getTopology(),mol.getTopology(),
                                  atomLabeler,bondLabeler,pms);
#else
    bool found=boost::vf2_all(query.getTopology(),mol.getTopology(),
                              atomLabeler,bondLabeler,matchChecker,pms);
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

#ifdef RDK_THREADSAFE_SSS
    if(recursionPossible){
      ROMol::ConstAtomIterator atIt;
      for(atIt=query.beginAtoms();atIt!=query.endAtoms();atIt++){
        if((*atIt)->getQuery()){
	  detail::ClearSubqueryLocks((*atIt)->getQuery());
        }
      }
    }
#endif
    return res;
  }

  namespace detail {
    unsigned int RecursiveMatcher(const ROMol &mol,const ROMol &query,
				  std::vector< int > &matches,bool useChirality,
				  SUBQUERY_MAP &subqueryMap,bool useQueryQueryMatches)
    {
      ROMol::ConstAtomIterator atIt;
      for(atIt=query.beginAtoms();atIt!=query.endAtoms();atIt++){
	if((*atIt)->getQuery()){
	  MatchSubqueries(mol,(*atIt)->getQuery(),useChirality,subqueryMap,
                          useQueryQueryMatches);
	}
      }
 
      detail::AtomLabelFunctor atomLabeler(query,mol,useChirality,useQueryQueryMatches);
      detail::BondLabelFunctor bondLabeler(query,mol,useChirality,useQueryQueryMatches);
      detail::MolMatchFinalCheckFunctor matchChecker(query,mol,useChirality);

      matches.clear();
      matches.resize(0);
      std::list<detail::ssPairType> pms;
#if 0
      bool found=boost::ullmann_all(query.getTopology(),mol.getTopology(),
				    atomLabeler,bondLabeler,pms);
#else
      bool found=boost::vf2_all(query.getTopology(),mol.getTopology(),
				atomLabeler,bondLabeler,matchChecker,pms);
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
	      if(pairIter->first==static_cast<unsigned int>(rootIdx)){
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
			 SUBQUERY_MAP &subqueryMap,bool useQueryQueryMatches){
      PRECONDITION(query,"bad query");
      //std::cout << "*-*-* MS: " << (int)query << std::endl;
      //std::cout << "\t\t" << typeid(*query).name() << std::endl;
      if(query->getDescription()=="RecursiveStructure"){
	RecursiveStructureQuery *rsq=(RecursiveStructureQuery *)query;
        //std::cerr<<"lock: "<<rsq<<std::endl;
#ifdef RDK_THREADSAFE_SSS
        rsq->d_mutex.lock();
#endif
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
						subqueryMap,useQueryQueryMatches);
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
	MatchSubqueries(mol,childIt->get(),useChirality,subqueryMap,useQueryQueryMatches);
      }
      //std::cout << "<<- back " << (int)query << std::endl;
    }

#ifdef RDK_THREADSAFE_SSS
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
#endif
  } // end of namespace detail
}

