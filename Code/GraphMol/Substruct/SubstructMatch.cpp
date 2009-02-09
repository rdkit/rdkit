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
          res=atomCompat(d_query.getAtomWithIdx(i),d_mol.getAtomWithIdx(j));
        } else {
          res=chiralAtomCompat(d_query.getAtomWithIdx(i),d_mol.getAtomWithIdx(j));
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
        bool res=bondCompat(d_query.getBondPMap()[i],d_mol.getBondPMap()[j]);
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
    bool res=boost::ullmann(*query.getTopology(),*mol.getTopology(),
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
    bool found=boost::ullmann_all(*query.getTopology(),*mol.getTopology(),
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
    bool found=boost::ullmann_all(*query.getTopology(),*mol.getTopology(),
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
