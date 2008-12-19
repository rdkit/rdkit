// $Id$
//
//  Copyright (c) 2008, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met: 
//
//     * Redistributions of source code must retain the above copyright 
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following 
//       disclaimer in the documentation and/or other materials provided 
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
//       nor the names of its contributors may be used to endorse or promote 
//       products derived from this software without specific prior
//       written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  Created by Greg Landrum, Sept. 2006
//

#include <GraphMol/SLNParse/SLNParse.h>
#include <GraphMol/SLNParse/SLNAttribs.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/Invariant.h>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

namespace RDKit{
  namespace SLNParse {
    namespace {
      int parseIntAttribVal(std::string attribName,std::string attribVal,int (*defaultFunc)(Atom const * at)=NULL,
                            Atom *atom=NULL){
        PRECONDITION((!defaultFunc)||atom,"If a default func is provided, an atom must be as well."); 
        int iVal;
        boost::to_lower(attribVal);    

        if(defaultFunc && attribVal=="f"){
          iVal = defaultFunc(atom);
        } else {
          try{
            iVal=boost::lexical_cast<int>(attribVal);
          } catch (boost::bad_lexical_cast &){
            std::stringstream err;
            err << "SLN Parser error: bad integer value (" << attribVal << ") provided for property: " << attribName;
            throw SLNParseException(err.str());
          }
        }
        return iVal;
      } 
    } // end of anonymous namespace

    QueryAtom::QUERYATOM_QUERY *makeQueryFromOp(std::string op,int val,int (*func)(Atom const * at),
                                                std::string description){
      PRECONDITION(func,"bad query function");
      QueryAtom::QUERYATOM_QUERY *res=0;  
      if(op=="="){
        ATOM_EQUALS_QUERY *tmp=new ATOM_EQUALS_QUERY;  
        tmp->setVal(val);
        tmp->setDataFunc(func);
        tmp->setDescription(description);
        res=tmp;
      } else if(op=="!="){
        ATOM_EQUALS_QUERY *tmp=new ATOM_EQUALS_QUERY;  
        tmp->setVal(val);
        tmp->setDataFunc(func);
        tmp->setDescription(description);
        tmp->setNegation(true);
        res=tmp;
      } else if(op==">"){
        // don't be alarmed by this use of the LessEqual query for >, it's not a bug.
        // The RD GreaterQuery(tgt) returns true if tgt is greater than the thing you
        // compare to. In this case we need to reverse that because we're interested in
        // seeing if the value is greater than the target; this is equiv to asking if
        // the target is < the value.  
        ATOM_LESS_QUERY *tmp=new ATOM_LESS_QUERY;  
        tmp->setVal(val);
        tmp->setDataFunc(func);
        tmp->setDescription(description);
        res=tmp;
      } else if(op==">="){
        ATOM_LESSEQUAL_QUERY *tmp=new ATOM_LESSEQUAL_QUERY;  
        tmp->setVal(val);
        tmp->setDataFunc(func);
        tmp->setDescription(description);
        res=tmp;
      } else if(op=="<"){
        ATOM_GREATER_QUERY *tmp=new ATOM_GREATER_QUERY;  
        tmp->setVal(val);
        tmp->setDataFunc(func);
        tmp->setDescription(description);
        res=tmp;
      } else if(op=="<="){
        ATOM_GREATEREQUAL_QUERY *tmp=new ATOM_GREATEREQUAL_QUERY;  
        tmp->setVal(val);
        tmp->setDataFunc(func);
        tmp->setDescription(description);
        res=tmp;
      } else {
        std::stringstream err;
        err << "SLN Parser error: bad attribute operator (" << op << ") provided.";
        throw SLNParseException(err.str());
      }
      
      POSTCONDITION(res,"bad query");
      return res;
    }

    void parseAtomAttribs(Atom *atom,AttribListType attribs,bool doingQuery){
      for(AttribListType::const_iterator it=attribs.begin();
          it!=attribs.end();++it){
        Queries::CompositeQueryType how;
        switch(it->first) {
        case AttribAnd:
          how=Queries::COMPOSITE_AND;  break;
        case AttribOr:
          how=Queries::COMPOSITE_OR;  break;
        case AttribLowPriAnd:
          how=Queries::COMPOSITE_AND;  break;
        default:
          throw SLNParseException("unrecognized query composition operator");
        }

        boost::shared_ptr<AttribType> attribPtr=it->second;
        std::string attribName=attribPtr->first;
        boost::to_lower(attribName);    
        std::string attribVal=attribPtr->second;
        if(attribName=="charge"){
          int chg=0;
          if(attribVal=="-"){
            chg=-1;
          } else if(attribVal=="+"){
            chg = +1;
          } else {
            chg = parseIntAttribVal(attribName,attribVal);
          }
          if(!doingQuery){
            atom->setFormalCharge(chg);
          } else {
            QueryAtom::QUERYATOM_QUERY *query=makeQueryFromOp(attribPtr->op,chg,queryAtomFormalCharge,
                                                              "AtomFormalCharge");
            if(attribPtr->negated) query->setNegation(!query->getNegation());
            atom->expandQuery(query,how);
          }
        } else if(attribName=="i"){
          int val=parseIntAttribVal(attribName,attribVal);
          if(!doingQuery){
            atom->setMass(static_cast<double>(val));
          } else {
            QueryAtom::QUERYATOM_QUERY *query=makeQueryFromOp(attribPtr->op,val,queryAtomMass,
                                                              "AtomMass");
            if(attribPtr->negated) query->setNegation(!query->getNegation());
            atom->expandQuery(query,how);
          }
        } else if(attribName=="r"){
          if(attribVal!=""){
            BOOST_LOG(rdWarningLog) << "Query value " << attribVal <<" ignored for r query\n";
          }
          if(!doingQuery){
            BOOST_LOG(rdWarningLog) << "Query property " << attribName <<" ignored on non-query atom\n";
          } else {
            QueryAtom::QUERYATOM_QUERY *query=makeAtomInRingQuery();
            if(attribPtr->negated) query->setNegation(true);
            atom->expandQuery(query,how);
          }
        } else if(attribName=="is"){
          // recursive queries:
          if(!attribPtr->structQuery) {
            throw SLNParseException("failed recursive query");
          } 
          atom->expandQuery(static_cast<QueryAtom::QUERYATOM_QUERY *>(attribPtr->structQuery),how);
        } else if(attribName=="s"){
          if(attribPtr->op!="="){
            std::stringstream err;
            err << "SLN Parser error: comparison operator '" << attribPtr->op << "' not supported for chiral attributes.\n";
            throw SLNParseException(err.str());
          }
          boost::to_lower(attribVal);
    
          if(attribVal[0]=='i' || attribVal[0]=='n'){
            if(attribVal.size()>1 && attribVal[1]=='*'){
              BOOST_LOG(rdWarningLog) << "Chiral modifier * ignored, chiral spec " << attribVal[0] <<" will be used\n";
            }
            if(attribVal.size()>1 && attribVal[1]=='m'){
              BOOST_LOG(rdWarningLog) << "Chiral modifier m ignored, chiral spec " << attribVal[0] <<" will be used\n";
            }
#if 0
          } else if(attribVal=="r"){
          } else if(attribVal=="r*"){
          } else if(attribVal=="rm"){
          } else if(attribVal=="s"){
          } else if(attribVal=="s*"){
          } else if(attribVal=="sm"){
#endif
          } else {
            BOOST_LOG(rdWarningLog) << "Unsupported stereochemistry specifier " << attribVal <<" ignored.\n";
          }              
          atom->setProp("_SLN_"+attribName,attribVal);
          atom->setProp("_OP_SLN_"+attribName,"=");
          atom->setProp("_HOW_SLN_"+attribName,static_cast<unsigned int>(how));
          atom->setProp("_NEG_SLN_"+attribName,(int)attribPtr->negated);
        } else {
          // anything we don't know how to deal now with we'll just store for later processing:
          atom->setProp("_SLN_"+attribName,attribVal);
          atom->setProp("_OP_SLN_"+attribName,attribPtr->op);
          atom->setProp("_HOW_SLN_"+attribName,static_cast<unsigned int>(how));
          atom->setProp("_NEG_SLN_"+attribName,(int)attribPtr->negated);
        }
      }
    }

    void parseFinalAtomAttribs(Atom *atom,bool doingQuery){
      STR_VECT propList=atom->getPropList();

      while(propList.size()){
        std::string attribName =propList.back();
        propList.pop_back(); 
        if(attribName.substr(0,5)=="_SLN_") {
          std::string attribVal;
          atom->getProp(attribName,attribVal);
          atom->clearProp(attribName);

          boost::erase_head(attribName,5);
          std::string attribOp;
          atom->getProp("_OP_SLN_"+attribName,attribOp);
          propList.erase(std::find(propList.begin(),propList.end(),"_OP_SLN_"+attribName));
          atom->clearProp("_OP_SLN_"+attribName);
          unsigned int howTmp;
          atom->getProp("_HOW_SLN_"+attribName,howTmp);
          propList.erase(std::find(propList.begin(),propList.end(),"_HOW_SLN_"+attribName));
          atom->clearProp("_HOW_SLN_"+attribName);

          int attribNeg;
          atom->getProp("_NEG_SLN_"+attribName,attribNeg);
          atom->clearProp("_NEG_SLN_"+attribName);
                  
          Queries::CompositeQueryType how=static_cast<Queries::CompositeQueryType>(howTmp);  
          boost::to_lower(attribName);    
          QueryAtom::QUERYATOM_QUERY *query=0;          
          if(attribName=="rbc"){
            int val=parseIntAttribVal(attribName,attribVal,queryAtomRingBondCount,atom);
            if(!doingQuery){
              BOOST_LOG(rdWarningLog) << "Query property " << attribName <<" ignored on non-query atom\n";
            } else {
              query=makeQueryFromOp(attribOp,val,queryAtomRingBondCount,
                                                                "AtomRingBondCount");
            }
          } else if(attribName=="tbo"){
            int val=parseIntAttribVal(attribName,attribVal,queryAtomTotalValence,atom);
            if(!doingQuery){
              BOOST_LOG(rdWarningLog) << "Query property " << attribName <<" ignored on non-query atom\n";
            } else {
              query=makeQueryFromOp(attribOp,val,queryAtomTotalValence,
                                                                "AtomTotalValence");
            }
          } else if(attribName=="tac"){
            int val=parseIntAttribVal(attribName,attribVal,queryAtomTotalDegree,atom);
            if(!doingQuery){
              BOOST_LOG(rdWarningLog) << "Query property " << attribName <<" ignored on non-query atom\n";
            } else {
              query=makeQueryFromOp(attribOp,val,queryAtomTotalDegree,
                                                                "AtomTotalDegree");
            }
          } else if(attribName=="hc"){
            int val=parseIntAttribVal(attribName,attribVal,queryAtomHCount,atom);
            if(!doingQuery){
              BOOST_LOG(rdWarningLog) << "Query property " << attribName <<" ignored on non-query atom\n";
            } else {
              query=makeQueryFromOp(attribOp,val,queryAtomHCount,
                                                                "AtomHCount");
            }
          } else if(attribName=="hac"){
            int val=parseIntAttribVal(attribName,attribVal,queryAtomHeavyAtomDegree,atom);
            if(!doingQuery){
              BOOST_LOG(rdWarningLog) << "Query property " << attribName <<" ignored on non-query atom\n";
            } else {
              query=makeQueryFromOp(attribOp,val,queryAtomHeavyAtomDegree,
                                                                "AtomHeavyAtomDegree");
            }
          } else if(attribName=="f"){
            int val=queryAtomAllBondProduct(atom);
            if(!doingQuery){
              BOOST_LOG(rdWarningLog) << "Query property " << attribName <<" ignored on non-query atom\n";
            } else {
              attribOp="=";
              query=makeQueryFromOp(attribOp,val,(int (*)(const RDKit::Atom*))(queryAtomAllBondProduct),
                                    "AtomBondEnvironment");
            }
          } else {
            // anything we don't know how to deal with we'll just store in raw form:
            atom->setProp(attribName,attribVal);
          }
          if(query){
            if(attribNeg) query->setNegation(!query->getNegation());
            atom->expandQuery(query,how);
          }
        }
      }
    }

    void parseBondAttribs(Bond *bond,AttribListType attribs,bool doingQuery){
      bool seenTypeQuery=false;
      for(AttribListType::const_iterator it=attribs.begin();
          it!=attribs.end();++it){
        Queries::CompositeQueryType how;
        switch(it->first) {
        case AttribAnd:
          how=Queries::COMPOSITE_AND;  break;
        case AttribOr:
          how=Queries::COMPOSITE_OR;  break;
        case AttribLowPriAnd:
          how=Queries::COMPOSITE_AND;  break;
        default:
          throw SLNParseException("unrecognized query composition operator");
        }

        boost::shared_ptr<AttribType> attribPtr=it->second;
        std::string attribName=attribPtr->first;
        boost::to_lower(attribName);    
        std::string attribVal=attribPtr->second;
        if(attribName=="type"){
          boost::to_lower(attribVal);
          Bond::BondType bondType;
          if(attribVal=="-" || attribVal=="1"){
            bondType = Bond::SINGLE;
          } else if(attribVal=="=" || attribVal=="2"){
            bondType = Bond::DOUBLE;
          } else if(attribVal=="#" || attribVal=="3"){
            bondType = Bond::TRIPLE;
          } else if(attribVal==":" || attribVal=="aromatic"){
            bondType = Bond::AROMATIC;
          } else {
            bondType = Bond::OTHER;
            bond->setProp("SLN_Type",attribVal);
          }
          if(!doingQuery){
            bond->setBondType(bondType);
          } else {
            QueryBond::QUERYBOND_QUERY *query=makeBondOrderEqualsQuery(bondType);
            if(attribPtr->negated) query->setNegation(!query->getNegation());
            if(seenTypeQuery){
              static_cast<RDKit::QueryBond *>(bond)->expandQuery(query,how,true);
            } else {
              // if this is the first type query, we need to replace any existing bond order queries:
              // FIX: this replaces tooo much, ring queries also get blown out
              bond->setQuery(query);
            }
            seenTypeQuery=true;
          }
        } else if(attribName=="r"){
          if(attribVal!=""){
            BOOST_LOG(rdWarningLog) << "Query value " << attribVal <<" ignored for r query\n";
          }
          if(!doingQuery){
            BOOST_LOG(rdWarningLog) << "Query property " << attribName <<" ignored on non-query bond\n";
          } else {
            QueryBond::QUERYBOND_QUERY *query=makeBondIsInRingQuery();
            if(attribPtr->negated) query->setNegation(true);
            static_cast<QueryBond *>(bond)->expandQuery(query,how);
          }
        } else {
          // anything we don't know how to deal with we'll just store in raw form:
          bond->setProp(attribName,attribVal);
        }
      } 
    }

    void parseMolAttribs(ROMol *mol,AttribListType attribs){
      for(AttribListType::const_iterator it=attribs.begin();
          it!=attribs.end();++it){
        CHECK_INVARIANT(it->first==AttribAnd,"bad attrib type");

        boost::shared_ptr<AttribType> attribPtr=it->second;
        std::string attribName=attribPtr->first;
        boost::to_lower(attribName);    
        std::string attribVal=attribPtr->second;
        if(attribVal.begin()!=attribVal.end() &&
           *(attribVal.begin())=='"' &&
           *(attribVal.begin())==*(attribVal.rbegin()) ){
          attribVal.erase(attribVal.begin());
          attribVal.erase(--(attribVal.end()));
        }
        if(attribName=="name"){
          mol->setProp("_Name",attribVal);
        } else {
          mol->setProp(attribName,attribVal);
        }
      } 
    }


    void adjustAtomChiralities(RWMol *mol){
      for(RWMol::AtomIterator atomIt=mol->beginAtoms();
          atomIt != mol->endAtoms();
          atomIt++){
        if((*atomIt)->hasProp("_SLN_s")){
          // the atom is marked as chiral, translate the sln chirality into
          // RDKit chirality
          std::string attribVal;
          (*atomIt)->getProp("_SLN_s",attribVal);
          
          // start with a straight map of the chirality value:
          // as a reminder, here are some SLN <-> SMILES pairs
          //      C[s=n]H(Cl)(F)Br  <->  [C@@H](Cl)(F)Br          (CHI_TETRAHEDRAL_CW)
          //      ClC[s=n]H(F)Br  <->  Cl[C@H](F)Br               (CHI_TETRAHEDRAL_CCW)
          //      FC[1:s=n](Cl)OCH2@1   <->  F[C@@]1(Cl)OC1       (CHI_TETRAHEDRAL_CW)
          if(attribVal[0]=='n'){
            (*atomIt)->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);  
          } else if(attribVal[0]=='i'){
            (*atomIt)->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);  
          }
          unsigned int numPrecedingAtoms=0;
          std::list< std::pair<int,int> > neighbors;
          RWMol::ADJ_ITER nbrIdx,endNbrs;
          boost::tie(nbrIdx,endNbrs) = mol->getAtomNeighbors(*atomIt);
          while(nbrIdx != endNbrs){
            Bond *nbrBond=mol->getBondBetweenAtoms((*atomIt)->getIdx(),*nbrIdx);
            neighbors.push_back(std::make_pair(*nbrIdx,nbrBond->getIdx()));
            if(*nbrIdx<(*atomIt)->getIdx()) ++numPrecedingAtoms;
            ++nbrIdx;
          }          

          //std::cerr << "CHIRAL " << (*atomIt)->getIdx();

          // sort by neighbor idx:
          neighbors.sort();
          // figure out the bond ordering:
          std::list<int> bondOrdering;
          for(std::list< std::pair<int,int> >::const_iterator nbrIt=neighbors.begin();
              nbrIt!= neighbors.end();++nbrIt){
            bondOrdering.push_back(nbrIt->second);
            //std::cerr << " " << nbrIt->second;
          }

          // ok, we now have the ordering of the bonds (used for RDKit chirality), 
          // figure out the permutation order relative to the atom numbering
          // (sln chirality):
          int nSwaps=(*atomIt)->getPerturbationOrder(bondOrdering);
          //nSwaps += numPrecedingAtoms;
          
          // if the atom has an explicit hydrogen, that means it was int
          //std::cerr << " swaps: " << nSwaps << std::endl;
          
          
          if(nSwaps%2){
            (*atomIt)->setChiralTag((*atomIt)->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW ?
                                    Atom::CHI_TETRAHEDRAL_CCW :
                                    Atom::CHI_TETRAHEDRAL_CW);
          }
        }
#if 0    
        Atom::ChiralType chiralType=(*atomIt)->getChiralTag();
        if(chiralType==Atom::CHI_TETRAHEDRAL_CW ||
           chiralType==Atom::CHI_TETRAHEDRAL_CCW ){
          //
          // The atom is marked as chiral, set the SMILES-order of the
          // atom's bonds.  This is easy for non-ring-closure bonds,
          // because the SMILES order is determined solely by the atom
          // indices.  Things are trickier for ring-closure bonds, which we
          // need to insert into the list in a particular order
          //
          INT_VECT ringClosures;
          if((*atomIt)->hasProp("_RingClosures"))
            (*atomIt)->getProp("_RingClosures",ringClosures);
          std::list<INT_PAIR> neighbors;
          // push this atom onto the list of neighbors (we'll use this
          // to find our place later):
          neighbors.push_back(std::make_pair((*atomIt)->getIdx(),-1));
          std::list<int> bondOrder;
          RWMol::ADJ_ITER nbrIdx,endNbrs;
          boost::tie(nbrIdx,endNbrs) = mol->getAtomNeighbors(*atomIt);
          while(nbrIdx != endNbrs){
            Bond *nbrBond=mol->getBondBetweenAtoms((*atomIt)->getIdx(),*nbrIdx);
            if(std::find(ringClosures.begin(),
                         ringClosures.end(),
                         nbrBond->getIdx())== ringClosures.end()){
              neighbors.push_back(std::make_pair(*nbrIdx,nbrBond->getIdx()));
            }
            nbrIdx++;
          }
          // sort the list of non-ring-closure bonds:
          neighbors.sort();
  
  
          // find the location of this atom.  it pretty much has to be
          // first in the list, e.g for smiles like [C@](F)(Cl)(Br)I, or
          // second (everything else).
          std::list<INT_PAIR>::iterator selfPos=neighbors.begin();
          if(selfPos->first != (*atomIt)->getIdx()){
            selfPos++;
          }
          CHECK_INVARIANT(selfPos->first==(*atomIt)->getIdx(),"weird atom ordering");
          
          // copy over the bond ids:
          INT_LIST bondOrdering;
          for(std::list<INT_PAIR>::iterator neighborIt=neighbors.begin();
              neighborIt != neighbors.end(); neighborIt++){
            if(neighborIt != selfPos){
              bondOrdering.push_back(neighborIt->second);
            } else {
              // we are not going to add the atom itself, but we will push on
              // ring closure bonds at this point (if required):
              for(INT_VECT::const_iterator closureIt=ringClosures.begin();
                  closureIt != ringClosures.end();closureIt++){
                bondOrdering.push_back(*closureIt);
              }
            }
          }
  
          // ok, we now have the SMILES ordering of the bonds, figure out the
          // permutation order:
          int nSwaps=(*atomIt)->getPerturbationOrder(bondOrdering);
          if(nSwaps%2){
            (*atomIt)->setChiralTag(chiralType==Atom::CHI_TETRAHEDRAL_CW ?
                                    Atom::CHI_TETRAHEDRAL_CCW :
                                    Atom::CHI_TETRAHEDRAL_CW);
          }
        }
#endif        
      }
    }
  } // end of SLNParse namespace
} // end of RDKit namespace
