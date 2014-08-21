// $Id$
//
//  Copyright (c) 2014, Novartis Institutes for BioMedical Research Inc.
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
//       products derived from this software without specific prior written permission.
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

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/QueryOps.h>
#include <boost/dynamic_bitset.hpp>
#include <boost/foreach.hpp>
#include <map>
#include <algorithm>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include "GraphMol/ChemReactions/ReactionRunner.h"
#include <RDGeneral/Invariant.h>

namespace RDKit {
  typedef std::vector<MatchVectType> VectMatchVectType;
  typedef std::vector< VectMatchVectType > VectVectMatchVectType;
    
  namespace ReactionRunnerUtils {
    //! returns whether or not all reactants matched
    bool getReactantMatches(const MOL_SPTR_VECT &reactants,
                            const ChemicalReaction &rxn,
                            VectVectMatchVectType &matchesByReactant)
    {
      PRECONDITION(reactants.size()==rxn.getNumReactantTemplates(),"reactant size mismatch");
      
      matchesByReactant.clear();
      matchesByReactant.resize(reactants.size());

      bool res=true;
      unsigned int i=0;
      for(MOL_SPTR_VECT::const_iterator iter = rxn.beginReactantTemplates();
            		  iter != rxn.endReactantTemplates(); ++iter, i++){
        std::vector< MatchVectType > matchesHere;
        // NOTE that we are *not* uniquifying the results.
        //   This is because we need multiple matches in reactions. For example, 
        //   The ring-closure coded as:
        //     [C:1]=[C:2] + [C:3]=[C:4][C:5]=[C:6] -> [C:1]1[C:2][C:3][C:4]=[C:5][C:6]1
        //   should give 4 products here:
        //     [Cl]C=C + [Br]C=CC=C ->
        //       [Cl]C1C([Br])C=CCC1
        //       [Cl]C1CC(Br)C=CC1
        //       C1C([Br])C=CCC1[Cl]
        //       C1CC([Br])C=CC1[Cl]
        //   Yes, in this case there are only 2 unique products, but that's
        //   a factor of the reactants' symmetry.
        //   
        //   There's no particularly straightforward way of solving this problem of recognizing cases
        //   where we should give all matches and cases where we shouldn't; it's safer to just
        //   produce everything and let the client deal with uniquifying their results.
        int matchCount=SubstructMatch(*(reactants[i]),*iter->get(), matchesHere, false, true, false);
        BOOST_FOREACH(const MatchVectType &match, matchesHere){
          bool keep=true;
          int pIdx,mIdx;
          BOOST_FOREACH(boost::tie(pIdx,mIdx),match){
            if(reactants[i]->getAtomWithIdx(mIdx)->hasProp("_protected")){
              keep=false;
              break;
            }
          }
          if(keep){
            matchesByReactant[i].push_back(match);
          } else {
            --matchCount;
          }
        }
        if(!matchCount){
          // no point continuing if we don't match one of the reactants:
          res=false;
          break;
        }
      }
      return res;  
    } // end of getReactantMatches()
    
    void recurseOverReactantCombinations(const VectVectMatchVectType &matchesByReactant,
                                         VectVectMatchVectType &matchesPerProduct,
                                         unsigned int level,VectMatchVectType combination){
      unsigned int nReactants=matchesByReactant.size();
      RANGE_CHECK(0,level,nReactants-1);
      PRECONDITION(combination.size()==nReactants,"bad combination size");
      for(VectMatchVectType::const_iterator reactIt=matchesByReactant[level].begin();
          reactIt != matchesByReactant[level].end(); ++reactIt){
        VectMatchVectType prod=combination;
        prod[level] = *reactIt;
        if(level==nReactants-1){
          // this is the bottom of the recursion:
          matchesPerProduct.push_back(prod);
        } else {
          recurseOverReactantCombinations(matchesByReactant,matchesPerProduct,level+1,prod);  
        }
      }
    } //end of recurseOverReactantCombinations

    void updateImplicitAtomProperties(Atom *prodAtom,const Atom *reactAtom){
      PRECONDITION(prodAtom,"no product atom");
      PRECONDITION(reactAtom,"no reactant atom");
      if(prodAtom->getAtomicNum()!=reactAtom->getAtomicNum()){
        // if we changed atom identity all bets are off, just
        // return
        return;
      }
      if(!prodAtom->hasProp("_QueryFormalCharge")){
        prodAtom->setFormalCharge(reactAtom->getFormalCharge());
      }
      if(!prodAtom->hasProp("_QueryIsotope")){
        prodAtom->setIsotope(reactAtom->getIsotope());
      }
      if(!prodAtom->hasProp("_ReactionDegreeChanged")){
        if(!prodAtom->hasProp("_QueryHCount")){
          prodAtom->setNumExplicitHs(reactAtom->getNumExplicitHs());
        }
        prodAtom->setNoImplicit(reactAtom->getNoImplicit());
      }
    }    

    void generateReactantCombinations(const VectVectMatchVectType &matchesByReactant,
                                      VectVectMatchVectType &matchesPerProduct){
      matchesPerProduct.clear();
      VectMatchVectType tmp;
      tmp.clear();
      tmp.resize(matchesByReactant.size());
      recurseOverReactantCombinations(matchesByReactant,matchesPerProduct,0,tmp);
    } // end of generateReactantCombinations()

    RWMOL_SPTR initProduct(const ROMOL_SPTR prodTemplateSptr){
      const ROMol *prodTemplate=prodTemplateSptr.get();
      RWMol *res=new RWMol();

      // --------- --------- --------- --------- --------- --------- 
      // Initialize by making a copy of the product template as a normal molecule.
      // NOTE that we can't just use a normal copy because we do not want to end up
      // with query atoms or bonds in the product.
      
      // copy in the atoms:
      ROMol::ATOM_ITER_PAIR atItP = prodTemplate->getVertices();
      while(atItP.first != atItP.second ){
        Atom *oAtom=(*prodTemplate)[*(atItP.first++)].get();
        Atom *newAtom=new Atom(*oAtom);
        res->addAtom(newAtom,false,true);
        if(newAtom->hasProp("molAtomMapNumber")){
          // set bookmarks for the mapped atoms:
          int mapNum;
          newAtom->getProp("molAtomMapNumber",mapNum);
          res->setAtomBookmark(newAtom,mapNum);
          // now clear the molAtomMapNumber property so that it doesn't
          // end up in the products (this was bug 3140490):
          newAtom->clearProp("molAtomMapNumber");
        }

        newAtom->setChiralTag(Atom::CHI_UNSPECIFIED);
        // if the product-template atom has the inversion flag set
        // to 4 (=SET), then bring its stereochem over, otherwise we'll
        // ignore it:
        if(oAtom->hasProp("molInversionFlag")){
          int iFlag;
          oAtom->getProp("molInversionFlag",iFlag);
          if(iFlag==4) newAtom->setChiralTag(oAtom->getChiralTag());
        }

        // check for properties we need to set:
        if(newAtom->hasProp("_QueryFormalCharge")){
          int val;
          newAtom->getProp("_QueryFormalCharge",val);
          newAtom->setFormalCharge(val);
        }
        if(newAtom->hasProp("_QueryHCount")){
          int val;
          newAtom->getProp("_QueryHCount",val);
          newAtom->setNumExplicitHs(val);
        }
        if(newAtom->hasProp("_QueryMass")){
          int val;
          newAtom->getProp("_QueryMass",val);
          newAtom->setMass(val);
        }
        if(newAtom->hasProp("_QueryIsotope")){
          int val;
          newAtom->getProp("_QueryIsotope",val);
          newAtom->setIsotope(val);
        }
      }
      // and the bonds:
      ROMol::BOND_ITER_PAIR bondItP = prodTemplate->getEdges();
      while(bondItP.first != bondItP.second ){
        const BOND_SPTR oldB=(*prodTemplate)[*(bondItP.first++)];
        unsigned int bondIdx;
        bondIdx=res->addBond(oldB->getBeginAtomIdx(),oldB->getEndAtomIdx(),oldB->getBondType())-1;
        // make sure we don't lose the bond dir information:
        Bond *newB=res->getBondWithIdx(bondIdx);
        newB->setBondDir(oldB->getBondDir());
        // Special case/hack:
        //  The product has been processed by the SMARTS parser.
        //  The SMARTS parser tags unspecified bonds as single, but then adds
        //  a query so that they match single or double
        //  This caused Issue 1748846
        //   http://sourceforge.net/tracker/index.php?func=detail&aid=1748846&group_id=160139&atid=814650
        //  We need to fix that little problem now:
        if( oldB->hasQuery()){
          //  remember that the product has been processed by the SMARTS parser.
          std::string queryDescription=oldB->getQuery()->getDescription();
          if(queryDescription=="BondOr" &&
             oldB->getBondType()==Bond::SINGLE){
            //  We need to fix that little problem now:
            if(newB->getBeginAtom()->getIsAromatic() && newB->getEndAtom()->getIsAromatic()){
              newB->setBondType(Bond::AROMATIC);
              newB->setIsAromatic(true);
            } else {
              newB->setBondType(Bond::SINGLE);
              newB->setIsAromatic(false);
            }
          } else if(queryDescription=="BondNull") {
            newB->setProp("NullBond",1);
          }
        }
      }
      return RWMOL_SPTR(res);  
    } // end of initProduct()

    struct ReactantProductAtomMapping{

    	ReactantProductAtomMapping(unsigned lenghtBitSet)
    	{
    		mappedAtoms.resize(lenghtBitSet);
    		skippedAtoms.resize(lenghtBitSet);
    	}

		boost::dynamic_bitset<> mappedAtoms;
		boost::dynamic_bitset<> skippedAtoms;
		std::map<unsigned int,std::vector<unsigned int> > reactProdAtomMap;
		std::map<unsigned int,unsigned int> prodReactAtomMap;
		std::map<unsigned int,unsigned int> prodAtomBondMap;
    };

    ReactantProductAtomMapping* getAtomMappingsReactantProduct(const MatchVectType &match,
    		const ROMol& reactantTemplate, RWMOL_SPTR product, unsigned numReactAtoms)
    {
      ReactantProductAtomMapping *mapping = new ReactantProductAtomMapping(numReactAtoms);
      for(unsigned int i=0;i<match.size();i++){
        const Atom *templateAtom=reactantTemplate.getAtomWithIdx(match[i].first);
        if(templateAtom->hasProp("molAtomMapNumber")){
          int molAtomMapNumber;
          templateAtom->getProp("molAtomMapNumber",molAtomMapNumber);

          if(product->hasAtomBookmark(molAtomMapNumber)){
        	RWMol::ATOM_PTR_LIST atomIdxs = product->getAllAtomsWithBookmark(molAtomMapNumber);
        	for(RWMol::ATOM_PTR_LIST::iterator iter = atomIdxs.begin();
        		iter != atomIdxs.end(); ++iter){
        	  Atom * a = *iter;
              unsigned int pIdx = a->getIdx();
              mapping->reactProdAtomMap[match[i].second].push_back(pIdx);
              mapping->mappedAtoms[match[i].second]=1;
              CHECK_INVARIANT(pIdx<product->getNumAtoms(),"yikes!");
              mapping->prodReactAtomMap[pIdx]=match[i].second;
        	}
          }
          else {
            // this skippedAtom has an atomMapNumber, but it's not in this product
            // (it's either in another product or it's not mapped at all).
          	mapping->skippedAtoms[match[i].second]=1;
          }
        }
        else {
          // This skippedAtom appears in the match, but not in a product:
      	  mapping->skippedAtoms[match[i].second]=1;
        }
      }
      return mapping;
    }

    void setReactantBondPropertiesToProduct(RWMOL_SPTR product, const ROMol& reactant,
    		ReactantProductAtomMapping* mapping)
    {
      ROMol::BOND_ITER_PAIR bondItP = product->getEdges();
      while(bondItP.first != bondItP.second ){
        BOND_SPTR pBond=(*product)[*(bondItP.first)];
        ++bondItP.first;
        if(pBond->hasProp("NullBond")){
          if(mapping->prodReactAtomMap.find(pBond->getBeginAtomIdx())!=mapping->prodReactAtomMap.end() &&
        	  mapping->prodReactAtomMap.find(pBond->getEndAtomIdx())!=mapping->prodReactAtomMap.end() ){
            // the bond is between two mapped atoms from this reactant:
        	unsigned begIdx = mapping->prodReactAtomMap[pBond->getBeginAtomIdx()];
        	unsigned endIdx = mapping->prodReactAtomMap[pBond->getEndAtomIdx()];
            const Bond *rBond=reactant.getBondBetweenAtoms(begIdx,endIdx);
            if(!rBond) continue;
            pBond->setBondType(rBond->getBondType());
            pBond->setBondDir(rBond->getBondDir());
            pBond->setIsAromatic(rBond->getIsAromatic());
            pBond->clearProp("NullBond");
          }
        }
      }
    }

    void checkProductChirality(Atom::ChiralType reactantChirality, Atom *productAtom)
    {
      int flagVal;
      productAtom->getProp("molInversionFlag",flagVal);

      switch(flagVal){
      case 0:
        // FIX: should we clear the chirality or leave it alone? for now we leave it alone
        //productAtom->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
        productAtom->setChiralTag(reactantChirality);
        break;
      case 1:
        // inversion
        if(reactantChirality!=Atom::CHI_TETRAHEDRAL_CW && reactantChirality!=Atom::CHI_TETRAHEDRAL_CCW){
          BOOST_LOG(rdWarningLog) << "unsupported chiral type on reactant atom ignored\n";
        }
        else {
          productAtom->setChiralTag(reactantChirality);
          productAtom->invertChirality();
        }
        break;
      case 2:
        // retention: just set to the reactant
        productAtom->setChiralTag(reactantChirality);
        break;
      case 3:
        // remove stereo
        productAtom->setChiralTag(Atom::CHI_UNSPECIFIED);
        break;
      case 4:
        // set stereo, so leave it the way it was in the product template
        break;
      default:
        BOOST_LOG(rdWarningLog) << "unrecognized chiral inversion/retention flag on product atom ignored\n";
      }
    }

    void setReactantAtomPropertiesToProduct(
    		Atom *productAtom, const Atom& reactantAtom,
    		bool setImplicitProperties)
    {
        // which properties need to be set from the reactant?
      if(productAtom->getAtomicNum()<=0){
        productAtom->setAtomicNum(reactantAtom.getAtomicNum());
        productAtom->setIsAromatic(reactantAtom.getIsAromatic());
        // don't copy isotope information over from dummy atoms
        productAtom->setIsotope(reactantAtom.getIsotope());

        // remove dummy labels (if present)
        if(productAtom->hasProp("dummyLabel")) productAtom->clearProp("dummyLabel");
        if(productAtom->hasProp("_MolFileRLabel")) productAtom->clearProp("_MolFileRLabel");
      }
      if(setImplicitProperties){
        updateImplicitAtomProperties(productAtom,&reactantAtom);
      }
      // One might be tempted to copy over the reactant atom's chirality into the
      // product atom if chirality is not specified on the product. This would be a
      // very bad idea because the order of bonds will almost certainly change on the
      // atom and the chirality is referenced to bond order.

      // --------- --------- --------- --------- --------- ---------
      // While we're here, set the stereochemistry
      // FIX: this should be free-standing, not in this function.
      if(reactantAtom.getChiralTag()!=Atom::CHI_UNSPECIFIED &&
         reactantAtom.getChiralTag()!=Atom::CHI_OTHER &&
         productAtom->hasProp("molInversionFlag")){
        checkProductChirality(reactantAtom.getChiralTag(), productAtom);
      }
    }

    void setNewProductBond(const Bond& origB, RWMOL_SPTR product,
    		unsigned bondBeginIdx, unsigned bondEndIdx)
    {
      unsigned bondIdx = product->addBond(bondBeginIdx, bondEndIdx, origB.getBondType())-1;
      Bond *newB = product->getBondWithIdx(bondIdx);
      newB->setBondDir(origB.getBondDir());
    }

    void addMissingProductBonds(const Bond& origB, RWMOL_SPTR product,
    		ReactantProductAtomMapping* mapping)
    {
      unsigned int begIdx=origB.getBeginAtomIdx();
      unsigned int endIdx=origB.getEndAtomIdx();

      std::vector<unsigned> prodBeginIdxs = mapping->reactProdAtomMap[begIdx];
      std::vector<unsigned> prodEndIdxs = mapping->reactProdAtomMap[endIdx];
      CHECK_INVARIANT(prodBeginIdxs.size() == prodEndIdxs.size(),
    	             "Different number of start-end points for product bonds.");
      for(unsigned i = 0; i < prodBeginIdxs.size(); i++){
        setNewProductBond(origB, product, prodBeginIdxs.at(i), prodEndIdxs.at(i));
      }
    }

    void addMissingProductAtom(const Atom& reactAtom, unsigned reactNeighborIdx,
    		unsigned prodNeighborIdx, RWMOL_SPTR product, const ROMol& reactant,
        	ReactantProductAtomMapping* mapping)
    {
      Atom *newAtom = new Atom(reactAtom);
      unsigned reactAtomIdx = reactAtom.getIdx();
      unsigned productIdx = product->addAtom(newAtom,false,true);
      mapping->reactProdAtomMap[reactAtomIdx].push_back(productIdx);
      mapping->prodReactAtomMap[productIdx] = reactAtomIdx;
      // add the bonds
      const Bond *origB = reactant.getBondBetweenAtoms(reactNeighborIdx,reactAtomIdx);
      unsigned int begIdx=origB->getBeginAtomIdx();
      unsigned int endIdx=origB->getEndAtomIdx();
      if(begIdx == reactNeighborIdx){
        setNewProductBond(*origB, product, prodNeighborIdx, productIdx);
      }
      else{
        setNewProductBond(*origB, product, productIdx, prodNeighborIdx);
      }
    }


    void addReactantNeighborsToProduct(const ROMol& reactant, const Atom& reactantAtom,
    		RWMOL_SPTR product, boost::dynamic_bitset<>& visitedAtoms,
    		std::vector<const Atom*>& chiralAtomsToCheck,
    		ReactantProductAtomMapping* mapping)
    {
      std::list<const Atom*> atomStack;
      atomStack.push_back(&reactantAtom);
      while(!atomStack.empty()){
        const Atom *lReactantAtom = atomStack.front();
        atomStack.pop_front();

        // each atom in the stack is guaranteed to already be in the product:
        CHECK_INVARIANT(mapping->reactProdAtomMap.find(lReactantAtom->getIdx()) !=
        		        mapping->reactProdAtomMap.end(),
                        "reactant atom on traversal stack not present in product.");

        std::vector<unsigned> lReactantAtomProductIndex =
        		mapping->reactProdAtomMap[lReactantAtom->getIdx()];
        unsigned lreactIdx = lReactantAtom->getIdx();
        visitedAtoms[lreactIdx] = 1;
        // Check our neighbors:
        ROMol::ADJ_ITER nbrIdx,endNbrs;
        boost::tie(nbrIdx,endNbrs) = reactant.getAtomNeighbors(lReactantAtom);
        while(nbrIdx!=endNbrs){
          // Four possibilities here. The neighbor:
          //  0) has been visited already: do nothing
          //  1) is part of the match (thus already in the product): set a bond to it
          //  2) has been added: set a bond to it
          //  3) has not yet been added: add it, set a bond to it, and push it
          //     onto the stack
          if(!visitedAtoms[*nbrIdx] && !mapping->skippedAtoms[*nbrIdx]){
            if(mapping->mappedAtoms[*nbrIdx]){
              // this is case 1 (neighbor in match); set a bond to the neighbor if this atom
              // is not also in the match (match-match bonds were set when the product template was
              // copied in to start things off).;
              if(!mapping->mappedAtoms[lreactIdx]){
                CHECK_INVARIANT(mapping->reactProdAtomMap.find(*nbrIdx)!=mapping->reactProdAtomMap.end(),
                              "reactant atom not present in product.");
                const Bond *origB=reactant.getBondBetweenAtoms(lreactIdx,*nbrIdx);
                addMissingProductBonds(*origB, product,mapping);
              }
            }
            else if(mapping->reactProdAtomMap.find(*nbrIdx)!= mapping->reactProdAtomMap.end()){
              // case 2, the neighbor has been added and we just need to set a bond to it:
              const Bond *origB=reactant.getBondBetweenAtoms(lreactIdx,*nbrIdx);
              addMissingProductBonds(*origB, product,mapping);
            }
            else {
              // case 3, add the atom, a bond to it, and push the atom onto the stack
              const Atom *neighbor=reactant.getAtomWithIdx(*nbrIdx);
              for(unsigned i = 0; i < lReactantAtomProductIndex.size(); i++){
            	addMissingProductAtom(*neighbor, lreactIdx,lReactantAtomProductIndex[i],
            	      		product, reactant, mapping);
              }
              // update the stack:
              atomStack.push_back(neighbor);
              // if the atom is chiral, we need to check its bond ordering later:
              if(neighbor->getChiralTag()!=Atom::CHI_UNSPECIFIED){
                chiralAtomsToCheck.push_back(neighbor);
              }
            }
          }
          nbrIdx++;
        }
      } // end of atomStack traversal
    }

    void checkAndCorrectChiralityOfMatchingAtomsInProduct(
    		const ROMol& reactant, unsigned reactantAtomIdx,
    		const Atom& reactantAtom, RWMOL_SPTR product,
    		ReactantProductAtomMapping* mapping)
    {
      for(unsigned i = 0; i < mapping->reactProdAtomMap[reactantAtomIdx].size(); i++){
        unsigned productAtomIdx = mapping->reactProdAtomMap[reactantAtomIdx][i];
        Atom *productAtom = product->getAtomWithIdx(productAtomIdx);

        if(productAtom->getChiralTag() != Atom::CHI_UNSPECIFIED ||
           reactantAtom.getChiralTag() == Atom::CHI_UNSPECIFIED ||
           reactantAtom.getChiralTag() == Atom::CHI_OTHER ||
           productAtom->hasProp("molInversionFlag")){
        	continue;
        }
        // we can only do something sensible here if we have the same number of bonds
        // in the reactants and the products:
        if(reactantAtom.getDegree() != productAtom->getDegree()){
          continue;
        }
        unsigned int nUnknown=0;
        INT_LIST pOrder;
        ROMol::ADJ_ITER nbrIdx,endNbrs;
        boost::tie(nbrIdx,endNbrs) = product->getAtomNeighbors(productAtom);
        while(nbrIdx != endNbrs){
          if(mapping->prodReactAtomMap.find(*nbrIdx) == mapping->prodReactAtomMap.end()){
            ++nUnknown;
            // if there's more than one bond in the product that doesn't correspond to
            // anything in the reactant, we're also doomed
            if(nUnknown > 1) break;
            // otherwise, add a -1 to the bond order that we'll fill in later
            pOrder.push_back(-1);
          }
          else {
            const Bond *rBond=reactant.getBondBetweenAtoms(reactantAtom.getIdx(),mapping->prodReactAtomMap[*nbrIdx]);
            CHECK_INVARIANT(rBond,"expected reactant bond not found");
            pOrder.push_back(rBond->getIdx());
          }
          ++nbrIdx;
        }
        if(nUnknown == 1){
          // find the reactant bond that hasn't yet been accounted for:
          int unmatchedBond = -1;
          boost::tie(nbrIdx,endNbrs) = reactant.getAtomNeighbors(&reactantAtom);
          while(nbrIdx!=endNbrs){
            const Bond *rBond=reactant.getBondBetweenAtoms(reactantAtom.getIdx(),*nbrIdx);
            if(std::find(pOrder.begin(),pOrder.end(),rBond->getIdx()) == pOrder.end()){
              unmatchedBond = rBond->getIdx();
              break;
            }
            ++nbrIdx;
          }
          // what must be true at this point:
          //  1) there's a -1 in pOrder that we'll substitute for
          //  2) unmatchedBond contains the index of the substitution
          INT_LIST::iterator bPos=std::find(pOrder.begin(),pOrder.end(),-1);
          if(unmatchedBond >= 0 && bPos != pOrder.end()){
            *bPos = unmatchedBond;
          }
          if(std::find(pOrder.begin(),pOrder.end(),-1) == pOrder.end()){
            nUnknown = 0;
          }
        }
        if(!nUnknown){
          productAtom->setChiralTag(reactantAtom.getChiralTag());
          int nSwaps = reactantAtom.getPerturbationOrder(pOrder);
          if(nSwaps%2){
            productAtom->invertChirality();
          }
        }
      }
    }

    void checkAndCorrectChiralityOfProduct(const std::vector<const Atom*>& chiralAtomsToCheck,
    		RWMOL_SPTR product, ReactantProductAtomMapping* mapping)
    {
      for(std::vector<const Atom *>::const_iterator atomIt = chiralAtomsToCheck.begin();
          atomIt != chiralAtomsToCheck.end(); ++atomIt){
        const Atom *reactantAtom = *atomIt;
        CHECK_INVARIANT(reactantAtom->getChiralTag()!=Atom::CHI_UNSPECIFIED,
        		"missing atom chirality.");
        unsigned int reactAtomDegree = reactantAtom->getOwningMol().getAtomDegree(reactantAtom);
        for(unsigned i = 0; i < mapping->reactProdAtomMap[reactantAtom->getIdx()].size(); i++){
          unsigned productAtomIdx = mapping->reactProdAtomMap[reactantAtom->getIdx()][i];
          Atom *productAtom = product->getAtomWithIdx(productAtomIdx);
          CHECK_INVARIANT(reactantAtom->getChiralTag() == productAtom->getChiralTag(),
                        "invalid product chirality.");

          if(reactAtomDegree != product->getAtomDegree(productAtom) ){
          // If the number of bonds to the atom has changed in the course of the
          // reaction we're lost, so remove chirality.
          //  A word of explanation here: the atoms in the chiralAtomsToCheck set are
          //  not explicitly mapped atoms of the reaction, so we really have no idea what
          //  to do with this case. At the moment I'm not even really sure how this
          //  could happen, but better safe than sorry.
            productAtom->setChiralTag(Atom::CHI_UNSPECIFIED);
          }
          else if(reactantAtom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
                    reactantAtom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW){
            // this will contain the indices of product bonds in the
            // reactant order:
            INT_LIST newOrder;
            ROMol::OEDGE_ITER beg,end;
            boost::tie(beg,end) = reactantAtom->getOwningMol().getAtomBonds(reactantAtom);
            while(beg!=end){
              const BOND_SPTR reactantBond = reactantAtom->getOwningMol()[*beg];
              unsigned int oAtomIdx = reactantBond->getOtherAtomIdx(reactantAtom->getIdx());
              CHECK_INVARIANT(mapping->reactProdAtomMap.find(oAtomIdx) != mapping->reactProdAtomMap.end(),
                            "other atom from bond not mapped.");
              const Bond *productBond;
              unsigned neighborBondIdx = mapping->reactProdAtomMap[oAtomIdx][i];
              productBond = product->getBondBetweenAtoms(productAtom->getIdx(),neighborBondIdx);
              CHECK_INVARIANT(productBond,"no matching bond found in product");
              newOrder.push_back(productBond->getIdx());
              ++beg;
            }
            int nSwaps=productAtom->getPerturbationOrder(newOrder);
            if(nSwaps%2){
              productAtom->invertChirality();
            }
          }
          else {
            // not tetrahedral chirality, don't do anything.
          }
        }
      } // end of loop over chiralAtomsToCheck
    }

    void generateProductConformers(Conformer *productConf, const ROMol& reactant,
    		ReactantProductAtomMapping* mapping)
    {
      if(!reactant.getNumConformers()){
    	  return;
      }
      const Conformer &reactConf=reactant.getConformer();
      if(reactConf.is3D()){
      	productConf->set3D(true);
      }
      for(std::map<unsigned int, std::vector<unsigned int> >::const_iterator pr = mapping->reactProdAtomMap.begin();
        pr!=mapping->reactProdAtomMap.end();++pr){
        std::vector<unsigned> prodIdxs =  pr->second;
    	if(prodIdxs.size()>1){
    	  BOOST_LOG(rdWarningLog) << "reactant atom match more than one product atom, coordinates need to be revised\n";
    	}
    	// is this reliable when multiple product atom mapping occurs????
    	for(unsigned i = 0; i < prodIdxs.size(); i++){
          productConf->setAtomPos(prodIdxs[i],reactConf.getAtomPos(pr->first));
    	}
      }
    }


    void addReactantAtomsAndBonds(const ChemicalReaction& rxn,
                                  RWMOL_SPTR product,const ROMOL_SPTR reactantSptr,
                                  const MatchVectType &match,
                                  const ROMOL_SPTR reactantTemplate,
                                  Conformer *productConf)
    {
      // start by looping over all matches and marking the reactant atoms that
      // have already been "added" by virtue of being in the product. We'll also
      // mark "skipped" atoms: those that are in the match, but not in this
      // particular product (or, perhaps, not in any product)
      // At the same time we'll set up a map between the indices of those
      // atoms and their index in the product.
      ReactantProductAtomMapping* mapping =
    		  getAtomMappingsReactantProduct(match,*reactantTemplate,product,reactantSptr->getNumAtoms());

      boost::dynamic_bitset<> visitedAtoms(reactantSptr->getNumAtoms());

      const ROMol *reactant=reactantSptr.get();

      // ---------- ---------- ---------- ---------- ---------- ----------
      // Loop over the bonds in the product and look for those that have
      // the NullBond property set. These are bonds for which no information
      // (other than their existance) was provided in the template:
      setReactantBondPropertiesToProduct(product, *reactant, mapping);

      // ---------- ---------- ---------- ---------- ---------- ----------
      // Loop over the atoms in the match that were added to the product
      // From the corresponding atom in the reactant, do a graph traversal
      // to find other connected atoms that should be added:

      std::vector<const Atom *> chiralAtomsToCheck;
      for(unsigned matchIdx=0; matchIdx<match.size(); matchIdx++){

    	int reactantAtomIdx = match[matchIdx].second;
    	if(mapping->mappedAtoms[reactantAtomIdx]){
          CHECK_INVARIANT(mapping->reactProdAtomMap.find(reactantAtomIdx)!=mapping->reactProdAtomMap.end(),
                          "mapped reactant atom not present in product.");

          const Atom *reactantAtom=reactant->getAtomWithIdx(reactantAtomIdx);
          for(unsigned i = 0; i < mapping->reactProdAtomMap[reactantAtomIdx].size(); i++){
            // here's a pointer to the atom in the product:
        	unsigned productAtomIdx = mapping->reactProdAtomMap[reactantAtomIdx][i];
            Atom *productAtom = product->getAtomWithIdx(productAtomIdx);
            setReactantAtomPropertiesToProduct(productAtom, *reactantAtom,rxn.getImplicitPropertiesFlag());
          }
          // now traverse:
          addReactantNeighborsToProduct(*reactant, *reactantAtom, product, visitedAtoms,
              		chiralAtomsToCheck, mapping);

          // now that we've added all the reactant's neighbors, check to see if
          // it is chiral in the reactant but is not in the reaction. If so
          // we need to worry about its chirality
          checkAndCorrectChiralityOfMatchingAtomsInProduct(*reactant, reactantAtomIdx, *reactantAtom, product, mapping);
    	}
      } // end of loop over matched atoms

      // ---------- ---------- ---------- ---------- ---------- ----------
      // now we need to loop over atoms from the reactants that were chiral but not
      // directly involved in the reaction in order to make sure their chirality hasn't
      // been disturbed
      checkAndCorrectChiralityOfProduct(chiralAtomsToCheck, product, mapping);

      // ---------- ---------- ---------- ---------- ---------- ----------
      // finally we may need to set the coordinates in the product conformer:
      if(productConf){
        productConf->resize(product->getNumAtoms());
        generateProductConformers(productConf, *reactant, mapping);
      }
    } // end of addReactantAtomsAndBonds


    MOL_SPTR_VECT generateOneProductSet(const ChemicalReaction& rxn,
    		const MOL_SPTR_VECT& reactants,
            const std::vector<MatchVectType> &reactantsMatch)
    {
      PRECONDITION(reactants.size() == reactantsMatch.size(),"vector size mismatch");

      // if any of the reactants have a conformer, we'll go ahead and
      // generate conformers for the products:
      bool doConfs=false;
      BOOST_FOREACH(ROMOL_SPTR reactant,reactants){
        if(reactant->getNumConformers()){
          doConfs=true;
          break;
        }
      }

      MOL_SPTR_VECT res;
      res.resize(rxn.getNumProductTemplates());
      unsigned int prodId=0;
      for(MOL_SPTR_VECT::const_iterator pTemplIt = rxn.beginProductTemplates();
          pTemplIt != rxn.endProductTemplates(); ++pTemplIt){
    	// copy product template and its properties to a new product RWMol
        RWMOL_SPTR product = initProduct(*pTemplIt);
        Conformer *conf = 0;
        if(doConfs){
          conf = new Conformer();
          conf->set3D(false);
        }
        
        unsigned int reactantId = 0;
        for(MOL_SPTR_VECT::const_iterator iter = rxn.beginReactantTemplates();
      		  iter != rxn.endReactantTemplates(); ++iter, reactantId++){
          addReactantAtomsAndBonds(rxn, product, reactants.at(reactantId),
                                                  reactantsMatch.at(reactantId),
                                                  *iter, conf);
        }
        product->clearAllAtomBookmarks();
        if(doConfs){
          product->addConformer(conf,true);
        }
        res[prodId] = product;
        ++prodId;
      }
      return res;
    }
  }// end namespace RectionRunnerUtils

  std::vector<MOL_SPTR_VECT> run_Reactants(const ChemicalReaction& rxn, const MOL_SPTR_VECT& reactants) {
    if(!rxn.isInitialized()) {
     throw ChemicalReactionException("initMatchers() must be called before runReactants()");
    }
    if(reactants.size() != rxn.getNumReactantTemplates()){
      throw ChemicalReactionException("Number of reactants provided does not match number of reactant templates.");
    }
    BOOST_FOREACH(ROMOL_SPTR msptr,reactants){
      CHECK_INVARIANT(msptr,"bad molecule in reactants");
    }

    std::vector<MOL_SPTR_VECT> productMols;
    productMols.clear();

    // if we have no products, return now:
    if(!rxn.getNumProductTemplates()){
      return productMols;
    }

    // find the matches for each reactant:
    VectVectMatchVectType matchesByReactant;
    if(!ReactionRunnerUtils::getReactantMatches(reactants, rxn, matchesByReactant)){
      // some reactants didn't find a match, return an empty product list:
      return productMols;
    }
    // -------------------------------------------------------
    // we now have matches for each reactant, so we can start creating products:
    // start by doing the combinatorics on the matches:
    VectVectMatchVectType reactantMatchesPerProduct;
    ReactionRunnerUtils::generateReactantCombinations(matchesByReactant,reactantMatchesPerProduct);
    productMols.resize(reactantMatchesPerProduct.size());

    for(unsigned int productId=0;productId!=productMols.size();++productId){
      MOL_SPTR_VECT lProds = ReactionRunnerUtils::generateOneProductSet(rxn, reactants,reactantMatchesPerProduct[productId]);
      productMols[productId]=lProds;
    }

    return productMols;
  } // end of ChemicalReaction::runReactants()

    
} // end of RDKit namespace
