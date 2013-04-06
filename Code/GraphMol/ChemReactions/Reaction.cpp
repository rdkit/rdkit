// $Id$
//
//  Copyright (c) 2007, Novartis Institutes for BioMedical Research Inc.
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
#include <GraphMol/ChemReactions/ReactionPickler.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/QueryOps.h>
#include <boost/dynamic_bitset.hpp>
#include <boost/foreach.hpp>
#include <map>
#include <algorithm>
#include <GraphMol/ChemTransforms/ChemTransforms.h>

namespace RDKit {
  typedef std::vector<MatchVectType> VectMatchVectType;
  typedef std::vector< VectMatchVectType > VectVectMatchVectType;
    
  namespace ReactionUtils {
    //! returns whether or not all reactants matched
    bool getReactantMatches(const MOL_SPTR_VECT &reactants,
                            const MOL_SPTR_VECT &reactantTemplates,
                            VectVectMatchVectType &matchesByReactant){
      PRECONDITION(reactants.size()==reactantTemplates.size(),"reactant size mismatch");
      
      matchesByReactant.clear();
      matchesByReactant.resize(reactants.size());

      bool res=true;
      for(unsigned int i=0;i<reactants.size();++i){
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
        int matchCount=SubstructMatch(*(reactants[i]),*(reactantTemplates[i]),matchesHere,false,true,false);
        BOOST_FOREACH(const MatchVectType &match,matchesHere){
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
    
    void addReactantAtomsAndBonds(const ChemicalReaction *rxn,
                                  RWMOL_SPTR product,const ROMOL_SPTR reactantSptr,
                                  const MatchVectType &match,
                                  const ROMOL_SPTR reactantTemplate,
                                  Conformer *productConf){
      PRECONDITION(rxn,"bad reaction");
      // start by looping over all matches and marking the reactant atoms that
      // have already been "added" by virtue of being in the product. We'll also
      // mark "skipped" atoms: those that are in the match, but not in this
      // particular product (or, perhaps, not in any product)
      // At the same time we'll set up a map between the indices of those
      // atoms and their index in the product.
      boost::dynamic_bitset<> mappedAtoms(reactantSptr->getNumAtoms());
      boost::dynamic_bitset<> skippedAtoms(reactantSptr->getNumAtoms());
      std::map<unsigned int,unsigned int> reactProdAtomMap; // this maps atom indices from reactant->product
      std::map<unsigned int,unsigned int> prodReactAtomMap; // this maps atom indices from product->reactant
      std::vector<const Atom *> chiralAtomsToCheck; 
      for(unsigned int i=0;i<match.size();i++){
        const Atom *templateAtom=reactantTemplate->getAtomWithIdx(match[i].first);
        if(templateAtom->hasProp("molAtomMapNumber")){
          int molAtomMapNumber;
          templateAtom->getProp("molAtomMapNumber",molAtomMapNumber);
          
          if(product->hasAtomBookmark(molAtomMapNumber)){
            unsigned int pIdx=product->getAtomWithBookmark(molAtomMapNumber)->getIdx();
            reactProdAtomMap[match[i].second] = pIdx;
            mappedAtoms[match[i].second]=1;
            CHECK_INVARIANT(pIdx<product->getNumAtoms(),"yikes!");
            prodReactAtomMap[pIdx]=match[i].second;
          } else {
            // this skippedAtom has an atomMapNumber, but it's not in this product 
            // (it's either in another product or it's not mapped at all).
            skippedAtoms[match[i].second]=1;  
          }
        } else {
          // This skippedAtom appears in the match, but not in a product:
          skippedAtoms[match[i].second]=1;  
        }
      }

      boost::dynamic_bitset<> visitedAtoms(reactantSptr->getNumAtoms());

      const ROMol *reactant=reactantSptr.get();

      // ---------- ---------- ---------- ---------- ---------- ---------- 
      // Loop over the bonds in the product and look for those that have
      // the NullBond property set. These are bonds for which no information
      // (other than their existance) was provided in the template:
      ROMol::BOND_ITER_PAIR bondItP = product->getEdges();
      while(bondItP.first != bondItP.second ){
        BOND_SPTR pBond=(*product)[*(bondItP.first)];
        ++bondItP.first;
        if(pBond->hasProp("NullBond")){
          if(prodReactAtomMap.find(pBond->getBeginAtomIdx())!=prodReactAtomMap.end() &&
             prodReactAtomMap.find(pBond->getEndAtomIdx())!=prodReactAtomMap.end() ){
            // the bond is between two mapped atoms from this reactant:
            const Bond *rBond=reactant->getBondBetweenAtoms(prodReactAtomMap[pBond->getBeginAtomIdx()],
                                                            prodReactAtomMap[pBond->getEndAtomIdx()]);
            if(!rBond) continue;
            pBond->setBondType(rBond->getBondType());
            pBond->setBondDir(rBond->getBondDir());
            pBond->setIsAromatic(rBond->getIsAromatic());
            pBond->clearProp("NullBond");
          }
        }
      }

      // ---------- ---------- ---------- ---------- ---------- ---------- 
      // Loop over the atoms in the match that were added to the product
      // From the corresponding atom in the reactant, do a graph traversal
      //  to find other connected atoms that should be added: 
      for(unsigned int matchIdx=0;matchIdx<match.size();matchIdx++){
        int reactantAtomIdx=match[matchIdx].second;
        if(mappedAtoms[reactantAtomIdx]){
          CHECK_INVARIANT(reactProdAtomMap.find(reactantAtomIdx)!=reactProdAtomMap.end(),
                          "mapped reactant atom not present in product.");
          
          // here's a pointer to the atom in the product:
          Atom *productAtom = product->getAtomWithIdx(reactProdAtomMap[reactantAtomIdx]);
          // and this is the corresponding atom in the reactant.
          const Atom *reactantAtom=reactant->getAtomWithIdx(reactantAtomIdx);

          if(rxn->getImplicitPropertiesFlag()){
            // --------- --------- --------- --------- --------- ---------
            // which properties need to be set from the reactant?
            if(productAtom->getAtomicNum()<=0){
              // If the product atom is a dummy, set everything
              productAtom->setAtomicNum(reactantAtom->getAtomicNum());
              productAtom->setIsAromatic(reactantAtom->getIsAromatic());
              // now that the atomic number is set, we need
              // to reset the isotope so that the mass is also correct:
              if(productAtom->getIsotope())
                productAtom->setIsotope(productAtom->getIsotope());              
            }
            updateImplicitAtomProperties(productAtom,reactantAtom);
          }
          // One might be tempted to copy over the reactant atom's chirality into the
          // product atom if chirality is not specified on the product. This would be a
          // very bad idea because the order of bonds will almost certainly change on the
          // atom and the chirality is referenced to bond order.

          // --------- --------- --------- --------- --------- --------- 
          // While we're here, set the stereochemistry 
          // FIX: this should be free-standing, not in this function.
          if(reactantAtom->getChiralTag()!=Atom::CHI_UNSPECIFIED &&
             reactantAtom->getChiralTag()!=Atom::CHI_OTHER &&
             productAtom->hasProp("molInversionFlag")){
            int flagVal;
            productAtom->getProp("molInversionFlag",flagVal);
            switch(flagVal){
            case 0:
              // FIX: should we clear the chirality or leave it alone? for now we leave it alone 
              //productAtom->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
              productAtom->setChiralTag(reactantAtom->getChiralTag());
              break;
            case 1:
              // inversion
              if(reactantAtom->getChiralTag()!=Atom::CHI_TETRAHEDRAL_CW &&
                 reactantAtom->getChiralTag()!=Atom::CHI_TETRAHEDRAL_CCW){
                BOOST_LOG(rdWarningLog) << "unsupported chiral type on reactant atom ignored\n";
              } else {
                productAtom->setChiralTag(reactantAtom->getChiralTag());
                productAtom->invertChirality();
              }
              break;
            case 2:
              // retention: just set to the reactant
              productAtom->setChiralTag(reactantAtom->getChiralTag());
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

          // now traverse:
          std::list< const Atom * > atomStack;
          atomStack.push_back(reactantAtom);
          while(!atomStack.empty()){
            const Atom *lReactantAtom = atomStack.front();
            atomStack.pop_front();

            // each atom in the stack is guaranteed to already be in the product:
            CHECK_INVARIANT(reactProdAtomMap.find(lReactantAtom->getIdx())!=reactProdAtomMap.end(),
                            "reactant atom on traversal stack not present in product.");
            int lReactantAtomProductIndex=reactProdAtomMap[lReactantAtom->getIdx()];
            productAtom = product->getAtomWithIdx(lReactantAtomProductIndex); 
            visitedAtoms[lReactantAtom->getIdx()]=1;

            // Check our neighbors:
            ROMol::ADJ_ITER nbrIdx,endNbrs;
            boost::tie(nbrIdx,endNbrs) = reactant->getAtomNeighbors(lReactantAtom);
            while(nbrIdx!=endNbrs){
              // Four possibilities here. The neighbor:
              //  0) has been visited already: do nothing
              //  1) is part of the match (thus already in the product): set a bond to it
              //  2) has been added: set a bond to it
              //  3) has not yet been added: add it, set a bond to it, and push it
              //     onto the stack
              if(!visitedAtoms[*nbrIdx] && !skippedAtoms[*nbrIdx]){
                unsigned int productIdx;
                bool addBond=false;
                if(mappedAtoms[*nbrIdx]){
                  // this is case 1 (neighbor in match); set a bond to the neighbor if this atom
                  // is not also in the match (match-match bonds were set when the product template was
                  // copied in to start things off).;
                  if(!mappedAtoms[lReactantAtom->getIdx()]){
                    CHECK_INVARIANT(reactProdAtomMap.find(*nbrIdx)!=reactProdAtomMap.end(),
                                  "reactant atom not present in product.");
                    addBond=true;
                  }                
                } else if(reactProdAtomMap.find(*nbrIdx)!=reactProdAtomMap.end()){
                  // case 2, the neighbor has been added and we just need to set a bond to it:
                  addBond=true;
                } else {
                  // case 3, add the atom, a bond to it, and push the atom onto the stack
                  const Atom *lReactantAtom=reactant->getAtomWithIdx(*nbrIdx);
                  Atom *newAtom = new Atom(*lReactantAtom);
                  productIdx=product->addAtom(newAtom,false,true);
                  reactProdAtomMap[*nbrIdx]=productIdx;
                  prodReactAtomMap[productIdx]=*nbrIdx;
                  addBond=true;
                  // update the stack:
                  atomStack.push_back(lReactantAtom);
                  // if the atom is chiral, we need to check its bond ordering later:
                  if(lReactantAtom->getChiralTag()!=Atom::CHI_UNSPECIFIED){
                    chiralAtomsToCheck.push_back(lReactantAtom);
                  }
                }
                if(addBond){
                  const Bond *origB=reactant->getBondBetweenAtoms(lReactantAtom->getIdx(),*nbrIdx);
                  unsigned int begIdx=origB->getBeginAtomIdx();
                  unsigned int endIdx=origB->getEndAtomIdx();
                  unsigned int bondIdx;
                  // add the bond, but make sure it has the same begin and end
                  // atom indices as the original:
                  bondIdx=product->addBond(reactProdAtomMap[begIdx],
                                           reactProdAtomMap[endIdx],
                                           origB->getBondType())-1;
                  //bondIdx=product->addBond(reactProdAtomMap[*nbrIdx],lReactantAtomProductIndex,
                  //                         origB->getBondType())-1;
                  Bond *newB=product->getBondWithIdx(bondIdx);
                  newB->setBondDir(origB->getBondDir());
                }
              }
              nbrIdx++;
            }
          } // end of atomStack traversal

          // now that we've added all the reactant's neighbors, check to see if 
          // it is chiral in the reactant but is not in the reaction. If so
          // we need to worry about its chirality
          productAtom = product->getAtomWithIdx(reactProdAtomMap[reactantAtomIdx]);
          if(productAtom->getChiralTag()==Atom::CHI_UNSPECIFIED &&
             reactantAtom->getChiralTag()!=Atom::CHI_UNSPECIFIED &&
             reactantAtom->getChiralTag()!=Atom::CHI_OTHER &&
             !productAtom->hasProp("molInversionFlag") ){
            // we can only do something sensible here if we have the same number of bonds
            // in the reactants and the products:
            if(reactantAtom->getDegree()==productAtom->getDegree()){
              unsigned int nUnknown=0;
              INT_LIST pOrder;

              ROMol::ADJ_ITER nbrIdx,endNbrs;
              boost::tie(nbrIdx,endNbrs) = product->getAtomNeighbors(productAtom);
              while(nbrIdx!=endNbrs){
                if(prodReactAtomMap.find(*nbrIdx)==prodReactAtomMap.end()){
                  ++nUnknown;
                  // if there's more than one bond in the product that doesn't correspond to
                  // anything in the reactant, we're also doomed
                  if(nUnknown>1) break;

                  // otherwise, add a -1 to the bond order that we'll fill in later
                  pOrder.push_back(-1);
                } else {
                  const Bond *rBond=reactant->getBondBetweenAtoms(reactantAtom->getIdx(),prodReactAtomMap[*nbrIdx]);
                  CHECK_INVARIANT(rBond,"expected reactant bond not found");
                  pOrder.push_back(rBond->getIdx());
                }
                ++nbrIdx;
              }
              if(nUnknown==1){
                // find the reactant bond that hasn't yet been accounted for:
                int unmatchedBond=-1;
                boost::tie(nbrIdx,endNbrs) = reactant->getAtomNeighbors(reactantAtom);
                while(nbrIdx!=endNbrs){
                  const Bond *rBond=reactant->getBondBetweenAtoms(reactantAtom->getIdx(),*nbrIdx);
                  if(std::find(pOrder.begin(),pOrder.end(),rBond->getIdx())==pOrder.end()){
                    unmatchedBond=rBond->getIdx();
                    break;
                  }
                  ++nbrIdx;
                }
                // what must be true at this point:
                //  1) there's a -1 in pOrder that we'll substitute for
                //  2) unmatchedBond contains the index of the substitution
                INT_LIST::iterator bPos=std::find(pOrder.begin(),pOrder.end(),-1);
                if(unmatchedBond>=0 && bPos!=pOrder.end()){
                  *bPos=unmatchedBond;
                }
                if(std::find(pOrder.begin(),pOrder.end(),-1)==pOrder.end()){
                  nUnknown=0;
                }
              }
              if(!nUnknown){
                productAtom->setChiralTag(reactantAtom->getChiralTag());
                int nSwaps = reactantAtom->getPerturbationOrder(pOrder);
                if(nSwaps%2){
                  productAtom->invertChirality();
                }
              }
            }
          }

        }
      } // end of loop over matched atoms

      // ---------- ---------- ---------- ---------- ---------- ---------- 
      // now we need to loop over atoms from the reactants that were chiral but not
      // directly involved in the reaction in order to make sure their chirality hasn't
      // been disturbed
      for(std::vector<const Atom *>::const_iterator atomIt=chiralAtomsToCheck.begin();
          atomIt!=chiralAtomsToCheck.end();++atomIt){
        const Atom *reactantAtom=*atomIt;
        Atom *productAtom=product->getAtomWithIdx(reactProdAtomMap[reactantAtom->getIdx()]);
        CHECK_INVARIANT(reactantAtom->getChiralTag()!=Atom::CHI_UNSPECIFIED,
                        "missing atom chirality.");
        CHECK_INVARIANT(reactantAtom->getChiralTag()==productAtom->getChiralTag(),
                        "invalid product chirality.");
        
        if( reactantAtom->getOwningMol().getAtomDegree(reactantAtom) !=
            product->getAtomDegree(productAtom) ){
          // If the number of bonds to the atom has changed in the course of the
          // reaction we're lost, so remove chirality.
          //  A word of explanation here: the atoms in the chiralAtomsToCheck set are
          //  not explicitly mapped atoms of the reaction, so we really have no idea what
          //  to do with this case. At the moment I'm not even really sure how this
          //  could happen, but better safe than sorry.
          productAtom->setChiralTag(Atom::CHI_UNSPECIFIED);
        } else if(reactantAtom->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW ||
                  reactantAtom->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW){
          // this will contain the indices of product bonds in the
          // reactant order:
          INT_LIST newOrder; 

          ROMol::OEDGE_ITER beg,end;
          boost::tie(beg,end) = reactantAtom->getOwningMol().getAtomBonds(reactantAtom);
          while(beg!=end){
            const BOND_SPTR reactantBond=reactantAtom->getOwningMol()[*beg];
            unsigned int oAtomIdx=reactantBond->getOtherAtomIdx(reactantAtom->getIdx());
            CHECK_INVARIANT(reactProdAtomMap.find(oAtomIdx)!=reactProdAtomMap.end(),
                            "other atom from bond not mapped.");
            const Bond *productBond;
            productBond=product->getBondBetweenAtoms(productAtom->getIdx(),
                                                     reactProdAtomMap[oAtomIdx]);
            CHECK_INVARIANT(productBond,"no matching bond found in product");
            newOrder.push_back(productBond->getIdx());
            ++beg;
          }
          int nSwaps=productAtom->getPerturbationOrder(newOrder);
          if(nSwaps%2){
            productAtom->invertChirality();
          }
        } else {
          // not tetrahedral chirality, don't do anything.
        }
      } // end of loop over chiralAtomsToCheck

      // ---------- ---------- ---------- ---------- ---------- ---------- 
      // finally we may need to set the coordinates in the product conformer:
      if(productConf){
        productConf->resize(product->getNumAtoms());
        if(reactantSptr->getNumConformers()){
          const Conformer &reactConf=reactantSptr->getConformer();
          if(reactConf.is3D()) productConf->set3D(true);
          for(std::map<unsigned int,unsigned int>::const_iterator pr=reactProdAtomMap.begin();
              pr!=reactProdAtomMap.end();++pr){
            productConf->setAtomPos(pr->second,reactConf.getAtomPos(pr->first));
          }
        }
      } // end of conformer update loop

    } // end of addReactantAtomsAndBonds()
  } // End of namespace ReactionUtils
  

  MOL_SPTR_VECT ChemicalReaction::generateOneProductSet(const MOL_SPTR_VECT &reactants,
                                                        const std::vector<MatchVectType> &reactantsMatch) const {
    PRECONDITION(reactants.size()==reactantsMatch.size(),"vector size mismatch");                                                          
    MOL_SPTR_VECT res;
    res.resize(this->getNumProductTemplates());

    // if any of the reactants have a conformer, we'll go ahead and
    // generate conformers for the products:
    bool doConfs=false;
    BOOST_FOREACH(ROMOL_SPTR reactant,reactants){
      if(reactant->getNumConformers()){
        doConfs=true;
        break;
      }
    }
    
    unsigned int prodId=0;
    for(MOL_SPTR_VECT::const_iterator pTemplIt=this->beginProductTemplates();
        pTemplIt!=this->endProductTemplates();++pTemplIt){
      RWMOL_SPTR product=ReactionUtils::initProduct(*pTemplIt);
      Conformer *conf=0;
      if(doConfs){
        conf = new Conformer();
        conf->set3D(false);
      }
      
      for(unsigned int reactantId=0;reactantId<reactants.size();++reactantId){
        ReactionUtils::addReactantAtomsAndBonds(this,
                                                product,reactants[reactantId],
                                                reactantsMatch[reactantId],
                                                this->m_reactantTemplates[reactantId],
                                                conf);  
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
  
  std::vector<MOL_SPTR_VECT> ChemicalReaction::runReactants(const MOL_SPTR_VECT reactants) const {
    if(this->df_needsInit) {
     throw ChemicalReactionException("initMatchers() must be called before runReactants()");
    }
    if(reactants.size() != this->getNumReactantTemplates()){
      throw ChemicalReactionException("Number of reactants provided does not match number of reactant templates.");    
    }
    BOOST_FOREACH(ROMOL_SPTR msptr,reactants){
      CHECK_INVARIANT(msptr,"bad molecule in reactants");
    }
    
    std::vector<MOL_SPTR_VECT> productMols;
    productMols.clear();
    
    // if we have no products, return now:
    if(!this->getNumProductTemplates()){ 
      return productMols;
    }
    
    // find the matches for each reactant:
    VectVectMatchVectType matchesByReactant;
    if(!ReactionUtils::getReactantMatches(reactants,this->m_reactantTemplates,matchesByReactant)){
      // some reactants didn't find a match, return an empty product list:
      return productMols;
    }
    
    // -------------------------------------------------------
    // we now have matches for each reactant, so we can start creating products:
    
    // start by doing the combinatorics on the matches:
    VectVectMatchVectType reactantMatchesPerProduct;
    ReactionUtils::generateReactantCombinations(matchesByReactant,reactantMatchesPerProduct);
    productMols.resize(reactantMatchesPerProduct.size());                                
    
    for(unsigned int productId=0;productId!=productMols.size();++productId){
      MOL_SPTR_VECT lProds=this->generateOneProductSet(reactants,reactantMatchesPerProduct[productId]);
      productMols[productId]=lProds;
    }
    
    return productMols;
  } // end of ChemicalReaction::runReactants()

  ChemicalReaction::ChemicalReaction(const std::string &pickle) {
    ReactionPickler::reactionFromPickle(pickle,this);
  }
  
  void ChemicalReaction::initReactantMatchers() {
    unsigned int nWarnings,nErrors;
    if(!this->validate(nWarnings,nErrors)){
      BOOST_LOG(rdErrorLog)<<"initialization failed\n";
      this->df_needsInit=true;
    } else {
      this->df_needsInit=false;
    }
  }

  bool ChemicalReaction::validate(unsigned int &numWarnings,
                                  unsigned int &numErrors,
                                  bool silent) const {
    bool res=true;
    numWarnings=0;
    numErrors=0;
    
    if(!this->getNumReactantTemplates()){
      if(!silent){
        BOOST_LOG(rdErrorLog)<<"reaction has no reactants\n";
      }
      numErrors++;
      res=false;
    }
    
    if(!this->getNumProductTemplates()){
      if(!silent){
        BOOST_LOG(rdErrorLog)<<"reaction has no products\n";
      }
      numErrors++;
      res=false;
    }
    
    std::vector<int> mapNumbersSeen;
    std::map<int,const Atom *> reactingAtoms;
    unsigned int molIdx=0;
    for(MOL_SPTR_VECT::const_iterator molIter=this->beginReactantTemplates();
        molIter!=this->endReactantTemplates();++molIter){
      bool thisMolMapped=false;
      for(ROMol::AtomIterator atomIt=(*molIter)->beginAtoms();
          atomIt!=(*molIter)->endAtoms();++atomIt){
        if((*atomIt)->hasProp("molAtomMapNumber")){
          thisMolMapped=true;
          int mapNum;
          (*atomIt)->getProp("molAtomMapNumber",mapNum);
          if(std::find(mapNumbersSeen.begin(),mapNumbersSeen.end(),mapNum)!=mapNumbersSeen.end()){
            if(!silent){
              BOOST_LOG(rdErrorLog)<<"reactant atom-mapping number "<<mapNum<<" found multiple times.\n";
            }
            numErrors++;
            res=false;
          } else {
            mapNumbersSeen.push_back(mapNum);
            reactingAtoms[mapNum]=*atomIt;
          }
        }
      }
      if(!thisMolMapped){
        if(!silent){
          BOOST_LOG(rdWarningLog)<<"reactant "<<molIdx<<" has no mapped atoms.\n";
        }
        numWarnings++;
      }
      molIdx++;
    }

    std::vector<int> productNumbersSeen;  
    molIdx=0;
    for(MOL_SPTR_VECT::const_iterator molIter=this->beginProductTemplates();
        molIter!=this->endProductTemplates();++molIter){

      // clear out some possible cached properties to prevent
      // misleading warnings
      for(ROMol::AtomIterator atomIt=(*molIter)->beginAtoms();
          atomIt!=(*molIter)->endAtoms();++atomIt){
        if((*atomIt)->hasProp("_QueryFormalCharge"))
          (*atomIt)->clearProp("_QueryFormalCharge");
        if((*atomIt)->hasProp("_QueryHCount"))
          (*atomIt)->clearProp("_QueryHCount");
        if((*atomIt)->hasProp("_QueryMass"))
          (*atomIt)->clearProp("_QueryMass");
        if((*atomIt)->hasProp("_QueryIsotope"))
          (*atomIt)->clearProp("_QueryIsotope");
      }
      bool thisMolMapped=false;
      for(ROMol::AtomIterator atomIt=(*molIter)->beginAtoms();
          atomIt!=(*molIter)->endAtoms();++atomIt){
        if((*atomIt)->hasProp("molAtomMapNumber")){
          thisMolMapped=true;
          int mapNum;
          (*atomIt)->getProp("molAtomMapNumber",mapNum);
          bool seenAlready=std::find(productNumbersSeen.begin(),
                                     productNumbersSeen.end(),mapNum)!=productNumbersSeen.end();
          if(seenAlready){
            if(!silent){
              BOOST_LOG(rdErrorLog)<<"product atom-mapping number "<<mapNum<<" found multiple times.\n";
            }
            numErrors++;
            res=false;
          } else {
            productNumbersSeen.push_back(mapNum);
          }
          std::vector<int>::iterator ivIt=std::find(mapNumbersSeen.begin(),
                                                    mapNumbersSeen.end(),mapNum);
          if(ivIt==mapNumbersSeen.end()){
            if(!seenAlready){
              if(!silent){
                BOOST_LOG(rdWarningLog)<<"product atom-mapping number "<<mapNum<<" not found in reactants.\n";
              }
              numWarnings++;
              //res=false;
            }
          } else {
            mapNumbersSeen.erase(ivIt); 

            // ------------
            //   The atom is mapped, check to see if its connectivity changes
            // ------------
            const Atom *rAtom=reactingAtoms[mapNum];
            CHECK_INVARIANT(rAtom,"missing atom");
            if(rAtom->getDegree()!=(*atomIt)->getDegree()){
              (*atomIt)->setProp("_ReactionDegreeChanged",1);
            }
          }
        }

        // ------------
        //    Deal with queries
        // ------------
        if((*atomIt)->hasQuery()){
          std::list<const Atom::QUERYATOM_QUERY *>queries;
          queries.push_back((*atomIt)->getQuery());
          while(!queries.empty()){
            const Atom::QUERYATOM_QUERY *query=queries.front();
            queries.pop_front();
            for(Atom::QUERYATOM_QUERY::CHILD_VECT_CI qIter=query->beginChildren();
                qIter!=query->endChildren();++qIter){
              queries.push_back((*qIter).get());
            }
            if(query->getDescription()=="AtomFormalCharge"){
              if((*atomIt)->hasProp("_QueryFormalCharge")){
                if(!silent){
                  BOOST_LOG(rdWarningLog)<<"atom "<<(*atomIt)->getIdx()<<" in product " 
                                         << molIdx << " has multiple charge specifications.\n";
                }
                numWarnings++;
              } else {
                (*atomIt)->setProp("_QueryFormalCharge",
                                   ((const ATOM_EQUALS_QUERY *)query)->getVal());
              }
            } else if(query->getDescription()=="AtomHCount"){
              if((*atomIt)->hasProp("_QueryHCount")){
                if(!silent){
                  BOOST_LOG(rdWarningLog)<<"atom "<<(*atomIt)->getIdx()<<" in product " 
                                         << molIdx << " has multiple H count specifications.\n";
                }
                numWarnings++;
              } else {
                (*atomIt)->setProp("_QueryHCount",
                                   ((const ATOM_EQUALS_QUERY *)query)->getVal());
              }
            } else if(query->getDescription()=="AtomMass"){
              if((*atomIt)->hasProp("_QueryMass")){
                if(!silent) {
                  BOOST_LOG(rdWarningLog)<<"atom "<<(*atomIt)->getIdx()<<" in product " 
                                         << molIdx << " has multiple mass specifications.\n";
                }
                numWarnings++;
              } else {
                (*atomIt)->setProp("_QueryMass",
                                   ((const ATOM_EQUALS_QUERY *)query)->getVal()/massIntegerConversionFactor);
              }
            } else if(query->getDescription()=="AtomIsotope"){
              if((*atomIt)->hasProp("_QueryIsotope")){
                if(!silent) {
                  BOOST_LOG(rdWarningLog)<<"atom "<<(*atomIt)->getIdx()<<" in product " 
                                         << molIdx << " has multiple isotope specifications.\n";
                }
                numWarnings++;
              } else {
                (*atomIt)->setProp("_QueryIsotope",
                                   ((const ATOM_EQUALS_QUERY *)query)->getVal());
              }
            }
          }
        }
      }
      if(!thisMolMapped){
        if(!silent){
          BOOST_LOG(rdWarningLog)<<"product "<<molIdx<<" has no mapped atoms.\n";
        }
        numWarnings++;
      }
      molIdx++;
    }
    if(!mapNumbersSeen.empty()){
      if(!silent){
        std::ostringstream ostr;
        ostr<<"mapped atoms in the reactants were not mapped in the products.\n";
        ostr<<"  unmapped numbers are: ";
        for(std::vector<int>::const_iterator ivIt=mapNumbersSeen.begin();
            ivIt!=mapNumbersSeen.end();++ivIt){
          ostr<< *ivIt << " ";
        }
        ostr<< "\n";
        BOOST_LOG(rdWarningLog)<<ostr.str();
      }
      numWarnings++;
    }
    
    return res;
  }

  bool isMoleculeReactantOfReaction(const ChemicalReaction &rxn,const ROMol &mol,
                                      unsigned int &which){
    if(!rxn.isInitialized()){
     throw ChemicalReactionException("initMatchers() must be called first");
    }
    which=0;
    for(MOL_SPTR_VECT::const_iterator iter=rxn.beginReactantTemplates();
        iter!=rxn.endReactantTemplates();++iter,++which){
      MatchVectType tvect;
      if(SubstructMatch(mol,**iter,tvect)){
        return true;
      }
    }
    return false;
  }

  bool isMoleculeProductOfReaction(const ChemicalReaction &rxn,const ROMol &mol,
                                      unsigned int &which){
    if(!rxn.isInitialized()){
     throw ChemicalReactionException("initMatchers() must be called first");
    }
    which=0;
    for(MOL_SPTR_VECT::const_iterator iter=rxn.beginProductTemplates();
        iter!=rxn.endProductTemplates();++iter,++which){
      MatchVectType tvect;
      if(SubstructMatch(mol,**iter,tvect)){
        return true;
      }
    }
    return false;
  }

  void addRecursiveQueriesToReaction(ChemicalReaction &rxn,
                  const std::map<std::string,ROMOL_SPTR> &queries,
                  std::string propName,
                  std::vector<std::vector<std::pair<unsigned int,std::string> > > *reactantLabels) {
    if(!rxn.isInitialized()){
      throw ChemicalReactionException("initMatchers() must be called first");
    }

    if (reactantLabels!=NULL) {
        (*reactantLabels).resize(0);
    }

    for(MOL_SPTR_VECT::const_iterator rIt=rxn.beginReactantTemplates();
            rIt!=rxn.endReactantTemplates();++rIt){
        if (reactantLabels!=NULL) {
          std::vector<std::pair<unsigned int, std::string> > labels;
          addRecursiveQueries(**rIt, queries, propName, &labels);
          (*reactantLabels).push_back(labels);
        } else {
          addRecursiveQueries(**rIt, queries, propName);
        }
    }
  }

  namespace {
    // recursively looks for atomic number queries anywhere in this set of children
    // or its children
    int numComplexQueries(Queries::Query<int,Atom const *,true>::CHILD_VECT_CI childIt,
                          Queries::Query<int,Atom const *,true>::CHILD_VECT_CI endChildren){
      int res=0;
      while(childIt!=endChildren){
        std::string descr=(*childIt)->getDescription();
        if(descr=="AtomAtomicNum"||descr=="AtomNull"){
          ++res;
        } else {
          res += numComplexQueries((*childIt)->beginChildren(),
                                   (*childIt)->endChildren());
        }
        ++childIt;
      }
      return res;
    }
    // FIX: this is adapted from Fingerprints.cpp and we really should have code
    // like this centralized
    bool isComplexQuery(const Atom &a){
      if( !a.hasQuery()) return false;
      // negated things are always complex:
      if( a.getQuery()->getNegation()) return true;
      std::string descr=a.getQuery()->getDescription();
      if(descr=="AtomAtomicNum") return false;
      if(descr=="AtomOr" || descr=="AtomXor") return true;
      if(descr=="AtomAnd"){
        Queries::Query<int,Atom const *,true>::CHILD_VECT_CI childIt=a.getQuery()->beginChildren();
        int ncq=numComplexQueries(childIt,a.getQuery()->endChildren());
        if(ncq==1){
          return false;
        }
      }
      return true;
    }


    bool isChangedAtom(const Atom &rAtom,const Atom &pAtom,int mapNum,
                       const std::map<int,const Atom *> &mappedProductAtoms) {
      PRECONDITION(mappedProductAtoms.find(mapNum)!=mappedProductAtoms.end(),"atom not mapped in products");

      if(rAtom.getAtomicNum()!=pAtom.getAtomicNum() &&
         pAtom.getAtomicNum()>0 ){
        // the atomic number changed and the product wasn't a dummy
        return true;
      } else if(rAtom.getDegree() != pAtom.getDegree()){
        // the degree changed
        return true;
      } else if(pAtom.getAtomicNum()>0 && isComplexQuery(rAtom)){
        // more than a simple query
        return true;
      }

      // now check bond layout:
      std::map<unsigned int,const Bond *> reactantBonds;
      ROMol::ADJ_ITER nbrIdx,endNbrs;
      boost::tie(nbrIdx,endNbrs) = rAtom.getOwningMol().getAtomNeighbors(&rAtom);
      while(nbrIdx!=endNbrs){
        const ATOM_SPTR nbr=rAtom.getOwningMol()[*nbrIdx];
        if(nbr->hasProp("molAtomMapNumber")){
          int mapNum;
          nbr->getProp("molAtomMapNumber",mapNum);
          reactantBonds[mapNum]=rAtom.getOwningMol().getBondBetweenAtoms(rAtom.getIdx(),
                                                                         nbr->getIdx());
        } else {
          // if we have an un-mapped neighbor, we are automatically a reacting atom:
          return true;
        }
        ++nbrIdx;
      }
      boost::tie(nbrIdx,endNbrs) = pAtom.getOwningMol().getAtomNeighbors(&pAtom);
      while(nbrIdx!=endNbrs){
        const ATOM_SPTR nbr=pAtom.getOwningMol()[*nbrIdx];
        if(nbr->hasProp("molAtomMapNumber")){
          int mapNum;
          nbr->getProp("molAtomMapNumber",mapNum);
          // if we don't have a bond to a similarly mapped atom in the reactant,
          // we're done:
          if(reactantBonds.find(mapNum)==reactantBonds.end()){
            return true;
          }
          const Bond *rBond=reactantBonds[mapNum];
          const Bond *pBond=pAtom.getOwningMol().getBondBetweenAtoms(pAtom.getIdx(),
                                                                     nbr->getIdx());

          // bond comparison logic:
          if(rBond->hasQuery()){
            if(!pBond->hasQuery()) {
              // reactant query, product not query: always a change
              return true;
            } else {
              if( pBond->getQuery()->getDescription()=="BondNull" ){
                // null queries are trump, they match everything
              } else if( rBond->getBondType()==Bond::SINGLE && pBond->getBondType()==Bond::SINGLE &&
                  rBond->getQuery()->getDescription()=="BondOr" &&
                  pBond->getQuery()->getDescription()=="BondOr" ) {
                // The SMARTS parser tags unspecified bonds as single, but then adds
                // a query so that they match single or double.
                // these cases match
              } else {
                if(  rBond->getBondType()==pBond->getBondType() &&
                     rBond->getQuery()->getDescription()=="BondOrder" &&
                     pBond->getQuery()->getDescription()=="BondOrder" &&
                     static_cast<BOND_EQUALS_QUERY *>(rBond->getQuery())->getVal()==
                     static_cast<BOND_EQUALS_QUERY *>(pBond->getQuery())->getVal()
                     ) {
                  // bond order queries with equal orders also match
                } else {
                  // anything else does not match
                  return true;
                }
              }
            }
          } else if(pBond->hasQuery()) {
            // reactant not query, product query
            // if product is anything other than the null query
            // it's a change:
            if(pBond->getQuery()->getDescription()!="BondNull"){
              return true;
            }
            
          } else {
            // neither has a query, just compare the types
            if(rBond->getBondType()!=pBond->getBondType()){
              return true;
            }
          }
        }
        ++nbrIdx;
      }

      // haven't found anything to say that we are changed, so we must
      // not be
      return false;
    }

    template <class T>
    bool getMappedAtoms(T &rIt,
                        std::map<int,const Atom *> &mappedAtoms){
      ROMol::ATOM_ITER_PAIR atItP = rIt->getVertices();
      while(atItP.first != atItP.second ){
        const Atom *oAtom=(*rIt)[*(atItP.first++)].get();
        // we only worry about mapped atoms:
        if(oAtom->hasProp("molAtomMapNumber")){
          int mapNum;
          oAtom->getProp("molAtomMapNumber",mapNum);
          mappedAtoms[mapNum]=oAtom;
        } else {
          // unmapped atom, return it
          return false;
        }
      }
      return true;
    }
  } // end of anonymous namespace
  
  VECT_INT_VECT getReactingAtoms(const ChemicalReaction &rxn,bool mappedAtomsOnly){
    if(!rxn.isInitialized()){
     throw ChemicalReactionException("initMatchers() must be called first");
    }
    VECT_INT_VECT res;
    res.resize(rxn.getNumReactantTemplates());

    // find mapped atoms in the products :
    std::map<int,const Atom *> mappedProductAtoms;
    for(MOL_SPTR_VECT::const_iterator rIt=rxn.beginProductTemplates();
        rIt!=rxn.endProductTemplates();++rIt){
      getMappedAtoms(*rIt,mappedProductAtoms);
    }

    // now loop over mapped atoms in the reactants, keeping track of
    // which reactant they are associated with, and check for changes.
    VECT_INT_VECT::iterator resIt=res.begin();
    for(MOL_SPTR_VECT::const_iterator rIt=rxn.beginReactantTemplates();
        rIt!=rxn.endReactantTemplates();++rIt,++resIt){
      ROMol::ATOM_ITER_PAIR atItP = (*rIt)->getVertices();
      while(atItP.first != atItP.second ){
        const Atom *oAtom=(**rIt)[*(atItP.first++)].get();
        // unmapped atoms are definitely changing:
        if(!oAtom->hasProp("molAtomMapNumber")){
          if(!mappedAtomsOnly){
            resIt->push_back(oAtom->getIdx());
          }
        } else {
          // but mapped ones require more careful consideration
          int mapNum;
          oAtom->getProp("molAtomMapNumber",mapNum);
          // if this is found in a reactant:
          if(mappedProductAtoms.find(mapNum)!=mappedProductAtoms.end()){
            if( isChangedAtom(*oAtom,*(mappedProductAtoms[mapNum]),
                              mapNum,mappedProductAtoms) ){
              resIt->push_back(oAtom->getIdx());
            }
          }
        }
      }
    }
    return res;
  }

} // end of RDKit namespace
