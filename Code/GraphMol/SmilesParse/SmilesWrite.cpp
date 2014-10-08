// $Id$
//
//  Copyright (C) 2002-2012 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "SmilesWrite.h"
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>
#include <GraphMol/Canon.h>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/dynamic_bitset.hpp>
#include <sstream>
#include <map>
#include <list>

//#define VERBOSE_CANON 1

namespace RDKit{

  namespace SmilesWrite{
    const int atomicSmiles[] = {5,6,7,8,9,15,16,17,35,53,-1};
    bool inOrganicSubset(int atomicNumber){
      unsigned int idx=0;
      while( atomicSmiles[idx]<atomicNumber &&
             atomicSmiles[idx]!=-1){
        ++idx;
      }
      if(atomicSmiles[idx]==atomicNumber){
        return true;
      }
      return false;
    }


    std::string GetAtomSmiles(const Atom *atom,bool doKekule,const Bond *bondIn){
      PRECONDITION(atom,"bad atom");
      INT_VECT atomicSmilesVect(atomicSmiles,
                                atomicSmiles+(sizeof(atomicSmiles)-1)/sizeof(atomicSmiles[0]));
      std::stringstream res;
      int fc = atom->getFormalCharge();
      int num = atom->getAtomicNum();
      int isotope = atom->getIsotope();

      bool needsBracket=false;
      std::string symb;
      if(atom->hasProp("smilesSymbol")){
        atom->getProp("smilesSymbol",symb);
      } else {
        symb=PeriodicTable::getTable()->getElementSymbol(num);
      }
      //symb = atom->getSymbol();
      if(inOrganicSubset(num)){
        // it's a member of the organic subset
        //if(!doKekule && atom->getIsAromatic() && symb[0] < 'a') symb[0] -= ('A'-'a');

        // -----
        // figure out if we need to put a bracket around the atom,
        // the conditions for this are:
        //   - formal charge specified
        //   - the atom has a nonstandard valence
        //   - chirality present and writing isomeric smiles
        //   - non-default isotope and writing isomeric smiles
        //   - atom-map information present
        const INT_VECT &defaultVs=PeriodicTable::getTable()->getValenceList(num);
        int totalValence= atom->getTotalValence();
        bool nonStandard=false;

        if(atom->getNumRadicalElectrons()){
          nonStandard=true;
        } else if((num==7||num==15) && atom->getIsAromatic() && atom->getNumExplicitHs()){
          // another type of "nonstandard" valence is an aromatic N or P with
          // explicit Hs indicated:
          nonStandard=true;
        } else {
          nonStandard = (totalValence!=defaultVs.front() && atom->getTotalNumHs());
        }

        if(fc || nonStandard){
          needsBracket=true;
        }
        if(atom->getOwningMol().hasProp("_doIsoSmiles")){
          if( atom->getChiralTag()!=Atom::CHI_UNSPECIFIED ){
            needsBracket = true;
          } else if(isotope){
            needsBracket=true;
          }
        }
        if(atom->hasProp("molAtomMapNumber")){
          needsBracket=true;
        }
      } else {
        needsBracket = true;
      }
      if( needsBracket ) res << "[";

      if(isotope && atom->getOwningMol().hasProp("_doIsoSmiles")){
        res <<isotope;
      }
      // this was originally only done for the organic subset,
      // applying it to other atom-types is a fix for Issue 3152751: 
      if(!doKekule && atom->getIsAromatic() && symb[0]>='A' && symb[0] <= 'Z'){
        symb[0] -= ('A'-'a');
      }
      res << symb;

      if(atom->getOwningMol().hasProp("_doIsoSmiles") &&
         atom->getChiralTag()!=Atom::CHI_UNSPECIFIED ){
        INT_LIST trueOrder;
        atom->getProp("_TraversalBondIndexOrder",trueOrder);
        int nSwaps=  atom->getPerturbationOrder(trueOrder);
        if(atom->getDegree()==3 && !bondIn){
          // This is a special case. Here's an example:
          //   Our internal representation of a chiral center is equivalent to:
          //     [C@](F)(O)(C)[H]
          //   we'll be dumping it without the H, which entails a reordering:
          //     [C@@H](F)(O)C
          ++nSwaps;
        }
        //BOOST_LOG(rdErrorLog)<<">>>> "<<atom->getIdx()<<" "<<nSwaps<<" "<<atom->getChiralTag()<<std::endl;
        std::string atStr="";
        switch(atom->getChiralTag()){
        case Atom::CHI_TETRAHEDRAL_CW:
          if(!(nSwaps%2))
            atStr = "@@";
          else
            atStr = "@";
          break;
        case Atom::CHI_TETRAHEDRAL_CCW:
          if(!(nSwaps%2))
            atStr = "@";
          else
            atStr = "@@";
          break;
        default:
          break;
        }
        res << atStr;
      }

      if(needsBracket){
        unsigned int totNumHs=atom->getTotalNumHs();
        if(totNumHs > 0){
          res << "H";
          if(totNumHs > 1) res << totNumHs;
        }
        if(fc > 0){
          res << "+";
          if(fc > 1) res << fc;
        } else if(fc < 0) {
          res << "-";
          if(fc < -1) res << -fc;
        }
    
        if(atom->hasProp("molAtomMapNumber")){
          int mapNum;
          atom->getProp("molAtomMapNumber",mapNum);
          res<<":"<<mapNum;
        }
        res << "]";
      }

      // If the atom has this property, the contained string will
      // be inserted directly in the SMILES:
      if(atom->hasProp("_supplementalSmilesLabel")){
        std::string label;
        atom->getProp("_supplementalSmilesLabel",label);
        res << label;
      }

      return res.str();
    }

    std::string GetBondSmiles(const Bond *bond,int atomToLeftIdx,bool doKekule,bool allBondsExplicit){
      PRECONDITION(bond,"bad bond");
      if(atomToLeftIdx<0) atomToLeftIdx=bond->getBeginAtomIdx();

      std::stringstream res;
      bool aromatic=false;
      if( !doKekule &&
          (bond->getBondType() == Bond::SINGLE ||
           bond->getBondType() == Bond::DOUBLE ||
           bond->getBondType() == Bond::AROMATIC) ){
        Atom *a1,*a2;
        a1 = bond->getOwningMol().getAtomWithIdx(atomToLeftIdx);
        a2 = bond->getOwningMol().getAtomWithIdx(bond->getOtherAtomIdx(atomToLeftIdx));
        if((a1->getIsAromatic() && a2->getIsAromatic()) &&
           (a1->getAtomicNum()||a2->getAtomicNum())) aromatic=true;
      }

      Bond::BondDir dir= bond->getBondDir();

      if(bond->hasProp("_TraversalRingClosureBond")){
        //std::cerr<<"FLIP: "<<bond->getIdx()<<" "<<bond->getBeginAtomIdx()<<"-"<<bond->getEndAtomIdx()<<std::endl;
        //if(dir==Bond::ENDDOWNRIGHT) dir=Bond::ENDUPRIGHT;
        //else if(dir==Bond::ENDUPRIGHT) dir=Bond::ENDDOWNRIGHT;
        bond->clearProp("_TraversalRingClosureBond");
      }
  
      switch(bond->getBondType()){
      case Bond::SINGLE:
        if( dir != Bond::NONE && dir != Bond::UNKNOWN ){
          switch(dir){
          case Bond::ENDDOWNRIGHT:
            if(bond->getOwningMol().hasProp("_doIsoSmiles"))  res << "\\";
            break;
          case Bond::ENDUPRIGHT:
            if(bond->getOwningMol().hasProp("_doIsoSmiles"))  res << "/";
            break;
          default:
            break;
          }
        } else {
          // if the bond is marked as aromatic and the two atoms
          //  are aromatic, we need no marker (this arises in kekulized
          //  molecules).
          // FIX: we should be able to dump kekulized smiles
          //   currently this is possible by removing all
          //   isAromatic flags, but there should maybe be another way
          if(allBondsExplicit) res<<"-";
          else if( aromatic && !bond->getIsAromatic() ) res << "-";
        }
        break;
      case Bond::DOUBLE:
        // see note above
        if( !aromatic || !bond->getIsAromatic() ) res << "=";
        break;
      case Bond::TRIPLE: res << "#"; break;
      case Bond::AROMATIC:
        if ( dir != Bond::NONE && dir != Bond::UNKNOWN ){
          switch(dir){
          case Bond::ENDDOWNRIGHT:
            if(bond->getOwningMol().hasProp("_doIsoSmiles"))  res << "\\";
            break;
          case Bond::ENDUPRIGHT:
            if(bond->getOwningMol().hasProp("_doIsoSmiles"))  res << "/";
            break;
          default:
            break;
          }
        } else if(allBondsExplicit || !aromatic ){
          res << ":";
        }
        break;
      case Bond::DATIVE:
        if(atomToLeftIdx>=0 &&
           bond->getBeginAtomIdx()==static_cast<unsigned int>(atomToLeftIdx) ) res << ">";
        else res << "<";
        break;
      default:
        res << "~";
      }
      return res.str();
    }

    std::string FragmentSmilesConstruct(ROMol &mol,int atomIdx,
                                        std::vector<Canon::AtomColors> &colors,
                                        INT_VECT &ranks,bool doKekule,bool canonical,
                                        bool allBondsExplicit,
                                        std::vector<unsigned int> &atomOrdering,
                                        const boost::dynamic_bitset<> *bondsInPlay=0,
                                        const std::vector<std::string> *atomSymbols=0,
                                        const std::vector<std::string> *bondSymbols=0
                                        ){
      PRECONDITION(!bondsInPlay||bondsInPlay->size()>=mol.getNumBonds(),"bad bondsInPlay");
      PRECONDITION(!atomSymbols||atomSymbols->size()>=mol.getNumAtoms(),"bad atomSymbols");
      PRECONDITION(!bondSymbols||bondSymbols->size()>=mol.getNumBonds(),"bad bondSymbols");
      Canon::MolStack molStack;
      // try to prevent excessive reallocation
      molStack.reserve(mol.getNumAtoms()+
                       mol.getNumBonds());
      std::stringstream res;

      std::map<int,int> ringClosureMap;
      int ringIdx,closureVal;
      if(!canonical) mol.setProp("_StereochemDone",1);
      std::list<unsigned int> ringClosuresToErase;

      Canon::canonicalizeFragment(mol,atomIdx,colors,ranks,
                                  molStack,bondsInPlay,bondSymbols);
      Bond *bond=0;
      BOOST_FOREACH(Canon::MolStackElem mSE,molStack){
        switch(mSE.type){
        case Canon::MOL_STACK_ATOM:
          if(!ringClosuresToErase.empty()){
            BOOST_FOREACH(unsigned int rclosure,ringClosuresToErase){
              ringClosureMap.erase(rclosure);
            }
            ringClosuresToErase.clear();
          }
          //std::cout<<"\t\tAtom: "<<mSE.obj.atom->getIdx()<<std::endl;
          if(!atomSymbols){
            res << GetAtomSmiles(mSE.obj.atom,doKekule,bond);
          } else {
            res << (*atomSymbols)[mSE.obj.atom->getIdx()];
          }
          atomOrdering.push_back(mSE.obj.atom->getIdx());
          
          break;
        case Canon::MOL_STACK_BOND:
          bond = mSE.obj.bond;
          //std::cout<<"\t\tBond: "<<bond->getIdx()<<std::endl;
          if(!bondSymbols){
            res << GetBondSmiles(bond,mSE.number,doKekule,allBondsExplicit);
          } else {
            res << (*bondSymbols)[bond->getIdx()];
          }
          break;
        case Canon::MOL_STACK_RING:
          ringIdx = mSE.number;
          //std::cout<<"\t\tRing: "<<ringIdx;
          if(ringClosureMap.count(ringIdx)){
            // the index is already in the map ->
            //   we're closing a ring, so grab
            //   the index and then delete the value:
            closureVal = ringClosureMap[ringIdx];
            //ringClosureMap.erase(ringIdx);
            ringClosuresToErase.push_back(ringIdx);
          } else {
            // we're opening a new ring, find the index for it:
            closureVal = 1;
            bool done=false;
            // EFF: there's got to be a more efficient way to do this
            while(!done){
              std::map<int,int>::iterator mapIt;
              for(mapIt=ringClosureMap.begin();
                  mapIt!=ringClosureMap.end();
                  mapIt++){
                if(mapIt->second==closureVal) break;
              }
              if(mapIt==ringClosureMap.end()){
                done=true;
              } else {
                closureVal+=1;
              }
            }
            ringClosureMap[ringIdx]=closureVal;
          }
          if(closureVal >= 10){
            res << "%";
          }
          //std::cerr << " > " << closureVal <<std::endl;
          res << closureVal;
          break;
        case Canon::MOL_STACK_BRANCH_OPEN:
          res << "(";
          break;
        case Canon::MOL_STACK_BRANCH_CLOSE:
          res << ")";
          break;
        default:
          break;
        }
      }
      return res.str();
    }

  } // end of namespace SmilesWrite


  std::string MolToSmiles(const ROMol &mol,bool doIsomericSmiles,
                          bool doKekule,int rootedAtAtom,bool canonical,
                          bool allBondsExplicit){
    if(!mol.getNumAtoms()) return "";
    PRECONDITION(rootedAtAtom<0||static_cast<unsigned int>(rootedAtAtom)<mol.getNumAtoms(),
                 "rootedAtomAtom must be less than the number of atoms");

    ROMol tmol(mol,true);
    if(doIsomericSmiles){
      tmol.setProp("_doIsoSmiles",1);
    }
#if 0
    std::cout << "----------------------------" << std::endl;
    std::cout << "MolToSmiles:"<< std::endl;
    tmol.debugMol(std::cout);
    std::cout << "----------------------------" << std::endl;
#endif  
    std::string res;

    for(ROMol::AtomIterator atIt=tmol.beginAtoms();atIt!=tmol.endAtoms();atIt++){
      (*atIt)->updatePropertyCache(false);
    }

    unsigned int nAtoms=tmol.getNumAtoms();
    INT_VECT ranks(nAtoms,-1);

    std::vector<unsigned int> atomOrdering;

    // clean up the chirality on any atom that is marked as chiral,
    // but that should not be:
    if(doIsomericSmiles){
      if(!mol.hasProp("_StereochemDone")){
        MolOps::assignStereochemistry(tmol,true);
      } else {
        tmol.setProp("_StereochemDone",1);
        // we need the CIP codes:
        for(unsigned int aidx=0;aidx<tmol.getNumAtoms();++aidx){
          const Atom *oAt=mol.getAtomWithIdx(aidx);
          if(oAt->hasProp("_CIPCode")){
            std::string cipCode;
            oAt->getProp("_CIPCode",cipCode);
            tmol.getAtomWithIdx(aidx)->setProp("_CIPCode",cipCode);
          }
        }
      }
    }
    if(canonical){
      MolOps::rankAtoms(tmol,ranks,true,doIsomericSmiles,doIsomericSmiles);
    } else {
      for(unsigned int i=0;i<tmol.getNumAtoms();++i) ranks[i]=i;
    }
#ifdef VERBOSE_CANON
    for(unsigned int tmpI=0;tmpI<ranks.size();tmpI++){
      std::cout << tmpI << " " << ranks[tmpI] << " " << *(tmol.getAtomWithIdx(tmpI)) << std::endl;
    }
#endif

    std::vector<Canon::AtomColors> colors(nAtoms,Canon::WHITE_NODE);
    std::vector<Canon::AtomColors>::iterator colorIt;
    colorIt = colors.begin();
    // loop to deal with the possibility that there might be disconnected fragments
    while(colorIt != colors.end()){
      int nextAtomIdx=-1;
      std::string subSmi;

      // find the next atom for a traverse
      if(rootedAtAtom>=0){
        nextAtomIdx=rootedAtAtom;
        rootedAtAtom=-1;
      } else {
        int nextRank = nAtoms+1;
        for(unsigned int i=0;i<nAtoms;i++){
          if( colors[i] == Canon::WHITE_NODE && ranks[i] < nextRank ){
            nextRank = ranks[i];
            nextAtomIdx = i;
          }
        }
      }
      CHECK_INVARIANT(nextAtomIdx>=0,"no start atom found");

      subSmi = SmilesWrite::FragmentSmilesConstruct(tmol, nextAtomIdx, colors,
                                                    ranks,doKekule,canonical,allBondsExplicit,
                                                    atomOrdering);

      res += subSmi;
      colorIt = std::find(colors.begin(),colors.end(),Canon::WHITE_NODE);
      if(colorIt != colors.end()){
        res += ".";
      }
    }
    mol.setProp("_smilesAtomOutputOrder",atomOrdering,true);
    return res;
  } // end of MolToSmiles()

  std::string MolFragmentToSmiles(const ROMol &mol,
                                  const std::vector<int> &atomsToUse,
                                  const std::vector<int> *bondsToUse,
                                  const std::vector<std::string> *atomSymbols,
                                  const std::vector<std::string> *bondSymbols,
                                  bool doIsomericSmiles,
                                  bool doKekule,
                                  int rootedAtAtom,
                                  bool canonical,
                                  bool allBondsExplicit){
    PRECONDITION(atomsToUse.size(),
                 "no atoms provided");
    PRECONDITION(rootedAtAtom<0||static_cast<unsigned int>(rootedAtAtom)<mol.getNumAtoms(),
                 "rootedAtomAtom must be less than the number of atoms");
    PRECONDITION(rootedAtAtom<0||std::find(atomsToUse.begin(),atomsToUse.end(),rootedAtAtom)!=atomsToUse.end(),
                 "rootedAtomAtom not found in atomsToUse");
    PRECONDITION(!atomSymbols || atomSymbols->size()>=mol.getNumAtoms(),
                 "bad atomSymbols vector");
    PRECONDITION(!bondSymbols || bondSymbols->size()>=mol.getNumBonds(),
                 "bad bondSymbols vector");
    if(!mol.getNumAtoms()) return "";

    ROMol tmol(mol,true);
    if(doIsomericSmiles){
      tmol.setProp("_doIsoSmiles",1);
    }
    std::string res;

    boost::dynamic_bitset<> atomsInPlay(mol.getNumAtoms(),0);
    BOOST_FOREACH(int aidx,atomsToUse){
      atomsInPlay.set(aidx);
    }
    // figure out which bonds are actually in play:
    boost::dynamic_bitset<> bondsInPlay(mol.getNumBonds(),0);
    if(bondsToUse){
      BOOST_FOREACH(int bidx,*bondsToUse){
        bondsInPlay.set(bidx);
      }
    } else {
      BOOST_FOREACH(int aidx,atomsToUse){
        ROMol::OEDGE_ITER beg,end;
        boost::tie(beg,end) = mol.getAtomBonds(mol.getAtomWithIdx(aidx));
        while(beg!=end){
          const BOND_SPTR bond=mol[*beg];
          if(atomsInPlay[bond->getOtherAtomIdx(aidx)])
            bondsInPlay.set(bond->getIdx());
          ++beg;
        }
        
      }
    }

    // copy over the rings that only involve atoms/bonds in this fragment:
    tmol.getRingInfo()->reset();
    tmol.getRingInfo()->initialize();
    for(unsigned int ridx=0;ridx<mol.getRingInfo()->numRings();++ridx){
      const INT_VECT &aring=mol.getRingInfo()->atomRings()[ridx];
      const INT_VECT &bring=mol.getRingInfo()->bondRings()[ridx];
      bool keepIt=true;
      BOOST_FOREACH(int aidx,aring){
        if(!atomsInPlay[aidx]){
          keepIt=false;
          break;
        }
      }
      if(keepIt){
        BOOST_FOREACH(int bidx,bring){
          if(!bondsInPlay[bidx]){
            keepIt=false;
            break;
          }
        }
      }
      if(keepIt){
        tmol.getRingInfo()->addRing(aring,bring);
      }
    }
    
    for(ROMol::AtomIterator atIt=tmol.beginAtoms();atIt!=tmol.endAtoms();atIt++){
      (*atIt)->updatePropertyCache(false);
    }

    INT_VECT ranks(tmol.getNumAtoms(),-1);

    std::vector<unsigned int> atomOrdering;

    // clean up the chirality on any atom that is marked as chiral,
    // but that should not be:
    if(doIsomericSmiles){
      if(!mol.hasProp("_StereochemDone")){
        MolOps::assignStereochemistry(tmol,true);
      } else {
        tmol.setProp("_StereochemDone",1);
        // we need the CIP codes:
        BOOST_FOREACH(int aidx,atomsToUse){
          const Atom *oAt=mol.getAtomWithIdx(aidx);
          if(oAt->hasProp("_CIPCode")){
            std::string cipCode;
            oAt->getProp("_CIPCode",cipCode);
            tmol.getAtomWithIdx(aidx)->setProp("_CIPCode",cipCode);
          }
        }
      }
    }
    if(canonical){
      MolOps::rankAtomsInFragment(tmol,ranks,atomsInPlay,bondsInPlay,atomSymbols,bondSymbols);
    } else {
      for(unsigned int i=0;i<tmol.getNumAtoms();++i) ranks[i]=i;
    }
#ifdef VERBOSE_CANON
    for(unsigned int tmpI=0;tmpI<ranks.size();tmpI++){
      std::cout << tmpI << " " << ranks[tmpI] << " " << *(tmol.getAtomWithIdx(tmpI)) << std::endl;
    }
#endif

    std::vector<Canon::AtomColors> colors(tmol.getNumAtoms(),Canon::BLACK_NODE);
    BOOST_FOREACH(int aidx,atomsToUse){
      colors[aidx]=Canon::WHITE_NODE;
    }
    std::vector<Canon::AtomColors>::iterator colorIt;
    colorIt = colors.begin();
    // loop to deal with the possibility that there might be disconnected fragments
    while(colorIt != colors.end()){
      int nextAtomIdx=-1;
      std::string subSmi;

      // find the next atom for a traverse
      if(rootedAtAtom>=0){
        nextAtomIdx=rootedAtAtom;
        rootedAtAtom=-1;
      } else {
        int nextRank = tmol.getNumAtoms()+1;
        BOOST_FOREACH(int i,atomsToUse){
          if( colors[i] == Canon::WHITE_NODE && ranks[i] < nextRank ){
            nextRank = ranks[i];
            nextAtomIdx = i;
          }
        }
      }
      CHECK_INVARIANT(nextAtomIdx>=0,"no start atom found");

      subSmi = SmilesWrite::FragmentSmilesConstruct(tmol, nextAtomIdx, colors,
                                                    ranks,doKekule,canonical,allBondsExplicit,
                                                    atomOrdering,
                                                    &bondsInPlay,
                                                    atomSymbols,bondSymbols);

      res += subSmi;
      colorIt = std::find(colors.begin(),colors.end(),Canon::WHITE_NODE);
      if(colorIt != colors.end()){
        res += ".";
      }
    }
    mol.setProp("_smilesAtomOutputOrder",atomOrdering,true);
    return res;
  } // end of MolFragmentToSmiles()
}
