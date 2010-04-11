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
// Created by Greg Landrum, September 2006
//
#ifndef __RD_SLNPARSEOPS_H__
#define __RD_SLNPARSEOPS_H__

#include <vector>
#include <GraphMol/SLNParse/SLNParse.h>
#include <GraphMol/SLNParse/SLNAttribs.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <boost/lexical_cast.hpp>

namespace RDKit{
  namespace SLNParse{
    namespace {
      //!  set a bookmark in the molecule if the atom has an associated ID:
      void bookmarkAtomID(RWMol *mp,Atom *atom){
        PRECONDITION(mp,"bad molecule");
        PRECONDITION(atom,"bad atom");
        if(atom->hasProp("_AtomID")){
          unsigned int label;
          atom->getProp("_AtomID",label);
          if(mp->hasAtomBookmark(label)){
            std::stringstream err;
            err << "SLN Parser error: Atom ID " << label << " used a second time.";
            throw SLNParseException(err.str());
          }
          if(mp->hasBondBookmark(label)){
            std::stringstream err;
            err << "SLN Parser error: Atom ID " << label << " appears *after* its ring closure.";
            throw SLNParseException(err.str());
          }
          mp->setAtomBookmark(atom,label);
        }
      }

      //! adds a bond, being careful to handle aromaticity properly
      template<typename BondType>
      void addBondToMol(RWMol *mp,BondType *bond){
        PRECONDITION(mp,"null molecule");
        PRECONDITION(bond,"null bond");
        mp->addBond(bond,true);
        if(bond->getBondType()==Bond::AROMATIC){
          // SLN doesn't have aromatic atom types, aromaticity is a property
          // of the bonds themselves, so we need to set the atom types:
          bond->setIsAromatic(true);
          bond->getBeginAtom()->setIsAromatic(true);
          bond->getEndAtom()->setIsAromatic(true);
        }
      }
    }// end of anonymous namespace

    // ------------------------------------------------------------------------------------
    //! initialize a molecule
    template <typename AtomType>
    int startMol(std::vector<RWMol *> &molList,AtomType *firstAtom,bool doingQuery){
      PRECONDITION(firstAtom,"empty atom");
      RWMol *mp = new RWMol();
      mp->addAtom(firstAtom,true,true);
      bookmarkAtomID(mp,firstAtom);

      if(!doingQuery){
        // add any hydrogens that are set on the atom, otherwise getting the numbering right
        // is just too hard:
        for(unsigned int i=0;i<firstAtom->getNumExplicitHs();++i){
          int hIdx=mp->addAtom(new Atom(1),false,true);
          mp->addBond(0,hIdx,Bond::SINGLE);
        }
        firstAtom->setNumExplicitHs(0);
      }

      int sz = molList.size();
      molList.push_back(mp);
      return sz;
    };

    // ------------------------------------------------------------------------------------
    //! adds an atom to a molecule
    template<typename AtomType,typename BondType>
    void addAtomToMol(std::vector<RWMol *> &molList,unsigned int idx,AtomType *atom,
                      BondType *bond,bool doingQuery){
      PRECONDITION(idx<molList.size(),"bad index");
      RWMol *mp=molList[idx];
      PRECONDITION(mp,"null molecule");
      PRECONDITION(atom,"empty atom");
      PRECONDITION(bond,"null bond");

      Atom *a1 = mp->getActiveAtom();
      int atomIdx1=a1->getIdx();
      int atomIdx2=mp->addAtom(atom,true,true);
      bookmarkAtomID(mp,atom);
      bond->setOwningMol(mp);
      bond->setBeginAtomIdx(atomIdx1);
      bond->setEndAtomIdx(atomIdx2);
      addBondToMol(mp,bond);
    
      if(!doingQuery){
        // add any hydrogens that are set on the atom, otherwise getting the numbering right
        // is just too hard:
        for(unsigned int i=0;i<atom->getNumExplicitHs();++i){
          int hIdx=mp->addAtom(new Atom(1),false,true);
          mp->addBond(atomIdx2,hIdx,Bond::SINGLE);
        }
        atom->setNumExplicitHs(0);
      }      
    }
    //! \overload
    template<typename AtomType>
    void addAtomToMol(std::vector<RWMol *> &molList,unsigned int idx,AtomType *atom,bool doingQuery){
      addAtomToMol(molList,idx,atom,new Bond(Bond::SINGLE),doingQuery);
    }

    // ------------------------------------------------------------------------------------
    //! closes an indexed ring in a molecule using the bond provided
    // The bond is formed from the atom in the molecule with the
    // corresponding bookmark to the active atom
    //
    template <typename BondType>
    void closeRingBond(std::vector<RWMol *> &molList,unsigned int molIdx,
                       unsigned int ringIdx,BondType *bond,
                       bool postponeAllowed=true){
      PRECONDITION(molIdx<molList.size(),"bad index");
      RWMol *mp=molList[molIdx];
      PRECONDITION(mp,"null molecule");
      PRECONDITION(bond,"Null bond");

      if(!mp->hasAtomBookmark(ringIdx)){
        if(postponeAllowed){
          // save it for later:
          bond->setOwningMol(mp);
          bond->setEndAtomIdx(mp->getActiveAtom()->getIdx());
          mp->setBondBookmark(bond,ringIdx);
          return;
        } else {
          std::stringstream err;
          err << "SLN Parser error: Ring closure " << ringIdx << " does not have a corresponding opener.";
          throw SLNParseException(err.str());
        }
      }
      Atom *opener=mp->getAtomWithBookmark(ringIdx);
      CHECK_INVARIANT(opener,"invalid atom");

      Atom *closer=mp->getActiveAtom();
      bond->setOwningMol(mp);
      bond->setBeginAtom(opener);
      bond->setEndAtom(closer);
      addBondToMol(mp,bond);
    };
    //! \overload
    void closeRingBond(std::vector<RWMol *> &molList,unsigned int molIdx,unsigned int ringIdx){
      closeRingBond(molList,molIdx,ringIdx,new Bond(Bond::SINGLE));
    };

    // ------------------------------------------------------------------------------------
    // NOTE: this takes over responsibility for the bond
    template <typename BondType>
    int addBranchToMol(std::vector<RWMol *> &molList,unsigned int molIdx,
                       unsigned int branchIdx,BondType *&bond){
      PRECONDITION(molIdx<molList.size(),"bad index");
      RWMol *mp=molList[molIdx];
      PRECONDITION(mp,"null molecule");
      PRECONDITION(branchIdx<molList.size(),"bad index");
      RWMol *branch=molList[branchIdx];
      PRECONDITION(branch,"null branch");
      PRECONDITION(bond,"null bond");

      unsigned int activeAtomIdx=mp->getActiveAtom()->getIdx();
      unsigned int nOrigAtoms=mp->getNumAtoms();

      //
      // Add the fragment's atoms and bonds to the molecule:
      //
      mp->insertMol(*branch);

      // copy in any atom bookmarks from the branch:
      for(ROMol::ATOM_BOOKMARK_MAP::const_iterator bmIt=branch->getAtomBookmarks()->begin();
          bmIt != branch->getAtomBookmarks()->end();++bmIt){
        if(bmIt->first<0) continue;
        if(mp->hasAtomBookmark(bmIt->first)){
          std::stringstream err;
          err << "SLN Parser error: Atom ID " << bmIt->first << " used a second time.";
          throw SLNParseException(err.str());
        } else if(mp->hasBondBookmark(bmIt->first)){
            std::stringstream err;
            err << "SLN Parser error: Atom ID " << bmIt->first << " appears *after* its ring closure.";
            throw SLNParseException(err.str());
          }
        else {
          CHECK_INVARIANT(bmIt->second.size()==1,"bad atom bookmark list on branch");
          Atom *tgtAtom=mp->getAtomWithIdx((*bmIt->second.begin())->getIdx()+nOrigAtoms);
          mp->setAtomBookmark(tgtAtom,bmIt->first);
        }
      }
      
      // loop over bond bookmarks in the branch and close the corresponding rings
      for(ROMol::BOND_BOOKMARK_MAP::const_iterator bmIt=branch->getBondBookmarks()->begin();
          bmIt != branch->getBondBookmarks()->end();++bmIt){
        CHECK_INVARIANT(bmIt->second.size()>=1,"bad bond bookmark list on branch");
        for(ROMol::BOND_PTR_LIST::const_iterator bondIt=bmIt->second.begin();
            bondIt!=bmIt->second.end();++bondIt){
          Bond *tgtBond=*bondIt;
          if(bmIt->first>0 && mp->hasAtomBookmark(bmIt->first)){
            Atom *tmpAtom=mp->getActiveAtom();
            mp->setActiveAtom(mp->getAtomWithIdx(tgtBond->getEndAtomIdx()+nOrigAtoms));
            closeRingBond(molList,molIdx,bmIt->first,tgtBond,false);
            mp->setActiveAtom(tmpAtom);
          } else {
            // no partner found yet, copy into this mol:
            tgtBond->setOwningMol(mp);
            tgtBond->setEndAtomIdx(tgtBond->getEndAtomIdx()+nOrigAtoms);
            mp->setBondBookmark(tgtBond,bmIt->first);
          }
        }
      }
      
      // set the connecting bond:
      if(bond->getBondType()!=Bond::IONIC){
        bond->setOwningMol(mp);
        bond->setBeginAtomIdx(activeAtomIdx);
        bond->setEndAtomIdx(nOrigAtoms);
        addBondToMol(mp,bond);
      } else {
        delete bond;
      }
      bond=0;
        
      delete branch;
      unsigned int sz = molList.size();
      if ( sz==branchIdx+1) {
        molList.resize( sz-1 );
      }
      return molIdx;
    };
    //! \overload
    int addBranchToMol(std::vector<RWMol *> &molList,unsigned int molIdx,unsigned int branchIdx){
      Bond *newBond=new Bond(Bond::SINGLE);
      return addBranchToMol(molList,molIdx,branchIdx,newBond);
    };

    // ------------------------------------------------------------------------------------
    //! adds the atoms and bonds from a fragment to the molecule, sets no bond between them
    int addFragToMol(std::vector<RWMol *> &molList,unsigned int molIdx,unsigned int fragIdx){
      Bond *newBond=new Bond(Bond::IONIC);
      return addBranchToMol(molList,molIdx,fragIdx,newBond);
    }

    //! convenience function to convert the argument to a string
    template <typename T>
    std::string convertToString(T val){
      std::string res=boost::lexical_cast<std::string>(val);
      return  res;
    }

    void CleanupAfterParseError(RWMol *mol){
      PRECONDITION(mol,"no molecule");
      // blow out any partial bonds:
      RWMol::BOND_BOOKMARK_MAP *marks = mol->getBondBookmarks();
      RWMol::BOND_BOOKMARK_MAP::iterator markI=marks->begin();
      while(markI != marks->end()){
        RWMol::BOND_PTR_LIST &bonds=markI->second;
        for(RWMol::BOND_PTR_LIST::iterator bondIt=bonds.begin();
            bondIt!=bonds.end();++bondIt){
          delete *bondIt;
        }
        ++markI;
      }
    }
  } // end of namespace SLNParse
} // end of namespace RDKit
#endif
