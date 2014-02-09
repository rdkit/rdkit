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
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include "SmilesParse.h"
#include "SmilesParseOps.h"
#include <list>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <boost/foreach.hpp>
#include <RDGeneral/RDLog.h>

namespace SmilesParseOps{
  using namespace RDKit;

  void ReportParseError(const char *message,bool throwIt){
    if(!throwIt) BOOST_LOG(rdErrorLog) <<"SMILES Parse Error: "<< message << std::endl;
    else throw SmilesParseException(message);
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

  //
  // set bondOrder to Bond::IONIC to skip the formation of a bond
  //  between the fragment and the molecule
  //
  void AddFragToMol(RWMol *mol,RWMol *frag,Bond::BondType bondOrder,
      Bond::BondDir bondDir,bool closeRings,bool doingQuery ){
    PRECONDITION(mol,"no molecule");
    PRECONDITION(frag,"no fragment");
    PRECONDITION(mol->getActiveAtom(),"no active atom");
    Atom *lastAt = mol->getActiveAtom();
    int nOrigAtoms = mol->getNumAtoms();
    int nOrigBonds = mol->getNumBonds();

    //
    // close any rings we can in the fragment
    //
    if(closeRings){
      CloseMolRings(frag,true);
    }

    //
    // Add the fragment's atoms and bonds to the molecule:
    //
    mol->insertMol(*frag);

    //
    // update ring-closure order information on the added atoms:
    //
    for(RWMol::AtomIterator atomIt=frag->beginAtoms();
        atomIt!=frag->endAtoms();atomIt++){
      if((*atomIt)->hasProp("_RingClosures")){
        INT_VECT tmpVect;
        (*atomIt)->getProp("_RingClosures",tmpVect);
        BOOST_FOREACH(int &v, tmpVect){
          // if the ring closure is not already a bond, don't touch it:
          if(v>=0) v += nOrigBonds;
        }
        Atom *newAtom = mol->getAtomWithIdx(nOrigAtoms+(*atomIt)->getIdx());
        newAtom->setProp("_RingClosures",tmpVect);
      }
    }


    //
    //  ses up the bond between the mol and the branch
    //
    if( bondOrder != Bond::IONIC ){
      // FIX: this is not so much with the elegance...
      Atom *firstAt = mol->getAtomWithIdx(nOrigAtoms);
      Bond::BondType bo;
      int atomIdx1,atomIdx2;
      atomIdx1 = firstAt->getIdx();
      atomIdx2 = lastAt->getIdx();
      if(frag->hasBondBookmark(ci_LEADING_BOND)){
        //std::cout << "found it" << std::endl;
        const ROMol::BOND_PTR_LIST &leadingBonds=frag->getAllBondsWithBookmark(ci_LEADING_BOND);
        BOOST_FOREACH(Bond *leadingBond,leadingBonds){
          // we've already got a bond, so just set its local info
          // and then add it to the molecule intact (no sense doing
          // any extra work).
          leadingBond->setOwningMol(mol);
          leadingBond->setEndAtomIdx(leadingBond->getBeginAtomIdx()+nOrigAtoms);
          leadingBond->setBeginAtomIdx(atomIdx2);
          mol->addBond(leadingBond,true);
        }
        mol->clearBondBookmark(ci_LEADING_BOND);
      } else {
        if(!doingQuery){
          if(bondOrder == Bond::UNSPECIFIED){
            // no bond order provided, figure it out ourselves
            if(lastAt->getIsAromatic() && firstAt->getIsAromatic() ){
              bo = Bond::AROMATIC;
            }
            else{
              bo = Bond::SINGLE;
            }
          }
          else{
            bo = bondOrder;
          }
          if(bo==Bond::DATIVEL){
            int tmp=atomIdx2;
            atomIdx2 = atomIdx1;
            atomIdx1 = tmp;
            bo = Bond::DATIVE;
          } else if(bo == Bond::DATIVER){
            bo = Bond::DATIVE;
          }
          int idx = mol->addBond(atomIdx2,atomIdx1,bo)-1;
          mol->getBondWithIdx(idx)->setBondDir(bondDir);
        } else {
          // semantics are different in SMARTS, unspecified bonds can be single or aromatic:
          if(bondOrder == Bond::UNSPECIFIED){
            QueryBond *newB=new QueryBond(Bond::SINGLE);
            newB->expandQuery(makeBondOrderEqualsQuery(Bond::AROMATIC),
                Queries::COMPOSITE_OR,
                true);
            newB->setOwningMol(mol);
            newB->setBeginAtomIdx(atomIdx1);
            newB->setEndAtomIdx(atomIdx2);
            mol->addBond(newB);
            delete newB;
          }
          else{
            bo=bondOrder;
            if(bo==Bond::DATIVEL){
              int tmp=atomIdx2;
              atomIdx2 = atomIdx1;
              atomIdx1 = tmp;
              bo = Bond::DATIVE;
            } else if(bo == Bond::DATIVER){
              bo = Bond::DATIVE;
            }
            int idx = mol->addBond(atomIdx2,atomIdx1,bo)-1;
            mol->getBondWithIdx(idx)->setBondDir(bondDir);
          }
        }
      }
    }

    //
    // okay, the next thing we have to worry about is the possibility
    // that there might be ring opening/closing in the fragment we just
    // dealt with e.g. for things like C1C(C1) and C1C.C1
    // We deal with this by copying in the bookmarks and partial bonds
    // that exist in the fragment
    //
    RWMol::ATOM_BOOKMARK_MAP::iterator atIt;
    for(atIt=frag->getAtomBookmarks()->begin();
        atIt!=frag->getAtomBookmarks()->end();
        ++atIt){
      // don't bother even considering bookmarks outside
      // the range used for loops
      if(atIt->first < 100 && atIt->first > 0){
        RWMol::ATOM_PTR_LIST::iterator otherAt;
        for(otherAt=atIt->second.begin();
            otherAt!=atIt->second.end();otherAt++){
          Atom *at2 = *otherAt;
          int newIdx = at2->getIdx()+nOrigAtoms;
          mol->setAtomBookmark(mol->getAtomWithIdx(newIdx),atIt->first);
          //frag->clearAtomBookmark(atIt->first,at2);
          while( frag->hasBondBookmark(atIt->first) ){
            Bond *b = frag->getBondWithBookmark(atIt->first);
            int atomIdx1 = b->getBeginAtomIdx()+nOrigAtoms;
            b->setOwningMol(mol);
            b->setBeginAtomIdx(atomIdx1);
            mol->setBondBookmark(b,atIt->first);
            frag->clearBondBookmark(atIt->first,b);
          }
        }
      }
    }

    frag->clearAllAtomBookmarks();
    frag->clearAllBondBookmarks();
  };


  void _invChiralRingAtomWithHs(Atom *atom) {
    PRECONDITION(atom,"bad atom");
    // we will assume that this function is called on a ring atom with a 
    // ring closure bond 
    if (atom->getNumExplicitHs() == 1) {
      atom->invertChirality();
    }
  }
  typedef std::pair<int,int> INT_PAIR;
  bool operator<(const INT_PAIR &p1,const INT_PAIR &p2){
    return p1.first<p2.first;
  }

  void AdjustAtomChiralityFlags(RWMol *mol){
    PRECONDITION(mol,"no molecule");
    for(RWMol::AtomIterator atomIt=mol->beginAtoms();
        atomIt != mol->endAtoms();
        ++atomIt){
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
#if 0
        std::cout << "CLOSURES: ";
        std::copy(ringClosures.begin(),ringClosures.end(),
            std::ostream_iterator<int>(std::cout," "));
        std::cout << std::endl;
#endif	
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
              static_cast<int>(nbrBond->getIdx()))== ringClosures.end()){
            neighbors.push_back(std::make_pair(*nbrIdx,nbrBond->getIdx()));
          }
          ++nbrIdx;
        }
        // sort the list of non-ring-closure bonds:
        neighbors.sort();


        // find the location of this atom.  it pretty much has to be
        // first in the list, e.g for smiles like [C@](F)(Cl)(Br)I, or
        // second (everything else).
        std::list<INT_PAIR>::iterator selfPos=neighbors.begin();
        if(selfPos->first != static_cast<int>((*atomIt)->getIdx())){
          ++selfPos;
        }
        CHECK_INVARIANT(selfPos->first==static_cast<int>((*atomIt)->getIdx()),"weird atom ordering");

        // copy over the bond ids:
        INT_LIST bondOrdering;
        for(std::list<INT_PAIR>::iterator neighborIt=neighbors.begin();
            neighborIt != neighbors.end(); ++neighborIt){
          if(neighborIt != selfPos){
            bondOrdering.push_back(neighborIt->second);
          } else {
            // we are not going to add the atom itself, but we will push on
            // ring closure bonds at this point (if required):
            BOOST_FOREACH(int closure,ringClosures){
              bondOrdering.push_back(closure);
            }
          }
        }

        // ok, we now have the SMILES ordering of the bonds, figure out the
        // permutation order.
        //
        //  This whole thing is necessary because the ring-closure bonds
        //  in the SMILES come before the bonds to the other neighbors, but
        //  they come after the neighbors in the molecule we build.
        //  A crude example:
        //   in F[C@](Cl)(Br)I the C-Cl bond is index 1 in both SMILES
        //         and as built
        //   in F[C@]1(Br)I.Cl1 the C-Cl bond is index 1 in the SMILES
        //         and index 3 as built.
        //
        int nSwaps=(*atomIt)->getPerturbationOrder(bondOrdering);
        // FIX: explain this one:
        if((*atomIt)->getDegree()==3 && (*atomIt)->hasProp("_SmilesStart")) ++nSwaps;
        if(nSwaps%2){
          (*atomIt)->invertChirality();
        }
      }
    }
  }
  
  Bond::BondType GetUnspecifiedBondType(const RWMol *mol,const Atom *atom1,const Atom *atom2){
    PRECONDITION(mol,"no molecule");
    PRECONDITION(atom1,"no atom1");
    PRECONDITION(atom2,"no atom2");
    Bond::BondType res;
    if(atom1->getIsAromatic() && atom2->getIsAromatic()) {
      res = Bond::AROMATIC;
    } else {
      res = Bond::SINGLE;
    }
    return res;
  }
  void CloseMolRings(RWMol *mol,bool toleratePartials){
    //  Here's what we want to do here:
    //    loop through the molecule's atom bookmarks
    //    for each bookmark:
    //       connect pairs of atoms sharing that bookmark
    //          left to right (in the order in which they were
    //          inserted into the molecule).
    //       whilst doing this, we have to be cognizant of the fact that
    //          there may well be partial bonds in the molecule which need
    //          to be tied in as well.  WOO HOO! IT'S A BIG MESS!
    PRECONDITION(mol,"no molecule");
    RWMol::ATOM_BOOKMARK_MAP::iterator bookmarkIt;
    bookmarkIt=mol->getAtomBookmarks()->begin();
    while(bookmarkIt!=mol->getAtomBookmarks()->end()){
      // don't bother even considering bookmarks outside
      // the range used for loops
      if(bookmarkIt->first < 100 && bookmarkIt->first >= 0){
        RWMol::ATOM_PTR_LIST::iterator atomIt,atomsEnd;
        RWMol::ATOM_PTR_LIST bookmarkedAtomsToRemove;
        atomIt = bookmarkIt->second.begin();
        atomsEnd = bookmarkIt->second.end();
        while(atomIt != atomsEnd){
          Atom *atom1 = *atomIt;
          ++atomIt;
          if(!toleratePartials && atomIt==atomsEnd){
            ReportParseError("unclosed ring");
          } else if(atomIt!=atomsEnd && *atomIt==atom1){
            // make sure we don't try to connect an atom to itself,
            // this was bug 3145697:
            ++atomIt;
          } else if(atomIt!=atomsEnd) {
            // we actually found an atom, so connect it to the first
            Atom *atom2 = *atomIt;
            ++atomIt;
            
            int bondIdx=-1;
            // We're guaranteed two partial bonds, one for each time
            // the ring index was used.  We give the first specification
            // priority.
            CHECK_INVARIANT(mol->hasBondBookmark(bookmarkIt->first),"Missing bond bookmark");

            // now use the info from the partial bond:
            // The partial bond itself will have a proper order and directionality
            // (with a minor caveat documented below) and will have its beginning
            // atom set already:
            RWMol::BOND_PTR_LIST bonds=mol->getAllBondsWithBookmark(bookmarkIt->first);
            RWMol::BOND_PTR_LIST::iterator bondIt=bonds.begin();
            CHECK_INVARIANT(bonds.size()>=2,"Missing bond");

            // get pointers to the two bonds:
            Bond *bond1=*bondIt;
            ++bondIt;
            Bond *bond2=*bondIt;

            // remove those bonds from the bookmarks:
            mol->clearBondBookmark(bookmarkIt->first,bond1);
            mol->clearBondBookmark(bookmarkIt->first,bond2);

            // Make sure the bonds have the correct starting atoms:
            CHECK_INVARIANT(bond1->getBeginAtomIdx()==atom1->getIdx(),"bad begin atom");
            CHECK_INVARIANT(bond2->getBeginAtomIdx()==atom2->getIdx(),"bad begin atom");

            Bond *matchedBond;

            // figure out which (if either) bond has a specified type, we'll
            // keep that one.  We also need to update the end atom index to match
            // FIX: daylight barfs when you give it multiple specs for the closure
            //   bond, we'll just take the first one and ignore others
            //   NOTE: we used to do this the other way (take the last specification),
            //   but that turned out to be troublesome in odd cases like C1CC11CC1.
            if(!bond1->hasProp("_unspecifiedOrder")){
              matchedBond = bond1;
              matchedBond->setEndAtomIdx(atom2->getIdx());
              delete bond2;
            } else {
              matchedBond = bond2;
              matchedBond->setEndAtomIdx(atom1->getIdx());
              delete bond1;
            }
            if(matchedBond->getBondType()==Bond::UNSPECIFIED){
              Bond::BondType bondT=GetUnspecifiedBondType(mol,atom1,atom2);
              matchedBond->setBondType(bondT);
            }
            matchedBond->setOwningMol(mol);
            if(matchedBond->getBondType()==Bond::AROMATIC){
              matchedBond->setIsAromatic(true);
            }

#if 0
            //
            //  In cases like this: Cl\C=C1.F/1, we need to
            //  reverse the directionality on the added bond
            //  (because the bond is added from C -> F, but the
            //  directionality is for F->C.  We recognize these
            //  cases because the matched bond direction isn't
            //  the same as the added bond direction (i.e. atom1
            //  isn't the begin atom for the matched bond).
            //
            //  This was Issue 175
            if(atom1->getIdx()!=matchedBond->getBeginAtomIdx()){
              switch(matchedBond->getBondDir()){
              case Bond::ENDUPRIGHT:
                matchedBond->setBondDir(Bond::ENDDOWNRIGHT);
                break;
              case Bond::ENDDOWNRIGHT:
                matchedBond->setBondDir(Bond::ENDUPRIGHT);
                break;
              default:
                break;
              }
            }
#endif
            // add the bond:
            bondIdx=mol->addBond(matchedBond,true);

            // we found a bond, so update the atom's _RingClosures
            // property:
            if(bondIdx>-1){
              CHECK_INVARIANT(atom1->hasProp("_RingClosures") &&
                  atom2->hasProp("_RingClosures"),
                  "somehow atom doesn't have _RingClosures property.");
              INT_VECT closures;
              atom1->getProp("_RingClosures",closures);
              INT_VECT::iterator closurePos= std::find(closures.begin(),
                  closures.end(),
                  -(bookmarkIt->first+1));
              CHECK_INVARIANT(closurePos!=closures.end(),
                  "could not find bookmark in atom _RingClosures");
              *closurePos = bondIdx-1;
              atom1->setProp("_RingClosures",closures);

              atom2->getProp("_RingClosures",closures);
              closurePos= std::find(closures.begin(),
                  closures.end(),
                  -(bookmarkIt->first+1));
              CHECK_INVARIANT(closurePos!=closures.end(),
                  "could not find bookmark in atom _RingClosures");
              *closurePos = bondIdx-1;
              atom2->setProp("_RingClosures",closures);
            }
            bookmarkedAtomsToRemove.push_back(atom1);
            bookmarkedAtomsToRemove.push_back(atom2);
          }
        }
        //
        // increment the bookmark before calling erase.  Otherwise we
        // get a seg fault under MSVC++
        //
        int mark=bookmarkIt->first;
        bookmarkIt++;
        RWMol::ATOM_PTR_LIST::const_iterator aplci;
        BOOST_FOREACH(Atom *atom,bookmarkedAtomsToRemove){
          mol->clearAtomBookmark(mark,atom);
        }
      } else {
        ++bookmarkIt;
      }
    }
  };

} // end of namespace SmilesParseOps
