// $Id$
//
//  Copyright (C) 2001-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Canon.h>
#include <RDBoost/Exceptions.h>

namespace Canon {
  using namespace RDKit;
  PossibleType makePossible(int rank,int atomIdx,Bond *bond) {
    return std::make_pair(rank,std::make_pair(atomIdx,bond));
  };

  int _possibleComp(const PossibleType &arg1,const PossibleType &arg2) {
    return (arg1.first < arg2.first);
  };

  void switchBondDir(Bond *bond){
    PRECONDITION(bond,"bad bond");
    PRECONDITION(bond->getBondType()==Bond::SINGLE,"bad bond type");
    switch(bond->getBondDir()){
    case Bond::ENDUPRIGHT:
      bond->setBondDir(Bond::ENDDOWNRIGHT);
      break;
    case Bond::ENDDOWNRIGHT:
      bond->setBondDir(Bond::ENDUPRIGHT);
      break;
    default:
      break;
    }
  }
  
  // FIX: this may only be of interest from the SmilesWriter, should we
  // move it there?
  //
  //
  void canonicalizeDoubleBond(Bond *dblBond,
                              INT_VECT &bondVisitOrders,
                              INT_VECT &atomVisitOrders,
                              INT_VECT &bondDirCounts){
    PRECONDITION(dblBond,"bad bond");
    PRECONDITION(dblBond->getBondType() == Bond::DOUBLE,"bad bond order");
    PRECONDITION(dblBond->getStereo() > Bond::STEREOANY,"bad bond stereo");
    PRECONDITION(dblBond->getStereoAtoms().size() >= 2,"bad bond stereo atoms");
    PRECONDITION(atomVisitOrders[dblBond->getBeginAtomIdx()]>0 ||
                 atomVisitOrders[dblBond->getEndAtomIdx()]>0,
                 "neither end atom traversed");

    // atom1 is the lower numbered atom of the double bond (the one traversed
    // first)
    Atom *atom1,*atom2;
    if(atomVisitOrders[dblBond->getBeginAtomIdx()] <
       atomVisitOrders[dblBond->getEndAtomIdx()] ){
      atom1 = dblBond->getBeginAtom();
      atom2 = dblBond->getEndAtom();
    } else {
      atom1 = dblBond->getEndAtom();
      atom2 = dblBond->getBeginAtom();
    }

    // we only worry about double bonds that begin and end at atoms
    // of degree 2 or 3:
    if( (atom1->getDegree() != 2 && atom1->getDegree() != 3) ||
        (atom2->getDegree() != 2 && atom2->getDegree() != 3) ) {
      return;
    }
    
    Bond *firstFromAtom1=NULL,*secondFromAtom1=NULL;
    Bond *firstFromAtom2=NULL,*secondFromAtom2=NULL;

    int firstVisitOrder=100000;
    
    ROMol::OBOND_ITER_PAIR atomBonds;
    ROMol::GRAPH_MOL_BOND_PMAP::type pMap = dblBond->getOwningMol().getBondPMap();

    // -------------------------------------------------------
    // find the lowest visit order bonds from each end and determine
    // if anything is already constraining our choice of directions:
    bool dir1Set=false,dir2Set=false;
    atomBonds = dblBond->getOwningMol().getAtomBonds(atom1);
    while( atomBonds.first != atomBonds.second ){
      if(pMap[*atomBonds.first] != dblBond){
        int bondIdx = pMap[*atomBonds.first]->getIdx();      
        if( bondDirCounts[bondIdx] > 0 ){
          dir1Set = true;
        }
        if(!firstFromAtom1 || bondVisitOrders[bondIdx] < firstVisitOrder ){
          if(firstFromAtom1) secondFromAtom1 = firstFromAtom1;
          firstFromAtom1 = pMap[*atomBonds.first];
          firstVisitOrder = bondVisitOrders[bondIdx];
        } else {
          secondFromAtom1 = pMap[*atomBonds.first];
        }
      }
      atomBonds.first++;
    }
    atomBonds = dblBond->getOwningMol().getAtomBonds(atom2);
    firstVisitOrder = 10000;
    while( atomBonds.first != atomBonds.second ){
      if(pMap[*atomBonds.first] != dblBond){
        int bondIdx = pMap[*atomBonds.first]->getIdx();      
        if( bondDirCounts[bondIdx] > 0 ){
          dir2Set = true;
        }
        if(!firstFromAtom2 || bondVisitOrders[bondIdx] < firstVisitOrder){
          if(firstFromAtom2) secondFromAtom2 = firstFromAtom2;
          firstFromAtom2 = pMap[*atomBonds.first];
          firstVisitOrder = bondVisitOrders[bondIdx];
        } else {
          secondFromAtom2 = pMap[*atomBonds.first];
        }
      }
      atomBonds.first++;
    }

    if(dir2Set){
      // I am pretty sure the only way this can happen is if we hit a double
      // bond that closes a ring and that has stereochem info encoded.
      // At the moment we are explicitly not supporting this situation:
      CHECK_INVARIANT(0,"Stereochemistry specification on ring double bonds is not supported");

      // if this restriction is lifted, additional code will have to be added below.
    }
    
    // make sure we found everything we need to find:
    CHECK_INVARIANT(firstFromAtom1,"could not find atom1");
    CHECK_INVARIANT(firstFromAtom2,"could not find atom2");
    CHECK_INVARIANT(atom1->getDegree()==2 || secondFromAtom1,"inconsistency at atom1");
    CHECK_INVARIANT(atom2->getDegree()==2 || secondFromAtom2,"inconsistency at atom2");

    Bond::BondDir atom1Dir=Bond::NONE;
    Bond *atom1ControllingBond=firstFromAtom1;
    if(!dir1Set && !dir2Set){
      // ----------------------------------
      // nothing has touched our bonds so far, so set the
      // directions to "arbitrary" values:

      // the bond we came in on becomes ENDUPRIGHT:
      atom1Dir = Bond::ENDUPRIGHT;
      firstFromAtom1->setBondDir(atom1Dir);
      bondDirCounts[firstFromAtom1->getIdx()] += 1;

    } else if(!dir2Set){
      // at least one of the bonds on atom1 has its directionality set already:
      if(bondDirCounts[firstFromAtom1->getIdx()]>0){
        // The first bond's direction has been set at some earlier point:
        atom1Dir = firstFromAtom1->getBondDir();
        bondDirCounts[firstFromAtom1->getIdx()] += 1;
        if(secondFromAtom1){
          // both bonds have their directionalities set, make sure
          // they are compatible:
          CHECK_INVARIANT(firstFromAtom1->getBondDir() !=
                          secondFromAtom1->getBondDir(),"inconsistent state");
        }
      } else {
        // the second bond must be present and setting the direction:
        CHECK_INVARIANT(secondFromAtom1,"inconsistent state");
        CHECK_INVARIANT(bondDirCounts[secondFromAtom1->getIdx()]>0,"inconsistent state");
        // It must be the second bond setting the direction.
        // This happens when the bond dir is set in a branch:
        //        v- this double bond
        //   CC(/C=P/N)=N/O
        //      ^- the second bond sets the direction
        // or when the first bond is a ring closure from an
        // earlier traversed atom:
        //             v- this double bond
        //   NC1=NOC/C1=N\O
        //     ^- this closure ends up being the first bond,
        //        and it does not set the direction.
        //
        // This addresses parts of Issue 185 and sf.net Issue 1842174
        // 
        atom1Dir = secondFromAtom1->getBondDir();

        firstFromAtom1->setBondDir(atom1Dir);
        bondDirCounts[firstFromAtom1->getIdx()]+=1;
        bondDirCounts[secondFromAtom1->getIdx()]+=1;
        atom1ControllingBond=secondFromAtom1;
      }
    } else {
      CHECK_INVARIANT(0,"Stereochemistry specification on ring double bonds is not supported");
    }

    // now set the directionality on the other side:
    Bond::BondDir atom2Dir=Bond::NONE;
    if( dblBond->getStereo() == Bond::STEREOE ){
      atom2Dir = atom1Dir;
    } else if( dblBond->getStereo() == Bond::STEREOZ ){
      atom2Dir = (atom1Dir==Bond::ENDUPRIGHT) ? Bond::ENDDOWNRIGHT : Bond::ENDUPRIGHT;
    }
    CHECK_INVARIANT(atom2Dir != Bond::NONE,"stereo not set");

    // If we're not looking at the bonds used to determine the
    // stereochemistry, we need to flip the setting on the other bond:
    const INT_VECT &stereoAtoms=dblBond->getStereoAtoms();
    if(atom1->getDegree()==3 &&
       std::find(stereoAtoms.begin(),stereoAtoms.end(),
                 static_cast<int>(atom1ControllingBond->getOtherAtomIdx(atom1->getIdx())))
       == stereoAtoms.end() ){
      atom2Dir = (atom2Dir == Bond::ENDUPRIGHT) ? Bond::ENDDOWNRIGHT : Bond::ENDUPRIGHT;
    }
    if(atom2->getDegree()==3 &&
       std::find(stereoAtoms.begin(),stereoAtoms.end(),
                 static_cast<int>(firstFromAtom2->getOtherAtomIdx(atom2->getIdx()))) == stereoAtoms.end() ){
      atom2Dir = (atom2Dir == Bond::ENDUPRIGHT) ? Bond::ENDDOWNRIGHT : Bond::ENDUPRIGHT;
    }
        
    firstFromAtom2->setBondDir(atom2Dir);
    bondDirCounts[firstFromAtom2->getIdx()] += 1;
    
    // -----------------------------------
    //
    // Check if there are other bonds from atoms 1 and 2 that need
    // to have their directionalities set:
    ///
    if(atom1->getDegree() == 3 && !bondDirCounts[secondFromAtom1->getIdx()] ){
      // This bond (the second bond from the starting atom of the double bond)
      // is a special case.  It's going to appear in a branch in the smiles:
      //     X\C(\Y)=C/Z
      //         ^
      //         |- here
      // so it actually needs to go down with the *same* direction as the
      // bond that's already been set (because "pulling the bond out of the
      // branch reverses its direction).
      // A quick example.  This SMILES:
      //     F/C(\Cl)=C/F
      // is *wrong*. This is the correct form:
      //     F/C(/Cl)=C/F
      // So, since we want this bond to have the opposite direction to the
      // other one, we put it in with the same direction.
      // This was Issue 183
      secondFromAtom1->setBondDir(firstFromAtom1->getBondDir());
      bondDirCounts[secondFromAtom1->getIdx()] += 1;
    }

    if(atom2->getDegree() == 3 && !bondDirCounts[secondFromAtom2->getIdx()] ){
      // Here we set the bond direction to be opposite the other one (since
      // both come after the atom connected to the double bond).
      Bond::BondDir otherDir;
      otherDir = (firstFromAtom2->getBondDir()==Bond::ENDUPRIGHT) ? Bond::ENDDOWNRIGHT : Bond::ENDUPRIGHT;
      secondFromAtom2->setBondDir(otherDir);
      bondDirCounts[secondFromAtom2->getIdx()] += 1;
    }

    // This is an odd case... The bonds off the beginning atom are
    // after the start atom in the traversal stack.  These need to
    // have their directions reversed.  An example SMILES (unlikely
    // to actually traverse this way is:
    //   C(=C/O)/F    or C(/F)=C/O
    // That bond is Z, without the reversal, this would come out:
    //   C(=C/O)\F    or C(\F)=C/O
    // which is E.
    //
    // In the case of three-coordinate atoms, we don't need to flip
    // the second bond because the Issue 183 fix (above) already got
    // that one.
    //
    // This was Issue 191 and continued into sf.net issue 1842174
    if( bondVisitOrders[atom1ControllingBond->getIdx()] >
        atomVisitOrders[atom1->getIdx()]){
      if(bondDirCounts[atom1ControllingBond->getIdx()]==1){
        switchBondDir(atom1ControllingBond);
      } else if(bondDirCounts[firstFromAtom2->getIdx()]==1){
        // the controlling bond at atom1 is being set by someone else, flip the direction
        // on the atom2 bond instead:
        switchBondDir(firstFromAtom2);
      }
    }
  }
  
  void canonicalDFSTraversal(ROMol &mol,int atomIdx,int inBondIdx,
                             std::vector<AtomColors> &colors,
                             VECT_INT_VECT &cycles,
                             INT_VECT &ranks,
                             INT_VECT &cyclesAvailable,
                             MolStack &molStack,
                             INT_VECT &atomOrders,
                             INT_VECT &bondVisitOrders){
    //ROMol *mol = molProps.getMol();
    int nAttached=0;

    Atom *atom = mol.getAtomWithIdx(atomIdx);
    // the atom will keep track of the order in which it sees bonds
    // using its _TraversalBondIndexOrder list:
    INT_LIST directTravList,cycleEndList;
    INT_VECT ringClosures(0);
    atom->setProp("_CanonRingClosureBondIndices",ringClosures,true);

    molStack.push_back(MolStackElem(atom));
    atomOrders[atom->getIdx()] = molStack.size();
    //atom->setProp("_CanonTravOrder",molStack.size(),1);
    colors[atomIdx] = GREY_NODE;

    // ---------------------
    //
    //  Build the list of possible destinations from here
    //
    // ---------------------
    std::vector< PossibleType > possibles;
    possibles.resize(0);
    ROMol::OBOND_ITER_PAIR bondsPair = mol.getAtomBonds(atom);
    possibles.reserve(bondsPair.second-bondsPair.first);
    ROMol::GRAPH_MOL_BOND_PMAP::type pMap = mol.getBondPMap();

    while(bondsPair.first != bondsPair.second){
      Bond *theBond = pMap[*(bondsPair.first)];
      if(inBondIdx<0 || theBond->getIdx() != static_cast<unsigned int>(inBondIdx)){
        int otherIdx = theBond->getOtherAtomIdx(atomIdx);
        long rank=ranks[otherIdx];
        // ---------------------
        //
        // things are a bit more complicated if we are sitting on a
        // ring atom we would like to traverse first to the
        // ring-closure atoms, then to atoms outside the ring first,
        // then to atoms in the ring that haven't already been visited
        // (non-ring-closure atoms).
        // 
        //  Here's how the black magic works:
        //   - non-ring atom neighbors have their original ranks
        //   - ring atom neighbors have this added to their ranks:
        //       (Bond::OTHER - bondOrder)*MAX_NATOMS*MAX_NATOMS
        //   - ring-closure neighbors lose a factor of:
        //       (Bond::OTHER+1)*MAX_NATOMS*MAX_NATOMS
        //
        //  This tactic biases us to traverse to non-ring neighbors first,
        //  original ordering if bond orders are all equal... crafty, neh?
        //  
        // ---------------------
        if( colors[otherIdx] == GREY_NODE ) {
          rank -= static_cast<int>(Bond::OTHER+1) *
            MAX_NATOMS*MAX_NATOMS;
          rank += static_cast<int>(Bond::OTHER - theBond->getBondType()) *
            MAX_NATOMS;
        } else if( theBond->getOwningMol().getRingInfo()->numBondRings(theBond->getIdx()) ){
          rank += static_cast<int>(Bond::OTHER - theBond->getBondType()) *
            MAX_NATOMS*MAX_NATOMS;
        }
        possibles.push_back(makePossible(rank,otherIdx,theBond));
      }
      bondsPair.first++;
    }

    // ---------------------
    //
    //  Sort on ranks
    //
    // ---------------------
    std::sort(possibles.begin(),possibles.end(),_possibleComp);


    // ---------------------
    //
    //  Now work the children
    //
    // ---------------------
    std::vector<MolStack> subStacks;
    for(std::vector<PossibleType>::iterator possiblesIt=possibles.begin();
        possiblesIt!=possibles.end();
        possiblesIt++){
      MolStack subStack;
      int possibleIdx = possiblesIt->second.first;
      Bond *bond = possiblesIt->second.second;
      Atom *otherAtom=mol.getAtomWithIdx(possibleIdx);
      INT_LIST otherTravList;
      unsigned int lowestRingIdx;
      switch(colors[possibleIdx]){
      case WHITE_NODE:
        // -----
        // we haven't seen this node at all before
        // -----

        // it might have some residual data from earlier calls, clean that up:
        if(otherAtom->hasProp("_TraversalBondIndexOrder")){
          otherAtom->clearProp("_TraversalBondIndexOrder");
        }

        directTravList.push_back(bond->getIdx());
        subStack.push_back(MolStackElem(bond,atomIdx));
        canonicalDFSTraversal(mol,possibleIdx,bond->getIdx(),colors,
                              cycles,ranks,cyclesAvailable,subStack,
                              atomOrders,bondVisitOrders);
        subStacks.push_back(subStack);
        nAttached += 1;
        break;
      case GREY_NODE:
        // -----
        // we've seen this, but haven't finished it (we're finishing a ring)
        // -----
        cycleEndList.push_back(bond->getIdx());
        lowestRingIdx = std::find(cyclesAvailable.begin(),
                                  cyclesAvailable.end(),1) -
          cyclesAvailable.begin();
        cyclesAvailable[lowestRingIdx] = 0;
        cycles[possibleIdx].push_back(lowestRingIdx);
        ++lowestRingIdx;

        molStack.push_back(MolStackElem(lowestRingIdx));

        // we need to add this bond (which closes the ring) to the traversal list for the
        // other atom as well:
        if(otherAtom->hasProp("_TraversalBondIndexOrder")){
          otherAtom->getProp("_TraversalBondIndexOrder",otherTravList);
        } else {
          otherTravList.clear();
        }
        otherTravList.push_back(bond->getIdx());
        otherAtom->setProp("_TraversalBondIndexOrder",otherTravList,true);

        otherAtom->getProp("_CanonRingClosureBondIndices",ringClosures);
        ringClosures.push_back(bond->getIdx());
        otherAtom->setProp("_CanonRingClosureBondIndices",ringClosures,true);

        break;
      default:
        // -----
        // this node has been finished. don't do anything.
        // -----
        break;
      }
    }
    

    atom->getProp("_CanonRingClosureBondIndices",ringClosures);
    
    CHECK_INVARIANT(ringClosures.size()==cycles[atomIdx].size(),
                    "ring closure mismatch");
    for(unsigned int i=0;i<ringClosures.size();i++){
      int ringIdx=cycles[atomIdx][i];
      ringIdx += 1;
      molStack.push_back(MolStackElem(mol.getBondWithIdx(ringClosures[i]),
                                      atom->getIdx()));
      molStack.push_back(MolStackElem(ringIdx));
    }
    cycles[atomIdx].resize(0);
  
    MolStack::const_iterator ciMS;
    for(int i=0;i<nAttached;i++){
      if(i<nAttached-1){
        int branchIdx=0;
        if(subStacks[i].begin()->type==MOL_STACK_ATOM){
          branchIdx=subStacks[i].begin()->obj.atom->getIdx();
        } else if(subStacks[i].begin()->type==MOL_STACK_BOND){
          branchIdx=-1*subStacks[i].begin()->obj.bond->getIdx();
        } else {
          ASSERT_INVARIANT(0,"branch started with something other than an atom or bond");
        }
        molStack.push_back(MolStackElem("(",branchIdx));
        for(ciMS=subStacks[i].begin();ciMS!=subStacks[i].end();ciMS++){
          molStack.push_back(*ciMS);
          switch(ciMS->type){
          case MOL_STACK_ATOM:
            atomOrders[ciMS->obj.atom->getIdx()] = molStack.size();
            break;
          case MOL_STACK_BOND:
            bondVisitOrders[ciMS->obj.bond->getIdx()] = molStack.size();
            break;
          default:
            break;
          }
        }
        molStack.push_back(MolStackElem(")",branchIdx));
      } else {
        for(ciMS=subStacks[i].begin();ciMS!=subStacks[i].end();ciMS++){
          molStack.push_back(*ciMS);
          switch(ciMS->type){
          case MOL_STACK_ATOM:
            atomOrders[ciMS->obj.atom->getIdx()] = molStack.size();
            break;
          case MOL_STACK_BOND:
            bondVisitOrders[ciMS->obj.bond->getIdx()] = molStack.size();
            //ciMS->obj.bond->setProp("_CanonTravOrder",molStack.size(),1);
            break;
          default:
            break;
          }
        }
      }
    }

    //std::cerr<<"*****>>>>>> Traversal results for atom: "<<atom->getIdx()<<"> ";
    INT_LIST travList;
    // first push on the incoming bond:
    if(inBondIdx >= 0){
      //std::cerr<<" "<<inBondIdx;
      travList.push_back(inBondIdx);
    }

    // ... ring closures that end here:
    for(INT_LIST_CI ilci=cycleEndList.begin();ilci!=cycleEndList.end();++ilci){
      //std::cerr<<" ["<<*ilci<<"]";
      travList.push_back(*ilci);
    }


    // ... ring closures that start here:
    if(atom->hasProp("_TraversalBondIndexOrder")){
      INT_LIST indirectTravList;
      atom->getProp("_TraversalBondIndexOrder",indirectTravList);
      for(INT_LIST_CI ilci=indirectTravList.begin();ilci!=indirectTravList.end();++ilci){
        //std::cerr<<" ("<<*ilci<<")";
        travList.push_back(*ilci);
      }
    }
    // and finally the bonds we directly traverse:
    for(INT_LIST_CI ilci=directTravList.begin();ilci!=directTravList.end();++ilci){
      //std::cerr<<" "<<*ilci;
      travList.push_back(*ilci);
    }
    //std::cerr<<"\n";
    atom->setProp("_TraversalBondIndexOrder",travList);
    colors[atomIdx] = BLACK_NODE;
  }

  bool canHaveDirection(const Bond *bond){
    PRECONDITION(bond,"bad bond");
    Bond::BondType bondType= bond->getBondType();
    return (bondType==Bond::SINGLE || bondType==Bond::AROMATIC);
  }

  void clearBondDirs(ROMol &mol,Bond *refBond,Atom *fromAtom,
                     INT_VECT &bondDirCounts){
    PRECONDITION(bondDirCounts.size()>=mol.getNumBonds(),"bad dirCount size");
    PRECONDITION(refBond,"bad bond");
    PRECONDITION(&refBond->getOwningMol()==&mol,"bad bond");
    PRECONDITION(fromAtom,"bad atom");
    PRECONDITION(&fromAtom->getOwningMol()==&mol,"bad bond");
    
    ROMol::GRAPH_MOL_BOND_PMAP::type pMap = mol.getBondPMap();
    ROMol::OEDGE_ITER beg,end;
    boost::tie(beg,end) = mol.getAtomBonds(fromAtom);
    bool nbrPossible=false,cleared=false;
    while(beg!=end){
      if( pMap[*beg] != refBond && canHaveDirection(pMap[*beg])){
        nbrPossible=true;
        if(bondDirCounts[pMap[*beg]->getIdx()] <= bondDirCounts[refBond->getIdx()]){
          bondDirCounts[pMap[*beg]->getIdx()] = 0;
          // no one is setting the direction here:
          pMap[*beg]->setBondDir(Bond::NONE);
          cleared=true;
        }
      }
      beg++;
    }
    if(nbrPossible && !cleared){
      // we found a neighbor that could have directionality set,
      // but it had a higher bondDirCount that us, so we must
      // need to be cleared:
      refBond->setBondDir(Bond::NONE);
    }
  }

  void removeRedundantBondDirSpecs(ROMol &mol,MolStack &molStack,INT_VECT bondDirCounts){
    PRECONDITION(bondDirCounts.size()>=mol.getNumBonds(),"bad dirCount size");
    // find bonds that have directions indicated that are redundant:
    for(MolStack::iterator msI=molStack.begin();
        msI!=molStack.end(); msI++) {
      if( msI->type == MOL_STACK_BOND ){
        Bond *tBond = msI->obj.bond;
        if(canHaveDirection(tBond) &&
           bondDirCounts[tBond->getIdx()]>=1 ) {
          ROMol::GRAPH_MOL_BOND_PMAP::type pMap = mol.getBondPMap();
          ROMol::OEDGE_ITER beg,end;

          // start by finding the double bond that sets tBond's direction:
          Atom *dblBondAtom=NULL;
          boost::tie(beg,end) = mol.getAtomBonds(tBond->getBeginAtom());
          while(beg!=end){
            if( pMap[*beg] != tBond && pMap[*beg]->getBondType()==Bond::DOUBLE &&
                pMap[*beg]->getStereo() > Bond::STEREOANY ){
              dblBondAtom = tBond->getBeginAtom();
              break;
            }
            beg++;
          }
          if(dblBondAtom != NULL){
            clearBondDirs(mol,tBond,dblBondAtom,bondDirCounts);
          }
          dblBondAtom = NULL;
          boost::tie(beg,end) = mol.getAtomBonds(tBond->getEndAtom());
          while(beg!=end){
            if( pMap[*beg] != tBond && pMap[*beg]->getBondType()==Bond::DOUBLE &&
                pMap[*beg]->getStereo() > Bond::STEREOANY ){
              dblBondAtom = tBond->getEndAtom();
              break;
            }
            beg++;
          }
          if(dblBondAtom != NULL){
            clearBondDirs(mol,tBond,dblBondAtom,bondDirCounts);
          }
        }
      }
    }
  }

  void canonicalizeFragment(ROMol &mol,int atomIdx,
                            std::vector<AtomColors> &colors,
                            INT_VECT &ranks,
                            MolStack &molStack){
    int nAtoms=mol.getNumAtoms();
    
    INT_VECT atomVisitOrders(nAtoms,0);
    INT_VECT bondVisitOrders(mol.getNumBonds(),0);
    INT_VECT bondDirCounts(mol.getNumBonds(),0);

    INT_VECT cyclesAvailable;
    INT_VECT_I viIt;
    cyclesAvailable.resize(MAX_CYCLES);
    for(viIt=cyclesAvailable.begin();viIt!=cyclesAvailable.end();viIt++) *viIt=1;

    VECT_INT_VECT  cycles;
    cycles.resize(nAtoms);
    VECT_INT_VECT_I vviIt;
    for(vviIt=cycles.begin();vviIt!=cycles.end();vviIt++) vviIt->resize(0);

    // make sure that we've done the bond stereo perception:
    MolOps::assignBondStereoCodes(mol,false);

    // we need ring information make sure findSSSR has been called before
    // if not call now
    if ( !mol.getRingInfo()->isInitialized() ) {
      MolOps::findSSSR(mol);
    }
    mol.getAtomWithIdx(atomIdx)->setProp("_TraversalStartPoint",true);
    if(mol.getAtomWithIdx(atomIdx)->hasProp("_TraversalBondIndexOrder")){
      mol.getAtomWithIdx(atomIdx)->clearProp("_TraversalBondIndexOrder");
    }
    Canon::canonicalDFSTraversal(mol,atomIdx,-1,colors,cycles,
                                 ranks,cyclesAvailable,molStack,atomVisitOrders,
                                 bondVisitOrders);
    
    // remove the current directions on single bonds around double bonds:
    for(ROMol::BondIterator bondIt=mol.beginBonds();
        bondIt!=mol.endBonds();
        bondIt++){
      Bond::BondDir dir = (*bondIt)->getBondDir();
      if (dir == Bond::ENDDOWNRIGHT || dir == Bond::ENDUPRIGHT) {
        (*bondIt)->setBondDir(Bond::NONE);
      }
    }

    // traverse the stack and clean up double bonds
    for(MolStack::iterator msI=molStack.begin();
        msI!=molStack.end(); msI++){
      if(msI->type == MOL_STACK_BOND &&
         msI->obj.bond->getBondType() == Bond::DOUBLE &&
         msI->obj.bond->getStereo() > Bond::STEREOANY){
        Canon::canonicalizeDoubleBond(msI->obj.bond,bondVisitOrders,atomVisitOrders,
                                      bondDirCounts);
      }
    }

    Canon::removeRedundantBondDirSpecs(mol,molStack,bondDirCounts);
  }
};




