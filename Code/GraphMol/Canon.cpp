// $Id$
//
//  Copyright (C) 2001-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Canon.h>

namespace Canon {
  using namespace RDKit;
  PossibleType makePossible(int rank,int atomIdx,Bond *bond) {
    return std::make_pair(rank,std::make_pair(atomIdx,bond));
  };

  int _possibleComp(const PossibleType arg1,const PossibleType arg2) {
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
	if(!firstFromAtom1 || bondVisitOrders[bondIdx] < firstVisitOrder){
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
      // bond that closes a ring and that has sterechem info encoded.
      // At the moment we are explicitly not supporting this situation:
      CHECK_INVARIANT(0,"Stereochemistry specification on ring double bonds is not supported");

      // if this restriction is lifted, additional code will have to be added below.
    }
    
    // make sure we found everything we need to find:
    CHECK_INVARIANT(firstFromAtom1,"could not find first atom");
    CHECK_INVARIANT(firstFromAtom2,"could not find first atom");
    CHECK_INVARIANT(atom1->getDegree()==2 || secondFromAtom1,"could not find second atom");
    CHECK_INVARIANT(atom2->getDegree()==2 || secondFromAtom2,"could not find second atom");

    Bond::BondDir atom1Dir=Bond::NONE;
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
	// This happens when the bond dir is set in a branch, e.g. in
	// molecules like: 

	// This addresses part of Issue 185
	// 
	//atom1Dir = (secondFromAtom1->getBondDir()==Bond::ENDUPRIGHT) ?
	//  Bond::ENDDOWNRIGHT : Bond::ENDUPRIGHT;
	atom1Dir = secondFromAtom1->getBondDir();

	firstFromAtom1->setBondDir(atom1Dir);
	bondDirCounts[firstFromAtom1->getIdx()]+=1;
	bondDirCounts[secondFromAtom1->getIdx()]+=1;
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
		 firstFromAtom1->getOtherAtomIdx(atom1->getIdx())) == stereoAtoms.end() ){
      atom2Dir = (atom2Dir == Bond::ENDUPRIGHT) ? Bond::ENDDOWNRIGHT : Bond::ENDUPRIGHT;
    }
    if(atom2->getDegree()==3 &&
       std::find(stereoAtoms.begin(),stereoAtoms.end(),
		 firstFromAtom2->getOtherAtomIdx(atom2->getIdx())) == stereoAtoms.end() ){
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
      Bond::BondDir otherDir;
      otherDir = firstFromAtom1->getBondDir();
      secondFromAtom1->setBondDir(otherDir);
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
    // after the double bond in the traversal stack.  These need to
    // have their directions reversed.  An example SMILES (unlikely
    // to actually traverse this way is:
    //   C(=C/O)/F
    // That bond is Z, without the reversal, this would come out:
    //   C(=C/O)\F
    // which is E.
    //
    // In the case of three-coordinate atoms, we don't need to flip
    // the second bond because the Issue 183 fix (above) already got
    // that one.
    //
    // This was Issue 191
    if( bondVisitOrders[firstFromAtom1->getIdx()] >
	bondVisitOrders[dblBond->getIdx()] &&
	bondDirCounts[firstFromAtom1->getIdx()]==1){
      switchBondDir(firstFromAtom1);
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
    INT_LIST travList,directTravList;
    if(inBondIdx >= 0){
      travList.push_back(inBondIdx);
    }
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

    //std::cerr<<" TRAVERSE: (" << atomIdx<<") ";
    while(bondsPair.first != bondsPair.second){
      Bond *theBond = pMap[*(bondsPair.first)];
      if(theBond->getIdx() != inBondIdx){
	int otherIdx = theBond->getOtherAtomIdx(atomIdx);
	long rank=ranks[otherIdx];
	// ---------------------
	//
	// things are a bit more complicated if we are sitting on a ring atom
	// we would like to traverse to atoms outside the ring first, then
	// to non-ring-closure atoms (atoms that haven't already been visited),
	// then, finally, to ring-closure atoms.  This business with ring-closure
	// visitation is to ensure that the canonical smiles we generate is
	// actually CORRECT.  (Someday I'll come back and explain this remark).
	// 
	//
	//  Here's how the black magic works:
	//   - non-ring atom neighbors have their original ranks
	//   - ring atom neighbors have this added to their ranks:
	//       (Bond::OTHER - bondOrder)*MAX_NATOMS*MAX_NATOMS
	//   - ring-closure neighbors have an additional factor of:
	//       (Bond::OTHER+1)*MAX_NATOMS*MAX_NATOMS
	//     added.
	//
	//  This tactic biases us to traverse to non-ring neighbors first
	//  original ordering if bond orders are all equal... crafty, neh?
	//  
	// ---------------------
	if( colors[otherIdx] == GREY_NODE ) {
	  rank += static_cast<int>(Bond::OTHER+1) *
	    MAX_NATOMS*MAX_NATOMS;
	}
	// FIX: this will not work
	if( theBond->getOwningMol().getRingInfo()->numBondRings(theBond->getIdx()) ){
	  rank += static_cast<int>(Bond::OTHER - theBond->getBondType()) *
	    MAX_NATOMS*MAX_NATOMS;
	}
	possibles.push_back(makePossible(rank,otherIdx,theBond));
	//std::cerr << otherIdx << "_" << ranks[otherIdx]<<"_"<<rank<<" ";
      }
      bondsPair.first++;
    }
    //std::cerr<<std::endl;

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
      int lowestRingIdx;
      Bond *bond = possiblesIt->second.second;
      Atom *otherAtom;
      INT_LIST otherTravList;

      switch(colors[possibleIdx]){
      case WHITE_NODE:
	// mmm, fresh baked node
	subStack.push_back(MolStackElem(bond,atomIdx));
	// we may need to update this atom's traversal order whilst
	// dealing with the rest of the molecule, so set travList
	// now and grab it back after the recursion:
	directTravList.push_back(bond->getIdx());
	atom->setProp("_TraversalBondIndexOrder",travList);
	canonicalDFSTraversal(mol,possibleIdx,bond->getIdx(),colors,
			      cycles,ranks,cyclesAvailable,subStack,
			      atomOrders,bondVisitOrders);
	atom->getProp("_TraversalBondIndexOrder",travList);
	subStacks.push_back(subStack);
	nAttached += 1;
	break;
      case GREY_NODE:
	// kind of stale, it closes a ring
	lowestRingIdx = std::find(cyclesAvailable.begin(),cyclesAvailable.end(),1) -
	  cyclesAvailable.begin();
	cyclesAvailable[lowestRingIdx] = 0;
	cycles[possibleIdx].push_back(lowestRingIdx);
	lowestRingIdx += 1;
	// we're not going to push the bond on here, but save it until later
	//molStack.push_back(MolStackElem(bond,atomIdx));
	//bondVisitOrders[bond->getIdx()] = molStack.size();
	//bond->setProp("_CanonTravOrder",molStack.size(),1);
	
	molStack.push_back(MolStackElem(lowestRingIdx));
	travList.push_back(bond->getIdx());

	// we need to add this bond to the traversal list for the
	// other atom as well:
	otherAtom=mol.getAtomWithIdx(possibleIdx);
	otherAtom->getProp("_TraversalBondIndexOrder",otherTravList);
	otherTravList.push_back(bond->getIdx());
	otherAtom->setProp("_TraversalBondIndexOrder",otherTravList);

	otherAtom->getProp("_CanonRingClosureBondIndices",ringClosures);
	ringClosures.push_back(bond->getIdx());
	otherAtom->setProp("_CanonRingClosureBondIndices",ringClosures,true);
	
	break;
      default:
	// whoa! this is hard as a rock! We will skip it!
	break;
      }
    }
    
    if(atom->hasProp("_CanonRingClosureBondIndices")){
      atom->getProp("_CanonRingClosureBondIndices",ringClosures);
    } else {
      ringClosures.resize(0);
    }
    
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
	molStack.push_back(MolStackElem("("));
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
      
	molStack.push_back(MolStackElem(")"));
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

    for(INT_LIST_CI ilci=directTravList.begin();ilci!=directTravList.end();ilci++){
      travList.push_back(*ilci);
    }
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
    while(beg!=end){
      if( pMap[*beg] != refBond && canHaveDirection(pMap[*beg]) ){
	if(bondDirCounts[pMap[*beg]->getIdx()]==1){
	  bondDirCounts[pMap[*beg]->getIdx()] = 0;
	  // no one is setting the direction here:
	  pMap[*beg]->setBondDir(Bond::NONE);
	}
      }
      beg++;
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




