// $Id$
//
//  Copyright (C) 2002-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include "SmilesWrite.h"
#include <GraphMol/Canon.h>

#include <sstream>
#include <map>


namespace SmilesWrite{
  using namespace RDKit;

  const int atomicSmiles[] = {5,6,7,8,15,16,9,17,35,53,-1};

  std::string GetAtomSmiles(const Atom *atom,bool doKekule){
    int i;
    static bool firstCall=true;
    static INT_VECT atomicSmilesVect;
    if(firstCall){
      i = 0;
      while(atomicSmiles[i] != -1){
	atomicSmilesVect.push_back(atomicSmiles[i]);
	i++;
      }
      firstCall = false;
    }
    std::stringstream res;
    int fc = atom->getFormalCharge();
    int num = atom->getAtomicNum();
    int numExplicit = atom->getNumExplicitHs();
    std::string symb = atom->getSymbol();
    bool needsBracket;
    if(std::find(atomicSmilesVect.begin(),atomicSmilesVect.end(),num) !=
       atomicSmilesVect.end()){
      // it's a "normal" atom:
      if(!doKekule && atom->getIsAromatic() && symb[0] < 'a') symb[0] -= ('A'-'a');
      needsBracket = (fc != 0) || (numExplicit != 0);
      if(atom->getOwningMol().hasProp("_doIsoSmiles") &&
	 atom->getChiralTag()!=Atom::CHI_UNSPECIFIED ){
	needsBracket = true;
      }
    } else {
      needsBracket = true;
    }
    if( needsBracket ) res << "[";

    res << symb;

    if(atom->getOwningMol().hasProp("_doIsoSmiles") &&
       atom->getChiralTag()!=Atom::CHI_UNSPECIFIED ){
      INT_LIST trueOrder;
      atom->getProp("_TraversalBondIndexOrder",trueOrder);
#ifdef VERBOSE_CANON
      std::copy(trueOrder.begin(),trueOrder.end(),
		std::ostream_iterator<int>(std::cout," "));
      std::cout << std::endl;
#endif    
      int nSwaps =  atom->getPerturbationOrder(trueOrder);

    
      //std::cout << "\t\tnSwaps: " << nSwaps << std::endl;
      std::string atStr;
      switch(atom->getChiralTag()){
      case Atom::CHI_TETRAHEDRAL_CW:
	//std::cout << "\tcw" << std::endl;
	if(!(nSwaps%2))
	  atStr = "@@";
	else
	  atStr = "@";
	break;
      case Atom::CHI_TETRAHEDRAL_CCW:
	//std::cout << "\tccw" << std::endl;
	if(!(nSwaps%2))
	  atStr = "@";
	else
	  atStr = "@@";
	break;
      }
      //std::cout << "\tats: " << atStr << std::endl;
      res << atStr;
    }

    if(needsBracket){
      if(numExplicit > 0) res << "H";
      if(numExplicit > 1) res << numExplicit;
      if(fc > 0){
	res << "+";
	if(fc > 1) res << fc;
      } else if(fc < 0) {
	res << "-";
	if(fc < -1) res << -fc;
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

  std::string GetBondSmiles(const Bond *bond,int atomToLeftIdx,bool doKekule){
    std::stringstream res;
    bool aromatic=false;
    if( !doKekule &&
	(bond->getBondType() == Bond::SINGLE ||
	 bond->getBondType() == Bond::DOUBLE ||
	 bond->getBondType() == Bond::AROMATIC) ){
      Atom *a1,*a2;
      a1 = bond->getOwningMol().getAtomWithIdx(atomToLeftIdx);
      a2 = bond->getOwningMol().getAtomWithIdx(bond->getOtherAtomIdx(atomToLeftIdx));
      if(a1->getIsAromatic() && a2->getIsAromatic()) aromatic=true;
    }

    Bond::BondDir dir= bond->getBondDir();
    //if(bond->getBeginAtomIdx() != atomToLeftIdx) {
    //  if(dir==Bond::ENDDOWNRIGHT) dir = Bond::ENDUPRIGHT;
    //  else if(dir==Bond::ENDUPRIGHT) dir = Bond::ENDDOWNRIGHT;
    //}
    Bond::BondType btype = bond->getBondType();
    int oaid = bond->getOtherAtomIdx(atomToLeftIdx);
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
	if( aromatic && !bond->getIsAromatic() ) res << "-";
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
      }
      if(!aromatic) res << ":"; 
      break;
    case Bond::DATIVE:
      if(bond->getBeginAtomIdx() == atomToLeftIdx) res << ">";
      else res << "<";
      break;
    default:
      res << "?";
    }
    return res.str();
  }
  std::string FragmentSmilesConstruct(ROMol &mol,int atomIdx,
				      std::vector<Canon::AtomColors> &colors,
				      INT_VECT &ranks,bool doKekule){

    Canon::MolStack molStack;
    // try to prevent excessive reallocation
    molStack.reserve(mol.getNumAtoms()+
		     mol.getNumBonds());
    std::stringstream res;

    std::map<int,int> ringClosureMap;
    int ringIdx,closureVal;
    Canon::canonicalizeFragment(mol,atomIdx,colors,ranks,
				 molStack);
    Canon::MolStack::const_iterator msCI,tmpIt;
    Bond *bond;
    for(msCI=molStack.begin();msCI!=molStack.end();msCI++){
      switch(msCI->type){
      case Canon::MOL_STACK_ATOM:
	res << GetAtomSmiles(msCI->obj.atom,doKekule);
	break;
      case Canon::MOL_STACK_BOND:
	bond = msCI->obj.bond;

	// FIX:  This is kind of a hacky fix to a problem.
	//   What's going on here is that the traversal code currently
	//   (June 2004) puts ring-closure bonds with the final atom in
	//   the ring (choices are first and last atom).  The bond
	//   directions calculated in the canonicalization code are
	//   configured to be correct if the bond dir goes from the
	//   first to the last atom.  So if we have a bond closing a
	//   ring that has direction indicated, that direction will end
	//   up being wrong.  The "correct" solution would be to go and 
	//   change the traversal code to insert the bond into the stack
	//   after the initial atom, but that will break loads of
	//   existing canonical smiles and, potentially, open another
	//   can of bugs... so what we'll do instead here is look ahead
	//   in the stack and reverse the direction on ring-closure
	//   bonds at the time of writing.  It's a hack, but a simple
	//   and functional hack.
	//   This fixes Issue 184.
#if 0
	tmpIt = msCI+1;
	if(tmpIt!=molStack.end()&&tmpIt->type==MOL_STACK_RING){
	  tmpDir = bond->getBondDir();
	  switch(tmpDir){
	  case Bond::ENDUPRIGHT:
	    bond->setBondDir(Bond::ENDDOWNRIGHT);
	    break;
	  case Bond::ENDDOWNRIGHT:
	    bond->setBondDir(Bond::ENDUPRIGHT);
	    break;
	  }
	  res << GetBondSmiles(bond,msCI->number);
	  bond->setBondDir(tmpDir);
	} else {
	  res << GetBondSmiles(bond,msCI->number);
	}	
#else
	res << GetBondSmiles(bond,msCI->number,doKekule);
#endif
	break;
      case Canon::MOL_STACK_RING:
	ringIdx = msCI->number;
	if(ringClosureMap.count(ringIdx)){
	  // the index is already in the map ->
	  //   we're closing a ring, so grab
	  //   the index and then delete the value:
	  closureVal = ringClosureMap[ringIdx];
	  ringClosureMap.erase(ringIdx);
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

namespace RDKit{
 
  // NOTE: I did not forget the const here... Producing SMILES for
  // a molecule actually can change the molecule.  Specifically,
  // things like the directionality of bonds may be changed when
  // the molecule is canonicalized.
  // Odds are good that this may be one of those foot-shooting
  // decisions and I'm gonna want to smack myself for doing this,
  // but we'll try anyway.
  std::string MolToSmiles(ROMol &mol,bool doIsomericSmiles,
			  bool doKekule){
    if(doIsomericSmiles){
      mol.setProp("_doIsoSmiles",1);
    } else if(mol.hasProp("_doIsoSmiles")){
      mol.clearProp("_doIsoSmiles");
    }

#if 0
    std::cout << "----------------------------" << std::endl;
    std::cout << "MolToSmiles:"<< std::endl;
    mol.debugMol(std::cout);
    std::cout << "----------------------------" << std::endl;
#endif  
    std::string res;

    int nAtoms=mol.getNumAtoms();
    for(ROMol::AtomIterator atIt=mol.beginAtoms();atIt!=mol.endAtoms();atIt++)
      (*atIt)->updatePropertyCache();

    INT_VECT ranks;
    ranks.resize(nAtoms);

    // clean up the chirality on any atom that is marked as chiral,
    // but that should not be:
    if(doIsomericSmiles){
      MolOps::assignAtomChiralCodes(mol,true);
      //MolOps::assignBondStereoCodes(mol,true);
      MolOps::rankAtoms(mol,ranks,ComprehensiveInvariants,
			true,0,true);
    } else {
      MolOps::rankAtoms(mol,ranks);
    }
#ifdef VERBOSE_CANON
    for(int tmpI=0;tmpI<ranks.size();tmpI++){
      std::cout << tmpI << " " << ranks[tmpI] << " " << *(mol.getAtomWithIdx(tmpI)) << std::endl;
    }
#endif
  
    std::vector<Canon::AtomColors> colors(nAtoms,Canon::WHITE_NODE);
    std::vector<Canon::AtomColors>::iterator colorIt;
    colorIt = colors.begin();
    // loop to deal with the possibility that there might be disconnected fragments
    while(colorIt != colors.end()){
      int i,nextAtomIdx,nextRank;
      std::string subSmi;
      // find the next atom for a traverse
      nextRank = nAtoms+1;
      for(i=0;i<nAtoms;i++){
	if( colors[i] == Canon::WHITE_NODE && ranks[i] < nextRank ){
	  nextRank = ranks[i];
	  nextAtomIdx = i;
	}
      }

      subSmi = SmilesWrite::FragmentSmilesConstruct(mol, nextAtomIdx, colors,
						    ranks,doKekule);

      res += subSmi;
      colorIt = std::find(colors.begin(),colors.end(),Canon::WHITE_NODE);
      if(colorIt != colors.end()){
	res += ".";
      }
    }

    return res;
  }
}
