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
#include "Bond.h"
#include "Atom.h"
#include "ROMol.h"
#include <RDGeneral/Invariant.h>

namespace RDKit {

Bond::Bond() {
  initBond();
  dp_props = new Dict();
};

Bond::Bond(BondType bT) {
  initBond();
  d_bondType = bT;
  dp_props = new Dict();
};

Bond::Bond(const Bond &other){
  // NOTE: we do *not* copy ownership!
  dp_mol = 0;
  d_bondType = other.d_bondType;
  d_beginAtomIdx = other.d_beginAtomIdx;
  d_endAtomIdx = other.d_endAtomIdx;
  d_dirTag = other.d_dirTag;
  d_stereo = other.d_stereo;
  d_stereoAtoms = other.d_stereoAtoms;
  df_isAromatic = other.df_isAromatic;
  df_isConjugated = other.df_isConjugated;
  d_index = other.d_index;
  if(other.dp_props){
    dp_props = new Dict(*other.dp_props);
  } else {
    dp_props = new Dict();
  }
}

Bond::~Bond()
{
  delete dp_props;
}

Bond &Bond::operator=(const Bond &other){
  dp_mol = other.dp_mol;
  d_bondType = other.d_bondType;
  d_beginAtomIdx = other.d_beginAtomIdx;
  d_endAtomIdx = other.d_endAtomIdx;
  d_dirTag = other.d_dirTag;
  d_stereoAtoms = other.d_stereoAtoms;
  df_isAromatic = other.df_isAromatic;
  df_isConjugated = other.df_isConjugated;
  d_index = other.d_index;
  if(other.dp_props){
    dp_props = new Dict(*other.dp_props);
  }else{
    dp_props = new Dict();
  }

    
  return *this;
}

Bond *Bond::copy() const {
  Bond *res = new Bond(*this);
  return res;
}


void Bond::setOwningMol(ROMol *other)
{
  // FIX: doesn't update topology
  dp_mol = other;
}

unsigned int Bond::getOtherAtomIdx(const unsigned int thisIdx) const
{
  PRECONDITION(d_beginAtomIdx == thisIdx ||
	       d_endAtomIdx == thisIdx, "bad index");
  if( d_beginAtomIdx == thisIdx ) return d_endAtomIdx;
  else if (d_endAtomIdx == thisIdx) return d_beginAtomIdx;
  // we cannot actually get down here
  return 0;
}

void Bond::setBeginAtomIdx(unsigned int what) {
  if(dp_mol) RANGE_CHECK(0,what,getOwningMol().getNumAtoms()-1);
  d_beginAtomIdx = what;
};

void Bond::setEndAtomIdx(unsigned int what) {
  if(dp_mol) RANGE_CHECK(0,what,getOwningMol().getNumAtoms()-1);
  d_endAtomIdx = what;
};


void Bond::setBeginAtom(Atom *at) {
  PRECONDITION( dp_mol != 0, "no owning molecule for bond");
  setBeginAtomIdx(at->getIdx());
}
void Bond::setBeginAtom(Atom::ATOM_SPTR at) {
  PRECONDITION( dp_mol != 0, "no owning molecule for bond");
  setBeginAtomIdx(at->getIdx());
}

void Bond::setEndAtom(Atom *at) {
  PRECONDITION( dp_mol != 0, "no owning molecule for bond");
  setEndAtomIdx(at->getIdx());
}
void Bond::setEndAtom(Atom::ATOM_SPTR at) {
  PRECONDITION( dp_mol != 0, "no owning molecule for bond");
  setEndAtomIdx(at->getIdx());
}

  Atom *Bond::getBeginAtom() const {
    PRECONDITION( dp_mol != 0, "no owning molecule for bond");
    return dp_mol->getAtomWithIdx(d_beginAtomIdx);
  };
  Atom *Bond::getEndAtom() const {
    PRECONDITION( dp_mol != 0, "no owning molecule for bond");
    return dp_mol->getAtomWithIdx(d_endAtomIdx);
  };
  Atom *Bond::getOtherAtom(Atom const *what) const {
    PRECONDITION( dp_mol != 0, "no owning molecule for bond");

    return dp_mol->getAtomWithIdx(getOtherAtomIdx(what->getIdx()));
  };


double Bond::getBondTypeAsDouble() const {
  switch(getBondType()){
  case UNSPECIFIED: return 0; break;
  case IONIC: return 0; break;
  case SINGLE: return 1; break;
  case DOUBLE: return 2; break;
  case TRIPLE: return 3; break;
  case QUADRUPLE: return 4; break;
  case QUINTUPLE: return 5; break;
  case HEXTUPLE: return 6; break;
  case ONEANDAHALF: return 1.5; break;
  case TWOANDAHALF: return 2.5; break;
  case THREEANDAHALF: return 3.5; break;
  case FOURANDAHALF: return 4.5; break;
  case FIVEANDAHALF: return 5.5; break;
  case AROMATIC: return 1.5; break;
  case DATIVEONE: return 1.0; break; // FIX: this should probably be different
  case DATIVE: return 1.0; break; //FIX: again probably wrong
  case ZERO: return 0; break; 
  default:
    UNDER_CONSTRUCTION("Bad bond type");
  }
}

double Bond::getValenceContrib(Atom::ATOM_SPTR at) const {
  return getValenceContrib(at.get());
}

double Bond::getValenceContrib(const Atom *atom) const {
  switch(getBondType()){
  case UNSPECIFIED: return 0; break;
  case IONIC: return 0; break;
  case SINGLE: return 1; break;
  case DOUBLE: return 2; break;
  case TRIPLE: return 3; break;
  case QUADRUPLE: return 4; break;
  case QUINTUPLE: return 5; break;
  case HEXTUPLE: return 6; break;
  case ONEANDAHALF: return 1.5; break;
  case TWOANDAHALF: return 2.5; break;
  case THREEANDAHALF: return 3.5; break;
  case FOURANDAHALF: return 4.5; break;
  case FIVEANDAHALF: return 5.5; break;
  case AROMATIC: return 1.5; break;
  case DATIVEONE:
    if(atom->getIdx()==getEndAtomIdx())return 1.0;
    else return 0.0;
    break;
  case DATIVE:
    if(atom->getIdx()==getEndAtomIdx())return 1.0;
    else return 0.0;
    break;
  case ZERO: return 0; break; 
  default:
    UNDER_CONSTRUCTION("Bad bond type");

  }
}

void Bond::setQuery(QUERYBOND_QUERY *what) {
  //  Bonds don't have queries at the moment because I have not
  //  yet figured out what a good base query should be.
  //  It would be nice to be able to do substructure searches
  //  using molecules alone, so it'd be nice if we got this
  //  issue resolved ASAP. 
  PRECONDITION(0,"plain bonds have no Query");
}

Bond::QUERYBOND_QUERY *Bond::getQuery() const {
  PRECONDITION(0,"plain bonds have no Query");
  return NULL;
};

bool Bond::Match(Bond const *what) const{
  bool res;
  if(getBondType()==Bond::UNSPECIFIED ||
     what->getBondType()==Bond::UNSPECIFIED){    
    res = true;
  } else {
    res = getBondType() == what->getBondType();
  }
  return res;
};

bool Bond::Match(const Bond::BOND_SPTR what) const {
  return Match(what.get());
};
  
void Bond::expandQuery(Bond::QUERYBOND_QUERY *what,
		       Queries::CompositeQueryType how,
		       bool maintainOrder) {
  PRECONDITION(0,"plain bonds have no query");
};

void Bond::initBond(){
  d_bondType = UNSPECIFIED;
  d_dirTag = NONE;
  d_stereo = STEREONONE;
  d_stereoAtoms.clear();
  dp_mol = 0;
  d_beginAtomIdx = 0;
  d_endAtomIdx = 0;
  df_isAromatic = 0;
  d_index = 0;
  df_isConjugated = 0;
};

}; // end o' namespace


std::ostream & operator<<(std::ostream& target, const RDKit::Bond &bond){
  target << bond.getIdx() << " ";
  target << bond.getBeginAtomIdx() << "->" << bond.getEndAtomIdx();
  target << " order: " << bond.getBondType();
  if(bond.getBondDir())
    target << " dir: " << bond.getBondDir();
  if(bond.getStereo())
    target << " stereo: " << bond.getStereo();
  target << " conj?: " << bond.getIsConjugated();
  target << " aromatic?: " << bond.getIsAromatic();

  return target;
}
