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
#include <math.h>

#include "ROMol.h"
#include "Atom.h"
#include "PeriodicTable.h"
#include "SanitException.h"
#include "MDLValence.h"
#include "QueryOps.h"
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/types.h>

namespace RDKit {
  namespace {
    // Determine whether or not a molecule is to the left of Carbon
    bool isEarlyAtom(int atomicNum){
      return (4 - PeriodicTable::getTable()->getNouterElecs(atomicNum)) > 0;
    }
  }
Atom::Atom(){
  d_atomicNum=0;
  initAtom();
}

Atom::Atom(unsigned int num) {
  d_atomicNum = num;
  initAtom();
};

Atom::Atom(std::string what) {
  d_atomicNum = PeriodicTable::getTable()->getAtomicNumber(what);
  initAtom();
};

Atom::Atom( const Atom & other){
  // NOTE: we do *not* copy ownership!
  d_atomicNum = other.d_atomicNum;
  dp_mol = 0;
  d_index=0;
  d_formalCharge = other.d_formalCharge;
  df_noImplicit = other.df_noImplicit;
  df_isAromatic = other.df_isAromatic;
  d_dativeFlag = other.d_dativeFlag;
  d_numRadicalElectrons = other.d_numRadicalElectrons;
  d_mass = other.d_mass;
  d_isotope = other.d_isotope;
  //d_pos = other.d_pos;
  d_chiralTag=other.d_chiralTag;
  d_hybrid = other.d_hybrid;
  d_implicitValence=other.d_implicitValence;
  d_explicitValence=other.d_explicitValence;
  if(other.dp_props){
    dp_props = new Dict(*other.dp_props);
  } else {
    dp_props = new Dict();
    STR_VECT computed;
    dp_props->setVal("__computedProps", computed);
  }
}
void Atom::initAtom(){
  df_isAromatic = false;
  df_noImplicit = false;
  d_dativeFlag=0;
  d_numRadicalElectrons=0;
  d_formalCharge = 0;
  d_index = 0;
  if(d_atomicNum){
    d_mass = PeriodicTable::getTable()->getAtomicWeight(d_atomicNum);
  } else{
    d_mass = 0.0;
  }
  d_isotope=0;
  d_chiralTag=CHI_UNSPECIFIED;
  d_hybrid = UNSPECIFIED;
  dp_mol = 0;
  dp_props = new Dict();

  d_implicitValence=-1;
  d_explicitValence=-1;

  // ok every Atom contains a property entry called "__computedProps"
  // which provides list of property keys that correspond to value
  // that have been computed this can used to blow out all computed
  // properties while leaving the rest along initialize this list to
  // an empty vector of strings
  //STR_VECT computed;
  //dp_props->setVal("__computedProps", computed);
}

Atom::~Atom()
{
  if(dp_props){
    delete dp_props;
    dp_props = 0;
  }
}

Atom *Atom::copy() const {
  Atom *res = new Atom(*this);
  return res;
}


void Atom::setOwningMol(ROMol *other)
{
  // NOTE: this operation does not update the topology of the owning
  // molecule (i.e. this atom is not added to the graph).  Only
  // molecules can add atoms to themselves.
  dp_mol = other;
}


std::string Atom::getSymbol() const {
  std::string res;
  // handle dummies differently:
  if(d_atomicNum != 0 || !hasProp("dummyLabel") ){
    res = PeriodicTable::getTable()->getElementSymbol(d_atomicNum);
  } else {
    getProp("dummyLabel",res);
  }
  return res;
}

unsigned int Atom::getDegree() const {
  PRECONDITION(dp_mol,"degree not defined for atoms not associated with molecules");
  return getOwningMol().getAtomDegree(this);
}

unsigned int Atom::getTotalDegree() const {
  PRECONDITION(dp_mol,"degree not defined for atoms not associated with molecules");
  unsigned int res=this->getNumImplicitHs()+this->getDegree();
  return res;
}

//
//  If includeNeighbors is set, we'll loop over our neighbors
//   and include any of them that are Hs in the count here
//
unsigned int Atom::getTotalNumHs() const {
  PRECONDITION(dp_mol,"valence not defined for atoms not associated with molecules")
  int res = getNumImplicitHs();
  ROMol::ADJ_ITER begin,end;
  const ROMol *parent = &getOwningMol();
  boost::tie(begin,end) = parent->getAtomNeighbors(this);
  while(begin!=end){
    const Atom *at = parent->getAtomWithIdx(*begin);
    if(at->getAtomicNum()==1) res++;
    ++begin;
  }
  return res;
}

unsigned int Atom::getNumImplicitHs() const {
  PRECONDITION(d_implicitValence>-1,
               "getNumImplicitHs() called without preceding call to calcImplicitValence()");
  return getImplicitValence();
}

int Atom::getExplicitValence() const {
  PRECONDITION(dp_mol,"valence not defined for atoms not associated with molecules");
  PRECONDITION(d_explicitValence>-1,
               "getExplicitValence() called without call to calcExplicitValence()");
  return d_explicitValence;
}

int Atom::calcExplicitValence(bool strict) {
  PRECONDITION(dp_mol,"valence not defined for atoms not associated with molecules");
  unsigned int res=0;
  unsigned int arom=0;
  ROMol::OEDGE_ITER beg,end;
  boost::tie(beg,end) = getOwningMol().getAtomBonds(this);
  while(beg!=end){
    switch (getOwningMol()[*beg]->getBondType()) {
    case Bond::UNSPECIFIED:
    case Bond::IONIC:
      break;
    case Bond::DOUBLE:
      res += 2;
      break;
    case Bond::TRIPLE:
      res += 3;
      break;
    case Bond::QUADRUPLE:
      res += 4;
      break;
    case Bond::QUINTUPLE:
      res += 5;
      break;
    case Bond::HEXTUPLE:
      res += 6;
      break;
    case Bond::AROMATIC:
      arom = 1;
      /* Fall thru */
    case Bond::SINGLE:
    default:
      res++;
    }
    ++beg;
  }
  res += arom;

  // special case for when we call this on an aromatic atom
  // which could kekulize to having no double bonds to it (we're not
  // doing the actual kekulization step here, just checking that
  // it could happen):
  if(getIsAromatic()&&d_implicitValence>-1){
    int explicitPlusRadV = res + d_numRadicalElectrons;
    explicitPlusRadV+=d_implicitValence;
    unsigned int mdlvalence1 = MDLValence(d_atomicNum,d_formalCharge,explicitPlusRadV);
    unsigned int mdlvalence2 = MDLValence(d_atomicNum,d_formalCharge,explicitPlusRadV-1);
    if((mdlvalence1==explicitPlusRadV+1) ||
       (mdlvalence1==explicitPlusRadV && mdlvalence2==explicitPlusRadV-1)
       ){
      --res;
    }
  }
  
  d_explicitValence = res;
  return res;
}

int Atom::getImplicitValence() const {
  return d_implicitValence;
}

// NOTE: this uses the explicitValence, so it will call
// calcExplicitValence() if it hasn't already been called
int Atom::calcImplicitValence(bool strict) {
  PRECONDITION(dp_mol,"valence not defined for atoms not associated with molecules");
  if(df_noImplicit) {
    if(d_implicitValence<0) d_implicitValence=0;
  } else {
    if(d_explicitValence==-1) this->calcExplicitValence();

    int chg = getFormalCharge();
    int explicitPlusRadV = d_explicitValence + getNumRadicalElectrons();
    unsigned int mdlvalence = MDLValence(d_atomicNum,chg,explicitPlusRadV);
    d_implicitValence = mdlvalence - explicitPlusRadV;
  }
  return d_implicitValence;
}

void Atom::setIsotope(unsigned int what){
  d_isotope=what;
  if(d_isotope){
    d_mass=PeriodicTable::getTable()->getMassForIsotope(d_atomicNum,d_isotope);
  } else {
    d_mass = PeriodicTable::getTable()->getAtomicWeight(d_atomicNum);
  }
}

void Atom::setQuery(Atom::QUERYATOM_QUERY *what) {
  //  Atoms don't have complex queries so this has to fail
  PRECONDITION(0,"plain atoms have no Query");
}
Atom::QUERYATOM_QUERY *Atom::getQuery() const {
  return NULL;
};
void Atom::expandQuery(Atom::QUERYATOM_QUERY *what,
                       Queries::CompositeQueryType how,
                       bool maintainOrder) {
  PRECONDITION(0,"plain atoms have no Query");
}

bool Atom::Match(Atom const *what) const {
  PRECONDITION(what,"bad query atom");
  bool res = getAtomicNum() == what->getAtomicNum();
  // special dummy--dummy match case:
  //   [*] matches [*],[1*],[2*],etc.
  //   [1*] only matches [*] and [1*]
  if(res){
    if(!getAtomicNum()){
      // this is the new behavior, based on the isotopes:
      int tgt=this->getIsotope();
      int test=what->getIsotope();
      if(tgt && test && tgt!=test){
        res = false;
      }
    } else {
      // standard atom-atom match: The general rule here is that if this atom has a property that
      // deviates from the default, then the other atom should match that value.
      if( (this->getFormalCharge() && this->getFormalCharge()!=what->getFormalCharge()) ||
          (this->getIsotope() && this->getIsotope()!=what->getIsotope()) ){
        res=false;
      }
    }
  }
  return res;
}
bool Atom::Match(const Atom::ATOM_SPTR what) const {
  return Match(what.get());
}

void Atom::updatePropertyCache(bool strict) {
  calcExplicitValence();
  if(!hasQuery()){
    calcImplicitValence();
  } else {
    d_implicitValence=0;
  }
}

// returns the number of swaps required to convert the ordering
// of the probe list to match the order of our incoming bonds:
//
//  e.g. if our incoming bond order is: [0,1,2,3]:
//   getPerturbationOrder([1,0,2,3]) = 1
//   getPerturbationOrder([1,2,3,0]) = 3
//   getPerturbationOrder([1,2,0,3]) = 2
int Atom::getPerturbationOrder(INT_LIST probe) const{
  PRECONDITION(dp_mol,"perturbation order not defined for atoms not associated with molecules")
  INT_LIST ref;
  ROMol::OEDGE_ITER beg,end;
  boost::tie(beg,end) = getOwningMol().getAtomBonds(this);
  while(beg!=end){
    ref.push_back(getOwningMol()[*beg]->getIdx());
    ++beg;
  }
  int nSwaps=static_cast<int>(countSwapsToInterconvert(ref,probe));
  return nSwaps;
}

void Atom::invertChirality(){
  switch(getChiralTag()){
  case CHI_TETRAHEDRAL_CW:
    setChiralTag(CHI_TETRAHEDRAL_CCW);
    break;
  case CHI_TETRAHEDRAL_CCW:
    setChiralTag(CHI_TETRAHEDRAL_CW);
    break;
  case CHI_OTHER:
  case CHI_UNSPECIFIED:
    break;
  }
}


} // end o' namespace

std::ostream & operator<<(std::ostream& target, const RDKit::Atom &at){
  target << at.getIdx() << " " << at.getAtomicNum() << " " << at.getSymbol();
  target << " chg: " << at.getFormalCharge();
  target << "  deg: " << at.getDegree();
  target << " exp: ";
  try {
    int explicitValence = at.getExplicitValence();
    target << explicitValence;
  } catch (...){
    target << "N/A";
  }
  target << " imp: ";
  try {
    int implicitValence = at.getImplicitValence();
    target << implicitValence;
  } catch (...){
    target << "N/A";
  }
  target << " no_imp: "<<at.getNoImplicit();
  target << " hyb: " << at.getHybridization();
  target << " arom?: " << at.getIsAromatic();
  target << " chi: " << at.getChiralTag();
  if(at.getNumRadicalElectrons()){
    target << " rad: " << at.getNumRadicalElectrons();
  }
  if(at.getIsotope()){
    target << " iso: " << at.getIsotope();
  }
  return target;
};
