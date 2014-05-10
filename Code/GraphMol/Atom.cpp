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
#include "QueryOps.h"
#include "MonomerInfo.h"

#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/types.h>
#include <RDGeneral/Dict.h>

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
  d_numExplicitHs = other.d_numExplicitHs;
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
    dp_props->setVal(detail::computedPropName, computed);
  }
  if(other.dp_monomerInfo){
    dp_monomerInfo = other.dp_monomerInfo->copy();
  } else {
    dp_monomerInfo=0;
  }

}
void Atom::initAtom(){
  df_isAromatic = false;
  df_noImplicit = false;
  d_dativeFlag=0;
  d_numExplicitHs = 0;
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
  dp_monomerInfo = 0;

  d_implicitValence=-1;
  d_explicitValence=-1;

}

Atom::~Atom()
{
  if(dp_props){
    delete dp_props;
    dp_props = 0;
  }
  if(dp_monomerInfo){
    delete dp_monomerInfo;
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
    res=getProp<std::string>("dummyLabel");
  }
  return res;
}

unsigned int Atom::getDegree() const {
  PRECONDITION(dp_mol,"degree not defined for atoms not associated with molecules");
  return getOwningMol().getAtomDegree(this);
}

unsigned int Atom::getTotalDegree() const {
  PRECONDITION(dp_mol,"degree not defined for atoms not associated with molecules");
  unsigned int res=this->getTotalNumHs(false)+this->getDegree();
  return res;
}

//
//  If includeNeighbors is set, we'll loop over our neighbors
//   and include any of them that are Hs in the count here
//
unsigned int Atom::getTotalNumHs(bool includeNeighbors) const {
  PRECONDITION(dp_mol,"valence not defined for atoms not associated with molecules")
  int res = getNumExplicitHs() + getNumImplicitHs();
  if(includeNeighbors){
    ROMol::ADJ_ITER begin,end;
    const ROMol *parent = &getOwningMol();
    boost::tie(begin,end) = parent->getAtomNeighbors(this);
    while(begin!=end){
      const Atom *at = parent->getAtomWithIdx(*begin);
      if(at->getAtomicNum()==1) res++;
      ++begin;
    }
  }
  return res;
}

unsigned int Atom::getNumImplicitHs() const {
  if(df_noImplicit) return 0;
  
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

unsigned int Atom::getTotalValence() const {
  PRECONDITION(dp_mol,"valence not defined for atoms not associated with molecules");
  return getExplicitValence()+getImplicitValence();
}
  
int Atom::calcExplicitValence(bool strict) {
  PRECONDITION(dp_mol,"valence not defined for atoms not associated with molecules");
  unsigned int res;
  // FIX: contributions of bonds to valence are being done at best
  // approximately
  double accum=0;
  ROMol::OEDGE_ITER beg,end;
  boost::tie(beg,end) = getOwningMol().getAtomBonds(this);
  while(beg!=end){
    accum += getOwningMol()[*beg]->getValenceContrib(this);
    ++beg;
  }
  accum += getNumExplicitHs();

  // check accum is greater than the default valence
  unsigned int dv = PeriodicTable::getTable()->getDefaultValence(d_atomicNum);
  int chr = getFormalCharge();
  if(isEarlyAtom(d_atomicNum)) chr*=-1;  // <- the usual correction for early atoms
  if (accum > (dv + chr) && this->getIsAromatic()){
    // this needs some explanation : if the atom is aromatic and
    // accum > (dv + chr) we assume that no hydrogen can be added
    // to this atom.  We set x = (v + chr) such that x is the
    // closest possible integer to "accum" but less than
    // "accum".
    //
    // "v" here is one of the allowed valences. For example:
    //    sulfur here : O=c1ccs(=O)cc1
    //    nitrogen here : c1cccn1C
    
    int pval = dv + chr;
    const INT_VECT &valens = PeriodicTable::getTable()->getValenceList(d_atomicNum);
    for (INT_VECT_CI vi = valens.begin(); vi != valens.end() && *vi!=-1; ++vi) {
      int val = (*vi) + chr;
      if (val > accum) {
        break;
      } else {
        pval = val;
      }
    }
    accum = pval;
  }
  // despite promising to not to blame it on him - this a trick Greg
  // came up with: if we have a bond order sum of x.5 (i.e. 1.5, 2.5
  // etc) we would like it to round to the higher integer value -- 
  // 2.5 to 3 instead of 2 -- so we will add 0.1 to accum.
  // this plays a role in the number of hydrogen that are implicitly
  // added. This will only happen when the accum is a non-integer
  // value and less than the default valence (otherwise the above if
  // statement should have caught it). An example of where this can
  // happen is the following smiles:
  //    C1ccccC1
  // Daylight accepts this smiles and we should be able to Kekulize
  // correctly.
  accum += 0.1;

  res = static_cast<int>(round(accum));

  if(strict){
    int effectiveValence;
    if(PeriodicTable::getTable()->getNouterElecs(d_atomicNum)>=4){
      effectiveValence=res-getFormalCharge();
    } else {
      // for boron and co, we move to the right in the PT, so adding
      // extra valences means adding negative charge
      effectiveValence=res+getFormalCharge();
    }
    const INT_VECT &valens = PeriodicTable::getTable()->getValenceList(d_atomicNum);
    int maxValence=*(valens.rbegin());
    // maxValence == -1 signifies that we'll take anything at the high end
    if( maxValence>0 &&effectiveValence>maxValence){
      // the explicit valence is greater than any
      // allowed valence for the atoms - raise an error
      std::ostringstream errout;
      errout << "Explicit valence for atom # " << getIdx() 
             << " " << PeriodicTable::getTable()->getElementSymbol(d_atomicNum)
             << ", " << effectiveValence <<", is greater than permitted";
      std::string msg = errout.str();
      BOOST_LOG(rdErrorLog) << msg << std::endl;
      throw MolSanitizeException(msg);
    }
  }
  d_explicitValence = res;

  return res;
}

int Atom::getImplicitValence() const {
  PRECONDITION(dp_mol,"valence not defined for atoms not associated with molecules");
  if(df_noImplicit) return 0;
  return d_implicitValence;
}

// NOTE: this uses the explicitValence, so it will call
// calcExplictValence() if it hasn't already been called
int Atom::calcImplicitValence(bool strict) {
  PRECONDITION(dp_mol,"valence not defined for atoms not associated with molecules");
  if(df_noImplicit) return 0;
  if(d_explicitValence==-1) this->calcExplicitValence(strict);
  // this is basically the difference between the allowed valence of
  // the atom and the explicit valence already specified - tells how
  // many Hs to add
  // 
  int res;

  // The d-block and f-block of the periodic table (i.e. transition metals,
  // lanthanoids and actinoids) have no default valence.
  int dv = PeriodicTable::getTable()->getDefaultValence(d_atomicNum);
  if (dv==-1) {
    d_implicitValence = 0;
    return 0;
  }

  // here is how we are going to deal with the possibility of
  // multiple valences
  // - check the explicit valence "ev"
  // - if it is already equal to one of the allowed valences for the
  //    atom return 0
  // - otherwise take return difference between next larger allowed
  //   valence and "ev"
  // if "ev" is greater than all allowed valences for the atom raise an exception
  // finally aromatic cases are dealt with differently - these atoms are allowed
  // only default valences
  const INT_VECT &valens = PeriodicTable::getTable()->getValenceList(d_atomicNum);
  int explicitPlusRadV = getExplicitValence() + getNumRadicalElectrons();
  int chg = getFormalCharge();

  // NOTE: this is here to take care of the difference in element on
  // the right side of the carbon vs left side of carbon
  // For elements on the right side of the periodic table
  // (electronegative elements):
  //     NHYD = V - SBO + CHG
  // For elements on the left side of the periodic table
  // (electropositive elements):
  //      NHYD = V - SBO - CHG
  // This reflects that hydrogen adds to, for example, O as H+ while
  // it adds to Na as H-.

  // V = valence
  // SBO = Sum of bond orders
  // CHG = Formal charge

  //  It seems reasonable that the line is drawn at Carbon (in Group
  //  IV), but we must assume on which side of the line C
  //  falls... an assumption which will not always be correct.  For
  //  example:
  //  - Electropositive Carbon: a C with three singly-bonded
  //    neighbors (DV = 4, SBO = 3, CHG = 1) and a positive charge (a
  //    'stable' carbocation) should not have any hydrogens added.
  //  - Electronegative Carbon: C in isonitrile, R[N+]#[C-] (DV = 4, SBO = 3, 
  //    CHG = -1), also should not have any hydrogens added.
  //  Because isonitrile seems more relevant to pharma problems, we'll be
  //  making the second assumption:  *Carbon is electronegative*.
  //
  // So assuming you read all the above stuff - you know why we are
  // changing signs for "chg" here
  if ( isEarlyAtom(d_atomicNum) ) {
    chg *= -1;
  }

  // if we have an aromatic case treat it differently
  if (getIsAromatic()) {
    if (explicitPlusRadV <= (static_cast<int>(dv) + chg)) {
      res = dv + chg - explicitPlusRadV;
    }
    else {
      // As we assume when finding the explicitPlusRadValence if we are
      // aromatic we should not be adding any hydrogen and already
      // be at an accepted valence state,

      // FIX: this is just ERROR checking and probably moot - the
      // explicitPlusRadValence function called above should assure us that
      // we satisfy one of the accepted valence states for the
      // atom. The only diff I can think of is in the way we handle
      // formal charge here vs the explicit valence function.
      bool satis = false;
      for (INT_VECT_CI vi = valens.begin();
	   vi!=valens.end() && *vi>0; ++vi) {
        if (explicitPlusRadV == ((*vi) + chg)) {
          satis = true;
          break;
        }
      }
      if (strict && !satis) {
        std::ostringstream errout;
        errout << "Explicit valence for aromatic atom # " << getIdx() 
               << " not equal to any accepted valence\n";
        std::string msg = errout.str();
        BOOST_LOG(rdErrorLog) << msg << std::endl;
        throw MolSanitizeException(msg);
      }
      res = 0;
    }
  }
  else {
    // non-aromatic case we are allowed to have non default valences
    // and be able to add hydrogens
    res = -1;
    for (INT_VECT_CI vi = valens.begin();
	 vi != valens.end() && *vi>=0; ++vi) {
      int tot = (*vi) + chg;
      if (explicitPlusRadV <= tot) {
        res = tot - explicitPlusRadV;
        break;
      }
    }
    if (res < 0) {
      if(strict){
        // this means that the explicit valence is greater than any
        // allowed valence for the atoms - raise an error
        std::ostringstream errout;
        errout << "Explicit valence for atom # " << getIdx() 
               << " " << PeriodicTable::getTable()->getElementSymbol(d_atomicNum)
               << " greater than permitted";
        std::string msg = errout.str();
        BOOST_LOG(rdErrorLog) << msg << std::endl;
        throw MolSanitizeException(msg);
      } else {
        res = 0;
      }
    }
  }
      
  d_implicitValence = res;
  return res;
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
    if(this->dp_mol && what->dp_mol &&
       this->getOwningMol().getRingInfo()->isInitialized() &&
       what->getOwningMol().getRingInfo()->isInitialized() &&
       this->getOwningMol().getRingInfo()->numAtomRings(d_index) >
       what->getOwningMol().getRingInfo()->numAtomRings(what->d_index)){
      res=false;
    } else if(!this->getAtomicNum()){
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
          (this->getIsotope() && this->getIsotope()!=what->getIsotope()) ||
          (this->getNumRadicalElectrons() && this->getNumRadicalElectrons()!=what->getNumRadicalElectrons())
          ){
        res=false;
      }
    }
  }
  return res;
}
void Atom::updatePropertyCache(bool strict) {
  calcExplicitValence(strict);
  calcImplicitValence(strict);
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
} // end o' namespace RDKit

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
