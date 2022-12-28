//
//  Copyright (C) 2001-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <cmath>

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

bool isAromaticAtom(const Atom &atom) {
  if (atom.getIsAromatic()) {
    return true;
  }
  if (atom.hasOwningMol()) {
    for (const auto &bond : atom.getOwningMol().atomBonds(&atom)) {
      if (bond->getIsAromatic() ||
          bond->getBondType() == Bond::BondType::AROMATIC) {
        return true;
      }
    }
  }
  return false;
}

// Determine whether or not an element is to the left of carbon.
bool isEarlyAtom(int atomicNum) {
  static const bool table[119] = {
    false, // #0 *
    false, // #1 H
    false,  // #2 He
    true,  // #3 Li
    true,  // #4 Be
    true,  // #5 B
    false, // #6 C
    false, // #7 N
    false, // #8 O
    false, // #9 F
    false, // #10 Ne
    true,  // #11 Na
    true,  // #12 Mg
    true,  // #13 Al
    false, // #14 Si
    false, // #15 P
    false, // #16 S
    false, // #17 Cl
    false, // #18 Ar
    true,  // #19 K
    true,  // #20 Ca
    true,  // #21 Sc
    true,  // #22 Ti
    false, // #23 V
    false, // #24 Cr
    false, // #25 Mn
    false, // #26 Fe
    false, // #27 Co
    false, // #28 Ni
    false, // #29 Cu
    true,  // #30 Zn
    true,  // #31 Ga
    true,  // #32 Ge  see github #2606
    false, // #33 As
    false, // #34 Se
    false, // #35 Br
    false, // #36 Kr
    true,  // #37 Rb
    true,  // #38 Sr
    true,  // #39 Y
    true,  // #40 Zr
    true,  // #41 Nb
    false, // #42 Mo
    false, // #43 Tc
    false, // #44 Ru
    false, // #45 Rh
    false, // #46 Pd
    false, // #47 Ag
    true,  // #48 Cd
    true,  // #49 In
    true,  // #50 Sn  see github #2606
    true,  // #51 Sb  see github #2775
    false, // #52 Te
    false, // #53 I
    false, // #54 Xe
    true,  // #55 Cs
    true,  // #56 Ba
    true,  // #57 La
    true,  // #58 Ce
    true,  // #59 Pr
    true,  // #60 Nd
    true,  // #61 Pm
    false, // #62 Sm
    false, // #63 Eu
    false, // #64 Gd
    false, // #65 Tb
    false, // #66 Dy
    false, // #67 Ho
    false, // #68 Er
    false, // #69 Tm
    false, // #70 Yb
    false, // #71 Lu
    true,  // #72 Hf
    true,  // #73 Ta
    false, // #74 W
    false, // #75 Re
    false, // #76 Os
    false, // #77 Ir
    false, // #78 Pt
    false, // #79 Au
    true,  // #80 Hg
    true,  // #81 Tl
    true,  // #82 Pb  see github #2606
    true,  // #83 Bi  see github #2775
    false, // #84 Po
    false, // #85 At
    false, // #86 Rn
    true,  // #87 Fr
    true,  // #88 Ra
    true,  // #89 Ac
    true,  // #90 Th
    true,  // #91 Pa
    true,  // #92 U
    true,  // #93 Np
    false, // #94 Pu
    false, // #95 Am
    false, // #96 Cm
    false, // #97 Bk
    false, // #98 Cf
    false, // #99 Es
    false, // #100 Fm
    false, // #101 Md
    false, // #102 No
    false, // #103 Lr
    true,  // #104 Rf
    true,  // #105 Db
    true,  // #106 Sg
    true,  // #107 Bh
    true,  // #108 Hs
    true,  // #109 Mt
    true,  // #110 Ds
    true,  // #111 Rg
    true,  // #112 Cn
    true,  // #113 Nh
    true,  // #114 Fl
    true,  // #115 Mc
    true,  // #116 Lv
    true,  // #117 Ts
    true,  // #118 Og
  };
  return ((unsigned int)atomicNum < 119) && table[atomicNum];
}

Atom::Atom() : RDProps() {
  d_atomicNum = 0;
  initAtom();
}

Atom::Atom(unsigned int num) : RDProps() {
  d_atomicNum = num;
  initAtom();
};

Atom::Atom(const std::string &what) : RDProps() {
  d_atomicNum = PeriodicTable::getTable()->getAtomicNumber(what);
  initAtom();
};

void Atom::initFromOther(const Atom &other) {
  RDProps::operator=(other);
  // NOTE: we do *not* copy ownership!
  dp_mol = nullptr;
  d_atomicNum = other.d_atomicNum;
  d_index = 0;
  d_formalCharge = other.d_formalCharge;
  df_noImplicit = other.df_noImplicit;
  df_isAromatic = other.df_isAromatic;
  d_numExplicitHs = other.d_numExplicitHs;
  d_numRadicalElectrons = other.d_numRadicalElectrons;
  d_isotope = other.d_isotope;
  // d_pos = other.d_pos;
  d_chiralTag = other.d_chiralTag;
  d_hybrid = other.d_hybrid;
  d_implicitValence = other.d_implicitValence;
  d_explicitValence = other.d_explicitValence;
  if (other.dp_monomerInfo) {
    dp_monomerInfo = other.dp_monomerInfo->copy();
  } else {
    dp_monomerInfo = nullptr;
  }
}

Atom::Atom(const Atom &other) : RDProps() { initFromOther(other); }

Atom &Atom::operator=(const Atom &other) {
  if (this == &other) {
    return *this;
  }
  initFromOther(other);
  return *this;
}

void Atom::initAtom() {
  df_isAromatic = false;
  df_noImplicit = false;
  d_numExplicitHs = 0;
  d_numRadicalElectrons = 0;
  d_formalCharge = 0;
  d_index = 0;
  d_isotope = 0;
  d_chiralTag = CHI_UNSPECIFIED;
  d_hybrid = UNSPECIFIED;
  dp_mol = nullptr;
  dp_monomerInfo = nullptr;

  d_implicitValence = -1;
  d_explicitValence = -1;
}

Atom::~Atom() {
  if (dp_monomerInfo) {
    delete dp_monomerInfo;
  }
}

Atom *Atom::copy() const {
  auto *res = new Atom(*this);
  return res;
}

void Atom::setOwningMol(ROMol *other) {
  // NOTE: this operation does not update the topology of the owning
  // molecule (i.e. this atom is not added to the graph).  Only
  // molecules can add atoms to themselves.
  dp_mol = other;
}

std::string Atom::getSymbol() const {
  std::string res;
  // handle dummies differently:
  if (d_atomicNum != 0 ||
      !getPropIfPresent<std::string>(common_properties::dummyLabel, res)) {
    res = PeriodicTable::getTable()->getElementSymbol(d_atomicNum);
  }
  return res;
}

unsigned int Atom::getDegree() const {
  PRECONDITION(dp_mol,
               "degree not defined for atoms not associated with molecules");
  return getOwningMol().getAtomDegree(this);
}

unsigned int Atom::getTotalDegree() const {
  PRECONDITION(dp_mol,
               "degree not defined for atoms not associated with molecules");
  unsigned int res = this->getTotalNumHs(false) + this->getDegree();
  return res;
}

//
//  If includeNeighbors is set, we'll loop over our neighbors
//   and include any of them that are Hs in the count here
//
unsigned int Atom::getTotalNumHs(bool includeNeighbors) const {
  PRECONDITION(dp_mol,
               "valence not defined for atoms not associated with molecules")
  int res = getNumExplicitHs() + getNumImplicitHs();
  if (includeNeighbors) {
    for (auto nbr : getOwningMol().atomNeighbors(this)) {
      if (nbr->getAtomicNum() == 1) {
        ++res;
      }
    }
  }
  return res;
}

unsigned int Atom::getNumImplicitHs() const {
  if (df_noImplicit) {
    return 0;
  }

  PRECONDITION(d_implicitValence > -1,
               "getNumImplicitHs() called without preceding call to "
               "calcImplicitValence()");
  return getImplicitValence();
}

int Atom::getExplicitValence() const {
  PRECONDITION(dp_mol,
               "valence not defined for atoms not associated with molecules");
  PRECONDITION(
      d_explicitValence > -1,
      "getExplicitValence() called without call to calcExplicitValence()");
  return d_explicitValence;
}

unsigned int Atom::getTotalValence() const {
  PRECONDITION(dp_mol,
               "valence not defined for atoms not associated with molecules");
  return getExplicitValence() + getImplicitValence();
}

int Atom::calcExplicitValence(bool strict) {
  PRECONDITION(dp_mol,
               "valence not defined for atoms not associated with molecules");
  unsigned int res;
  // FIX: contributions of bonds to valence are being done at best
  // approximately
  double accum = 0;
  for (const auto bnd : getOwningMol().atomBonds(this)) {
    accum += bnd->getValenceContrib(this);
  }
  accum += getNumExplicitHs();

  // check accum is greater than the default valence
  unsigned int dv = PeriodicTable::getTable()->getDefaultValence(d_atomicNum);
  int chr = getFormalCharge();
  if (isEarlyAtom(d_atomicNum)) {
    chr *= -1;  // <- the usual correction for early atoms
  }
  // special case for carbon - see GitHub #539
  if (d_atomicNum == 6 && chr > 0) {
    chr = -chr;
  }
  if (accum > (dv + chr) && isAromaticAtom(*this)) {
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
    const INT_VECT &valens =
        PeriodicTable::getTable()->getValenceList(d_atomicNum);
    for (auto val : valens) {
      if (val == -1) {
        break;
      }
      val += chr;
      if (val > accum) {
        break;
      } else {
        pval = val;
      }
    }
    // if we're within 1.5 of the allowed valence, go ahead and take it.
    // this reflects things like the N in c1cccn1C, which starts with
    // accum of 4, but which can be kekulized to C1=CC=CN1C, where
    // the valence is 3 or the bridging N in c1ccn2cncc2c1, which starts
    // with a valence of 4.5, but can be happily kekulized down to a valence
    // of 3
    if (accum - pval <= 1.5) {
      accum = pval;
    }
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

  res = static_cast<int>(std::round(accum));

  if (strict) {
    int effectiveValence;
    if (PeriodicTable::getTable()->getNouterElecs(d_atomicNum) >= 4) {
      effectiveValence = res - getFormalCharge();
    } else {
      // for boron and co, we move to the right in the PT, so adding
      // extra valences means adding negative charge
      effectiveValence = res + getFormalCharge();
    }
    const INT_VECT &valens =
        PeriodicTable::getTable()->getValenceList(d_atomicNum);

    int maxValence = valens.back();
    // maxValence == -1 signifies that we'll take anything at the high end
    if (maxValence > 0 && effectiveValence > maxValence) {
      // the explicit valence is greater than any
      // allowed valence for the atoms - raise an error
      std::ostringstream errout;
      errout << "Explicit valence for atom # " << getIdx() << " "
             << PeriodicTable::getTable()->getElementSymbol(d_atomicNum) << ", "
             << effectiveValence << ", is greater than permitted";
      std::string msg = errout.str();
      BOOST_LOG(rdErrorLog) << msg << std::endl;
      throw AtomValenceException(msg, getIdx());
    }
  }
  d_explicitValence = res;

  return res;
}

int Atom::getImplicitValence() const {
  PRECONDITION(dp_mol,
               "valence not defined for atoms not associated with molecules");
  if (df_noImplicit) {
    return 0;
  }
  return d_implicitValence;
}

// NOTE: this uses the explicitValence, so it will call
// calcExplictValence() if it hasn't already been called
int Atom::calcImplicitValence(bool strict) {
  PRECONDITION(dp_mol,
               "valence not defined for atoms not associated with molecules");
  if (df_noImplicit) {
    return 0;
  }
  if (d_explicitValence == -1) {
    this->calcExplicitValence(strict);
  }
  // special cases
  if (d_atomicNum == 0) {
    d_implicitValence = 0;
    return 0;
  }
  for (const auto &nbri :
       boost::make_iterator_range(getOwningMol().getAtomBonds(this))) {
    const auto bnd = getOwningMol()[nbri];
    if (QueryOps::hasComplexBondTypeQuery(*bnd)) {
      d_implicitValence = 0;
      return 0;
    }
  }
  if (d_explicitValence == 0 && d_atomicNum == 1 &&
      d_numRadicalElectrons == 0) {
    if (d_formalCharge == 1 || d_formalCharge == -1) {
      d_implicitValence = 0;
      return 0;
    } else if (d_formalCharge == 0) {
      d_implicitValence = 1;
      return 1;
    } else {
      if (strict) {
        std::ostringstream errout;
        errout << "Unreasonable formal charge on hydrogen # " << getIdx()
               << ".";
        std::string msg = errout.str();
        BOOST_LOG(rdErrorLog) << msg << std::endl;
        throw AtomValenceException(msg, getIdx());
      } else {
        d_implicitValence = 0;
        return 0;
      }
    }
  }

  // this is basically the difference between the allowed valence of
  // the atom and the explicit valence already specified - tells how
  // many Hs to add
  //
  int res;

  // The d-block and f-block of the periodic table (i.e. transition metals,
  // lanthanoids and actinoids) have no default valence.
  int dv = PeriodicTable::getTable()->getDefaultValence(d_atomicNum);
  if (dv == -1) {
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
  // if "ev" is greater than all allowed valences for the atom raise an
  // exception
  // finally aromatic cases are dealt with differently - these atoms are allowed
  // only default valences
  const INT_VECT &valens =
      PeriodicTable::getTable()->getValenceList(d_atomicNum);

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
  if (isEarlyAtom(d_atomicNum)) {
    chg *= -1;
  }
  // special case for carbon - see GitHub #539
  if (d_atomicNum == 6 && chg > 0) {
    chg = -chg;
  }

  // if we have an aromatic case treat it differently
  if (isAromaticAtom(*this)) {
    if (explicitPlusRadV <= (static_cast<int>(dv) + chg)) {
      res = dv + chg - explicitPlusRadV;
    } else {
      // As we assume when finding the explicitPlusRadValence if we are
      // aromatic we should not be adding any hydrogen and already
      // be at an accepted valence state,

      // FIX: this is just ERROR checking and probably moot - the
      // explicitPlusRadValence function called above should assure us that
      // we satisfy one of the accepted valence states for the
      // atom. The only diff I can think of is in the way we handle
      // formal charge here vs the explicit valence function.
      bool satis = false;
      for (auto vi = valens.begin(); vi != valens.end() && *vi > 0; ++vi) {
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
        throw AtomValenceException(msg, getIdx());
      }
      res = 0;
    }
  } else {
    // non-aromatic case we are allowed to have non default valences
    // and be able to add hydrogens
    res = -1;
    for (auto vi = valens.begin(); vi != valens.end() && *vi >= 0; ++vi) {
      int tot = (*vi) + chg;
      if (explicitPlusRadV <= tot) {
        res = tot - explicitPlusRadV;
        break;
      }
    }
    if (res < 0) {
      if (strict && valens.back() != -1) {
        // this means that the explicit valence is greater than any
        // allowed valence for the atoms - raise an error
        std::ostringstream errout;
        errout << "Explicit valence for atom # " << getIdx() << " "
               << PeriodicTable::getTable()->getElementSymbol(d_atomicNum)
               << " greater than permitted";
        std::string msg = errout.str();
        BOOST_LOG(rdErrorLog) << msg << std::endl;
        throw AtomValenceException(msg, getIdx());
      } else {
        res = 0;
      }
    }
  }

  d_implicitValence = res;
  return res;
}

void Atom::setIsotope(unsigned int what) { d_isotope = what; }

double Atom::getMass() const {
  if (d_isotope) {
    double res =
        PeriodicTable::getTable()->getMassForIsotope(d_atomicNum, d_isotope);
    if (d_atomicNum != 0 && res == 0.0) {
      res = d_isotope;
    }
    return res;
  } else {
    return PeriodicTable::getTable()->getAtomicWeight(d_atomicNum);
  }
}

void Atom::setQuery(Atom::QUERYATOM_QUERY *) {
  //  Atoms don't have complex queries so this has to fail
  PRECONDITION(0, "plain atoms have no Query");
}
Atom::QUERYATOM_QUERY *Atom::getQuery() const { return nullptr; };
void Atom::expandQuery(Atom::QUERYATOM_QUERY *, Queries::CompositeQueryType,
                       bool) {
  PRECONDITION(0, "plain atoms have no Query");
}

bool Atom::Match(Atom const *what) const {
  PRECONDITION(what, "bad query atom");
  bool res = getAtomicNum() == what->getAtomicNum();

  // special dummy--dummy match case:
  //   [*] matches [*],[1*],[2*],etc.
  //   [1*] only matches [*] and [1*]
  if (res) {
    if (this->dp_mol && what->dp_mol &&
        this->getOwningMol().getRingInfo()->isInitialized() &&
        what->getOwningMol().getRingInfo()->isInitialized() &&
        this->getOwningMol().getRingInfo()->numAtomRings(d_index) >
            what->getOwningMol().getRingInfo()->numAtomRings(what->d_index)) {
      res = false;
    } else if (!this->getAtomicNum()) {
      // this is the new behavior, based on the isotopes:
      int tgt = this->getIsotope();
      int test = what->getIsotope();
      if (tgt && test && tgt != test) {
        res = false;
      }
    } else {
      // standard atom-atom match: The general rule here is that if this atom
      // has a property that
      // deviates from the default, then the other atom should match that value.
      if ((this->getFormalCharge() &&
           this->getFormalCharge() != what->getFormalCharge()) ||
          (this->getIsotope() && this->getIsotope() != what->getIsotope()) ||
          (this->getNumRadicalElectrons() &&
           this->getNumRadicalElectrons() != what->getNumRadicalElectrons())) {
        res = false;
      }
    }
  }
  return res;
}
void Atom::updatePropertyCache(bool strict) {
  calcExplicitValence(strict);
  calcImplicitValence(strict);
}

bool Atom::needsUpdatePropertyCache() const {
  return !(this->d_explicitValence >= 0 &&
           (this->df_noImplicit || this->d_implicitValence >= 0));
}

// returns the number of swaps required to convert the ordering
// of the probe list to match the order of our incoming bonds:
//
//  e.g. if our incoming bond order is: [0,1,2,3]:
//   getPerturbationOrder([1,0,2,3]) = 1
//   getPerturbationOrder([1,2,3,0]) = 3
//   getPerturbationOrder([1,2,0,3]) = 2
int Atom::getPerturbationOrder(const INT_LIST &probe) const {
  PRECONDITION(
      dp_mol,
      "perturbation order not defined for atoms not associated with molecules")
  INT_LIST ref;
  ROMol::OEDGE_ITER beg, end;
  boost::tie(beg, end) = getOwningMol().getAtomBonds(this);
  while (beg != end) {
    ref.push_back(getOwningMol()[*beg]->getIdx());
    ++beg;
  }
  int nSwaps = static_cast<int>(countSwapsToInterconvert(probe, ref));
  return nSwaps;
}

static const unsigned char octahedral_invert[31] = {
    0,   //  0 -> 0
    2,   //  1 -> 2
    1,   //  2 -> 1
    16,  //  3 -> 16
    14,  //  4 -> 14
    15,  //  5 -> 15
    18,  //  6 -> 18
    17,  //  7 -> 17
    10,  //  8 -> 10
    11,  //  9 -> 11
    8,   // 10 -> 8
    9,   // 11 -> 9
    13,  // 12 -> 13
    12,  // 13 -> 12
    4,   // 14 -> 4
    5,   // 15 -> 5
    3,   // 16 -> 3
    7,   // 17 -> 7
    6,   // 18 -> 6
    24,  // 19 -> 24
    23,  // 20 -> 23
    22,  // 21 -> 22
    21,  // 22 -> 21
    20,  // 23 -> 20
    19,  // 24 -> 19
    30,  // 25 -> 30
    29,  // 26 -> 29
    28,  // 27 -> 28
    27,  // 28 -> 27
    26,  // 29 -> 26
    25   // 30 -> 25
};

static const unsigned char trigonalbipyramidal_invert[21] = {
    0,   //  0 -> 0
    2,   //  1 -> 2
    1,   //  2 -> 1
    4,   //  3 -> 4
    3,   //  4 -> 3
    6,   //  5 -> 6
    5,   //  6 -> 5
    8,   //  7 -> 8
    7,   //  8 -> 7
    11,  //  9 -> 11
    12,  // 10 -> 12
    9,   // 11 -> 9
    10,  // 12 -> 10
    14,  // 13 -> 14
    13,  // 14 -> 13
    20,  // 15 -> 20
    19,  // 16 -> 19
    18,  // 17 -> 28
    17,  // 18 -> 17
    16,  // 19 -> 16
    15   // 20 -> 15
};

bool Atom::invertChirality() {
  unsigned int perm;
  switch (getChiralTag()) {
    case CHI_TETRAHEDRAL_CW:
      setChiralTag(CHI_TETRAHEDRAL_CCW);
      return true;
    case CHI_TETRAHEDRAL_CCW:
      setChiralTag(CHI_TETRAHEDRAL_CW);
      return true;
    case CHI_TETRAHEDRAL:
      if (getPropIfPresent(common_properties::_chiralPermutation, perm)) {
        if (perm == 1) {
          perm = 2;
        } else if (perm == 2) {
          perm = 1;
        } else {
          perm = 0;
        }
        setProp(common_properties::_chiralPermutation, perm);
        return perm != 0;
      }
      break;
    case CHI_TRIGONALBIPYRAMIDAL:
      if (getPropIfPresent(common_properties::_chiralPermutation, perm)) {
        perm = (perm <= 20) ? trigonalbipyramidal_invert[perm] : 0;
        setProp(common_properties::_chiralPermutation, perm);
        return perm != 0;
      }
      break;
    case CHI_OCTAHEDRAL:
      if (getPropIfPresent(common_properties::_chiralPermutation, perm)) {
        perm = (perm <= 30) ? octahedral_invert[perm] : 0;
        setProp(common_properties::_chiralPermutation, perm);
        return perm != 0;
      }
      break;
    default:
      break;
  }
  return false;
}

void setAtomRLabel(Atom *atm, int rlabel) {
  PRECONDITION(atm, "bad atom");
  // rlabel ==> n2 => 0..99
  PRECONDITION(rlabel >= 0 && rlabel < 100,
               "rlabel out of range for MDL files");
  if (rlabel) {
    atm->setProp(common_properties::_MolFileRLabel,
                 static_cast<unsigned int>(rlabel));
  } else if (atm->hasProp(common_properties::_MolFileRLabel)) {
    atm->clearProp(common_properties::_MolFileRLabel);
  }
}
//! Gets the atom's RLabel
int getAtomRLabel(const Atom *atom) {
  PRECONDITION(atom, "bad atom");
  unsigned int rlabel = 0;
  atom->getPropIfPresent(common_properties::_MolFileRLabel, rlabel);
  return static_cast<int>(rlabel);
}

void setAtomAlias(Atom *atom, const std::string &alias) {
  PRECONDITION(atom, "bad atom");
  if (alias != "") {
    atom->setProp(common_properties::molFileAlias, alias);
  } else if (atom->hasProp(common_properties::molFileAlias)) {
    atom->clearProp(common_properties::molFileAlias);
  }
}

std::string getAtomAlias(const Atom *atom) {
  PRECONDITION(atom, "bad atom");
  std::string alias;
  atom->getPropIfPresent(common_properties::molFileAlias, alias);
  return alias;
}

void setAtomValue(Atom *atom, const std::string &value) {
  PRECONDITION(atom, "bad atom");
  if (value != "") {
    atom->setProp(common_properties::molFileValue, value);
  } else if (atom->hasProp(common_properties::molFileValue)) {
    atom->clearProp(common_properties::molFileValue);
  }
}

std::string getAtomValue(const Atom *atom) {
  PRECONDITION(atom, "bad atom");
  std::string value;
  atom->getPropIfPresent(common_properties::molFileValue, value);
  return value;
}

void setSupplementalSmilesLabel(Atom *atom, const std::string &label) {
  PRECONDITION(atom, "bad atom");
  if (label != "") {
    atom->setProp(common_properties::_supplementalSmilesLabel, label);
  } else if (atom->hasProp(common_properties::_supplementalSmilesLabel)) {
    atom->clearProp(common_properties::_supplementalSmilesLabel);
  }
}

std::string getSupplementalSmilesLabel(const Atom *atom) {
  PRECONDITION(atom, "bad atom");
  std::string label;
  atom->getPropIfPresent(common_properties::_supplementalSmilesLabel, label);
  return label;
}

}  // namespace RDKit

std::ostream &operator<<(std::ostream &target, const RDKit::Atom &at) {
  target << at.getIdx() << " " << at.getAtomicNum() << " " << at.getSymbol();
  target << " chg: " << at.getFormalCharge();
  target << "  deg: " << at.getDegree();
  target << " exp: ";
  try {
    int explicitValence = at.getExplicitValence();
    target << explicitValence;
  } catch (...) {
    target << "N/A";
  }
  target << " imp: ";
  try {
    int implicitValence = at.getImplicitValence();
    target << implicitValence;
  } catch (...) {
    target << "N/A";
  }
  target << " hyb: " << at.getHybridization();
  target << " arom?: " << at.getIsAromatic();
  target << " chi: " << at.getChiralTag();
  if (at.getNumRadicalElectrons()) {
    target << " rad: " << at.getNumRadicalElectrons();
  }
  if (at.getIsotope()) {
    target << " iso: " << at.getIsotope();
  }
  if (at.getAtomMapNum()) {
    target << " mapno: " << at.getAtomMapNum();
  }
  return target;
};
