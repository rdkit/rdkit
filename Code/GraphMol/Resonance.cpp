//
//
//  Copyright (C) 2015 Paolo Tosco
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Resonance.h>
#include <RDGeneral/hash/hash.hpp>
#include <RDGeneral/RDThreads.h>
#ifdef RDK_BUILD_THREADSAFE_SSS
#include <thread>
#include <future>
#endif
#include <algorithm>

namespace RDKit {
// class definitions that do not need being exposed in Resonance.h

typedef std::set<std::size_t> CESet;
typedef std::unordered_map<std::size_t, unsigned int> CEDegCount;

class CEVect2 {
 public:
  CEVect2(const CEMap &ceMap);
  ConjElectrons *getCE(unsigned int depth, unsigned int width);
  size_t ceCount() { return d_ceVect.size(); }
  size_t depth() { return d_degVect.size(); }
  void resize(size_t size);
  void idxToDepthWidth(size_t idx, unsigned int &d, unsigned int &w);
  unsigned int ceCountAtDepth(unsigned int depth);
  unsigned int ceCountUntilDepth(unsigned int depth);

 private:
  static bool resonanceStructureCompare(const ConjElectrons *a,
                                        const ConjElectrons *b);
  CEVect d_ceVect;
  std::vector<unsigned int> d_degVect;
};

class CEMetrics {
  friend class ConjElectrons;

 public:
  CEMetrics(){};
  bool operator==(const CEMetrics &other) const {
    return (d_hash == other.d_hash);
  }
  bool operator!=(const CEMetrics &other) const {
    return (d_hash != other.d_hash);
  }

 private:
  int d_absFormalCharges{0};
  int d_fcSameSignDist{0};
  int d_fcOppSignDist{0};
  int d_nbMissing{0};
  int d_wtdFormalCharges{0};
  int d_sumFormalChargeIdxs{0};
  int d_sumMultipleBondIdxs{0};
  std::size_t d_hash{0};
};

struct CEStats {
  CEStats()
      : d_init(false),
        d_haveNoCationsRightOfN(false),
        d_haveNoChargeSeparation(false),
        d_minNbMissing(0) {}
  bool d_init;
  bool d_haveNoCationsRightOfN;
  bool d_haveNoChargeSeparation;
  unsigned int d_minNbMissing;
};

class ConjElectrons {
 public:
  typedef enum {
    HAVE_CATION_RIGHT_OF_N = (1 << 0),
    HAVE_CATION = (1 << 1),
    HAVE_ANION = (1 << 2)
  } ConjElectronsFlags;
  typedef enum { FP_BONDS = (1 << 0), FP_ATOMS = (1 << 1) } FPFlags;
  ConjElectrons(ResonanceMolSupplier *parent, unsigned int groupIdx);
  ConjElectrons(const ConjElectrons &ce);
  ~ConjElectrons();
  unsigned int groupIdx() const { return d_groupIdx; };
  unsigned int currElectrons() const { return d_currElectrons; };
  unsigned int totalElectrons() const { return d_totalElectrons; };
  void decrCurrElectrons(unsigned int d);
  AtomElectrons *getAtomElectronsWithIdx(unsigned int ai);
  BondElectrons *getBondElectronsWithIdx(unsigned int bi);
  void pushToBeginStack(unsigned int ai);
  bool popFromBeginStack(unsigned int &ai);
  bool isBeginStackEmpty() { return d_beginAIStack.empty(); };
  int allowedChgLeftOfN() const { return d_allowedChgLeftOfN; };
  void decrAllowedChgLeftOfN(int d) { d_allowedChgLeftOfN -= d; };
  int totalFormalCharge() const { return d_totalFormalCharge; };
  bool hasCationRightOfN() const {
    return static_cast<bool>(d_flags & HAVE_CATION_RIGHT_OF_N);
  };
  bool hasChargeSeparation() const {
    return static_cast<bool>((d_flags & HAVE_CATION) && (d_flags & HAVE_ANION));
  };
  unsigned int absFormalCharges() const {
    return d_ceMetrics.d_absFormalCharges;
  };
  unsigned int fcSameSignDist() const { return d_ceMetrics.d_fcSameSignDist; };
  unsigned int fcOppSignDist() const { return d_ceMetrics.d_fcOppSignDist; };
  unsigned int nbMissing() const { return d_ceMetrics.d_nbMissing; }
  unsigned int sumFormalChargeIdxs() const {
    return d_ceMetrics.d_sumFormalChargeIdxs;
  }
  unsigned int sumMultipleBondIdxs() const {
    return d_ceMetrics.d_sumMultipleBondIdxs;
  }
  std::size_t hash() { return d_ceMetrics.d_hash; }
  int wtdFormalCharges() const { return d_ceMetrics.d_wtdFormalCharges; };
  void enumerateNonBonded(CEMap &ceMap, CEDegCount &ceDegCount,
                          CEStats &ceStats);
  void initCeFromMol();
  void assignNonBonded();
  void assignFormalCharge();
  bool assignFormalChargesAndStore(CEMap &ceMap, CEDegCount &ceDegCount,
                                   CEStats &ceStats, unsigned int fpFlags);
  void assignBondsFormalChargesToMol(ROMol &mol);
  bool checkChargesAndBondOrders();
  void computeMetrics();
  bool checkMetrics(CEStats &ceStats, bool &ok) const;
  bool purgeMaps(CEMap &ceMap, CEDegCount &ceDegCount, CEStats &ceStats) const;
  void updateDegCount(CEDegCount &ceDegCount);
  std::size_t computeFP(unsigned int flags);
  inline bool haveFP(CESet &ceSet, unsigned int flags);
  inline bool storeFP(CEMap &ceMap, unsigned int flags);
  ResonanceMolSupplier *parent() const { return d_parent; };

 private:
  unsigned int d_groupIdx;
  unsigned int d_totalElectrons;
  unsigned int d_currElectrons;
  unsigned int d_numFormalCharges;
  int d_totalFormalCharge;
  int d_allowedChgLeftOfN;
  std::uint8_t d_flags;
  CEMetrics d_ceMetrics;
  ConjBondMap d_conjBondMap;
  ConjAtomMap d_conjAtomMap;
  std::stack<unsigned int> d_beginAIStack;
  ResonanceMolSupplier *d_parent;
  ConjElectrons &operator=(const ConjElectrons &);
  unsigned int countTotalElectrons();
  void computeDistFormalCharges();
  void computeSumFormalChargeIdxs();
  void computeSumMultipleBondIdxs();
  void checkOctets();
};

class CEVect2Store {
 public:
  CEVect2Store(CEVect &ceVect);
  unsigned int ceCount() { return d_ceCount; }
  ConjElectrons *getCE(unsigned int i, unsigned int j);

 private:
  unsigned int d_ceCount;
  CEVect2 d_ceVect2;
};

class AtomElectrons {
 public:
  typedef enum {
    LAST_BOND = (1 << 0),
    DEFINITIVE = (1 << 1),
    STACKED = (1 << 2),
    HAS_MULTIPLE_BOND = (1 << 3),
  } AtomElectronsFlags;
  typedef enum { NEED_CHARGE_BIT = 1 } AllowedBondFlag;
  AtomElectrons(ConjElectrons *parent, const Atom *a);
  AtomElectrons(ConjElectrons *parent, const AtomElectrons &ae);
  ~AtomElectrons() = default;
  std::uint8_t findAllowedBonds(unsigned int bi);
  bool hasOctet() const { return ((d_nb + d_tv * 2) == 8); };
  bool isLastBond() const { return (d_flags & LAST_BOND); };
  void setLastBond() { d_flags |= LAST_BOND; };
  bool isDefinitive() const { return (d_flags & DEFINITIVE); };
  void setDefinitive() { d_flags |= DEFINITIVE; };
  bool hasMultipleBond() const { return (d_flags & HAS_MULTIPLE_BOND); };
  void setHasMultipleBond() { d_flags |= HAS_MULTIPLE_BOND; };
  bool isStacked() const { return (d_flags & STACKED); };
  void setStacked() { d_flags |= STACKED; };
  void clearStacked() { d_flags &= ~STACKED; };
  int conjGrpIdx() const {
    return d_parent->parent()->getAtomConjGrpIdx(d_atom->getIdx());
  };
  void finalizeAtom();
  unsigned int nb() const { return d_nb; };
  unsigned int tv() const { return d_tv; };
  unsigned int oe() const {
    return PeriodicTable::getTable()->getNouterElecs(d_atom->getAtomicNum());
  };
  int fc() const { return d_fc; };
  void tvIncr(unsigned int i) { d_tv += i; };
  unsigned int neededNbForOctet() const { return (8 - (2 * d_tv + d_nb)); }
  const Atom *atom() { return d_atom; }
  void initTvNbFcFromAtom();
  void assignNonBonded(unsigned int nb) { d_nb = nb; }
  void assignFormalCharge() { d_fc = oe() - (d_nb + d_tv); }
  bool isNbrCharged(unsigned int bo, unsigned int oeConstraint = 0);
  AtomElectrons &operator=(const AtomElectrons &) = delete;
  AtomElectrons(const AtomElectrons &) = delete;

 private:
  std::uint8_t d_nb;
  std::uint8_t d_tv;
  std::int8_t d_fc;
  std::uint8_t d_flags;
  const Atom *d_atom;
  ConjElectrons *d_parent;
  std::uint8_t canAddBondWithOrder(unsigned int bo);
  void allConjBondsDefinitiveBut(unsigned int bi);
};

class BondElectrons {
 public:
  typedef enum { DEFINITIVE = (1 << 0) } BondElectronsFlags;
  BondElectrons(ConjElectrons *parent, const Bond *b);
  BondElectrons(ConjElectrons *parent, const BondElectrons &be);
  ~BondElectrons() = default;
  bool isDefinitive() const { return (d_flags & DEFINITIVE); };
  void setDefinitive() { d_flags |= DEFINITIVE; };
  int conjGrpIdx() const {
    return d_parent->parent()->getBondConjGrpIdx(d_bond->getIdx());
  };
  void setOrder(unsigned int bo);
  unsigned int order() const { return d_bo; };
  unsigned int orderFromBond();
  void initOrderFromBond() { d_bo = orderFromBond(); };
  const Bond *bond() { return d_bond; };
  BondElectrons &operator=(const BondElectrons &) = delete;
  BondElectrons(const BondElectrons &) = delete;

 private:
  std::uint8_t d_bo;
  std::uint8_t d_flags;
  const Bond *d_bond;
  ConjElectrons *d_parent;
};

namespace ResonanceUtils {
// depending on the number of atoms which won't have a complete
// octet (nCandSlots) and the number of those which need non-bonded
// electrons (nTotalSlots), the number of permutations (numComb) and
// a binary code (v) which indicates which of the atom indices in
// aiVec will be octet-unsatisfied for each permutation are computed
void getNumCombStartV(unsigned int nCandSlots, unsigned int nTotalSlots,
                      unsigned int &numComb, unsigned int &v) {
  numComb = 1;
  v = 0;
  for (unsigned int i = 0; i < nCandSlots; ++i) {
    numComb *= (nTotalSlots - i);
    numComb /= (i + 1);
    v |= (1 << i);
  }
}

// get the next permutation
void updateV(unsigned int &v) {
  unsigned int t = (v | (v - 1)) + 1;
  v = t | ((((t & (~t + 1)) / (v & (~v + 1))) >> 1) - 1);
}

// sanitize the resonance structure which has been assembled
void sanitizeMol(RWMol &mol) {
  unsigned int opFailed;
  MolOps::sanitizeMol(
      mol, opFailed, MolOps::SANITIZE_FINDRADICALS | MolOps::SANITIZE_ADJUSTHS);
}

// fix the number of explicit and implicit Hs in the
// resonance structure which has been assembled
void fixExplicitImplicitHs(ROMol &mol) {
  mol.clearComputedProps(false);
  for (ROMol::AtomIterator ai = mol.beginAtoms(); ai != mol.endAtoms(); ++ai) {
    (*ai)->clearComputedProps();
    (*ai)->setNumExplicitHs((*ai)->getNumImplicitHs() +
                            (*ai)->getNumExplicitHs());
    (*ai)->updatePropertyCache();
  }
}
}  // end of namespace ResonanceUtils

// object constructor
AtomElectrons::AtomElectrons(ConjElectrons *parent, const Atom *a)
    : d_nb(0), d_fc(0), d_flags(0), d_atom(a), d_parent(parent) {
  PRECONDITION(d_atom, "d_atom cannot be NULL");
  d_tv = static_cast<std::uint8_t>(a->getTotalDegree());
  const ROMol &mol = d_atom->getOwningMol();
  for (const auto &bNbri :
       boost::make_iterator_range(mol.getAtomBonds(d_atom))) {
    if (d_tv && mol[bNbri]->getBondType() >= Bond::DATIVEONE) {
      --d_tv;
    }
  }
}

// copy constructor
AtomElectrons::AtomElectrons(ConjElectrons *parent, const AtomElectrons &ae)
    : d_nb(ae.d_nb),
      d_tv(ae.d_tv),
      d_fc(ae.d_fc),
      d_flags(ae.d_flags),
      d_atom(ae.d_atom),
      d_parent(parent) {}

// assign total valence, formal charge and non-bonded electrons
// from the original atom
void AtomElectrons::initTvNbFcFromAtom() {
  d_tv = d_atom->getTotalValence();
  d_fc = d_atom->getFormalCharge();
  d_nb = oe() - d_tv - d_fc;
}

std::uint8_t AtomElectrons::findAllowedBonds(unsigned int bi) {
  // AtomElectrons::findAllowedBonds returns a 6-bit result
  // encoded as follows:
  // +-------------------------------------+
  // | BIT | MEANING                       |
  // |   1 | can accept single bond        |
  // |   2 | needs charge if single bonded |
  // |   3 | can accept double bond        |
  // |   4 | needs charge if double bonded |
  // |   5 | can accept triple bond        |
  // |   6 | needs charge if triple bonded |
  // +-------------------------------------+

  allConjBondsDefinitiveBut(bi);
  std::uint8_t res = 0;
  for (unsigned int i = 0; i < 3; ++i) {
    res |= (canAddBondWithOrder(i + 1) << (i * 2));
  }
  return res;
}

// returns true if any conjugated neighbor is charged (and
// has oe() == oeConstraint, if the latter is non-zero)
bool AtomElectrons::isNbrCharged(unsigned int bo, unsigned int oeConstraint) {
  bool res = false;
  const ROMol &mol = d_parent->parent()->mol();
  ROMol::OEDGE_ITER nbrIdx, endNbrs;
  boost::tie(nbrIdx, endNbrs) = mol.getAtomBonds(d_atom);
  for (; !res && (nbrIdx != endNbrs); ++nbrIdx) {
    const Bond *bondNbr = mol[*nbrIdx];
    unsigned int biNbr = bondNbr->getIdx();
    if (d_parent->parent()->getBondConjGrpIdx(biNbr) != conjGrpIdx()) {
      continue;
    }
    BondElectrons *beNbr = d_parent->getBondElectronsWithIdx(biNbr);
    if (!beNbr) {
      continue;
    }
    const Atom *atomNbr = bondNbr->getOtherAtom(d_atom);
    unsigned int aiNbr = atomNbr->getIdx();
    AtomElectrons *aeNbr = d_parent->getAtomElectronsWithIdx(aiNbr);
    if (!aeNbr) {
      continue;
    }
    res = (((beNbr->isDefinitive() && !aeNbr->hasOctet()) ||
            (!beNbr->isDefinitive() && aeNbr->isDefinitive() &&
             (aeNbr->oe() < (5 - bo)))) &&
           (!oeConstraint || (aeNbr->oe() == oeConstraint)));
  }
  return res;
}

// returns a 2-bit value, where the least significant bit is true
// if the atom can add a bond with order bo, and the most significant
// bit is true if the atom needs be charged
std::uint8_t AtomElectrons::canAddBondWithOrder(unsigned int bo) {
  std::uint8_t canAdd = !isDefinitive();
  if (canAdd) {
    canAdd = (d_tv <= (5 - bo));
  }
  // if canAdd is true, loop over neighboring conjugated bonds
  // and check their definitive flag; if all neighbors are
  // definitive, we need an additional check on total valence
  if (canAdd && isLastBond()) {
    // we allow a formal charge up to 2 on atoms right of
    // carbon, not more than 1 on atoms left of N
    unsigned int rightOfC = ((oe() > 4) ? 1 : 0);
    unsigned int fcInc = 0;
    if (rightOfC) {
      fcInc = (!isNbrCharged(bo, 4) ? 1 : 0);
    } else {
      // atoms left of N can be charged only if:
      // - the neighbor is uncharged and either the conjugate
      //   group does not bear a non-zero total formal charge
      //   or it does and there is no other element left of N
      //   which has already been assigned a formal charge
      // - the neighbor is charged, and triple-bonded to
      //   this atom, which is left of N (e.g., isonitrile)
      bool isnc = isNbrCharged(bo);
      fcInc = (((!isnc && !(!d_parent->allowedChgLeftOfN() &&
                            d_parent->totalFormalCharge())) ||
                (isnc && (bo == 3) && (oe() < 5)))
                   ? 1
                   : 0);
    }
    unsigned int e = oe() + d_tv - 1 + bo;
    canAdd = ((e + fcInc + rightOfC) >= 8);
    if (canAdd && (e < 8)) {
      if (d_parent->totalFormalCharge() != 0 ||
          (d_parent->parent()->flags() &
           (ResonanceMolSupplier::UNCONSTRAINED_CATIONS |
            ResonanceMolSupplier::UNCONSTRAINED_ANIONS)) ||
          rightOfC) {
        canAdd |= (1 << NEED_CHARGE_BIT);
      } else {
        canAdd = 0;
      }
    }
  }
  return canAdd;
}

// sets the LAST_BOND flag on this atom if there is only one
// non-definitive bond left
void AtomElectrons::allConjBondsDefinitiveBut(unsigned int bi) {
  bool allDefinitive = true;
  ROMol &mol = d_atom->getOwningMol();
  ROMol::OEDGE_ITER nbrIdx, endNbrs;
  boost::tie(nbrIdx, endNbrs) = mol.getAtomBonds(d_atom);
  for (; allDefinitive && (nbrIdx != endNbrs); ++nbrIdx) {
    unsigned int nbi = mol[*nbrIdx]->getIdx();
    if ((nbi != bi) &&
        (d_parent->parent()->getBondConjGrpIdx(nbi) == conjGrpIdx())) {
      allDefinitive = d_parent->getBondElectronsWithIdx(nbi)->isDefinitive();
    }
  }
  if (allDefinitive) {
    setLastBond();
  }
}

// called after all bonds to the atom have been marked as definitive
void AtomElectrons::finalizeAtom() {
  // if the atom is left of N and needs non-bonded electrons
  // to achieve the octet, it is the total formal charge
  // of the conjugated group which determines if it is
  // going to be a cation or an anion; once the formal
  // charge is assigned to an atom left of N, the counter of allowed
  // charges on atoms left of N (signed) is decremented by one
  if (oe() < 5) {
    unsigned int nb = neededNbForOctet();
    if (nb) {
      int fc = nb / 2;
      if (d_parent->allowedChgLeftOfN()) {
        if (d_parent->allowedChgLeftOfN() < 0) {
          fc = -fc;
        }
        d_parent->decrAllowedChgLeftOfN(fc);
      }
    }
  }
}

// object constructor
BondElectrons::BondElectrons(ConjElectrons *parent, const Bond *b)
    : d_bo(1), d_flags(0), d_bond(b), d_parent(parent) {
  PRECONDITION(d_bond, "d_bond cannot be NULL");
}

// copy constructor
BondElectrons::BondElectrons(ConjElectrons *parent, const BondElectrons &be)
    : d_bo(be.d_bo), d_flags(be.d_flags), d_bond(be.d_bond), d_parent(parent) {}

// returns the bond order given the bond type
unsigned int BondElectrons::orderFromBond() {
  unsigned int bo = 0;
  switch (d_bond->getBondType()) {
    case Bond::SINGLE:
      bo = 1;
      break;
    case Bond::DOUBLE:
      bo = 2;
      break;
    case Bond::TRIPLE:
      bo = 3;
      break;
    default:
      std::stringstream ss;
      ss << "Bond idx " << d_bond->getIdx() << " between atoms "
         << d_bond->getBeginAtomIdx() << " and " << d_bond->getEndAtomIdx()
         << " has an invalid bond type";
      throw std::runtime_error(ss.str());
  }
  return bo;
}

// set bond order, update total valence on the atoms involved
// in the bond and update count of current available electrons
void BondElectrons::setOrder(unsigned int bo) {
  d_parent->getAtomElectronsWithIdx(d_bond->getBeginAtomIdx())->tvIncr(bo - 1);
  d_parent->getAtomElectronsWithIdx(d_bond->getEndAtomIdx())->tvIncr(bo - 1);
  setDefinitive();
  d_parent->decrCurrElectrons(bo * 2);
  d_bo = bo;
}

// object constructor
ConjElectrons::ConjElectrons(ResonanceMolSupplier *parent,
                             unsigned int groupIdx)
    : d_groupIdx(groupIdx),
      d_totalElectrons(0),
      d_numFormalCharges(0),
      d_totalFormalCharge(0),
      d_flags(0),
      d_parent(parent) {
  const ROMol &mol = d_parent->mol();
  unsigned int nb = mol.getNumBonds();
  unsigned int na = mol.getNumAtoms();
  for (unsigned int ai = 0; ai < na; ++ai) {
    if (d_parent->getAtomConjGrpIdx(ai) != -1) {
      d_totalFormalCharge += mol.getAtomWithIdx(ai)->getFormalCharge();
    }
  }
  d_allowedChgLeftOfN = d_totalFormalCharge;
  for (unsigned int bi = 0; bi < nb; ++bi) {
    if (d_parent->getBondConjGrpIdx(bi) != static_cast<int>(groupIdx)) {
      continue;
    }
    const Bond *bond = mol.getBondWithIdx(bi);
    // store the pointers to BondElectrons objects in a map
    d_conjBondMap[bi] = new BondElectrons(this, bond);
    // store the pointers to AtomElectrons objects in a map
    const Atom *atom[2] = {bond->getBeginAtom(), bond->getEndAtom()};
    for (auto &i : atom) {
      unsigned int ai = i->getIdx();
      if (d_conjAtomMap.find(ai) == d_conjAtomMap.end()) {
        d_conjAtomMap[ai] = new AtomElectrons(this, i);
      }
    }
  }
  // count total number of valence electrons in conjugated group
  d_currElectrons = countTotalElectrons();
}

// copy constructor
ConjElectrons::ConjElectrons(const ConjElectrons &ce)
    : d_groupIdx(ce.d_groupIdx),
      d_totalElectrons(ce.d_totalElectrons),
      d_currElectrons(ce.d_currElectrons),
      d_numFormalCharges(ce.d_numFormalCharges),
      d_totalFormalCharge(ce.d_totalFormalCharge),
      d_allowedChgLeftOfN(ce.d_allowedChgLeftOfN),
      d_flags(ce.d_flags),
      d_ceMetrics(ce.d_ceMetrics),
      d_beginAIStack(ce.d_beginAIStack),
      d_parent(ce.d_parent) {
  for (const auto &it : ce.d_conjAtomMap) {
    d_conjAtomMap[it.first] = new AtomElectrons(this, *(it.second));
  }
  for (const auto &it : ce.d_conjBondMap) {
    d_conjBondMap[it.first] = new BondElectrons(this, *(it.second));
  }
}

// object destructor
ConjElectrons::~ConjElectrons() {
  for (ConjAtomMap::const_iterator it = d_conjAtomMap.begin();
       it != d_conjAtomMap.end(); ++it) {
    delete it->second;
  }
  for (ConjBondMap::const_iterator it = d_conjBondMap.begin();
       it != d_conjBondMap.end(); ++it) {
    delete it->second;
  }
}

std::size_t ConjElectrons::computeFP(unsigned int flags) {
  std::uint8_t byte;
  ConjFP fp;
  size_t fpSize = 0;
  if (flags & FP_ATOMS) {
    fpSize += d_conjAtomMap.size();
  }
  if (flags & FP_BONDS) {
    fpSize += (d_conjBondMap.size() - 1) / 4 + 1;
  }
  fp.reserve(fpSize);
  if (flags & FP_ATOMS) {
    // for each atom, we push a byte to the FP vector whose
    // 4 least significant bits are total valence and the
    // 4 most significant bits are non-bonded electrons
    for (ConjAtomMap::const_iterator it = d_conjAtomMap.begin();
         it != d_conjAtomMap.end(); ++it) {
      byte = it->second->tv() | (it->second->nb() << 4);
      fp.push_back(byte);
    }
  }
  if (flags & FP_BONDS) {
    unsigned int i = 0;
    byte = 0;
    for (ConjBondMap::const_iterator it = d_conjBondMap.begin();
         it != d_conjBondMap.end(); ++it) {
      // for each bond, we push 2 bits to the FP vector which
      // represent the bond order; the FP vector is byte-aligned
      // anyway
      if (i && !(i % 4)) {
        fp.push_back(byte);
        byte = 0;
        i = 0;
      }
      byte |= (static_cast<std::uint8_t>(it->second->order()) << (i * 2));
      ++i;
    }
    if (i) {
      fp.push_back(byte);
    }
  }
  // convert the FP vector to a hash
  return boost::hash_range(fp.begin(), fp.end());
}

// store fingerprints for this ConjElectrons object in ceSet
// return true if the FP did not already exist in the set, false
// if they did (so the ConjElectrons object can be deleted)
bool ConjElectrons::haveFP(CESet &ceSet, unsigned int flags) {
  // return true if the FP did not already exist in ceMap,
  // false if it did
  return ceSet.insert(computeFP(flags)).second;
}

// store fingerprints for this ConjElectrons object in ceMap
// return true if the FP did not already exist in the map, false
// if they did (so the ConjElectrons object can be deleted)
bool ConjElectrons::storeFP(CEMap &ceMap, unsigned int flags) {
  // return true if the FP did not already exist in ceMap,
  // false if it did
  return ceMap.insert(std::make_pair(computeFP(flags), this)).second;
}

// assign bond orders and formal charges as encoded in the
// ConjElectrons object to the ROMol passed as reference
void ConjElectrons::assignBondsFormalChargesToMol(ROMol &mol) {
  const Bond::BondType bondType[3] = {Bond::SINGLE, Bond::DOUBLE, Bond::TRIPLE};
  for (ConjAtomMap::const_iterator it = d_conjAtomMap.begin();
       it != d_conjAtomMap.end(); ++it) {
    unsigned int ai = it->first;
    AtomElectrons *ae = it->second;
    mol.getAtomWithIdx(ai)->setFormalCharge(ae->fc());
  }
  for (ConjBondMap::const_iterator it = d_conjBondMap.begin();
       it != d_conjBondMap.end(); ++it) {
    unsigned int bi = it->first;
    BondElectrons *be = it->second;
    if ((be->order() < 1) || (be->order() > 3)) {
      std::stringstream ss;
      ss << "bond order for bond with index " << bi << " is " << be->order()
         << "; it should be between 1 and 3";
      throw std::runtime_error(ss.str());
    }
    mol.getBondWithIdx(bi)->setBondType(bondType[be->order() - 1]);
  }
}

// init atom total valences and bond orders from the
// respective atoms and bonds
void ConjElectrons::initCeFromMol() {
  for (ConjAtomMap::const_iterator it = d_conjAtomMap.begin();
       it != d_conjAtomMap.end(); ++it) {
    it->second->initTvNbFcFromAtom();
  }
  for (ConjBondMap::const_iterator it = d_conjBondMap.begin();
       it != d_conjBondMap.end(); ++it) {
    it->second->initOrderFromBond();
  }
  d_currElectrons = 0;
}

// assign non-bonded electrons to atoms
void ConjElectrons::assignNonBonded() {
  for (ConjAtomMap::const_iterator it = d_conjAtomMap.begin();
       it != d_conjAtomMap.end(); ++it) {
    AtomElectrons *ae = it->second;
    unsigned int nb = std::min(ae->neededNbForOctet(), d_currElectrons);
    decrCurrElectrons(nb);
    ae->assignNonBonded(nb);
  }
}

// assign formal charges to atoms
void ConjElectrons::assignFormalCharge() {
  for (ConjAtomMap::const_iterator it = d_conjAtomMap.begin();
       it != d_conjAtomMap.end(); ++it) {
    it->second->assignFormalCharge();
  }
}

// return true if bond orders and formal charges for this ConjElectrons
// object are OK, false if they aren't
bool ConjElectrons::checkChargesAndBondOrders() {
  bool areAcceptable = true;
  bool haveIncompleteOctetRightOfC = false;
  bool haveNitrogenGroup = false;
  bool havePosLeftOfN = false;
  bool haveNegLeftOfN = false;
  for (ConjAtomMap::const_iterator it = d_conjAtomMap.begin();
       areAcceptable && (it != d_conjAtomMap.end()); ++it) {
    AtomElectrons *ae = it->second;
    // formal charges should be between -2 and +1
    areAcceptable = ((ae->fc() < 2) && (ae->fc() > -3));
    if (areAcceptable) {
      if (ae->fc() > 0) {
        d_flags |= HAVE_CATION;
      } else if (ae->fc() < 0) {
        d_flags |= HAVE_ANION;
      }
      if (ae->oe() > 4) {
        if (ae->oe() == 5) {
          haveNitrogenGroup = true;
        }
        if (!ae->hasOctet()) {
          haveIncompleteOctetRightOfC = true;
        }
        if ((ae->fc() > 0) && (ae->oe() > 5)) {
          d_flags |= HAVE_CATION_RIGHT_OF_N;
        }
      } else {
        if (ae->fc() < 0) {
          haveNegLeftOfN = true;
        } else if (ae->fc() > 0) {
          havePosLeftOfN = true;
        }
      }
      // no carbanions should coexist with incomplete octets
      areAcceptable = !(haveIncompleteOctetRightOfC && haveNegLeftOfN);
    }
  }
  if (areAcceptable && havePosLeftOfN &&
      !(d_parent->flags() & ResonanceMolSupplier::UNCONSTRAINED_CATIONS)) {
    // if the UNCONSTRAINED_CATIONS flag is not set, positively charged
    // atoms left of N with an incomplete octet are acceptable
    // only if the conjugated group has a positive total formal charge,
    // no atoms bearing a negative charge and no nitrogens
    areAcceptable = (d_totalFormalCharge > 0 && !(d_flags & HAVE_ANION) &&
                     !haveNitrogenGroup);
  }
  if (areAcceptable && haveNegLeftOfN &&
      !(d_parent->flags() & ResonanceMolSupplier::UNCONSTRAINED_ANIONS)) {
    // if the UNCONSTRAINED_ANIONS flag is not set, negatively charged
    // atoms left of N are acceptable only if the conjugated group has
    // a negative total formal charge and no atoms bearing a positive
    // charge
    areAcceptable = (d_totalFormalCharge < 0 && !(d_flags & HAVE_CATION));
  }
  // if the ALLOW_INCOMPLETE_OCTETS flag is not set, incomplete octets
  // right of C are only allowed in the absence of anions
  if (areAcceptable &&
      !(d_parent->flags() & ResonanceMolSupplier::ALLOW_INCOMPLETE_OCTETS)) {
    areAcceptable = !(haveIncompleteOctetRightOfC && (d_flags & HAVE_ANION));
  }
  for (ConjBondMap::const_iterator it = d_conjBondMap.begin();
       areAcceptable && (it != d_conjBondMap.end()); ++it) {
    BondElectrons *be = it->second;
    AtomElectrons *ae[2] = {d_conjAtomMap[be->bond()->getBeginAtomIdx()],
                            d_conjAtomMap[be->bond()->getEndAtomIdx()]};
    for (unsigned int i = 0; areAcceptable && (i < 2); ++i) {
      // cumulated bonds are disallowed in aromatic rings
      if (be->order() > 1 && ae[i]->atom()->getIsAromatic()) {
        if (ae[i]->hasMultipleBond()) {
          areAcceptable = false;
        } else {
          ae[i]->setHasMultipleBond();
        }
      }
      if (areAcceptable && ae[i]->oe() < 5) {
        // no carbocations allowed on carbons bearing double
        // or triple bonds, carbanions allowed on all carbons
        // no double charged allowed left of N
        areAcceptable =
            ((be->order() == 1) || ((be->order() > 1) && (ae[i]->fc() < 1)));
      }
    }
    if (areAcceptable) {
      // charged neighboring atoms left of N are not acceptable
      areAcceptable = !((ae[0]->oe() < 5) && ae[0]->fc() && (ae[1]->oe() < 5) &&
                        ae[1]->fc());
    }
  }
  return areAcceptable;
}

bool ConjElectrons::checkMetrics(CEStats &ceStats, bool &changed) const {
  bool ok = true;
  changed = false;
  if (!ceStats.d_init || (nbMissing() < ceStats.d_minNbMissing)) {
    changed = ceStats.d_init;
    ceStats.d_init = true;
    ceStats.d_minNbMissing = nbMissing();
  }
  if ((d_parent->flags() & ResonanceMolSupplier::ALLOW_INCOMPLETE_OCTETS) ||
      (nbMissing() <= ceStats.d_minNbMissing)) {
    if (!hasCationRightOfN() && !ceStats.d_haveNoCationsRightOfN) {
      ceStats.d_haveNoCationsRightOfN = true;
      changed = true;
    }
    if (!hasChargeSeparation() && !ceStats.d_haveNoChargeSeparation) {
      ceStats.d_haveNoChargeSeparation = true;
      changed = true;
    }
  }
  // if the flag ALLOW_INCOMPLETE_OCTETS is not set, ConjElectrons
  // objects having less electrons than the most electron-complete
  // structure (which most often will have all complete octets) are
  // discarded
  if ((!(d_parent->flags() & ResonanceMolSupplier::ALLOW_INCOMPLETE_OCTETS) &&
       (nbMissing() > ceStats.d_minNbMissing)) ||
      (!(d_parent->flags() & ResonanceMolSupplier::UNCONSTRAINED_CATIONS) &&
       hasCationRightOfN() && ceStats.d_haveNoCationsRightOfN) ||
      (!(d_parent->flags() & ResonanceMolSupplier::ALLOW_CHARGE_SEPARATION) &&
       hasChargeSeparation() && ceStats.d_haveNoChargeSeparation)) {
    ok = false;
  }
  return ok;
}

bool ConjElectrons::purgeMaps(CEMap &ceMap, CEDegCount &ceDegCount,
                              CEStats &ceStats) const {
  bool ok = true;
  bool changed = true;
  while (changed) {
    for (auto it = ceMap.begin(); it != ceMap.end();) {
      if (!it->second->checkMetrics(ceStats, changed)) {
        auto it2 = ceDegCount.find(it->second->hash());
        if (it2 != ceDegCount.end()) {
          --it2->second;
          if (!it2->second) {
            ceDegCount.erase(it2);
          }
        }
        if (it->second == this) {
          // postpone self deletion
          ok = false;
        } else {
          delete it->second;
        }
        it = ceMap.erase(it);
        changed = true;
      } else {
        ++it;
      }
      if (changed) {
        break;
      }
    }
  }
  return ok;
}

// assign formal charges and, if they are acceptable, store
// return true if FPs did not already exist in ceMap, false if they did
bool ConjElectrons::assignFormalChargesAndStore(CEMap &ceMap,
                                                CEDegCount &ceDegCount,
                                                CEStats &ceStats,
                                                unsigned int fpFlags) {
  assignFormalCharge();
  bool ok = checkChargesAndBondOrders();
  bool changed = false;
  if (ok) {
    computeMetrics();
    ok = checkMetrics(ceStats, changed);
  }
  if (ok) {
    ok = storeFP(ceMap, fpFlags);
  }
  if (changed) {
    ok = purgeMaps(ceMap, ceDegCount, ceStats) && ok;
  }
  if (ok) {
    updateDegCount(ceDegCount);
  }
  return ok;
}

// enumerate all possible permutations of non-bonded electrons
void ConjElectrons::enumerateNonBonded(CEMap &ceMap, CEDegCount &ceDegCount,
                                       CEStats &ceStats) {
  // the way we compute FPs for a resonance structure depends
  // on whether we want to enumerate all Kekule structures
  // or not; in the first case, we also include bond orders
  // in the FP computation in addition to atom valences
  const unsigned int fpFlags =
      FP_ATOMS |
      ((d_parent->flags() & ResonanceMolSupplier::KEKULE_ALL) ? FP_BONDS : 0);
  // count how many atoms need non-bonded electrons to complete
  // their octet ant store their indices in aiVec
  std::vector<unsigned int> aiVec;
  unsigned int nbTotal = 0;
  for (ConjAtomMap::const_iterator it = d_conjAtomMap.begin();
       it != d_conjAtomMap.end(); ++it) {
    unsigned int nb = it->second->neededNbForOctet();
    if (nb) {
      nbTotal += nb;
      aiVec.push_back(it->first);
    }
  }
  if (nbTotal > currElectrons()) {
    // if  the electrons required to satisfy all octets
    // are more than those currently available, some atoms will
    // be satisfied and some won't: we enumerate all permutations
    // of the unsatisfied atoms
    unsigned int missingElectrons = nbTotal - currElectrons();
    // number of atoms which won't have a complete octet
    unsigned int numCand = (missingElectrons - 1) / 2 + 1;
    unsigned int numComb;
    unsigned int v;
    // depending on the number of atoms which won't have a complete
    // octet and the number of those which need non-bonded electrons
    // we compute the number of permutations (numComb) and a
    // binary code (v) which indicates which of the atom indices in
    // aiVec will be octet-unsatisfied for each permutation
    ResonanceUtils::getNumCombStartV(numCand, aiVec.size(), numComb, v);
    // if there are multiple permutations, make a copy of the original
    // ConjElectrons object, since the latter will be modified
    ConjElectrons *ceCopy = new ConjElectrons(*this);
    // enumerate all permutations
    for (unsigned int c = 0; c < numComb; ++c) {
      ConjElectrons *ce = new ConjElectrons(*ceCopy);
      unsigned int vc = v;
      for (unsigned int i : aiVec) {
        AtomElectrons *ae = ce->getAtomElectronsWithIdx(i);
        unsigned int e = ae->neededNbForOctet();
        // if this atom was chosen to be octet-unsatisfied in
        // this permutation, give it one electron pair less than
        // its need (which most likely means no electrons at all)
        if (vc & 1) {
          e -= 2;
        }
        ce->decrCurrElectrons(e);
        ae->assignNonBonded(e);
        vc >>= 1;
      }
      // delete this candidate if it fails the formal charge check
      if (!ce->assignFormalChargesAndStore(ceMap, ceDegCount, ceStats,
                                           fpFlags)) {
        delete ce;
      }

      // get the next binary code
      ResonanceUtils::updateV(v);
    }
    delete ceCopy;
  } else if (nbTotal == currElectrons()) {
    ConjElectrons *ce = new ConjElectrons(*this);
    // if the electrons required to satisfy all octets
    // are as many as those currently available, assignment
    // is univocal
    ce->assignNonBonded();
    // delete this candidate if it fails the formal charge check
    if (!ce->assignFormalChargesAndStore(ceMap, ceDegCount, ceStats, fpFlags)) {
      delete ce;
    }
  }
  // if the electrons required to satisfy all octets are less
  // than those currently available, we must have failed the bond
  // assignment
}

void ConjElectrons::computeMetrics() {
  // 1000 * Electronegativity according to the Allen scale
  // (Allen, L.C. J. Am. Chem. Soc. 1989, 111, 9003-9014)
  static const std::vector<unsigned int> en{
      1000, 2300, 4160, 912,  1576, 2051, 2544, 3066, 3610, 4193, 4789, 869,
      1293, 1613, 1916, 2253, 2589, 2869, 3242, 734,  1034, 1190, 1380, 1530,
      1650, 1750, 1800, 1840, 1880, 1850, 1590, 1756, 1994, 2211, 2434, 2685,
      2966, 706,  963,  1120, 1320, 1410, 1470, 1510, 1540, 1560, 1590, 1870,
      1520, 1656, 1824, 1984, 2158, 2359, 2582, 659,  881,  1000, 1000, 1000,
      1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1090,
      1160, 1340, 1470, 1600, 1650, 1680, 1720, 1920, 1760, 1789, 1854, 2010,
      2190, 2390, 2600, 670,  890};
  for (ConjAtomMap::const_iterator it = d_conjAtomMap.begin();
       it != d_conjAtomMap.end(); ++it) {
    d_ceMetrics.d_absFormalCharges += abs(it->second->fc());
    size_t anIdx = it->second->atom()->getAtomicNum();
    d_ceMetrics.d_wtdFormalCharges +=
        (it->second->fc() * ((anIdx >= en.size()) ? 1000 : en[anIdx]));
    d_ceMetrics.d_nbMissing += it->second->neededNbForOctet();
  }
  computeDistFormalCharges();
  computeSumFormalChargeIdxs();
  computeSumMultipleBondIdxs();
}

void ConjElectrons::updateDegCount(CEDegCount &ceDegCount) {
  static const size_t ceMetricsArrayLen = 5;
  const auto metricsIt = reinterpret_cast<const int *>(&d_ceMetrics);
  d_ceMetrics.d_hash =
      boost::hash_range(metricsIt, metricsIt + ceMetricsArrayLen);
  ++ceDegCount[d_ceMetrics.d_hash];
}

// compute sum of shortest path distances between all pairs of
// formal charges of the same sign and of opposite signs
void ConjElectrons::computeDistFormalCharges() {
  unsigned int n = d_parent->mol().getNumAtoms();
  for (ConjAtomMap::const_iterator it1 = d_conjAtomMap.begin();
       it1 != d_conjAtomMap.end(); ++it1) {
    if (!it1->second->fc()) {
      continue;
    }
    size_t y = it1->first * n;
    for (auto it2 = it1; it2 != d_conjAtomMap.end(); ++it2) {
      if ((it1 == it2) || !it2->second->fc()) {
        continue;
      }
      unsigned int dist = static_cast<unsigned int>(
          MolOps::getDistanceMat(d_parent->mol())[y + it2->first] + 0.1);
      if ((it1->second->fc() * it2->second->fc()) > 0) {
        d_ceMetrics.d_fcSameSignDist += dist;
      } else {
        d_ceMetrics.d_fcOppSignDist += dist;
      }
    }
  }
}

// compute the sum of indices of atoms bearing a formal charge
void ConjElectrons::computeSumFormalChargeIdxs() {
  for (ConjAtomMap::const_iterator it = d_conjAtomMap.begin();
       it != d_conjAtomMap.end(); ++it) {
    if (it->second->fc()) {
      d_ceMetrics.d_sumFormalChargeIdxs += it->first;
    }
  }
}

// compute the sum of indices of multiple bonds
void ConjElectrons::computeSumMultipleBondIdxs() {
  for (ConjBondMap::const_iterator it = d_conjBondMap.begin();
       it != d_conjBondMap.end(); ++it) {
    if (it->second->order() > 1) {
      d_ceMetrics.d_sumMultipleBondIdxs += it->first;
    }
  }
}

// decrement the count of current electrons by d
void ConjElectrons::decrCurrElectrons(unsigned int d) {
  if (d_currElectrons < d) {
    std::stringstream ss;
    ss << "d_currElectrons = " << d_currElectrons << ", d = " << d;
    throw std::runtime_error(ss.str());
  }
  d_currElectrons -= d;
}

// push atom index ai to the begin stack if it is not already there
void ConjElectrons::pushToBeginStack(unsigned int ai) {
  if (!d_conjAtomMap[ai]->isStacked()) {
    d_conjAtomMap[ai]->setStacked();
    d_beginAIStack.push(ai);
  }
}

// pop an atom from the begin stack and put its index in ai
bool ConjElectrons::popFromBeginStack(unsigned int &ai) {
  bool ok = false;
  while (!d_beginAIStack.empty() && !ok) {
    ai = d_beginAIStack.top();
    d_beginAIStack.pop();
    AtomElectrons *ae = d_conjAtomMap[ai];
    ae->clearStacked();
    ok = (!ae->isDefinitive());
  }
  return ok;
}

// function used to sort ConjElectrons objects based on their
// importance to describe the structure; criteria in order of decreasing
// priority follow:
// 1) Number of unsatisfied octets
// 2) Number of formal charges
// 3) Number of formal charges weighted by atom electronegativity
// 4) Distance between formal charges with the same sign
// 5) Distance between formal charges with opposite signs
// 6) Sum of the indices of atoms bearing a formal charge
// 7) Sum of the indices of multiple bonds
bool CEVect2::resonanceStructureCompare(const ConjElectrons *a,
                                        const ConjElectrons *b) {
  return (
      (a->nbMissing() != b->nbMissing())
          ? (a->nbMissing() < b->nbMissing())
          : (a->absFormalCharges() != b->absFormalCharges())
                ? (a->absFormalCharges() < b->absFormalCharges())
                : (a->wtdFormalCharges() != b->wtdFormalCharges())
                      ? (a->wtdFormalCharges() < b->wtdFormalCharges())
                      : (a->fcSameSignDist() != b->fcSameSignDist())
                            ? (a->fcSameSignDist() > b->fcSameSignDist())
                            : (a->fcOppSignDist() != b->fcOppSignDist())
                                  ? (a->fcOppSignDist() > b->fcOppSignDist())
                                  : (a->sumFormalChargeIdxs() !=
                                     b->sumFormalChargeIdxs())
                                        ? (a->sumFormalChargeIdxs() <
                                           b->sumFormalChargeIdxs())
                                        : (a->sumMultipleBondIdxs() <
                                           b->sumMultipleBondIdxs()));
}

CEVect2::CEVect2(const CEMap &ceMap) {
  d_ceVect.reserve(ceMap.size());
  for (CEMap::const_iterator it = ceMap.begin(); it != ceMap.end(); ++it) {
    d_ceVect.push_back(it->second);
  }
  std::sort(d_ceVect.begin(), d_ceVect.end(), resonanceStructureCompare);
  bool first = true;
  std::size_t hashPrev;
  for (CEVect::const_iterator it = d_ceVect.begin(); it != d_ceVect.end();
       ++it) {
    if (first || ((*it)->hash() != hashPrev)) {
      first = false;
      hashPrev = (*it)->hash();
      d_degVect.push_back(1);
    } else {
      ++d_degVect.back();
    }
  }
}

ConjElectrons *CEVect2::getCE(unsigned int depth, unsigned int width) {
  if (depth >= d_degVect.size()) {
    std::stringstream ss;
    ss << "depth = " << depth << ", d_degVect.size() = " << d_degVect.size();
    throw std::runtime_error(ss.str());
  }
  if (width >= d_degVect[depth]) {
    std::stringstream ss;
    ss << "width = " << width << ", d_degVect[" << depth
       << "] = " << d_degVect[depth];
    throw std::runtime_error(ss.str());
  }
  size_t i =
      std::accumulate(d_degVect.begin(), d_degVect.begin() + depth, 0U) + width;
  return d_ceVect[i];
}

void CEVect2::resize(size_t size) {
  d_ceVect.resize(size ? ceCountUntilDepth(size - 1) : 0);
  d_degVect.resize(size);
}

unsigned int CEVect2::ceCountAtDepth(unsigned int depth) {
  if (depth >= d_degVect.size()) {
    std::stringstream ss;
    ss << "depth = " << depth << ", d_degVect.size() = " << d_degVect.size();
    throw std::runtime_error(ss.str());
  }
  return d_degVect[depth];
}

unsigned int CEVect2::ceCountUntilDepth(unsigned int depth) {
  if (depth >= d_degVect.size()) {
    std::stringstream ss;
    ss << "depth = " << depth << ", d_degVect.size() = " << d_degVect.size();
    throw std::runtime_error(ss.str());
  }
  return std::accumulate(d_degVect.begin(), d_degVect.begin() + depth + 1, 0U);
}

void CEVect2::idxToDepthWidth(size_t idx, unsigned int &d, unsigned int &w) {
  if (idx >= d_ceVect.size()) {
    std::stringstream ss;
    ss << "idx = " << idx << ", d_ceVect.size() = " << d_ceVect.size();
    throw std::runtime_error(ss.str());
  }
  d = 0;
  while (idx >= d_degVect[d]) {
    idx -= d_degVect[d];
    ++d;
  }
  w = idx;
}

// get the pointer to the BondElectrons object for bond having index bi
BondElectrons *ConjElectrons::getBondElectronsWithIdx(unsigned int bi) {
  auto it = d_conjBondMap.find(bi);
  return (it != d_conjBondMap.end() ? it->second : nullptr);
}

// get the pointer to the AtomElectrons object for atom having index ai
AtomElectrons *ConjElectrons::getAtomElectronsWithIdx(unsigned int ai) {
  auto it = d_conjAtomMap.find(ai);
  return (it != d_conjAtomMap.end() ? it->second : nullptr);
}

// count number of total electrons
unsigned int ConjElectrons::countTotalElectrons() {
  // count total number of valence electrons in conjugated group
  for (ConjBondMap::const_iterator it = d_conjBondMap.begin();
       it != d_conjBondMap.end(); ++it) {
    d_totalElectrons += (2 * it->second->orderFromBond());
  }
  for (ConjAtomMap::const_iterator it = d_conjAtomMap.begin();
       it != d_conjAtomMap.end(); ++it) {
    const Atom *a = it->second->atom();
    d_totalElectrons +=
        it->second->oe() - a->getTotalValence() - a->getFormalCharge();
  }
  return d_totalElectrons;
}

size_t ResonanceMolSupplierCallback::getNumStructures(
    unsigned int conjGrpIdx) const {
  PRECONDITION(conjGrpIdx < d_progress.size(), "conjGrpIdx out of bounds");
  return d_progress.at(conjGrpIdx).d_totalStructs;
}

size_t ResonanceMolSupplierCallback::getNumDiverseStructures(
    unsigned int conjGrpIdx) const {
  PRECONDITION(conjGrpIdx < d_progress.size(), "conjGrpIdx out of bounds");
  return d_progress.at(conjGrpIdx).d_diverseStructs;
}

void ResonanceMolSupplier::setProgressCallback(
    ResonanceMolSupplierCallback *callback) {
  d_callback.reset(callback);
}

// if the total number of resonance structures exceeds d_maxStructs,
// trim the number of ConjElectrons objects for each conjugated group
// to exclude the less likely structures and save memory
// we are going to generate the complete resonance structures
// combining the most likely ConjElectrons objects in a breadth-first
// fashion
void ResonanceMolSupplier::trimCeVect2() {
  if (d_length == d_maxStructs) {
    std::vector<unsigned int> s(d_nConjGrp, 0);
    std::vector<unsigned int> t(d_nConjGrp, 0);
    unsigned int currSize = 0;
    while (currSize < d_length) {
      currSize = 1;
      for (unsigned int conjGrpIdx = 0;
           (currSize < d_length) && (conjGrpIdx < d_nConjGrp); ++conjGrpIdx) {
        if (s[conjGrpIdx] < d_ceVect3[conjGrpIdx]->depth()) {
          t[conjGrpIdx] += d_ceVect3[conjGrpIdx]->ceCountAtDepth(s[conjGrpIdx]);
          ++s[conjGrpIdx];
        }
        currSize *= t[conjGrpIdx];
      }
    }
    for (unsigned int conjGrpIdx = 0; conjGrpIdx < d_nConjGrp; ++conjGrpIdx) {
      for (unsigned int d = s[conjGrpIdx]; d < d_ceVect3[conjGrpIdx]->depth();
           ++d) {
        for (unsigned int w = 0; w < d_ceVect3[conjGrpIdx]->ceCountAtDepth(d);
             ++w) {
          delete d_ceVect3[conjGrpIdx]->getCE(d, w);
        }
      }
      d_ceVect3[conjGrpIdx]->resize(s[conjGrpIdx]);
    }
  }
}

// get the ConjElectrons indices to be combined given
// the complete resonance structure index
void ResonanceMolSupplier::idxToCEPerm(size_t idx,
                                       std::vector<unsigned int> &c) const {
  // the c vector holds a pair of values for each ConjGrp,
  // namely depth and width where the ConjElectrons
  // object lies in the CEVect2 vector
  c.resize(d_nConjGrp * 2);
  unsigned int g = d_nConjGrp;
  while (g) {
    --g;
    unsigned int gt2 = g * 2;
    unsigned int d = 1;
    for (unsigned int j = 0; j < g; ++j) {
      d *= d_ceVect3[j]->ceCount();
    }
    if (d_ceVect3[g]->ceCount()) {
      d_ceVect3[g]->idxToDepthWidth(idx / d, c[gt2], c[gt2 + 1]);
      idx %= d;
    }
  }
}

// sort the vectors of ConjElectrons indices in such a way that we
// generate resonance structures out of the most likely ConjElectrons
// objects in a breadth-first fashion
bool ResonanceMolSupplier::cePermCompare(const CEPerm *a, const CEPerm *b) {
  unsigned int aSum = 0;
  unsigned int bSum = 0;
  for (size_t i = 0; i < a->v.size(); i += 2) {
    aSum += a->v[i];
    bSum += b->v[i];
  }
  if (aSum != bSum) {
    return (aSum < bSum);
  }
  unsigned int aMax = 0;
  unsigned int bMax = 0;
  for (size_t i = 0; i < a->v.size(); i += 2) {
    if (!i || (a->v[i] > aMax)) {
      aMax = a->v[i];
    }
    if (!i || (b->v[i] > bMax)) {
      bMax = b->v[i];
    }
  }
  if (aMax != bMax) {
    return (aMax < bMax);
  }
  for (size_t i = 0; i < a->v.size(); i += 2) {
    if (a->v[i] != b->v[i]) {
      return (a->v[i] < b->v[i]);
    }
  }
  // if the other criteria didn't discriminate,
  // sort based on degenerate resonance structures
  for (size_t i = 1; i < a->v.size(); i += 2) {
    if (a->v[i] != b->v[i]) {
      return (a->v[i] < b->v[i]);
    }
  }
  // we'll never get here, this it is just to silence a warning
  return false;
}

// enumerate upfront all the d_length indices for complete resonance
// structures and sort them to privilege the most likely
// ConjElectrons index combinations in a breadth-first fashion
void ResonanceMolSupplier::prepEnumIdxVect() {
  d_enumIdx.resize(d_length);
  std::vector<CEPerm *> cePermVect(d_length);
  for (size_t i = 0; i < d_length; ++i) {
    cePermVect[i] = new CEPerm;
    cePermVect[i]->idx = i;
    idxToCEPerm(i, cePermVect[i]->v);
  }
  std::sort(cePermVect.begin(), cePermVect.end(), cePermCompare);
  for (unsigned int i = 0; i < d_length; ++i) {
    d_enumIdx[i] = cePermVect[i]->idx;
    delete cePermVect[i];
  }
}

// object constructor
ResonanceMolSupplier::ResonanceMolSupplier(ROMol &mol, unsigned int flags,
                                           unsigned int maxStructs)
    : d_nConjGrp(0),
      d_flags(flags),
      d_idx(0),
      d_numThreads(1),
      d_isEnumerated(false),
      d_wasCanceled(false) {
  const unsigned int MAX_STRUCTS = 1000000;
  if (d_flags & UNCONSTRAINED_CATIONS) {
    d_flags |= (ALLOW_CHARGE_SEPARATION | ALLOW_INCOMPLETE_OCTETS);
  }
  if (d_flags & UNCONSTRAINED_ANIONS) {
    d_flags |= ALLOW_CHARGE_SEPARATION;
  }
  d_maxStructs = std::min(maxStructs, MAX_STRUCTS);
  d_length = std::min(1U, d_maxStructs);
  d_mol = new ROMol(mol);
  // call here to avoid a race condition in threads
  MolOps::getDistanceMat(*d_mol);
  MolOps::Kekulize((RWMol &)*d_mol, false);
  // identify conjugate substructures
  assignConjGrpIdx();
}

// object destructor
ResonanceMolSupplier::~ResonanceMolSupplier() {
  for (CEVect3::const_iterator ceVect3It = d_ceVect3.begin();
       ceVect3It != d_ceVect3.end(); ++ceVect3It) {
    if (!(*ceVect3It)) {
      continue;
    }
    for (unsigned int d = 0; d < (*ceVect3It)->depth(); ++d) {
      for (unsigned int w = 0; w < (*ceVect3It)->ceCountAtDepth(d); ++w) {
        delete (*ceVect3It)->getCE(d, w);
      }
    }
    delete *ceVect3It;
  }
  if (d_mol) {
    delete d_mol;
  }
}

void ResonanceMolSupplier::setNumThreads(int numThreads) {
  d_numThreads = std::min(d_nConjGrp, getNumThreadsToUse(numThreads));
}

void ResonanceMolSupplier::enumerate() {
  if (d_isEnumerated) {
    return;
  }
  if (d_callback.get()) {
    d_callback->d_progress.resize(d_nConjGrp);
  }
  resizeCeVect();
  if (d_numThreads == 1) {
    mainLoop(0, 1);
  }
#ifdef RDK_BUILD_THREADSAFE_SSS
  else {
    std::vector<std::future<void>> tg;
    auto functor = [this](unsigned int ti, unsigned int d_numThreads) -> void {
      mainLoop(ti, d_numThreads);
    };
    for (unsigned int ti = 0; ti < d_numThreads; ++ti) {
      tg.emplace_back(
          std::async(std::launch::async, functor, ti, d_numThreads));
    }
    for (auto &fut : tg) {
      fut.get();
    }
  }
#endif
  setResonanceMolSupplierLength();
  trimCeVect2();
  prepEnumIdxVect();
  d_isEnumerated = true;
}

void ResonanceMolSupplier::mainLoop(unsigned int ti, unsigned int nt) {
  for (unsigned int conjGrpIdx = 0; conjGrpIdx < d_nConjGrp; ++conjGrpIdx) {
    if ((conjGrpIdx % nt) != ti) {
      continue;
    }
    CEMap ceMap;
    buildCEMap(ceMap, conjGrpIdx);
    storeCEMap(ceMap, conjGrpIdx);
  }
}

// each bond an atom is assigned an index representing the conjugated
// group it belongs to; such indices are stored in two vectors
// (d_bondConjGrpIdx and d_atomConjGrpIdx, respectively)
// atoms and bonds which do not belong to a conjugated group are given
// index -1
void ResonanceMolSupplier::assignConjGrpIdx() {
  unsigned int nb = d_mol->getNumBonds();
  if (!nb) {
    return;
  }
  std::stack<unsigned int> bondIndexStack;
  d_bondConjGrpIdx.resize(nb, -1);
  unsigned int na = d_mol->getNumAtoms();
  d_atomConjGrpIdx.resize(na, -1);
  for (d_nConjGrp = 0; true; ++d_nConjGrp) {
    for (const auto b : d_mol->bonds()) {
      unsigned int i = b->getIdx();
      if (b->getIsConjugated() && (d_bondConjGrpIdx[i] == -1)) {
        bondIndexStack.push(i);
        break;
      }
    }
    if (bondIndexStack.empty()) {
      break;
    }
    while (!bondIndexStack.empty()) {
      unsigned int i = bondIndexStack.top();
      bondIndexStack.pop();
      const Bond *bi = d_mol->getBondWithIdx(i);
      for (const Atom *bondAtom : {bi->getBeginAtom(), bi->getEndAtom()}) {
        // loop over bonds sprouting from bondAtom
        for (const auto &bNbri :
             boost::make_iterator_range(d_mol->getAtomBonds(bondAtom))) {
          const auto &bNbr = (*d_mol)[bNbri];
          unsigned int biNbr = bNbr->getIdx();
          if (bNbr->getIsConjugated() && (d_bondConjGrpIdx[biNbr] == -1)) {
            d_bondConjGrpIdx[biNbr] = d_nConjGrp;
            bondIndexStack.push(biNbr);
          }
        }
      }
    }
  }
  for (unsigned int i = 0; i < nb; ++i) {
    if (d_bondConjGrpIdx[i] != -1) {
      const Bond *bi = d_mol->getBondWithIdx(i);
      unsigned int biBeginIdx = bi->getBeginAtomIdx();
      unsigned int biEndIdx = bi->getEndAtomIdx();
      d_atomConjGrpIdx[biBeginIdx] = d_bondConjGrpIdx[i];
      d_atomConjGrpIdx[biEndIdx] = d_bondConjGrpIdx[i];
    }
  }
}

// function which enumerates all possible multiple bond arrangements
// for each conjugated group, stores each arrangement in a ConjElectrons
// object and stores the latter in a map (ceMapTmp), keyed with its
// fingerprints. Depending on whether the KEKULE_ALL flag is set, the FP
// computation will be based either on the bond arrangement or on atom
// valences; this provides a convenient mechanism to remove degenerate
// Kekule structures upfront. In a subsequent step, ConjElectrons objects
// stored in ceMapTmp are further enumerated for their non-bonded
// electron arrangements, and the final ConjElectrons objects are stored
// in ceMap. Depending on whether the KEKULE_ALL flag is set, the FP
// computation will involve or not bond arrangement in addition to atom
// valences
void ResonanceMolSupplier::buildCEMap(CEMap &ceMap, unsigned int conjGrpIdx) {
  const unsigned int BEGIN_POS = 0;
  const unsigned int END_POS = 1;
  const unsigned int fpFlags =
      ConjElectrons::FP_ATOMS |
      ((d_flags & KEKULE_ALL) ? ConjElectrons::FP_BONDS : 0);
  const unsigned int fpFlagsTmp =
      ((d_flags & KEKULE_ALL) ? ConjElectrons::FP_BONDS
                              : ConjElectrons::FP_ATOMS);
  unsigned int nb = d_mol->getNumBonds();
  unsigned int na = d_mol->getNumAtoms();
  CESet ceSetTmp;
  CEDegCount ceDegCount;
  CEStats ceStats;
  // There are two stacks in this algorithm:
  // 1) ConjElectrons stack (ceStack, unique). It stores partial bond
  //    arrangements whenever multiple paths arise. Partial bond
  //    arrangements are then popped from the stack until the latter is
  //    empty
  // 2) atom index stack (d_beginAIStack, associated to each
  //    ConjElectrons object). It stores a stack of atom indices to
  //    start from to complete the bond arrangement for each
  //    ConjElectrons object
  std::stack<ConjElectrons *> ceStack;
  auto *ce = new ConjElectrons(this, conjGrpIdx);
  auto *ceCopy = new ConjElectrons(*ce);
  // the first ConjElectrons object has the user-supplied bond
  // and formal charge arrangement and is stored as such
  ce->initCeFromMol();
  // we ignore the result of the call to checkChargesAndBondOrders()
  // but we need to call it so that the HAVE_CATION_RIGHT_OF_N flag
  // may eventually be set
  ce->checkChargesAndBondOrders();
  ce->computeMetrics();
  ce->updateDegCount(ceDegCount);
  ce->storeFP(ceMap, fpFlags);
  ce = ceCopy;
  // initialize ceStack
  ceStack.push(ce);
  // loop until ceStack is empty
  while (!ceStack.empty()) {
    ce = ceStack.top();
    ceStack.pop();
    unsigned int aiBegin = 0;
    // if the atom index stack is empty, initialize it with a primer;
    // any atom index belonging to this conjugated group will do
    if (ce->isBeginStackEmpty()) {
      bool aiFound = false;
      while (aiBegin < na) {
        aiFound = (d_atomConjGrpIdx[aiBegin] == static_cast<int>(conjGrpIdx));
        if (aiFound) {
          break;
        }
        ++aiBegin;
      }
      if (!aiFound) {
        continue;
      }
      ce->pushToBeginStack(aiBegin);
    }
    // loop until the atom index stack is empty
    while (ce && ce->popFromBeginStack(aiBegin)) {
      // aiBegin holds the atom index just popped from stack
      unsigned int ai[2] = {aiBegin, na};
      AtomElectrons *ae[2] = {ce->getAtomElectronsWithIdx(ai[BEGIN_POS]),
                              nullptr};
      unsigned int bi = nb;
      BondElectrons *be = nullptr;
      // loop over neighbors of the atom popped from the
      // atom index stack
      ROMol::ADJ_ITER nbrIdx, endNbrs;
      boost::tie(nbrIdx, endNbrs) =
          d_mol->getAtomNeighbors(ae[BEGIN_POS]->atom());
      for (; nbrIdx != endNbrs; ++nbrIdx) {
        unsigned int aiNbr = *nbrIdx;
        // if this neighbor is not part of the conjugated group,
        // ignore it
        if (ce->parent()->getAtomConjGrpIdx(aiNbr) !=
            static_cast<int>(conjGrpIdx)) {
          continue;
        }
        AtomElectrons *aeNbr = ce->getAtomElectronsWithIdx(aiNbr);
        // if we've already dealt with this neighbor before, ignore it
        if (!aeNbr || aeNbr->isDefinitive()) {
          continue;
        }
        unsigned int biNbr =
            d_mol->getBondBetweenAtoms(ai[BEGIN_POS], aiNbr)->getIdx();
        BondElectrons *beNbr = ce->getBondElectronsWithIdx(biNbr);
        // if we have already assigned the bond order to this bond,
        // ignore it
        if (!beNbr || beNbr->isDefinitive()) {
          continue;
        }
        // if this is the first neighbor we find, process it
        if (!ae[END_POS]) {
          bi = biNbr;
          be = beNbr;
          ai[END_POS] = aiNbr;
          ae[END_POS] = aeNbr;
        }
        // otherwise, if there are multiple neighbors, process the first
        // and store the rest to the atom index stack
        else {
          ce->pushToBeginStack(aiNbr);
        }
      }
      // if no neighbors were found, move on to the next atom index
      // in the stack
      if (!be) {
        continue;
      }
      std::uint8_t allowedBondsPerAtom[2];
      // loop over the two atoms that need be bonded
      for (unsigned int i = 0; i < 2; ++i) {
        // for each atom, find which bond orders are allowed
        allowedBondsPerAtom[i] = ae[i]->findAllowedBonds(bi);
        // if for either atom this is the last bond, mark the
        // atom as definitive, otherwise push the atom index
        // to the stack
        if (ae[i]->isLastBond()) {
          ae[i]->setDefinitive();
        } else {
          ce->pushToBeginStack(ai[i]);
        }
      }
      // logical AND between allowed bond masks for the two atoms
      std::uint8_t allowedBonds =
          allowedBondsPerAtom[BEGIN_POS] & allowedBondsPerAtom[END_POS];
      bool isAnyBondAllowed = false;
      bool needToFork = false;
      // make a copy of the current ConjElectrons object
      ceCopy = new ConjElectrons(*ce);
      // consider single, double and triple bond alternatives in turn
      for (unsigned int i = 0; i < 3; ++i) {
        unsigned int t = i * 2;
        ConjElectrons *ceToSet = nullptr;
        std::uint8_t orderMask = (1 << t);
        std::uint8_t chgMask = (1 << (t + 1));
        unsigned int bo = i + 1;
        // if the currently available electrons are enough for this
        // bond type and both atoms can accept this bond type
        // and we are not in a situation where both atoms are left
        // of N and both need a charge to accept this bond type,
        // then this bond type is feasible
        if ((ce->currElectrons() >= (bo * 2)) && (allowedBonds & orderMask) &&
            !((allowedBonds & chgMask) &&
              ((ae[BEGIN_POS]->oe() < 5) || (ae[END_POS]->oe() < 5)))) {
          isAnyBondAllowed = true;
          // if another bond type will be possible, then we'll
          // need to save that to ceStack and keep on with the current
          // ConjElectrons object
          if (!needToFork) {
            needToFork = true;
            ceToSet = ce;
          } else {
            auto *ceFork = new ConjElectrons(*ceCopy);
            ceStack.push(ceFork);
            ceToSet = ceFork;
          }
        }
        if (ceToSet) {
          // set bond orders and finalize atoms as needed
          ceToSet->getBondElectronsWithIdx(bi)->setOrder(bo);
          for (unsigned int j = 0; j < 2; ++j) {
            if (ae[j]->isLastBond()) {
              ceToSet->getAtomElectronsWithIdx(ai[j])->finalizeAtom();
            }
          }
        }
      }
      delete ceCopy;
      // if a dead end was hit, discard this ConjElectrons object
      if (!isAnyBondAllowed) {
        delete ce;
        ce = nullptr;
      }
    }
    if (ce) {
      // if this bond arrangement was already stored previously,
      // discard it
      if (!ce->haveFP(ceSetTmp, fpFlagsTmp)) {
        delete ce;
      } else {
        // enumerate possible non-bonded electron
        // arrangements, and store them in ceMap
        ce->enumerateNonBonded(ceMap, ceDegCount, ceStats);
        delete ce;
        // quit the loop early if the number of non-degenerate
        // structures already exceeds d_maxStructs
        if (d_callback.get()) {
          d_callback->d_progress[conjGrpIdx].d_totalStructs = ceMap.size();
          d_callback->d_progress[conjGrpIdx].d_diverseStructs =
              ceDegCount.size();
          if (!(*d_callback)()) {
            d_wasCanceled = true;
            break;
          }
        }
        if (ceDegCount.size() >= d_maxStructs) {
          break;
        }
      }
    }
  }
  while (!ceStack.empty()) {
    delete ceStack.top();
    ceStack.pop();
  }
}

// getter function which returns the bondConjGrpIdx for a given
// bond index, or -1 if the bond is not conjugated
int ResonanceMolSupplier::getBondConjGrpIdx(unsigned int bi) const {
  if (bi >= d_bondConjGrpIdx.size()) {
    std::stringstream ss;
    ss << "d_bondConjGrpIdx.size() = " << d_bondConjGrpIdx.size()
       << ", bi = " << bi;
    throw std::runtime_error(ss.str());
  }
  return d_bondConjGrpIdx[bi];
}

// getter function which returns the atomConjGrpIdx for a given
// atom index, or -1 if the atom is not conjugated
int ResonanceMolSupplier::getAtomConjGrpIdx(unsigned int ai) const {
  if (ai >= d_atomConjGrpIdx.size()) {
    std::stringstream ss;
    ss << "d_atomConjGrpIdx.size() = " << d_atomConjGrpIdx.size()
       << ", ai = " << ai;
    throw std::runtime_error(ss.str());
  }
  return d_atomConjGrpIdx[ai];
}

// resizes d_ceVect3 vector
inline void ResonanceMolSupplier::resizeCeVect() {
  d_ceVect3.resize(d_nConjGrp, nullptr);
}

// stores the ConjElectrons pointers currently stored in ceMap
// in the d_ceVect3 vector
inline void ResonanceMolSupplier::storeCEMap(const CEMap &ceMap,
                                             unsigned int conjGrpIdx) {
  d_ceVect3[conjGrpIdx] = new CEVect2(ceMap);
}

void ResonanceMolSupplier::setResonanceMolSupplierLength() {
  for (unsigned int i = 0; (d_length < d_maxStructs) && (i < d_ceVect3.size());
       ++i) {
    boost::uint64_t p = d_length * d_ceVect3[i]->ceCount();
    d_length =
        ((p < d_maxStructs) ? static_cast<unsigned int>(p) : d_maxStructs);
  }
}

// Returns the number of resonance structures in the
// ResonanceMolSupplier
unsigned int ResonanceMolSupplier::length() {
  enumerate();
  return d_length;
}

// Resets the ResonanceMolSupplier index
void ResonanceMolSupplier::reset() {
  enumerate();
  d_idx = 0;
}

// Returns true if there are no more resonance structures left
bool ResonanceMolSupplier::atEnd() {
  enumerate();
  return (d_idx == d_length);
}

// Returns a pointer to the next resonance structure as a ROMol,
// or NULL if there are no more resonance structures left.
// The caller is responsible for freeing memory associated to
// the pointer
ROMol *ResonanceMolSupplier::next() {
  enumerate();
  return (atEnd() ? nullptr : (*this)[d_idx++]);
}

// sets the ResonanceMolSupplier index to idx
void ResonanceMolSupplier::moveTo(unsigned int idx) {
  enumerate();
  if (idx >= d_length) {
    std::stringstream ss;
    ss << "d_length = " << d_length << ", idx = " << idx;
    throw std::runtime_error(ss.str());
  }
  d_idx = idx;
}

// returns the resonance structure with index idx as a ROMol *
// the index returns resonance structures combining ConjElectrons
// objects in a breadth-first fashion, in order to return the most
// likely complete resonance structures first
ROMol *ResonanceMolSupplier::operator[](unsigned int idx) {
  enumerate();
  if (idx >= d_length) {
    std::stringstream ss;
    ss << "d_length = " << d_length << ", idx = " << idx;
    throw std::runtime_error(ss.str());
  }
  std::vector<unsigned int> c;
  idxToCEPerm(d_enumIdx[idx], c);
  return assignBondsFormalCharges(c);
}

// helper function to assign bond orders and formal charges to the
// mol passed as reference out of the ConjElectrons objects whose
// indices are passed with the c vector
void ResonanceMolSupplier::assignBondsFormalChargesHelper(
    ROMol &mol, std::vector<unsigned int> &c) const {
  for (unsigned int conjGrpIdx = 0; conjGrpIdx < d_nConjGrp; ++conjGrpIdx) {
    unsigned int i = conjGrpIdx * 2;
    ConjElectrons *ce = d_ceVect3[conjGrpIdx]->getCE(c[i], c[i + 1]);
    ce->assignBondsFormalChargesToMol(mol);
  }
}

// returns a pointer to a ROMol with bond orders and formal charges
// assigned out of the ConjElectrons objects whose indices are passed
// with the c vector
ROMol *ResonanceMolSupplier::assignBondsFormalCharges(
    std::vector<unsigned int> &c) const {
  auto *mol = new ROMol(this->mol());
  assignBondsFormalChargesHelper(*mol, c);
  ResonanceUtils::fixExplicitImplicitHs(*mol);
  ResonanceUtils::sanitizeMol((RWMol &)*mol);
  return mol;
}
}  // namespace RDKit
