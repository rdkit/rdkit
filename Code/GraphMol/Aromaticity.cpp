//
//  Copyright (C) 2003-2017 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/Rings.h>
#include <RDGeneral/types.h>
#include <boost/dynamic_bitset.hpp>
#include <set>

// introduced for the sake of efficiency
// this is the maximum ring size that will be considered
// as a candidate for fused-ring aromaticity. This is picked to
// be a bit bigger than the outer ring in a porphyrin
// This came up while fixing sf.net issue249
const unsigned int maxFusedAromaticRingSize = 24;

/****************************************************************
Here are some molecules that have troubled us aromaticity wise

1. We get 2 aromatic rings for this molecule, but daylight says none. Also
   if we replace [O+] with the N daylight says we have two aromatic rings
   [O-][N+]1=CC2=CC=C[O+]2[Cu]13[O+]4C=CC=C4C=[N+]3[O-]

2. We the fused ring system with all the rings in it is considered aroamtic by
our code.
   This is because we count electrons from buried atoms as well when
   we are dealing with fused rings

3. Here's a weird fused ring case, Pattern NAN, A overall
   O=C3C2=CC1=CC=COC1=CC2=CC=C3
 ************************************************************/
namespace RingUtils {
using namespace RDKit;

void pickFusedRings(int curr, const INT_INT_VECT_MAP &neighMap, INT_VECT &res,
                    boost::dynamic_bitset<> &done, int depth) {
  auto pos = neighMap.find(curr);
  PRECONDITION(pos != neighMap.end(), "bad argument");
  done[curr] = 1;
  res.push_back(curr);

  const INT_VECT &neighs = pos->second;
#if 0
    std::cerr<<"depth: "<<depth<<" ring: "<<curr<<" size: "<<res.size()<<" neighs: "<<neighs.size()<<std::endl;
    std::cerr<<"   ";
    std::copy(neighs.begin(),neighs.end(),std::ostream_iterator<int>(std::cerr," "));
    std::cerr<<"\n";
#endif
  for (int neigh : neighs) {
    if (!done[neigh]) {
      pickFusedRings(neigh, neighMap, res, done, depth + 1);
    }
  }
}

bool checkFused(const INT_VECT &rids, INT_INT_VECT_MAP &ringNeighs) {
  INT_INT_VECT_MAP_CI nci;
  int nrings = rdcast<int>(ringNeighs.size());
  boost::dynamic_bitset<> done(nrings);
  int rid;
  INT_VECT fused;

  // mark all rings in the system other than those in rids as done
  for (nci = ringNeighs.begin(); nci != ringNeighs.end(); nci++) {
    rid = (*nci).first;
    if (std::find(rids.begin(), rids.end(), rid) == rids.end()) {
      done[rid] = 1;
    }
  }

  // then pick a fused system from the remaining (i.e. rids)
  // If the rings in rids are fused we should get back all of them
  // in fused
  // if we get a smaller number in fused then rids are not fused
  pickFusedRings(rids.front(), ringNeighs, fused, done);

  CHECK_INVARIANT(fused.size() <= rids.size(), "");
  return (fused.size() == rids.size());
}

void makeRingNeighborMap(const VECT_INT_VECT &brings,
                         INT_INT_VECT_MAP &neighMap, unsigned int maxSize,
                         unsigned int maxOverlapSize) {
  int nrings = rdcast<int>(brings.size());
  int i, j;
  INT_VECT ring1;
  for (i = 0; i < nrings; i++) {
    INT_VECT neighs;
    neighMap[i] = neighs;
  }

  for (i = 0; i < nrings; i++) {
    if (maxSize && brings[i].size() > maxSize) {
      continue;
    }
    ring1 = brings[i];
    for (j = i + 1; j < nrings; j++) {
      if (maxSize && brings[j].size() > maxSize) {
        continue;
      }
      INT_VECT inter;
      Intersect(ring1, brings[j], inter);
      if (inter.size() > 0 &&
          (!maxOverlapSize || inter.size() <= maxOverlapSize)) {
        neighMap[i].push_back(j);
        neighMap[j].push_back(i);
      }
    }
  }
#if 0
    for (i = 0; i < nrings; i++) {
      std::cerr<<"**************\n    "<<i<<"\n*************\n";
      std::copy(neighMap[i].begin(),neighMap[i].end(),std::ostream_iterator<int>(std::cerr," "));
      std::cerr<<"\n";
    }
#endif
}

}  // end of namespace RingUtils

// local utility namespace
namespace {
using namespace RDKit;

typedef enum {
  VacantElectronDonorType,
  OneElectronDonorType,
  TwoElectronDonorType,
  OneOrTwoElectronDonorType,
  AnyElectronDonorType,
  NoElectronDonorType
} ElectronDonorType;  // used in setting aromaticity
typedef std::vector<ElectronDonorType> VECT_EDON_TYPE;
typedef VECT_EDON_TYPE::iterator VECT_EDON_TYPE_I;
typedef VECT_EDON_TYPE::const_iterator VECT_EDON_TYPE_CI;

/******************************************************************************
 * SUMMARY:
 *  Apply Huckel rule to the specified ring and mark the bonds and atoms in the
 *  ring accordingly.
 *
 * ARGUMENTS:
 *  mol - molecule of interest
 *  ring - list of atoms that form the simple ring
 *  edon - list of electron donor type of the atoms in the molecule
 *
 * RETURN :
 *  none
 *
 * ASSUMPTIONS:
 *  the ring has to a simple ring and the electron donor type of each atom
 *  in the ring have been computed
 ******************************************************************************/
static bool applyHuckel(ROMol &mol, const INT_VECT &ring,
                        const VECT_EDON_TYPE &edon);

static void applyHuckelToFused(
    ROMol &mol,                   // molecule of interest
    const VECT_INT_VECT &srings,  // list of all ring as atom IDS
    const VECT_INT_VECT &brings,  // list of all rings as bond ids
    const INT_VECT &fused,       // list of ring ids in the current fused system
    const VECT_EDON_TYPE &edon,  // electron donor state for each atom
    INT_INT_VECT_MAP &ringNeighs,
    int &narom,  // number of aromatic ring so far
    unsigned int maxNumFusedRings, const std::vector<Bond*>& bondsByIdx, unsigned int minRingSize= 0);

void markAtomsBondsArom(ROMol &mol, const VECT_INT_VECT &srings,
                        const VECT_INT_VECT &brings, const INT_VECT &ringIds,
                        std::set<unsigned int> &doneBonds,
                        const std::vector<Bond*>& bondsByIdx) {
  INT_VECT aring, bring;
  INT_VECT_CI ri, ai, bi;

  // first mark the atoms in the rings
  for (ri = ringIds.begin(); ri != ringIds.end(); ri++) {
    aring = srings[*ri];

    // first mark the atoms in the ring
    for (ai = aring.begin(); ai != aring.end(); ai++) {
      mol.getAtomWithIdx(*ai)->setIsAromatic(true);
    }
  }

  // mark the bonds
  // here we want to be more careful. We don't want to mark the fusing bonds
  // as aromatic - only the outside bonds in a fused system are marked aromatic.
  // - loop through the rings and count the number of times each bond appears in
  //   all the fused rings.
  // - bonds that appears only once are marked aromatic
  INT_MAP_INT bndCntr;
  INT_MAP_INT_I bci;

  for (ri = ringIds.begin(); ri != ringIds.end(); ri++) {
    bring = brings[*ri];
    for (bi = bring.begin(); bi != bring.end(); bi++) {
      if (bndCntr.find(*bi) == bndCntr.end()) {
        bndCntr[*bi] = 1;
      } else {
        bndCntr[*bi] += 1;
      }
    }
  }
  // now mark bonds that have a count of 1 to be aromatic;
  // std::cerr << "bring:";
  for (bci = bndCntr.begin(); bci != bndCntr.end(); ++bci) {
    // std::cerr << " " << bci->first << "(" << bci->second << ")";
    if (bci->second == 1) {
      auto bond = bondsByIdx[bci->first];
      // Bond *bond = mol.get BondWithIdx(bci->first);
      bond->setIsAromatic(true);
      switch (bond->getBondType()) {
        case Bond::SINGLE:
        case Bond::DOUBLE:
          bond->setBondType(Bond::AROMATIC);
          break;
        default:
          break;
      }
      doneBonds.insert(bond->getIdx());
    }
  }
  // std::cerr << std::endl;
}

void getMinMaxAtomElecs(ElectronDonorType dtype, int &atlw, int &atup) {
  switch (dtype) {
    case AnyElectronDonorType:
      atlw = 0;
      atup = 2;
      break;
    case OneOrTwoElectronDonorType:
      atlw = 1;
      atup = 2;
      break;
    case OneElectronDonorType:
      atlw = atup = 1;
      break;
    case TwoElectronDonorType:
      atlw = atup = 2;
      break;
    case NoElectronDonorType:
    case VacantElectronDonorType:
    default:
      atlw = atup = 0;
      break;
  }
}

bool incidentNonCyclicMultipleBond(const Atom *at, int &who) {
  PRECONDITION(at, "bad atom");
  // check if "at" has an non-cyclic multiple bond on it
  // if yes check which atom this bond goes to
  // and record the atomID in who
  const ROMol &mol = at->getOwningMol();
  for(auto *bond : at->bonds()) {
    if (!mol.getRingInfo()->numBondRings(bond->getIdx())) {
      if (bond->getValenceContrib(at) >= 2.0) {
        who = bond->getOtherAtomIdx(at->getIdx());
        return true;
      }
    }
  }
  return false;
}

bool incidentCyclicMultipleBond(const Atom *at) {
  PRECONDITION(at, "bad atom");
  const ROMol &mol = at->getOwningMol();
  for(auto *bond : at->bonds()) {
    if (mol.getRingInfo()->numBondRings(bond->getIdx())) {
      if (bond->getValenceContrib(at) >= 2.0) {
        return true;
      }
    }
  }
  return false;
}

bool incidentMultipleBond(const Atom *at) {
  PRECONDITION(at, "bad atom");
  int deg = at->getDegree() + at->getNumExplicitHs();
  for(auto *bond : at->bonds()) {
    if (!std::lround(bond->getValenceContrib(at))) {
      --deg;
    }
  }
  return at->getExplicitValence() != static_cast<int>(deg);
}

bool applyHuckel(ROMol &mol, const INT_VECT &ring, const VECT_EDON_TYPE &edon,
                 unsigned int minRingSize) {
  RDUNUSED_PARAM(mol);
  if (ring.size() < minRingSize) {
    return false;
  }
  int atlw, atup, rlw, rup, rie;
  bool aromatic = false;
  rlw = 0;
  rup = 0;
  for (auto idx : ring) {
    ElectronDonorType edonType = edon[idx];
    getMinMaxAtomElecs(edonType, atlw, atup);
    rlw += atlw;
    rup += atup;
  }

  if (rup >= 6) {
    for (rie = rlw; rie <= rup; rie++) {
      if ((rie - 2) % 4 == 0) {
        aromatic = true;
        break;
      }
    }
  } else if (rup == 2) {
    aromatic = true;
  }
#if 0
    std::cerr <<" ring: ";
    std::copy(ring.begin(),ring.end(),std::ostream_iterator<int>(std::cerr," "));
    std::cerr <<" rlw: "<<rlw<<" rup: "<<rup<<" aromatic? "<<aromatic<<std::endl;
#endif
  return aromatic;
}

void applyHuckelToFused(
    ROMol &mol,                   // molecule of interest
    const VECT_INT_VECT &srings,  // list of all ring as atom IDS
    const VECT_INT_VECT &brings,  // list of all rings as bond ids
    const INT_VECT &fused,       // list of ring ids in the current fused system
    const VECT_EDON_TYPE &edon,  // electron donor state for each atom
    INT_INT_VECT_MAP &ringNeighs,  // list of neighbors for each candidate ring
    int &narom,                    // number of aromatic ring so far
    unsigned int maxNumFusedRings,
    const std::vector<Bond*>& bondsByIdx,
    unsigned int minRingSize) {
  // this function check huckel rule on a fused system it starts
  // with the individual rings in the system and then proceeds to
  // check for larger system i.e. if we have a 3 ring fused system,
  // huckel rule checked first on all the 1 ring subsystems then 2
  // rung subsystems etc.

  INT_VECT aromRings;
  aromRings.resize(0);
  auto nrings = rdcast<unsigned int>(fused.size());
  INT_VECT curRs;
  INT_VECT_CI mri;
  curRs.push_back(fused.front());
  int pos;
  unsigned int i, curSize = 0;
  INT_VECT comb;
  pos = -1;

  size_t nRingBonds;
  {
    boost::dynamic_bitset<> fusedBonds(mol.getNumBonds());
    for (auto ridx: fused) {
      for (auto bidx: brings[ridx]) {
        fusedBonds[bidx] = true;
      }
    }
    nRingBonds = rdcast<unsigned int>(fusedBonds.count());
  }
  std::set<unsigned int> doneBonds;
  while (1) {
    if (pos == -1) {
      ++curSize;
      // check is we are done with all the atoms in the fused
      // system, if so quit. This is a fix for Issue252 REVIEW: is
      // this check sufficient or should we add an additional
      // constraint on the number of combinations of rings in a
      // fused system that we will try. The number of combinations
      // can obviously be quite large when the number of rings in
      // the fused system is large
      if (curSize > std::min(nrings, maxNumFusedRings) || doneBonds.size() >= nRingBonds) {
        break;
      }
      comb.resize(curSize);
      pos = 0;
      for (i = 0; i < curSize; i++) {
        comb[i] = i;
      }
    } else {
      pos = nextCombination(comb, nrings);
    }

    if (pos == -1) {
      continue;
    }

    curRs.resize(0);
    for (i = 0; i < comb.size(); i++) {
      curRs.push_back(fused[comb[i]]);
    }

    // check if the picked subsystem is fused
    if (ringNeighs.size() && !RingUtils::checkFused(curRs, ringNeighs)) {
      continue;
    }

    // check aromaticity on the current fused system
    INT_VECT atsInRingSystem(mol.getNumAtoms(), 0);
    for (auto ridx : curRs) {
      auto sring = srings[ridx];
      for (const auto rid : sring) {
        atsInRingSystem[rid]++;
      }
    }
    INT_VECT unon;
    for (i = 0; i < atsInRingSystem.size(); ++i) {
      // condition for inclusion of an atom in the aromaticity of a fused ring system
      // is that it's present in one or two of the rings.
      // this was #2895: the central atom in acepentalene was being included in
      // the count of aromatic atoms
      if (atsInRingSystem[i] == 1 || atsInRingSystem[i] == 2) {
        unon.push_back(i);
      }
    }

    if (applyHuckel(mol, unon, edon, minRingSize)) {
      // mark the atoms and bonds in these rings to be aromatic
      markAtomsBondsArom(mol, srings, brings, curRs, doneBonds, bondsByIdx);

      // add the ring IDs to the aromatic rings found so far
      // avoid duplicates
      for (mri = curRs.begin(); mri != curRs.end(); mri++) {
        if (std::find(aromRings.begin(), aromRings.end(), (*mri)) ==
            aromRings.end()) {
          aromRings.push_back(*mri);
        }
      }
    }  // end check huckel rule
  }    // end while(1)
  narom += rdcast<int>(aromRings.size());
}

bool isAtomCandForArom(const Atom *at, const ElectronDonorType edon,
                       bool allowThirdRow = true, bool allowTripleBonds = true,
                       bool allowHigherExceptions = true, bool onlyCorN = false,
                       bool allowExocyclicMultipleBonds = true) {
  PRECONDITION(at, "bad atom");
  if (onlyCorN && at->getAtomicNum() != 6 && at->getAtomicNum() != 7) {
    return false;
  }
  if (!allowThirdRow && at->getAtomicNum() > 10) {
    return false;
  }

  // limit aromaticity to:
  //   - the first two rows of the periodic table
  //   - Se and Te
  if (at->getAtomicNum() > 18 &&
      (!allowHigherExceptions ||
       (at->getAtomicNum() != 34 && at->getAtomicNum() != 52))) {
    return false;
  }
  switch (edon) {
    case VacantElectronDonorType:
    case OneElectronDonorType:
    case TwoElectronDonorType:
    case OneOrTwoElectronDonorType:
    case AnyElectronDonorType:
      break;
    default:
      return (false);
  }

  // atoms that aren't in their default valence state also get shut out
  int defVal = PeriodicTable::getTable()->getDefaultValence(at->getAtomicNum());
  if (defVal > 0 && rdcast<int>(at->getTotalValence()) >
                        (PeriodicTable::getTable()->getDefaultValence(
                            at->getAtomicNum() - at->getFormalCharge()))) {
    return false;
  }

  // heteroatoms or charged carbons with radicals also disqualify us from being
  // considered. This was github issue 432 (heteroatoms) and 1936 (charged
  // carbons)
  if (at->getNumRadicalElectrons() &&
      (at->getAtomicNum() != 6 || at->getFormalCharge())) {
    return false;
  }

  // We are going to explicitly disallow atoms that have more
  // than one double or triple bond. This is to handle
  // the situation:
  //   C1=C=NC=N1 (sf.net bug 1934360)
  int nUnsaturations = at->getExplicitValence() - at->getDegree();
  if (nUnsaturations > 1) {
    unsigned int nMult = 0;
    for(auto *bond : at->bonds()) {
      switch (bond->getBondType()) {
        case Bond::SINGLE:
        case Bond::AROMATIC:
          break;
        case Bond::DOUBLE:
          ++nMult;
          break;
        case Bond::TRIPLE:
          if (!allowTripleBonds) {
            return false;
          }
          ++nMult;
          break;
        default:
          // hopefully we had enough sense that we don't even
          // get here with these bonds... If we do land here,
          // just bail... I have no good answer for them.
          break;
      }
      if (nMult > 1) {
        break;
      }
    }
    if (nMult > 1) {
      return (false);
    }
  }

  if (!allowExocyclicMultipleBonds) {
    for(auto *bnd : at->bonds()) {
      if ((bnd->getBondType() == Bond::DOUBLE ||
           bnd->getBondType() == Bond::TRIPLE) &&
          !queryIsBondInRing(bnd)) {
        return false;
      }
    }
  }

  return (true);
}

ElectronDonorType getAtomDonorTypeArom(
    const Atom *at, bool exocyclicBondsStealElectrons = true) {
  PRECONDITION(at, "bad atom");
  if (at->getAtomicNum() == 0) {
    // dummies can be anything:
    return AnyElectronDonorType;
  }

  ElectronDonorType res = NoElectronDonorType;
  int nelec = MolOps::countAtomElec(at);
  int who = -1;
  const ROMol &mol = at->getOwningMol();
  if (nelec < 0) {
    res = NoElectronDonorType;
  } else if (nelec == 0) {
    if (incidentNonCyclicMultipleBond(at, who)) {
      // This is borderline:  no electron to spare but may have an empty
      // p-orbital
      // Not sure if this should return vacantElectronDonorType
      // FIX: explicitly doing this as a note for potential problems
      //
      res = VacantElectronDonorType;
    } else if (incidentCyclicMultipleBond(at)) {
      // no electron but has one in a in cycle multiple bond
      res = OneElectronDonorType;
    } else {
      // no multiple bonds no electrons
      res = NoElectronDonorType;
    }
  } else if (nelec == 1) {
    if (incidentNonCyclicMultipleBond(at, who)) {
      // the only available electron is going to be from the
      // external multiple bond this electron will not be available
      // for aromaticity if this atom is bonded to a more electro
      // negative atom
      const Atom *at2 = mol.getAtomWithIdx(who);
      if (exocyclicBondsStealElectrons &&
          PeriodicTable::getTable()->moreElectroNegative(at2->getAtomicNum(),
                                                         at->getAtomicNum())) {
        res = VacantElectronDonorType;
      } else {
        res = OneElectronDonorType;
      }
    } else {
      // require that the atom have at least one multiple bond
      if (incidentMultipleBond(at)) {
        res = OneElectronDonorType;
      }
      // account for the tropylium and cyclopropenyl cation cases
      else if (at->getFormalCharge() == 1) {
        res = VacantElectronDonorType;
      }
    }
  } else {
    if (incidentNonCyclicMultipleBond(at, who)) {
      // for cases with more than one electron :
      // if there is an incident multiple bond with an element that
      // is more electronegative than the this atom, count one less
      // electron
      const Atom *at2 = mol.getAtomWithIdx(who);
      if (exocyclicBondsStealElectrons &&
          PeriodicTable::getTable()->moreElectroNegative(at2->getAtomicNum(),
                                                         at->getAtomicNum())) {
        nelec--;
      }
    }
    if (nelec % 2 == 1) {
      res = OneElectronDonorType;
    } else {
      res = TwoElectronDonorType;
    }
  }
  return (res);
}
}  // namespace

namespace RDKit {
namespace MolOps {
int countAtomElec(const Atom *at) {
  PRECONDITION(at, "bad atom");

  // default valence :
  int dv = PeriodicTable::getTable()->getDefaultValence(at->getAtomicNum());
  if (dv <= 1) {
    // univalent elements can't be either aromatic or conjugated
    return -1;
  }

  // total atom degree:
  int degree = at->getDegree() + at->getTotalNumHs();

  for(auto *bond : at->bonds()) {
    // don't count bonds that aren't actually contributing to the valence here:
    if (!std::lround(bond->getValenceContrib(at))) {
      --degree;
    }
  }

  // if we are more than 3 coordinated we should not be aromatic
  if (degree > 3) {
    return -1;
  }

  // number of lone pair electrons = (outer shell elecs) - (default valence)
  int nlp;
  nlp = PeriodicTable::getTable()->getNouterElecs(at->getAtomicNum()) - dv;

  // subtract the charge to get the true number of lone pair electrons:
  nlp = std::max(nlp - at->getFormalCharge(), 0);

  int nRadicals = at->getNumRadicalElectrons();

  // num electrons available for donation into the pi system:
  int res = (dv - degree) + nlp - nRadicals;

  if (res > 1) {
    // if we have an incident bond with order higher than 2,
    // (e.g. triple or higher), we only want to return 1 electron
    // we detect this using the total unsaturation, because we
    // know that there aren't multiple unsaturations (detected
    // above in isAtomCandForArom())
    int nUnsaturations = at->getExplicitValence() - at->getDegree();
    if (nUnsaturations > 1) {
      res = 1;
    }
  }

  return res;
}

namespace {
int mdlAromaticityHelper(RWMol &mol, const VECT_INT_VECT &srings) {
  int narom = 0;
  // loop over all the atoms in the rings that can be candidates
  // for aromaticity
  // Atoms are candidates if
  //   - it is part of ring
  //   - has one or more electron to donate or has empty p-orbitals
  int natoms = mol.getNumAtoms();
  boost::dynamic_bitset<> acands(natoms);
  boost::dynamic_bitset<> aseen(natoms);
  VECT_EDON_TYPE edon(natoms);

  VECT_INT_VECT cRings;  // holder for rings that are candidates for aromaticity
  for (auto &sring : srings) {
    bool allAromatic = true;
    bool allDummy = true;
    for (auto firstIdx : sring) {
      Atom *at = mol.getAtomWithIdx(firstIdx);

      if (allDummy && at->getAtomicNum() != 0) {
        allDummy = false;
      }

      if (aseen[firstIdx]) {
        if (!acands[firstIdx]) {
          allAromatic = false;
        }
        continue;
      }
      aseen[firstIdx] = 1;

      // now that the atom is part of ring check if it can donate
      // electron or has empty orbitals. Record the donor type
      // information in 'edon' - we will need it when we get to
      // the Huckel rule later
      edon[firstIdx] = getAtomDonorTypeArom(at, false);
      // we only accept one electron donors?
      if (edon[firstIdx] != OneElectronDonorType) {
        allAromatic = false;
        continue;
      }
      const bool allowThirdRow = false;
      const bool allowTripleBonds = false;
      const bool allowHigherExceptions = false;
      const bool onlyCorN = true;
      const bool allowExocyclicMultipleBonds = false;
      acands[firstIdx] = isAtomCandForArom(
          at, edon[firstIdx], allowThirdRow, allowTripleBonds,
          allowHigherExceptions, onlyCorN, allowExocyclicMultipleBonds);
      if (!acands[firstIdx]) {
        allAromatic = false;
      }
    }
    if (allAromatic && !allDummy) {
      cRings.push_back(sring);
    }
  }

  // first convert all rings to bonds ids
  VECT_INT_VECT brings;
  RingUtils::convertToBonds(cRings, brings, mol);

  // make the neighbor map for the rings
  // i.e. a ring is a neighbor a another candidate ring if
  // shares at least one bond
  // useful to figure out fused systems
  INT_INT_VECT_MAP neighMap;
  RingUtils::makeRingNeighborMap(brings, neighMap, maxFusedAromaticRingSize, 1);

  // now loop over all the candidate rings and check the
  // huckel rule - of course paying attention to fused systems.
  INT_VECT doneRs;
  int curr = 0;
  int cnrs = rdcast<int>(cRings.size());
  boost::dynamic_bitset<> fusDone(cnrs);
  INT_VECT fused;

  std::vector<Bond*> bondsByIdx;
  bondsByIdx.reserve(mol.getNumBonds());
  for (auto b: mol.bonds()) {
    bondsByIdx.push_back(b);
  }

  while (curr < cnrs) {
    fused.resize(0);
    RingUtils::pickFusedRings(curr, neighMap, fused, fusDone);
    const unsigned int maxFused = 6;
    const unsigned int minRingSize = 6;
    applyHuckelToFused(mol, cRings, brings, fused, edon, neighMap, narom,
                       maxFused, bondsByIdx, minRingSize);

    int rix;
    for (rix = 0; rix < cnrs; rix++) {
      if (!fusDone[rix]) {
        curr = rix;
        break;
      }
    }
    if (rix == cnrs) {
      break;
    }
  }

  mol.setProp(common_properties::numArom, narom, true);

  return narom;
}

// use minRingSize=0 or maxRingSize=0 to ignore these constraints
int aromaticityHelper(RWMol &mol, const VECT_INT_VECT &srings,
                      unsigned int minRingSize, unsigned int maxRingSize,
                      bool includeFused) {
  int narom = 0;
  // loop over all the atoms in the rings that can be candidates
  // for aromaticity
  // Atoms are candidates if
  //   - it is part of ring
  //   - has one or more electron to donate or has empty p-orbitals
  int natoms = mol.getNumAtoms();
  boost::dynamic_bitset<> acands(natoms);
  boost::dynamic_bitset<> aseen(natoms);
  VECT_EDON_TYPE edon(natoms);

  VECT_INT_VECT cRings;  // holder for rings that are candidates for aromaticity
  for (auto &sring : srings) {
    size_t ringSz = sring.size();
    // test ring size:
    if ((minRingSize && ringSz < minRingSize) ||
        (maxRingSize && ringSz > maxRingSize)) {
      continue;
    }

    bool allAromatic = true;
    bool allDummy = true;
    for (auto firstIdx : sring) {
      Atom *at = mol.getAtomWithIdx(firstIdx);

      if (allDummy && at->getAtomicNum() != 0) {
        allDummy = false;
      }

      if (aseen[firstIdx]) {
        if (!acands[firstIdx]) {
          allAromatic = false;
        }
        continue;
      }
      aseen[firstIdx] = 1;

      // now that the atom is part of ring check if it can donate
      // electron or has empty orbitals. Record the donor type
      // information in 'edon' - we will need it when we get to
      // the Huckel rule later
      edon[firstIdx] = getAtomDonorTypeArom(at);
      acands[firstIdx] = isAtomCandForArom(at, edon[firstIdx]);
      if (!acands[firstIdx]) {
        allAromatic = false;
      }
    }
    if (allAromatic && !allDummy) {
      cRings.push_back(sring);
    }
  }

  // first convert all rings to bonds ids
  VECT_INT_VECT brings;
  RingUtils::convertToBonds(cRings, brings, mol);


  std::vector<Bond*> bondsByIdx;
  bondsByIdx.reserve(mol.getNumBonds());
  for (auto b: mol.bonds()) {
    bondsByIdx.push_back(b);
  }

  if (!includeFused) {
    // now loop over all the candidate rings and check the
    // huckel rule - skipping fused systems
    INT_INT_VECT_MAP neighMap;
    for (size_t ri = 0; ri < cRings.size(); ++ri) {
      INT_VECT fused;
      fused.push_back(ri);
      const unsigned int maxFused = 6;
      const unsigned int minRingSize = 0;
      applyHuckelToFused(mol, cRings, brings, fused, edon, neighMap, narom,
                         maxFused, bondsByIdx, minRingSize);
    }
  } else {
    // make the neighbor map for the rings
    // i.e. a ring is a neighbor a another candidate ring if
    // shares at least one bond
    // useful to figure out fused systems
    INT_INT_VECT_MAP neighMap;
    RingUtils::makeRingNeighborMap(brings, neighMap, maxFusedAromaticRingSize, 1);

    // now loop over all the candidate rings and check the
    // huckel rule - of course paying attention to fused systems.
    INT_VECT doneRs;
    int curr = 0;
    int cnrs = rdcast<int>(cRings.size());
    boost::dynamic_bitset<> fusDone(cnrs);
    INT_VECT fused;
    while (curr < cnrs) {
      fused.resize(0);
      RingUtils::pickFusedRings(curr, neighMap, fused, fusDone);
      applyHuckelToFused(mol, cRings, brings, fused, edon, neighMap, narom, 6, bondsByIdx);

      int rix;
      for (rix = 0; rix < cnrs; rix++) {
        if (!fusDone[rix]) {
          curr = rix;
          break;
        }
      }
      if (rix == cnrs) {
        break;
      }
    }
  }

  mol.setProp(common_properties::numArom, narom, true);

  return narom;
}

}  // end of anonymous namespace

int setAromaticity(RWMol &mol, AromaticityModel model, int (*func)(RWMol &)) {
  // This function used to check if the input molecule came
  // with aromaticity information, assumed it is correct and
  // did not touch it. Now it ignores that information entirely.

  // first find the all the simple rings in the molecule
  VECT_INT_VECT srings;
  if (mol.getRingInfo()->isInitialized()) {
    srings = mol.getRingInfo()->atomRings();
  } else {
    MolOps::symmetrizeSSSR(mol, srings);
  }

  int res;
  switch (model) {
    case AROMATICITY_DEFAULT:
    case AROMATICITY_RDKIT:
      res = aromaticityHelper(mol, srings, 0, 0, true);
      break;
    case AROMATICITY_SIMPLE:
      res = aromaticityHelper(mol, srings, 5, 6, false);
      break;
    case AROMATICITY_MDL:
      res = mdlAromaticityHelper(mol, srings);
      break;
    case AROMATICITY_CUSTOM:
      PRECONDITION(
          func,
          "function must be set when aromaticity model is AROMATICITY_CUSTOM");
      res = func(mol);
      break;
    default:
      throw ValueErrorException("Bad AromaticityModel");
  }
  return res;
}

};  // end of namespace MolOps
};  // end of namespace RDKit
