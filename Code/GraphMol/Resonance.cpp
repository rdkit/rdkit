//
//  Copyright (C) 2015 Paolo Tosco
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <iostream>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include "Resonance.h"
#include <boost/functional/hash.hpp>

namespace RDKit {
  namespace ResonanceUtils {
    // depending on the number of atoms which won't have a complete
    // octet (nCandSlots) and the number of those which need non-bonded
    // electrons (nTotalSlots), the number of permutations (numComb) and
    // a binary code (v) which indicates which of the atom indices in
    // aiVec will be octet-unsatisfied for each permutation are computed
    void getNumCombStartV(unsigned int nCandSlots,
      unsigned int nTotalSlots, unsigned int &numComb, unsigned int &v)
    {
      numComb = 1;
      v = 0;
      for (unsigned int i = 0; i < nCandSlots; ++i) {
        numComb *= (nTotalSlots - i);
        numComb /= (i + 1);
        v |= (1 << i);
      }
    }
    
    // get the next permutation
    void updateV(unsigned int &v)
    {
      unsigned int t = (v | (v - 1)) + 1;
      v = t | ((((t & (~t + 1)) / (v & (~v + 1))) >> 1) - 1);
    }

    // sanitize the resonance structure which has been assembled
    void sanitizeMol(RWMol &mol)
    {
      unsigned int opFailed;
      MolOps::sanitizeMol(mol, opFailed, MolOps::SANITIZE_FINDRADICALS
        | MolOps::SANITIZE_SETAROMATICITY | MolOps::SANITIZE_ADJUSTHS);
      if (opFailed != MolOps::SANITIZE_NONE) {
        std::string s;
        switch (opFailed) {
          case MolOps::SANITIZE_FINDRADICALS:
            s = "SANITIZE_FINDRADICALS";
            break;
          case MolOps::SANITIZE_SETAROMATICITY:
            s = "SANITIZE_SETAROMATICITY";
            break;
          case MolOps::SANITIZE_ADJUSTHS:
            s = "SANITIZE_ADJUSTHS";
            break;
        }
        s += " failed";
        throw std::runtime_error(s);
      }
    }
    
    // fix the number of explicit and implicit Hs in the
    // resonance structure which has been assembled
    void fixExplicitImplicitHs(ROMol &mol)
    {
      mol.clearComputedProps(false);
      for (ROMol::AtomIterator ai = mol.beginAtoms();
           ai != mol.endAtoms(); ++ai) {
        (*ai)->clearComputedProps();
        (*ai)->setNumExplicitHs((*ai)->getNumImplicitHs()
             + (*ai)->getNumExplicitHs());
        (*ai)->updatePropertyCache();
      }
    }

  }

  // object constructor
  AtomElectrons::AtomElectrons(ConjElectrons *parent, const Atom *a) :
    d_nb(0),
    d_tv(static_cast<boost::uint8_t>(a->getTotalDegree())),
    d_fc(0),
    d_flags(0),
    d_atom(a),
    d_parent(parent)
  {
  }
  
  // copy constructor
  AtomElectrons::AtomElectrons(ConjElectrons *parent, const AtomElectrons &ae) :
    d_nb(ae.d_nb),
    d_tv(ae.d_tv),
    d_fc(ae.d_fc),
    d_flags(ae.d_flags),
    d_atom(ae.d_atom),
    d_parent(parent)
  {
  }
  
  // assign total valence, formal charge and non-bonded electrons
  // from the original atom
  void AtomElectrons::initTvNbFcFromAtom()
  {
    d_tv = d_atom->getTotalValence();
    d_fc = d_atom->getFormalCharge();
    d_nb = oe() - d_tv - d_fc;
  }
  
  boost::uint8_t AtomElectrons::findAllowedBonds(unsigned int bi)
  {
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
    boost::uint8_t res = 0;
    for (unsigned int i = 0; i < 3; ++i) {
      res |= (canAddBondWithOrder(bi, i + 1) << (i * 2));
    }
    return res;
  }

  // returns true if any conjugated neighbor is charged (and
  // has atomic number == atomicNum, if the latter is non-zero)
  bool AtomElectrons::isNbrCharged(unsigned int bo,
    unsigned int atomicNum)
  {
    bool res = false;
    const ROMol &mol = d_parent->parent()->mol();
    ROMol::OEDGE_ITER nbrIdx, endNbrs;
    boost::tie(nbrIdx, endNbrs) = mol.getAtomBonds(d_atom);
    for (; !res && (nbrIdx != endNbrs); ++nbrIdx) {
      const Bond *bondNbr = mol[*nbrIdx].get();
      unsigned int biNbr = bondNbr->getIdx();
      if (d_parent->parent()->getBondConjGrpIdx(biNbr) != conjGrpIdx())
        continue;
      BondElectrons *beNbr =
        d_parent->getBondElectronsWithIdx(biNbr);
      const Atom *atomNbr = ((bondNbr->getBeginAtom() == d_atom)
        ? bondNbr->getEndAtom() : bondNbr->getBeginAtom());
      unsigned int aiNbr = atomNbr->getIdx();
      AtomElectrons *aeNbr =
        d_parent->getAtomElectronsWithIdx(aiNbr);
      res = (((beNbr->isDefinitive() && !aeNbr->hasOctet())
        || (!beNbr->isDefinitive() && aeNbr->isDefinitive()
        && (aeNbr->oe() < (5 - bo))))
        && (!atomicNum || (atomNbr->getAtomicNum() == atomicNum)));
    }
    return res;
  }
  
  // returns a 2-bit value, where the least significant bit is true
  // if the atom can add a bond with order bo, and the most significant
  // bit is true if the atom needs be charged
  boost::uint8_t AtomElectrons::canAddBondWithOrder(unsigned int bi,
    unsigned int bo)
  {
    boost::uint8_t canAdd = !isDefinitive();
    if (canAdd)
      canAdd = (d_tv <= (5 - bo));
    // if canAdd is true, loop over neighboring conjugated bonds
    // and check their definitive flag; if all neighbors are
    // definitive, we need an additional check on total valence
    if (canAdd && isLastBond()) {
      // we allow a formal charge up to 2 on atoms right of
      // carbon, not more than 1 on atoms left of N
      unsigned int rightOfC = ((d_atom->getAtomicNum() > 6) ? 1 : 0);
      unsigned int fcInc = 0;
      if (rightOfC)
       fcInc = (!isNbrCharged(bo, 6) ? 1 : 0);
      else {
        // atoms left of N can be charged only if:
        // - the neighbor is uncharged and either the conjugate
        //   group does not bear a non-zero total formal charge
        //   or it does and there is no other element left of N
        //   which has already been assigned a formal charge
        // - the neighbor is charged, and triple-bonded to
        //   this atom, which is left of N (e.g., isonitrile)
        bool isnc = isNbrCharged(bo);
        fcInc = (((!isnc && !(!d_parent->allowedChgLeftOfN()
          && d_parent->totalFormalCharge()))
          || (isnc && (bo == 3) && (d_atom->getAtomicNum() < 7)))
          ? 1 : 0);
      }
      unsigned int e = oe() + d_tv - 1 + bo;
      canAdd = ((e + fcInc + rightOfC) >= 8);
      if (canAdd && (e < 8))
        canAdd |= (1 << NEED_CHARGE_BIT);
    }
    return canAdd;
  }

  // sets the LAST_BOND flag on this atom if there is only one
  // non-definitive bond left
  void AtomElectrons::allConjBondsDefinitiveBut(unsigned int bi)
  {
    bool allDefinitive = true;
    ROMol &mol = d_atom->getOwningMol();
    ROMol::OEDGE_ITER nbrIdx, endNbrs;
    boost::tie(nbrIdx, endNbrs) = mol.getAtomBonds(d_atom);
    for (; allDefinitive && (nbrIdx != endNbrs); ++nbrIdx) {
      unsigned int nbi = mol[*nbrIdx].get()->getIdx();
      if ((nbi != bi)
        && (d_parent->parent()->getBondConjGrpIdx(nbi) == conjGrpIdx()))
        allDefinitive =
          d_parent->getBondElectronsWithIdx(nbi)->isDefinitive();
    }
    if (allDefinitive)
      setLastBond();
  }
  
  // called after all bonds to the atom have been marked as definitive
  void AtomElectrons::finalizeAtom()
  {
    // if the atom is left of N and needs non-bonded electrons
    // to achieve the octet, it is the total formal charge
    // of the conjugated group which determines if it is
    // going to be a cation or an anion; once the formal
    // charge is assigned to an atom left of N, the counter of allowed
    // charges on atoms left of N (signed) is decremented by one
    if (d_atom->getAtomicNum() < 7) {
      unsigned int nb = neededNbForOctet();
      if (nb) {
        int fc = nb / 2;
        if (d_parent->allowedChgLeftOfN()) {
          if (d_parent->allowedChgLeftOfN() < 0)
            fc = -fc;
          d_parent->decrAllowedChgLeftOfN(fc);
        }
      }
    }
  }
  
  // object constructor
  BondElectrons::BondElectrons(ConjElectrons *parent, const Bond *b):
    d_bo(1),
    d_flags(0),
    d_bond(b),
    d_parent(parent)
  {
  }
  
  // copy constructor
  BondElectrons::BondElectrons(ConjElectrons *parent, const BondElectrons &be) :
    d_bo(be.d_bo),
    d_flags(be.d_flags),
    d_bond(be.d_bond),
    d_parent(parent)
  {
  }

  // returns the bond order given the bond type
  unsigned int BondElectrons::orderFromBond()
  {
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
        ss << "Bond idx " << d_bond->getIdx()
          << " between atoms " << d_bond->getBeginAtomIdx()
          << " and " << d_bond->getEndAtomIdx()
          << " has an invalid bond type";
        throw std::runtime_error(ss.str());
    }
    return bo;
  }
  
  // set bond order, update total valence on the atoms involved
  // in the bond and update count of current available electrons
  void BondElectrons::setOrder(unsigned int bo)
  {
    d_parent->getAtomElectronsWithIdx
      (d_bond->getBeginAtomIdx())->tvIncr(bo - 1);
    d_parent->getAtomElectronsWithIdx
      (d_bond->getEndAtomIdx())->tvIncr(bo - 1);
    setDefinitive();
    d_parent->decrCurrElectrons(bo * 2);
    d_bo = bo;
  }
  
  // object constructor
  ConjElectrons::ConjElectrons(ResonanceMolSupplier *parent,
    unsigned int groupIdx) :
    d_groupIdx(groupIdx),
    d_totalElectrons(0),
    d_numFormalCharges(0),
    d_totalFormalCharge(0),
    d_absFormalCharges(0),
    d_wtdFormalCharges(0),
    d_fcSameSignDist(0),
    d_fcOppSignDist(0),
    d_nbMissing(0),
    d_flags(0),
    d_parent(parent)
  {
    const ROMol &mol = d_parent->mol();
    unsigned int nb = mol.getNumBonds();
    unsigned int na = mol.getNumAtoms();
    for (unsigned int ai = 0; ai < na; ++ai) {
      if (d_parent->getAtomConjGrpIdx(ai) != na)
        d_totalFormalCharge +=
          mol.getAtomWithIdx(ai)->getFormalCharge();
    }
    d_allowedChgLeftOfN = d_totalFormalCharge;
    for (unsigned int bi = 0; bi < nb; ++bi) {
      if (d_parent->getBondConjGrpIdx(bi) != groupIdx)
        continue;
      const Bond *bond = mol.getBondWithIdx(bi);
      // store the pointers to BondElectrons objects in a map
      d_conjBondMap[bi] = new BondElectrons(this, bond);
      // store the pointers to AtomElectrons objects in a map
      const Atom *atom[2] = {
        bond->getBeginAtom(),
        bond->getEndAtom()
      };
      for (unsigned int i = 0; i < 2; ++i) {
        unsigned int ai = atom[i]->getIdx();
        if (d_conjAtomMap.find(ai) == d_conjAtomMap.end())
          d_conjAtomMap[ai] = new AtomElectrons(this, atom[i]);
      }
    }
    // count total number of valence electrons in conjugated group
    d_currElectrons = countTotalElectrons();
  }

  // copy constructor
  ConjElectrons::ConjElectrons(const ConjElectrons &ce) :
    d_groupIdx(ce.d_groupIdx),
    d_totalElectrons(ce.d_totalElectrons),
    d_currElectrons(ce.d_currElectrons),
    d_numFormalCharges(ce.d_numFormalCharges),
    d_totalFormalCharge(ce.d_totalFormalCharge),
    d_allowedChgLeftOfN(ce.d_allowedChgLeftOfN),
    d_absFormalCharges(ce.d_absFormalCharges),
    d_wtdFormalCharges(ce.d_wtdFormalCharges),
    d_fcSameSignDist(ce.d_fcSameSignDist),
    d_fcOppSignDist(ce.d_fcOppSignDist),
    d_nbMissing(ce.d_nbMissing),
    d_flags(ce.d_flags),
    d_beginAIStack(ce.d_beginAIStack),
    d_parent(ce.d_parent)
  {
    for (ConjAtomMap::const_iterator it = ce.d_conjAtomMap.begin();
      it != ce.d_conjAtomMap.end(); ++it)
      d_conjAtomMap[it->first] = new AtomElectrons(this, *(it->second));
    for (ConjBondMap::const_iterator it = ce.d_conjBondMap.begin();
      it != ce.d_conjBondMap.end(); ++it)
      d_conjBondMap[it->first] = new BondElectrons(this, *(it->second));
  }

  // object destructor
  ConjElectrons::~ConjElectrons()
  {
    for (ConjAtomMap::const_iterator it = d_conjAtomMap.begin();
      it != d_conjAtomMap.end(); ++it) {
      if (it->second)
        delete it->second;
    }
    for (ConjBondMap::const_iterator it = d_conjBondMap.begin();
      it != d_conjBondMap.end(); ++it) {
      if (it->second)
        delete it->second;
    }
  }

  // store fingerprints for this ConjElectrons object in ceMap
  // return true if the FP did not already exist in the map, false
  // if they did (so the ConjElectrons object can be deleted)
  bool ConjElectrons::storeFP(CEMap &ceMap, unsigned int flags)
  {
    boost::uint8_t byte;
    ConjFP fp;
    unsigned int fpSize = 0;
    if (flags & FP_ATOMS)
      fpSize += d_conjAtomMap.size();
    if (flags & FP_BONDS)
      fpSize += ((d_conjBondMap.size() - 1) / 4 + 1);
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
        byte |= (static_cast<boost::uint8_t>
          (it->second->order()) << (i * 2));
        ++i;
      }
      if (i)
        fp.push_back(byte);
    }
    // convert the FP vector to a hash
    std::size_t hash = boost::hash_range(fp.begin(), fp.end());
    // return true if the FP did not already exist in ceMap,
    // false if it did
    return ceMap.insert(std::make_pair(hash, this)).second;
  }
  
  // assign bond orders and formal charges as encoded in the
  // ConjElectrons object to the ROMol passed as reference
  void ConjElectrons::assignBondsFormalChargesToMol(ROMol &mol)
  {
    const Bond::BondType bondType[3] = {
      Bond::SINGLE, Bond::DOUBLE, Bond::TRIPLE
    };
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
        ss << "bond order for bond with index " << bi << " is "
          << be->order() << "; it should be between 1 and 3";
        throw std::runtime_error(ss.str());
      }
      mol.getBondWithIdx(bi)->setBondType(bondType[be->order() - 1]);
    }
  }

  // init atom total valences and bond orders from the
  // respective atoms and bonds
  void ConjElectrons::initCeFromMol()
  {
    for (ConjAtomMap::const_iterator it = d_conjAtomMap.begin();
      it != d_conjAtomMap.end(); ++it)
      it->second->initTvNbFcFromAtom();
    for (ConjBondMap::const_iterator it = d_conjBondMap.begin();
      it != d_conjBondMap.end(); ++it)
      it->second->initOrderFromBond();
    d_currElectrons = 0;
  }
  
  // assign non-bonded electrons to atoms
  void ConjElectrons::assignNonBonded()
  {
    for (ConjAtomMap::const_iterator it = d_conjAtomMap.begin();
      it != d_conjAtomMap.end(); ++it) {
      AtomElectrons *ae = it->second;
      unsigned int nb = std::min
        (ae->neededNbForOctet(), d_currElectrons);
      decrCurrElectrons(nb);
      ae->assignNonBonded(nb);
    }
  }
  
  // assign formal charges to atoms
  void ConjElectrons::assignFormalCharge()
  {
    for (ConjAtomMap::const_iterator it = d_conjAtomMap.begin();
      it != d_conjAtomMap.end(); ++it)
      it->second->assignFormalCharge();
  }
  
  // return true if formal charges for this ConjElectrons
  // object are OK, false if they aren't
  bool ConjElectrons::checkCharges()
  {
    bool areAcceptable = true;
    bool haveIncompleteOctetRightOfC = false;
    bool havePosLeftOfN = false;
    bool haveNegLeftOfN = false;
    bool havePosRightOfNNoOctet = false;
    for (ConjAtomMap::const_iterator it = d_conjAtomMap.begin();
      areAcceptable && (it != d_conjAtomMap.end()); ++it) {
      AtomElectrons *ae = it->second;
      // formal charges should be between -2 and +1
      areAcceptable = ((ae->fc() < 2) && (ae->fc() > -3));
      if (areAcceptable) {
        if (ae->atom()->getAtomicNum() > 6) {
          if (!ae->hasOctet())
            haveIncompleteOctetRightOfC = true;
          if ((ae->fc() > 0) && (ae->atom()->getAtomicNum() > 7)) {
            d_flags |= HAVE_CATION_RIGHT_OF_N;
            if (!ae->hasOctet())
              havePosRightOfNNoOctet = true;
          }
        }
        else {
          if (ae->fc() < 0)
            haveNegLeftOfN = true;
          else if (ae->fc() > 0)
            havePosLeftOfN = true;
        }
        // no carbanions should coexist with incomplete octets
        areAcceptable = !(haveIncompleteOctetRightOfC
          && haveNegLeftOfN);
      }
    }
    if (areAcceptable && !(d_parent->flags()
      & ResonanceMolSupplier::UNCONSTRAINED_CATIONS)) {
      // if the UNCONSTRAINED_CATIONS flag is not set, positively charged
      // atoms left and right of N with an incomplete octet are acceptable
      // only if the conjugated group has a positive total formal charge
      if (havePosLeftOfN || havePosRightOfNNoOctet)
        areAcceptable = (d_totalFormalCharge > 0);
    }
    if (areAcceptable && haveNegLeftOfN && !(d_parent->flags()
      & ResonanceMolSupplier::UNCONSTRAINED_ANIONS))
      // if the UNCONSTRAINED_ANIONS flag is not set, negatively charged
      // atoms left of N are acceptable only if the conjugated group has
      // a negative total formal charge
      areAcceptable = (d_totalFormalCharge < 0);
    for (ConjBondMap::const_iterator it = d_conjBondMap.begin();
      areAcceptable && (it != d_conjBondMap.end()); ++it) {
      BondElectrons *be = it->second;
      AtomElectrons *ae[2] = {
        d_conjAtomMap[be->bond()->getBeginAtomIdx()],
        d_conjAtomMap[be->bond()->getEndAtomIdx()]
      };
      for (unsigned int i = 0; areAcceptable && (i < 2); ++i) {
        if (ae[i]->atom()->getAtomicNum() < 7) {
          // no carbocations allowed on carbons bearing double
          // or triple bonds, carbanions allowed on all carbons
          // no double charged allowed left of N
          areAcceptable = ((be->order() == 1)
            || ((be->order() > 1) && (ae[i]->fc() < 1)));
        }
      }
      if (areAcceptable)
        // charged neighboring atoms left of N are not acceptable
        areAcceptable = !((ae[0]->atom()->getAtomicNum() < 7)
          && ae[0]->fc() && (ae[1]->atom()->getAtomicNum() < 7)
          && ae[1]->fc());
    }
    return areAcceptable;
  }
  
  // assign formal charges and, if the latter acceptable, store
  // return true if FPs did not already exist in ceMap, false if they did
  bool ConjElectrons::assignFormalChargesAndStore
    (CEMap &ceMap, unsigned int fpFlags)
  {
    assignFormalCharge();
    bool ok = checkCharges();
    if (ok)
      ok = storeFP(ceMap, fpFlags);
    if (ok)
      computeMetrics();
    return ok;
  }
  
  // enumerate all possible permutations of non-bonded electrons
  void ConjElectrons::enumerateNonBonded(CEMap &ceMap)
  {
    ConjElectrons *ce = this;
    // the way we compute FPs for a resonance structure depends
    // on whether we want to enumerate all Kekule structures
    // or not; in the first case, we also include bond orders
    // in the FP computation in addition to atom valences
    const unsigned int fpFlags = FP_ATOMS | ((d_parent->flags()
      & ResonanceMolSupplier::KEKULE_ALL) ? FP_BONDS : 0);
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
      ConjElectrons *ceCopy = ((numComb > 1)
        ? new ConjElectrons(*ce) : NULL);
      // enumerate all permutations
      for (unsigned int c = 0; c < numComb; ++c) {
        if (c)
          ce = new ConjElectrons(*ceCopy);
        unsigned int vc = v;
        for (unsigned int i = 0; i < aiVec.size(); ++i) {
          AtomElectrons *ae = ce->getAtomElectronsWithIdx(aiVec[i]);
          unsigned int e = ae->neededNbForOctet();
          // if this atom was chosen to be octet-unsatisfied in
          // this permutation, give it one electron pair less than
          // its need (which most likely means no electrons at all)
          if (vc & 1)
            e -= 2;
          ce->decrCurrElectrons(e);
          ae->assignNonBonded(e);
          vc >>= 1;
        }
        // delete this candidate if it fails the formal charge check
        if (!ce->assignFormalChargesAndStore(ceMap, fpFlags))
          delete ce;
        // get the next binary code
        ResonanceUtils::updateV(v);
      }
      if (ceCopy)
        delete ceCopy;
    }
    else if (nbTotal == currElectrons()) {
      // if the electrons required to satisfy all octets
      // are as many as those currently available, assignment
      // is univocal
      ce->assignNonBonded();
      // delete this candidate if it fails the formal charge check
      if (!ce->assignFormalChargesAndStore(ceMap, fpFlags))
        delete ce;
    }
    else
      // if the electrons required to satisfy all octets are less
      // than those currently available, we must have failed the bond
      // assignment, so the candidate must be deleted
      delete ce;
  }
  
  void ConjElectrons::computeMetrics()
  {
    for (ConjAtomMap::const_iterator it = d_conjAtomMap.begin();
      it != d_conjAtomMap.end(); ++it) {
      d_absFormalCharges += abs(it->second->fc());
      d_wtdFormalCharges += it->second->fc() * it->second->oe();
      d_nbMissing += it->second->neededNbForOctet();
    }
    computeDistFormalCharges();
  }

  // compute sum of shortest path distances between all pairs of
  // formal charges of the same sign and of opposite signs
  void ConjElectrons::computeDistFormalCharges()
  {
    for (ConjAtomMap::const_iterator it1 = d_conjAtomMap.begin();
      it1 != d_conjAtomMap.end(); ++it1) {
      if (!it1->second->fc())
        continue;
      for (ConjAtomMap::const_iterator it2 = it1;
        it2 != d_conjAtomMap.end(); ++it2) {
        if ((it1 == it2) || !it2->second->fc())
          continue;
        unsigned int dist = MolOps::getShortestPath(d_parent->mol(),
          it1->first, it2->first).size();
        if ((it1->second->fc() * it2->second->fc()) > 0)
          d_fcSameSignDist += dist;
        else
          d_fcOppSignDist += dist;
      }
    }
  }
  
  // return the lowest index of atom bearing a formal charge
  unsigned int ConjElectrons::lowestFcIndex() const
  {
    unsigned int na = d_parent->mol().getNumAtoms();
    unsigned int ai = na;
    for (ConjAtomMap::const_iterator it = d_conjAtomMap.begin();
      (ai == na) && (it != d_conjAtomMap.end()); ++it) {
      if (it->second->fc())
        ai = it->first;
    }
    return ai;
  }
  
  // return the lowest index of a multiple bond
  unsigned int ConjElectrons::lowestMultipleBondIndex() const
  {
    unsigned int nb = d_parent->mol().getNumBonds();
    unsigned int bi = nb;
    for (ConjBondMap::const_iterator it = d_conjBondMap.begin();
      (bi == nb) && (it != d_conjBondMap.end()); ++it) {
      if (it->second->order() > 1)
        bi = it->first;
    }
    return bi;
  }
  
  // decrement the count of current electrons by d
  void ConjElectrons::decrCurrElectrons(unsigned int d)
  {
    PRECONDITION((d_currElectrons >= d), "d_currElectrons < d");
    d_currElectrons -= d;
  }
  
  // push atom index ai to the begin stack if it is not already there
  void ConjElectrons::pushToBeginStack(unsigned int ai)
  {
    if (!d_conjAtomMap[ai]->isStacked()) {
      d_conjAtomMap[ai]->setStacked();
      d_beginAIStack.push(ai);
    }
  }
  
  // pop an atom from the begin stack and put its index in ai
  bool ConjElectrons::popFromBeginStack(unsigned int &ai)
  {
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

  // get the pointer to the BondElectrons object for bond having index bi
  BondElectrons *ConjElectrons::getBondElectronsWithIdx
    (unsigned int bi)
  {
    return d_conjBondMap[bi];
  }

// get the pointer to the AtomElectrons object for atom having index ai
  AtomElectrons *ConjElectrons::getAtomElectronsWithIdx
    (unsigned int ai)
  {
    return d_conjAtomMap[ai];
  }

  // count number of total electrons
  unsigned int ConjElectrons::countTotalElectrons()
  {
    // count total number of valence electrons in conjugated group
    PeriodicTable *pt = PeriodicTable::getTable();
    for (ConjBondMap::const_iterator it =
      d_conjBondMap.begin(); it != d_conjBondMap.end(); ++it)
      d_totalElectrons += (2 * it->second->orderFromBond());
    for (ConjAtomMap::const_iterator it =
      d_conjAtomMap.begin(); it != d_conjAtomMap.end(); ++it) {
      const Atom *a = it->second->atom();
      d_totalElectrons += it->second->oe()
        - a->getTotalValence() - a->getFormalCharge();
    }
    return d_totalElectrons;
  }

  // sort CEVectVect by increasing size
  bool ResonanceMolSupplier::vectSizeCompare
    (const std::vector<ConjElectrons *> *a,
    const std::vector<ConjElectrons *> *b)
  {
    if (a->size() != b->size())
      return (a->size() < b->size());
    if (a->size() && b->size())
      return ((*a)[0]->groupIdx() < (*b)[0]->groupIdx());
    return false;
  }

  // if the total number of resonance structures exceeds d_maxStructs,
  // trim the number of ConjElectrons object for each conjugated group
  // to exclude the less likely structures and save memory
  // we are going to generate the complete resonance structures
  // combining the most likely ConjElectrons objects in a breadth-first
  // fashion
  void ResonanceMolSupplier::trimCeVectVect()
  {
    if (d_length > d_maxStructs) {
      d_length = std::min(d_length, d_maxStructs);
      std::vector<unsigned int> s(d_nConjGrp, 0);
      unsigned int currSize = 0;
      while (currSize < d_length) {
        currSize = 1;
        for (unsigned int conjGrpIdx = 0;
          conjGrpIdx < d_nConjGrp; ++conjGrpIdx) {
          if (s[conjGrpIdx] < d_ceVectVect[conjGrpIdx]->size())
            ++s[conjGrpIdx];
          currSize *= s[conjGrpIdx];
        }
      }
      for (unsigned int conjGrpIdx = 0;
        conjGrpIdx < d_nConjGrp; ++conjGrpIdx) {
        for (unsigned int i = s[conjGrpIdx];
          i < d_ceVectVect[conjGrpIdx]->size(); ++i)
          delete (*(d_ceVectVect[conjGrpIdx]))[i];
        d_ceVectVect[conjGrpIdx]->resize(s[conjGrpIdx]);
      }
    }
    std::sort(d_ceVectVect.begin(), d_ceVectVect.end(), vectSizeCompare);
  }
  
  // get the ConjElectrons indices to be combined given
  // the complete resonance structure index
  void ResonanceMolSupplier::idxToCEPerm
    (unsigned int idx, std::vector<unsigned int> &c) const
  {
    unsigned int g = d_nConjGrp;
    while (g) {
      --g;
      unsigned int d = 1;
      for (unsigned int j = 0; j < g; ++j) {
        d *= d_ceVectVect[j]->size();
      }
      c[g] = idx / d;
      idx %= d;
    }
  }

  // sort the vectors of ConjElectrons indices in such a way that we
  // generate resonance structures out of the most likely ConjElectrons
  // objects in a breadth-first fashion
  bool ResonanceMolSupplier::cePermCompare
    (const CEPerm *a, const CEPerm *b)
  {
    unsigned int aSum = 0;
    unsigned int bSum = 0;
    for (unsigned int i = 0; i < a->v.size(); ++i) {
      aSum += a->v[i];
      bSum += b->v[i];
    }
    if (aSum != bSum)
      return (aSum < bSum);
    unsigned int aMax = 0;
    unsigned int bMax = 0;
    for (unsigned int i = 0; i < a->v.size(); ++i) {
      if (!i || (a->v[i] > aMax))
        aMax = a->v[i];
      if (!i || (b->v[i] > bMax))
        bMax = b->v[i];
    }
    if (aMax != bMax)
      return (aMax < bMax);
    for (unsigned int i = 0; i < a->v.size(); ++i) {
      if (a->v[i] != b->v[i])
        return (a->v[i] < b->v[i]);
    }
  }

  // enumerate upfront all the d_length indices for complete resonance
  // structures and sort them to privilege the most likely
  // ConjElectrons index combinations in a breadth-first fashion
  void ResonanceMolSupplier::prepEnumIdxVect()
  {
    d_enumIdx.resize(d_length);
    std::vector<CEPerm *> cePermVect(d_length);
    for (unsigned int i = 0; i < d_length; ++i) {
      cePermVect[i] = new CEPerm;
      cePermVect[i]->idx = i;
      cePermVect[i]->v.resize(d_nConjGrp);
      idxToCEPerm(i, cePermVect[i]->v);
    }
    std::sort(cePermVect.begin(), cePermVect.end(), cePermCompare);
    for (unsigned int i = 0; i < d_length; ++i) {
      d_enumIdx[i] = cePermVect[i]->idx;
      delete cePermVect[i];
    }
  }

  // object constructor
  ResonanceMolSupplier::ResonanceMolSupplier(ROMol &mol,
    unsigned int flags, unsigned int maxStructs) :
    d_nConjGrp(0),
    d_length(1),
    d_flags(flags),
    d_maxStructs(maxStructs),
    d_idx(0)
  {
    d_mol = new ROMol(mol);
    MolOps::Kekulize((RWMol &)*d_mol, false);
    // identify conjugate substructures
    assignConjGrpIdx();
    resizeCeVect();
    for (unsigned int conjGrpIdx = 0;
      conjGrpIdx < d_nConjGrp; ++conjGrpIdx) {
      buildCEMap(conjGrpIdx);
      storeCEMap(conjGrpIdx);
    }
    trimCeVectVect();
    prepEnumIdxVect();
  }
  
  // object destructor
  ResonanceMolSupplier::~ResonanceMolSupplier()
  {
    for (CEVectVect::const_iterator ceVectIt = d_ceVectVect.begin();
      ceVectIt != d_ceVectVect.end();
      ++ceVectIt) {
      for (std::vector<ConjElectrons *>::const_iterator
        ceIt = (*ceVectIt)->begin(); ceIt != (*ceVectIt)->end(); ++ceIt)
        delete(*ceIt);
      delete(*ceVectIt);
    }
    if (d_mol)
      delete d_mol;
  }

  // each bond an atom is assigned an index representing the conjugated
  // group it belongs to; such indices are stored in two vectors
  // (d_bondConjGrpIdx and d_atomConjGrpIdx, respectively)
  void ResonanceMolSupplier::assignConjGrpIdx()
  {
    unsigned int nb = d_mol->getNumBonds();
    d_bondConjGrpIdx.resize(nb, nb);
    unsigned int na = d_mol->getNumAtoms();
    d_atomConjGrpIdx.resize(na, na);
    for (unsigned int i = 0; i < nb; ++i) {
      const Bond *bi = d_mol->getBondWithIdx(i);
      unsigned int biBeginIdx = bi->getBeginAtomIdx();
      unsigned int biEndIdx = bi->getEndAtomIdx();
      if (bi->getIsConjugated() && (d_bondConjGrpIdx[i] == nb)) {
        // assign this conjugate bond to the matching group, if any
        for (unsigned int j = 0;
          (d_bondConjGrpIdx[i] == nb) && (j < nb); ++j) {
          if ((i == j) || (d_bondConjGrpIdx[j] == nb))
            continue;
          const Bond *bj = d_mol->getBondWithIdx(j);
          if ((bj->getBeginAtomIdx() == biBeginIdx)
            || (bj->getBeginAtomIdx() == biEndIdx)
            || (bj->getEndAtomIdx() == biBeginIdx)
            || (bj->getEndAtomIdx() == biEndIdx))
            d_bondConjGrpIdx[i] = d_bondConjGrpIdx[j];
        }
        // no existing group matches: create a new group
        if (d_bondConjGrpIdx[i] == nb)
          d_bondConjGrpIdx[i] = d_nConjGrp++;
      }
    }
    for (unsigned int i = 0; i < nb; ++i) {
      if (d_bondConjGrpIdx[i] != nb) {
        const Bond *bi = d_mol->getBondWithIdx(i);
        unsigned int biBeginIdx = bi->getBeginAtomIdx();
        unsigned int biEndIdx = bi->getEndAtomIdx();
        d_atomConjGrpIdx[biBeginIdx] = d_bondConjGrpIdx[i];
        d_atomConjGrpIdx[biEndIdx] = d_bondConjGrpIdx[i];
      }
    }
  }

  // enumerateNonBonded() is called for each ConjElectrons object
  // retrieved from d_ceMapTmp; the resulting ConjElectrons
  // objects are collected in d_ceMap
  void ResonanceMolSupplier::enumerateNbArrangements()
  {
    for (CEMap::iterator it = d_ceMapTmp.begin();
      it != d_ceMapTmp.end(); ++it)
      it->second->enumerateNonBonded(d_ceMap);
  }

  void ResonanceMolSupplier::pruneStructures()
  {
    unsigned int minNbMissing = 0;
    bool first = true;
    bool haveNoCationsRightOfN = false;
    for (CEMap::const_iterator it = d_ceMap.begin();
      (it != d_ceMap.end()); ++it) {
      if (first || (it->second->nbMissing() < minNbMissing)) {
        first = false;
        minNbMissing = it->second->nbMissing();
      }
      if (!it->second->hasCationRightOfN())
        haveNoCationsRightOfN = true;
    }
    for (CEMap::iterator it = d_ceMap.begin(); it != d_ceMap.end();) {
      // if the flag ALLOW_INCOMPLETE_OCTETS is not set, ConjElectrons
      // objects having less electrons than the most electron-complete
      // structure (which most often will have all complete octets) are
      // discarded
      if ((!(d_flags & ALLOW_INCOMPLETE_OCTETS)
        && (it->second->nbMissing() > minNbMissing))
        || (!(d_flags & UNCONSTRAINED_CATIONS)
        && it->second->hasCationRightOfN() && haveNoCationsRightOfN)) {
        CEMap::iterator toBeDeleted = it;
        ++it;
        delete(toBeDeleted->second);
        d_ceMap.erase(toBeDeleted);
      }
      else
        ++it;
    }
  }

  // function used to sort ConjElectrons objects based on their
  // importance to describe the structure; criteria in order of decreasing
  // priority follow:
  // 1) Number of unsatisfied octets
  // 2) Number of formal charges
  // 3) Number of formal charges weighted by atom electronegativity
  // 4) Distance between formal charges with the same sign
  // 5) Distance between formal charges with opposite signs
  // 6) Index of the first atom bearing a formal charge
  // 7) Index of the first multiple bond
  bool ResonanceMolSupplier::resonanceStructureCompare
    (const ConjElectrons *a, const ConjElectrons *b)
  {
    return ((a->nbMissing() != b->nbMissing())
      ? (a->nbMissing() < b->nbMissing())
      : (a->absFormalCharges() != b->absFormalCharges())
      ? (a->absFormalCharges() < b->absFormalCharges())
      : (a->wtdFormalCharges() != b->wtdFormalCharges())
      ? (a->wtdFormalCharges() < b->wtdFormalCharges())
      : (a->fcSameSignDist() != b->fcSameSignDist())
      ? (a->fcSameSignDist() > b->fcSameSignDist())
      : (a->fcOppSignDist() != b->fcOppSignDist())
      ? (a->fcOppSignDist() > b->fcOppSignDist())
      : (a->lowestFcIndex() != b->lowestFcIndex())
      ? (a->lowestFcIndex() < b->lowestFcIndex())
      : (a->lowestMultipleBondIndex() < b->lowestMultipleBondIndex())
      );
  }
  
  // function which enumerates all possible multiple bond arrangements
  // for each conjugated group, stores each arrangement in a ConjElectrons
  // object ans stores the latter in a map (d_ceMapTmp), keyed with its
  // fingerprints. Depending on whether the KEKULE_ALL flag is set, the FP
  // computation will be based either on the bond arrangement or on atom
  // valences; this provides a convenient mechanism to remove degenerate
  // Kekule structures upfront. In a subsequent step, ConjElectrons objects
  // stored in d_ceMapTmp are further enumerated for their non-bonded
  // electron arrangements, and the final ConjElectrons objects are stored
  // in d_ceMap. Depending on whether the KEKULE_ALL flag is set, the FP
  // computation will involve or not bond arrangement in addition to atom
  // valences
  void ResonanceMolSupplier::buildCEMap(unsigned int conjGrpIdx)
  {
    const unsigned int BEGIN_POS = 0;
    const unsigned int END_POS = 1;
    const unsigned int fpFlags = ConjElectrons::FP_ATOMS
      | ((d_flags & KEKULE_ALL) ? ConjElectrons::FP_BONDS : 0);
    const unsigned int fpFlagsTmp = ((d_flags & KEKULE_ALL)
      ? ConjElectrons::FP_BONDS : ConjElectrons::FP_ATOMS);
    unsigned int nb = d_mol->getNumBonds();
    unsigned int na = d_mol->getNumAtoms();
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
    ConjElectrons *ce = new ConjElectrons(this, conjGrpIdx);
    ConjElectrons *ceCopy = new ConjElectrons(*ce);
    // the first ConjElectrons object has the user-supplied bond
    // and formal charge arrangement and is stored as such
    ce->initCeFromMol();
    // we ignore the result of the call to checkCharges()
    // but we need to call it so that the HAVE_CATION_RIGHT_OF_N flag
    // may eventually be set
    ce->checkCharges();
    ce->computeMetrics();
    ce->storeFP(d_ceMap, fpFlags);
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
          aiFound = (d_atomConjGrpIdx[aiBegin] == conjGrpIdx);
          if (aiFound)
            break;
          ++aiBegin;
        }
        if (!aiFound)
          continue;
        ce->pushToBeginStack(aiBegin);
      }
      // loop until the atom index stack is empty
      while (ce && ce->popFromBeginStack(aiBegin)) {
        // aiBegin holds the atom index just popped from stack
        unsigned int ai[2] = {
          aiBegin, na
        };
        AtomElectrons *ae[2] = {
          ce->getAtomElectronsWithIdx(ai[BEGIN_POS]), NULL
        };
        unsigned int bi = nb;
        BondElectrons *be = NULL;
        // loop over neighbors of the atom popped from the
        // atom index stack
        ROMol::ADJ_ITER nbrIdx, endNbrs;
        boost::tie(nbrIdx, endNbrs) =
          d_mol->getAtomNeighbors(ae[BEGIN_POS]->atom());
        for (; nbrIdx != endNbrs; ++nbrIdx) {
          unsigned int aiNbr = (*d_mol)[*nbrIdx].get()->getIdx();
          // if this neighbor is not part of the conjugated group,
          // ignore it
          if ((ce->parent()->getAtomConjGrpIdx(aiNbr) != conjGrpIdx))
            continue;
          AtomElectrons *aeNbr = ce->getAtomElectronsWithIdx(aiNbr);
          // if we've already dealt with this neighbor before, ignore it
          if (aeNbr->isDefinitive())
            continue;
          unsigned int biNbr = d_mol->getBondBetweenAtoms(ai[BEGIN_POS],
            aiNbr)->getIdx();
          BondElectrons *beNbr = ce->getBondElectronsWithIdx(biNbr);
          // if we have already assigned the bond order to this bond,
          // ignore it
          if (beNbr->isDefinitive())
            continue;
          // if this is the first neighbor we find, process it
          if (!ae[END_POS]) {
            bi = biNbr;
            be = beNbr;
            ai[END_POS] = aiNbr;
            ae[END_POS] = aeNbr;
          }
          // otherwise, if there are multiple neighbors, process the first
          // and store the rest to the atom index stack
          else
            ce->pushToBeginStack(aiNbr);
        }
        // if no neighbors were found, move on to the next atom index
        // in the stack
        if (!be)
          continue;
        boost::uint8_t allowedBondsPerAtom[2];
        // loop over the two atoms that need be bonded
        for (unsigned int i = 0; i < 2; ++i) {
          // for each atom, find which bond orders are allowed
          allowedBondsPerAtom[i] = ae[i]->findAllowedBonds(bi);
          // if for either atom this is the last bond, mark the
          // atom as definitive, otherwise push the atom index
          // to the stack
          if (ae[i]->isLastBond()) {
            ae[i]->setDefinitive();
          }
          else
            ce->pushToBeginStack(ai[i]);
        }
        // logical AND between allowed bond masks for the two atoms
        boost::uint8_t allowedBonds = allowedBondsPerAtom[BEGIN_POS]
          & allowedBondsPerAtom[END_POS];
        bool isAnyBondAllowed = false;
        bool needToFork = false;
        // make a copy of the current ConjElectrons object
        ceCopy = new ConjElectrons(*ce);
        // consider single, double and triple bond alternatives in turn
        for (unsigned int i = 0; i < 3; ++i) {
          unsigned int t = i * 2;
          ConjElectrons *ceToSet = NULL;
          boost::uint8_t orderMask = (1 << t);
          boost::uint8_t chgMask = (1 << (t + 1));
          unsigned int bo = i + 1;
          // if the currently available electrons are enough for this
          // bond type and both atoms can accept this bond type
          // and we are not in a situation where both atoms are left
          // of N and both need a charge to accept this bond type,
          // then this bond type is feasible
          if ((ce->currElectrons() >= (bo * 2))
            && (allowedBonds & orderMask)
            && !((allowedBonds & chgMask)
            && ((ae[BEGIN_POS]->atom()->getAtomicNum() < 7)
            || (ae[END_POS]->atom()->getAtomicNum() < 7)))) {
            isAnyBondAllowed = true;
            // if another bond type will be possible, then we'll
            // need to save that to ceStack and keep on with the current
            // ConjElectrons object
            if (!needToFork) {
              needToFork = true;
              ceToSet = ce;
            }
            else {
              ConjElectrons *ceFork = new ConjElectrons(*ceCopy);
              ceStack.push(ceFork);
              ceToSet = ceFork;
            }
          }
          if (ceToSet) {
            // set bond orders and finalize atoms as needed
            ceToSet->getBondElectronsWithIdx(bi)->setOrder(bo);
            for (unsigned int j = 0; j < 2; ++j) {
              if (ae[j]->isLastBond())
                ceToSet->getAtomElectronsWithIdx(ai[j])->finalizeAtom();
            }
          }
        }
        delete ceCopy;
        // if a dead end was hit, discard this ConjElectrons object
        if (!isAnyBondAllowed) {
          delete ce;
          ce = NULL;
        }
      }
      if (ce) {
        // if this bond arrangement was already stored previously,
        // discard it
        if (!ce->storeFP(d_ceMapTmp, fpFlagsTmp))
          delete ce;
      }
    }
    // for each bond arrangement in d_ceMapTmp, enumerate possible
    // non-bonded electron arrangements, amd store them in d_ceMap
    enumerateNbArrangements();
    // prune structures depending on how flags were set
    pruneStructures();
  }

  // getter function which returns the bondConjGrpIdx for a given
  // bond index
  unsigned int ResonanceMolSupplier::getBondConjGrpIdx
    (unsigned int bi) const
  {
    PRECONDITION(bi < d_bondConjGrpIdx.size(),
      "bi >= d_bondConjGrpIdx.size()");
    return d_bondConjGrpIdx[bi];
  }

  // getter function which returns the atomConjGrpIdx for a given
  // atom index
  unsigned int ResonanceMolSupplier::getAtomConjGrpIdx
    (unsigned int ai) const
  {
    PRECONDITION(ai < d_atomConjGrpIdx.size(),
      "ai >= d_atomConjGrpIdx.size()");
    return d_atomConjGrpIdx[ai];
  }

  // resizes d_ceVectVect vector
  void ResonanceMolSupplier::resizeCeVect()
  {
    d_ceVectVect.resize(d_nConjGrp, NULL);
  }

  // stores the ConjElectrons pointers currently stored in d_ceMap
  // in the d_ceVectVect vector, and clears d_ceMapTmp and d_ceMap
  void ResonanceMolSupplier::storeCEMap(unsigned int conjGrpIdx)
  {
    d_ceVectVect[conjGrpIdx] = new std::vector<ConjElectrons *>();
    d_ceVectVect[conjGrpIdx]->reserve(d_ceMap.size());
    for (CEMap::const_iterator it = d_ceMap.begin();
      it != d_ceMap.end(); ++it)
      d_ceVectVect[conjGrpIdx]->push_back(it->second);
    std::sort(d_ceVectVect[conjGrpIdx]->begin(),
      d_ceVectVect[conjGrpIdx]->end(), resonanceStructureCompare);
    d_length *= d_ceMap.size();
    d_ceMapTmp.clear();
    d_ceMap.clear();
  }
  
  // resets the ResonanceMolSupplier index
  void ResonanceMolSupplier::reset()
  {
    d_idx = 0;
  }

  // returns the next resonance structure as a ROMol *, or NULL
  // if there are no more resonance structures left
  ROMol *ResonanceMolSupplier::next()
  {
    return (atEnd() ? NULL : (*this)[d_idx++]);
  }

  // returns true if there are no more resonance structures left
  bool ResonanceMolSupplier::atEnd() const
  {
    return (d_idx == d_length);
  }

  // sets the ResonanceMolSupplier index to idx
  void ResonanceMolSupplier::moveTo(unsigned int idx)
  {
    PRECONDITION(idx < d_length, "idx >= d_length");
    d_idx = idx;
  }

  // returns the resonance structure with index idx as a ROMol *
  // the index returns resonance structures combining ConjElectrons
  // objects in a breadth-first fashion, in order to return the most
  // likely complete resonance structures first
  ROMol *ResonanceMolSupplier::operator[](unsigned int idx) const
  {
    PRECONDITION(idx < d_length, "idx >= d_length");
    std::vector<unsigned int> c(d_nConjGrp);
    idxToCEPerm(d_enumIdx[idx], c);
    return assignBondsFormalCharges(c);
  }

  // helper function to assign bond orders and formal charges to the
  // mol passed as reference out of the ConjElectrons objects whose
  // indices are passed with the c vector
  void ResonanceMolSupplier::assignBondsFormalChargesHelper
    (ROMol &mol, std::vector<unsigned int> &c) const
  {
    for (unsigned int conjGrpIdx = 0;
      conjGrpIdx < d_nConjGrp; ++conjGrpIdx) {
      ConjElectrons *ce = (*d_ceVectVect[conjGrpIdx])[c[conjGrpIdx]];
      ce->assignBondsFormalChargesToMol(mol);
    }
  }

  // returns a pointer to a ROMol with bond orders and formal charges
  // assigned out of the ConjElectrons objects whose indices are passed
  // with the c vector
  ROMol *ResonanceMolSupplier::assignBondsFormalCharges
    (std::vector<unsigned int> &c) const
  {
    ROMol *mol = new ROMol(this->mol());
    assignBondsFormalChargesHelper(*mol, c);
    ResonanceUtils::fixExplicitImplicitHs(*mol);
    ResonanceUtils::sanitizeMol((RWMol &)*mol);
    return mol;
  }
}
