// $Id$
//
//  Copyright (C) 2013 Paolo Tosco
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include <ForceField/MMFF/Params.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <boost/dynamic_bitset.hpp>

#include <GraphMol/QueryOps.h>
#include "AtomTyper.h"
#include <cstdarg>


namespace RDKit {
  namespace MMFF {
    using namespace ForceFields::MMFF;

    // given the atomic num, this function returns the periodic
    // table row number, starting from 0 for hydrogen
    unsigned int getPeriodicTableRow(const int atomicNum)
    {
      unsigned int periodicTableRow = 0;
      
      if ((atomicNum >= 3) && (atomicNum <= 10)) {
        periodicTableRow = 1;
      }
      else if ((atomicNum >= 11) && (atomicNum <= 18)) {
        periodicTableRow = 2;
      }
      else if ((atomicNum >= 19) && (atomicNum <= 36)) {
        periodicTableRow = 3;
      }
      else if ((atomicNum >= 37) && (atomicNum <= 54)) {
        periodicTableRow = 4;
      }
      
      return periodicTableRow;
    }
    
    // given the atomic num, this function returns the periodic
    // table row number, starting from 1 for helium
    // Hydrogen has a special row number (0), while transition
    // metals have the row number multiplied by 10
    unsigned int getPeriodicTableRowHL(const int atomicNum)
    {
      unsigned int periodicTableRow = 0;
      
      if (atomicNum == 2) {
        periodicTableRow = 1;
      }
      else if ((atomicNum >= 3) && (atomicNum <= 10)) {
        periodicTableRow = 2;
      }
      else if ((atomicNum >= 11) && (atomicNum <= 18)) {
        periodicTableRow = 3;
      }
      else if ((atomicNum >= 19) && (atomicNum <= 36)) {
        periodicTableRow = 4;
      }
      else if ((atomicNum >= 37) && (atomicNum <= 54)) {
        periodicTableRow = 5;
      }
      if (((atomicNum >= 21) && (atomicNum <= 30))
        || ((atomicNum >= 39) && (atomicNum <= 48))
        || ((atomicNum >= 39) && (atomicNum <= 48))) {
        periodicTableRow *= 10;
      }
      
      return periodicTableRow;
    }
    
    // given the MMFF atom type, this function returns true
    // if it is aromatic
    bool isAromaticAtomType(const unsigned int atomType)
    {
      const unsigned int aromatic_array[] =
        { 37, 38, 39, 44, 58, 59, 63, 64, 65, 66, 69, 76, 78, 79, 80, 81, 82 };
      const std::vector<int> aromaticTypes(aromatic_array,
        aromatic_array + sizeof(aromatic_array) / sizeof(aromatic_array[0]));
      
      return ((std::find(aromaticTypes.begin(), aromaticTypes.end(),
        atomType) != aromaticTypes.end()) ? true : false);
    }
    
    
    // returns true if the atom is in a ring of size ringSize
    bool isAtomInAromaticRingOfSize(const Atom *atom, const unsigned int ringSize)
    {
      unsigned int i;
      unsigned int j;
      bool isAromatic = false;
      ROMol mol = atom->getOwningMol();
      VECT_INT_VECT atomRings = mol.getRingInfo()->atomRings();

      if (atom->getIsAromatic()) {
        for (i = 0; (!isAromatic) && (i < atomRings.size()); ++i) {
          if ((atomRings[i].size() != ringSize) || (std::find(atomRings[i].begin(),
            atomRings[i].end(), atom->getIdx()) == atomRings[i].end())) {
            continue;
          }
          for (j = 0, isAromatic = true; isAromatic && (j < atomRings[i].size() - 1); ++j) {
            isAromatic = (mol.getBondBetweenAtoms(atomRings[i][j],
              atomRings[i][j + 1])->getBondType() == Bond::AROMATIC);
          }
        }
      }
      
      return isAromatic;
    }


    // returns true if the atom is an N-oxide
    bool isAtomNOxide(const Atom *atom)
    {
      bool isNOxide = false;
      ROMol mol = atom->getOwningMol();
      ROMol::ADJ_ITER nbrIdx;
      ROMol::ADJ_ITER endNbrs;

      if ((atom->getAtomicNum() == 7) && (atom->getTotalDegree() >= 3)) {
        // loop over neighbors
        boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
        for (; (!isNOxide) && (nbrIdx != endNbrs); ++nbrIdx) {
          const Atom *nbrAtom = mol[*nbrIdx].get();
          isNOxide = ((nbrAtom->getAtomicNum() == 8) && (nbrAtom->getTotalDegree() == 1));
        }
      }
        
      return isNOxide;
    }


    // if the angle formed by atoms with indexes idx1, idx2, idx3
    // is in a ring of {3,4} atoms returns 3 or 4, respectively;
    // otherwise it returns 0
    unsigned int isAngleInRingOfSize3or4(const ROMol &mol, const unsigned int idx1,
      const unsigned int idx2, const unsigned int idx3)
    {
      unsigned int ringSize = 0;
      
      if (mol.getBondBetweenAtoms(idx1, idx2)
        && mol.getBondBetweenAtoms(idx2, idx3)) {
        if (mol.getBondBetweenAtoms(idx3, idx1)) {
          ringSize = 3;
        }
        else {
          std::set<unsigned int> s1;
          std::set<unsigned int> s2;
          std::vector<int> intersect;
          ROMol::ADJ_ITER nbrIdx;
          ROMol::ADJ_ITER endNbrs;
          unsigned int newIdx;
          boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(mol.getAtomWithIdx(idx1));
          for (; nbrIdx != endNbrs; ++nbrIdx) {
            newIdx = mol[*nbrIdx].get()->getIdx();
            if (newIdx != idx2) {
              s1.insert(newIdx);
            }
          }
          boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(mol.getAtomWithIdx(idx3));
          for (; nbrIdx != endNbrs; ++nbrIdx) {
            newIdx = mol[*nbrIdx].get()->getIdx();
            if (newIdx != idx2) {
              s2.insert(newIdx);
            }
          }
          std::set_intersection(s1.begin(), s1.end(),
            s2.begin(), s2.end(), std::back_inserter(intersect));
          if (intersect.size()) {
            ringSize = 4;
          }
        }
      }
      
      return ringSize;
    }
    
    
    // if the dihedral angle formed by atoms with indexes idx1,
    // idx2, idx3, idx4 is in a ring of {4,5} atoms returns 4 or 5,
    // respectively; otherwise it returns 0
    unsigned int isTorsionInRingOfSize4or5(const ROMol &mol, const unsigned int idx1,
      const unsigned int idx2, const unsigned int idx3, const unsigned int idx4)
    {
      unsigned int ringSize = 0;
      
      if (mol.getBondBetweenAtoms(idx1, idx2)
        && mol.getBondBetweenAtoms(idx2, idx3)
        && mol.getBondBetweenAtoms(idx3, idx4)) {
        if (mol.getBondBetweenAtoms(idx4, idx1)) {
          ringSize = 4;
        }
        else {
          std::set<unsigned int> s1;
          std::set<unsigned int> s2;
          std::vector<int> intersect;
          ROMol::ADJ_ITER nbrIdx;
          ROMol::ADJ_ITER endNbrs;
          unsigned int newIdx;
          boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(mol.getAtomWithIdx(idx1));
          for (; nbrIdx != endNbrs; ++nbrIdx) {
            newIdx = mol[*nbrIdx].get()->getIdx();
            if (newIdx != idx2) {
              s1.insert(newIdx);
            }
          }
          boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(mol.getAtomWithIdx(idx4));
          for (; nbrIdx != endNbrs; ++nbrIdx) {
            newIdx = mol[*nbrIdx].get()->getIdx();
            if (newIdx != idx3) {
              s2.insert(newIdx);
            }
          }
          std::set_intersection(s1.begin(), s1.end(),
            s2.begin(), s2.end(), std::back_inserter(intersect));
          if (intersect.size()) {
            ringSize = 5;
          }
        }
      }
      
      return ringSize;
    }

    
    // return true if atoms are in the same ring of size ringSize
    bool areAtomsInSameRingOfSize(const ROMol &mol,
      const unsigned int ringSize, const unsigned int numAtoms, ...)
    {
      unsigned int i;
      bool areInSameRingOfSize = false;
      VECT_INT_VECT atomRings = mol.getRingInfo()->atomRings();
      unsigned int idx;
      va_list atomIdxs;

      for (i = 0; (!areInSameRingOfSize) && (i < atomRings.size()); ++i) {
        if (atomRings[i].size() != ringSize) {
          continue;
        }
        areInSameRingOfSize = true;
        va_start(atomIdxs, numAtoms);
        for (unsigned int j = 0; areInSameRingOfSize && (j < numAtoms); ++j) {
          idx = va_arg(atomIdxs, unsigned int);
          areInSameRingOfSize = (std::find(atomRings[i].begin(),
            atomRings[i].end(), idx) != atomRings[i].end());
        }
        va_end(atomIdxs);
      }
      
      return areInSameRingOfSize;
    }


    // return true if atoms are in the same aromatic ring
    bool areAtomsInSameAromaticRing(const ROMol &mol,
      const unsigned int idx1, const unsigned int idx2)
    {
      unsigned int i;
      unsigned int j;
      bool areInSameAromatic = false;
      VECT_INT_VECT atomRings = mol.getRingInfo()->atomRings();

      if (mol.getAtomWithIdx(idx1)->getIsAromatic()
        && mol.getAtomWithIdx(idx2)->getIsAromatic()) {
        for (i = 0; (!areInSameAromatic) && (i < atomRings.size()); ++i) {
          if ((std::find(atomRings[i].begin(),
            atomRings[i].end(), idx1) != atomRings[i].end())
            && (std::find(atomRings[i].begin(),
            atomRings[i].end(), idx2) != atomRings[i].end())) {
            areInSameAromatic = true;
            for (j = 0; areInSameAromatic && (j < (atomRings[i].size() - 1)); ++j) {
              areInSameAromatic = (mol.getBondBetweenAtoms
                (atomRings[i][j], atomRings[i][j + 1])->getBondType() == Bond::AROMATIC);
            }
          }
        }
      }
      
      return areInSameAromatic;
    }


    // sets the aromaticity flags according to MMFF
    void setMMFFAromaticity(RWMol &mol)
    {
      bool moveToNextRing = false;
      bool isNOSinRing = false;
      bool aromRingsAllSet = false;
      bool exoDoubleBond = false;
      bool canBeAromatic = false;
      unsigned int i;
      unsigned int j;
      unsigned int nextInRing;
      unsigned int pi_e = 0;
      int nAromSet = 0;
      int old_nAromSet = -1;
      RingInfo *ringInfo = mol.getRingInfo();
      Atom *atom;
      Bond *bond;
      VECT_INT_VECT atomRings = ringInfo->atomRings();
      ROMol::ADJ_ITER nbrIdx;
      ROMol::ADJ_ITER endNbrs;
      boost::dynamic_bitset<> aromBitVect(mol.getNumAtoms());
      boost::dynamic_bitset<> aromRingBitVect(atomRings.size());
      
      while ((!aromRingsAllSet) && atomRings.size() && (nAromSet > old_nAromSet)) {
        // loop over all rings
        for (i = 0; i < atomRings.size(); ++i) {
          // add 2 pi electrons for each double bond in the ring
          for (j = 0, pi_e = 0, moveToNextRing = false, isNOSinRing = false,
            exoDoubleBond = false; (!moveToNextRing) && (j < atomRings[i].size()); ++j) {
            atom = mol.getAtomWithIdx(atomRings[i][j]);
            // remember if this atom is nitrogen, oxygen or divalent sulfur
            if ((atom->getAtomicNum() == 7) || (atom->getAtomicNum() == 8)
              || ((atom->getAtomicNum() == 16) && (atom->getDegree() == 2))) {
              isNOSinRing = true;
            }
            // check whether this atom is double-bonded to next one in the ring
            nextInRing = (j == (atomRings[i].size() - 1))
              ? atomRings[i][0] : atomRings[i][j + 1];
            if (mol.getBondBetweenAtoms(atomRings[i][j],
              nextInRing)->getBondType() == Bond::DOUBLE) {
              pi_e += 2;
            }
            // if this is not a double bond, check whether this is carbon
            // or nitrogen with total bond order = 4
            else {
              atom = mol.getAtomWithIdx(atomRings[i][j]);
              // if not, move on
              if ((atom->getAtomicNum() != 6) && (!((atom->getAtomicNum() == 7)
                && ((atom->getExplicitValence() + atom->getNumImplicitHs()) == 4)))) {
                continue;
              }
              // loop over neighbors
              boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
              for (; nbrIdx != endNbrs; ++nbrIdx) {
                const Atom *nbrAtom = mol[*nbrIdx].get();
                // if the neighbor is one of the ring atoms, skip it
                // since we are looking for exocyclic neighbors
                if (std::find(atomRings[i].begin(), atomRings[i].end(),
                  nbrAtom->getIdx()) != atomRings[i].end()) {
                  continue;
                }
                // it the neighbor is single-bonded, skip it
                if (mol.getBondBetweenAtoms(atomRings[i][j],
                  nbrAtom->getIdx())->getBondType() == Bond::SINGLE) {
                  continue;
                }
                // if the neighbor is in a ring and its aromaticity
                // bit has not yet been set, then move to the next ring
                // we'll take care of this later
                if (queryIsAtomInRing(nbrAtom)
                  && (!(aromBitVect[nbrAtom->getIdx()]))) {
                  moveToNextRing = true;
                  break;
                }
                // if the neighbor is in an aromatic ring and is
                // double-bonded to the current atom, add 1 pi electron
                if (mol.getBondBetweenAtoms(atomRings[i][j],
                  nbrAtom->getIdx())->getBondType() == Bond::DOUBLE) {
                  if (nbrAtom->getIsAromatic()) {
                    ++pi_e;
                  }
                  else {
                    exoDoubleBond = true;
                  }
                }
              }
            }
          }
          // if we quit the loop at an early stage because aromaticity
          // had not yet been set, then move to the next ring
          if (moveToNextRing) {
            continue;
          }
          // loop again over all ring atoms
          for (j = 0, canBeAromatic = true; j < atomRings[i].size(); ++j) {
            // set aromaticity as perceived
            aromBitVect[atomRings[i][j]] = 1;
            atom = mol.getAtomWithIdx(atomRings[i][j]);
            // if this is is a non-sp2 carbon or nitrogen
            // then this ring can't be aromatic
            if (((atom->getAtomicNum() == 6) || (atom->getAtomicNum() == 7))
              && (atom->getHybridization() != Atom::SP2)) {
              canBeAromatic = false;
            }
          }
          // if this ring can't be aromatic, move to the next one
          if (!canBeAromatic) {
            continue;
          }
          // if there is N, O, S; no exocyclic double bonds;
          // the ring has an odd number of terms: add 2 pi electrons
          if (isNOSinRing && (!exoDoubleBond) && (atomRings[i].size() % 2)) {
            pi_e += 2;
          }
          // if this ring satisfies the 4n+2 rule,
          // then mark its atoms as aromatic
          if ((pi_e > 2) && (!((pi_e - 2) % 4))) {
            aromRingBitVect[i] = 1;
            for (j = 0; j < atomRings[i].size(); ++j) {
              atom=mol.getAtomWithIdx(atomRings[i][j]);
              atom->setIsAromatic(true);
              if(atom->getAtomicNum()!=6) {
                //                std::cerr<<"   orig: "<<atom->getNumExplicitHs()<<std::endl;
#if 1
                atom->calcImplicitValence(false);
                int iv=atom->getImplicitValence();
                if(iv){
                  atom->setNumExplicitHs(iv);
                  atom->calcImplicitValence(false);
                }
#endif
              }
            }
          }
        }
        // termination criterion: if we did not manage to set any more
        // aromatic atoms compared to the previous iteration, then
        // stop looping
        old_nAromSet = nAromSet;
        nAromSet = 0;
        aromRingsAllSet = true;
        for (i = 0; i < atomRings.size(); ++i) {
          for (j = 0; j < atomRings[i].size(); ++j) {
            if (aromBitVect[atomRings[i][j]]) {
              ++nAromSet;
            }
            else {
              aromRingsAllSet = false;
            }
          }
        }
      }
      for (i = 0; i < atomRings.size(); ++i) {
        // if the ring is not aromatic, move to the next one
        if (!aromRingBitVect[i]) {
          continue;
        }
        for (j = 0; j < atomRings[i].size(); ++j) {
          // mark all ring bonds as aromatic
          nextInRing = (j == (atomRings[i].size() - 1))
            ? atomRings[i][0] : atomRings[i][j + 1];
          bond = mol.getBondBetweenAtoms(atomRings[i][j], nextInRing);
          bond->setBondType(Bond::AROMATIC);
          bond->setIsAromatic(true);
        }
      }
    }


    // sets the MMFF atomType for a heavy atom
    void MMFFMolProperties::setMMFFHeavyAtomType(const Atom *atom)
    {
      unsigned int atomType = 0;
      unsigned int i;
      unsigned int j;
      unsigned int nTermObondedToN = 0;
      bool alphaOrBetaInSameRing = false;
      bool isAlphaOS = false;
      bool isBetaOS = false;
      bool isNSO2orNSO3orNCN = false;
      ROMol mol = atom->getOwningMol();
      RingInfo *ringInfo = mol.getRingInfo();
      ROMol::ADJ_ITER nbrIdx;
      ROMol::ADJ_ITER endNbrs;
      ROMol::ADJ_ITER nbr2Idx;
      ROMol::ADJ_ITER end2Nbrs;
      ROMol::ADJ_ITER nbr3Idx;
      ROMol::ADJ_ITER end3Nbrs;
      std::vector<const Atom *> alphaHet;
      std::vector<const Atom *> betaHet;


      if (atom->getIsAromatic()) {
        if (isAtomInAromaticRingOfSize(atom, 5)) {
          // 5-membered aromatic rings
          // if ipso is carbon or nitrogen, find eventual alpha and beta heteroatoms
          if ((atom->getAtomicNum() == 6) || (atom->getAtomicNum() == 7)) {
            // loop over alpha neighbors
            boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
            for (; nbrIdx != endNbrs; ++nbrIdx) {
              const Atom *nbrAtom = mol[*nbrIdx].get();
              // if the alpha neighbor is not in a 5-membered aromatic
              // ring, skip to the next neighbor
              if (!isAtomInAromaticRingOfSize(nbrAtom, 5)) {
                continue;
              }
              // if the alpha neighbor belongs to the same ring of ipso atom
              // and it is either oxygen, sulfur, or non-N-oxide trivalent nitrogen,
              // add it to the alpha atom vector
              if (areAtomsInSameRingOfSize(mol, 5, 2, atom->getIdx(), nbrAtom->getIdx())
                && ((nbrAtom->getAtomicNum() == 8)
                || (nbrAtom->getAtomicNum() == 16) || ((nbrAtom->getAtomicNum() == 7)
                && (nbrAtom->getTotalDegree() == 3) && (!isAtomNOxide(nbrAtom))))) {
                alphaHet.push_back(nbrAtom);
              }
              // loop over beta neighbors
              boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(nbrAtom);
              for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
                const Atom *nbr2Atom = mol[*nbr2Idx].get();
                // if we have gone back to the ipso atom, move on
                if (nbr2Atom->getIdx() == atom->getIdx()) {
                  continue;
                }
                // if the beta neighbor is not in a 5-membered aromatic
                // ring, skip to the next neighbor
                if (!isAtomInAromaticRingOfSize(nbr2Atom, 5)) {
                  continue;
                }
                // if the beta neighbor belongs to the same ring of ipso atom
                // and it is either oxygen, sulfur, or non-N-oxide trivalent nitrogen,
                // add it to the beta atom vector
                if (areAtomsInSameRingOfSize(mol, 5, 2, atom->getIdx(), nbr2Atom->getIdx())
                  && ((nbr2Atom->getAtomicNum() == 8)
                  || (nbr2Atom->getAtomicNum() == 16) || ((nbr2Atom->getAtomicNum() == 7)
                  && (nbr2Atom->getTotalDegree() == 3) && (!isAtomNOxide(nbr2Atom))))) {
                  betaHet.push_back(nbr2Atom);
                }
              }
            }
            isAlphaOS = false;
            for (i = 0; (!isAlphaOS) && (i < alphaHet.size()); ++i) {
              isAlphaOS = ((alphaHet[i]->getAtomicNum() == 8)
                || (alphaHet[i]->getAtomicNum() == 16));
            }
            isBetaOS = false;
            for (i = 0; (!isBetaOS) && (i < betaHet.size()); ++i) {
              isBetaOS = ((betaHet[i]->getAtomicNum() == 8)
                || (betaHet[i]->getAtomicNum() == 16));
            }
            if (alphaHet.size() && betaHet.size()) {
              // do alpha and beta heteroatoms belong to the same ring?
              for (i = 0; (!alphaOrBetaInSameRing) && (i < alphaHet.size()); ++i) {
                for (j = 0; (!alphaOrBetaInSameRing) && (j < betaHet.size()); ++j) {
                  alphaOrBetaInSameRing = areAtomsInSameRingOfSize
                    (mol, 5, 2, alphaHet[i]->getIdx(), betaHet[j]->getIdx());
                }
              }
            }
          }

          switch (atom->getAtomicNum()) {
            // Carbon
            case 6:
            // if there are no beta heteroatoms
            if (!(betaHet.size())) {
              // count how many 3-neighbor nitrogens we have
              // to be CIM+, there must be at least two such nitrogens,
              // one of which in a 5-membered aromatic ring and none
              // in a 6-membered aromatic ring; additionally, one
              // one of the hydrogens must be protonated
              unsigned int nN = 0;
              unsigned int nFormalCharge = 0;
              unsigned int nInAromatic5Ring = 0;
              unsigned int nInAromatic6Ring = 0;
              boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
              for (; nbrIdx != endNbrs; ++nbrIdx) {
                const Atom *nbrAtom = mol[*nbrIdx].get();
                if ((nbrAtom->getAtomicNum() == 7)
                  && (nbrAtom->getTotalDegree() == 3)) {
                  ++nN;
                  if ((nbrAtom->getFormalCharge() > 0)
                    && (!isAtomNOxide(nbrAtom))) {
                    ++nFormalCharge;
                  }
                  if (isAtomInAromaticRingOfSize(nbrAtom, 5)) {
                    ++nInAromatic5Ring;
                  }
                  if (isAtomInAromaticRingOfSize(nbrAtom, 6)) {
                    ++nInAromatic6Ring;
                  }
                }
              }
              if ((((nN == 2) && nInAromatic5Ring) || ((nN == 3) && (nInAromatic5Ring == 2)))
                && nFormalCharge && (!nInAromatic6Ring)) {
                // CIM+
                // Aromatic carbon between N's in imidazolium
                atomType = 80;
                break;
              }
            }
            // if there are neither alpha nor beta heteroatoms
            // or if there are both, but they belong to different rings
            if (((!(alphaHet.size())) && (!(betaHet.size())))
              || (alphaHet.size() && betaHet.size())) {
              bool surroundedByBenzeneC = true;
              bool surroundedByArom = true;
              // loop over neighbors
              // are all neighbors aromatic?
              // are all neighbors benzene carbons?
              boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
              for (; nbrIdx != endNbrs; ++nbrIdx) {
                const Atom *nbrAtom = mol[*nbrIdx].get();
                
                if ((nbrAtom->getAtomicNum() != 6)
                  || (!(ringInfo->isAtomInRingOfSize(nbrAtom->getIdx(), 6)))) {
                  surroundedByBenzeneC = false;
                }
                if (areAtomsInSameRingOfSize(mol, 5, 2, atom->getIdx(), nbrAtom->getIdx())
                  && (!(nbrAtom->getIsAromatic()))) {
                  surroundedByArom = false;
                }
              }
              // if there are no alpha and beta heteroatoms and
              // all neighbors are aromatic but not all of them
              // benzene carbons, or if there are alpha and beta
              // atoms but they belong to different rings, or if
              // there are alpha and beta heteroatoms but no alpha
              // oxygen or sulfur, then it's C5
              if (((!(alphaHet.size())) && (!(betaHet.size()))
                && (!surroundedByBenzeneC) && surroundedByArom)
                || (alphaHet.size() && betaHet.size()
                && ((!alphaOrBetaInSameRing) || ((!isAlphaOS) && (!isBetaOS))))) {
                // C5
                // General carbon in 5-membered heteroaromatic ring
                atomType = 78;
                break;
              }
            }
            if (alphaHet.size() && ((!(betaHet.size())) || isAlphaOS)) {
              // C5A
              // Aromatic 5-ring C, alpha to N:, O: or S:
              atomType = 63;
              break;
            }
            if (betaHet.size() && ((!(alphaHet.size())) || isBetaOS)) {
              // C5B
              // Aromatic 5-ring C, alpha to N:, O: or S:
              atomType = 64;
              break;
            }
            break;
          
          // Nitrogen
            case 7:
            if (isAtomNOxide(atom)) {
              // N5AX
              // N-oxide nitrogen in 5-ring alpha position
              // N5BX
              // N-oxide nitrogen in 5-ring beta position
              // N5OX
              // N-oxide nitrogen in other 5-ring position
              atomType = 82;
              break;
            }
            // if there are neither alpha nor beta heteroatoms
            if ((!(alphaHet.size())) && (!(betaHet.size()))) {
              // if it is nitrogen
              // if valence is 3, it's pyrrole nitrogen
              if (atom->getTotalDegree() == 3) {
                // NPYL
                // Aromatic 5-ring nitrogen with pi lone pair
                atomType = 39;
                break;
              }
              // otherwise it is anionic
              // N5M
              // Nitrogen in 5-ring aromatic anion
              atomType = 76;
              break;
            }
            if ((atom->getTotalDegree() == 3)
              && ((alphaHet.size() && (!(betaHet.size())))
              || (betaHet.size() && (!(alphaHet.size()))))) {
              // NIM+
              // Aromatic nitrogen in imidazolium
              // N5A+
              // Positive nitrogen in 5-ring alpha position
              // N5B+
              // Positive nitrogen in 5-ring beta position
              // N5+
              // Positive nitrogen in other 5-ring position
              atomType = 81;
              break;
            }
            // if there are alpha heteroatoms and either no beta heteroatoms
            // or no alpha oxygen/sulfur
            if (alphaHet.size() && ((!(betaHet.size())) || isAlphaOS)) {
              // N5A
              // Aromatic 5-ring N, alpha to N:, O: or S:
              atomType = 65;
              break;
            }
            // if there are beta heteroatoms and either no alpha heteroatoms
            // or no beta oxygen/sulfur
            if (betaHet.size() && ((!(alphaHet.size())) || isBetaOS)) {
              // N5B
              // Aromatic 5-ring N, beta to N:, O: or S:
              atomType = 66;
              break;
            }
            // if there are both alpha and beta heteroatoms
            if (alphaHet.size() && betaHet.size()) {
              // N5
              // General nitrogen in 5-memebered heteroaromatic ring
              atomType = 79;
              break;
            }
            break;
          
          // Oxygen
            case 8:
              // OFUR
              // Aromatic 5-ring oxygen with pi lone pair
              atomType = 59;
            break;
          
          // Sulfur
            case 16:
              // STHI
              // Aromatic 5-ring sulfur with pi lone pair
              atomType = 44;
            break;
          }
        }
        
        if ((!atomType) && (isAtomInAromaticRingOfSize(atom, 6))) {
          // 6-membered aromatic rings
          switch (atom->getAtomicNum()) {
            // Carbon
            case 6:
              // CB
              // Aromatic carbon, e.g., in benzene
              atomType = 37;
            break;
          
            // Nitrogen
            case 7:
            if (isAtomNOxide(atom)) {
              // NPOX
              // Pyridinium N-oxide nitrogen
              atomType = 69;
              break;
            }
            if (atom->getTotalDegree() == 3) {
              // NPD+
              // Aromatic nitrogen in pyridinium
              atomType = 58;
              break;
            }
            // NPYD
            // Aromatic nitrogen with sigma lone pair
            atomType = 38;
            break;
          }
        }
      }

      if (!atomType) {

        // Aliphatic heavy atom types

        switch (atom->getAtomicNum()) {
          // Lithium
          case 3:
            if (atom->getDegree() == 0) {
              // LI+
              // Lithium cation
              atomType = 92;
              break;
            }
          break;

          // Carbon
          case 6:
            // 4 neighbors
            if (atom->getTotalDegree() == 4) {
              if (ringInfo->isAtomInRingOfSize(atom->getIdx(), 3)) {
                // CR3R
                // Aliphatic carbon in 3-membered ring
                atomType = 22;
                break;
              }
              if (ringInfo->isAtomInRingOfSize(atom->getIdx(), 4)) {
                // CR4R
                // Aliphatic carbon in 4-membered ring
                atomType = 20;
                break;
              }
              // CR
              // Alkyl carbon
              atomType = 1;
              break;
            }
            // 3 neighbors
            if (atom->getTotalDegree() == 3) {
              unsigned int nN2 = 0;
              unsigned int nN3 = 0;
              unsigned int nO = 0;
              unsigned int nS = 0;
              unsigned int doubleBondedElement = 0;
              boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
              for (; nbrIdx != endNbrs; ++nbrIdx) {
                const Atom *nbrAtom = mol[*nbrIdx].get();
                // find if there is a double-bonded element
                if ((mol.getBondBetweenAtoms(nbrAtom->getIdx(),
                  atom->getIdx()))->getBondType() == Bond::DOUBLE) {
                  doubleBondedElement = nbrAtom->getAtomicNum();
                }
                // count how many terminal oxygen/sulfur atoms
                // are bonded to ipso
                if (nbrAtom->getTotalDegree() == 1) {
                  if (nbrAtom->getAtomicNum() == 8) {
                    ++nO;
                  }
                  else if (nbrAtom->getAtomicNum() == 16) {
                    ++nS;
                  }
                }
                else if (nbrAtom->getAtomicNum() == 7) {
                  // count how many nitrogens with 3 neighbors
                  // are bonded to ipso
                  if (nbrAtom->getTotalDegree() == 3) {
                    ++nN3;
                  }
                  // count how many nitrogens with 2 neighbors
                  // are double-bonded to ipso
                  else if ((nbrAtom->getTotalDegree() == 2)
                    && ((mol.getBondBetweenAtoms(nbrAtom->getIdx(),
                    atom->getIdx()))->getBondType() == Bond::DOUBLE)) {
                    ++nN2;
                  }
                }
              }
              // if there are two or more nitrogens with 3 neighbors each,
              // and there are no nitrogens with two neighbors only,
              // and carbon is double-bonded to nitrogen
              if ((nN3 >= 2) && (!nN2) && (doubleBondedElement == 7)) {
                // CNN+
                // Carbon in +N=C-N: resonance structures
                // CGD+
                // Guanidinium carbon
                atomType = 57;
                break;
              }
              // if there are two terminal oxygen/sulfur atoms
              if ((nO == 2) || (nS == 2)) {
                // CO2M
                // Carbon in carboxylate anion
                // CS2M
                // Carbon in thiocarboxylate anion
                atomType = 41;
                break;
              }
              // if this carbon is in a 4-membered ring and
              // is double-bonded to another carbon
              if (ringInfo->isAtomInRingOfSize(atom->getIdx(), 4)
                && (doubleBondedElement == 6)) {
                // CR4E
                // Olefinic carbon in 4-membered ring
                atomType = 30;
                break;
              }
              // if this carbon is is double-bonded to nitrogen,
              // oxygen, phosphorus or sulfur
              if ((doubleBondedElement == 7) || (doubleBondedElement == 8)
                || (doubleBondedElement == 15) || (doubleBondedElement == 16)) {
                // C=N
                // Imine-atomType carbon
                // CGD
                // Guanidine carbon
                // C=O
                // Generic carbonyl carbon
                // C=OR
                // Ketone or aldehyde carbonyl carbon
                // C=ON
                // Amide carbonyl carbon
                // COO
                // Carboxylic acid or ester carbonyl carbon
                // COON
                // Carbamate carbonyl carbon
                // COOO
                // Carbonic acid or ester carbonyl function
                // C=OS
                // Thioester carbonyl carbon, double bonded to O
                // C=P
                // Carbon doubly bonded to P
                // C=S
                // Thioester carbon, double bonded to S
                // C=SN
                // Thioamide carbon, double bonded to S
                // CSO2
                // Carbon in >C=SO2
                // CS=O
                // Sulfinyl carbon in >C=S=O
                // CSS
                // Thiocarboxylic acid or ester carbon
                atomType = 3;
                break;
              }
              // otherwise it must be generic sp2 carbon
              // C=C
              // Vinylic carbon
              // CSP2
              // Generic sp2 carbon
              atomType = 2;
              break;
            }
            // 2 neighbors
            if (atom->getTotalDegree() == 2) {
              // CSP
              // Acetylenic carbon
              // =C=
              // Allenic carbon
              atomType = 4;
              break;
            }
            // 1 neighbor
            if (atom->getTotalDegree() == 1) {
              // C%-
              // Isonitrile carbon
              atomType = 60;
              break;
            }
          break;

          // Nitrogen
          case 7:
            // if the neighbor is phosphorus or sulfur
            // count the number of terminal oxygens bonded
            // to that phosphorus or sulfur atom
            boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
            for (; nbrIdx != endNbrs; ++nbrIdx) {
              const Atom *nbrAtom = mol[*nbrIdx].get();
              // count how many terminal oxygen atoms
              // are bonded to ipso
              if ((nbrAtom->getAtomicNum() == 8)
                && (nbrAtom->getTotalDegree() == 1)) {
                ++nTermObondedToN;
              }
              if (((atom->getExplicitValence() + atom->getNumImplicitHs()) >= 3)
                && ((nbrAtom->getAtomicNum() == 15) || (nbrAtom->getAtomicNum() == 16))) {
                unsigned int nObondedToSP = 0;
                boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(nbrAtom);
                for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
                  const Atom *nbr2Atom = mol[*nbr2Idx].get();
                  if ((nbr2Atom->getAtomicNum() == 8)
                    && (nbr2Atom->getTotalDegree() == 1)) {
                    ++nObondedToSP;
                  }
                }
                // if there are two or more oxygens, ipso is a sulfonamide nitrogen
                if (!isNSO2orNSO3orNCN) {
                  isNSO2orNSO3orNCN = (nObondedToSP >= 2);
                }
              }
            }
            // 4 neighbors
            if (atom->getTotalDegree() == 4) {
              if (isAtomNOxide(atom)) {
                // N3OX
                // sp3-hybridized N-oxide nitrogen
                atomType = 68;
                break;
              }
              // NR+
              // Quaternary nitrogen
              atomType = 34;
              break;
            }
            // 3 neighbors
            if (atom->getTotalDegree() == 3) {
              // total bond order >= 4
              if ((atom->getExplicitValence() + atom->getNumImplicitHs()) >= 4) {
                bool doubleBondedCN = false;
                boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
                for (; nbrIdx != endNbrs; ++nbrIdx) {
                  const Atom *nbrAtom = mol[*nbrIdx].get();
                  // find if there is a double-bonded nitrogen,
                  // or a carbon which is not bonded to other
                  // nitrogen atoms with 3 neighbors
                  if ((mol.getBondBetweenAtoms(nbrAtom->getIdx(),
                    atom->getIdx()))->getBondType() == Bond::DOUBLE) {
                    doubleBondedCN = ((nbrAtom->getAtomicNum() == 7)
                      || (nbrAtom->getAtomicNum() == 6));
                    if (nbrAtom->getAtomicNum() == 6) {
                      boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(nbrAtom);
                      for (; doubleBondedCN && (nbr2Idx != end2Nbrs); ++nbr2Idx) {
                        const Atom *nbr2Atom = mol[*nbr2Idx].get();
                        if (nbr2Atom->getIdx() == atom->getIdx()) {
                          continue;
                        }
                        doubleBondedCN = (!((nbr2Atom->getAtomicNum() == 7)
                          && (nbr2Atom->getTotalDegree() == 3)));
                      }
                    }
                  }
                }
                // if there is a single terminal oxygen
                if (nTermObondedToN == 1) {
                  // N2OX
                  // sp2-hybridized N-oxide nitrogen
                  atomType = 67;
                  break;
                }
                // if there are two or more terminal oxygens
                if (nTermObondedToN >= 2) {
                  // NO2
                  // Nitrogen in nitro group
                  // NO3
                  // Nitrogen in nitrate group
                  atomType = 45;
                  break;
                }
                // if the carbon bonded to ipso is bonded to 1 nitrogen
                // with 3 neighbors, that nitrogen is ipso (>N+=C)
                // alternatively, if there is no carbon but ipso is
                // double bonded to nitrogen, we have >N+=N
                if (doubleBondedCN) {
                  // N+=C
                  // Iminium nitrogen
                  // N+=N
                  // Positively charged nitrogen doubly bonded to N
                  atomType = 54;
                  break;
                }
              }
              // total bond order >= 3
              if ((atom->getExplicitValence() + atom->getNumImplicitHs()) >= 3) {
                bool isNCOorNCS = false;
                bool isNCNplus = false;
                bool isNGDplus = false;
                bool isNNNorNNC = false;
                bool isNbrC = false;
                bool isNbrBenzeneC = false;
                unsigned int elementDoubleBondedToC = 0;
                unsigned int elementTripleBondedToC = 0;
                unsigned int nN2bondedToC = 0;
                unsigned int nN3bondedToC = 0;
                unsigned int nObondedToC = 0;
                unsigned int nSbondedToC = 0;
                // loop over neighbors
                boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
                for (; nbrIdx != endNbrs; ++nbrIdx) {
                  const Atom *nbrAtom = mol[*nbrIdx].get();
                  // if the neighbor is carbon
                  if (nbrAtom->getAtomicNum() == 6) {
                    isNbrC = true;
                    // check if we have a benzene carbon close to ipso
                    if (nbrAtom->getIsAromatic()
                      && ringInfo->isAtomInRingOfSize(nbrAtom->getIdx(), 6)) {
                      isNbrBenzeneC = true;
                    }
                    nN2bondedToC = 0;
                    nN3bondedToC = 0;
                    nObondedToC = 0;
                    nSbondedToC = 0;
                    unsigned int nFormalCharge = 0;
                    unsigned int nInAromatic6Ring = 0;
                    // loop over carbon neighbors
                    boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(nbrAtom);
                    for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
                      const Atom *nbr2Atom = mol[*nbr2Idx].get();
                      const Bond *bond = mol.getBondBetweenAtoms
                        (nbrAtom->getIdx(), nbr2Atom->getIdx());
                      // check if we have oxygen or sulfur double-bonded to this carbon
                      if ((bond->getBondType() == Bond::DOUBLE)
                        && ((nbr2Atom->getAtomicNum() == 8)
                        || (nbr2Atom->getAtomicNum() == 16))) {
                        isNCOorNCS = true;
                      }
                      // check if there is an atom double-bonded to this carbon,
                      // and if so find which element; if it is carbon or
                      // nitrogen (provided that the latter does not belong to
                      // multiple rings), also an aromatic bond is acceptable
                      if ((bond->getBondType() == Bond::DOUBLE)
                        || ((bond->getBondType() == Bond::AROMATIC)
                        && ((nbr2Atom->getAtomicNum() == 6)
                        || ((nbr2Atom->getAtomicNum() == 7)
                        && (queryIsAtomInNRings(nbr2Atom) == 1))))) {
                        elementDoubleBondedToC = nbr2Atom->getAtomicNum();
                      }
                      // check there is an atom triple-bonded to this carbon,
                      // and if so find which element
                      if (bond->getBondType() == Bond::TRIPLE) {
                        elementTripleBondedToC = nbr2Atom->getAtomicNum();
                      }
                      // if this carbon is bonded to a nitrogen with 3 neighbors
                      if ((nbr2Atom->getAtomicNum() == 7)
                        && (nbr2Atom->getTotalDegree() == 3)) {
                        // count the number of +1 formal charges that we have
                        if (nbr2Atom->getFormalCharge() == 1) {
                          ++nFormalCharge;
                        }
                        if (isAtomInAromaticRingOfSize(nbrAtom, 6)) {
                          ++nInAromatic6Ring;
                        }
                        // count how many oxygens are bonded to this nitrogen
                        // with 3 neighbors
                        unsigned int nObondedToN3 = 0;
                        boost::tie(nbr3Idx, end3Nbrs) = mol.getAtomNeighbors(nbr2Atom);
                        for (; nbr3Idx != end3Nbrs; ++nbr3Idx) {
                          const Atom *nbr3Atom = mol[*nbr3Idx].get();
                          if (nbr3Atom->getAtomicNum() == 8) {
                            ++nObondedToN3;
                          }
                        }
                        // if there are less than 2 oxygens, this is neither
                        // a nitro group nor a nitrate, so increment the counter
                        // of nitrogens with 3 neighbors bonded to this carbon (C-N<)
                        if (nObondedToN3 < 2) {
                          ++nN3bondedToC;
                        }
                      }
                      // if this carbon is bonded to a nitrogen with 2 neighbors
                      // via a double or aromatic bond, increment the counter
                      // of nitrogens with 2 neighbors bonded to this carbon
                      // via a double or aromatic bond (C=N-)
                      if ((nbr2Atom->getAtomicNum() == 7)
                        && (nbr2Atom->getTotalDegree() == 2)
                        && ((bond->getBondType() == Bond::DOUBLE)
                        || (bond->getBondType() == Bond::AROMATIC))) {
                        ++nN2bondedToC;
                      }
                      // if this carbon is bonded to an aromatic atom
                      if (nbr2Atom->getIsAromatic()) {
                        // if it is oxygen, increment the counter of
                        // aromatic oxygen atoms bonded to this carbon
                        if (nbr2Atom->getAtomicNum() == 8) {
                          ++nObondedToC;
                        }
                        // if it is sulfur, increment the counter of
                        // aromatic sulfur atoms bonded to this carbon
                        if (nbr2Atom->getAtomicNum() == 16) {
                          ++nSbondedToC;
                        }
                      }
                    }
                    // if nitrogen is bonded to this carbon via a double or aromatic bond
                    if (elementDoubleBondedToC == 7) {
                      // if 2 nitrogens with 3 neighbors and no nitrogens with 2 neighbors
                      // are bonded to this carbon, and we have a formal charge,
                      // but not a 6-membered aromatic ring, and the carbon atom
                      // is not sp3, then this is an amidinium nitrogen (>N-C=N+<)
                      if ((nN3bondedToC == 2) && (!nN2bondedToC)
                        && nFormalCharge && (!nInAromatic6Ring)
                        && (nbrAtom->getTotalDegree() < 4)) {
                        isNCNplus = true;
                      }
                      // if 3 nitrogens with 3 neighbors are bonded
                      // to this carbon, then this is a guanidinium nitrogen ((>N-)2-C=N+<)
                      if (nN3bondedToC == 3) {
                        isNGDplus =true;
                      }
                    }
                  }
                  // if the neighbor is nitrogen
                  if (nbrAtom->getAtomicNum() == 7) {
                    unsigned int nNbondedToN = 0;
                    unsigned int nObondedToN = 0;
                    unsigned int nSbondedToN = 0;
                    // loop over nitrogen neighbors
                    boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(nbrAtom);
                    for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
                      const Atom *nbr2Atom = mol[*nbr2Idx].get();
                      const Bond *bond = mol.getBondBetweenAtoms
                        (nbrAtom->getIdx(), nbr2Atom->getIdx());
                      // if the bond to nitrogen is double
                      if (bond->getBondType() == Bond::DOUBLE) {
                        // if the neighbor is carbon (N=N-C)
                        if (nbr2Atom->getAtomicNum() == 6) {
                          // loop over carbon neighbors
                          boost::tie(nbr3Idx, end3Nbrs) = mol.getAtomNeighbors(nbr2Atom);
                          for (; nbr3Idx != end3Nbrs; ++nbr3Idx) {
                            const Atom *nbr3Atom = mol[*nbr3Idx].get();
                            // if the nitrogen neighbor to ipso is met, move on
                            if (nbr3Atom->getIdx() == nbrAtom->getIdx()) {
                              continue;
                            }
                            // count how many nitrogen, oxygen, sulfur atoms
                            // are bonded to this carbon
                            switch (nbr3Atom->getAtomicNum()) {
                              case 7:
                              ++nNbondedToN;
                              break;
                              case 8:
                              ++nObondedToN;
                              break;
                              case 16:
                              ++nSbondedToN;
                              break;
                            }
                          }
                          // if there are no more nitrogens, no oxygen, no sulfur bonded
                          // to carbon, and the latter is not a benzene carbon
                          // then it is N=N-C
                          if ((!nObondedToN) && (!nSbondedToN)
                            && (!nNbondedToN) && (!isNbrBenzeneC)) {
                            isNNNorNNC = true;
                          }
                        }
                        // if the neighbor is nitrogen (N=N-N) and ipso is not bonded
                        // to benzene carbon then it is N=N-N
                        if ((nbr2Atom->getAtomicNum() == 7) && (!isNbrBenzeneC)) {
                          isNNNorNNC = true;
                        }
                      }
                    }
                  }
                }
                // if ipso nitrogen is bonded to carbon
                if (isNbrC) {
                  // if neighbor carbon is triple-bonded to N, then ipso is N-C%N
                  if (elementTripleBondedToC == 7) {
                    isNSO2orNSO3orNCN = true;
                  }
                  // if neighbor carbon is amidinium
                  if (isNCNplus) {
                    // NCN+
                    // Either nitrogen in N+=C-N
                    atomType = 55;
                    break;
                  }
                  // if neighbor carbon is guanidinium
                  if (isNGDplus) {
                    // NGD+
                    // Guanidinium nitrogen
                    atomType = 56;
                    break;
                  }
                  // if neighbor carbon is not bonded to oxygen or sulfur
                  // and is not cyano, there two possibilities:
                  // 1) ipso nitrogen is bonded to benzene carbon while no oxygen
                  //    or sulfur are bonded to the latter: ipso is aniline nitrogen
                  // 2) ipso nitrogen is bonded to a carbon which is double-bonded to
                  //    carbon, nitrogen or phosphorus, or triple-bonded to carbon
                  if (((!isNCOorNCS) && (!isNSO2orNSO3orNCN))
                    && (((!nObondedToC) && (!nSbondedToC) && isNbrBenzeneC)
                    || ((elementDoubleBondedToC == 6) || (elementDoubleBondedToC == 7)
                    || (elementDoubleBondedToC == 15) || (elementTripleBondedToC == 6)))) {
                    // NC=C
                    // Enamine or aniline nitrogen, deloc. lp
                    // NC=N
                    // Nitrogen in N-C=N with deloc. lp
                    // NC=P
                    // Nitrogen in N-C=P with deloc. lp
                    // NC%C
                    // Nitrogen attached to C-C triple bond
                    atomType = 40;
                    break;
                  }
                }
                // if ipso is not sulfonamide while it is either amide/thioamide
                // or >N-N=N-/>N-N=C<
                if ((!isNSO2orNSO3orNCN) && (isNCOorNCS || isNNNorNNC)) {
                  // NC=O
                  // Amide nitrogen
                  // NC=S
                  // Thioamide nitrogen
                  // NN=C
                  // Nitrogen in N-N=C moiety with deloc. lp
                  // NN=N
                  // Nitrogen in N-N=N moiety with deloc. lp
                  atomType = 10;
                  break;
                }
              }
            }
            // 2 neighbors
            if (atom->getTotalDegree() == 2) {
              // total bond order = 4
              if ((atom->getExplicitValence() + atom->getNumImplicitHs()) == 4) {
                // loop over neighbors
                bool isIsonitrile = false;
                boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
                for (; (!isIsonitrile) && (nbrIdx != endNbrs); ++nbrIdx) {
                  const Atom *nbrAtom = mol[*nbrIdx].get();
                  // if neighbor is triple-bonded
                  isIsonitrile = ((mol.getBondBetweenAtoms(atom->getIdx(),
                    nbrAtom->getIdx()))->getBondType() == Bond::TRIPLE);
                }
                if (isIsonitrile) {
                  // NR%
                  // Isonitrile nitrogen
                  atomType = 61;
                  break;
                }
                // =N=
                // Central nitrogen in C=N=N or N=N=N
                atomType = 53;
                break;
              }
              // total bond order = 3
              if ((atom->getExplicitValence() + atom->getNumImplicitHs()) == 3) {
                // loop over neighbors
                bool isNitroso = false;
                bool isImineOrAzo = false;
                boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
                for (; nbrIdx != endNbrs; ++nbrIdx) {
                  const Atom *nbrAtom = mol[*nbrIdx].get();
                  // if the neighbor is double bonded (-N=)
                  if ((mol.getBondBetweenAtoms(atom->getIdx(),
                    nbrAtom->getIdx()))->getBondType() == Bond::DOUBLE) {
                    // if it is terminal oxygen (-N=O)
                    isNitroso = ((nbrAtom->getAtomicNum() == 8)
                      && (nTermObondedToN == 1));
                    // if it is carbon or nitrogen (-N=N-, -N=C<),
                    // ipso is imine or azo
                    isImineOrAzo = ((nbrAtom->getAtomicNum() == 6)
                      || (nbrAtom->getAtomicNum() == 7));
                  }
                }
                if (isNitroso && (!isImineOrAzo)) {
                  // N=O
                  // Nitrogen in nitroso group
                  atomType = 46;
                  break;
                }
                if (isImineOrAzo) {
                  // N=C
                  // Imine nitrogen
                  // N=N
                  // Azo-group nitrogen
                  atomType = 9;
                  break;
                }
              }
              // total bond order >= 2
              if ((atom->getExplicitValence() + atom->getNumImplicitHs()) >= 2) {
                // loop over neighbors
                bool isNSO = false;
                boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
                for (; (!isNSO) && (nbrIdx != endNbrs); ++nbrIdx) {
                  const Atom *nbrAtom = mol[*nbrIdx].get();
                  // if the neighbor is sulfur bonded to a single terminal oxygen
                  if (nbrAtom->getAtomicNum() == 16) {
                    // loop over neighbors and count how many
                    // terminal oxygens are bonded to sulfur
                    unsigned int nTermObondedToS = 0;
                    boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(nbrAtom);
                    for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
                      const Atom *nbr2Atom = mol[*nbr2Idx].get();
                      if ((nbr2Atom->getAtomicNum() == 8)
                        && (nbr2Atom->getTotalDegree() == 1)) {
                        ++nTermObondedToS;
                      }
                    }
                    isNSO = (nTermObondedToS == 1);
                  }
                }
                if (isNSO) {
                  // NSO
                  // Divalent nitrogen replacing monovalent O in SO2 group
                  atomType = 48;
                  break;
                }
                if (!isNSO2orNSO3orNCN) {
                  // If it is not sulfonamide deprotonated nitrogen,
                  // it is anionic nitrogen (>N::-)
                  // NM
                  // Anionic divalent nitrogen
                  atomType = 62;
                  break;
                }
              }
            }
            // if it is sulfonamide (3 neighbors) or cyano (2 neighbors)
            if (isNSO2orNSO3orNCN) {
              // NSO2
              // Sulfonamide nitrogen
              // NSO3
              // Sulfonamide nitrogen
              // NC%N
              // Nitrogen attached to cyano group
              atomType = 43;
              break;
            }
            // 1 neighbor
            if (atom->getTotalDegree() == 1) {
              bool isNSP = false;
              bool isNAZT = false;
              boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
              for (; (!isNSP) && (!isNAZT) && (nbrIdx != endNbrs); ++nbrIdx) {
                const Atom *nbrAtom = mol[*nbrIdx].get();
                // if ipso is triple-bonded to its only neighbor
                isNSP = ((mol.getBondBetweenAtoms(atom->getIdx(),
                  nbrAtom->getIdx()))->getBondType() == Bond::TRIPLE);
                // ipso is bonded to a nitrogen atom with 2 neighbors
                if ((nbrAtom->getAtomicNum() == 7)
                  && (nbrAtom->getTotalDegree() == 2)) {
                  // loop over nitrogen neighbors
                  boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(nbrAtom);
                  for (; (!isNAZT) && (nbr2Idx != end2Nbrs); ++nbr2Idx) {
                    const Atom *nbr2Atom = mol[*nbr2Idx].get();
                    // if another nitrogen with 2 neighbors, or a carbon
                    // with 3 neighbors is found, ipso is NAZT
                    isNAZT = (((nbr2Atom->getAtomicNum() == 7)
                      && (nbr2Atom->getTotalDegree() == 2))
                      || ((nbr2Atom->getAtomicNum() == 6)
                      && (nbr2Atom->getTotalDegree() == 3)));
                  }
                }
              }
              if (isNSP) {
                // NSP
                // Triply bonded nitrogen
                atomType = 42;
                break;
              }
              if (isNAZT) {
                // NAZT
                // Terminal nitrogen in azido or diazo group
                atomType = 47;
                break;
              }
            }
            // if nothing else was found
            // NR
            // Amine nitrogen
            atomType = 8;
          break;

          // Oxygen
          case 8:
            // 3 neighbors
            if (atom->getTotalDegree() == 3) {
              // O+
              // Oxonium oxygen
              atomType = 49;
              break;
            }
            // 2 neighbors
            if (atom->getTotalDegree() == 2) {
              if ((atom->getExplicitValence() + atom->getNumImplicitHs()) == 3) {
                // O=+
                // Oxenium oxygen
                atomType = 51;
                break;
              }
              // count how many hydrogens are bound to ipso
              unsigned int nHbondedToO = 0;
              boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
              for (; nbrIdx != endNbrs; ++nbrIdx) {
                const Atom *nbrAtom = mol[*nbrIdx].get();
                if (nbrAtom->getAtomicNum() == 1) {
                  ++nHbondedToO;
                }
              }
              if ((nHbondedToO + atom->getNumImplicitHs()) == 2) {
                // OH2
                // Oxygen in water
                atomType = 70;
                break;
              }
              // otherwise, ipso must be one of the following
              // OC=O
              // Carboxylic acid or ester oxygen
              // OC=C
              // Enolic or phenolic oxygen
              // OC=N
              // Oxygen in -O-C=N- moiety
              // OC=S
              // Divalent oxygen in thioacid or ester
              // ONO2
              // Divalent nitrate "ether" oxygen
              // ON=O
              // Divalent nitrate "ether" oxygen
              // OSO3
              // Divalent oxygen in sulfate group
              // OSO2
              // Divalent oxygen in sulfite group
              // OSO
              // One of two divalent oxygens attached to sulfur
              // OS=O
              // Divalent oxygen in R(RO)S=O
              // -OS
              // Other divalent oxygen attached to sulfur
              // OPO3
              // Divalent oxygen in phosphate group
              // OPO2
              // Divalent oxygen in phosphite group
              // OPO
              // Divalent oxygen, one of two oxygens attached to P
              // -OP
              // Other divalent oxygen attached to phosphorus
              atomType = 6;
              break;
            }
            // 1 neighbor
            if (atom->getDegree() <= 1) {
              unsigned int nNbondedToCorNorS = 0;
              unsigned int nObondedToCorNorS = 0;
              unsigned int nSbondedToCorNorS = 0;
              bool isOxideOBondedToH = atom->getNumExplicitHs() + atom->getNumImplicitHs();
              bool isCarboxylateO = false;
              bool isCarbonylO = false;
              bool isOxideOBondedToC = false;
              bool isNitrosoO = false;
              bool isOxideOBondedToN = false;
              bool isNOxideO = false;
              bool isNitroO = false;
              bool isThioSulfinateO = false;
              bool isSulfateO = false;
              bool isSulfoxideO = false;
              bool isPhosphateOrPerchlorateO = false;
              // loop over neighbors
              boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
              for (; (nbrIdx != endNbrs) && (!isOxideOBondedToC) && (!isOxideOBondedToN)
                && (!isOxideOBondedToH) && (!isCarboxylateO) && (!isNitroO) && (!isNOxideO)
                && (!isThioSulfinateO) && (!isSulfateO) && (!isPhosphateOrPerchlorateO)
                && (!isCarbonylO) && (!isNitrosoO) && (!isSulfoxideO); ++nbrIdx) {
                const Atom *nbrAtom = mol[*nbrIdx].get();
                const Bond *bond = mol.getBondBetweenAtoms
                  (atom->getIdx(), nbrAtom->getIdx());
                // if the neighbor is carbon, nitrogen or sulfur
                if ((nbrAtom->getAtomicNum() == 6)
                  || (nbrAtom->getAtomicNum() == 7)
                  || (nbrAtom->getAtomicNum() == 16)) {
                  // count how many terminal oxygen/sulfur atoms
                  // or secondary nitrogens
                  // are bonded to the carbon or nitrogen neighbor of ipso
                  boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(nbrAtom);
                  for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
                    const Atom *nbr2Atom = mol[*nbr2Idx].get();
                    if ((nbr2Atom->getAtomicNum() == 7)
                      && (nbr2Atom->getTotalDegree() == 2)) {
                      ++nNbondedToCorNorS;
                    }
                    if ((nbr2Atom->getAtomicNum() == 8)
                      && (nbr2Atom->getTotalDegree() == 1)) {
                      ++nObondedToCorNorS;
                    }
                    if ((nbr2Atom->getAtomicNum() == 16)
                      && (nbr2Atom->getTotalDegree() == 1)) {
                      ++nSbondedToCorNorS;
                    }
                  }
                }
                // if ipso neighbor is hydrogen
                isOxideOBondedToH = (nbrAtom->getAtomicNum() == 1);
                
                // if ipso neighbor is carbon
                if (nbrAtom->getAtomicNum() == 6) {
                  // if carbon neighbor is bonded to 2 oxygens,
                  // ipso is carboxylate oxygen
                  isCarboxylateO = (nObondedToCorNorS == 2);
                  // if ipso oxygen is bonded to carbon
                  // via a double bond, ipso is carbonyl oxygen
                  isCarbonylO = (bond->getBondType() == Bond::DOUBLE);
                  // if ipso oxygen is bonded to carbon via a
                  // single bond, and there are no other bonded oxygens,
                  // ipso is oxide oxygen
                  isOxideOBondedToC = ((bond->getBondType() == Bond::SINGLE)
                    && (nObondedToCorNorS == 1));
                }
                
                // if ipso neighbor is nitrogen
                if (nbrAtom->getAtomicNum() == 7) {
                  // if ipso oxygen is bonded to nitrogen
                  // via a double bond, ipso is nitroso oxygen
                  isNitrosoO = (bond->getBondType() == Bond::DOUBLE);
                  // if ipso oxygen is bonded to nitrogen via a single bond
                  // and there are no other bonded oxygens
                  if ((bond->getBondType() == Bond::SINGLE) && (nObondedToCorNorS == 1)) {
                    // if nitrogen has 2 neighbors or, if the neighbors are 3,
                    // the total bond order on nitrogen is 3, ipso is oxide oxygen
                    isOxideOBondedToN = ((nbrAtom->getTotalDegree() == 2)
                      || ((nbrAtom->getExplicitValence() + nbrAtom->getNumImplicitHs()) == 3));
                    // if the total bond order on nitrogen is 4, ipso is N-oxide oxygen
                    isNOxideO = ((nbrAtom->getExplicitValence() + nbrAtom->getNumImplicitHs()) == 4);
                  }
                  // if ipso oxygen is bonded to nitrogen which is bonded
                  // to multiple oxygens, ipso is nitro/nitrate oxygen
                  isNitroO = (nObondedToCorNorS >= 2);
                }
                
                // if ipso neighbor is sulfur
                if (nbrAtom->getAtomicNum() == 16) {
                  // if ipso oxygen is bonded to sulfur and
                  // the latter is bonded to another sulfur,
                  // ipso is thiosulfinate oxygen
                  isThioSulfinateO = (nSbondedToCorNorS == 1);
                  // if ipso oxygen is bonded to sulfur via a single
                  // bond or, if the bond is double, there are multiple
                  // oxygen/nitrogen atoms bonded to that sulfur,
                  // ipso is sulfate oxygen
                  isSulfateO = ((bond->getBondType() == Bond::SINGLE)
                    || ((bond->getBondType() == Bond::DOUBLE)
                    && ((nObondedToCorNorS + nNbondedToCorNorS) > 1)));
                  // if ipso oxygen is bonded to sulfur via a double
                  // bond and the sum of oxygen/nitrogen atoms bonded
                  // to that sulfur is 1, ipso is sulfoxide oxygen
                  isSulfoxideO = ((bond->getBondType() == Bond::DOUBLE)
                    && ((nObondedToCorNorS + nNbondedToCorNorS) == 1));
                }
                
                // if ipso neighbor is phosphorus or chlorine
                isPhosphateOrPerchlorateO = ((nbrAtom->getAtomicNum() == 15)
                  || (nbrAtom->getAtomicNum() == 17));
              }
              if (isOxideOBondedToC || isOxideOBondedToN || isOxideOBondedToH) {
                // OM
                // Oxide oxygen on sp3 carbon
                // OM2
                // Oxide oxygen on sp2 carbon
                // OM
                // Oxide oxygen on sp3 nitrogen (not in original MMFF.I Table III)
                // OM2
                // Oxide oxygen on sp2 nitrogen (not in original MMFF.I Table III)
                atomType = 35;
                break;
              }
              if (isCarboxylateO || isNitroO || isNOxideO || isThioSulfinateO
                || isSulfateO || isPhosphateOrPerchlorateO) {
                // O2CM
                // Oxygen in carboxylate group
                // ONX
                // Oxygen in N-oxides
                // O2N
                // Oxygen in nitro group
                // O2NO
                // Nitro-group oxygen in nitrate
                // O3N
                // Nitrate anion oxygen
                // OSMS
                // Terminal oxygen in thiosulfinate anion
                // O-S
                // Single terminal O on tetracoordinate sulfur
                // O2S
                // One of 2 terminal O's on sulfur
                // O3S
                // One of 3 terminal O's on sulfur
                // O4S
                // Terminal O in sulfate anion
                // OP
                // Oxygen in phosphine oxide
                // O2P
                // One of 2 terminal O's on P
                // O3P
                // One of 3 terminal O's on P
                // O4P
                // One of 4 terminal O's on P
                // O4Cl
                // Oxygen in perchlorate anion
                atomType = 32;
                break;
              }
              if (isCarbonylO || isNitrosoO || isSulfoxideO) {
                // O=C
                // Generic carbonyl oxygen
                // O=CN
                // Carbonyl oxygen in amides
                // O=CR
                // Carbonyl oxygen in aldehydes and ketones
                // O=CO
                // Carbonyl oxygen in acids and esters
                // O=N
                // Nitroso oxygen
                // O=S
                // Doubly bonded sulfoxide oxygen
                atomType = 7;
                break;
              }
            }
          break;

          // Fluorine
          case 9:
            // 1 neighbor
            if (atom->getDegree() == 1) {
              // F
              // Fluorine
              atomType = 11;
              break;
            }
            if (atom->getDegree() == 0) {
              // F-
              // Fluoride anion
              atomType = 89;
              break;
            }
          break;

          // Sodium
          case 11:
            if (atom->getDegree() == 0) {
              // NA+
              // Sodium cation
              atomType = 93;
              break;
            }
          break;
          
          // Magnesium
          case 12:
            if (atom->getDegree() == 0) {
              // MG+2
              // Dipositive magnesium cation
              atomType = 99;
              break;
            }
          break;
          
          // Silicon
          case 14:
            // SI
            // Silicon
            atomType = 19;
          break;
          
          // Phosphorus
          case 15:
            if (atom->getTotalDegree() == 4) {
              // PO4
              // Phosphate group phosphorus
              // PO3
              // Phosphorus with 3 attached oxygens
              // PO2
              // Phosphorus with 2 attached oxygens
              // PO
              // Phosphine oxide phosphorus
              // PTET
              // General tetracoordinate phosphorus
              atomType = 25;
              break;
            }
            if (atom->getTotalDegree() == 3) {
              // P
              // Phosphorus in phosphines
              atomType = 26;
              break;
            }
            if (atom->getTotalDegree() == 2) {
              // -P=C
              // Phosphorus doubly bonded to C
              atomType = 75;
              break;
            }
          break;
          
          // Sulfur
          case 16:
            // 3  or 4 neighbors
            if ((atom->getTotalDegree() == 3) || (atom->getTotalDegree() == 4)) {
              unsigned int nOorNbondedToS = 0;
              unsigned int nSbondedToS = 0;
              bool isCDoubleBondedToS = false;
              // loop over neighbors
              boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
              for (; nbrIdx != endNbrs; ++nbrIdx) {
                const Atom *nbrAtom = mol[*nbrIdx].get();
                // check if ipso sulfur is double-bonded to carbon
                if ((nbrAtom->getAtomicNum() == 6)
                  && ((mol.getBondBetweenAtoms(atom->getIdx(),
                  nbrAtom->getIdx()))->getBondType() == Bond::DOUBLE)) {
                  isCDoubleBondedToS = true;
                }
                // if the neighbor is terminal oxygen/sulfur
                // or secondary nitrogen, increment the respective counter
                if (((nbrAtom->getDegree() == 1) && (nbrAtom->getAtomicNum() == 8))
                  || ((nbrAtom->getTotalDegree() == 2) && (nbrAtom->getAtomicNum() == 7))) {
                  ++nOorNbondedToS;
                }
                if ((nbrAtom->getDegree() == 1) && (nbrAtom->getAtomicNum() == 16)) {
                  ++nSbondedToS;
                }
              }
              // if ipso sulfur has 3 neighbors and is bonded to
              // two atoms of oxygen/nitrogen and double-bonded
              // to carbon, or if it has 4 neighbors
              if (((atom->getTotalDegree() == 3) && (nOorNbondedToS == 2)
                && (isCDoubleBondedToS)) || (atom->getTotalDegree() == 4)) {
                // =SO2
                // Sulfone sulfur, doubly bonded to carbon
                atomType = 18;
                break;
              }
              // if ipso sulfur is bonded to both oxygen/nitrogen and sulfur
              if ((nOorNbondedToS && nSbondedToS)
                || ((nOorNbondedToS == 2) && (!isCDoubleBondedToS))) {
                // SSOM
                // Tricoordinate sulfur in anionic thiosulfinate group
                atomType = 73;
                break;
              }
              // otherwise ipso sulfur is double bonded to oxygen or nitrogen
              // S=O
              // Sulfoxide sulfur
              // >S=N
              // Tricoordinate sulfur doubly bonded to N
              atomType = 17;
              break;
            }
            // 2 neighbors
            if (atom->getTotalDegree() == 2) {
              // loop over neighbors
              bool isODoubleBondedToS = false;
              boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
              for (; nbrIdx != endNbrs; ++nbrIdx) {
                const Atom *nbrAtom = mol[*nbrIdx].get();
                // check if ipso sulfur is double-bonded to oxygen
                if ((nbrAtom->getAtomicNum() == 8)
                  && ((mol.getBondBetweenAtoms(atom->getIdx(),
                  nbrAtom->getIdx()))->getBondType() == Bond::DOUBLE)) {
                  isODoubleBondedToS = true;
                }
              }
              // if ipso sulfur is double-bonded to oxygen
              if (isODoubleBondedToS) {
                // =S=O
                // Sulfinyl sulfur, e.g., in C=S=O
                atomType = 74;
                break;
              }
              // otherwise it is a thiol, sulfide or disulfide
              // S
              // Thiol, sulfide, or disulfide sulfur
              atomType = 15;
              break;
            }
            // 1 neighbor
            if (atom->getDegree() == 1) {
              unsigned int nTermSbondedToNbr = 0;
              bool isCDoubleBondedToS = false;
              // find the neighbor and count how many terminal sulfur
              // atoms are there, including ipso
              boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
              for (; nbrIdx != endNbrs; ++nbrIdx) {
                const Atom *nbrAtom = mol[*nbrIdx].get();
                boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(nbrAtom);
                for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
                  const Atom *nbr2Atom = mol[*nbr2Idx].get();
                  if ((nbr2Atom->getAtomicNum() == 16)
                    && (nbr2Atom->getTotalDegree() == 1)) {
                    ++nTermSbondedToNbr;
                  }
                }
                // check if ipso sulfur is double-bonded to carbon
                if ((nbrAtom->getAtomicNum() == 6)
                  && ((mol.getBondBetweenAtoms(atom->getIdx(),
                  nbrAtom->getIdx()))->getBondType() == Bond::DOUBLE)) {
                  isCDoubleBondedToS = true;
                }
              }
              // if ipso sulfur is double bonded to carbon and the latter
              // is not bonded to other terminal sulfur atoms, then it is
              // not a dithiocarboxylate, but a thioketone, etc.
              if (isCDoubleBondedToS && (nTermSbondedToNbr != 2)) {
                // S=C
                // Sulfur doubly bonded to carbon
                atomType = 16;
                break;
              }
              // otherwise ipso must be one of these
              // S-P
              // Terminal sulfur bonded to P
              // SM
              // Anionic terminal sulfur
              // SSMO
              // Terminal sulfur in thiosulfinate group
              atomType = 72;
              break;
            }
          break;
          
          // Chlorine
          case 17:
            // 4 neighbors
            if (atom->getTotalDegree() == 4) {
              // loop over neighbors and count the number
              // of bonded oxygens
              unsigned int nObondedToCl = 0;
              boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
              for (; nbrIdx != endNbrs; ++nbrIdx) {
                const Atom *nbrAtom = mol[*nbrIdx].get();
                if (nbrAtom->getAtomicNum() == 8) {
                  ++nObondedToCl;
                }
              }
              // if there are 4 oxygens
              if (nObondedToCl == 4) {
                // CLO4
                // Perchlorate anione chlorine
                atomType = 77;
                break;
              }
            }
            // 1 neighbor
            if (atom->getTotalDegree() == 1) {
              // Cl
              // Chlorine
              atomType = 12;
              break;
            }
            // 0 neighbors
            if (atom->getDegree() == 0) {
              // Cl-
              // Chloride anion
              atomType = 90;
              break;
            }
          break;
          
          // Potassium
          case 19:
            if (atom->getDegree() == 0) {
              // K+
              // Potassium cation
              atomType = 94;
              break;
            }
          break;
          
          // Calcium
          case 20:
            if (atom->getDegree() == 0) {
              // CA+2
              // Dipositive calcium cation
              atomType = 96;
              break;
            }
          break;
          
          // Iron
          case 26:
            if (atom->getDegree() == 0) {
              if (atom->getFormalCharge() == 2) {
                // FE+2
                // Dipositive iron cation
                atomType = 87;
                break;
              }
              if (atom->getFormalCharge() == 3) {
                // FE+3
                // Tripositive iron cation
                atomType = 88;
                break;
              }
            }
          break;
          
          // Copper
          case 29:
            if (atom->getDegree() == 0) {
              if (atom->getFormalCharge() == 1) {
                // CU+1
                // Monopositive copper cation
                atomType = 97;
                break;
              }
              if (atom->getFormalCharge() == 2) {
                // CU+2
                // Dipositive copper cation
                atomType = 98;
                break;
              }
            }
          break;
          
          // Zinc
          case 30:
            if (atom->getDegree() == 0) {
              // ZN+2
              // Dipositive zinc cation
              atomType = 95;
              break;
            }
          break;
          
          // Bromine
          case 35:
            if (atom->getDegree() == 1) {
              // Br
              // Bromine
              atomType = 13;
              break;
            }
            if (atom->getDegree() == 0) {
              // BR-
              // Bromide anion
              atomType = 91;
              break;
            }
          break;
          
          // Iodine
          case 53:
            if (atom->getDegree() == 1) {
              // I
              // Iodine
              atomType = 14;
              break;
            }
          break;
        }
      }
      d_MMFFAtomPropertiesPtrVect[atom->getIdx()]->mmffAtomType = atomType;
      if (!atomType) {
        d_valid = false;
      }
    }


    // finds the MMFF atomType for a hydrogen atom
    void MMFFMolProperties::setMMFFHydrogenType(const Atom *atom)
    {
      unsigned int atomType=0;
      bool isHOCCorHOCN = false;
      bool isHOCO = false;
      bool isHOP = false;
      bool isHOS = false;
      ROMol mol = atom->getOwningMol();
      ROMol::ADJ_ITER nbrIdx;
      ROMol::ADJ_ITER endNbrs;
      ROMol::ADJ_ITER nbr2Idx;
      ROMol::ADJ_ITER end2Nbrs;
      ROMol::ADJ_ITER nbr3Idx;
      ROMol::ADJ_ITER end3Nbrs;
      
      
      // loop over neighbors (actually there can be only one)
      boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
      for (; nbrIdx != endNbrs; ++nbrIdx) {
        const Atom *nbrAtom = mol[*nbrIdx].get();
        switch (nbrAtom->getAtomicNum()) {
          // carbon, silicon
          case 6:
          case 14:
            // HC
            // Hydrogen attached to carbon
            // HSI
            // Hydrogen attached to silicon
            atomType = 5;
          break;
          
          // nitrogen
          case 7:
            switch (this->getMMFFAtomType(nbrAtom->getIdx())) {
              case 8:
                // HNR
                // Generic hydrogen on sp3 nitrogen, e.g. in amines
                // H3N
                // Hydrogen in ammonia
              case 39:
                // HPYL
                // Hydrogen on nitrogen in pyrrole
              case 62:
                // HNR
                // Generic hydrogen on sp3 nitrogen, e.g. in amines
              case 67:
              case 68:
                // HNOX
                // Hydrogen on N in a N-oxide
                atomType = 23;
              break;
              
              case 34:
                // NR+
                // Quaternary nitrogen
              case 54:
                // N+=C
                // Iminium nitrogen
                // N+=N
                // Positively charged nitrogen doubly bonded to N
              case 55:
                // HNN+
                // Hydrogen on amidinium nitrogen
              case 56:
                // HGD+
                // Hydrogen on guanidinium nitrogen
              case 58:
                // NPD+
                // Aromatic nitrogen in pyridinium
              case 81:
                // HIM+
                // Hydrogen on imidazolium nitrogen
                atomType = 36;
              break;
              
              case 9:
                // HN=N
                // Hydrogen on azo nitrogen
                // HN=C
                // Hydrogen on imine nitrogen
                atomType = 27;
              break;
              
              default:
                // HNCC
                // Hydrogen on enamine nitrogen
                // HNCN
                // Hydrogen in H-N-C=N moiety
                // HNCO
                // Hydrogen on amide nitrogen
                // HNCS
                // Hydrogen on thioamide nitrogen
                // HNNC
                // Hydrogen in H-N-N=C moiety
                // HNNN
                // Hydrogen in H-N-N=N moiety
                // HNSO
                // Hydrogen on NSO, NSO2, or NSO3 nitrogen
                // HNC%
                // Hydrogen on N triply bonded to C
                // HSP2
                // Generic hydrogen on sp2 nitrogen
                atomType = 28;
              break;
            }
          break;
          
          // oxygen
          case 8:
            switch (this->getMMFFAtomType(nbrAtom->getIdx())) {
              case 49:
                // HO+
                // Hydrogen on oxonium oxygen
                atomType = 50;
              break;
              
              case 51:
                // HO=+
                // Hydrogen on oxenium oxygen
                atomType = 52;
              break;
              
              case 70:
                // HOH
                // Hydroxyl hydrogen in water
                atomType = 31;
              break;
              
              case 6:
                // for hydrogen bonded to atomType 6 oxygen we need to distinguish
                // among acidic hydrogens belonging to carboxylic/phospho acids,
                // enolic/phenolic/hydroxamic hydrogens and hydrogens whose oxygen
                // partner is bonded to sulfur. If none of these is found
                // it is either an alcohol or a generic hydroxyl hydrogen
                // loop over oxygen neighbors
                boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(nbrAtom);
                for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
                  const Atom *nbr2Atom = mol[*nbr2Idx].get();
                  // if the neighbor of oxygen is carbon, loop over the carbon neighbors
                  if (nbr2Atom->getAtomicNum() == 6) {
                    boost::tie(nbr3Idx, end3Nbrs) = mol.getAtomNeighbors(nbr2Atom);
                    for (; nbr3Idx != end3Nbrs; ++nbr3Idx) {
                      const Atom *nbr3Atom = mol[*nbr3Idx].get();
                      const Bond *bond = mol.getBondBetweenAtoms(nbr2Atom->getIdx(),
                        nbr3Atom->getIdx());
                      // if the starting oxygen is met, move on
                      if (nbr3Atom->getIdx() == nbrAtom->getIdx()) {
                        continue;
                      }
                      // if the carbon neighbor is another carbon or nitrogen
                      // bonded via a double or aromatic bond, ipso is HOCC/HOCN
                      if (((nbr3Atom->getAtomicNum() == 6)
                        || (nbr3Atom->getAtomicNum() == 7))
                        && ((bond->getBondType() == Bond::DOUBLE)
                        || (bond->getBondType() == Bond::AROMATIC))) {
                        isHOCCorHOCN = true;
                      }
                      // if the carbon neighbor is an oxygen bonded
                      // via a double bond, ipso is HOCO
                      if ((nbr3Atom->getAtomicNum() == 8)
                        && (bond->getBondType() == Bond::DOUBLE)) {
                        isHOCO = true;
                      }
                    }
                  }
                  // if the neighbor of oxygen is phosphorus, ipso is HOCO
                  if (nbr2Atom->getAtomicNum() == 15) {
                    isHOP = true;
                  }
                  // if the neighbor of oxygen is sulfur, ipso is HOS
                  if (nbr2Atom->getAtomicNum() == 16) {
                    isHOS = true;
                  }
                }
                if (isHOCO || isHOP) {
                  // HOCO
                  // Hydroxyl hydrogen in carboxylic acids
                  atomType = 24;
                  break;
                }
                if (isHOCCorHOCN) {
                  // HOCC
                  // Enolic or phenolic hydroxyl hydrogen
                  // HOCN
                  // Hydroxyl hydrogen in HO-C=N moiety
                  atomType = 29;
                  break;
                }
                if (isHOS) {
                  // HOS
                  // Hydrogen on oxygen attached to sulfur
                  atomType = 33;
                  break;
                }
                
              default:
                // HO
                // Generic hydroxyl hydrogen
                // HOR
                // Hydroxyl hydrogen in alcohols
                atomType = 21;
              break;
            }
          break;
          
          // phosphorus and sulfur
          case 15:
          case 16:
            // HP
            // Hydrogen attached to phosphorus
            // HS
            // Hydrogen attached to sulfur
            // HS=N
            // Hydrogen attached to >S= sulfur doubly bonded to N
            atomType = 71;
           break;
        }
      }
      d_MMFFAtomPropertiesPtrVect[atom->getIdx()]->mmffAtomType = atomType;
      if (!atomType) {
        d_valid = false;
      }
    }


    // sanitizes molecule according to MMFF requirements
    // returns MolOps::SANITIZE_NONE on success, the flag
    // which caused trouble in case of failure
    unsigned int sanitizeMMFFMol(RWMol &mol)
    {
      unsigned int error = 0;
      
      try { 
        MolOps::sanitizeMol(mol, error,
                            (unsigned int)(MolOps::SANITIZE_CLEANUP
                                           | MolOps::SANITIZE_PROPERTIES
                                           | MolOps::SANITIZE_SYMMRINGS
                                           | MolOps::SANITIZE_KEKULIZE
                                           | MolOps::SANITIZE_FINDRADICALS
                                           | MolOps::SANITIZE_SETCONJUGATION
                                           | MolOps::SANITIZE_SETHYBRIDIZATION
                                           | MolOps::SANITIZE_CLEANUPCHIRALITY
                                           | MolOps::SANITIZE_ADJUSTHS));
        if (!(mol.hasProp("_MMFFSanitized"))) {
          mol.setProp("_MMFFSanitized",1,true);
        }
      } catch (MolSanitizeException &e){
      
      }
      
      return error;
    }
    
    
    // constructs a MMFFMolProperties object for ROMol mol filled
    // with MMFF atom types, formal and partial charges
    // in case atom types are missing, d_valid is set to false,
    // charges are set to 0.0 and the force-field is unusable
    MMFFMolProperties::MMFFMolProperties(ROMol &mol, 
      std::string mmffVariant, boost::uint8_t verbosity,
      std::ostream &oStream) :
      d_valid(true),
      d_mmffs(mmffVariant == "MMFF94s" ? true : false),
      d_bondTerm(true),
      d_angleTerm(true),
      d_stretchBendTerm(true),
      d_oopTerm(true),
      d_torsionTerm(true),
      d_vdWTerm(true),
      d_eleTerm(true),
      d_dielConst(1.0),
      d_dielModel(CONSTANT),
      d_verbosity(verbosity),
      d_oStream(&oStream),
      d_MMFFAtomPropertiesPtrVect(mol.getNumAtoms()) {
      ROMol::AtomIterator it;
      if (!(mol.hasProp("_MMFFSanitized"))) {
        bool isAromaticSet = false;
        for (it = mol.beginAtoms(); (!isAromaticSet) && (it != mol.endAtoms()); ++it) {
          isAromaticSet = (*it)->getIsAromatic();
        }
        if (isAromaticSet) {
          MolOps::Kekulize((RWMol &)mol, true);
        }
        mol.setProp("_MMFFSanitized",1,true);
      }      
      for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        d_MMFFAtomPropertiesPtrVect[i]
          = MMFFAtomPropertiesPtr(new MMFFAtomProperties());
      }
      unsigned int idx;
      boost::uint8_t atomType = 1;
      
      setMMFFAromaticity((RWMol &)mol);
      for (it = mol.beginAtoms(); it != mol.endAtoms(); ++it) {
        if ((*it)->getAtomicNum() != 1) {
          this->setMMFFHeavyAtomType(*it);
        }
      }
      for (it = mol.beginAtoms(); atomType && (it != mol.endAtoms()); ++it) {
        if ((*it)->getAtomicNum() == 1) {
          this->setMMFFHydrogenType(*it);
        }
      }
      if (this->isValid()) {
        this->computeMMFFCharges(mol);
      }
      if (verbosity == MMFF_VERBOSITY_HIGH) {
        oStream <<
          "\n"
          "A T O M   T Y P E S   A N D   C H A R G E S\n\n"
          "          ATOM    FORMAL   PARTIAL\n"
          " ATOM     TYPE    CHARGE    CHARGE\n"
          "-----------------------------------"
          << std::endl;
        for (idx = 0; idx < mol.getNumAtoms(); ++idx) {
          oStream
            << std::left << std::setw(2) << mol.getAtomWithIdx(idx)->getSymbol()
            << std::left << " #" << std::setw(5) << idx + 1
            << std::right << std::setw(5)
            << (unsigned int)(this->getMMFFAtomType(idx))
            << std::right << std::setw(10) << std::fixed << std::setprecision(3)
            << this->getMMFFFormalCharge(idx)
            << std::right << std::setw(10) << this->getMMFFPartialCharge(idx)
            << std::endl;
        }
        if (!(this->isValid())) {
          oStream << "\nMissing atom types - charges were not computed"
            << std::endl;
        }
      }
    }
    

    // returns the MMFF angle type of the angle formed
    // by atoms with indexes idx1, idx2, idx3
    unsigned int MMFFMolProperties::getMMFFAngleType
      (const ROMol &mol, const unsigned int idx1,
      const unsigned int idx2, const unsigned int idx3)
    {
      PRECONDITION(this->isValid(), "missing atom types - invalid force-field");
      
      // ftp://ftp.wiley.com/public/journals/jcc/suppmat/17/553/MMFF-III_AppendixA.html
      //
      // AT[IJK]    Structural significance
      //--------------------------------------------------------------------------
      //  0		      The angle i-j-k is a "normal" bond angle
      //  1 		    Either bond i-j or bond j-k has a bond type of 1
      //  2		      Bonds i-j and j-k each have bond types of 1; the sum is 2.
      //  3		      The angle occurs in a three-membered ring
      //  4		      The angle occurs in a four-membered ring
      //  5		      Is in a three-membered ring and the sum of the bond types is 1
      //  6		      Is in a three-membered ring and the sum of the bond types is 2
      //  7		      Is in a four-membered ring and the sum of the bond types is 1
      //  8		      Is in a four-membered ring and the sum of the bond types is 2
      
      unsigned int bondTypeSum =
        this->getMMFFBondType(mol.getBondBetweenAtoms(idx1, idx2))
        + this->getMMFFBondType(mol.getBondBetweenAtoms(idx2, idx3));
      unsigned int angleType = bondTypeSum;
      
      unsigned int size = isAngleInRingOfSize3or4(mol, idx1, idx2, idx3);
      if (size) {
        angleType = size;
        if (bondTypeSum) {
          angleType += (bondTypeSum + size - 2);
        }
      }
      
      return angleType;
    }
    

    // returns the MMFF bond type of the bond
    unsigned int MMFFMolProperties::getMMFFBondType(const Bond *bond)
    {
      PRECONDITION(this->isValid(), "missing atom types - invalid force-field");

      MMFFPropCollection *mmffProp = MMFFPropCollection::getMMFFProp();
      const ForceFields::MMFF::MMFFProp *mmffPropAtom1 =
        (*mmffProp)(this->getMMFFAtomType(bond->getBeginAtomIdx()));
      const ForceFields::MMFF::MMFFProp *mmffPropAtom2 =
        (*mmffProp)(this->getMMFFAtomType(bond->getEndAtomIdx()));

      // return 1 if the bond is single and the properties for this 
      // single bond match either those of sbmb or aromatic bonds
      // for this atom pair, 0 if they don't
      return (unsigned int)(((bond->getBondType() == Bond::SINGLE)
        && ((mmffPropAtom1->sbmb && mmffPropAtom2->sbmb)
        || (mmffPropAtom1->arom && mmffPropAtom2->arom))) ? 1 : 0);
    }
    
    
    // given the angle type and the two bond types of the bond
    // which compose the angle, it returns the MMFF stretch-bend
    // type of the angle
    unsigned int getMMFFStretchBendType(const unsigned int angleType,
      const unsigned int bondType1, const unsigned int bondType2)
    {
      unsigned int stretchBendType = 0;
      
      switch (angleType) {
        case 1:
          stretchBendType = ((bondType1 || (bondType1 == bondType2)) ? 1 : 2);
        break;
        
        case 2:
          stretchBendType = 3;
        break;

        case 4:
          stretchBendType = 4;
        break;
        
        case 3:
          stretchBendType = 5;
        break;
        
        case 5:
          stretchBendType = ((bondType1 || (bondType1 == bondType2)) ? 6 : 7);
        break;
        
        case 6:
          stretchBendType = 8;
        break;
        
        case 7:
          stretchBendType = ((bondType1 || (bondType1 == bondType2)) ? 9 : 10);
        break;
        
        case 8:
          stretchBendType = 11;
        break;
      }
      
      return stretchBendType;
    }

    
    // given a dihedral angle formed by 4 atoms with indexes
    // idx1, idx2, idx3, idx4, it returns a std::pair whose first element
    // is the principal torsion type, and the second is the secondary
    // torsion type, to be used only if parameters could not be found
    // (empirically found - this is not mentioned either in MMFF.IV
    // nor in MMFF.V)
    const std::pair<unsigned int, unsigned int> MMFFMolProperties::getMMFFTorsionType
      (const ROMol &mol, const unsigned int idx1, const unsigned int idx2,
      const unsigned int idx3, const unsigned int idx4)
    {
      PRECONDITION(this->isValid(), "missing atom types - invalid force-field");

      const Bond *bondJK = mol.getBondBetweenAtoms(idx2, idx3);
      unsigned int bondTypeIJ = this->getMMFFBondType
        (mol.getBondBetweenAtoms(idx1, idx2));
      unsigned int bondTypeJK = this->getMMFFBondType(bondJK);
      unsigned int bondTypeKL = this->getMMFFBondType
        (mol.getBondBetweenAtoms(idx3, idx4));
      unsigned int torsionType = bondTypeJK;
      unsigned int secondTorsionType = 0;
      
      // according to MMFF.IV page 609 the condition should be as simple as
      // if ((bondTypeJK == 0) && ((bondTypeIJ == 1) || (bondTypeKL == 1))) {
      // but CYGUAN01 fails the test, so the following condition was
      // empirically determined to be the correct one
      if ((bondTypeJK == 0) && (bondJK->getBondType() == Bond::SINGLE)
        && ((bondTypeIJ == 1) || (bondTypeKL == 1))) {
        torsionType = 2;
      }
      unsigned int size = isTorsionInRingOfSize4or5(mol, idx1, idx2, idx3, idx4);
      // the additional check on the existence of a bond between I and K or J and L
      // is to avoid assigning torsionType 4 to those torsions in a 4-membered ring
      // constituted by the fusion of two 3-membered rings, even though it would
      // be harmless for the energy calculation since parameters for 
      // 4,22,22,22,22 and 0,22,22,22,22 are identical
      if ((size == 4) && (!(mol.getBondBetweenAtoms(idx1, idx3)
        || mol.getBondBetweenAtoms(idx2, idx4)))) {
        secondTorsionType = torsionType;
        torsionType = 4;
      }
      else if ((size == 5)
        && ((this->getMMFFAtomType(idx1) == 1)|| (this->getMMFFAtomType(idx2) == 1)
        || (this->getMMFFAtomType(idx3) == 1) || (this->getMMFFAtomType(idx4) == 1))) {
        secondTorsionType = torsionType;
        torsionType = 5;
      }
      
      return std::make_pair(torsionType, secondTorsionType);
    }


    // empirical rule to compute bond stretching parameters if
    // tabulated parameters could not be found. The returned
    // pointer to a MMFFBond object must be freed by the caller
    const ForceFields::MMFF::MMFFBond *
      MMFFMolProperties::getMMFFBondStretchEmpiricalRuleParams
      (const ROMol &mol, const Bond *bond)
    {
      PRECONDITION(this->isValid(), "missing atom types - invalid force-field");

      const MMFFBond *mmffBndkParams;
      const MMFFHerschbachLaurie *mmffHerschbachLaurieParams;
      const MMFFProp *mmffAtomPropParams[2];
      const MMFFCovRadPauEle *mmffAtomCovRadPauEleParams[2];
      MMFFBndkCollection *mmffBndk = MMFFBndkCollection::getMMFFBndk();
      MMFFHerschbachLaurieCollection *mmffHerschbachLaurie =
        MMFFHerschbachLaurieCollection::getMMFFHerschbachLaurie();
      MMFFCovRadPauEleCollection *mmffCovRadPauEle =
        MMFFCovRadPauEleCollection::getMMFFCovRadPauEle();
      MMFFPropCollection *mmffProp = MMFFPropCollection::getMMFFProp();

      unsigned int atomicNum1 = bond->getBeginAtom()->getAtomicNum();
      unsigned int atomicNum2 = bond->getEndAtom()->getAtomicNum();
      mmffBndkParams = (*mmffBndk)(atomicNum1, atomicNum2);
      mmffAtomCovRadPauEleParams[0] = (*mmffCovRadPauEle)(atomicNum1);
      mmffAtomCovRadPauEleParams[1] = (*mmffCovRadPauEle)(atomicNum2);
      mmffAtomPropParams[0] = (*mmffProp)(this->getMMFFAtomType(bond->getBeginAtomIdx()));
      mmffAtomPropParams[1] = (*mmffProp)(this->getMMFFAtomType(bond->getEndAtomIdx()));

      PRECONDITION(mmffAtomCovRadPauEleParams[0],
        "covalent radius/Pauling electronegativity parameters for atom 1 not found");
      PRECONDITION(mmffAtomCovRadPauEleParams[1],
        "covalent radius/Pauling electronegativity parameters for atom 2 not found");
      PRECONDITION(mmffAtomPropParams[0],
        "property parameters for atom 1 not found");
      PRECONDITION(mmffAtomPropParams[1],
        "property parameters for atom 2 not found");

      ForceFields::MMFF::MMFFBond *mmffBondParams = new ForceFields::MMFF::MMFFBond();
      const double c = (((atomicNum1 == 1) || (atomicNum2 == 1)) ? 0.050 : 0.085);
      const double n = 1.4;
      #if 0
      const double delta = 0.008;
      #endif
      #if 1
      const double delta = 0.0;
      #endif
      double r0_i[2];

      // MMFF.V, page 625
      for (unsigned int i = 0; i < 2; ++i) {
        r0_i[i] = mmffAtomCovRadPauEleParams[i]->r0;
        // the part of the empirical rule concerning H
        // parameters appears not to be used - tests are
        // passed only in its absence, hence it is
        // currently excluded
        #if 0
        switch (mmffAtomPropParams[i]->mltb) {
          case 1:
          case 2:
            H_i[i] = 2;
          break;
          
          case 3:
            H_i[i] = 1;
          
          default:
            H_i[i] = 3;
        }
        #endif
      }
      // also the part of the empirical rule concerning BO
      // parameters appears not to be used - tests are
      // passed only in its absence, hence it is
      // currently excluded
      #if 0
      unsigned int BO_ij = (unsigned int)(bond->getBondTypeAsDouble());
      if ((mmffAtomPropParams[0]->mltb == 1)
        && (mmffAtomPropParams[1]->mltb == 1)) {
        BO_ij = 4;
      }
      if (((mmffAtomPropParams[0]->mltb == 1)
        && (mmffAtomPropParams[1]->mltb == 2))
        || ((mmffAtomPropParams[0]->mltb == 2)
        && (mmffAtomPropParams[1]->mltb == 1))) {
        BO_ij = 5;
      }
      if (areAtomsInSameAromaticRing(mol,
        bond->getBeginAtomIdx(), bond->getEndAtomIdx())) {
        BO_ij = (((mmffAtomPropParams[0]->pilp == 0)
          && (mmffAtomPropParams[1]->pilp == 0)) ? 4 : 5);
      }
      if (BO_ij == 1) {
        for (unsigned int i = 0; i < 2; ++i) {
          std::cout << "H" << i << "=" << H_i[i] << std::endl;
          switch (H_i[i]) {
            case 1:
              r0_i[i] -= 0.08;
            break;
            
            case 2:
              r0_i[i] -= 0.03;
            break;
          }
        }
      }
      else {
        double dec = 0.0;
        switch (BO_ij) {
          case 5:
            dec = 0.04;
          break;
          
          case 4:
            dec = 0.075;
          break;
          
          case 3:
            dec = 0.17;
          break;
          
          case 2:
            dec = 0.10;
          break;
        }
        r0_i[0] -= dec;
        r0_i[1] -= dec;
      }
      std::cout << "BO_ij=" << BO_ij << std::endl;
      #endif
      // equation (18) - MMFF.V, page 625
      mmffBondParams->r0 = (r0_i[0] + r0_i[1] - c * pow(fabs
        (mmffAtomCovRadPauEleParams[0]->chi
        - mmffAtomCovRadPauEleParams[1]->chi), n) - delta);
      if (mmffBndkParams) {
        // equation (19) - MMFF.V, page 625
        double coeff = mmffBndkParams->r0 / mmffBondParams->r0;
        double coeff2 = coeff * coeff;
        double coeff6 = coeff2 * coeff2 * coeff2;
        mmffBondParams->kb = mmffBndkParams->kb * coeff6;
      }
      else {
        // MMFF.V, page 627
        // Herschbach-Laurie version of Badger's rule
        // J. Chem. Phys. 35, 458 (1961); http://dx.doi.org/10.1063/1.1731952
        // equation (8), page 5
        mmffHerschbachLaurieParams = (*mmffHerschbachLaurie)
          (getPeriodicTableRowHL(atomicNum1), getPeriodicTableRowHL(atomicNum2));
        mmffBondParams->kb = pow(10.0, -(mmffBondParams->r0
          - mmffHerschbachLaurieParams->a_ij) / mmffHerschbachLaurieParams->d_ij);
      }
      
      return (const ForceFields::MMFF::MMFFBond *)mmffBondParams;
    }


    // empirical rule to compute angle bending parameters if
    // tabulated parameters could not be found. The returned
    // pointer to a MMFFAngle object must be freed by the caller
    const ForceFields::MMFF::MMFFAngle *getMMFFAngleBendEmpiricalRuleParams
      (const ROMol &mol, const ForceFields::MMFF::MMFFAngle *oldMMFFAngleParams,
      const ForceFields::MMFF::MMFFProp *mmffPropParamsCentralAtom,
      const ForceFields::MMFF::MMFFBond *mmffBondParams1,
      const ForceFields::MMFF::MMFFBond *mmffBondParams2,
      unsigned int idx1, unsigned int idx2, unsigned int idx3)
    {
      int atomicNum[3];
      atomicNum[0] = mol.getAtomWithIdx(idx1)->getAtomicNum();
      atomicNum[1] = mol.getAtomWithIdx(idx2)->getAtomicNum();
      atomicNum[2] = mol.getAtomWithIdx(idx3)->getAtomicNum();
      ForceFields::MMFF::MMFFAngle *mmffAngleParams = new ForceFields::MMFF::MMFFAngle();
      unsigned int ringSize = isAngleInRingOfSize3or4(mol, idx1, idx2, idx3);
      if (!oldMMFFAngleParams) {
        // angle rest value empirical rule
        mmffAngleParams->theta0 = 120.0;
        switch (mmffPropParamsCentralAtom->crd) {
          case 4:
            // if the central atom has crd = 4
            mmffAngleParams->theta0 = 109.45;
          break;
          
          case 2:
            // if the central atom is oxygen
            if (atomicNum[1] == 6) {
              mmffAngleParams->theta0 = 105.0;
            }
            // if the central atom is linear
            else if (mmffPropParamsCentralAtom->linh == 1) {
              mmffAngleParams->theta0 = 180.0;
            }
          break;
          
          case 3:
            if ((mmffPropParamsCentralAtom->val == 3)
              && (mmffPropParamsCentralAtom->mltb == 0)) {
              // if the central atom is nitrogen
              if (atomicNum[1] == 5) {
                mmffAngleParams->theta0 = 107.0;
              }
              else {
                mmffAngleParams->theta0 = 92.0;
              }
            }
          break;
        }
        if (ringSize == 3) {
          mmffAngleParams->theta0 = 60.0;
        }
        else if (ringSize == 4) {
          mmffAngleParams->theta0 = 90.0;
        }
      }
      else {
        mmffAngleParams->theta0 = oldMMFFAngleParams->theta0;
      }
      // angle force constant empirical rule
      double Z[3] = { 0.0, 0.0, 0.0 };
      double C[3] = { 0.0, 0.0, 0.0 };
      double beta = 1.75;
      for (unsigned int i = 0; i < 3; ++i) {
        // Table VI - MMFF.V, page 628
        switch (atomicNum[i]) {
          // Hydrogen
          case 1:
            Z[i] = 1.395;
          break;
          
          // Carbon
          case 6:
            Z[i] = 2.494;
            C[i] = 1.016;
          break;
          
          // Nitrogen
          case 7:
            Z[i] = 2.711;
            C[i] = 1.113;
          break;
          
          // Oxygen
          case 8:
            Z[i] = 3.045;
            C[i] = 1.337;
          break;
          
          // Fluorine
          case 9:
            Z[i] = 2.847;
          break;
          
          // Silicon
          case 14:
            Z[i] = 2.350;
            C[i] = 0.811;
          break;
          
          // Phosphorus
          case 15:
            Z[i] = 2.350;
            C[i] = 1.068;
          break;
          
          // Sulfur
          case 16:
            Z[i] = 2.980;
            C[i] = 1.249;
          break;
          
          // Chlorine
          case 17:
            Z[i] = 2.909;
            C[i] = 1.078;
          break;
          
          // Bromine
          case 35:
            Z[i] = 3.017;
          break;
          
          // Iodine
          case 53:
            Z[i] = 3.086;
          break;
        }
      }
      double r0_ij = mmffBondParams1->r0;
      double r0_jk = mmffBondParams2->r0;
      double D = (r0_ij - r0_jk) * (r0_ij - r0_jk) / ((r0_ij + r0_jk) * (r0_ij + r0_jk));
      double theta0_rad = DEG2RAD * mmffAngleParams->theta0;
      if (ringSize == 4) {
        beta *= 0.85;
      }
      else if (ringSize == 3) {
        beta *= 0.05;
      }
      // equation (20) - MMFF.V, page 628
      mmffAngleParams->ka = beta * Z[0] * C[1] * Z[2]
        / ((r0_ij + r0_jk) * theta0_rad * theta0_rad * exp(2.0 * D));
      
      return (const ForceFields::MMFF::MMFFAngle *)mmffAngleParams;
    }
    

    // empirical rule to compute torsional parameters if
    // tabulated parameters could not be found
    // the indexes of the two central atoms J and K
    // idx2 and idx3 must be supplied. The returned pointer
    // to a MMFFTor object must be freed by the caller
    const ForceFields::MMFF::MMFFTor *
      MMFFMolProperties::getMMFFTorsionEmpiricalRuleParams
      (const ROMol &mol, unsigned int idx2, unsigned int idx3)
    {
      PRECONDITION(this->isValid(), "missing atom types - invalid force-field");

      MMFFPropCollection *mmffProp = MMFFPropCollection::getMMFFProp();
      MMFFAromCollection *mmffArom = MMFFAromCollection::getMMFFArom();
      ForceFields::MMFF::MMFFTor *mmffTorParams = new ForceFields::MMFF::MMFFTor();
      unsigned int jAtomType = this->getMMFFAtomType(idx2);
      unsigned int kAtomType = this->getMMFFAtomType(idx3);
      const MMFFProp *jMMFFProp = (*mmffProp)(jAtomType);
      const MMFFProp *kMMFFProp = (*mmffProp)(kAtomType);
      const Bond *bond = mol.getBondBetweenAtoms(idx2, idx3);
      double U[2] = { 0.0, 0.0 };
      double V[2] = { 0.0, 0.0 };
      double W[2] = { 0.0, 0.0 };
      double beta = 0.0;
      double pi_jk = 0.0;
      const double N_jk = (double)((jMMFFProp->crd - 1)
        * (kMMFFProp->crd - 1));
      int atomicNum[2] = { mol.getAtomWithIdx(idx2)->getAtomicNum(),
        mol.getAtomWithIdx(idx3)->getAtomicNum() };
      
      for (unsigned int i = 0; i < 2; ++i) {
        switch (atomicNum[i]) {
          // carbon
          case 6:
            U[i] = 2.0;
            V[i] = 2.12;
          break;
          
          // nitrogen
          case 7:
            U[i] = 2.0;
            V[i] = 1.5;
          break;
          
          // oxygen
          case 8:
            U[i] = 2.0;
            V[i] = 0.2;
            W[i] = 2.0;
          break;
          
          // silicon
          case 14:
            U[i] = 1.25;
            V[i] = 1.22;
          break;
          
          // phosphorus
          case 15:
            U[i] = 1.25;
            V[i] = 2.40;
          break;
          
          // sulfur
          case 16:
            U[i] = 1.25;
            V[i] = 0.49;
            W[i] = 8.0;
          break;
        }
      }
      
      // rule (a)
      if (jMMFFProp->linh || kMMFFProp->linh) {
        mmffTorParams->V1 = 0.0;
        mmffTorParams->V2 = 0.0;
        mmffTorParams->V3 = 0.0;
      }
      
      // rule (b)
      else if (mmffArom->isMMFFAromatic(jAtomType)
        && mmffArom->isMMFFAromatic(kAtomType)
        && bond->getIsAromatic()) {
        beta = ((((jMMFFProp->val == 3) && (kMMFFProp->val == 4))
          || ((jMMFFProp->val == 4) && (kMMFFProp->val == 3))) ? 3.0 : 6.0);
        pi_jk = (((jMMFFProp->pilp == 0)
          && (kMMFFProp->pilp == 0)) ? 0.5 : 0.3);
        mmffTorParams->V2 = beta * pi_jk * sqrt(U[0] * U[1]);
      }
      
      // rule (c)
      else if (bond->getBondType() == Bond::DOUBLE) {
        beta = 6.0;
        pi_jk =(((jMMFFProp->mltb == 2)
          && (kMMFFProp->mltb == 2)) ? 1.0 : 0.4);
        mmffTorParams->V2 = beta * pi_jk * sqrt(U[0] * U[1]);
      }
      
      // rule (d)
      else if ((jMMFFProp->crd == 4) && (kMMFFProp->crd == 4)) {
        mmffTorParams->V3 = sqrt(V[0] * V[1]) / N_jk;
      }
      
      // rule (e)
      else if ((jMMFFProp->crd == 4) && (kMMFFProp->crd != 4)) {
        if (((kMMFFProp->crd == 3) && (((kMMFFProp->val == 4)
          || (kMMFFProp->val == 34)) || kMMFFProp->mltb))
          || ((kMMFFProp->crd == 2)
          && ((kMMFFProp->val == 3) || kMMFFProp->mltb))) {
          mmffTorParams->V1 = 0.0;
          mmffTorParams->V2 = 0.0;
          mmffTorParams->V3 = 0.0;
        }
        else {
          mmffTorParams->V3 = sqrt(V[0] * V[1]) / N_jk;
        }
      }

      // rule (f)
      else if ((kMMFFProp->crd == 4) && (jMMFFProp->crd != 4)) {
        if (((jMMFFProp->crd == 3) && (((jMMFFProp->val == 4)
          || (jMMFFProp->val == 34)) || jMMFFProp->mltb))
          || ((jMMFFProp->crd == 2)
          && ((jMMFFProp->val == 3) || jMMFFProp->mltb))) {
          mmffTorParams->V1 = 0.0;
          mmffTorParams->V2 = 0.0;
          mmffTorParams->V3 = 0.0;
        }
        else {
          mmffTorParams->V3 = sqrt(V[0] * V[1]) / N_jk;
        }
      }

      // rule (g)
      else if (((bond->getBondType() == Bond::SINGLE)
        && jMMFFProp->mltb && kMMFFProp->mltb)
        || (jMMFFProp->mltb && kMMFFProp->pilp)
        || (jMMFFProp->pilp && kMMFFProp->mltb)) {
        // case (1)
        if (jMMFFProp->pilp && kMMFFProp->pilp) {
          mmffTorParams->V1 = 0.0;
          mmffTorParams->V2 = 0.0;
          mmffTorParams->V3 = 0.0;
        }
        // case (2)
        else if (jMMFFProp->pilp && kMMFFProp->mltb) {
          beta = 6.0;
          if (jMMFFProp->mltb == 1) {
            pi_jk = 0.5;
          }
          else if ((getPeriodicTableRow(atomicNum[0]) == 2)
            && (getPeriodicTableRow(atomicNum[1]) == 2)) {
            pi_jk = 0.3;
          }
          else if ((getPeriodicTableRow(atomicNum[0]) != 2)
            || (getPeriodicTableRow(atomicNum[1]) != 2)) {
            pi_jk = 0.15;
          }
          mmffTorParams->V2 = beta * pi_jk * sqrt(U[0] * U[1]);
        }
        // case (3)
        else if (kMMFFProp->pilp && jMMFFProp->mltb) {
          beta = 6.0;
          if (kMMFFProp->mltb == 1) {
            pi_jk = 0.5;
          }
          else if ((getPeriodicTableRow(atomicNum[0]) == 2)
            && (getPeriodicTableRow(atomicNum[1]) == 2)) {
            pi_jk = 0.3;
          }
          else if ((getPeriodicTableRow(atomicNum[0]) != 2)
            || (getPeriodicTableRow(atomicNum[1]) != 2)) {
            pi_jk = 0.15;
          }
          mmffTorParams->V2 = beta * pi_jk * sqrt(U[0] * U[1]);
        }
        // case (4)
        else if (((jMMFFProp->mltb == 1) || (kMMFFProp->mltb == 1))
          && ((atomicNum[0] != 6) || (atomicNum[1] != 6))) {
          beta = 6.0;
          pi_jk = 0.4;
          mmffTorParams->V2 = beta * pi_jk * sqrt(U[0] * U[1]);
        }
        // case (5)
        else {
          beta = 6.0;
          pi_jk = 0.15;
          mmffTorParams->V2 = beta * pi_jk * sqrt(U[0] * U[1]);
        }
      }
      
      // rule (h)
      else {
        if (((atomicNum[0] == 8) || (atomicNum[0] == 16))
          && ((atomicNum[1] == 8) || (atomicNum[1] == 16))) {
          mmffTorParams->V2 = -sqrt(W[0] * W[1]);
        }
        else {
          mmffTorParams->V3 = sqrt(V[0] * V[1]) / N_jk;
        }
      }

      return (const MMFFTor *)mmffTorParams;
    }


    // populates the MMFFMolProperties object with MMFF
    // formal and partial charges
    void MMFFMolProperties::computeMMFFCharges(const ROMol &mol)
    {
      PRECONDITION(this->isValid(), "missing atom types - invalid force-field");

      unsigned int idx;
      unsigned int i;
      unsigned int j;
      unsigned int atomType;
      unsigned int nbrAtomType;
      unsigned int nConj = 0;
      unsigned int old_nConj = 0;
      std::pair<int, double> bci;
      double pChg = 0.0;
      double fChg = 0.0;
      boost::dynamic_bitset<> conjNBitVect(mol.getNumAtoms());
      VECT_INT_VECT atomRings = mol.getRingInfo()->atomRings();
      ROMol::ADJ_ITER nbrIdx;
      ROMol::ADJ_ITER endNbrs;
      ROMol::ADJ_ITER nbr2Idx;
      ROMol::ADJ_ITER end2Nbrs;
      MMFFPropCollection *mmffProp = MMFFPropCollection::getMMFFProp();
      MMFFPBCICollection *mmffPBCI = MMFFPBCICollection::getMMFFPBCI();
      MMFFChgCollection *mmffChg = MMFFChgCollection::getMMFFChg();


      // We need to set formal charges upfront
      for (idx = 0; idx < mol.getNumAtoms(); ++idx) {
        const Atom *atom = mol.getAtomWithIdx(idx);
        atomType = this->getMMFFAtomType(idx);
        fChg = 0.0;
        switch (atomType) {
          // special cases
          case 32:
            // O2CM
            // Oxygen in carboxylate group
          case 72:
            // SM
            // Anionic terminal sulfur
            // loop over neighbors
            boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
            for (; nbrIdx != endNbrs; ++nbrIdx) {
              const Atom *nbrAtom = mol[*nbrIdx].get();
              nbrAtomType = this->getMMFFAtomType(nbrAtom->getIdx());
              // loop over neighbors of the neighbor
              // count how many terminal oxygen/sulfur atoms
              // or secondary nitrogens
              // are bonded to the neighbor of ipso
              int nSecNbondedToNbr = 0;
              int nTermOSbondedToNbr = 0;
              boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(nbrAtom);
              for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
                const Atom *nbr2Atom = mol[*nbr2Idx].get();
                // if it's nitrogen with 2 neighbors and it is not aromatic,
                // increment the counter of secondary nitrogens
                if ((nbr2Atom->getAtomicNum() == 7) && (nbr2Atom->getDegree() == 2)
                  && (!(nbr2Atom->getIsAromatic()))) {
                    ++nSecNbondedToNbr;
                }
                // if it's terminal oxygen/sulfur,
                // increment the terminal oxygen/sulfur counter
                if (((nbr2Atom->getAtomicNum() == 8) || (nbr2Atom->getAtomicNum() == 16))
                  && (nbr2Atom->getDegree() == 1)) {
                  ++nTermOSbondedToNbr;
                }
              }
              // in case its sulfur with two terminal oxygen/sulfur atoms and one secondary
              // nitrogen, this is a deprotonated sulfonamide, so we should not consider
              // nitrogen as a replacement for oxygen/sulfur in a sulfone
              if ((nbrAtom->getAtomicNum() == 16)
                && (nTermOSbondedToNbr == 2) && (nSecNbondedToNbr == 1)) {
                nSecNbondedToNbr = 0;
              }
              // if the neighbor is carbon
              if ((nbrAtom->getAtomicNum() == 6) && nTermOSbondedToNbr) {
                // O2CM
                // Oxygen in (thio)carboxylate group: charge is shared
                // across 2 oxygens/sulfur atoms in (thio)carboxylate,
                // 3 oxygen/sulfur atoms in (thio)carbonate
                // SM
                // Anionic terminal sulfur: charge is localized
                fChg = ((nTermOSbondedToNbr == 1) ? -1.0
                  : -((double)(nTermOSbondedToNbr - 1) / (double)nTermOSbondedToNbr));
                break;
              }
              // if the neighbor is NO2 or NO3
              if ((nbrAtomType == 45) && (nTermOSbondedToNbr == 3)) {
                // O3N
                // Nitrate anion oxygen
                fChg = -1.0 / 3.0;
                break;
              }
              // if the neighbor is PO2, PO3, PO4
              if ((nbrAtomType == 25) && nTermOSbondedToNbr) {
                // OP
                // Oxygen in phosphine oxide
                // O2P
                // One of 2 terminal O's on P
                // O3P
                // One of 3 terminal O's on P
                // O4P
                // One of 4 terminal O's on P
                fChg = ((nTermOSbondedToNbr == 1) ? 0.0
                  : -((double)(nTermOSbondedToNbr - 1) / (double)nTermOSbondedToNbr));
                break;
              }
              // if the neighbor is SO2, SO2N, SO3, SO4, SO2M, SSOM
              if ((nbrAtomType == 18) && nTermOSbondedToNbr) {
                // SO2
                // Sulfone sulfur
                // SO2N
                // Sulfonamide sulfur
                // SO3
                // Sulfonate group sulfur
                // SO4
                // Sulfate group sulfur
                // SNO
                // Sulfur in nitrogen analog of a sulfone
                fChg = (((nSecNbondedToNbr + nTermOSbondedToNbr) == 2) ? 0.0
                  : -((double)((nSecNbondedToNbr + nTermOSbondedToNbr) - 2)
                  / (double)nTermOSbondedToNbr));
                break;
              }
              if ((nbrAtomType == 73)  && nTermOSbondedToNbr) {
                // SO2M
                // Sulfur in anionic sulfinate group
                // SSOM
                // Tricoordinate sulfur in anionic thiosulfinate group
                fChg = ((nTermOSbondedToNbr == 1) ? 0.0
                  : -((double)(nTermOSbondedToNbr - 1) / (double)nTermOSbondedToNbr));
                break;
              }
              if ((nbrAtomType == 77) && nTermOSbondedToNbr) {
                // O4Cl
                // Oxygen in perchlorate anion
                fChg = -(1.0 / (double)nTermOSbondedToNbr);
                break;
              }
            }
          break;
          
          case 76:
            // N5M
            // Nitrogen in 5-ring aromatic anion
            // we don't need to bother about the neighbors with N5M
            for (i = 0; i < atomRings.size(); ++i) {
              if ((std::find(atomRings[i].begin(), atomRings[i].end(),
                idx) != atomRings[i].end())) {
                break;
              }
            }
            // find how many nitrogens with atom type 76 we have
            // and share the formal charge accordingly
            if (i < atomRings.size()) {
              unsigned int nNitrogensIn5Ring = 0;
              for (j = 0; j < atomRings[i].size(); ++j) {
                if (this->getMMFFAtomType(atomRings[i][j]) == 76) {
                  ++nNitrogensIn5Ring;
                }
              }
              if (nNitrogensIn5Ring) {
                fChg = -(1.0 / (double)nNitrogensIn5Ring);
              }
            }
          break;
            
          case 55:
          case 56:
          case 81:
            // NIM+
            // Aromatic nitrogen in imidazolium
            // N5A+
            // Positive nitrogen in 5-ring alpha position
            // N5B+
            // Positive nitrogen in 5-ring beta position
            // N5+
            // Positive nitrogen in other 5-ring position
            // we need to loop over all molecule atoms
            // and find all those nitrogens with atom type
            // 81, 55 or 56, check whether they are conjugated
            // with ipso and keep on looping until no more
            // conjugated atoms can be found. Finally, we divide
            // the total formal charge that was found on the
            // conjugated system by the number of conjugated nitrogens
            // of types 81, 55 or 56 that were found.
            // This is not strictly what is described
            // in the MMFF papers, but it is the only way to get an
            // integer total formal charge, which makes sense to me
            // probably such conjugated systems are anyway out of the
            // scope of MMFF, but this is an attempt to correctly
            // deal with them somehow
            fChg = (double)(atom->getFormalCharge());
            nConj = 1;
            old_nConj = 0;
            conjNBitVect.reset();
            conjNBitVect[idx] = 1;
            while (nConj > old_nConj) {
              old_nConj = nConj;
              for (i = 0; i < mol.getNumAtoms(); ++i) {
                // if this atom is not marked as conj, move on
                if (!conjNBitVect[i]) {
                  continue;
                }
                // loop over neighbors
                boost::tie(nbrIdx, endNbrs) =
                  mol.getAtomNeighbors(mol.getAtomWithIdx(i));
                for (; nbrIdx != endNbrs; ++nbrIdx) {
                  const Atom *nbrAtom = mol[*nbrIdx].get();
                  nbrAtomType = this->getMMFFAtomType(nbrAtom->getIdx());
                  // if atom type is not 80 or 57, move on
                  if ((nbrAtomType != 57) && (nbrAtomType != 80)) {
                    continue;
                  }
                  // loop over neighbors of the neighbor
                  // if they are nitrogens of type 81, 55 or 56 and
                  // they are not not marked as conjugated yet, do it
                  // and increment the nConj counter by 1
                  boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(nbrAtom);
                  for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
                    const Atom *nbr2Atom = mol[*nbr2Idx].get();
                    // if atom type is not 81, 55 or 56, move on
                    nbrAtomType = this->getMMFFAtomType(nbr2Atom->getIdx());
                    if ((nbrAtomType != 55) && (nbrAtomType != 56)
                      && (nbrAtomType != 81)) {
                      continue;
                    }
                    j = nbr2Atom->getIdx();
                    // if this nitrogen is not yet marked as conjugated,
                    // mark it and increment the counter and eventually
                    // adjust the total formal charge of the conjugated system
                    if (!conjNBitVect[j]) {
                      conjNBitVect[j] = 1;
                      fChg += (double)(nbr2Atom->getFormalCharge());
                      ++nConj;
                    }
                  }
                }
              }
            }
            if (nConj) {
              fChg /= (double)nConj;
            }
          break;
          
          case 61:
            // loop over neighbors
            boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
            for (; nbrIdx != endNbrs; ++nbrIdx) {
              const Atom *nbrAtom = mol[*nbrIdx].get();
              // if it is diazonium, set a +1 formal charge on
              // the secondary nitrogen
              if (this->getMMFFAtomType(nbrAtom->getIdx()) == 42) {
                fChg = 1.0;
              }
            }
          break;

          // non-complicated +1 atom types
          case 34:
            // NR+
            // Quaternary nitrogen
          case 49:
            // O+
            // Oxonium oxygen
          case 51:
            // O=+
            // Oxenium oxygen
          case 54:
            // N+=C
            // Iminium nitrogen
            // N+=N
            // Positively charged nitrogen doubly bonded to N
          case 58:
            // NPD+
            // Aromatic nitrogen in pyridinium
          case 92:
            // LI+
            // Lithium cation
          case 93:
            // NA+
            // Sodium cation
          case 94:
            // K+
            // Potassium cation
          case 97:
            // CU+1
            // Monopositive copper cation
            fChg = 1.0;
          break;
            
          // non-complicated +2 atom types
          case 87:
            // FE+2
            // Dipositive iron cation
          case 95:
            // ZN+2
            // Dipositive zinc cation
          case 96:
            // CA+2
            // Dipositive calcium cation
          case 98:
            // CU+2
            // Dipositive copper cation
          case 99:
            // MG+2
            // Dipositive magnesium cation
            fChg = 2.0;
          break;
          
          // non-complicated +3 atom types
          case 88:
            // FE+3
            // Tripositive iron cation
            fChg = 3.0;
          break;
            
          // non-complicated -1 atom types
          case 35:
            // OM
            // Oxide oxygen on sp3 carbon
            // OM2
            // Oxide oxygen on sp2 carbon
            // OM
            // Oxide oxygen on sp3 nitrogen (not in original MMFF.I Table III)
            // OM2
            // Oxide oxygen on sp2 nitrogen (not in original MMFF.I Table III)
          case 62:
            // NM
            // Anionic divalent nitrogen
          case 89:
            // F-
            // Fluoride anion
          case 90:
            // Cl-
            // Chloride anion
          case 91:
            // BR-
            // Bromide anion
            fChg = -1.0;
          break;
        }
        this->setMMFFFormalCharge(idx, fChg);
      }
      // now we compute partial charges
      // See Halgren, T. MMFF.V, J. Comput. Chem. 1996, 17, 616-641
      // http://dx.doi.org/10.1002/(SICI)1096-987X(199604)17:5/6<616::AID-JCC5>3.0.CO;2-X
      for (idx = 0; idx < mol.getNumAtoms(); ++idx) {
        const Atom *atom = mol.getAtomWithIdx(idx);
        atomType = this->getMMFFAtomType(idx);
        double q0 = this->getMMFFFormalCharge(idx);
        double M = (double)((*mmffProp)(atomType)->crd);
        double v = (*mmffPBCI)(atomType)->fcadj;
        double sumFormalCharge = 0.0;
        double sumPartialCharge = 0.0;
        double nbrFormalCharge;
        std::pair<int, const MMFFChg *> mmffChgParams;

        if (isDoubleZero(v)) {
          // loop over neighbors
          boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
          for (; nbrIdx != endNbrs; ++nbrIdx) {
            const Atom *nbrAtom = mol[*nbrIdx].get();
            nbrFormalCharge = this->getMMFFFormalCharge(nbrAtom->getIdx());
            // if neighbors have a negative formal charge, the latter
            // influences the charge on ipso
            if (nbrFormalCharge < 0.0) {
              q0 += (nbrFormalCharge / (2.0 * (double)(nbrAtom->getDegree())));
            }
          }
        }
        // there is a special case for anionic divalent nitrogen
        // with positively charged neighbor
        if (atomType == 62) {
          boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
          for (; nbrIdx != endNbrs; ++nbrIdx) {
            const Atom *nbrAtom = mol[*nbrIdx].get();
            nbrFormalCharge = this->getMMFFFormalCharge(nbrAtom->getIdx());
            if (nbrFormalCharge > 0.0){
              q0 -= (nbrFormalCharge / 2.0);
            }
          }
        }
        // loop over neighbors
        boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
        for (; nbrIdx != endNbrs; ++nbrIdx) {
          const Atom *nbrAtom = mol[*nbrIdx].get();
          const Bond *bond = mol.getBondBetweenAtoms
            (atom->getIdx(), nbrAtom->getIdx());
          // we need to determine the sign of bond charge
          // increments depending on the bonding relationship
          // i.e. we have parameters for [a,b] bonds
          // but it depends whether ipso is a or b
          unsigned int nbrAtomType = this->getMMFFAtomType(nbrAtom->getIdx());
          unsigned int bondType = this->getMMFFBondType(bond);
          mmffChgParams = mmffChg->getMMFFChgParams(bondType, atomType, nbrAtomType);
          sumPartialCharge += (mmffChgParams.second
            ? (double)(mmffChgParams.first) * ((mmffChgParams.second)->bci)
            : ((*mmffPBCI)(atomType)->pbci - (*mmffPBCI)(nbrAtomType)->pbci));
          nbrFormalCharge = this->getMMFFFormalCharge(nbrAtom->getIdx());
          sumFormalCharge += nbrFormalCharge;
        }
        // we compute ipso partial charge according to
        // equation 15, page 622 MMFF.V paper
        pChg = (1.0 - M * v) * q0 + v * sumFormalCharge + sumPartialCharge;
        this->setMMFFPartialCharge(atom->getIdx(), pChg);
      }
    }
  }
}
