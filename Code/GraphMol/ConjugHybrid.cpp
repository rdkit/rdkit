//
//  Copyright (C) 2001-2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "ROMol.h"
#include "RWMol.h"
#include "Atom.h"
#include "Bond.h"
#include "MolOps.h"
#include "PeriodicTable.h"
#include "AtomIterators.h"
#include "BondIterators.h"

namespace RDKit {
// local utility namespace:
namespace {
bool isAtomConjugCand(const RDMol &mol, std::uint32_t atomIdx,
                      const AtomData &at) {
  // return false for neutral atoms where the current valence exceeds the
  // minimal valence for the atom. logic: if we're hypervalent we aren't
  // conjugated
  const auto &vals =
      PeriodicTable::getTable()->getValenceList(at.getAtomicNum());
  if (!at.getFormalCharge() && vals.front() >= 0 &&
      at.getTotalValence() > static_cast<unsigned int>(vals.front())) {
    return false;
  }
  // the second check here is for Issue211, where the c-P bonds in
  // Pc1ccccc1 were being marked as conjugated.  This caused the P atom
  // itself to be SP2 hybridized.  This is wrong.  For now we'll do a quick
  // hack and forbid this check from adding conjugation to anything out of
  // the first row of the periodic table.  (Conjugation in aromatic rings
  // has already been attended to, so this is safe.)
  int nouter = PeriodicTable::getTable()->getNouterElecs(at.getAtomicNum());
  auto res = ((at.getAtomicNum() <= 10) || (nouter != 5 && nouter != 6) ||
              (nouter == 6 && mol.getAtomTotalDegree(atomIdx) < 2u)) &&
             MolOps::countAtomElec(mol, atomIdx) > 0;

  return res;
}

void markConjAtomBonds(RDMol &mol, std::uint32_t atomIdx) {
  const AtomData &at = mol.getAtom(atomIdx);
  if (!isAtomConjugCand(mol, atomIdx, at)) {
    return;
  }

  // make sure that have either 2 or 3 substitutions on this atom
  int sbo = mol.getAtomTotalDegree(atomIdx);
  if ((sbo < 2) || (sbo > 3)) {
    return;
  }
  auto [beg1, end1] = mol.getAtomBonds(atomIdx);
  for (; beg1 != end1; ++beg1) {
    uint32_t bondIdx = *beg1;
    BondData &bnd1 = mol.getBond(bondIdx);
    if (bnd1.getTwiceValenceContrib(atomIdx) < 3 ||
        !isAtomConjugCand(mol, bnd1.getOtherAtomIdx(atomIdx),
                          mol.getAtom(bnd1.getOtherAtomIdx(atomIdx)))) {
      continue;
    }
    auto [beg2, end2] = mol.getAtomBonds(atomIdx);
    for (; beg2 != end2; ++beg2) {
      if (beg1 == beg2) {
        continue;
      }
      BondData &bnd2 = mol.getBond(*beg2);
      auto atomIdx2 = bnd2.getOtherAtomIdx(atomIdx);
      auto at2 = mol.getAtom(atomIdx2);
      sbo = mol.getAtomTotalDegree(atomIdx2);
      if (sbo > 3) {
        continue;
      }
      if (isAtomConjugCand(mol, atomIdx2, at2)) {
        bnd1.setIsConjugated(true);
        bnd2.setIsConjugated(true);
      }
    }
  }
}

int numBondsPlusLonePairs(const RDMol &mol, std::uint32_t atomIdx,
                          const AtomData &atom) {
  int deg = mol.getAtomTotalDegree(atomIdx);

  auto [beg, end] = mol.getAtomBonds(atomIdx);
  for (; beg != end; ++beg) {
    uint32_t bondIdx = *beg;
    const BondData &bond = mol.getBond(bondIdx);
    if (bond.getBondType() == BondEnums::BondType::ZERO ||
        (isDative(bond.getBondType()) && atomIdx != bond.getEndAtomIdx())) {
      --deg;
    }
  }

  if (atom.getAtomicNum() <= 1) {
    return deg;
  }
  int nouter = PeriodicTable::getTable()->getNouterElecs(atom.getAtomicNum());
  int totalValence = atom.getTotalValence();
  int chg = atom.getFormalCharge();

  int numFreeElectrons = nouter - (totalValence + chg);
  if (totalValence + nouter - chg < 8) {
    // we're below an octet, so we need to think
    // about radicals:
    int numRadicals = atom.getNumRadicalElectrons();
    int numLonePairs = (numFreeElectrons - numRadicals) / 2;
    return deg + numLonePairs + numRadicals;
  } else {
    int numLonePairs = numFreeElectrons / 2;
    return deg + numLonePairs;
  }
}
}  // namespace

namespace MolOps {
bool atomHasConjugatedBond(const Atom *at) {
  PRECONDITION(at, "bad atom");

  auto &mol = at->getOwningMol();
  for (const auto bnd : mol.atomBonds(at)) {
    if (bnd->getIsConjugated()) {
      return true;
    }
  }
  return false;
}
bool atomHasConjugatedBond(const RDMol &mol, std::uint32_t atomIdx) {
  auto [beg, end] = mol.getAtomBonds(atomIdx);
  for (; beg != end; ++beg) {
    uint32_t bondIdx = *beg;
    const BondData &bnd = mol.getBond(bondIdx);
    if (bnd.getIsConjugated()) {
      return true;
    }
  }
  return false;
}

void setConjugation(ROMol &mol) {
  // start with all bonds being marked unconjugated
  // except for aromatic bonds
  auto &rdmol = mol.asRDMol();

  for (uint32_t bondIdx = 0, numBonds = rdmol.getNumBonds(); bondIdx < numBonds;
       ++bondIdx) {
    BondData &bond = rdmol.getBond(bondIdx);
    bond.setIsConjugated(bond.getIsAromatic());
  }

  // loop over each atom and check if the bonds connecting to it can
  // be conjugated

  for (uint32_t atomIdx = 0, numAtoms = rdmol.getNumAtoms(); atomIdx < numAtoms;
       ++atomIdx) {
    markConjAtomBonds(rdmol, atomIdx);
  }
}

void setHybridization(ROMol &mol) {
  auto &rdmol = mol.asRDMol();

  for (uint32_t atomIdx = 0, numAtoms = rdmol.getNumAtoms(); atomIdx < numAtoms;
       ++atomIdx) {
    auto &atom = rdmol.getAtom(atomIdx);
    if (atom.getAtomicNum() == 0) {
      atom.setHybridization(AtomEnums::HybridizationType::UNSPECIFIED);
    } else {
      // if the stereo spec matches the coordination number, this is easy
      const auto totalDegree = rdmol.getAtomTotalDegree(atomIdx);
      switch (atom.getChiralTag()) {
        case Atom::ChiralType::CHI_TETRAHEDRAL:
        case Atom::ChiralType::CHI_TETRAHEDRAL_CW:
        case Atom::ChiralType::CHI_TETRAHEDRAL_CCW:
          if (totalDegree == 4) {
            atom.setHybridization(AtomEnums::HybridizationType::SP3);
            continue;
          }
          break;
        case Atom::ChiralType::CHI_SQUAREPLANAR:
          if (totalDegree <= 4 && totalDegree >= 2) {
            atom.setHybridization(AtomEnums::HybridizationType::SP2D);
            continue;
          }
          break;
        case Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL:
          if (totalDegree <= 5 && totalDegree >= 2) {
            atom.setHybridization(AtomEnums::HybridizationType::SP3D);
            continue;
          }
          break;
        case Atom::ChiralType::CHI_OCTAHEDRAL:
          if (totalDegree <= 6 && totalDegree >= 2) {
            atom.setHybridization(AtomEnums::HybridizationType::SP3D2);
            continue;
          }
          break;
        default:
          break;
      }
      // otherwise we have to do some work
      int norbs = 0;
      // try to be smart for early elements, but for later
      // ones just use the degree
      // FIX: we should probably also be using the degree for metals
      if (atom.getAtomicNum() < 89) {
        norbs = numBondsPlusLonePairs(rdmol, atomIdx, atom);
      } else {
        norbs = totalDegree;
      }
      switch (norbs) {
        case 0:
          // This occurs for things like Na+
          atom.setHybridization(AtomEnums::HybridizationType::S);
          break;
        case 1:
          atom.setHybridization(AtomEnums::HybridizationType::S);
          break;
        case 2:
          atom.setHybridization(AtomEnums::HybridizationType::SP);
          break;
        case 3:
          atom.setHybridization(AtomEnums::HybridizationType::SP2);
          break;
        case 4:
          // potentially SP3, but we'll set it down to SP2
          // if we have a conjugated bond (like the second O
          // in O=CO)
          // we'll also avoid setting the hybridization down to
          // SP2 in the case of an atom with degree higher than 3
          // (e.g. things like CP1(C)=CC=CN=C1C, where the P
          //   has norbs = 4, and a conjugated bond, but clearly should
          //   not be SP2)
          // This is Issue276
          if (totalDegree > 3 ||
              !MolOps::atomHasConjugatedBond(rdmol, atomIdx)) {
            atom.setHybridization(AtomEnums::HybridizationType::SP3);
          } else {
            atom.setHybridization(AtomEnums::HybridizationType::SP2);
          }
          break;
        case 5:
          atom.setHybridization(AtomEnums::HybridizationType::SP3D);
          break;
        case 6:
          atom.setHybridization(AtomEnums::HybridizationType::SP3D2);
          break;
        default:
          atom.setHybridization(AtomEnums::HybridizationType::UNSPECIFIED);
      }
    }
  }
}
}  // end of namespace MolOps
}  // end of namespace RDKit
