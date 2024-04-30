//
//  Copyright (C) 2001-2021 Greg Landrum and Rational Discovery LLC
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
bool isAtomConjugCand(const Atom *at) {
  PRECONDITION(at, "bad atom");
  // the second check here is for Issue211, where the c-P bonds in
  // Pc1ccccc1 were being marked as conjugated.  This caused the P atom
  // itself to be SP2 hybridized.  This is wrong.  For now we'll do a quick
  // hack and forbid this check from adding conjugation to anything out of
  // the first row of the periodic table.  (Conjugation in aromatic rings
  // has already been attended to, so this is safe.)
  int nouter = PeriodicTable::getTable()->getNouterElecs(at->getAtomicNum());
  return ((at->getAtomicNum() <= 10) || (nouter != 5 && nouter != 6) ||
          (nouter == 6 && at->getTotalDegree() < 2u)) &&
         (MolOps::countAtomElec(at) > 0);
}

void markConjAtomBonds(Atom *at) {
  PRECONDITION(at, "bad atom");
  if (!isAtomConjugCand(at)) {
    return;
  }
  auto &mol = at->getOwningMol();

  int atx = at->getIdx();
  // make sure that have either 2 or 3 substitutions on this atom
  int sbo = at->getDegree() + at->getTotalNumHs();
  if ((sbo < 2) || (sbo > 3)) {
    return;
  }

  for (const auto &nbri : boost::make_iterator_range(mol.getAtomBonds(at))) {
    auto bnd1 = mol[nbri];
    if (bnd1->getValenceContrib(at) < 1.5) {
      continue;
    }

    for (const auto &nbrj : boost::make_iterator_range(mol.getAtomBonds(at))) {
      auto bnd2 = mol[nbrj];
      if (bnd1 == bnd2) {
        continue;
      }
      auto at2 = mol.getAtomWithIdx(bnd2->getOtherAtomIdx(atx));
      sbo = at2->getDegree() + at2->getTotalNumHs();
      if (sbo > 3) {
        continue;
      }
      if (isAtomConjugCand(at2)) {
        bnd1->setIsConjugated(true);
        bnd2->setIsConjugated(true);
      }
    }
  }
}

int numBondsPlusLonePairs(Atom *at) {
  PRECONDITION(at, "bad atom");
  int deg = at->getTotalDegree();

  auto &mol = at->getOwningMol();
  for (const auto bond : mol.atomBonds(at)) {
    if (bond->getBondType() == Bond::ZERO ||
        (isDative(*bond) && at->getIdx() != bond->getEndAtomIdx())) {
      --deg;
    }
  }

  if (at->getAtomicNum() <= 1) {
    return deg;
  }
  int nouter = PeriodicTable::getTable()->getNouterElecs(at->getAtomicNum());
  int totalValence = at->getExplicitValence() + at->getImplicitValence();
  int chg = at->getFormalCharge();

  int numFreeElectrons = nouter - (totalValence + chg);
  if (totalValence + nouter - chg < 8) {
    // we're below an octet, so we need to think
    // about radicals:
    int numRadicals = at->getNumRadicalElectrons();
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
  for (const auto &nbri : boost::make_iterator_range(mol.getAtomBonds(at))) {
    auto bnd = mol[nbri];
    if (bnd->getIsConjugated()) {
      return true;
    }
  }
  return false;
}

void setConjugation(ROMol &mol) {
  // start with all bonds being marked unconjugated
  // except for aromatic bonds
  for (auto bond : mol.bonds()) {
    bond->setIsConjugated(bond->getIsAromatic());
  }

  // loop over each atom and check if the bonds connecting to it can
  // be conjugated
  for (auto atom : mol.atoms()) {
    markConjAtomBonds(atom);
  }
}

void setHybridization(ROMol &mol) {
  for (auto atom : mol.atoms()) {
    if (atom->getAtomicNum() == 0) {
      atom->setHybridization(Atom::UNSPECIFIED);
    } else {
      // if the stereo spec matches the coordination number, this is easy
      switch (atom->getChiralTag()) {
        case Atom::ChiralType::CHI_TETRAHEDRAL:
        case Atom::ChiralType::CHI_TETRAHEDRAL_CW:
        case Atom::ChiralType::CHI_TETRAHEDRAL_CCW:
          if (atom->getTotalDegree() == 4) {
            atom->setHybridization(Atom::HybridizationType::SP3);
            continue;
          }
          break;
        case Atom::ChiralType::CHI_SQUAREPLANAR:
          if (atom->getTotalDegree() <= 4 && atom->getTotalDegree() >= 2) {
            atom->setHybridization(Atom::HybridizationType::SP2D);
            continue;
          }
          break;
        case Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL:
          if (atom->getTotalDegree() <= 5 && atom->getTotalDegree() >= 2) {
            atom->setHybridization(Atom::HybridizationType::SP3D);
            continue;
          }
          break;
        case Atom::ChiralType::CHI_OCTAHEDRAL:
          if (atom->getTotalDegree() <= 6 && atom->getTotalDegree() >= 2) {
            atom->setHybridization(Atom::HybridizationType::SP3D2);
            continue;
          }
          break;
        default:
          break;
      }
      // otherwise we have to do some work
      int norbs;
      // try to be smart for early elements, but for later
      // ones just use the degree
      // FIX: we should probably also be using the degree for metals
      if (atom->getAtomicNum() < 89) {
        norbs = numBondsPlusLonePairs(atom);
      } else {
        norbs = atom->getTotalDegree();
      }
      switch (norbs) {
        case 0:
          // This occurs for things like Na+
          atom->setHybridization(Atom::S);
          break;
        case 1:
          atom->setHybridization(Atom::S);
          break;
        case 2:
          atom->setHybridization(Atom::SP);
          break;
        case 3:
          atom->setHybridization(Atom::SP2);
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
          if (atom->getTotalDegree() > 3 ||
              !MolOps::atomHasConjugatedBond(atom)) {
            atom->setHybridization(Atom::SP3);
          } else {
            atom->setHybridization(Atom::SP2);
          }
          break;
        case 5:
          atom->setHybridization(Atom::SP3D);
          break;
        case 6:
          atom->setHybridization(Atom::SP3D2);
          break;
        default:
          atom->setHybridization(Atom::UNSPECIFIED);
      }
    }
  }
}
}  // end of namespace MolOps
}  // end of namespace RDKit
