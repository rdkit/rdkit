//
//
//  Copyright (C) 2020 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

///
/// Utilities to handle atomic numbers in Mancude (maximum number
/// of noncumulative double bonds) rings
///
/// This guarantees that aromatic rings containing heteroatoms
/// are always resolved in the same way
///

#include <vector>
#include "Mancude.h"
#include "CIPMol.h"

namespace RDKit {

namespace CIPLabeler {

int Fraction::compare(int anum, int aden, int bnum, int bden) {
  unsigned long x = anum * bden;
  unsigned long y = bnum * aden;

  return x < y ? -1 : (x == y ? 0 : 1);
}

Fraction::Fraction(int num) : d_numerator{num} {}

Fraction::Fraction(int num, int den) : d_numerator{num}, d_denominator{den} {}

int Fraction::numerator() const { return d_numerator; }

int Fraction::denominator() const { return d_denominator; }

int Fraction::compareTo(const Fraction &o) const {
  return compare(this->d_numerator, this->d_denominator, o.d_numerator,
                 o.d_denominator);
}

namespace {

bool SeedTypes(std::vector<Type> &types, const CIPMol *mol) {
  bool result = false;
  for (const auto &atom : mol->atoms()) {
    const int aidx = atom->getIdx();

    // check ring
    int btypes = atom->getTotalNumHs();
    bool ring = false;
    for (const auto &bond : mol->getBonds(atom)) {
      switch (mol->getBondOrder(bond)) {
      case 1:
        btypes += 0x00000001;
        break;
      case 2:
        btypes += 0x00000100;
        break;
      case 3:
        btypes += 0x00010000;
        break;
      default:
        btypes += 0x01000000;
        break;
      }
      if (mol->isInRing(bond)) {
        ring = true;
      }
    }
    if (ring) {
      int q = atom->getFormalCharge();
      switch (atom->getAtomicNum()) {
      case 6:  // C
      case 14: // Si
      case 32: // Ge
        if (q == 0 && btypes == 0x0102) {
          types[aidx] = Type::Cv4D4;
        } else if (q == -1 && btypes == 0x0003) {
          types[aidx] = Type::Cv3D3Minus;
          result = true;
        }
        break;
      case 7:  // N
      case 15: // P
      case 33: // As
        if (q == 0 && btypes == 0x0101) {
          types[aidx] = Type::Nv3D2;
          result = true;
        } else if (q == -1 && btypes == 0x0002) {
          types[aidx] = Type::Nv2D2Minus;
          result = true;
        } else if (q == +1 && btypes == 0x0102) {
          types[aidx] = Type::Nv4D3Plus;
          result = true;
        }
        break;
      case 8: // O
        if (q == 1 && btypes == 0x0101) {
          types[aidx] = Type::Ov3D2Plus;
          result = true;
        }
        break;
      }
    }
  }
  return result;
}

void RelaxTypes(std::vector<Type> &types, const CIPMol *mol) {
  auto counts = std::vector<int>(mol->getNumAtoms());
  auto queue = std::vector<Atom *>();

  for (auto atom : mol->atoms()) {
    const auto aidx = atom->getIdx();
    for (const auto &nbr : mol->getNeighbors(atom)) {
      if (types[nbr->getIdx()] != Type::Other) {
        ++counts[aidx];
      }
    }

    if (counts[aidx] == 1) {
      queue.push_back(atom);
    }
  }

  for (auto pos = 0u; pos < queue.size(); ++pos) {
    auto atom = queue[0];

    auto aidx = atom->getIdx();
    if (types[aidx] != Type::Other) {
      types[aidx] = Type::Other;

      for (auto &nbr : mol->getNeighbors(atom)) {
        auto nbridx = nbr->getIdx();
        --counts[nbridx];
        if (counts[nbridx] == 1) {
          queue.push_back(nbr);
        }
      }
    }
  }
}

void VisitPart(std::vector<int> &parts, const std::vector<Type> &types,
               int part, Atom *atom, const CIPMol *mol) {
  Atom *next;
  do {
    next = nullptr;
    for (auto &bond : mol->getBonds(atom)) {
      if (!mol->isInRing(bond)) {
        continue;
      }

      auto nbr = bond->getOtherAtom(atom);
      int aidx = nbr->getIdx();

      if (parts[aidx] == 0 && types[aidx] != Type::Other) {
        parts[aidx] = part;
        if (next != nullptr) {
          VisitPart(parts, types, part, nbr, mol);
        } else {
          next = nbr;
        }
      }
    }
    atom = next;
  } while (atom != nullptr);
}

int VisitParts(std::vector<int> &parts, const std::vector<Type> &types,
               const CIPMol *mol) {
  int numparts = 0;
  for (auto &atom : mol->atoms()) {
    int aidx = atom->getIdx();
    if (parts[aidx] == 0 && types[aidx] != Type::Other) {
      parts[aidx] = ++numparts;
      VisitPart(parts, types, parts[aidx], atom, mol);
    }
  }
  return numparts;
}

} // namespace

std::vector<Fraction> calcFracAtomNums(const CIPMol *mol) {
  const auto num_atoms = mol->getNumAtoms();
  auto fractions = std::vector<Fraction>();
  fractions.reserve(num_atoms);

  for (const auto &atom : mol->atoms()) {
    fractions.emplace_back(atom->getAtomicNum(), 1);
  }

  auto types = std::vector<Type>(num_atoms, Type::Other);
  if (SeedTypes(types, mol)) {
    RelaxTypes(types, mol);

    auto parts = std::vector<int>(num_atoms);
    int numparts = VisitParts(parts, types, mol);

    auto resparts = std::vector<int>(numparts);
    int numres = 0;

    if (numparts > 0) {
      for (auto i = 0u; i < num_atoms; ++i) {
        if (parts[i] == 0) {
          continue;
        }
        auto atom = mol->getAtom(i);

        if (types[i] == Type::Cv3D3Minus || types[i] == Type::Nv2D2Minus) {
          int j = 0;
          for (; j < numres; ++j)
            if (resparts[j] == parts[i]) {
              break;
            }
          if (j >= numres) {
            resparts[numres] = parts[i];
            ++numres;
          }
        }

        fractions[i].d_numerator = 0;
        fractions[i].d_denominator = 0;

        for (const auto &nbr : mol->getNeighbors(atom)) {
          if (parts[nbr->getIdx()] == parts[i]) {
            fractions[i].d_numerator += nbr->getAtomicNum();
            ++fractions[i].d_denominator;
          }
        }
      }
    }

    if (numres > 0) {
      for (int j = 0; j < numres; ++j) {
        auto frac = Fraction(0, 0);
        int part = resparts[j];
        for (auto i = 0u; i < num_atoms; ++i) {
          if (parts[i] == part) {
            fractions[i] = frac;
            ++frac.d_denominator;
            auto atom = mol->getAtom(i);

            for (auto &bond : mol->getBonds(atom)) {
              auto nbr = bond->getOtherAtom(atom);
              int bord = mol->getBondOrder(bond);
              if (bord > 1 && parts[nbr->getIdx()] == part) {
                frac.d_numerator += (bord - 1) * nbr->getAtomicNum();
              }
            }
          }
        }
      }
    }
  }
  return fractions;
}

} // namespace CIPLabeler
} // namespace RDKit
