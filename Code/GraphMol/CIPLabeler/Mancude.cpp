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

#include <list>
#include <vector>

#include "Mancude.h"
#include "CIPMol.h"

namespace RDKit {

namespace CIPLabeler {

namespace {

enum class Type : uint8_t {
  Cv4D3,       // =C(X)-
  Nv3D2,       // =N-
  Nv4D3Plus,   // =[N+]<
  Nv2D2Minus,  // -[N-]-
  Cv3D3Minus,  // -[C(X)-]-
  Ov3D2Plus,   // -[O+]=
  Other
};

// Initialize atom types for ring atoms, looking at connectivity,
// bond orders, atomic number and charge.
bool SeedTypes(std::vector<Type> &types, const CIPMol &mol) {
  bool result = false;
  for (const auto &atom : mol.atoms()) {
    const int aidx = atom->getIdx();

    // check ring
    int btypes = atom->getTotalNumHs();
    bool ring = false;
    for (const auto &bond : mol.getBonds(atom)) {
      // Given the possible types we have, we only care
      // for single and double bonds which are in rings.
      switch (mol.getBondOrder(bond)) {
        case 1:
          btypes += 0x00000001;
          break;
        case 2:
          btypes += 0x00000100;
          break;
        default:
          btypes += 0x01000000;
          break;
      }
      if (mol.isInRing(bond)) {
        ring = true;
      }
    }
    if (ring) {
      int q = atom->getFormalCharge();
      switch (atom->getAtomicNum()) {
        case 6:   // C
        case 14:  // Si
        case 32:  // Ge
          if (q == 0 && btypes == 0x0102) {
            types[aidx] = Type::Cv4D3;
          } else if (q == -1 && btypes == 0x0003) {
            types[aidx] = Type::Cv3D3Minus;
            result = true;
          }
          break;
        case 7:   // N
        case 15:  // P
        case 33:  // As
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
        case 8:  // O
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

// Reset (to Type::Other) atom types for atoms that have been
//  given a type but cannot be part of a mancude system (more
//  than one typed neighbor is required for resonance to be possible)
void RelaxTypes(std::vector<Type> &types, const CIPMol &mol) {
  std::list<Atom *> queue;
  std::vector<int> counts(mol.getNumAtoms(), 0);
  for (auto atom : mol.atoms()) {
    const auto aidx = atom->getIdx();
    if (types[aidx] == Type::Other) {
      // This is already Type::Other, no need to reset it!
      continue;
    }
    for (const auto &bond : mol.getBonds(atom)) {
      if (!mol.isInRing(bond)) {
        // No non-ring bonds are Type::Other already
        continue;
      }
      const auto nbr = bond->getOtherAtom(atom);
      if (types[nbr->getIdx()] != Type::Other) {
        ++counts[aidx];
      }
    }

    // Atoms that do not have at least two typed neighbors
    // go to the queue to be reset to Type::Other.
    if (counts[aidx] <= 1) {
      queue.push_back(atom);
    }
  }

  while (!queue.empty()) {
    const auto atom = queue.front();
    queue.pop_front();

    const auto aidx = atom->getIdx();
    if (types[aidx] == Type::Other) {
      // shouldn't happen, we only enqueue
      // typed atoms
      continue;
    }

    // the actual reset
    types[aidx] = Type::Other;

    // Update the counts for the neighbors, and check
    // whether they need to be reset too.
    for (const auto &bond : mol.getBonds(atom)) {
      if (!mol.isInRing(bond)) {
        continue;
      }
      const auto nbr = bond->getOtherAtom(atom);
      const auto nbridx = nbr->getIdx();
      if (types[nbridx] == Type::Other) {
        continue;
      }

      --counts[nbridx];
      if (counts[nbridx] == 1) {
        queue.push_back(nbr);
      }
    }
  }
}

// Mark mol atoms in the same resonant part of the mol.
void VisitPart(std::vector<int> &parts, const std::vector<Type> &types,
               int part, Atom *atom, const CIPMol &mol) {
  Atom *next;
  do {
    next = nullptr;
    for (auto &bond : mol.getBonds(atom)) {
      if (!mol.isInRing(bond)) {
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

// Classify mol atoms into different resonant parts of the mol.
int VisitParts(std::vector<int> &parts, const std::vector<Type> &types,
               const CIPMol &mol) {
  int numparts = 0;
  for (auto &atom : mol.atoms()) {
    int aidx = atom->getIdx();
    if (parts[aidx] == 0 && types[aidx] != Type::Other) {
      parts[aidx] = ++numparts;
      VisitPart(parts, types, parts[aidx], atom, mol);
    }
  }
  return numparts;
}

}  // namespace

std::vector<FractionalAtomicNum> calcFracAtomNums(const CIPMol &mol) {
  const auto num_atoms = mol.getNumAtoms();
  std::vector<FractionalAtomicNum> fractions;
  fractions.reserve(num_atoms);
  for (const auto &atom : mol.atoms()) {
    fractions.emplace_back(atom->getAtomicNum(), 1);
  }

  // Mark all atoms which are potentially part of a resonance system.
  auto types = std::vector<Type>(num_atoms, Type::Other);
  if (!SeedTypes(types, mol)) {
    // No relevant rings.
    return fractions;
  }

  // Filter out atoms which cannot be resonant because
  // of not having the proper environment.
  RelaxTypes(types, mol);

  // Find resonant systems: parts stores the ids of the
  // systems each atom is involved in.
  std::vector<int> parts(num_atoms);
  int numparts = VisitParts(parts, types, mol);

  if (numparts <= 0) {
    // No rings (shouldn't happen, we'd seen it before)
    return fractions;
  }

  // parts are 1-based, so we need to allocate one extra element.
  // This means that resonant_parts[0] will never be used, and
  // we can use it as a shortcut to check whether any resonant
  // systems were found.
  std::vector<bool> resonant_parts(numparts + 1, false);

  for (auto i = 0u; i < num_atoms; ++i) {
    if (parts[i] == 0) {
      continue;
    }
    auto atom = mol.getAtom(i);

    // Find resonant structures caused by relocation of a negative charge.
    if (types[i] == Type::Cv3D3Minus || types[i] == Type::Nv2D2Minus) {
      resonant_parts[parts[i]] = true;
      resonant_parts[0] = true;
    }

    int numerator = 0;
    int denominator = 0;
    for (const auto &nbr : mol.getNeighbors(atom)) {
      if (parts[nbr->getIdx()] == parts[i]) {
        numerator += nbr->getAtomicNum();
        ++denominator;
      }
    }
    fractions[i] = FractionalAtomicNum(numerator, denominator);
  }

  if (!resonant_parts[0]) {
    // No resonant systems.
    return fractions;
  }

  // If there are any resonant parts due to negative charges,
  // update the average atomic number considering relocation
  // of the charge through higher order bonds.
  std::vector<int> numerators(numparts + 1);
  std::vector<int> denominators(numparts + 1);
  for (auto i = 0u; i < num_atoms; ++i) {
    const auto part = parts[i];
    if (part == 0 || !resonant_parts[part]) {
      continue;
    }

    ++denominators[part];
    const auto atom = mol.getAtom(i);
    for (const auto &bond : mol.getBonds(atom)) {
      const auto nbr = bond->getOtherAtom(atom);
      const auto bord = mol.getBondOrder(bond);
      if (bord > 1 && parts[nbr->getIdx()] == part) {
        numerators[part] += (bord - 1) * nbr->getAtomicNum();
      }
    }
  }

  // Once all the numerators and denominators for the different
  // resonant parts of the mol have been calculated, distribute
  // the values to the atoms in the corresponding parts.
  for (auto i = 0u; i < num_atoms; ++i) {
    const auto part = parts[i];
    if (part != 0 && resonant_parts[part] && denominators[part] != 0) {
      fractions[i] = FractionalAtomicNum(numerators[part], denominators[part]);
    }
  }

  return fractions;
}

}  // namespace CIPLabeler
}  // namespace RDKit
