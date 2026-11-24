//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_BRIDGEHEADS_H
#define RD_BRIDGEHEADS_H

#include <vector>

namespace RDKit {
class ROMol;

namespace MolOps {

// a bridgehead is defined by two atom rings that share at least two bonds and
// two bridgehead atoms that are members of those rings and which have at least
// three ring bonds each
struct RDKIT_GRAPHMOL_EXPORT Bridgehead {
  unsigned int atom1 = std::numeric_limits<unsigned int>::max();
  unsigned int atom2 = std::numeric_limits<unsigned int>::max();
  std::vector<unsigned int> ring1;
  std::vector<unsigned int> ring2;
  Bridgehead(unsigned int a1, unsigned int a2,
             const std::vector<unsigned int> &r1,
             const std::vector<unsigned int> &r2)
      : atom1(a1), atom2(a2), ring1(r1), ring2(r2) {}
  bool operator==(const Bridgehead &other) const {
    return atom1 == other.atom1 && atom2 == other.atom2 &&
           ring1 == other.ring1 && ring2 == other.ring2;
  }
};

//! finds bridgehead atoms in a molecule
/*!
 *
 * @param mol the molecule to inspect
 * @return a vector of atom indices corresponding to bridgehead atoms
 */
RDKIT_GRAPHMOL_EXPORT std::vector<unsigned int> findBridgeheadAtoms(
    const ROMol &mol);
RDKIT_GRAPHMOL_EXPORT std::vector<Bridgehead> findBridgeheads(const ROMol &mol);

}  // namespace MolOps
}  // namespace RDKit

#endif