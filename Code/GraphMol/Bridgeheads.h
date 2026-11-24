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

#include <iosfwd>
#include <vector>
#include <limits>

namespace RDKit {
class ROMol;

namespace MolOps {

// a bridgehead is defined by two atom rings that share at least two bonds and
// two bridgehead atoms that are members of those rings and which have at least
// three ring bonds from the shared rings
struct RDKIT_GRAPHMOL_EXPORT Bridgehead {
  unsigned int atom1 =
      std::numeric_limits<unsigned int>::max();  // index of the first
                                                 // bridgehead atom
  unsigned int atom2 =
      std::numeric_limits<unsigned int>::max();  // index of the second
                                                 // bridgehead atom
  std::vector<unsigned int> ring1;               // one of the atom rings
  std::vector<unsigned int> ring2;               // the other atom ring
  Bridgehead() = default;
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
 *  Bridgehead atoms are part of two rings that share at least two bonds and
 *  have at least three ring bonds from these shared rings.
 *
 * @param mol the molecule to inspect
 * @return a vector of atom indices corresponding to bridgehead atoms
 */
RDKIT_GRAPHMOL_EXPORT std::vector<unsigned int> findBridgeheadAtoms(
    const ROMol &mol);
//! finds bridgeheads in a molecule
/*!
 *
 * A bridgehead is defined by two atom rings that share at least two bonds and
 * two bridgehead atoms that are members of those rings and which have at least
 * three ring bonds from the shared rings.
 *
 * @param mol the molecule to inspect
 * @return a vector of Bridgehead objects
 */
RDKIT_GRAPHMOL_EXPORT std::vector<Bridgehead> findBridgeheads(const ROMol &mol);

namespace details {
// returns indices for all pairs of rings that share at least two bonds, along
// with the indices of the bonds shared between those rings
RDKIT_GRAPHMOL_EXPORT
std::vector<
    std::pair<std::pair<unsigned int, unsigned int>, std::vector<unsigned int>>>
findRingPairs(const ROMol &mol, unsigned int minCommonBonds = 2);

}  // namespace details
}  // namespace MolOps
}  // namespace RDKit

//! allows Bridgehead objects to be dumped to streams
RDKIT_GRAPHMOL_EXPORT std::ostream &operator<<(
    std::ostream &target, const RDKit::MolOps::Bridgehead &);

#endif