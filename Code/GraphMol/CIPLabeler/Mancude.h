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

#pragma once

#include <vector>

#include <boost/rational.hpp>

namespace RDKit {

namespace CIPLabeler {

class CIPMol;

enum class Type {
  Cv4D3,       // =C(X)-
  Nv3D2,       // =N-
  Nv4D3Plus,   // =[N+]<
  Nv2D2Minus,  // -[N-]-
  Cv3D3Minus,  // -[C(X)-]-
  Ov3D2Plus,   // -[O+]=
  Other
};

/**
 * Calculate fractional atomic numbers for all atoms in the mol.
 * Using fractional atomic numbers makes sure that atoms in rings
 * that have resonant structures are always considered with the same
 * priority.
 *
 */
std::vector<boost::rational<int>> calcFracAtomNums(const CIPMol &mol);

}  // namespace CIPLabeler
}  // namespace RDKit
