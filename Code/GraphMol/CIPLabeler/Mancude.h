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

#include <RDGeneral/BoostStartInclude.h>
#include <boost/rational.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {

namespace CIPLabeler {

class CIPMol;

// boost::rational reduces its numerator and denominator on construction. The
// unreduced denominator is also significant here: it records that an atomic
// number was averaged over a Mancude system, even when the average happens to
// be an integer (for example, 12/3). Keep both the value used for comparison
// and that metadata.
class FractionalAtomicNum {
 public:
  FractionalAtomicNum(int numerator, int denominator)
      : d_numerator{numerator}, d_denominator{denominator} {}

  boost::rational<int> value() const {
    return {d_numerator, d_denominator};
  }

  int numerator() const { return d_numerator; }
  int denominator() const { return d_denominator; }
  bool isAveraged() const { return d_denominator > 1; }

 private:
  int d_numerator;
  int d_denominator;
};

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
std::vector<FractionalAtomicNum> calcFracAtomNums(const CIPMol &mol);

}  // namespace CIPLabeler
}  // namespace RDKit
