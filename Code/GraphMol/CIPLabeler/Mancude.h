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

namespace RDKit {

namespace CIPLabeler {

class CIPMol;

class Fraction {
private:
  friend std::vector<Fraction> calcFracAtomNums(const CIPMol *mol);
  int d_numerator;
  int d_denominator = 1;

public:
  static int compare(int anum, int aden, int bnum, int bden);

  Fraction() = delete;
  explicit Fraction(int num);
  Fraction(int num, int den);

  int numerator() const;

  int denominator() const;

  int compareTo(const Fraction &o) const;
};

enum class Type {
  Cv4D4,      // =CH-
  Nv3D2,      // =N-
  Nv4D3Plus,  // =[N+]<
  Nv2D2Minus, // -[N-]-
  Cv3D3Minus, // -[CH-]-
  Ov3D2Plus,  // -[O+]=
  Other
};

std::vector<Fraction> calcFracAtomNums(const CIPMol *mol);

} // namespace CIPLabeler
} // namespace RDKit
