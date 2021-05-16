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
#pragma once

namespace RDKit {
namespace CIPLabeler {

/**
 * Holds some properties that are determined when sorting/prioritising
 * substituents.
 *
 */
class Priority {
 public:
  Priority() = delete;

  Priority(bool unique, bool pseudoAsym)
      : d_unique{unique}, d_pseudoAsym{pseudoAsym} {}

  /**
   * Indicates whether the substituents were unique (i.e. could be ordered)
   *
   * @return whether the substituents were unique
   */
  bool isUnique() const { return d_unique; }

  /**
   * Indicates the descriptor type used to. This allows methods that represent
   * pseudo-asymmetric molecules to indicate that the centre is
   * pseudo-asymmetric.
   *
   * @return The type of the descriptor that should be assigned
   */
  bool isPseudoAsymetric() const { return d_pseudoAsym; }

 private:
  bool d_unique;
  bool d_pseudoAsym;
};

}  // namespace CIPLabeler
}  // namespace RDKit
