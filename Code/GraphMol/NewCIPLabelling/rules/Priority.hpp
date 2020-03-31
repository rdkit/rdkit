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
namespace NewCIPLabelling {

/**
 * Holds some properties that are determined when sorting/prioritising ligands.
 *
 */
class Priority {
private:
  bool unique;
  bool pseudoAsym;
  int ruleIdx;

public:
  Priority() = delete;

  Priority(bool unique, int ruleIdx, bool pseudoAsym)
      : unique{unique}, pseudoAsym{pseudoAsym}, ruleIdx{ruleIdx} {}

  /**
   * Indicates whether the ligands were unique (i.e. could be ordered)
   *
   * @return whether the ligands were unique
   */
  bool isUnique() const { return unique; }

  int getRuleIdx() const { return ruleIdx; }

  /**
   * Indicates the descriptor type used to. This allows methods that represent
   * pseudo-asymmetric molecules to indicate that the centre is
   * pseudo-asymmetric.
   *
   * @return The type of the descriptor that should be assigned
   */
  bool isPseudoAsymetric() const { return pseudoAsym; }
};

} // namespace NewCIPLabelling
} // namespace RDKit