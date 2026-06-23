//
//  Copyright (C) 2026 Tad Hurst, Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RD_MACROMOL_H
#define RD_MACROMOL_H

#include <optional>
#include <string>

#include "RWMol.h"

namespace RDKit {

enum class MonomerClass {
  AA,
  NA,
  CHEM,
  OTHER
};

RDKIT_GRAPHMOL_EXPORT std::string monomerClassToString(
    MonomerClass monomerClass);

RDKIT_GRAPHMOL_EXPORT MonomerClass
stringToMonomerClass(const std::string &className);

class RDKIT_GRAPHMOL_EXPORT MacroMol : public RWMol {
 public:
  unsigned int addMacroAtom(
      MonomerClass monomerClass, std::string templateName,
      std::optional<unsigned int> residueNumber = std::nullopt,
      std::optional<std::string> chainId = std::nullopt,
      std::optional<std::string> insertionCode = std::nullopt);

  void addMacroBond(unsigned int fromAtomIdx, unsigned int toAtomIdx,
                    int fromConnectionPoint, int toConnectionPoint,
                    std::optional<Bond::BondType> bondType = std::nullopt);
};
}  // namespace RDKit

#endif
