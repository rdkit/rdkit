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

#include <string>
#include <vector>

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

struct RDKIT_GRAPHMOL_EXPORT MacroBondProps {
  int beginAttachPt;
  int endAttachPt;
  Bond::BondType bondType;
  bool isDirectional;
};

class RDKIT_GRAPHMOL_EXPORT MacroMol : public RWMol {
 public:
  using RWMol::addBond;

  unsigned int addMacroAtom(MonomerClass monomerClass, std::string symbol);

  unsigned int addMacroBond(unsigned int beginAtomIdx, unsigned int endAtomIdx,
                            int beginAttachPt, int endAttachPt,
                            Bond::BondType bondType = Bond::BondType::SINGLE);

  unsigned int addAtomToMacroAtomBond(
      unsigned int beginAtomIdx, unsigned int endMacroAtomIdx, int endAttachPt,
      Bond::BondType bondType = Bond::BondType::SINGLE);

  unsigned int addMacroAtomToAtomBond(
      unsigned int beginMacroAtomIdx, unsigned int endAtomIdx,
      int beginAttachPt, Bond::BondType bondType = Bond::BondType::SINGLE);

  unsigned int addBond(unsigned int beginAtomIdx, unsigned int endAtomIdx,
                       Bond::BondType bondType = Bond::BondType::SINGLE);
};
}  // namespace RDKit

#endif
