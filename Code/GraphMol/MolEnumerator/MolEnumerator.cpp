//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "MolEnumerator.h"
#include <RDGeneral/Exceptions.h>
#include <GraphMol/FileParsers/MolSGroupParsing.h>

namespace RDKit {
namespace MolEnumerator {

void PositionVariationOp::initFromMol() {
  d_variationPoints.clear();
  if (!dp_mol) {
    return;
  }
  for (const auto bond : dp_mol->bonds()) {
    std::string endpts;
    std::string attach;
    if (bond->getPropIfPresent("ENDPTS", endpts) &&
        bond->getPropIfPresent("ATTACH", attach) && attach == "ANY") {
      const auto atom = bond->getBeginAtom();
      if (atom->getAtomicNum() == 0) {
        std::stringstream iss(endpts);
        std::vector<unsigned int> oats =
            RDKit::SGroupParsing::ParseV3000Array<unsigned int>(iss);
        // decrement the indices and do error checking:
        for (auto &oat : oats) {
          --oat;
          if (oat < 0 || oat >= dp_mol->getNumAtoms()) {
            throw ValueErrorException("Bad variation point index");
          }
        }
        d_variationPoints.push_back(std::make_pair(atom->getIdx(), oats));
      }
    }
  }
}

}  // namespace MolEnumerator
}  // namespace RDKit
