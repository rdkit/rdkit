//
//  Copyright (C) 2026 Glysade LLC and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*! \file SpacialScore.h

  \brief C++ port of rdkit.Chem.SpacialScore

  Krzyzanowski, A.; Pahl, A.; Grigalunas, M.; Waldmann, H. Spacial
  Score - A Comprehensive Topological Indicator for Small-Molecule
  Complexity. J. Med. Chem. 2023.
  https://doi.org/10.1021/acs.jmedchem.3c00689

*/
#include <RDGeneral/export.h>
#ifndef RD_SPACIALSCORE_H
#define RD_SPACIALSCORE_H

#include <string>

namespace RDKit {
class ROMol;
namespace Descriptors {

//! calculates the spacial score (SPS) of a molecule
/*!
  \param mol        the molecule of interest
  \param normalize  if set (the default) the score is divided by the number
                    of heavy atoms (nSPS). Molecules without heavy atoms
                    return 0 in that case.
*/
RDKIT_DESCRIPTORS_EXPORT double calcSPS(const ROMol &mol,
                                        bool normalize = true);
const std::string SPSVersion = "1.0.0";

}  // namespace Descriptors
}  // namespace RDKit
#endif
