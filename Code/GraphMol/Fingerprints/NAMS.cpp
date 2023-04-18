//
//  Copyright (C) 2023 Rocco Moretti and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// This file implements the non-contiguous atom molecular structural similarity (NAMS)
// method of Teixeira & Falcao (https://pubs.acs.org/doi/abs/10.1021/ci400324u)
// for an RDKit context

#include <GraphMol/Fingerprints/NAMS.h>
#include <GraphMol/ROMol.h>

namespace RDKit {
namespace NAMS {

NAMSMolInfo * getNAMSMolInfo(const ROMol &mol) {
  return new NAMSMolInfo; 
}

double getNAMSSimilarity(const NAMSMolInfo & molinfo1, const NAMSMolInfo & molinfo2) {
  return 0;
}

} // namespace NAMS
} // namespace RDKit
