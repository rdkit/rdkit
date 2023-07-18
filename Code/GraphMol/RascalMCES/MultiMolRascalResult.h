//
// Copyright (C) David Cosgrove 2023
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// A class to hold the results from an MCES determination for multiple
// molecules.  Holds the SMARTS string of the MCES plus shared pointers
// to the molecules that contain it.

#ifndef MULTIMOLRASCALRESULT_H
#define MULTIMOLRASCALRESULT_H

#include <string>
#include <vector>

#include <GraphMol/ROMol.h>

namespace RDKit {

namespace RascalMCES {
class MultiMolRascalResult {
 public:
  MultiMolRascalResult(std::string smarts, size_t num_atoms, size_t num_bonds,
                       std::vector<std::shared_ptr<RDKit::ROMol>> &match_mols)
      : d_smarts(smarts),
        d_numAtoms(num_atoms),
        d_numBonds(num_bonds),
        d_mols(match_mols) {}

  bool operator<(const MultiMolRascalResult &rhs) {
    if (d_smarts == rhs.d_smarts) {
      if (d_numBonds == rhs.d_numBonds) {
        return d_numAtoms < rhs.d_numAtoms;
      } else {
        return d_numBonds < rhs.d_numBonds;
      }
    } else {
      return d_smarts < rhs.d_smarts;
    }
  }

  std::string d_smarts;
  size_t d_numAtoms, d_numBonds;
  std::vector<std::shared_ptr<RDKit::ROMol>> d_mols;
};
}  // namespace RascalMCES
}  // namespace RDKit

#endif  // MULTIMOLRASCALRESULT_H
