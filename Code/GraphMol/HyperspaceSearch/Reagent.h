//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// This class holds a Reagent that will be part of a ReagentSet.

#ifndef RDKIT_REAGENT_H
#define RDKIT_REAGENT_H

#include <string>

#include <DataStructs/ExplicitBitVect.h>

namespace RDKit {
class Atom;
class ROMol;

namespace HyperspaceSSSearch {

class Reagent {
 public:
  Reagent(const std::string &smi, const std::string &id)
      : d_smiles(smi), d_id(id) {}

  std::string d_smiles;
  std::string d_id;

  const std::unique_ptr<ROMol> &mol() const;
  const std::unique_ptr<ExplicitBitVect> &pattFP() const;
  const std::vector<std::shared_ptr<ROMol>> &connRegions() const;

 private:
  // All mutable because of lazy evaluation.
  mutable std::unique_ptr<ROMol> d_mol;
  mutable std::unique_ptr<ExplicitBitVect> d_pattFP;
  // SMILES strings of any connector regions.  Normally there will only
  // be 1 or 2.
  mutable std::vector<std::shared_ptr<ROMol>> d_connRegions;
};

// Return a molecule containing the portions of the molecule starting at
// each dummy atom and going out up to 3 bonds.  There may be more than
// 1 fragment if there are dummy atoms more than 3 bonds apart, and there
// may be fragments with more than 1 dummy atom if their fragments fall
// within 3 bonds of each other.  E.g. the molecule [1*]CN(C[2*])Cc1ccccc1
// will give [1*]CN(C)C[1*].  The 2 dummy atoms are 4 bonds apart, but the
// fragments overlap.  All dummy atoms given isotope 1 whatever they had before.
std::unique_ptr<ROMol> getConnRegion(const ROMol &mol);

}  // namespace HyperspaceSSSearch
}  // namespace RDKit

#endif  // RDKIT_REAGENT_H
