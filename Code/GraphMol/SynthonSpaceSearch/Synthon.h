//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RDKIT_REAGENT_H
#define RDKIT_REAGENT_H

#include <string>

#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/ROMol.h>

namespace RDKit {
class Atom;

namespace SynthonSpaceSearch {

// This class holds a Synthon that will be part of a SynthonSet.
class RDKIT_SYNTHONSPACESEARCH_EXPORT Synthon {
 public:
  Synthon() = default;
  Synthon(const std::string &smi, const std::string &id);
  Synthon(const Synthon &other);
  Synthon(Synthon &&other) = default;
  Synthon &operator=(const Synthon &other);
  Synthon &operator=(Synthon &&other) = default;

  [[nodiscard]] const std::string &getSmiles() const { return d_smiles; }
  [[nodiscard]] const std::string &getId() const { return d_id; }
  [[nodiscard]] const std::unique_ptr<ROMol> &getMol() const;
  [[nodiscard]] const std::unique_ptr<ExplicitBitVect> &getPattFP() const;
  [[nodiscard]] const std::vector<std::shared_ptr<ROMol>> &getConnRegions()
      const;

  // Writes to/reads from a binary stream.
  void writeToDBStream(std::ostream &os) const;
  void readFromDBStream(std::istream &is);

 private:
  std::string d_smiles;
  std::string d_id;

  std::unique_ptr<ROMol> dp_mol{nullptr};
  std::unique_ptr<ExplicitBitVect> dp_pattFP{nullptr};
  // SMILES strings of any connector regions.  Normally there will only
  // be 1 or 2.
  std::vector<std::shared_ptr<ROMol>> d_connRegions;
};

// Return a molecule containing the portions of the molecule starting at
// each dummy atom and going out up to 3 bonds.  There may be more than
// 1 fragment if there are dummy atoms more than 3 bonds apart, and there
// may be fragments with more than 1 dummy atom if their fragments fall
// within 3 bonds of each other.  E.g. the molecule [1*]CN(C[2*])Cc1ccccc1
// will give [1*]CN(C)C[1*].  The 2 dummy atoms are 4 bonds apart, but the
// fragments overlap.  All dummy atoms given isotope 1 whatever they had before.
RDKIT_SYNTHONSPACESEARCH_EXPORT std::unique_ptr<ROMol> getConnRegion(
    const ROMol &mol);

}  // namespace SynthonSpaceSearch
}  // namespace RDKit

#endif  // RDKIT_REAGENT_H
