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

#include <RDGeneral/export.h>

namespace RDKit {
class Atom;
class ROMol;
class RWMol;

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

  const std::string &getSmiles() const { return d_smiles; }
  const std::string &getId() const { return d_id; }
  const std::unique_ptr<ROMol> &getOrigMol() const;
  const std::unique_ptr<ROMol> &getSearchMol() const;
  const std::unique_ptr<ExplicitBitVect> &getPattFP() const;
  const std::vector<std::shared_ptr<ROMol>> &getConnRegions() const;
  void setSearchMol(std::unique_ptr<RWMol> mol);

  // Writes to/reads from a binary stream.
  void writeToDBStream(std::ostream &os) const;
  void readFromDBStream(std::istream &is);

  // Tag each atom and bond with the given molecule number and its index,
  // so we can find them again in a product if necessary.
  void tagAtomsAndBonds(int molNum) const;

 private:
  std::string d_smiles;
  std::string d_id;

  // Keep 2 copies of the molecule.  The first is as passed in, which
  // will be used for building products.  The second will have its
  // atoms and bonds fiddled with to make them match the product (the
  // aliphatic precursor to aromatic product issue).  The search mol
  // doesn't always work with product building.
  std::unique_ptr<ROMol> dp_origMol{nullptr};
  std::unique_ptr<ROMol> dp_searchMol{nullptr};
  std::unique_ptr<ExplicitBitVect> dp_pattFP{nullptr};
  // SMILES strings of any connector regions.  Normally there will only
  // be 1 or 2.  These are derived from the search mol.
  std::vector<std::shared_ptr<ROMol>> d_connRegions;

  // One the search molecule has been added, get the connector regions,
  // connector fingerprint etc.
  void finishInitialization();
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
