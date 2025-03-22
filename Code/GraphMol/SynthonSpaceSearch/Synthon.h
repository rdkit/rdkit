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

// These are the numbers of bits used in the internal fingerprints.
// The user is not restricted to these numbers for the search.
inline constexpr unsigned int PATT_FP_NUM_BITS = 1024;

// This class holds a Synthon that will be part of a SynthonSet.
class RDKIT_SYNTHONSPACESEARCH_EXPORT Synthon {
 public:
  Synthon() = default;
  Synthon(const std::string &smi);
  Synthon(const Synthon &other);
  Synthon(Synthon &&other) = default;
  Synthon &operator=(const Synthon &other);
  Synthon &operator=(Synthon &&other) = default;

  const std::string &getSmiles() const { return d_smiles; }
  const std::unique_ptr<ROMol> &getOrigMol() const;
  const std::unique_ptr<ROMol> &getSearchMol() const;
  const std::unique_ptr<ExplicitBitVect> &getPattFP() const;
  const std::unique_ptr<ExplicitBitVect> &getFP() const;
  const std::vector<std::shared_ptr<ROMol>> &getConnRegions() const;
  void setSearchMol(std::unique_ptr<ROMol> mol);
  void setFP(std::unique_ptr<ExplicitBitVect> fp);

  // Writes to/reads from a binary stream.
  void writeToDBStream(std::ostream &os) const;
  void readFromDBStream(std::istream &is);

 private:
  std::string d_smiles;

  // Keep 2 copies of the molecule.  The first is as passed in, which
  // will be used for building products.  The second will have its
  // atoms and bonds fiddled with to make them match the product (the
  // aliphatic precursor to aromatic product issue).  The search mol
  // doesn't always work with product building.
  std::unique_ptr<ROMol> dp_origMol{nullptr};
  std::unique_ptr<ROMol> dp_searchMol{nullptr};
  // The pattern fingerprint, used for substructure search screening.
  std::unique_ptr<ExplicitBitVect> dp_pattFP{nullptr};
  // The fingerprint of the dp_searchMol, used in fingerprint similarity
  // searching.  Its type is known by the SynthonSpace that holds the
  // Synthon.
  std::unique_ptr<ExplicitBitVect> dp_FP{nullptr};
  // SMILES strings of any connector regions.  Normally there will only
  // be 1 or 2.  These are derived from the search mol.
  std::vector<std::shared_ptr<ROMol>> d_connRegions;

  // Once the search molecule has been added, get the connector regions,
  // connector fingerprint etc.
  void finishInitialization();
};

}  // namespace SynthonSpaceSearch
}  // namespace RDKit

#endif  // RDKIT_REAGENT_H
