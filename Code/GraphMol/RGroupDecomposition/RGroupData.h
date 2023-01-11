//
//  Copyright (C) 2017 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RGROUP_DATA
#define RGROUP_DATA

#include "../RDKitBase.h"
#include "RGroupUtils.h"
#include <DataStructs/ExplicitBitVect.h>
#include <set>
#include <vector>

namespace RDKit {

//! A single rgroup attached to a given core.
struct RDKIT_RGROUPDECOMPOSITION_EXPORT RGroupData {
  boost::shared_ptr<RWMol> combinedMol;
  std::vector<boost::shared_ptr<ROMol>> mols;  // All the mols in the rgroup
  std::vector<std::string> smilesVect;         // used for rgroup equivalence
  std::string
      smiles;  // smiles for all the mols in the rgroup (with attachments)
  std::set<int> attachments;  // core attachment points
  std::unique_ptr<ExplicitBitVect>
      fingerprint;  // fingerprint for score calculations
  std::vector<int> fingerprintOnBits;
  bool is_hydrogen = false;
  bool single_fragment = true;
  bool is_linker = false;
  bool labelled = false;

 public:
  RGroupData() {}

  void add(boost::shared_ptr<ROMol> newMol,
           const std::vector<int> &rlabel_attachments);

  std::map<int, int> getNumBondsToRlabels() const;

  std::string toString() const;

 private:
  void computeIsHydrogen();

  bool isMolHydrogen(ROMol &mol);

  //! compute the canonical smiles for the attachments (bug: removes dupes since
  //! we are using a set...)
  std::string getSmiles() const;
};
}  // namespace RDKit

#endif
