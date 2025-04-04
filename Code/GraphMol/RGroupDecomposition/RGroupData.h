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
  RWMOL_SPTR combinedMol;
  std::vector<ROMOL_SPTR> mols;         // All the mols in the rgroup
  std::vector<std::string> smilesVect;  // used for rgroup equivalence
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

  void add(const ROMOL_SPTR &newMol,
           const std::vector<int> &rlabel_attachments);

  std::map<int, int> getNumBondsToRlabels() const;

  std::string toString() const;
  static std::string getRGroupLabel(int rlabel);
  static const std::string &getCoreLabel();
  static const std::string &getMolLabel();
  static bool isMolHydrogen(const ROMol &mol);

 private:
  void computeIsHydrogen();

  //! compute the canonical smiles for the attachments (bug: removes dupes since
  //! we are using a set...)
  std::string getSmiles() const;
  //! merge mol into combinedMol, including atom and bond highlights if present
  void mergeIntoCombinedMol(const ROMOL_SPTR &mol);
  std::map<int, std::vector<int>> rlabelAtomIndicesMap;
  std::map<int, std::vector<int>> rlabelBondIndicesMap;
};
}  // namespace RDKit

#endif
