//
//  Copyright (C) 2020-2021 Novartis Institutes for BioMedical Research and
//  other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RGROUP_CORE
#define RGROUP_CORE

#include <GraphMol/SmilesParse/SmartsWrite.h>
#include "../RDKitBase.h"
#include "RGroupUtils.h"
#include "GraphMol/Substruct/SubstructMatch.h"

// #define VERBOSE 1

namespace RDKit {

//! RCore is the core common to a series of molecules
struct RCore {
  boost::shared_ptr<RWMol> core;
  // core with terminal user R groups stripped for matching
  boost::shared_ptr<RWMol> matchingMol;
  boost::shared_ptr<RWMol> labelledCore;

  // Bitset: indices corresponding to atoms bearing user-defined labels are 1
  boost::dynamic_bitset<> core_atoms_with_user_labels;
  // Number of user labelled rgroups in the core
  size_t numberUserRGroups = 0;
  RCore() {}
  RCore(const RWMol &c) : core(new RWMol(c)) { init(); }

  void init();

  inline bool isCoreAtomUserLabelled(int idx) const {
    return core_atoms_with_user_labels.test(idx);
  }

  void countUserRGroups() {
    numberUserRGroups = core_atoms_with_user_labels.count();
  }

  // Find all the core atoms that have user
  // label and set their indices to 1 into core_atoms_with_user_labels
  void findIndicesWithRLabel();

  // Return a copy of core where dummy atoms are replaced by
  // the respective matching atom in mol, while other atoms have
  // their aromatic flag and formal charge copied from 
  // the respective matching atom in mol
  ROMOL_SPTR replaceCoreAtomsWithMolMatches(bool &hasCoreDummies,
                                            const ROMol &mol,
                                            const MatchVectType &match) const;

  // Final core returned to user with dummy atoms and bonds set to those in the
  // match
  RWMOL_SPTR coreWithMatches(const ROMol &coreReplacedAtoms) const;

  std::vector<MatchVectType> matchTerminalUserRGroups(
      const RWMol &target, MatchVectType match,
      const SubstructMatchParameters &sssParams) const;

  inline bool isTerminalRGroupWithUserLabel(const int idx) const {
    return terminalRGroupDummyAtoms.find(idx) != terminalRGroupDummyAtoms.end();
  }

  /*
   * For when onlyMatchAtRGroups = true.  Checks the query core can satisfy all
   * attachment points. Including when two user defined attachment points can
   * match the same target atom.
   */
  bool checkAllBondsToAttachmentPointPresent(
      const ROMol &mol, const int attachmentIdx,
      const MatchVectType &mapping) const;

 private:
  // The set of atom indices in the core for terminal R groups with user label
  std::set<int> terminalRGroupDummyAtoms;
  // The set of atom indices in the core for terminal R groups with atom indices with or without user labels
  std::set<int> terminalRGroupAtoms;
  // An atom index map of terminal R groups to their heavy atom neighbor
  std::map<int, int> terminalRGroupAtomToNeighbor;

  void replaceCoreAtom(RWMol &mol, Atom &atom, const Atom &other) const;

  // Convert a matching molecule index to a core index
  int matchingIndexToCoreIndex(int matchingIndex) const;

  // Build the matching molecule (core minus user R groups)
  void buildMatchingMol();
  
  // Add attachement points to unlabelled R Groups
  void addDummyAtomsToUnlabelledCoreAtoms();
};

}  // namespace RDKit
#endif
