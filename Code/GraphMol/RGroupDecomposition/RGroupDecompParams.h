//
//  Copyright (c) 2017-2023, Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RDKIT_RGROUPDECOMPPARAMS_H
#define RDKIT_RGROUPDECOMPPARAMS_H

#include "../RDKitBase.h"
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/BetterEnums.h>

namespace RDKit {

BETTER_ENUM(RGroupLabels, unsigned int,
  IsotopeLabels = 0x01,
  AtomMapLabels = 0x02,
  AtomIndexLabels = 0x04,
  RelabelDuplicateLabels = 0x08,
  MDLRGroupLabels = 0x10,
  DummyAtomLabels = 0x20,  // These are rgroups but will get relabelled
  AutoDetect = 0xFF
);

BETTER_ENUM(RGroupMatching, unsigned int,
  Greedy = 0x01,
  GreedyChunks = 0x02,
  Exhaustive = 0x04,  // not really useful for large sets
  NoSymmetrization = 0x08,
  GA = 0x10
);

BETTER_ENUM(
  RGroupLabelling, unsigned int,
  AtomMap = 0x01,
  Isotope = 0x02,
  MDLRGroup = 0x04
);

BETTER_ENUM(RGroupCoreAlignment, unsigned int,
  NoAlignment = 0x0,
  MCS = 0x01
);

BETTER_ENUM(RGroupScore, unsigned int,
  Match = 0x1,
  FingerprintVariance = 0x4
);

struct RDKIT_RGROUPDECOMPOSITION_EXPORT RGroupDecompositionParameters {
  unsigned int labels = RGroupLabels::AutoDetect;
  unsigned int matchingStrategy = RGroupMatching::GreedyChunks;
  unsigned int scoreMethod = RGroupScore::Match;
  unsigned int rgroupLabelling =
      RGroupLabelling::AtomMap | RGroupLabelling::MDLRGroup;
  unsigned int alignment = RGroupCoreAlignment::MCS;

  unsigned int chunkSize = 5;
  //! only allow rgroup decomposition at the specified rgroups
  bool onlyMatchAtRGroups = false;
  //! remove all user-defined rgroups that only have hydrogens
  bool removeAllHydrogenRGroups = true;
  //! remove all user-defined rgroups that only have hydrogens,
  //! and also remove the corresponding labels from the core
  bool removeAllHydrogenRGroupsAndLabels = true;
  //! remove all hydrogens from the output molecules
  bool removeHydrogensPostMatch = true;
  //! allow labelled Rgroups of degree 2 or more
  bool allowNonTerminalRGroups = false;
  //! unlabelled core atoms can have multiple rgroups
  bool allowMultipleRGroupsOnUnlabelled = false;
  // extended query settings for core matching
  bool doTautomers = false;
  bool doEnumeration = false;
  //! include target molecule (featuring explicit hydrogens where they
  //! coincide with R groups in the core) into RGD results,
  //! and set _rgroupTargetAtoms and _rgroupTargetBonds properties
  //! on R groups and core as vectors of target atom and bond indices
  //! to enable highlighting for SAR analysis (see
  //! https://greglandrum.github.io/rdkit-blog/posts/2021-08-07-rgd-and-highlighting.html)
  bool includeTargetMolInResults = false;

  double timeout = -1.0;  ///< timeout in seconds. <=0 indicates no timeout

  // Determine how to assign the rgroup labels from the given core
  unsigned int autoGetLabels(const RWMol &);

  // Prepare the core for substructure searching and rgroup assignment
  bool prepareCore(RWMol &, const RWMol *alignCore);

  // Add r groups to unlabelled atoms if allowMultipleRGroupsOnUnlabelled is set
  void addDummyAtomsToUnlabelledCoreAtoms(RWMol &core);

  // Parameters specific to GA

  // GA population size or -1 to use best guess
  int gaPopulationSize = -1;
  // GA maximum number of operations or -1 to use best guess
  int gaMaximumOperations = -1;
  // GA number of operations permitted without improvement before exiting (-1
  // for best guess)
  int gaNumberOperationsWithoutImprovement = -1;
  // GA random number seed (-1 for default, -2 for random seed)
  int gaRandomSeed = -1;
  // Number of runs
  int gaNumberRuns = 1;
  // Sequential or parallel runs?
#ifdef RDK_BUILD_THREADSAFE_SSS
  bool gaParallelRuns = true;
#else
  bool gaParallelRuns = false;
#endif
  // Controls the way substructure matching with the core is done
  SubstructMatchParameters substructmatchParams;

  RGroupDecompositionParameters() { substructmatchParams.useChirality = true; }

 private:
  int indexOffset{-1};
  void checkNonTerminal(const Atom &atom) const;
};

}  // namespace RDKit

#endif  // RDKIT_RGROUPDECOMPPARAMS_H
