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

namespace RDKit {

#define RGROUPLABELS_ENUM_ITEMS                                         \
  RGD_ENUM_ITEM(IsotopeLabels, 1 << 0)                                  \
  RGD_ENUM_ITEM(AtomMapLabels, 1 << 1)                                  \
  RGD_ENUM_ITEM(AtomIndexLabels, 1 << 2)                                \
  RGD_ENUM_ITEM(RelabelDuplicateLabels, 1 << 3)                         \
  RGD_ENUM_ITEM(MDLRGroupLabels, 1 << 4)                                \
  RGD_ENUM_ITEM(DummyAtomLabels,                                        \
                1 << 5) /* These are rgroups but will get relabelled */ \
  RGD_ENUM_ITEM(AutoDetect, 0xFF)

#define RGROUPMATCHING_ENUM_ITEMS                                          \
  RGD_ENUM_ITEM(Greedy, 1 << 0)                                            \
  RGD_ENUM_ITEM(GreedyChunks, 1 << 1)                                      \
  RGD_ENUM_ITEM(Exhaustive, 1 << 2) /* not really useful for large sets */ \
  RGD_ENUM_ITEM(NoSymmetrization, 1 << 3)                                  \
  RGD_ENUM_ITEM(GA, 1 << 4)

#define RGROUPLABELLING_ENUM_ITEMS \
  RGD_ENUM_ITEM(AtomMap, 1 << 0)   \
  RGD_ENUM_ITEM(Isotope, 1 << 1)   \
  RGD_ENUM_ITEM(MDLRGroup, 1 << 2)

#define RGROUPCOREALIGNMENT_ENUM_ITEMS \
  RGD_ENUM_ITEM(NoAlignment, 0)        \
  RGD_ENUM_ITEM(MCS, 1 << 0)

#define RGROUPSCORE_ENUM_ITEMS \
  RGD_ENUM_ITEM(Match, 1 << 0) \
  RGD_ENUM_ITEM(FingerprintVariance, 1 << 2)

#define RGD_ENUM_ITEM(k, v) k = v,
typedef enum { RGROUPLABELS_ENUM_ITEMS } RGroupLabels;

typedef enum { RGROUPMATCHING_ENUM_ITEMS } RGroupMatching;

typedef enum { RGROUPLABELLING_ENUM_ITEMS } RGroupLabelling;

typedef enum { RGROUPCOREALIGNMENT_ENUM_ITEMS } RGroupCoreAlignment;

typedef enum { RGROUPSCORE_ENUM_ITEMS } RGroupScore;
#undef RGD_ENUM_ITEM
#define RGD_STD_MAP_ITEM(k) {#k, k},
#define RGD_ENUM_ITEM(k, v) RGD_STD_MAP_ITEM(k)

struct RDKIT_RGROUPDECOMPOSITION_EXPORT RGroupDecompositionParameters {
  unsigned int labels = AutoDetect;
  unsigned int matchingStrategy = GreedyChunks;
  unsigned int scoreMethod = Match;
  unsigned int rgroupLabelling = AtomMap | MDLRGroup;
  unsigned int alignment = MCS;

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

void updateRGroupDecompositionParametersFromJSON(
    RGroupDecompositionParameters &params, const std::string &details_json);
void updateRGroupDecompositionParametersFromJSON(
    RGroupDecompositionParameters &params, const char *details_json);

}  // namespace RDKit

#endif  // RDKIT_RGROUPDECOMPPARAMS_H
