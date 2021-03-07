//
//  Copyright (C) 2017 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RDKIT_RGROUPDECOMP_H
#define RDKIT_RGROUPDECOMP_H

#include "../RDKitBase.h"
#include <GraphMol/Substruct/SubstructMatch.h>
#include <chrono>

namespace RDKit {

//! Compute the isomorphic degenerative points in the
//!  molecule.  These are points that are symmetrically
//!  equivalent.
/*!
   \param mol     Molecule to compute the degenerative points

   \return the set of degenerative points (set<unsigned int>)
*/

typedef enum {
  IsotopeLabels = 0x01,
  AtomMapLabels = 0x02,
  AtomIndexLabels = 0x04,
  RelabelDuplicateLabels = 0x08,
  MDLRGroupLabels = 0x10,
  DummyAtomLabels = 0x20,  // These are rgroups but will get relabelled
  AutoDetect = 0xFF,
} RGroupLabels;

typedef enum {
  Greedy = 0x01,
  GreedyChunks = 0x02,
  Exhaustive = 0x04,  // not really useful for large sets
  NoSymmetrization = 0x08,
  GA = 0x10,
} RGroupMatching;

typedef enum {
  AtomMap = 0x01,
  Isotope = 0x02,
  MDLRGroup = 0x04,
} RGroupLabelling;

typedef enum {
  // DEPRECATED, remove the following line in release 2021.03
  None = 0x0,
  NoAlignment = 0x0,
  MCS = 0x01,
} RGroupCoreAlignment;

typedef enum {
  Match = 0x1,
  FingerprintVariance = 0x4,
} RGroupScore;

struct RDKIT_RGROUPDECOMPOSITION_EXPORT RGroupDecompositionProcessResult {
  const bool success;
  const double score;
  RGroupDecompositionProcessResult(const bool success, const double score)
      : success(success), score(score) {}
};

class RGroupMatch;

struct RDKIT_RGROUPDECOMPOSITION_EXPORT RGroupDecompositionParameters {
  unsigned int labels = AutoDetect;
  unsigned int matchingStrategy = GreedyChunks;
  unsigned int scoreMethod = Match;
  unsigned int rgroupLabelling = AtomMap | MDLRGroup;
  unsigned int alignment = MCS;

  unsigned int chunkSize = 5;
  bool onlyMatchAtRGroups = false;
  bool removeAllHydrogenRGroups = true;
  bool removeHydrogensPostMatch = true;
  double timeout = -1.0;  ///< timeout in seconds. <=0 indicates no timeout

  // Determine how to assign the rgroup labels from the given core
  unsigned int autoGetLabels(const RWMol &);

  // Prepare the core for substructure searching and rgroup assignment
  bool prepareCore(RWMol &, const RWMol *alignCore);

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
  bool gaParallelRuns = true;

 private:
  int indexOffset{-1};
};

typedef std::map<std::string, ROMOL_SPTR> RGroupRow;
typedef std::vector<ROMOL_SPTR> RGroupColumn;

typedef std::vector<RGroupRow> RGroupRows;
typedef std::map<std::string, RGroupColumn> RGroupColumns;

struct RGroupDecompData;
class RDKIT_RGROUPDECOMPOSITION_EXPORT RGroupDecomposition {
 private:
  RGroupDecompData *data;                            // implementation details
  RGroupDecomposition(const RGroupDecomposition &);  // no copy construct
  RGroupDecomposition &operator=(
      const RGroupDecomposition &);  // Prevent assignment
  RWMOL_SPTR outputCoreMolecule(const RGroupMatch &match, const std::map<int, bool> &usedRGroupMap = std::map<int, bool>()) const;
  std::map<int, bool> getBlankRGroupMap() const;

 public:
  RGroupDecomposition(const ROMol &core,
                      const RGroupDecompositionParameters &params =
                          RGroupDecompositionParameters());
  RGroupDecomposition(const std::vector<ROMOL_SPTR> &cores,
                      const RGroupDecompositionParameters &params =
                          RGroupDecompositionParameters());

  ~RGroupDecomposition();

  int add(const ROMol &mol);
  RGroupDecompositionProcessResult processAndScore();
  bool process();

  const RGroupDecompositionParameters &params() const;
  //! return the current group labels
  std::vector<std::string> getRGroupLabels() const;

  //! return rgroups in row order group[row][attachment_point] = ROMol
  RGroupRows getRGroupsAsRows() const;
  //! return rgroups in column order group[attachment_point][row] = ROMol
  RGroupColumns getRGroupsAsColumns() const;
};

RDKIT_RGROUPDECOMPOSITION_EXPORT unsigned int RGroupDecompose(
    const std::vector<ROMOL_SPTR> &cores, const std::vector<ROMOL_SPTR> &mols,
    RGroupRows &rows, std::vector<unsigned int> *unmatched = nullptr,
    const RGroupDecompositionParameters &options =
        RGroupDecompositionParameters());

RDKIT_RGROUPDECOMPOSITION_EXPORT unsigned int RGroupDecompose(
    const std::vector<ROMOL_SPTR> &cores, const std::vector<ROMOL_SPTR> &mols,
    RGroupColumns &columns, std::vector<unsigned int> *unmatched = nullptr,
    const RGroupDecompositionParameters &options =
        RGroupDecompositionParameters());

inline bool checkForTimeout(const std::chrono::steady_clock::time_point &t0,
                            double timeout, bool throwOnTimeout = true) {
  if (timeout <= 0) return false;
  auto t1 = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed = t1 - t0;
  if (elapsed.count() >= timeout) {
    if (throwOnTimeout) {
      throw std::runtime_error("operation timed out");
    }
    return true;
  }
  return false;
}

}  // namespace RDKit

#endif
