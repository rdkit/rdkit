//
//  Copyright (c) 2017-2021, Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
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
#include "RGroupDecompParams.h"
#include <GraphMol/Substruct/SubstructMatch.h>
#include <chrono>

namespace RDKit {

struct RDKIT_RGROUPDECOMPOSITION_EXPORT RGroupDecompositionProcessResult {
  const bool success;
  const double score;
  RGroupDecompositionProcessResult(const bool success, const double score)
      : success(success), score(score) {}
};

struct RGroupMatch;

typedef std::map<std::string, ROMOL_SPTR> RGroupRow;
typedef std::vector<ROMOL_SPTR> RGroupColumn;

typedef std::vector<RGroupRow> RGroupRows;
typedef std::map<std::string, RGroupColumn> RGroupColumns;

class UsedLabelMap {
 public:
  UsedLabelMap(const std::map<int, int> &mapping) {
    for (const auto &rl : mapping) {
      d_map[rl.second] = std::make_pair(false, (rl.first > 0));
    }
  }
  bool has(int label) const { return d_map.find(label) != d_map.end(); }
  bool getIsUsed(int label) const { return d_map.at(label).first; }
  void setIsUsed(int label) { d_map[label].first = true; }
  bool isUserDefined(int label) const { return d_map.at(label).second; }

 private:
  std::map<int, std::pair<bool, bool>> d_map;
};

struct RGroupDecompData;
class RDKIT_RGROUPDECOMPOSITION_EXPORT RGroupDecomposition {
 private:
  RGroupDecompData *data;                            // implementation details
  RGroupDecomposition(const RGroupDecomposition &);  // no copy construct
  RGroupDecomposition &operator=(
      const RGroupDecomposition &);  // Prevent assignment
  RWMOL_SPTR outputCoreMolecule(const RGroupMatch &match,
                                const UsedLabelMap &usedRGroupMap) const;

 public:
  RGroupDecomposition(const ROMol &core,
                      const RGroupDecompositionParameters &params =
                          RGroupDecompositionParameters());
  RGroupDecomposition(const std::vector<ROMOL_SPTR> &cores,
                      const RGroupDecompositionParameters &params =
                          RGroupDecompositionParameters());

  ~RGroupDecomposition();

  //! Returns the index of the added molecule in the RGroupDecomposition
  ///  or a negative error code
  /// :param mol: Molecule to add to the decomposition
  /// :result: index of the molecle or
  ///             -1 if none of the core matches
  ///             -2 if the matched molecule has no sidechains, i.e. is the
  ///                same as the scaffold
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
  if (timeout <= 0) {
    return false;
  }
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
