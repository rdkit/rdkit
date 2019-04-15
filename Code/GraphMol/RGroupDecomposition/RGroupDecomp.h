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

namespace RDKit {

//! Compute the isomporphic degenerative points in the
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
  AutoDetect = 0x0F,
} RGroupLabels;

typedef enum {
  Greedy = 0x01,
  GreedyChunks = 0x02,
  Exhaustive = 0x04,  // not really useful for large sets
} RGroupMatching;

typedef enum {
  AtomMap = 0x01,
  Isotope = 0x02,
  MDLRGroup = 0x04,
} RGroupLabelling;

typedef enum {
  None = 0x0,
  MCS = 0x01,
} RGroupCoreAlignment;

struct RDKIT_RGROUPDECOMPOSITION_EXPORT RGroupDecompositionParameters {
  unsigned int labels;
  unsigned int matchingStrategy;
  unsigned int rgroupLabelling;
  unsigned int alignment;

  unsigned int chunkSize;
  bool onlyMatchAtRGroups;
  bool removeAllHydrogenRGroups;
  bool removeHydrogensPostMatch;

  RGroupDecompositionParameters(unsigned int labels = AutoDetect,
                                unsigned int strategy = GreedyChunks,
                                unsigned int labelling = AtomMap | MDLRGroup,
                                unsigned int alignment = MCS,
                                unsigned int chunkSize = 5,
                                bool matchOnlyAtRGroups = false,
                                bool removeHydrogenOnlyGroups = true,
                                bool removeHydrogensPostMatch = false)
      : labels(labels),
        matchingStrategy(strategy),
        rgroupLabelling(labelling),
        alignment(alignment),
        chunkSize(chunkSize),
        onlyMatchAtRGroups(matchOnlyAtRGroups),
        removeAllHydrogenRGroups(removeHydrogenOnlyGroups),
        removeHydrogensPostMatch(removeHydrogensPostMatch),
        indexOffset(-1) {}
  bool prepareCore(RWMol &, const RWMol *alignCore);

 private:
  int indexOffset;
};

typedef std::map<std::string, boost::shared_ptr<ROMol>> RGroupRow;
typedef std::vector<boost::shared_ptr<ROMol>> RGroupColumn;

typedef std::vector<RGroupRow> RGroupRows;
typedef std::map<std::string, RGroupColumn> RGroupColumns;

struct RGroupDecompData;
class RDKIT_RGROUPDECOMPOSITION_EXPORT RGroupDecomposition {
  RGroupDecompData *data;                            // implementation details
  RGroupDecomposition(const RGroupDecomposition &);  // no copy construct
  RGroupDecomposition &operator=(
      const RGroupDecomposition &);  // Prevent assignment

 public:
  RGroupDecomposition(const ROMol &core,
                      const RGroupDecompositionParameters &params =
                          RGroupDecompositionParameters());
  RGroupDecomposition(const std::vector<ROMOL_SPTR> &cores,
                      const RGroupDecompositionParameters &params =
                          RGroupDecompositionParameters());

  ~RGroupDecomposition();

  int add(const ROMol &mol);
  bool process();

  //! return rgroups in row order group[row][attachment_point] = ROMol
  RGroupRows getRGroupsAsRows() const;
  //! return rgroups in column order group[attachment_point][row] = ROMol
  RGroupColumns getRGroupsAsColumns() const;
};

RDKIT_RGROUPDECOMPOSITION_EXPORT unsigned int RGroupDecompose(
    const std::vector<ROMOL_SPTR> &cores, const std::vector<ROMOL_SPTR> &mols,
    RGroupRows &rows, std::vector<unsigned int> *unmatched = 0,
    const RGroupDecompositionParameters &options =
        RGroupDecompositionParameters());

RDKIT_RGROUPDECOMPOSITION_EXPORT unsigned int RGroupDecompose(
    const std::vector<ROMOL_SPTR> &cores, const std::vector<ROMOL_SPTR> &mols,
    RGroupColumns &columns, std::vector<unsigned int> *unmatched = 0,
    const RGroupDecompositionParameters &options =
        RGroupDecompositionParameters());
}  // namespace RDKit

#endif
