//
//  Copyright (C) 2017 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
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

typedef enum  {
  IsotopeLabels          = 0x01,
  AtomMapLabels          = 0x02,
  AtomIndexLabels        = 0x04,
  RelabelDuplicateLabels = 0x08,
  AutoDetect             = 0x0F,
} RGroupLabels;

typedef enum  {
  Greedy        = 0x01,
  GreedyChunks  = 0x02,
  Exhaustive    = 0x04, // not really useful for large sets
} RGroupMatching;

typedef enum  {
  AtomNum     = 0x01,
  Isotope     = 0x02,
  MDLRGroup   = 0x04,
} RGroupLabelling;

typedef enum  {
  None    = 0x0,
  MCS     = 0x01,
} RGroupCoreAlignment;

struct RGroupDecompositionParameters {
  RGroupLabels labels;
  RGroupMatching matchingStrategy;
  RGroupLabelling rgroupLabelling;
  RGroupCoreAlignment alignment;
  
  unsigned int chunkSize;
  bool onlyMatchAtRGroups;
  bool removeAllHydrogenRGroups;

 RGroupDecompositionParameters(RGroupLabels labels=AutoDetect,
                               RGroupMatching strategy=GreedyChunks,
                               RGroupLabelling labelling=MDLRGroup,
                               RGroupCoreAlignment alignment=MCS,
                               unsigned int chunkSize=5,
                               bool matchOnlyAtRGroups=false,
                               bool removeHydrogenOnlyGroups=true) :
  labels(labels),
    matchingStrategy(strategy),
    rgroupLabelling(labelling),
    alignment(alignment),
    chunkSize(chunkSize),
    onlyMatchAtRGroups(matchOnlyAtRGroups),
    removeAllHydrogenRGroups(removeHydrogenOnlyGroups),
    indexOffset(-1) {
  }
  bool prepareCore(RWMol &, const RWMol *alignCore);
  void SetRGroupLabels(RGroupLabels rgroupLabels) {
    labels = rgroupLabels;
  }
  void SetRGroupLabelling(RGroupLabelling labelling) {
    rgroupLabelling = labelling;
  }
  void SetRGroupMatching(RGroupMatching matching) {
    matchingStrategy = matching;
  }
  void SetRGroupCoreAlignment(RGroupCoreAlignment coreAlignment) {
    alignment = coreAlignment;
  }
  void SetChunkSize(unsigned int size) { chunkSize = size; }
  void SetOnlyMatchAtRGroups(bool matchAtRGroups) { onlyMatchAtRGroups = matchAtRGroups; }
  void SetRemoveRGroupsThatAreAllHydrogen( bool remove ) { removeAllHydrogenRGroups = remove; }

  RGroupLabels GetRGroupLabels() const { return labels; }
  RGroupMatching GetRGroupLabelling() const { return matchingStrategy; }
  RGroupLabelling GetRGroupMatching() const { return rgroupLabelling; }
  RGroupCoreAlignment GetRGroupCoreAlignment() const { return alignment; }
  unsigned int  GetChunkSize() const { return chunkSize; }
  bool GetOnlyMatchAtRGroups() const { return onlyMatchAtRGroups; }
  bool GetRemoveRGroupsThatAreAllHydrogen() const { return removeAllHydrogenRGroups; }

  

 private:
  int indexOffset;
};

typedef std::map<std::string,  boost::shared_ptr<ROMol> >  RGroupRow;
typedef std::vector<boost::shared_ptr<ROMol> > RGroupColumn;

typedef std::vector<RGroupRow> RGroupRows;
typedef std::map<std::string, RGroupColumn> RGroupColumns;

struct RGroupDecompData;
class RGroupDecomposition {
  RGroupDecompData *data; // implementation details
  RGroupDecomposition(const RGroupDecomposition &);             // no copy construct
  RGroupDecomposition& operator=(const RGroupDecomposition&);   // Prevent assignment

public:
  RGroupDecomposition(const ROMol &core,
                      const RGroupDecompositionParameters &params=RGroupDecompositionParameters());
  RGroupDecomposition(const std::vector<ROMOL_SPTR> &cores,
                      const RGroupDecompositionParameters &params=RGroupDecompositionParameters());
  
 ~RGroupDecomposition();
 
  int  add(const ROMol &mol);
  bool process();
    
  //! return rgroups in row order group[row][attachment_point] = ROMol
  RGroupRows getRGroupsAsRows() const;
  //! return rgroups in column order group[attachment_point][row] = ROMol
  RGroupColumns getRGroupsAsColumns() const;
};

}

#endif
