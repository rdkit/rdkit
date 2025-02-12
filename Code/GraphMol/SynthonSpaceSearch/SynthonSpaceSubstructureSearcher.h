//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// This file declares a concrete class derived from SynthonSpaceSearcher
// that does substructure searching of the SynthonSpace.

#ifndef SYNTHONSPACESUBSTRUCTURESEARCHER_H
#define SYNTHONSPACESUBSTRUCTURESEARCHER_H

#include <RDGeneral/export.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearcher.h>

namespace RDKit::SynthonSpaceSearch {

class SynthonSpaceSubstructureSearcher : public SynthonSpaceSearcher {
 public:
  SynthonSpaceSubstructureSearcher() = delete;
  SynthonSpaceSubstructureSearcher(const ROMol &query,
                                   const SynthonSpaceSearchParams &params,
                                   SynthonSpace &space)
      : SynthonSpaceSearcher(query, params, space) {}

 private:
  // These are the pattern fingerprints for the fragments in this
  // search.  They are used for screening the fragments prior to
  // a substructure search.
  // It's thread-safe because each search creates its own Searcher
  // object so multiple searches in different threads will be in
  // different Searcher objects.  The fingerprints are keyed on
  // the address of the corresponding fragment.
  std::map<void *, std::unique_ptr<ExplicitBitVect>> d_pattFPs;
  // Likewise, the connector regions and connector region
  // fingerprints.
  std::map<void *, std::vector<std::unique_ptr<ROMol>>> d_connRegs;
  std::map<void *, std::vector<std::string>> d_connRegSmis;
  std::map<void *, std::vector<std::unique_ptr<ExplicitBitVect>>> d_connRegFPs;

  void extraSearchSetup(
      std::vector<std::vector<std::unique_ptr<ROMol>>> &fragSets) override;

  std::vector<std::unique_ptr<SynthonSpaceHitSet>> searchFragSet(
      std::vector<std::unique_ptr<ROMol>> &fragSet,
      const SynthonSet &reaction) const override;
  bool verifyHit(const ROMol &hit) const override;

  void getConnectorRegions(
      const std::vector<std::unique_ptr<ROMol>> &molFrags,
      std::vector<std::vector<ROMol *>> &connRegs,
      std::vector<std::vector<const std::string *>> &connRegSmis,
      std::vector<std::vector<ExplicitBitVect *>> &connRegFPs) const;
};

}  // namespace RDKit::SynthonSpaceSearch

#endif  // SYNTHONSPACESUBSTRUCTURESEARCHER_H
