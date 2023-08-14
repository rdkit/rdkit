//
// Copyright (C) David Cosgrove 2023
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef RDKIT_RASCAL_DETAILS_H
#define RDKIT_RASCAL_DETAILS_H

#include <map>

#include <GraphMol/RascalMCES/RascalOptions.h>
#include <GraphMol/RascalMCES/RascalResult.h>
namespace RDKit {
class ROMol;

namespace RascalMCES {

namespace details {

struct ClusNode {
  std::shared_ptr<RascalResult> d_res;
  double d_sim;
  unsigned int d_mol1Num, d_mol2Num;
};

RDKIT_RASCALMCES_EXPORT double tier1Sim(
    const RDKit::ROMol &mol1, const RDKit::ROMol &mol2,
    std::map<int, std::vector<std::pair<int, int>>> &degSeqs1,
    std::map<int, std::vector<std::pair<int, int>>> &degSeqs2);

RDKIT_RASCALMCES_EXPORT double tier2Sim(
    const ROMol &mol1, const ROMol &mol2,
    const std::map<int, std::vector<std::pair<int, int>>> &degSeqs1,
    const std::map<int, std::vector<std::pair<int, int>>> &degSeqs2,
    const std::vector<unsigned int> &bondLabels1,
    const std::vector<unsigned int> &bondLabels2);

RDKIT_RASCALMCES_EXPORT void getBondLabels(
    const RDKit::ROMol &mol1, const RDKit::ROMol &mol2,
    const RascalOptions &opts, std::vector<unsigned int> &bondLabels1,
    std::vector<unsigned int> &bondLabels2);

std::vector<std::vector<ClusNode>> buildProximityGraph(
    const std::vector<std::shared_ptr<ROMol>> &mols,
    const RascalClusterOptions &clusOpts);
}  // namespace details

}  // namespace RascalMCES
}  // namespace RDKit
#endif  // RDKIT_RASCAL_MCES_H
