//
//  Copyright (C) 2019 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_SCAFFOLDNETWORK_H
#define RD_SCAFFOLDNETWORK_H

#include <vector>
#include <map>
#include <string>
#include <memory>

namespace RDKit {
class ROMol;
class ChemicalReaction;

namespace ScaffoldNetwork {

struct RDKIT_SCAFFOLDNETWORK_EXPORT ScaffoldNetworkParams {
  bool includeGenericScaffolds =
      true;  ///< include scaffolds with all atoms replaced by dummies
  bool includeGenericBondScaffolds =
      false;  ///< include scaffolds with all bonds replaced by single bonds
  bool includeScaffoldsWithoutAttachments =
      true;  ///< remove attachment points from scaffolds and include the result
  bool keepOnlyFirstFragment =
      true;  ///<  keep only the first fragment from the bond breaking rule
  bool pruneBeforeFragmenting =
      true;  ///<  Do a pruning/flattening step before starting fragmenting
  bool flattenIsotopes = true;  ///< remove isotopes when flattening
  bool flattenChirality =
      true;  ///< remove chirality and bond stereo when flattening
  bool flattenKeepLargest =
      true;  ///< keep only the largest fragment when doing flattening
  std::vector<std::shared_ptr<ChemicalReaction>>
      bondBreakersRxns;  ///< the reaction(s) used to fragment. Should expect a
                         ///< single reactant and produce two products
  ScaffoldNetworkParams()
      : ScaffoldNetworkParams{
            {"[!#0;R:1]-!@[!#0:2]>>[*:1]-[#0].[#0]-[*:2]"}} {};
  ScaffoldNetworkParams(const std::vector<std::string> &bondBreakersSmarts);
};

enum class RDKIT_SCAFFOLDNETWORK_EXPORT EdgeType {
  Fragment = 1,     ///< molecule -> fragment
  Generic = 2,      ///< molecule -> generic molecule (all atoms are dummies)
  GenericBond = 3,  ///< molecule -> generic bond molecule (all bonds single)
  RemoveAttachment = 4,  ///< molecule -> molecule with no attachment points
  Initialize = 5         ///< molecule -> flattened molecule
};

struct RDKIT_SCAFFOLDNETWORK_EXPORT NetworkEdge {
  size_t beginIdx;
  size_t endIdx;
  EdgeType type;
  NetworkEdge(size_t bi, size_t ei, EdgeType typ)
      : beginIdx(bi), endIdx(ei), type(typ){};
  bool operator==(const RDKit::ScaffoldNetwork::NetworkEdge &o) const {
    return (beginIdx == o.beginIdx) && (endIdx == o.endIdx) && (type == o.type);
  }
};

struct RDKIT_SCAFFOLDNETWORK_EXPORT ScaffoldNetwork {
  std::vector<std::string> nodes;  ///< SMILES for the scaffolds
  std::vector<unsigned>
      counts;  ///< number of times each scaffold was encountered
  std::vector<NetworkEdge> edges;  ///< edges in the network
};

//! update an existing ScaffoldNetwork using a set of molecules
template <typename T>
void updateScaffoldNetwork(const T &mols, ScaffoldNetwork &network,
                           const ScaffoldNetworkParams &params);

//! create a new ScaffoldNetwork for a set of molecules
template <typename T>
ScaffoldNetwork createScaffoldNetwork(const T &mols,
                                      const ScaffoldNetworkParams &params) {
  ScaffoldNetwork res;
  updateScaffoldNetwork(mols, res, params);
  return res;
}
}  // namespace ScaffoldNetwork
}  // namespace RDKit
#endif
