//
//  Copyright (C) 2019 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/ScaffoldNetwork/ScaffoldNetwork.h>

// declarations of stuff we want to test that isn't in the public API
namespace RDKit {
namespace ScaffoldNetwork {
namespace detail {
RDKIT_SCAFFOLDNETWORK_EXPORT std::vector<std::pair<std::string, ROMOL_SPTR>>
getMolFragments(const ROMol &mol, const ScaffoldNetworkParams &params);
RDKIT_SCAFFOLDNETWORK_EXPORT ROMol *makeScaffoldGeneric(const ROMol &mol,
                                                        bool doAtoms,
                                                        bool doBonds);
RDKIT_SCAFFOLDNETWORK_EXPORT ROMol *removeAttachmentPoints(
    const ROMol &mol, const ScaffoldNetworkParams &params);
RDKIT_SCAFFOLDNETWORK_EXPORT ROMol *pruneMol(
    const ROMol &mol, const ScaffoldNetworkParams &params);
RDKIT_SCAFFOLDNETWORK_EXPORT ROMol *flattenMol(
    const ROMol &mol, const ScaffoldNetworkParams &params);
RDKIT_SCAFFOLDNETWORK_EXPORT void addMolToNetwork(
    const ROMol &mol, ScaffoldNetwork &network,
    const ScaffoldNetworkParams &params);
}  // namespace detail
}  // namespace ScaffoldNetwork
}  // namespace RDKit
