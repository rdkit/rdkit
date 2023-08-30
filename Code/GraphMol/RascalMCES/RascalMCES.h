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
#ifndef RDKIT_RASCAL_MCES_H
#define RDKIT_RASCAL_MCES_H

#include <vector>

#include <GraphMol/RascalMCES/RascalClusterOptions.h>
#include <GraphMol/RascalMCES/RascalOptions.h>
#include <GraphMol/RascalMCES/RascalResult.h>
namespace RDKit {
class ROMol;

namespace RascalMCES {

// Find one or more MCESs between the two molecules.  The MCES is the
// Maximum Common Edge Substructure, and is the largest set of bonds
// common to the 2 molecules.
/*!
 *
 * @param mol1 : first molecule
 * @param mol2 : second molecule for MCES determination.
 * @param opts : (optional) set of options controlling the MCES determination
 * @return : vector of RascalResult objects.
 */
RDKIT_RASCALMCES_EXPORT std::vector<RascalResult> rascalMCES(
    const ROMol &mol1, const ROMol &mol2,
    const RascalOptions &opts = RascalOptions());

// Cluster the molecules using the Johnson similarity from rascalMCES
// and the algorithm of
// 'A Line Graph Algorithm for Clustering Chemical Structures Based
// on Common Substructural Cores', JW Raymond, PW Willett.
// https://match.pmf.kg.ac.rs/electronic_versions/Match48/match48_197-207.pdf
// https://eprints.whiterose.ac.uk/77598/
// This is a fuzzy clustering algorithm, so a molecule may appear in more than
// one cluster.  The final cluster is all the molecules that didn't fit into
// another cluster (the singletons).
/*!
 *
 * @param mols : molecules to cluster
 * @param clusOpts : (optional) cluster options
 * @return clusters as vector of vectors of unsigned ints - indices into the
 *         input mols vector
 */
RDKIT_RASCALMCES_EXPORT std::vector<std::vector<unsigned int>> rascalCluster(
    const std::vector<std::shared_ptr<ROMol>> &mols,
    const RascalClusterOptions &clusOpts = RascalClusterOptions());
// Cluster the molecules using the Johnson similarity from rascalMCES and
// the Butina algorithm.  Butina JCICS 39 747-750 (1999).
/*!
 *
 * @param mols : molecules to cluster
 * @param clusOpts : (optional) cluster options
 * @return clusters as vector of vectors of unsigned ints - indices into the
 *         input mols vector
 */
RDKIT_RASCALMCES_EXPORT std::vector<std::vector<unsigned int>>
rascalButinaCluster(
    const std::vector<std::shared_ptr<ROMol>> &mols,
    const RascalClusterOptions &clusOpts = RascalClusterOptions());
}  // namespace RascalMCES
}  // namespace RDKit
#endif  // RDKIT_RASCAL_MCES_H
