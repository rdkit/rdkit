//
//  Copyright (C) 2025 Hussein Faara, Brian Kelley and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RD_SUBSET_H
#define RD_SUBSET_H

#include <RDGeneral/export.h>
#include <boost/dynamic_bitset.hpp>

#include <memory>
#include <vector>

namespace RDKit {
class RWMol;

//! Subsetting Methods for copyMolSubset
/*
 * These control what the path structures mean in copyMolSubset
 *
 *  NOTE: when the atoms and bonds are fully specified, the subset method is
 * ignored.
 */
enum class SubsetMethod {
  BONDS_BETWEEN_ATOMS =
      0,     //<! extract the specified atoms and the bonds between them
  BONDS = 1  //<! extract the specified bonds and their attached atoms
};

//! Subsetting Options for copyMolSubset
/*
 * These control what is copied over from the original molecule
 *
 *  NOTE: when the atoms and bonds are fully specified, the subset method is
 * ignored.
 */
struct RDKIT_GRAPHMOL_EXPORT SubsetOptions {
  bool sanitize =
      false; /**< perform sanitization automatically on the subset */
  bool clearComputedProps =
      false; /**< clear all computed props on the subsetted molecule */
  bool copyAsQuery = false; /**< Return the subset as a query */
  bool copyCoordinates =
      true; /**< *Copy the active coordinates from the molecule */
  unsigned int conformerIdx =
      -1u; /**< What conformer idx to use for the coordinates default is -1 */

  SubsetMethod method =
      SubsetMethod::BONDS_BETWEEN_ATOMS; /**< Subsetting method to use.  *Note*
                                          * if atoms and bonds are fully
                                          * specified the method is ignored.*/
};

struct RDKIT_GRAPHMOL_EXPORT SubsetInfo {
  boost::dynamic_bitset<> selectedAtoms; /**< atoms selected in the subset */
  boost::dynamic_bitset<> selectedBonds; /**< bonds selected in the subset */
  std::map<unsigned int, unsigned int>
      atomMapping; /**< mapping of original atom indices to new ones */
  std::map<unsigned int, unsigned int>
      bondMapping; /**< mapping of original bond indices to new ones */
};

//!
/*
 * Extract a subgraph from an ROMol. Bonds, atoms, substance groups and
 * stereo groups are only extracted to the subgraph if all participant entities
 * are contained within the given atoms and bonds.
 *
 * \param mol - starting mol
 * \param atoms - atoms to extract
 * \param bonds - bonds to extract
 * \param subsetInfo - optional subsetInfo to record the atoms and bonds used
 * \param options - subset options, note the method is ignored since all the
 * atoms and bonds are specified
 *
 *
 * NOTE: Bookmarks are currently copied, StereoGroups that are not entirely
 *        included in the subset are not copied.
 *
 * NOTE: when using this method the SubsetOptions.method is ignored
 *   for resolving the atoms and bonds to use to subset
 */

RDKIT_GRAPHMOL_EXPORT
std::unique_ptr<RDKit::RWMol> copyMolSubset(
    const RDKit::ROMol &mol, const std::vector<unsigned int> &atoms,
    const std::vector<unsigned int> &bonds,
    const SubsetOptions &options = SubsetOptions());

RDKIT_GRAPHMOL_EXPORT
std::unique_ptr<RDKit::RWMol> copyMolSubset(
    const RDKit::ROMol &mol, const std::vector<unsigned int> &atoms,
    const std::vector<unsigned int> &bonds, SubsetInfo &subsetInfo,
    const SubsetOptions &options = SubsetOptions());

//!
/*
 * Extract a subgraph from an ROMol. Bonds, substance groups and
 * stereo groups are only extracted to the subgraph if all participant entities
 * are selected by the `path` parameter.
 *
 * \param mol - starting mol
 * \param path - the indices of atoms or bonds to extract. If an index falls
 *             outside of the acceptable indices, it is ignored.
 *             use SubsetMethod::BONDS to indicate a bond path and
 *                 SubsetMethod::BONDS_BETWEEN_ATOMS to indicate an atom path
 *                 and any bond that includes both atoms in the path
 * \param subsetInfo - optional subsetInfo to record the atoms and bonds used
 * \param option - optional SubsetOptions to control how the subset is created
 *
 * NOTE: Bookmarks are currently copied, StereoGroups that are not entirely
 *        included in the subset are not copied.
 *
 */
RDKIT_GRAPHMOL_EXPORT std::unique_ptr<RDKit::RWMol> copyMolSubset(
    const RDKit::ROMol &mol, const std::vector<unsigned int> &path,
    const SubsetOptions &options = SubsetOptions());

RDKIT_GRAPHMOL_EXPORT std::unique_ptr<RDKit::RWMol> copyMolSubset(
    const RDKit::ROMol &mol, const std::vector<unsigned int> &path,
    SubsetInfo &subsetInfo, const SubsetOptions &options = SubsetOptions());

}  // namespace RDKit
#endif
