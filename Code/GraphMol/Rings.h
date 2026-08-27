//
//  Copyright (C) 2004-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

//! \file Rings.h
//! \brief utility functionality for working with ring systems

#include <RDGeneral/export.h>
#ifndef RDRINGS_H
#define RDRINGS_H

#include <vector>
#include <map>
#include <boost/dynamic_bitset_fwd.hpp>

struct RDL_cycle;

namespace RDKit {
class ROMol;
};

namespace RingUtils {
typedef std::vector<int> INT_VECT;
typedef std::vector<std::vector<int>> VECT_INT_VECT;
typedef std::map<int, std::vector<int>> INT_INT_VECT_MAP;

//! Pick a set of rings that are fused together and contain a specified ring
/*!

   \param curr      the ID for the ring that should be in the fused system
   \param neighMap  adjacency lists for for all rings in the molecule.
           See documentation for makeNeighMap
   \param res       used to return the results: a list of rings that are fused
           with curr in them
   \param done      a bit vector recording the rings that are already dealt with
           this also can be used to avoid any rings that should not be included
           in the fused system
   \param depth retained for API compatibility; no longer used

*/
RDKIT_GRAPHMOL_EXPORT void pickFusedRings(int curr,
                                          const INT_INT_VECT_MAP &neighMap,
                                          INT_VECT &res,
                                          boost::dynamic_bitset<> &done,
                                          int depth = 0);

//! \brief For each ring in bring compute and store the ring that are fused
//! (share at least one bond with it).
/*!
  Useful both for the kekulization stuff and aromaticity perception.

  \param brings   list of rings - each ring is specified as a list of bond IDs
  \param neighMap an STL map into which the results are stored. Each entry in
  the
              map is indexed by the ring ID and the contents are the list
              rings (rather their IDs) that are fused with this ring
  \param maxSize if this is >0, rings that are larger than the threshold
                 will not be considered as candidates to be neighbors
  \param maxOverlapSize if this is >0, rings that overlap by more bonds than
                        this will not be considered to be neighbors

*/
RDKIT_GRAPHMOL_EXPORT void makeRingNeighborMap(const VECT_INT_VECT &brings,
                                               INT_INT_VECT_MAP &neighMap,
                                               unsigned int maxSize = 0,
                                               unsigned int maxOverlapSize = 0);

//! converts a list of atom indices into a list of bond indices
/*!

   \param res    list of ring - each ring is a list of atom ids
   \param brings reference to a list of rings to the write the results to
                 each ring here is list of bonds ids
   \param mol    the molecule of interest

  <b>Assumptions:</b>
   - each list of atom ids in "res" form a legitimate ring
   - each of these list of ordered such that a ring can be traversed
*/
RDKIT_GRAPHMOL_EXPORT void convertToBonds(const VECT_INT_VECT &res,
                                          VECT_INT_VECT &brings,
                                          const RDKit::ROMol &mol);

// normalizes a ring by rotating/reversing it so that the first atom
// is the one with the smallest index, and the second atom is the neighbor
// to the first one that again has the smallest index.
// This change should have a small performance footprint while it helps
// keeping test results consistent when making changes to ring detection.
void normalizeRing(std::vector<int> &ring);
std::vector<int> rdlCycleToAtomRing(RDL_cycle *cycle);

}  // namespace RingUtils

#endif
