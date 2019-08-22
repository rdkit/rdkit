//
//  Copyright (C) 2008 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file Chirality.h

  \brief Not intended for client-code use.

*/
#include <RDGeneral/export.h>
#ifndef _RD_CHIRALITY_20AUG2008_H_
#define _RD_CHIRALITY_20AUG2008_H_
#include <RDGeneral/types.h>
#include <GraphMol/Bond.h>

/// @cond
namespace RDKit {
class Atom;
class Bond;
class ROMol;

namespace Chirality {
/*!
  \param mol the molecule to be altered
  \param ranks  used to return the set of ranks.
                Should be at least mol.getNumAtoms() long.

  <b>Notes:</b>
     - All atoms gain a property common_properties::_CIPRank with their overall
       CIP ranking.

*/
RDKIT_GRAPHMOL_EXPORT void assignAtomCIPRanks(const ROMol &mol,
                                              UINT_VECT &ranks);

RDKIT_GRAPHMOL_EXPORT bool hasStereoBondDir(const Bond *bond);

/**
 *  Returns the first neighboring bond that can be found which has a stereo
 * bond direction set. If no such bond can be found, it returns null. No
 * checks are made to ensure there aren't any other conflicting directed bonds.
 */
RDKIT_GRAPHMOL_EXPORT const Bond *getNeighboringDirectedBond(const ROMol &mol,
                                                             const Atom *atom);

/**
 *  This just translates the labels, setting/translating StereoAtoms or the
 * label is not the responsibility of this function. If the passed label is not
 * E/Z, it will be returned unchanged.
 */
RDKIT_GRAPHMOL_EXPORT Bond::BondStereo translateEZLabelToCisTrans(
    Bond::BondStereo label);

}  // namespace Chirality
}  // namespace RDKit
/// @endcond
#endif
