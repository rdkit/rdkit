//
//  Copyright (C) 2008-2020 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file Chirality.h

*/
#include <RDGeneral/export.h>
#ifndef RD_CHIRALITY_20AUG2008_H
#define RD_CHIRALITY_20AUG2008_H
#include <RDGeneral/types.h>
#include <GraphMol/Bond.h>
#include <boost/dynamic_bitset.hpp>
#include <limits>

namespace RDKit {
class Atom;
class Bond;
class ROMol;

namespace Chirality {
/// @cond
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
/// @endcond

enum class StereoType {
  Unspecified,
  Atom_Tetrahedral,
  Bond_Double,         // single double bond and odd-numbered cumulenes
  Bond_Cumulene_Even,  // even-numbered cumulenes
  Bond_Atropisomer
};

enum class StereoDescriptor { None, Tet_CW, Tet_CCW, Bond_Cis, Bond_Trans };

enum class StereoSpecified {
  Unspecified,  // no information provided
  Specified,
  Unknown  // deliberately marked as unknown
};

struct RDKIT_GRAPHMOL_EXPORT StereoInfo {
  // REVIEW: absolute stereo data member?
#ifdef _MSC_VER
  static const unsigned NOATOM =
      std::numeric_limits<unsigned>::max();  // used to mark missing atoms
#else
  static const unsigned NOATOM;  // used to mark missing atoms
#endif
  StereoType type = StereoType::Unspecified;
  StereoSpecified specified = StereoSpecified::Unspecified;
  unsigned centeredOn = NOATOM;
  StereoDescriptor descriptor = StereoDescriptor::None;
  std::vector<unsigned> controllingAtoms;  // all atoms around the atom or bond.
  // Order is important
  bool operator==(const StereoInfo &other) const {
    return type == other.type && specified == other.specified &&
           centeredOn == other.centeredOn && descriptor == other.descriptor &&
           controllingAtoms == other.controllingAtoms;
  }
};

//! identifies potential stereoatoms and stereobonds in a molecule
/*!
  Note that this function is still somewhat experimental and the API
  and results may change in a future release.

  \param mol the molecule to look for stereo in
  \param cleanIt remove chirality/stereo specifications from atoms/bonds that
     cannot be chiral/stereo
*/
RDKIT_GRAPHMOL_EXPORT std::vector<StereoInfo> findPotentialStereo(
    ROMol &mol, bool cleanIt, bool flagPossible = true);
//! overload
RDKIT_GRAPHMOL_EXPORT std::vector<StereoInfo> findPotentialStereo(
    const ROMol &mol);

/// @cond
namespace detail {
RDKIT_GRAPHMOL_EXPORT bool isAtomPotentialTetrahedralCenter(const Atom *atom);
RDKIT_GRAPHMOL_EXPORT bool isAtomPotentialStereoAtom(const Atom *atom);
RDKIT_GRAPHMOL_EXPORT bool isBondPotentialStereoBond(const Bond *bond);
RDKIT_GRAPHMOL_EXPORT StereoInfo getStereoInfo(const Bond *bond);
RDKIT_GRAPHMOL_EXPORT StereoInfo getStereoInfo(const Atom *atom);

}  // namespace detail
/// @endcond

RDKIT_GRAPHMOL_EXPORT INT_VECT findStereoAtoms(const Bond *bond);

}  // namespace Chirality
}  // namespace RDKit
#endif
