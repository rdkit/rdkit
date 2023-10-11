//
//  Copyright (C) 2008-2022 Greg Landrum and other RDKit contributors
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
class Conformer;

namespace Chirality {

//! double bond stereo will be ignored/removed for rings smaller than this:
constexpr unsigned int minRingSizeForDoubleBondStereo = 8;

constexpr auto nonTetrahedralStereoEnvVar = "RDK_ENABLE_NONTETRAHEDRAL_STEREO";
constexpr auto useLegacyStereoEnvVar = "RDK_USE_LEGACY_STEREO_PERCEPTION";
constexpr bool nonTetrahedralStereoDefaultVal =
    true;  //!< whether or not nontetrahedral stereo is perceived by default
constexpr bool useLegacyStereoDefaultVal =
    true;  //!< whether or not the legacy stereo perception code is used by
           //!< default

RDKIT_GRAPHMOL_EXPORT extern void setAllowNontetrahedralChirality(bool val);
RDKIT_GRAPHMOL_EXPORT extern bool getAllowNontetrahedralChirality();

RDKIT_GRAPHMOL_EXPORT extern void setUseLegacyStereoPerception(bool val);
RDKIT_GRAPHMOL_EXPORT extern bool getUseLegacyStereoPerception();

RDKIT_GRAPHMOL_EXPORT extern bool
    useLegacyStereoPerception;  //!< Toggle usage of the legacy stereo
                                //!< perception code

RDKIT_GRAPHMOL_EXPORT extern bool
    useLegacyStereoPerception;  //!< Toggle usage of the legacy stereo
                                //!< perception code

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
  Atom_SquarePlanar,
  Atom_TrigonalBipyramidal,
  Atom_Octahedral,
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
  unsigned permutation = 0;  // for the non-tetrahedral stereo cases
  std::vector<unsigned> controllingAtoms;  // all atoms around the atom or bond.
  // Order is important
  bool operator==(const StereoInfo &other) const {
    return type == other.type && specified == other.specified &&
           centeredOn == other.centeredOn && descriptor == other.descriptor &&
           permutation == other.permutation &&
           controllingAtoms == other.controllingAtoms;
  }
  bool operator!=(const StereoInfo &other) const { return !(*this == other); }
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

//! removes atoms without specified chirality from stereo groups
RDKIT_GRAPHMOL_EXPORT void cleanupStereoGroups(ROMol &mol);

//! calls the approximate legacy code for assigning CIP labels
RDKIT_GRAPHMOL_EXPORT void assignLegacyCIPLabels(
    ROMol &mol, bool flagPossibleStereoCenters = false);

/// @cond
namespace detail {
RDKIT_GRAPHMOL_EXPORT bool isAtomPotentialNontetrahedralCenter(
    const Atom *atom);
RDKIT_GRAPHMOL_EXPORT bool isAtomPotentialTetrahedralCenter(const Atom *atom);
RDKIT_GRAPHMOL_EXPORT bool isAtomPotentialStereoAtom(const Atom *atom);
RDKIT_GRAPHMOL_EXPORT bool isBondPotentialStereoBond(const Bond *bond);
RDKIT_GRAPHMOL_EXPORT StereoInfo getStereoInfo(const Bond *bond);
RDKIT_GRAPHMOL_EXPORT StereoInfo getStereoInfo(const Atom *atom);
RDKIT_GRAPHMOL_EXPORT bool bondAffectsAtomChirality(const Bond *bond,
                                                    const Atom *atom);
RDKIT_GRAPHMOL_EXPORT unsigned int getAtomNonzeroDegree(const Atom *atom);

RDKIT_GRAPHMOL_EXPORT bool has_protium_neighbor(const ROMol &mol,
                                                const Atom *atom);

}  // namespace detail
/// @endcond

RDKIT_GRAPHMOL_EXPORT INT_VECT findStereoAtoms(const Bond *bond);

//! \name Non-tetrahedral stereochemistry
//! @{
RDKIT_GRAPHMOL_EXPORT bool hasNonTetrahedralStereo(const Atom *center);
RDKIT_GRAPHMOL_EXPORT Bond *getChiralAcrossBond(const Atom *center,
                                                const Bond *qry);
RDKIT_GRAPHMOL_EXPORT Bond *getChiralAcrossBond(const Atom *center,
                                                const Atom *qry);
RDKIT_GRAPHMOL_EXPORT Atom *getChiralAcrossAtom(const Atom *center,
                                                const Bond *qry);
RDKIT_GRAPHMOL_EXPORT Atom *getChiralAcrossAtom(const Atom *center,
                                                const Atom *qry);
//! \param which: if this is -1 then the second axial bond will be returned,
//! otherwise the first
RDKIT_GRAPHMOL_EXPORT Bond *getTrigonalBipyramidalAxialBond(const Atom *center,
                                                            int which = 0);
RDKIT_GRAPHMOL_EXPORT Atom *getTrigonalBipyramidalAxialAtom(const Atom *center,
                                                            int which = 0);

//! \returns 1 if it's the first axial atom, -1 if it's the second
RDKIT_GRAPHMOL_EXPORT int isTrigonalBipyramidalAxialBond(const Atom *center,
                                                         const Bond *qry);
RDKIT_GRAPHMOL_EXPORT int isTrigonalBipyramidalAxialAtom(const Atom *center,
                                                         const Atom *qry);

RDKIT_GRAPHMOL_EXPORT double getIdealAngleBetweenLigands(const Atom *center,
                                                         const Atom *lig1,
                                                         const Atom *lig2);

RDKIT_GRAPHMOL_EXPORT unsigned int getChiralPermutation(const Atom *center,
                                                        const INT_LIST &probe);
//! @}

RDKIT_GRAPHMOL_EXPORT std::ostream &operator<<(std::ostream &oss,
                                               const StereoSpecified &s);
RDKIT_GRAPHMOL_EXPORT std::ostream &operator<<(std::ostream &oss,
                                               const StereoType &s);

struct RDKIT_GRAPHMOL_EXPORT BondWedgingParameters {
  bool wedgeTwoBondsIfPossible =
      false;  //!< If this is enabled then two bonds will be wedged at chiral
              //!< centers subject to the following constraints:
              //!<   1. ring bonds will not be wedged
              //!<   2. bonds to chiral centers will not be wedged
              //!<   3. bonds separated by more than 120 degrees will not be
              //!<      wedged
};

namespace detail {
RDKIT_GRAPHMOL_EXPORT Bond::BondDir determineBondWedgeState(
    const Bond *bond, unsigned int fromAtomIdx, const Conformer *conf);
RDKIT_GRAPHMOL_EXPORT Bond::BondDir determineBondWedgeState(
    const Bond *bond, const INT_MAP_INT &wedgeBonds, const Conformer *conf);
RDKIT_GRAPHMOL_EXPORT std::pair<bool, INT_VECT> countChiralNbrs(
    const ROMol &mol, int noNbrs);
RDKIT_GRAPHMOL_EXPORT int pickBondToWedge(const Atom *atom, const ROMol &mol,
                                          const INT_VECT &nChiralNbrs,
                                          const INT_MAP_INT &resSoFar,
                                          int noNbrs);
}  // namespace detail

//! picks the bonds which should be wedged
/// \returns a map from bond idx -> controlling atom idx
RDKIT_GRAPHMOL_EXPORT INT_MAP_INT pickBondsToWedge(
    const ROMol &mol, const BondWedgingParameters *params = nullptr);

RDKIT_GRAPHMOL_EXPORT void wedgeMolBonds(
    ROMol &mol, const Conformer *conf = nullptr,
    const BondWedgingParameters *params = nullptr);
RDKIT_GRAPHMOL_EXPORT void wedgeBond(Bond *bond, unsigned int fromAtomIdx,
                                     const Conformer *conf);

}  // namespace Chirality
}  // namespace RDKit
#endif
