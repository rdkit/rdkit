//
//  Copyright (C) 2018 Boran Adas, Google Summer of Code
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef RD_FINGERPRINTUTIL_H_2018_07
#define RD_FINGERPRINTUTIL_H_2018_07

#include <GraphMol/RDKitBase.h>
#include <DataStructs/SparseIntVect.h>
#include <DataStructs/BitVects.h>
#include <cstdint>
#include <tuple>
#include <vector>
#include <map>
#include <DataStructs/ExplicitBitVect.h>

#include <GraphMol/Subgraphs/Subgraphs.h>

namespace RDKit {
namespace AtomPairs {
const unsigned int numTypeBits = 4;
const unsigned int atomNumberTypes[1 << numTypeBits] = {
    5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 51, 52, 43};
const unsigned int numPiBits = 2;
const unsigned int maxNumPi = (1 << numPiBits) - 1;
const unsigned int numBranchBits = 3;
const unsigned int maxNumBranches = (1 << numBranchBits) - 1;
const unsigned int numChiralBits = 2;
const unsigned int codeSize = numTypeBits + numPiBits + numBranchBits;
const unsigned int numPathBits = 5;
const unsigned int maxPathLen = (1 << numPathBits) - 1;
const unsigned int numAtomPairFingerprintBits =
    numPathBits + 2 * codeSize;  // note that this is only accurate if chirality
                                 // is not included

//! returns a numeric code for the atom (the atom's hash in the
//! atom-pair scheme)
/*!
  \param atom            the atom to be considered
  \param branchSubtract  (optional) a constant to subtract from
  the number of neighbors when the hash
  is calculated (used in the topological
  torsions code)
  \param includeChirality toggles the inclusions of bits indicating R/S
  chirality
*/
RDKIT_FINGERPRINTS_EXPORT std::uint32_t getAtomCode(
    const Atom *atom, unsigned int branchSubtract = 0,
    bool includeChirality = false);

//! returns an atom pair hash based on two atom hashes and the
//! distance between the atoms.
/*!
  \param codeI  the hash for the first atom
  \param codeJ  the hash for the second atom
  \param dist   the distance (number of bonds) between the two
  atoms
  \param includeChirality toggles the inclusions of bits indicating R/S
  chirality
*/
RDKIT_FINGERPRINTS_EXPORT std::uint32_t getAtomPairCode(
    std::uint32_t codeI, std::uint32_t codeJ, unsigned int dist,
    bool includeChirality = false);

//! returns an topological torsion hash based on the atom hashes
//! passed in
/*!
  \param atomCodes  the vector of atom hashes
*/
RDKIT_FINGERPRINTS_EXPORT std::uint64_t getTopologicalTorsionCode(
    const std::vector<std::uint32_t> &atomCodes, bool includeChirality = false);

RDKIT_FINGERPRINTS_EXPORT std::uint32_t getTopologicalTorsionHash(
    const std::vector<std::uint32_t> &pathCodes);

}  // namespace AtomPairs

namespace MorganFingerprints {

class RDKIT_FINGERPRINTS_EXPORT ss_matcher {
 public:
  ss_matcher();
  ss_matcher(const std::string &pattern);

  // const RDKit::ROMOL_SPTR &getMatcher() const { return m_matcher; }
  const RDKit::ROMol *getMatcher() const;

 private:
  RDKit::ROMOL_SPTR m_matcher;
};

typedef std::tuple<boost::dynamic_bitset<>, uint32_t, unsigned int> AccumTuple;

RDKIT_FINGERPRINTS_EXPORT extern std::vector<std::string> defaultFeatureSmarts;

//! returns the connectivity invariants for a molecule
/*!

  \param mol :    the molecule to be considered
  \param invars : used to return the results
  \param includeRingMembership : if set, whether or not the atom is in
             a ring will be used in the invariant list.
*/
RDKIT_FINGERPRINTS_EXPORT void getConnectivityInvariants(
    const ROMol &mol, std::vector<std::uint32_t> &invars,
    bool includeRingMembership = true);
const std::string morganConnectivityInvariantVersion = "1.0.0";

//! returns the feature invariants for a molecule
/*!

  \param mol:    the molecule to be considered
  \param invars : used to return the results
  \param patterns: if provided should contain the queries used to assign
  atom-types.
                   if not provided, feature definitions adapted from reference:
                   Gobbi and Poppinger, Biotech. Bioeng. _61_ 47-54 (1998)
                   will be used for Donor, Acceptor, Aromatic, Halogen, Basic,
  Acidic

*/
RDKIT_FINGERPRINTS_EXPORT void getFeatureInvariants(
    const ROMol &mol, std::vector<std::uint32_t> &invars,
    std::vector<const ROMol *> *patterns = nullptr);
const std::string morganFeatureInvariantVersion = "0.1.0";

}  // namespace MorganFingerprints

namespace RDKitFPUtils {

RDKIT_FINGERPRINTS_EXPORT void buildDefaultRDKitFingerprintAtomInvariants(
    const ROMol &mol, std::vector<std::uint32_t> &lAtomInvariants);

RDKIT_FINGERPRINTS_EXPORT void enumerateAllPaths(
    const ROMol &mol, std::map<int, std::list<std::vector<int>>> &allPaths,
    const std::vector<std::uint32_t> *fromAtoms, bool branchedPaths, bool useHs,
    unsigned int minPath, unsigned int maxPath);

RDKIT_FINGERPRINTS_EXPORT void identifyQueryBonds(
    const ROMol &mol, std::vector<const Bond *> &bondCache,
    std::vector<short> &isQueryBond);

RDKIT_FINGERPRINTS_EXPORT std::vector<unsigned int> generateBondHashes(
    const ROMol &mol, boost::dynamic_bitset<> &atomsInPath,
    const std::vector<const Bond *> &bondCache,
    const std::vector<short> &isQueryBond, const std::vector<int> &path,
    bool useBondOrder, const std::vector<std::uint32_t> *atomInvariants);

}  // namespace RDKitFPUtils

}  // namespace RDKit

#endif
