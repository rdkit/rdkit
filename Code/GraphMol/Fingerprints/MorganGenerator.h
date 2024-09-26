//
//  Copyright (C) 2018-2022 Boran Adas and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef RD_MORGANGEN_H_2018_07
#define RD_MORGANGEN_H_2018_07

#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <cstdint>

namespace RDKit {

namespace MorganFingerprint {

/**
 \brief Default atom invariants generator for Morgan fingerprint, generates
 ECFP-type invariants

 */
class RDKIT_FINGERPRINTS_EXPORT MorganAtomInvGenerator
    : public AtomInvariantsGenerator {
  const bool df_includeRingMembership;

 public:
  /**
   \brief Construct a new MorganAtomInvGenerator object

   \param includeRingMembership : if set, whether or not the atom is in a ring
   will be used in the invariant list.
   */
  MorganAtomInvGenerator(const bool includeRingMembership = true);

  std::vector<std::uint32_t> *getAtomInvariants(
      const ROMol &mol) const override;

  std::string infoString() const override;
  MorganAtomInvGenerator *clone() const override;
};

/**
 \brief Alternative atom invariants generator for Morgan fingerprint, generate
 FCFP-type invariants

 */
class RDKIT_FINGERPRINTS_EXPORT MorganFeatureAtomInvGenerator
    : public AtomInvariantsGenerator {
  std::vector<const ROMol *> *dp_patterns;

 public:
  /**
   \brief Construct a new MorganFeatureAtomInvGenerator object

   \param patterns : if provided should contain the queries used to assign
   atom-types. if not provided, feature definitions adapted from reference:
   Gobbi and Poppinger, Biotech. Bioeng. _61_ 47-54 (1998) will be used for
   Donor, Acceptor, Aromatic, Halogen, Basic, Acidic.
   */
  MorganFeatureAtomInvGenerator(std::vector<const ROMol *> *patterns = nullptr);

  std::vector<std::uint32_t> *getAtomInvariants(
      const ROMol &mol) const override;

  std::string infoString() const override;
  MorganFeatureAtomInvGenerator *clone() const override;
};

/**
 \brief Bond invariants generator for Morgan fingerprint

 */
class RDKIT_FINGERPRINTS_EXPORT MorganBondInvGenerator
    : public BondInvariantsGenerator {
  const bool df_useBondTypes;
  const bool df_useChirality;

 public:
  /**
   \brief Construct a new MorganBondInvGenerator object

   \param useBondTypes : if set, bond types will be included as a part of the
   bond invariants
   \param useChirality : if set, chirality information will be included as a
   part of the bond invariants
   */
  MorganBondInvGenerator(const bool useBondTypes = true,
                         const bool useChirality = false);

  std::vector<std::uint32_t> *getBondInvariants(
      const ROMol &mol) const override;

  std::string infoString() const override;
  MorganBondInvGenerator *clone() const override;
  ~MorganBondInvGenerator() override = default;
};

/**
 \brief Class for holding Morgan fingerprint specific arguments

 */
class RDKIT_FINGERPRINTS_EXPORT MorganArguments : public FingerprintArguments {
 public:
  bool df_onlyNonzeroInvariants = false;
  unsigned int d_radius = 3;
  bool df_includeRedundantEnvironments = false;

  std::string infoString() const override;

  /**
   \brief Construct a new MorganArguments object

   \param radius the number of iterations to grow the fingerprint
   \param countSimulation if set, use count simulation while generating the
    fingerprint
   \param includeChirality if set, chirality information will be added to the
   generated bit id, independently from bond invariants
   \param onlyNonzeroInvariants if set, bits will only be set from atoms that
   have a nonzero invariant
   \param countBounds boundaries for count simulation, corresponding bit will
   be set if the count is higher than the number provided for that spot
   \param fpSize size of the generated fingerprint, does not affect the sparse
   versions
   \param includeRedundantEnvironments if set redundant environments will be
   included in the fingerprint
  */
  MorganArguments(unsigned int radius, bool countSimulation = false,
                  bool includeChirality = false,
                  bool onlyNonzeroInvariants = false,
                  std::vector<std::uint32_t> countBounds = {1, 2, 4, 8},
                  std::uint32_t fpSize = 2048,
                  bool includeRedundantEnvironments = false)
      : FingerprintArguments(countSimulation, countBounds, fpSize, 1,
                             includeChirality),
        df_onlyNonzeroInvariants(onlyNonzeroInvariants),
        d_radius(radius),
        df_includeRedundantEnvironments(includeRedundantEnvironments) {};
};

/**
 \brief Class for holding the bit-id created from Morgan fingerprint
 environments and the additional data necessary extra outputs

 */
template <typename OutputType>
class RDKIT_FINGERPRINTS_EXPORT MorganAtomEnv
    : public AtomEnvironment<OutputType> {
  const OutputType d_code;
  const unsigned int d_atomId;
  const unsigned int d_layer;

 public:
  OutputType getBitId(
      FingerprintArguments *arguments,                   // unused
      const std::vector<std::uint32_t> *atomInvariants,  // unused
      const std::vector<std::uint32_t> *bondInvariants,  // unused
      AdditionalOutput *additionalOutput,                // unused
      const bool hashResults = false,                    // unused
      const std::uint64_t fpSize = 0                     // unused
  ) const override;
  void updateAdditionalOutput(AdditionalOutput *output,
                              size_t bitId) const override;

  /**
   \brief Construct a new MorganAtomEnv object

   \param code bit id generated from this environment
   \param atomId atom id of the atom at the center of this environment
   \param layer radius of this environment
   */
  MorganAtomEnv(const std::uint32_t code, const unsigned int atomId,
                const unsigned int layer);
};

/**
 \brief Class that generates atom environments for Morgan fingerprint

 */
template <typename OutputType>
class RDKIT_FINGERPRINTS_EXPORT MorganEnvGenerator
    : public AtomEnvironmentGenerator<OutputType> {
 public:
  std::vector<AtomEnvironment<OutputType> *> getEnvironments(
      const ROMol &mol, FingerprintArguments *arguments,
      const std::vector<std::uint32_t> *fromAtoms,
      const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
      const AdditionalOutput *additionalOutput,
      const std::vector<std::uint32_t> *atomInvariants,
      const std::vector<std::uint32_t> *bondInvariants,
      const bool hashResults = false) const override;

  std::string infoString() const override;
  OutputType getResultSize() const override;
};

/**
 \brief Get a fingerprint generator for Morgan fingerprint

 \tparam OutputType determines the size of the bitIds and the result, can be 32
 or 64 bit unsigned integer

 \param radius the number of iterations to grow the fingerprint

 \param countSimulation if set, use count simulation while generating the
 fingerprint

 \param includeChirality if set, chirality information will be added to the
 generated bit id, independently from bond invariants

 \param onlyNonzeroInvariants if set, bits will only be set from atoms that
 have a nonzero invariant

 \param countBounds boundaries for count simulation, corresponding bit will be
 set if the count is higher than the number provided for that spot

 \param fpSize size of the generated fingerprint, does not affect the sparse
 versions
 \param countSimulation if set, use count simulation while generating the
 fingerprint
 \param includeChirality sets includeChirality flag for both MorganArguments
 and the default bond generator MorganBondInvGenerator
 \param useBondTypes if set, bond types will be included as a part of the
 default bond invariants
 \param onlyNonzeroInvariants if set, bits will only be set from atoms that
 have a nonzero invariant
 \param includeRedundantEnvironments if set redundant environments will be
 included in the fingerprint
 \param atomInvariantsGenerator custom atom invariants generator to use
 \param bondInvariantsGenerator custom  bond invariants generator to use
 \param ownsAtomInvGen  if set atom invariants  generator is destroyed with the
 fingerprint generator
 \param ownsBondInvGen  if set bond invariants generator is destroyed with the
 fingerprint generator

 \return FingerprintGenerator<OutputType>* that generates Morgan fingerprints

This generator supports the following \c AdditionalOutput types:
  - \c atomToBits : which bits each atom is the central atom for
  - \c atomCounts : how many bits each atom sets
  - \c bitInfoMap : map from bitId to (atomId, radius) pairs

 */
template <typename OutputType>
RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<OutputType> *getMorganGenerator(
    unsigned int radius, bool countSimulation, bool includeChirality,
    bool useBondTypes, bool onlyNonzeroInvariants,
    bool includeRedundantEnvironments,
    AtomInvariantsGenerator *atomInvariantsGenerator = nullptr,
    BondInvariantsGenerator *bondInvariantsGenerator = nullptr,
    std::uint32_t fpSize = 2048,
    std::vector<std::uint32_t> countBounds = {1, 2, 4, 8},
    bool ownsAtomInvGen = false, bool ownsBondInvGen = false);

/**
 \brief Get a fingerprint generator for Morgan fingerprint

 \tparam OutputType determines the size of the bitIds and the result, can be 32
 or 64 bit unsigned integer

 \param radius the number of iterations to grow the fingerprint

 \param countSimulation if set, use count simulation while generating the
 fingerprint

 \param includeChirality if set, chirality information will be added to the
 generated bit id, independently from bond invariants

 \param onlyNonzeroInvariants if set, bits will only be set from atoms that
 have a nonzero invariant

 \param countBounds boundaries for count simulation, corresponding bit will be
 set if the count is higher than the number provided for that spot

 \param fpSize size of the generated fingerprint, does not affect the sparse
 versions
 \param countSimulation if set, use count simulation while generating the
 fingerprint
 \param includeChirality sets includeChirality flag for both MorganArguments
 and the default bond generator MorganBondInvGenerator
 \param useBondTypes if set, bond types will be included as a part of the
 default bond invariants
 \param onlyNonzeroInvariants if set, bits will only be set from atoms that
 have a nonzero invariant
 \param atomInvariantsGenerator custom atom invariants generator to use
 \param bondInvariantsGenerator custom  bond invariants generator to use
 \param ownsAtomInvGen  if set atom invariants  generator is destroyed with the
 fingerprint generator
 \param ownsBondInvGen  if set bond invariants generator is destroyed with the
 fingerprint generator

 \return FingerprintGenerator<OutputType>* that generates Morgan fingerprints

This generator supports the following \c AdditionalOutput types:
  - \c atomToBits : which bits each atom is the central atom for
  - \c atomCounts : how many bits each atom sets
  - \c bitInfoMap : map from bitId to (atomId, radius) pairs

 */
template <typename OutputType>
FingerprintGenerator<OutputType> *getMorganGenerator(
    unsigned int radius, bool countSimulation = false,
    bool includeChirality = false, bool useBondTypes = true,
    bool onlyNonzeroInvariants = false,
    AtomInvariantsGenerator *atomInvariantsGenerator = nullptr,
    BondInvariantsGenerator *bondInvariantsGenerator = nullptr,
    std::uint32_t fpSize = 2048,
    std::vector<std::uint32_t> countBounds = {1, 2, 4, 8},
    bool ownsAtomInvGen = false, bool ownsBondInvGen = false) {
  return getMorganGenerator<OutputType>(
      radius, countSimulation, includeChirality, useBondTypes,
      onlyNonzeroInvariants, false, atomInvariantsGenerator,
      bondInvariantsGenerator, fpSize, countBounds, ownsAtomInvGen,
      ownsBondInvGen);
};

}  // namespace MorganFingerprint
}  // namespace RDKit

#endif
