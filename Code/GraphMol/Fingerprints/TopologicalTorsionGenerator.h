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
#ifndef RD_TOPOLOGICALTORSIONGEN_H_2018_07
#define RD_TOPOLOGICALTORSIONGEN_H_2018_07

#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/FingerprintUtil.h>

namespace RDKit {
namespace TopologicalTorsion {

class RDKIT_FINGERPRINTS_EXPORT TopologicalTorsionArguments
    : public FingerprintArguments {
 public:
  uint32_t d_torsionAtomCount = 4;
  bool df_onlyShortestPaths = false;

  std::string infoString() const override;

  /**
   \brief Construct a new Topological Torsion Arguments object

   \param includeChirality if set, chirality will be used in sparse result
   \param torsionAtomCount the number of atoms to include in the "torsions"
   \param useCountSimulation         if set, use count simulation while
   generating the fingerprint
   \param countBounds  boundaries for count simulation, corresponding bit will
   be set if the count is higher than the number provided for that spot
   \param fpSize size of the generated fingerprint, does not affect the sparse
   versions
   */
  TopologicalTorsionArguments(
      const bool includeChirality = false, const uint32_t torsionAtomCount = 4,
      const bool countSimulation = true,
      const std::vector<std::uint32_t> countBounds = {1, 2, 4, 8},
      const std::uint32_t fpSize = 2048);
};

template <typename OutputType>
class RDKIT_FINGERPRINTS_EXPORT TopologicalTorsionAtomEnv
    : public AtomEnvironment<OutputType> {
  const OutputType d_bitId;
  const INT_VECT d_atomPath;

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
   \brief Construct a new Topological Torsion Atom Env object

   \param bitId bitId generated for this environment
   */
  TopologicalTorsionAtomEnv(OutputType bitId, INT_VECT atomPath)
      : d_bitId(bitId), d_atomPath(std::move(atomPath)) {}
};

template <typename OutputType>
class RDKIT_FINGERPRINTS_EXPORT TopologicalTorsionEnvGenerator
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
 \brief Get the Topological Torsion Generator object

 \tparam OutputType determines the size of the bitIds and the result, can only
 be 64 bit unsigned integer for this type
 \param includeChirality includeChirality argument for both the default atom
 invariants generator and the fingerprint arguments
 \param torsionAtomCount the number of atoms to include in the "torsions"
 \param atomInvariantsGenerator custom atom invariants generator to use
 \param useCountSimulation         if set, use count simulation while
 generating the fingerprint
 \param countBounds  boundaries for count simulation, corresponding bit will
 be set if the count is higher than the number provided for that spot
 \param fpSize size of the generated fingerprint, does not affect the sparse
 versions
 \param ownsAtomInvGen  if set atom invariants generator is destroyed with the
 fingerprint generator

 /return FingerprintGenerator<OutputType>* that generates topological-torsion
 fingerprints

 This generator supports the following \c AdditionalOutput types:
  - \c atomToBits : which bits each atom is involved in
  - \c atomCounts : how many bits each atom sets
  - \c bitPaths : map from bitId to vectors of atom indices

 */
template <typename OutputType>
RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<OutputType> *
getTopologicalTorsionGenerator(
    bool includeChirality = false, uint32_t torsionAtomCount = 4,
    AtomInvariantsGenerator *atomInvariantsGenerator = nullptr,
    bool countSimulation = true, std::uint32_t fpSize = 2048,
    std::vector<std::uint32_t> countBounds = {1, 2, 4, 8},
    bool ownsAtomInvGen = false);
//! overload
template <typename OutputType>
FingerprintGenerator<OutputType> *getTopologicalTorsionGenerator(
    const TopologicalTorsionArguments &args,
    AtomInvariantsGenerator *atomInvariantsGenerator = nullptr,
    const bool ownsAtomInvGen = false);
}  // namespace TopologicalTorsion
}  // namespace RDKit

#endif
