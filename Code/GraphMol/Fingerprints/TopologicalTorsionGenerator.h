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
#ifndef RD_TOPOLOGICALTORSIONGEN_H_2018_07
#define RD_TOPOLOGICALTORSIONGEN_H_2018_07

#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/FingerprintUtil.h>

namespace RDKit {
namespace TopologicalTorsion {

template <typename OutputType>
class RDKIT_FINGERPRINTS_EXPORT TopologicalTorsionArguments
    : public FingerprintArguments<OutputType> {
 public:
  const bool df_includeChirality;
  const uint32_t d_torsionAtomCount;

  OutputType getResultSize() const;

  std::string infoString() const;

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
  TopologicalTorsionArguments(const bool includeChirality,
                              const uint32_t torsionAtomCount,
                              const bool countSimulation,
                              const std::vector<std::uint32_t> countBounds,
                              const std::uint32_t fpSize);
};

template <typename OutputType>
class RDKIT_FINGERPRINTS_EXPORT TopologicalTorsionAtomEnv
    : public AtomEnvironment<OutputType> {
  const OutputType d_bitId;

 public:
  OutputType getBitId(FingerprintArguments<OutputType> *arguments,
                      const std::vector<std::uint32_t> *atomInvariants,
                      const std::vector<std::uint32_t> *bondInvariants,
                      const AdditionalOutput *additionalOutput,
                      const bool hashResults = false) const;
  /**
   \brief Construct a new Topological Torsion Atom Env object

   \param bitId bitId generated for this environment
   */
  TopologicalTorsionAtomEnv(OutputType bitId);
};

template <typename OutputType>
class RDKIT_FINGERPRINTS_EXPORT TopologicalTorsionEnvGenerator
    : public AtomEnvironmentGenerator<OutputType> {
 public:
  std::vector<AtomEnvironment<OutputType> *> getEnvironments(
      const ROMol &mol, FingerprintArguments<OutputType> *arguments,
      const std::vector<std::uint32_t> *fromAtoms,
      const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
      const AdditionalOutput *additionalOutput,
      const std::vector<std::uint32_t> *atomInvariants,
      const std::vector<std::uint32_t> *bondInvariants,
      const bool hashResults = false) const;

  std::string infoString() const;
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
 */
template <typename OutputType>
RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<OutputType> *
getTopologicalTorsionGenerator(
    const bool includeChirality = false, const uint32_t torsionAtomCount = 4,
    AtomInvariantsGenerator *atomInvariantsGenerator = nullptr,
    const bool countSimulation = true,
    const std::vector<std::uint32_t> countBounds = {1, 2, 4, 8},
    const std::uint32_t fpSize = 2048, const bool ownsAtomInvGen = false);
}  // namespace TopologicalTorsion
}  // namespace RDKit

#endif
