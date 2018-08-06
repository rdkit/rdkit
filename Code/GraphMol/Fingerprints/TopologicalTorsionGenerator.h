#ifndef RD_TOPOLOGICALTORSIONGEN_H_2018_07
#define RD_TOPOLOGICALTORSIONGEN_H_2018_07

#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/FingerprintUtil.h>

namespace RDKit {
namespace TopologicalTorsion {

template <typename OutputType>
class TopologicalTorsionArguments : public FingerprintArguments<OutputType> {
 public:
  const bool df_includeChirality;
  const uint32_t d_torsionAtomCount;

  OutputType getResultSize() const;

  std::string infoString() const;

  TopologicalTorsionArguments(const bool includeChirality,
                              const uint32_t torsionAtomCount,
                              const bool countSimulation,
                              const std::vector<std::uint32_t> countBounds,
                              const std::uint32_t foldedSize);
};

template <typename OutputType>
class TopologicalTorsionAtomEnv : public AtomEnvironment<OutputType> {
  const OutputType d_bitId;

 public:
  OutputType getBitId(FingerprintArguments<OutputType> *arguments,
                      const std::vector<std::uint32_t> *atomInvariants,
                      const std::vector<std::uint32_t> *bondInvariants,
                      const AdditionalOutput *additionalOutput,
                      const bool hashResults = false) const;

  TopologicalTorsionAtomEnv(OutputType bitId);
};

template <typename OutputType>
class TopologicalTorsionEnvGenerator
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

template <typename OutputType>
FingerprintGenerator<OutputType> *getTopologicalTorsionGenerator(
    const bool includeChirality = false, const uint32_t torsionAtomCount = 4,
    AtomInvariantsGenerator *atomInvariantsGenerator = nullptr,
    const bool countSimulation = true,
    const std::vector<std::uint32_t> countBounds = {1, 2, 4, 8},
    const std::uint32_t foldedSize = 2048, const bool ownsAtomInvGen = false);
}  // namespace TopologicalTorsion
}  // namespace RDKit

#endif
