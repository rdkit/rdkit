#ifndef RD_MORGANGEN_H_2018_07
#define RD_MORGANGEN_H_2018_07

#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <cstdint>

namespace RDKit {

namespace MorganFingerprint {

class MorganAtomInvGenerator : public AtomInvariantsGenerator {
  const bool df_includeRingMembership;

 public:
  MorganAtomInvGenerator(const bool includeRingMembership = true);

  std::vector<std::uint32_t> *getAtomInvariants(const ROMol &mol) const;

  std::string infoString() const;
};

class MorganArguments : public FingerprintArguments {
 public:
  const bool df_includeChirality;
  const bool df_useBondTypes;
  const bool df_onlyNonzeroInvariants;
  const unsigned int d_radius;

  std::uint64_t getResultSize() const;

  std::string infoString() const;

  MorganArguments(const unsigned int radius, const bool countSimulation = true,
                  const bool includeChirality = false,
                  const bool useBondTypes = true,
                  const bool onlyNonzeroInvariants = false,
                  const std::vector<std::uint32_t> countBounds = {1, 2, 4, 8},
                  const std::uint32_t foldedSize = 2048);
};

class MorganAtomEnv : public AtomEnvironment {
  const std::uint32_t d_code;
  const unsigned int d_atomId;
  const unsigned int d_layer;

 public:
  std::uint32_t getBitId(FingerprintArguments *arguments,
                         const std::vector<std::uint32_t> *atomInvariants,
                         const std::vector<std::uint32_t> *bondInvariants,
                         const AdditionalOutput *additionalOutput) const;

  MorganAtomEnv(const std::uint32_t code, const unsigned int atomId,
                const unsigned int layer);
};

class MorganEnvGenerator : public AtomEnvironmentGenerator {
 public:
  std::vector<AtomEnvironment *> getEnvironments(
      const ROMol &mol, FingerprintArguments *arguments,
      const std::vector<std::uint32_t> *fromAtoms,
      const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
      const AdditionalOutput *additionalOutput,
      const std::vector<std::uint32_t> *atomInvariants,
      const std::vector<std::uint32_t> *bondInvariants) const;

  std::string infoString() const;
};

FingerprintGenerator *getMorganGenerator(
    const unsigned int radius, const bool countSimulation = true,
    const bool includeChirality = false, const bool useBondTypes = true,
    const bool onlyNonzeroInvariants = false,
    const std::vector<std::uint32_t> countBounds = {1, 2, 4, 8},
    const std::uint32_t foldedSize = 2048,
    AtomInvariantsGenerator *atomInvariantsGenerator = nullptr,
    BondInvariantsGenerator *bondInvariantsGenerator = nullptr);

}  // namespace MorganFingerprint
}  // namespace RDKit

#endif