#ifndef RD_RDFINGERPRINTGEN_H_2018_07
#define RD_RDFINGERPRINTGEN_H_2018_07

#include <GraphMol/Fingerprints/FingerprintGenerator.h>

namespace RDKit {
namespace RDKitFP {

class RDKitFPArguments : public FingerprintArguments {
 public:
  const unsigned int d_minPath;
  const unsigned int d_maxPath;
  const bool df_useHs;
  const bool df_branchedPaths;
  const bool df_useBondOrder;

  std::uint64_t getResultSize() const;

  std::string infoString() const;

  RDKitFPArguments(unsigned int minPath, unsigned int maxPath, bool useHs,
                   bool branchedPaths, bool useBondOrder,
                   const bool countSimulation,
                   const std::vector<std::uint32_t> countBounds,
                   const std::uint32_t foldedSize);
};

class RDKitFPAtomInvGenerator : public AtomInvariantsGenerator {
 public:
  std::vector<std::uint32_t> *getAtomInvariants(const ROMol &mol) const;

  std::string infoString() const;
};

class RDKitFPAtomEnv : public AtomEnvironment {
  const std::uint32_t d_bitId;
  const boost::dynamic_bitset<> d_atomsInPath;

 public:
  std::uint32_t getBitId(FingerprintArguments *arguments,
                         const std::vector<std::uint32_t> *atomInvariants,
                         const std::vector<std::uint32_t> *bondInvariants,
                         const AdditionalOutput *additionalOutput) const;

  RDKitFPAtomEnv(const std::uint32_t bitId,
                 const boost::dynamic_bitset<> atomsInPath);
};

class RDKitFPEnvGenerator : public AtomEnvironmentGenerator {
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

FingerprintGenerator *getRDKitFPGenerator(
    const unsigned int minPath = 1, const unsigned int maxPath = 7,
    const bool useHs = true, const bool branchedPaths = true,
    const bool useBondOrder = true,
    AtomInvariantsGenerator *atomInvariantsGenerator = nullptr,
    BondInvariantsGenerator *bondInvariantsGenerator = nullptr,
    const bool countSimulation = true,
    const std::vector<std::uint32_t> countBounds = {1, 2, 4, 8},
    const std::uint32_t foldedSize = 2048);

}  // namespace RDKitFP
}  // namespace RDKit

#endif