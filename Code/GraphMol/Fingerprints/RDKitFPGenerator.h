#ifndef RD_RDFINGERPRINTGEN_H_2018_07
#define RD_RDFINGERPRINTGEN_H_2018_07

#include <GraphMol/Fingerprints/FingerprintGenerator.h>

namespace RDKit {
namespace RDKitFP {

template <typename OutputType>
class RDKitFPArguments : public FingerprintArguments<OutputType> {
 public:
  const unsigned int d_minPath;
  const unsigned int d_maxPath;
  const bool df_useHs;
  const bool df_branchedPaths;
  const bool df_useBondOrder;

  OutputType getResultSize() const;

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
  RDKitFPAtomInvGenerator *clone() const;
};

template <typename OutputType>
class RDKitFPAtomEnv : public AtomEnvironment<OutputType> {
  const OutputType d_bitId;
  const boost::dynamic_bitset<> d_atomsInPath;

 public:
  OutputType getBitId(FingerprintArguments<OutputType> *arguments,
                      const std::vector<std::uint32_t> *atomInvariants,
                      const std::vector<std::uint32_t> *bondInvariants,
                      const AdditionalOutput *additionalOutput) const;

  RDKitFPAtomEnv(const OutputType bitId,
                 const boost::dynamic_bitset<> atomsInPath);
};

template <typename OutputType>
class RDKitFPEnvGenerator : public AtomEnvironmentGenerator<OutputType> {
 public:
  std::vector<AtomEnvironment<OutputType> *> getEnvironments(
      const ROMol &mol, FingerprintArguments<OutputType> *arguments,
      const std::vector<std::uint32_t> *fromAtoms,
      const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
      const AdditionalOutput *additionalOutput,
      const std::vector<std::uint32_t> *atomInvariants,
      const std::vector<std::uint32_t> *bondInvariants) const;

  std::string infoString() const;
};

template <typename OutputType>
FingerprintGenerator<OutputType> *getRDKitFPGenerator(
    const unsigned int minPath = 1, const unsigned int maxPath = 7,
    const bool useHs = true, const bool branchedPaths = true,
    const bool useBondOrder = true,
    AtomInvariantsGenerator *atomInvariantsGenerator = nullptr,
    const bool countSimulation = true,
    const std::vector<std::uint32_t> countBounds = {1, 2, 4, 8},
    const std::uint32_t foldedSize = 2048, const bool ownsAtomInvGen = false);

}  // namespace RDKitFP
}  // namespace RDKit

#endif