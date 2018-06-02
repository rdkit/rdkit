#ifndef RD_ATOMPAIRGEN_H_2018_06
#define RD_ATOMPAIRGEN_H_2018_06

#include <GraphMol/Fingerprints/FingerprintGenerator.h>

namespace RDKit {
namespace AtomPair {

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
const unsigned int numAtomPairFingerprintBits = numPathBits + 2 * codeSize;

class AtomPairArguments : public FingerprintArguments {
 public:
  const bool includeChirality;
  const bool use2D;
  const unsigned int minDistance;
  const unsigned int maxDistance;

  unsigned int getResultSize() const;

  AtomPairArguments(const bool countSimulation, const bool includeChirality,
                    const bool use2D, const unsigned int minDistance,
                    const unsigned int maxDistance);
};

class AtomPairAtomEnv : public AtomEnvironment<AtomPairArguments> {
  const std::vector<boost::uint32_t> *atomCodeCache;
  const unsigned int atomIdFirst;
  const unsigned int atomIdSecond;
  const unsigned int distance;

 public:
  boost::uint32_t getBitId(
      AtomPairArguments arguments,
      const std::vector<boost::uint32_t> *atomInvariants,
      const std::vector<boost::uint32_t> *bondInvariants) const;

  const std::vector<boost::uint32_t> *getAtomCodeCache() const;

  AtomPairAtomEnv(const std::vector<boost::uint32_t> *atomCodeCache,
                  const unsigned int atomIdFirst,
                  const unsigned int atomIdSecond, const unsigned int distance);
};

class AtomPairEnvGenerator
    : public AtomEnvironmentGenerator<AtomPairAtomEnv, AtomPairArguments> {
 public:
  std::vector<AtomPairAtomEnv> getEnvironments(
      const ROMol &mol, const AtomPairArguments arguments,
      const std::vector<boost::uint32_t> *fromAtoms = 0,
      const std::vector<boost::uint32_t> *ignoreAtoms = 0,
      const int confId = -1, const AdditionalOutput *additionalOutput = 0,
      const std::vector<boost::uint32_t> *atomInvariants = 0,
      const std::vector<boost::uint32_t> *bondInvariants = 0) const;

  void cleanUpEnvironments(
      const std::vector<AtomPairAtomEnv> &atomEnvironments) const;
};

}  // namespace AtomPair
}  // namespace RDKit

#endif
