#ifndef RD_ATOMPAIRGEN_H_2018_06
#define RD_ATOMPAIRGEN_H_2018_06

#include <GraphMol/Fingerprints/FingerprintGenerator.h>

namespace RDKit {
namespace AtomPair {

// constants taken from existing atom pairs implementation
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

/*
Class that holds atom pair fingerprinting arguments
includeChirality: if set, chirality will be used in the atom invariants (note:
this is ignored if atomInvariantsGenerator is present for the fingerprint
generator that uses this)
use2D: if set, the 2D (topological) distance matrix is used
minDistance: minimum distance between atoms to be considered in a pair. Default
is 1 bond.
maxDistance: maximum distance between atoms to be considered in a pair. Default
is maxPathLen-1 bonds
*/
class AtomPairArguments : public FingerprintArguments {
 public:
  const bool includeChirality;
  const bool use2D;
  const unsigned int minDistance;
  const unsigned int maxDistance;

  unsigned int getResultSize() const;

  AtomPairArguments(const bool countSimulation, const bool includeChirality,
                    const bool use2D, const unsigned int minDistance = 1,
                    const unsigned int maxDistance = (maxPathLen - 1));
};

/*
Class to hold atom environment data for atom pair fingerprinting
atomCodeCache: list of atom codes for all atoms in the molecule, used to avoid
generating atom code for the same atom multiple times. This is deleted by
calling AtomPairEnvGenerator::cleanUpEnvironments
atomIdFirst: index of the first atom of the pair
atomIdSecond: index of the second atom of the pair
distance: distance between the atoms
*/
class AtomPairAtomEnv : public AtomEnvironment {
  const std::vector<boost::uint32_t> *atomCodeCache;
  const unsigned int atomIdFirst;
  const unsigned int atomIdSecond;
  const unsigned int distance;

 public:
  boost::uint32_t getBitId(
      FingerprintArguments *arguments,
      const std::vector<boost::uint32_t> *atomInvariants,
      const std::vector<boost::uint32_t> *bondInvariants) const;

  /*
  returns a pointer to the atom code cache so it can be deleted
  */
  const std::vector<boost::uint32_t> *getAtomCodeCache() const;

  AtomPairAtomEnv(const std::vector<boost::uint32_t> *atomCodeCache,
                  const unsigned int atomIdFirst,
                  const unsigned int atomIdSecond, const unsigned int distance);
};

/*
Class to generate atom environments from a molecule to be used for atom pair
fingerprinting
*/
class AtomPairEnvGenerator : public AtomEnvironmentGenerator {
 public:
  std::vector<AtomEnvironment *> getEnvironments(
      const ROMol &mol, FingerprintArguments *arguments,
      const std::vector<boost::uint32_t> *fromAtoms = 0,
      const std::vector<boost::uint32_t> *ignoreAtoms = 0,
      const int confId = -1, const AdditionalOutput *additionalOutput = 0,
      const std::vector<boost::uint32_t> *atomInvariants = 0,
      const std::vector<boost::uint32_t> *bondInvariants = 0) const;

  void cleanUpEnvironments(
      std::vector<AtomEnvironment *> atomEnvironments) const;
};

FingerprintGenerator getAtomPairGenerator(
    const unsigned int minDistance = 1,
    const unsigned int maxDistance = maxPathLen - 1,
    const bool includeChirality = false, const bool use2D = true,
    const bool useCountSimulation = true,
    AtomInvariantsGenerator *atomInvariantsGenerator = 0,
    BondInvariantsGenerator *bondInvariantsGenerator = 0);
}  // namespace AtomPair
}  // namespace RDKit

#endif
