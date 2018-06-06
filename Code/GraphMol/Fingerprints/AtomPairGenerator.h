#ifndef RD_ATOMPAIRGEN_H_2018_06
#define RD_ATOMPAIRGEN_H_2018_06

#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <cstdint>

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

unsigned int numPiElectrons(const Atom *atom);
std::uint32_t getAtomCode(const Atom *atom, unsigned int branchSubtract = 0,
                          bool includeChirality = false);
std::uint32_t getAtomPairCode(std::uint32_t codeI, std::uint32_t codeJ,
                              unsigned int dist, bool includeChirality = false);

/*!
  /brief class that holds atom-pair fingerprint specific arguments

 */
class AtomPairArguments : public FingerprintArguments {
 public:
  const bool df_includeChirality;
  const bool df_use2D;
  const unsigned int d_minDistance;
  const unsigned int d_maxDistance;

  std::uint64_t getResultSize() const;

  /*!
    /brief construct a new AtomPairArguments object

    /param countSimulation  if set, use count simulation while generating the
    fingerprint
    /param includeChirality if set, chirality will be used in the atom
    invariants, this is ignored if atomInvariantsGenerator is present for
    the /c FingerprintGenerator that uses this
    /param use2D            if set, the 2D (topological) distance matrix will be
    used
    /param minDistance      minimum distance between atoms to be considered in a
    pair, default is 1 bond
    /param maxDistance      maximum distance between atoms to be considered in a
    pair, default is maxPathLen-1 bonds
   */
  AtomPairArguments(const bool countSimulation, const bool includeChirality,
                    const bool use2D, const unsigned int minDistance = 1,
                    const unsigned int maxDistance = (maxPathLen - 1));
};

/*!
  /brief class that holds atom-environment data needed for atom-pair fingerprint
  generation

 */
class AtomPairAtomEnv : public AtomEnvironment {
  const unsigned int d_atomIdFirst;
  const unsigned int d_atomIdSecond;
  const unsigned int d_distance;
  const std::uint32_t d_atomCodeFirst;
  const std::uint32_t d_atomCodeSecond;

 public:
  std::uint32_t getBitId(
      FingerprintArguments *arguments,
      const std::vector<std::uint32_t> *atomInvariants,
      const std::vector<std::uint32_t> *bondInvariants) const;

  /*!
    /brief construct a new AtomPairAtomEnv object

    /param atomIdFirst      id of the first atom of the atom-pair
    /param atomIdSecond     id of the second atom of the atom-pair
    /param distance         distance between the atoms
    /param atomCodeFirst    hashed atom code for the first atom
    /param atomCodeSecond   hashed atom code for the second atom
   */
  AtomPairAtomEnv(const unsigned int atomIdFirst,
                  const unsigned int atomIdSecond, const unsigned int distance,
                  const std::uint32_t atomCodeFirst,
                  const std::uint32_t atomCodeSecond);
};

/*!
  /brief class that generates atom-environments for atom-pair fingerprint

 */
class AtomPairEnvGenerator : public AtomEnvironmentGenerator {
 public:
  std::vector<AtomEnvironment *> getEnvironments(
      const ROMol &mol, FingerprintArguments *arguments,
      const std::vector<std::uint32_t> *fromAtoms = nullptr,
      const std::vector<std::uint32_t> *ignoreAtoms = nullptr,
      const int confId = -1, const AdditionalOutput *additionalOutput = nullptr,
      const std::vector<std::uint32_t> *atomInvariants = nullptr,
      const std::vector<std::uint32_t> *bondInvariants = nullptr) const;

  void cleanUpEnvironments(
      std::vector<AtomEnvironment *> atomEnvironments) const;
};

/*!
  /brief helper function that generates a /c FingerprintGenerator that generates
  atom-pair fingerprints


  /param minDistance            minimum distance between atoms to be considered
  in a pair, default is 1 bond
  /param maxDistance            maximum distance between atoms to be considered
  in a pair, default is maxPathLen-1 bonds
  /param includeChirality       if set, chirality will be used in the atom
  invariants, this is ignored if atomInvariantsGenerator is provided
  /param use2D                  if set, the 2D (topological) distance matrix
  will be used
  /param useCountSimulation         if set, use count simulation while
  generating the fingerprint
  /param atomInvariantsGenerator    atom invariants to be used during
  fingerprint generation
  /param bondInvariantsGenerator    bond invariants to be used during
  fingerprint generation

  /return FingerprintGenerator that generates atom-pair fingerprints
 */
FingerprintGenerator *getAtomPairGenerator(
    const unsigned int minDistance = 1,
    const unsigned int maxDistance = maxPathLen - 1,
    const bool includeChirality = false, const bool use2D = true,
    const bool useCountSimulation = true,
    AtomInvariantsGenerator *atomInvariantsGenerator = nullptr,
    BondInvariantsGenerator *bondInvariantsGenerator = nullptr);
}  // namespace AtomPair
}  // namespace RDKit

#endif
