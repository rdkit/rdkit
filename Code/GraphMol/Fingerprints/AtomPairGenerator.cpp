#include <GraphMol/RDKitBase.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/AtomPairGenerator.h>
#include <cstdint>

namespace RDKit {
namespace AtomPair {

//! taken from the existing implementation
unsigned int numPiElectrons(const Atom *atom) {
  PRECONDITION(atom, "no atom");
  unsigned int res = 0;
  if (atom->getIsAromatic()) {
    res = 1;
  } else if (atom->getHybridization() != Atom::SP3) {
    unsigned int val = static_cast<unsigned int>(atom->getExplicitValence());
    val -= atom->getNumExplicitHs();
    CHECK_INVARIANT(val >= atom->getDegree(),
                    "explicit valence exceeds atom degree");
    res = val - atom->getDegree();
  }
  return res;
}

//! taken from the existing implementation
std::uint32_t getAtomCode(const Atom *atom, unsigned int branchSubtract,
                          bool includeChirality) {
  PRECONDITION(atom, "no atom");

  std::uint32_t numBranches = 0;
  if (atom->getDegree() > branchSubtract) {
    numBranches = atom->getDegree() - branchSubtract;
  }

  std::uint32_t branchBits = numBranches % maxNumBranches;
  std::uint32_t nPi = numPiElectrons(atom) % maxNumPi;

  std::uint32_t typeIdx = 0;
  std::uint32_t nTypes = 1 << numTypeBits;
  while (typeIdx < nTypes) {
    if (atomNumberTypes[typeIdx] ==
        static_cast<unsigned int>(atom->getAtomicNum())) {
      break;
    } else if (atomNumberTypes[typeIdx] >
               static_cast<unsigned int>(atom->getAtomicNum())) {
      typeIdx = nTypes;
      break;
    }
    ++typeIdx;
  }
  if (typeIdx == nTypes) --typeIdx;

  std::uint32_t cipCodeBits = 0;
  if (includeChirality) {
    std::string cipCode;
    if (atom->getPropIfPresent(common_properties::_CIPCode, cipCode)) {
      if (cipCode == "R") {
        cipCodeBits = 1;
      } else if (cipCode == "S") {
        cipCodeBits = 2;
      }
    }
  }

  std::uint32_t code;
  code = branchBits;
  code |= nPi << numBranchBits;
  code |= typeIdx << (numBranchBits + numPiBits);
  code |= cipCodeBits << (numBranchBits + numPiBits + numTypeBits);

  POSTCONDITION(
      code < static_cast<std::uint32_t>(
                 1 << (codeSize + (includeChirality ? numChiralBits : 0))),
      "code exceeds number of bits");
  return code;
};

std::uint32_t getAtomPairCode(std::uint32_t codeI, std::uint32_t codeJ,
                              unsigned int dist, bool includeChirality) {
  PRECONDITION(dist < maxPathLen, "dist too long");
  std::uint32_t res = dist;
  res |= std::min(codeI, codeJ) << numPathBits;
  res |= std::max(codeI, codeJ)
         << (numPathBits + codeSize + (includeChirality ? numChiralBits : 0));
  return res;
}

std::uint64_t AtomPairArguments::getResultSize() const {
  return (1 << (numAtomPairFingerprintBits +
                2 * (df_includeChirality ? numChiralBits : 0)));
}

AtomPairArguments::AtomPairArguments(const bool countSimulation,
                                     const bool includeChirality,
                                     const bool use2D,
                                     const unsigned int minDistance,
                                     const unsigned int maxDistance)
    : FingerprintArguments(countSimulation),
      df_includeChirality(includeChirality),
      df_use2D(use2D),
      d_minDistance(minDistance),
      d_maxDistance(maxDistance) {
  PRECONDITION(minDistance <= maxDistance, "bad distances provided");
}

std::uint32_t AtomPairAtomEnv::getBitId(
    FingerprintArguments *arguments,
    const std::vector<std::uint32_t> *atomInvariants,
    const std::vector<std::uint32_t> *bondInvariants,
    const AdditionalOutput *additionalOutput) const {
  AtomPairArguments *atomPairArguments =
      dynamic_cast<AtomPairArguments *>(arguments);

  std::uint32_t bitId =
      getAtomPairCode(d_atomCodeFirst, d_atomCodeSecond, d_distance,
                      atomPairArguments->df_includeChirality);

  if (additionalOutput && additionalOutput->atomToBits) {
    additionalOutput->atomToBits->at(d_atomIdFirst).push_back(bitId);
    additionalOutput->atomToBits->at(d_atomIdSecond).push_back(bitId);
  }
  return bitId;
}

AtomPairAtomEnv::AtomPairAtomEnv(const unsigned int atomIdFirst,
                                 const unsigned int atomIdSecond,
                                 const unsigned int distance,
                                 const std::uint32_t atomCodeFirst,
                                 const std::uint32_t atomCodeSecond)
    : d_atomIdFirst(atomIdFirst),
      d_atomIdSecond(atomIdSecond),
      d_distance(distance),
      d_atomCodeFirst(atomCodeFirst),
      d_atomCodeSecond(atomCodeSecond) {}

std::vector<AtomEnvironment *> AtomPairEnvGenerator::getEnvironments(
    const ROMol &mol, FingerprintArguments *arguments,
    const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput,
    const std::vector<std::uint32_t> *atomInvariants,
    const std::vector<std::uint32_t> *bondInvariants) const {
  PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
               "bad atomInvariants size");
  const unsigned int atomCount = mol.getNumAtoms();
  PRECONDITION(!additionalOutput || !additionalOutput->atomToBits ||
                   additionalOutput->atomToBits->size() == atomCount,
               "bad atomToBits size in AdditionalOutput");

  AtomPairArguments *atomPairArguments =
      dynamic_cast<AtomPairArguments *>(arguments);
  std::vector<AtomEnvironment *> result = std::vector<AtomEnvironment *>();
  const double *distanceMatrix;
  if (atomPairArguments->df_use2D) {
    distanceMatrix = MolOps::getDistanceMat(mol);
  } else {
    distanceMatrix = MolOps::get3DDistanceMat(mol, confId);
  }

  std::vector<std::uint32_t> atomCodeCache = std::vector<std::uint32_t>();
  for (ROMol::ConstAtomIterator atomItI = mol.beginAtoms();
       atomItI != mol.endAtoms(); ++atomItI) {
    if (!atomInvariants) {
      atomCodeCache.push_back(
          getAtomCode(*atomItI, 0, atomPairArguments->df_includeChirality));
    } else {
      atomCodeCache.push_back((*atomInvariants)[(*atomItI)->getIdx()] %
                              ((1 << codeSize) - 1));
    }
  }

  for (ROMol::ConstAtomIterator atomItI = mol.beginAtoms();
       atomItI != mol.endAtoms(); ++atomItI) {
    unsigned int i = (*atomItI)->getIdx();
    if (ignoreAtoms && std::find(ignoreAtoms->begin(), ignoreAtoms->end(), i) !=
                           ignoreAtoms->end()) {
      continue;
    }

    for (ROMol::ConstAtomIterator atomItJ = atomItI + 1;
         atomItJ != mol.endAtoms(); ++atomItJ) {
      unsigned int j = (*atomItJ)->getIdx();
      if (ignoreAtoms && std::find(ignoreAtoms->begin(), ignoreAtoms->end(),
                                   j) != ignoreAtoms->end()) {
        continue;
      }

      if (fromAtoms &&
          (std::find(fromAtoms->begin(), fromAtoms->end(), i) ==
           fromAtoms->end()) &&
          (std::find(fromAtoms->begin(), fromAtoms->end(), j) ==
           fromAtoms->end())) {
        continue;
      }
      unsigned int distance =
          static_cast<unsigned int>(floor(distanceMatrix[i * atomCount + j]));

      if (distance >= atomPairArguments->d_minDistance &&
          distance <= atomPairArguments->d_maxDistance) {
        result.push_back(new AtomPairAtomEnv(i, j, distance, atomCodeCache[i],
                                             atomCodeCache[j]));
      }
    }
  }

  return result;
}

FingerprintGenerator *getAtomPairGenerator(
    const unsigned int minDistance, const unsigned int maxDistance,
    const bool includeChirality, const bool use2D,
    const bool useCountSimulation,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    BondInvariantsGenerator *bondInvariantsGenerator) {
  AtomEnvironmentGenerator *atomPairEnvGenerator =
      new AtomPair::AtomPairEnvGenerator();
  FingerprintArguments *atomPairArguments = new AtomPair::AtomPairArguments(
      useCountSimulation, includeChirality, use2D, minDistance, maxDistance);

  return new FingerprintGenerator(atomPairEnvGenerator, atomPairArguments,
                                  atomInvariantsGenerator,
                                  bondInvariantsGenerator);
}

}  // namespace AtomPair
}  // namespace RDKit
