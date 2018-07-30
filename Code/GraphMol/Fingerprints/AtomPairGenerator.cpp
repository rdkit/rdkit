#include <GraphMol/RDKitBase.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/AtomPairGenerator.h>
#include <GraphMol/Fingerprints/FingerprintUtil.h>

namespace RDKit {
namespace AtomPair {
using namespace AtomPairs;

AtomPairAtomInvGenerator::AtomPairAtomInvGenerator(bool includeChirality)
    : df_includeChirality(includeChirality) {}

std::vector<std::uint32_t> *AtomPairAtomInvGenerator::getAtomInvariants(
    const ROMol &mol) const {
  std::vector<std::uint32_t> *atomInvariants =
      new std::vector<std::uint32_t>(mol.getNumAtoms());

  for (ROMol::ConstAtomIterator atomItI = mol.beginAtoms();
       atomItI != mol.endAtoms(); ++atomItI) {
    (*atomInvariants)[(*atomItI)->getIdx()] =
        getAtomCode(*atomItI, 0, df_includeChirality);
  }

  return atomInvariants;
}

std::string AtomPairAtomInvGenerator::infoString() const {
  return "AtomPairInvariantGenerator includeChirality=" +
         std::to_string(df_includeChirality);
}

template <>
std::uint64_t AtomPairArguments<std::uint64_t>::getResultSize() const {
  return (1 << (numAtomPairFingerprintBits +
                2 * (df_includeChirality ? numChiralBits : 0)));
}

template <>
std::uint32_t AtomPairArguments<std::uint32_t>::getResultSize() const {
  return (1 << (numAtomPairFingerprintBits +
                2 * (df_includeChirality ? numChiralBits : 0)));
}

template <typename OutputType>
AtomPairArguments<OutputType>::AtomPairArguments(
    const bool countSimulation, const bool includeChirality, const bool use2D,
    const unsigned int minDistance, const unsigned int maxDistance,
    const std::vector<std::uint32_t> countBounds,
    const std::uint32_t foldedSize)
    : FingerprintArguments<OutputType>(countSimulation, countBounds,
                                       foldedSize),
      df_includeChirality(includeChirality),
      df_use2D(use2D),
      d_minDistance(minDistance),
      d_maxDistance(maxDistance) {
  PRECONDITION(minDistance <= maxDistance, "bad distances provided");
}

template <typename OutputType>
std::string AtomPairArguments<OutputType>::infoString() const {
  return "AtomPairArguments includeChirality=" +
         std::to_string(df_includeChirality) +
         " use2D=" + std::to_string(df_use2D) +
         " minDistance=" + std::to_string(d_minDistance) +
         " maxDistance=" + std::to_string(d_maxDistance);
}

template <typename OutputType>
std::uint32_t AtomPairAtomEnv<OutputType>::getBitId(
    FingerprintArguments<OutputType> *arguments,
    const std::vector<std::uint32_t> *atomInvariants,
    const std::vector<std::uint32_t> *bondInvariants,
    const AdditionalOutput *additionalOutput) const {
  PRECONDITION((atomInvariants->size() >= d_atomIdFirst) &&
                   (atomInvariants->size() >= d_atomIdSecond),
               "bad atom invariants size");

  AtomPairArguments<OutputType> *atomPairArguments =
      dynamic_cast<AtomPairArguments<OutputType> *>(arguments);

  std::uint32_t codeSizeLimit =
      (1 << (codeSize +
             (atomPairArguments->df_includeChirality ? numChiralBits : 0))) -
      1;

  std::uint32_t atomCodeFirst =
      (*atomInvariants)[d_atomIdFirst] % codeSizeLimit;

  std::uint32_t atomCodeSecond =
      (*atomInvariants)[d_atomIdSecond] % codeSizeLimit;

  std::uint32_t bitId =
      getAtomPairCode(atomCodeFirst, atomCodeSecond, d_distance,
                      atomPairArguments->df_includeChirality);

  if (additionalOutput && additionalOutput->atomToBits) {
    additionalOutput->atomToBits->at(d_atomIdFirst).push_back(bitId);
    additionalOutput->atomToBits->at(d_atomIdSecond).push_back(bitId);
  }
  return bitId;
}

template <typename OutputType>
AtomPairAtomEnv<OutputType>::AtomPairAtomEnv(const unsigned int atomIdFirst,
                                             const unsigned int atomIdSecond,
                                             const unsigned int distance)
    : d_atomIdFirst(atomIdFirst),
      d_atomIdSecond(atomIdSecond),
      d_distance(distance) {}

template <typename OutputType>
std::vector<AtomEnvironment<OutputType> *>
AtomPairEnvGenerator<OutputType>::getEnvironments(
    const ROMol &mol, FingerprintArguments<OutputType> *arguments,
    const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput,
    const std::vector<std::uint32_t> *atomInvariants,
    const std::vector<std::uint32_t> *bondInvariants) const {
  const unsigned int atomCount = mol.getNumAtoms();
  PRECONDITION(!additionalOutput || !additionalOutput->atomToBits ||
                   additionalOutput->atomToBits->size() == atomCount,
               "bad atomToBits size in AdditionalOutput");

  AtomPairArguments<OutputType> *atomPairArguments =
      dynamic_cast<AtomPairArguments<OutputType> *>(arguments);
  std::vector<AtomEnvironment<OutputType> *> result =
      std::vector<AtomEnvironment<OutputType> *>();
  const double *distanceMatrix;
  if (atomPairArguments->df_use2D) {
    distanceMatrix = MolOps::getDistanceMat(mol);
  } else {
    distanceMatrix = MolOps::get3DDistanceMat(mol, confId);
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
        result.push_back(new AtomPairAtomEnv<OutputType>(i, j, distance));
      }
    }
  }

  return result;
}

template <typename OutputType>
std::string AtomPairEnvGenerator<OutputType>::infoString() const {
  return "AtomPairEnvironmentGenerator";
}

FingerprintGenerator<std::uint32_t> *getAtomPairGenerator(
    const unsigned int minDistance, const unsigned int maxDistance,
    const bool includeChirality, const bool use2D,
    const bool useCountSimulation,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    BondInvariantsGenerator *bondInvariantsGenerator) {
  AtomEnvironmentGenerator<std::uint32_t> *atomPairEnvGenerator =
      new AtomPair::AtomPairEnvGenerator<std::uint32_t>();
  FingerprintArguments<std::uint32_t> *atomPairArguments =
      new AtomPair::AtomPairArguments<std::uint32_t>(useCountSimulation,
                                                     includeChirality, use2D,
                                                     minDistance, maxDistance);

  bool ownsAtomInvGenerator = false;
  if (!atomInvariantsGenerator) {
    atomInvariantsGenerator = new AtomPairAtomInvGenerator(includeChirality);
    ownsAtomInvGenerator = true;
  }

  return new FingerprintGenerator<std::uint32_t>(
      atomPairEnvGenerator, atomPairArguments, atomInvariantsGenerator,
      bondInvariantsGenerator, ownsAtomInvGenerator, false);
}

}  // namespace AtomPair
}  // namespace RDKit
