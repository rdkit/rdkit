//
//  Copyright (C) 2018 Boran Adas, Google Summer of Code
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/AtomPairGenerator.h>
#include <GraphMol/Fingerprints/FingerprintUtil.h>
#include <RDGeneral/hash/hash.hpp>

namespace RDKit {
namespace AtomPair {
using namespace AtomPairs;

AtomPairAtomInvGenerator::AtomPairAtomInvGenerator(
    bool includeChirality, bool topologicalTorsionCorrection)
    : df_includeChirality(includeChirality),
      df_topologicalTorsionCorrection(topologicalTorsionCorrection) {}

std::vector<std::uint32_t> *AtomPairAtomInvGenerator::getAtomInvariants(
    const ROMol &mol) const {
  std::vector<std::uint32_t> *atomInvariants =
      new std::vector<std::uint32_t>(mol.getNumAtoms());

  for (ROMol::ConstAtomIterator atomItI = mol.beginAtoms();
       atomItI != mol.endAtoms(); ++atomItI) {
    (*atomInvariants)[(*atomItI)->getIdx()] =
        getAtomCode(*atomItI, 0, df_includeChirality) -
        (df_topologicalTorsionCorrection ? 2 : 0);
  }

  return atomInvariants;
}

std::string AtomPairAtomInvGenerator::infoString() const {
  return "AtomPairInvariantGenerator includeChirality=" +
         std::to_string(df_includeChirality) +
         " topologicalTorsionCorrection=" +
         std::to_string(df_topologicalTorsionCorrection);
}

AtomPairAtomInvGenerator *AtomPairAtomInvGenerator::clone() const {
  return new AtomPairAtomInvGenerator(df_includeChirality,
                                      df_topologicalTorsionCorrection);
}

template <typename OutputType>
OutputType AtomPairArguments<OutputType>::getResultSize() const {
  OutputType result = 1;
  return (result << (numAtomPairFingerprintBits +
                     2 * (df_includeChirality ? numChiralBits : 0)));
}

template <typename OutputType>
AtomPairArguments<OutputType>::AtomPairArguments(
    const bool countSimulation, const bool includeChirality, const bool use2D,
    const unsigned int minDistance, const unsigned int maxDistance,
    const std::vector<std::uint32_t> countBounds, const std::uint32_t fpSize)
    : FingerprintArguments<OutputType>(countSimulation, countBounds, fpSize),
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
OutputType AtomPairAtomEnv<OutputType>::getBitId(
    FingerprintArguments<OutputType> *arguments,
    const std::vector<std::uint32_t> *atomInvariants,
    const std::vector<std::uint32_t> *,  // bondInvariants
    const AdditionalOutput *additionalOutput, const bool hashResults,
    const std::uint64_t fpSize) const {
  PRECONDITION((atomInvariants->size() >= d_atomIdFirst) &&
                   (atomInvariants->size() >= d_atomIdSecond),
               "bad atom invariants size");

  auto *atomPairArguments =
      dynamic_cast<AtomPairArguments<OutputType> *>(arguments);

  std::uint32_t codeSizeLimit =
      (1 << (codeSize +
             (atomPairArguments->df_includeChirality ? numChiralBits : 0))) -
      1;

  std::uint32_t atomCodeFirst =
      (*atomInvariants)[d_atomIdFirst] % codeSizeLimit;

  std::uint32_t atomCodeSecond =
      (*atomInvariants)[d_atomIdSecond] % codeSizeLimit;

  std::uint32_t bitId = 0;
  if (hashResults) {
    gboost::hash_combine(bitId, std::min(atomCodeFirst, atomCodeSecond));
    gboost::hash_combine(bitId, d_distance);
    gboost::hash_combine(bitId, std::max(atomCodeFirst, atomCodeSecond));
  } else {
    bitId = getAtomPairCode(atomCodeFirst, atomCodeSecond, d_distance,
                            atomPairArguments->df_includeChirality);
  }

  if (additionalOutput) {
    std::uint32_t tBitId = bitId;
    if (fpSize) {
      tBitId = tBitId % fpSize;
    }
    if (additionalOutput->bitInfoMap) {
      (*additionalOutput->bitInfoMap)[tBitId].emplace_back(d_atomIdFirst,
                                                           d_atomIdSecond);
    }
    if (additionalOutput->atomToBits) {
      additionalOutput->atomToBits->at(d_atomIdFirst).push_back(tBitId);
      additionalOutput->atomToBits->at(d_atomIdSecond).push_back(tBitId);
    }
    if (additionalOutput->atomCounts) {
      additionalOutput->atomCounts->at(d_atomIdFirst)++;
      additionalOutput->atomCounts->at(d_atomIdSecond)++;
    }
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
    const std::vector<std::uint32_t> *,  // atomInvariants
    const std::vector<std::uint32_t> *,  // bondInvariants,
    const bool                           // hashResults
) const {
  const unsigned int atomCount = mol.getNumAtoms();
  PRECONDITION(!additionalOutput || !additionalOutput->atomToBits ||
                   additionalOutput->atomToBits->size() == atomCount,
               "bad atomToBits size in AdditionalOutput");

  auto *atomPairArguments =
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
      auto distance =
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

template <typename OutputType>
FingerprintGenerator<OutputType> *getAtomPairGenerator(
    const unsigned int minDistance, const unsigned int maxDistance,
    const bool includeChirality, const bool use2D,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    const bool useCountSimulation, const std::uint32_t fpSize,
    const std::vector<std::uint32_t> countBounds, const bool ownsAtomInvGen) {
  AtomEnvironmentGenerator<OutputType> *atomPairEnvGenerator =
      new AtomPair::AtomPairEnvGenerator<OutputType>();
  FingerprintArguments<OutputType> *atomPairArguments =
      new AtomPair::AtomPairArguments<OutputType>(
          useCountSimulation, includeChirality, use2D, minDistance, maxDistance,
          countBounds, fpSize);

  bool ownsAtomInvGenerator = ownsAtomInvGen;
  if (!atomInvariantsGenerator) {
    atomInvariantsGenerator = new AtomPairAtomInvGenerator(includeChirality);
    ownsAtomInvGenerator = true;
  }

  return new FingerprintGenerator<OutputType>(
      atomPairEnvGenerator, atomPairArguments, atomInvariantsGenerator, nullptr,
      ownsAtomInvGenerator, false);
}

template RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<std::uint32_t> *
getAtomPairGenerator(const unsigned int minDistance,
                     const unsigned int maxDistance,
                     const bool includeChirality, const bool use2D,
                     AtomInvariantsGenerator *atomInvariantsGenerator,
                     const bool useCountSimulation, const std::uint32_t fpSize,
                     const std::vector<std::uint32_t> countBounds,
                     const bool ownsAtomInvGen);

template RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<std::uint64_t> *
getAtomPairGenerator(const unsigned int minDistance,
                     const unsigned int maxDistance,
                     const bool includeChirality, const bool use2D,
                     AtomInvariantsGenerator *atomInvariantsGenerator,
                     const bool useCountSimulation, const std::uint32_t fpSize,
                     const std::vector<std::uint32_t> countBounds,
                     const bool ownsAtomInvGen);
}  // namespace AtomPair
}  // namespace RDKit
