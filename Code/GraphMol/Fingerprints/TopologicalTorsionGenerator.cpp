//
//  Copyright (C) 2018 Boran Adas, Google Summer of Code
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/Fingerprints/TopologicalTorsionGenerator.h>
#include <GraphMol/Fingerprints/FingerprintUtil.h>
#include <GraphMol/Fingerprints/AtomPairGenerator.h>

namespace RDKit {
namespace TopologicalTorsion {

using namespace AtomPairs;

template <typename OutputType>
TopologicalTorsionArguments<OutputType>::TopologicalTorsionArguments(
    const bool includeChirality, const uint32_t torsionAtomCount,
    const bool countSimulation, const std::vector<std::uint32_t> countBounds,
    const std::uint32_t fpSize)
    : FingerprintArguments<OutputType>(countSimulation, countBounds, fpSize),
      df_includeChirality(includeChirality),
      d_torsionAtomCount(torsionAtomCount){};

template <typename OutputType>
OutputType TopologicalTorsionArguments<OutputType>::getResultSize() const {
  OutputType result = 1;
  return (result << (d_torsionAtomCount *
                     (codeSize + (df_includeChirality ? numChiralBits : 0))));
};

template <typename OutputType>
std::string TopologicalTorsionArguments<OutputType>::infoString() const {
  return "TopologicalTorsionArguments includeChirality=" +
         std::to_string(df_includeChirality) +
         " torsionAtomCount=" + std::to_string(d_torsionAtomCount);
};
template <typename OutputType>
OutputType TopologicalTorsionAtomEnv<OutputType>::getBitId(
    FingerprintArguments<OutputType> *,  // arguments
    const std::vector<std::uint32_t> *,  // atomInvariants
    const std::vector<std::uint32_t> *,  // bondInvariants
    const AdditionalOutput *additionalOutput,
    const bool,  // hashResults
    const std::uint64_t fpSize) const {
  if (additionalOutput) {
    OutputType bitId = d_bitId;
    if (fpSize) {
      bitId %= fpSize;
    }
    if (additionalOutput->atomToBits || additionalOutput->atomCounts) {
      for (auto aid : d_atomPath) {
        if (additionalOutput->atomToBits) {
          additionalOutput->atomToBits->at(aid).push_back(bitId);
        }
        if (additionalOutput->atomCounts) {
          additionalOutput->atomCounts->at(aid)++;
        }
      }
    }
    if (additionalOutput->bitPaths) {
      (*additionalOutput->bitPaths)[bitId].push_back(d_atomPath);
    }
  }

  return d_bitId;
};

template <typename OutputType>
std::vector<AtomEnvironment<OutputType> *>
TopologicalTorsionEnvGenerator<OutputType>::getEnvironments(
    const ROMol &mol, FingerprintArguments<OutputType> *arguments,
    const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms,
    const int,                 // confId
    const AdditionalOutput *,  // additionalOutput
    const std::vector<std::uint32_t> *atomInvariants,
    const std::vector<std::uint32_t> *,  // bondInvariants
    const bool hashResults) const {
  auto *topologicalTorsionArguments =
      dynamic_cast<TopologicalTorsionArguments<OutputType> *>(arguments);

  std::vector<AtomEnvironment<OutputType> *> result;

  boost::dynamic_bitset<> *fromAtomsBV = nullptr;
  if (fromAtoms) {
    fromAtomsBV = new boost::dynamic_bitset<>(mol.getNumAtoms());
    for (auto fAt : *fromAtoms) {
      fromAtomsBV->set(fAt);
    }
  }
  boost::dynamic_bitset<> *ignoreAtomsBV = nullptr;
  if (ignoreAtoms) {
    ignoreAtomsBV = new boost::dynamic_bitset<>(mol.getNumAtoms());
    for (auto fAt : *ignoreAtoms) {
      ignoreAtomsBV->set(fAt);
    }
  }
  boost::dynamic_bitset<> pAtoms(mol.getNumAtoms());
  PATH_LIST paths = findAllPathsOfLengthN(
      mol, topologicalTorsionArguments->d_torsionAtomCount, false);
  for (PATH_LIST::const_iterator pathIt = paths.begin(); pathIt != paths.end();
       ++pathIt) {
    bool keepIt = true;
    if (fromAtomsBV) {
      keepIt = false;
    }
    std::vector<std::uint32_t> pathCodes;
    const PATH_TYPE &path = *pathIt;
    if (fromAtomsBV) {
      if (fromAtomsBV->test(static_cast<std::uint32_t>(path.front())) ||
          fromAtomsBV->test(static_cast<std::uint32_t>(path.back()))) {
        keepIt = true;
      }
    }
    if (keepIt && ignoreAtomsBV) {
      for (int pElem : path) {
        if (ignoreAtomsBV->test(pElem)) {
          keepIt = false;
          break;
        }
      }
    }
    if (keepIt) {
      pAtoms.reset();
      for (auto pIt = path.begin(); pIt < path.end(); ++pIt) {
        // look for a cycle that doesn't start at the first atom
        // we can't effectively canonicalize these at the moment
        // (was github #811)
        if (pIt != path.begin() && *pIt != *(path.begin()) && pAtoms[*pIt]) {
          pathCodes.clear();
          break;
        }
        pAtoms.set(*pIt);
        unsigned int code = (*atomInvariants)[*pIt] % ((1 << codeSize) - 1) + 1;
        // subtract off the branching number:
        if (pIt != path.begin() && pIt + 1 != path.end()) {
          --code;
        }
        pathCodes.push_back(code);
      }
      if (pathCodes.size()) {
        OutputType code;
        if (hashResults) {
          code = getTopologicalTorsionHash(pathCodes);
        } else {
          code = getTopologicalTorsionCode(
              pathCodes, topologicalTorsionArguments->df_includeChirality);
        }
        result.push_back(new TopologicalTorsionAtomEnv<OutputType>(code, path));
      }
    }
  }
  delete fromAtomsBV;
  delete ignoreAtomsBV;

  return result;
};

template <typename OutputType>
std::string TopologicalTorsionEnvGenerator<OutputType>::infoString() const {
  return "TopologicalTorsionEnvGenerator";
};

template <typename OutputType>
FingerprintGenerator<OutputType> *getTopologicalTorsionGenerator(
    const bool includeChirality, const uint32_t torsionAtomCount,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    const bool countSimulation, const std::vector<std::uint32_t> countBounds,
    const std::uint32_t fpSize, const bool ownsAtomInvGen) {
  auto *envGenerator = new TopologicalTorsionEnvGenerator<OutputType>();

  auto *arguments = new TopologicalTorsionArguments<OutputType>(
      includeChirality, torsionAtomCount, countSimulation, countBounds, fpSize);

  bool ownsAtomInvGenerator = ownsAtomInvGen;
  if (!atomInvariantsGenerator) {
    atomInvariantsGenerator =
        new AtomPair::AtomPairAtomInvGenerator(includeChirality, true);
    ownsAtomInvGenerator = true;
  }

  return new FingerprintGenerator<OutputType>(envGenerator, arguments,
                                              atomInvariantsGenerator, nullptr,
                                              ownsAtomInvGenerator, false);
};

// Topological torsion fingerprint does not support 32 bit output yet

template RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<std::uint64_t> *
getTopologicalTorsionGenerator(const bool includeChirality,
                               const uint32_t torsionAtomCount,
                               AtomInvariantsGenerator *atomInvariantsGenerator,
                               const bool countSimulation,
                               const std::vector<std::uint32_t> countBounds,
                               const std::uint32_t fpSize,
                               const bool ownsAtomInvGen);

}  // namespace TopologicalTorsion
}  // namespace RDKit
