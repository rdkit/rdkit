//
//  Copyright (C) 2018-2022 Boran Adas and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/Fingerprints/RDKitFPGenerator.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <RDGeneral/hash/hash.hpp>
#include <cstdint>

#include <GraphMol/QueryOps.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/random.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <climits>
#include <RDGeneral/hash/hash.hpp>
#include <RDGeneral/types.h>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>

#include <GraphMol/Fingerprints/FingerprintUtil.h>

namespace RDKit {
namespace RDKitFP {

std::vector<std::uint32_t> *RDKitFPAtomInvGenerator::getAtomInvariants(
    const ROMol &mol) const {
  auto *result = new std::vector<std::uint32_t>();
  result->reserve(mol.getNumAtoms());
  for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
       atomIt != mol.endAtoms(); ++atomIt) {
    unsigned int aHash = ((*atomIt)->getAtomicNum() % 128) << 1 |
                         static_cast<unsigned int>((*atomIt)->getIsAromatic());
    result->push_back(aHash);
  }
  return result;
}

std::string RDKitFPAtomInvGenerator::infoString() const {
  return "RDKitFPAtomInvGenerator";
}

RDKitFPAtomInvGenerator *RDKitFPAtomInvGenerator::clone() const {
  return new RDKitFPAtomInvGenerator();
}

template <typename OutputType>
OutputType RDKitFPEnvGenerator<OutputType>::getResultSize() const {
  return std::numeric_limits<OutputType>::max();
}

std::string RDKitFPArguments::infoString() const {
  return "RDKitFPArguments minPath=" + std::to_string(d_minPath) +
         " maxPath=" + std::to_string(d_maxPath) +
         " useHs=" + std::to_string(df_useHs) +
         " branchedPaths=" + std::to_string(df_branchedPaths) +
         " useBondOrder=" + std::to_string(df_useBondOrder);
}

RDKitFPArguments::RDKitFPArguments(unsigned int minPath, unsigned int maxPath,
                                   bool useHs, bool branchedPaths,
                                   bool useBondOrder, bool countSimulation,
                                   const std::vector<std::uint32_t> countBounds,
                                   std::uint32_t fpSize,
                                   std::uint32_t numBitsPerFeature)
    : FingerprintArguments(countSimulation, countBounds, fpSize,
                           numBitsPerFeature),
      d_minPath(minPath),
      d_maxPath(maxPath),
      df_useHs(useHs),
      df_branchedPaths(branchedPaths),
      df_useBondOrder(useBondOrder) {
  PRECONDITION(minPath != 0, "minPath==0");
  PRECONDITION(maxPath >= minPath, "maxPath<minPath");
}

template <typename OutputType>
void RDKitFPAtomEnv<OutputType>::updateAdditionalOutput(
    AdditionalOutput *additionalOutput, size_t bitId) const {
  PRECONDITION(additionalOutput, "bad output pointer");
  if (additionalOutput->bitPaths) {
    (*additionalOutput->bitPaths)[bitId].push_back(d_bondPath);
  }
  if (additionalOutput->atomToBits || additionalOutput->atomCounts) {
    for (size_t i = 0; i < d_atomsInPath.size(); ++i) {
      if (d_atomsInPath[i]) {
        if (additionalOutput->atomToBits) {
          auto &alist = additionalOutput->atomToBits->at(i);
          if (std::find(alist.begin(), alist.end(), bitId) == alist.end()) {
            alist.push_back(bitId);
          }
        }
        if (additionalOutput->atomCounts) {
          additionalOutput->atomCounts->at(i)++;
        }
      }
    }
  }
}

template <typename OutputType>
OutputType RDKitFPAtomEnv<OutputType>::getBitId(
    FingerprintArguments *,              // arguments
    const std::vector<std::uint32_t> *,  // atomInvariants
    const std::vector<std::uint32_t> *,  // bondInvariants
    AdditionalOutput *,                  // additional Output
    const bool,                          // hashResults
    const std::uint64_t                  // fpSize
) const {
  return d_bitId;
}

template <typename OutputType>
std::string RDKitFPEnvGenerator<OutputType>::infoString() const {
  return "RDKitFPEnvGenerator";
}

template <typename OutputType>
std::vector<AtomEnvironment<OutputType> *>
RDKitFPEnvGenerator<OutputType>::getEnvironments(
    const ROMol &mol, FingerprintArguments *arguments,
    const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *,  // ignoreAtoms
    const int,                           // confId
    const AdditionalOutput *,            // additionalOutput
    const std::vector<std::uint32_t> *atomInvariants,
    const std::vector<std::uint32_t> *,  // bondInvariants
    const bool                           // hashResults
) const {
  PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
               "bad atomInvariants size");

  auto *fpArguments = dynamic_cast<RDKitFPArguments *>(arguments);

  std::vector<AtomEnvironment<OutputType> *> result;

  // get all paths
  INT_PATH_LIST_MAP allPaths;
  RDKitFPUtils::enumerateAllPaths(
      mol, allPaths, fromAtoms, fpArguments->df_branchedPaths,
      fpArguments->df_useHs, fpArguments->d_minPath, fpArguments->d_maxPath);

  // identify query bonds
  std::vector<short> isQueryBond(mol.getNumBonds(), 0);
  std::vector<const Bond *> bondCache;
  RDKitFPUtils::identifyQueryBonds(mol, bondCache, isQueryBond);

  boost::dynamic_bitset<> atomsInPath(mol.getNumAtoms());
  for (INT_PATH_LIST_MAP_CI paths = allPaths.begin(); paths != allPaths.end();
       paths++) {
    for (const auto &path : paths->second) {
      // the bond hashes of the path
      std::vector<std::uint32_t> bondHashes = RDKitFPUtils::generateBondHashes(
          mol, atomsInPath, bondCache, isQueryBond, path,
          fpArguments->df_useBondOrder, atomInvariants);
      if (!bondHashes.size()) {
        continue;
      }

      // hash the path to generate a seed:
      unsigned long seed;
      if (path.size() > 1) {
        std::sort(bondHashes.begin(), bondHashes.end());

        // finally, we will add the number of distinct atoms in the path at the
        // end
        // of the vect. This allows us to distinguish C1CC1 from CC(C)C
        bondHashes.push_back(static_cast<std::uint32_t>(atomsInPath.count()));
        seed = gboost::hash_range(bondHashes.begin(), bondHashes.end());
      } else {
        seed = bondHashes[0];
      }

      result.push_back(new RDKitFPAtomEnv<OutputType>(
          static_cast<OutputType>(seed), atomsInPath, path));
    }
  }

  return result;
}

template <typename OutputType>
FingerprintGenerator<OutputType> *getRDKitFPGenerator(
    unsigned int minPath, unsigned int maxPath, bool useHs, bool branchedPaths,
    bool useBondOrder, AtomInvariantsGenerator *atomInvariantsGenerator,
    bool countSimulation, const std::vector<std::uint32_t> countBounds,
    std::uint32_t fpSize, std::uint32_t numBitsPerFeature,
    bool ownsAtomInvGen) {
  auto *envGenerator = new RDKitFPEnvGenerator<OutputType>();
  auto *arguments = new RDKitFPArguments(
      minPath, maxPath, useHs, branchedPaths, useBondOrder, countSimulation,
      countBounds, fpSize, numBitsPerFeature);

  bool ownsAtomInvGenerator = ownsAtomInvGen;
  if (!atomInvariantsGenerator) {
    atomInvariantsGenerator = new RDKitFPAtomInvGenerator();
    ownsAtomInvGenerator = true;
  }

  return new FingerprintGenerator<OutputType>(envGenerator, arguments,
                                              atomInvariantsGenerator, nullptr,
                                              ownsAtomInvGenerator, false);
}

template RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<std::uint32_t>
    *getRDKitFPGenerator(unsigned int minPath, unsigned int maxPath, bool useHs,
                         bool branchedPaths, bool useBondOrder,
                         AtomInvariantsGenerator *atomInvariantsGenerator,
                         bool countSimulation,
                         const std::vector<std::uint32_t> countBounds,
                         std::uint32_t fpSize, std::uint32_t numBitsPerFeature,
                         bool ownsAtomInvGen);

template RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<std::uint64_t>
    *getRDKitFPGenerator(unsigned int minPath, unsigned int maxPath, bool useHs,
                         bool branchedPaths, bool useBondOrder,
                         AtomInvariantsGenerator *atomInvariantsGenerator,
                         bool countSimulation,
                         const std::vector<std::uint32_t> countBounds,
                         std::uint32_t fpSize, std::uint32_t numBitsPerFeature,
                         bool ownsAtomInvGen);

}  // namespace RDKitFP

}  // namespace RDKit
