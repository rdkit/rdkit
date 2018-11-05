//
//  Copyright (C) 2018 Boran Adas, Google Summer of Code
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
#include <limits.h>
#include <RDGeneral/hash/hash.hpp>
#include <RDGeneral/types.h>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>

#include <GraphMol/Fingerprints/FingerprintUtil.h>

namespace RDKit {
namespace RDKitFP {

std::vector<std::uint32_t> *RDKitFPAtomInvGenerator::getAtomInvariants(
    const ROMol &mol) const {
  std::vector<std::uint32_t> *result = new std::vector<std::uint32_t>();
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
OutputType RDKitFPArguments<OutputType>::getResultSize() const {
  return std::numeric_limits<OutputType>::max();
}

template <typename OutputType>
std::string RDKitFPArguments<OutputType>::infoString() const {
  return "RDKitFPArguments minPath=" + std::to_string(d_minPath) +
         " maxPath=" + std::to_string(d_maxPath) +
         " useHs=" + std::to_string(df_useHs) +
         " branchedPaths=" + std::to_string(df_branchedPaths) +
         " useBondOrder=" + std::to_string(df_useBondOrder);
}

template <typename OutputType>
RDKitFPArguments<OutputType>::RDKitFPArguments(
    unsigned int minPath, unsigned int maxPath, bool useHs, bool branchedPaths,
    bool useBondOrder, const bool countSimulation,
    const std::vector<std::uint32_t> countBounds, const std::uint32_t fpSize)
    : FingerprintArguments<OutputType>(countSimulation, countBounds, fpSize),
      d_minPath(minPath),
      d_maxPath(maxPath),
      df_useHs(useHs),
      df_branchedPaths(branchedPaths),
      df_useBondOrder(useBondOrder) {
  PRECONDITION(minPath != 0, "minPath==0");
  PRECONDITION(maxPath >= minPath, "maxPath<minPath");
}

template <typename OutputType>
OutputType RDKitFPAtomEnv<OutputType>::getBitId(
    FingerprintArguments<OutputType> *, // arguments
    const std::vector<std::uint32_t> *, // atomInvariants
    const std::vector<std::uint32_t> *, // bondInvariants
    const AdditionalOutput *, // additionalOutput
    const bool // hashResults
) const {
  // todo set additional outputs
  return d_bitId;
}

template <typename OutputType>
RDKitFPAtomEnv<OutputType>::RDKitFPAtomEnv(
    const OutputType bitId, const boost::dynamic_bitset<> atomsInPath)
    : d_bitId(bitId), d_atomsInPath(atomsInPath) {}

template <typename OutputType>
std::string RDKitFPEnvGenerator<OutputType>::infoString() const {
  return "RDKitFPEnvGenerator";
}

template <typename OutputType>
std::vector<AtomEnvironment<OutputType> *>
RDKitFPEnvGenerator<OutputType>::getEnvironments(
    const ROMol &mol, FingerprintArguments<OutputType> *arguments,
    const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *, // ignoreAtoms
    const int, // confId
    const AdditionalOutput *, // additionalOutput
    const std::vector<std::uint32_t> *atomInvariants,
    const std::vector<std::uint32_t> *, // bondInvariants
    const bool // hashResults
) const {
  PRECONDITION(!atomInvariants || atomInvariants->size() >= mol.getNumAtoms(),
               "bad atomInvariants size");

  RDKitFPArguments<OutputType> *rDKitFPArguments =
      dynamic_cast<RDKitFPArguments<OutputType> *>(arguments);

  std::vector<AtomEnvironment<OutputType> *> result;

  // get all paths
  INT_PATH_LIST_MAP allPaths;
  RDKitFPUtils::enumerateAllPaths(
      mol, allPaths, fromAtoms, rDKitFPArguments->df_branchedPaths,
      rDKitFPArguments->df_useHs, rDKitFPArguments->d_minPath,
      rDKitFPArguments->d_maxPath);

  // identify query bonds
  std::vector<short> isQueryBond(mol.getNumBonds(), 0);
  std::vector<const Bond *> bondCache;
  RDKitFPUtils::identifyQueryBonds(mol, bondCache, isQueryBond);

  boost::dynamic_bitset<> atomsInPath(mol.getNumAtoms());
  for (INT_PATH_LIST_MAP_CI paths = allPaths.begin(); paths != allPaths.end();
       paths++) {
    BOOST_FOREACH (const PATH_TYPE &path, paths->second) {
      // the bond hashes of the path
      std::vector<std::uint32_t> bondHashes = RDKitFPUtils::generateBondHashes(
          mol, atomsInPath, bondCache, isQueryBond, path,
          rDKitFPArguments->df_useBondOrder, atomInvariants);
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
          static_cast<OutputType>(seed), atomsInPath));
    }
  }

  return result;
}

template <typename OutputType>
FingerprintGenerator<OutputType> *getRDKitFPGenerator(
    const unsigned int minPath, const unsigned int maxPath, const bool useHs,
    const bool branchedPaths, const bool useBondOrder,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    const bool countSimulation, const std::vector<std::uint32_t> countBounds,
    const std::uint32_t fpSize, const bool ownsAtomInvGen) {
  AtomEnvironmentGenerator<OutputType> *envGenerator =
      new RDKitFPEnvGenerator<OutputType>();
  FingerprintArguments<OutputType> *arguments =
      new RDKitFPArguments<OutputType>(minPath, maxPath, useHs, branchedPaths,
                                       useBondOrder, countSimulation,
                                       countBounds, fpSize);

  bool ownsAtomInvGenerator = ownsAtomInvGen;
  if (!atomInvariantsGenerator) {
    atomInvariantsGenerator = new RDKitFPAtomInvGenerator();
    ownsAtomInvGenerator = true;
  }

  return new FingerprintGenerator<OutputType>(envGenerator, arguments,
                                              atomInvariantsGenerator, nullptr,
                                              ownsAtomInvGenerator, false);
}

template RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<std::uint32_t> *getRDKitFPGenerator(
    const unsigned int minPath, const unsigned int maxPath, const bool useHs,
    const bool branchedPaths, const bool useBondOrder,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    const bool countSimulation, const std::vector<std::uint32_t> countBounds,
    const std::uint32_t fpSize, const bool ownsAtomInvGen);

template RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<std::uint64_t> *getRDKitFPGenerator(
    const unsigned int minPath, const unsigned int maxPath, const bool useHs,
    const bool branchedPaths, const bool useBondOrder,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    const bool countSimulation, const std::vector<std::uint32_t> countBounds,
    const std::uint32_t fpSize, const bool ownsAtomInvGen);
}  // namespace RDKitFP
}  // namespace RDKit
