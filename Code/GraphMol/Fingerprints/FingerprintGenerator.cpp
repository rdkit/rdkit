//
//  Copyright (C) 2018-2022 Boran Adas and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <DataStructs/SparseIntVect.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/SparseBitVect.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <RDGeneral/hash/hash.hpp>
#include <cstdint>

#include <GraphMol/Fingerprints/AtomPairGenerator.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/Fingerprints/RDKitFPGenerator.h>
#include <GraphMol/Fingerprints/TopologicalTorsionGenerator.h>

#include <RDGeneral/RDThreads.h>
#ifdef RDK_BUILD_THREADSAFE_SSS
#include <thread>
#include <future>
#endif

namespace RDKit {

FingerprintArguments::FingerprintArguments(
    const bool countSimulation, const std::vector<std::uint32_t> countBounds,
    std::uint32_t fpSize, std::uint32_t numBitsPerFeature,
    bool includeChirality)
    : df_countSimulation(countSimulation),
      df_includeChirality(includeChirality),
      d_countBounds(countBounds),
      d_fpSize(fpSize),
      d_numBitsPerFeature(numBitsPerFeature) {
  PRECONDITION(!countSimulation || !countBounds.empty(),
               "bad count bounds provided");
  PRECONDITION(d_numBitsPerFeature > 0, "numBitsPerFeature must be >0");
}

std::string FingerprintArguments::commonArgumentsString() const {
  return "Common arguments : countSimulation=" +
         std::to_string(df_countSimulation) +
         " fpSize=" + std::to_string(d_fpSize) +
         " bitsPerFeature=" + std::to_string(d_numBitsPerFeature) +
         " includeChirality=" + std::to_string(df_includeChirality);
}

template <typename OutputType>
FingerprintGenerator<OutputType>::FingerprintGenerator(
    AtomEnvironmentGenerator<OutputType> *atomEnvironmentGenerator,
    FingerprintArguments *fingerprintArguments,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    BondInvariantsGenerator *bondInvariantsGenerator, bool ownsAtomInvGenerator,
    bool ownsBondInvGenerator)
    : df_ownsAtomInvGenerator(ownsAtomInvGenerator),
      df_ownsBondInvGenerator(ownsBondInvGenerator) {
  this->dp_atomEnvironmentGenerator = atomEnvironmentGenerator;
  this->dp_atomEnvironmentGenerator->dp_fingerprintArguments =
      fingerprintArguments;

  this->dp_fingerprintArguments = fingerprintArguments;
  this->dp_atomInvariantsGenerator = atomInvariantsGenerator;
  this->dp_bondInvariantsGenerator = bondInvariantsGenerator;
}

template FingerprintGenerator<std::uint32_t>::FingerprintGenerator(
    AtomEnvironmentGenerator<std::uint32_t> *atomEnvironmentGenerator,
    FingerprintArguments *fingerprintArguments,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    BondInvariantsGenerator *bondInvariantsGenerator, bool ownsAtomInvGenerator,
    bool ownsBondInvGenerator);

template FingerprintGenerator<std::uint64_t>::FingerprintGenerator(
    AtomEnvironmentGenerator<std::uint64_t> *atomEnvironmentGenerator,
    FingerprintArguments *fingerprintArguments,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    BondInvariantsGenerator *bondInvariantsGenerator, bool ownsAtomInvGenerator,
    bool ownsBondInvGenerator);

template <typename OutputType>
FingerprintGenerator<OutputType>::~FingerprintGenerator() {
  delete dp_atomEnvironmentGenerator;
  delete dp_fingerprintArguments;
  if (df_ownsAtomInvGenerator) {
    delete dp_atomInvariantsGenerator;
  }
  if (df_ownsBondInvGenerator) {
    delete dp_bondInvariantsGenerator;
  }
}

namespace {
void reinitAdditionalOutput(AdditionalOutput &ao, size_t numAtoms) {
  if (ao.atomCounts) {
    ao.atomCounts->resize(numAtoms);
    std::fill(ao.atomCounts->begin(), ao.atomCounts->end(), 0);
  }
  if (ao.atomToBits) {
    ao.atomToBits->resize(numAtoms);
    std::fill(ao.atomToBits->begin(), ao.atomToBits->end(),
              std::vector<std::uint64_t>());
  }
  if (ao.bitInfoMap) {
    ao.bitInfoMap->clear();
  }
  if (ao.bitPaths) {
    ao.bitPaths->clear();
  }
}
}  // namespace

template FingerprintGenerator<std::uint32_t>::~FingerprintGenerator();

template FingerprintGenerator<std::uint64_t>::~FingerprintGenerator();

template std::string FingerprintGenerator<std::uint32_t>::infoString() const;

template std::string FingerprintGenerator<std::uint64_t>::infoString() const;

template <typename OutputType>
std::string FingerprintGenerator<OutputType>::infoString() const {
  std::string separator = " --- ";
  return dp_fingerprintArguments->commonArgumentsString() + separator +
         dp_fingerprintArguments->infoString() + separator +
         dp_atomEnvironmentGenerator->infoString() + separator +
         (dp_atomInvariantsGenerator
              ? (dp_atomInvariantsGenerator->infoString() + separator)
              : ("No atom invariants generator" + separator)) +
         (dp_bondInvariantsGenerator
              ? (dp_bondInvariantsGenerator->infoString())
              : "No bond invariants generator");
}

template <typename OutputType>
std::unique_ptr<SparseIntVect<OutputType>>
FingerprintGenerator<OutputType>::getFingerprintHelper(
    const ROMol &mol, FingerprintFuncArguments &args,
    const std::uint64_t fpSize) const {
  const ROMol *lmol = &mol;
  std::unique_ptr<ROMol> tmol;
  if (dp_fingerprintArguments->df_includeChirality &&
      !mol.hasProp(common_properties::_StereochemDone)) {
    tmol = std::unique_ptr<ROMol>(new ROMol(mol));
    MolOps::assignStereochemistry(*tmol);
    lmol = tmol.get();
  }

  if (args.additionalOutput) {
    reinitAdditionalOutput(*args.additionalOutput, mol.getNumAtoms());
  }

  bool hashResults = false;
  if (fpSize != 0) {
    hashResults = true;
  }

  std::unique_ptr<std::vector<std::uint32_t>> atomInvariants = nullptr;
  if (args.customAtomInvariants) {
    atomInvariants.reset(
        new std::vector<std::uint32_t>(*args.customAtomInvariants));
  } else if (dp_atomInvariantsGenerator) {
    atomInvariants.reset(dp_atomInvariantsGenerator->getAtomInvariants(mol));
  }

  std::unique_ptr<std::vector<std::uint32_t>> bondInvariants = nullptr;
  if (args.customBondInvariants) {
    bondInvariants.reset(
        new std::vector<std::uint32_t>(*args.customBondInvariants));
  } else if (dp_bondInvariantsGenerator) {
    bondInvariants.reset(dp_bondInvariantsGenerator->getBondInvariants(mol));
  }

  // create all atom environments that will generate the bit-ids that will make
  // up the fingerprint
  auto atomEnvironments = dp_atomEnvironmentGenerator->getEnvironments(
      *lmol, dp_fingerprintArguments, args.fromAtoms, args.ignoreAtoms,
      args.confId, args.additionalOutput, atomInvariants.get(),
      bondInvariants.get(), hashResults);

  // allocate the result
  auto res = std::make_unique<SparseIntVect<OutputType>>(
      fpSize ? fpSize : dp_atomEnvironmentGenerator->getResultSize());

  // define a mersenne twister with customized parameters.
  // The standard parameters (used to create boost::mt19937)
  // result in an RNG that's much too computationally intensive
  // to seed.
  // These are the parameters that have been used for the RDKit fingerprint.
  typedef boost::random::mersenne_twister<std::uint32_t, 32, 4, 2, 31,
                                          0x9908b0df, 11, 7, 0x9d2c5680, 15,
                                          0xefc60000, 18, 3346425566U>
      rng_type;
  typedef boost::uniform_int<> distrib_type;
  typedef boost::variate_generator<rng_type &, distrib_type> source_type;
  std::unique_ptr<rng_type> generator;
  //
  // if we generate arbitrarily sized ints then mod them down to the
  // appropriate size, we can guarantee that a fingerprint of
  // size x has the same bits set as one of size 2x that's been folded
  // in half.  This is a nice guarantee to have.
  //
  std::unique_ptr<distrib_type> dist;
  std::unique_ptr<source_type> randomSource;
  if (dp_fingerprintArguments->d_numBitsPerFeature > 1) {
    // we will only create the RNG if we're going to need it
    generator.reset(new rng_type(42u));
    dist.reset(new distrib_type(0, INT_MAX));
    randomSource.reset(new source_type(*generator, *dist));
  }

  // iterate over every atom environment and generate bit-ids that will make up
  // the fingerprint
  for (const auto env : atomEnvironments) {
    OutputType seed = env->getBitId(dp_fingerprintArguments,
                                    atomInvariants.get(), bondInvariants.get(),
                                    args.additionalOutput, hashResults, fpSize);

    auto bitId = seed;
    if (fpSize != 0) {
      bitId %= fpSize;
    }
    res->setVal(bitId, res->getVal(bitId) + 1);
    if (args.additionalOutput) {
      env->updateAdditionalOutput(args.additionalOutput, bitId);
    }
    // do the additional bits if required:
    if (dp_fingerprintArguments->d_numBitsPerFeature > 1) {
      generator->seed(static_cast<rng_type::result_type>(seed));

      for (boost::uint32_t bitN = 1;
           bitN < dp_fingerprintArguments->d_numBitsPerFeature; ++bitN) {
        bitId = (*randomSource)();
        if (fpSize != 0) {
          bitId %= fpSize;
        }
        res->setVal(bitId, res->getVal(bitId) + 1);
        if (args.additionalOutput) {
          env->updateAdditionalOutput(args.additionalOutput, bitId);
        }
      }
    }
    delete env;
  }

  return res;
}
namespace {
template <typename OutputType>
void duplicateAdditionalOutputBit(AdditionalOutput &oldAO,
                                  AdditionalOutput &newAO, OutputType origBitId,
                                  OutputType newBitId) {
  PRECONDITION(!((oldAO.bitInfoMap != nullptr) ^ (newAO.bitInfoMap != nullptr)),
               "bitInfoMap not allocated");
  PRECONDITION(!((oldAO.atomToBits != nullptr) ^ (newAO.atomToBits != nullptr)),
               "atomToBits not allocated");
  PRECONDITION(!((oldAO.bitPaths != nullptr) ^ (newAO.bitPaths != nullptr)),
               "bitPaths not allocated");

  // we don't need to do anything with atomCounts

  if (oldAO.atomToBits) {
    if (newAO.atomToBits->empty()) {
      newAO.atomToBits->resize(oldAO.atomToBits->size());
    }
    for (unsigned int i = 0; i < oldAO.atomToBits->size(); ++i) {
      const auto &nv = oldAO.atomToBits->at(i);
      if (std::find(nv.begin(), nv.end(), origBitId) != nv.end()) {
        newAO.atomToBits->at(i).push_back(newBitId);
      }
    }
  }
  if (oldAO.bitInfoMap) {
    const auto v = oldAO.bitInfoMap->find(origBitId);
    if (v != oldAO.bitInfoMap->end()) {
      (*newAO.bitInfoMap)[newBitId] = v->second;
    }
  }
  if (oldAO.bitPaths) {
    const auto v = oldAO.bitPaths->find(origBitId);
    if (v != oldAO.bitPaths->end()) {
      (*newAO.bitPaths)[newBitId] = v->second;
    }
  }
}

void setupTempAdditionalOutput(RDKit::FingerprintFuncArguments &args,
                               AdditionalOutput &countSimulationOutput,
                               size_t numAtoms) {
  if (args.additionalOutput->atomToBits) {
    countSimulationOutput.allocateAtomToBits();
  }
  if (args.additionalOutput->atomCounts) {
    countSimulationOutput.allocateAtomCounts();
  }
  if (args.additionalOutput->bitInfoMap) {
    countSimulationOutput.allocateBitInfoMap();
  }
  if (args.additionalOutput->bitPaths) {
    countSimulationOutput.allocateBitPaths();
  }
  reinitAdditionalOutput(*args.additionalOutput, numAtoms);
}
}  // namespace

template <typename OutputType>
std::unique_ptr<SparseIntVect<OutputType>>
FingerprintGenerator<OutputType>::getSparseCountFingerprint(
    const ROMol &mol, FingerprintFuncArguments &args) const {
  return getFingerprintHelper(mol, args);
}

// todo getSparseFingerprint does not completely produce the same output as
// getSparseCountFingerprint. Count simulation and potential 64 bit outputs
// makes size limiting necessary for getSparseFingerprint. This can be
// changed if there is another way to avoid the size limitation of SparseBitVect
template <typename OutputType>
std::unique_ptr<SparseBitVect>
FingerprintGenerator<OutputType>::getSparseFingerprint(
    const ROMol &mol, FingerprintFuncArguments &args) const {
  // make sure the result will fit into SparseBitVect
  std::uint32_t resultSize =
      std::min((std::uint64_t)std::numeric_limits<std::uint32_t>::max(),
               (std::uint64_t)dp_atomEnvironmentGenerator->getResultSize());

  std::uint32_t effectiveSize = resultSize;
  if (dp_fingerprintArguments->df_countSimulation) {
    // effective size needs to be smaller than result size to compansate for
    // count simulation
    effectiveSize /= dp_fingerprintArguments->d_countBounds.size();
  }

  AdditionalOutput countSimulationOutput;
  AdditionalOutput *origAO = nullptr;
  if (dp_fingerprintArguments->df_countSimulation && args.additionalOutput) {
    setupTempAdditionalOutput(args, countSimulationOutput, mol.getNumAtoms());
    origAO = args.additionalOutput;
    args.additionalOutput = &countSimulationOutput;
  }

  auto tempResult = getFingerprintHelper(mol, args, effectiveSize);

  auto result = std::make_unique<SparseBitVect>(resultSize);

  for (auto val : tempResult->getNonzeroElements()) {
    if (dp_fingerprintArguments->df_countSimulation) {
      for (unsigned int i = 0;
           i < dp_fingerprintArguments->d_countBounds.size(); ++i) {
        // for every bound in the d_countBounds in dp_fingerprintArguments, set
        // a bit if the occurrence count is equal or higher than the bound for
        // that bit
        const auto &bounds_count = dp_fingerprintArguments->d_countBounds;
        if (val.second >= static_cast<int>(bounds_count[i])) {
          OutputType nBitId = val.first * bounds_count.size() + i;
          result->setBit(nBitId);
          if (args.additionalOutput) {
            duplicateAdditionalOutputBit(*args.additionalOutput, *origAO,
                                         static_cast<OutputType>(val.first),
                                         nBitId);
          }
        }
      }
    } else {
      result->setBit(val.first);
    }
  }
  if (origAO) {
    if (origAO->atomCounts) {
      *origAO->atomCounts = *countSimulationOutput.atomCounts;
    }
    args.additionalOutput = origAO;
  }

  return result;
}

template <typename OutputType>
std::unique_ptr<SparseIntVect<std::uint32_t>>
FingerprintGenerator<OutputType>::getCountFingerprint(
    const ROMol &mol, FingerprintFuncArguments &args) const {
  auto tempResult =
      getFingerprintHelper(mol, args, dp_fingerprintArguments->d_fpSize);

  auto result = std::make_unique<SparseIntVect<std::uint32_t>>(
      dp_fingerprintArguments->d_fpSize);
  for (auto val : tempResult->getNonzeroElements()) {
    result->setVal(val.first, val.second);
  }

  return result;
}

template <typename OutputType>
std::unique_ptr<ExplicitBitVect>
FingerprintGenerator<OutputType>::getFingerprint(
    const ROMol &mol, FingerprintFuncArguments &args) const {
  std::uint32_t effectiveSize = dp_fingerprintArguments->d_fpSize;
  if (dp_fingerprintArguments->df_countSimulation) {
    if (dp_fingerprintArguments->d_countBounds.empty()) {
      throw ValueErrorException("Count bounds are empty");
    }

    if (dp_fingerprintArguments->d_countBounds.size() >= effectiveSize) {
      throw ValueErrorException("Count bounds size is >= fingerprint size");
    }

    // effective size needs to be smaller than result size to compensate for
    // count simulation
    effectiveSize /= dp_fingerprintArguments->d_countBounds.size();
  }

  AdditionalOutput countSimulationOutput;
  AdditionalOutput *origAO = nullptr;
  if (dp_fingerprintArguments->df_countSimulation && args.additionalOutput) {
    setupTempAdditionalOutput(args, countSimulationOutput, mol.getNumAtoms());
    origAO = args.additionalOutput;
    args.additionalOutput = &countSimulationOutput;
  }
  auto tempResult = getFingerprintHelper(mol, args, effectiveSize);

  auto result =
      std::make_unique<ExplicitBitVect>(dp_fingerprintArguments->d_fpSize);
  for (auto val : tempResult->getNonzeroElements()) {
    if (dp_fingerprintArguments->df_countSimulation) {
      for (unsigned int i = 0;
           i < dp_fingerprintArguments->d_countBounds.size(); ++i) {
        // for every bound in the d_countBounds in dp_fingerprintArguments,
        // set a bit if the occurrence count is equal or higher than the bound
        // for that bit
        const auto &bounds_count = dp_fingerprintArguments->d_countBounds;
        if (val.second >= static_cast<int>(bounds_count[i])) {
          OutputType nBitId = val.first * bounds_count.size() + i;
          result->setBit(nBitId);
          if (args.additionalOutput) {
            duplicateAdditionalOutputBit(*args.additionalOutput, *origAO,
                                         static_cast<OutputType>(val.first),
                                         nBitId);
          }
        }
      }
    } else {
      result->setBit(val.first);
    }
  }

  if (origAO) {
    if (origAO->atomCounts) {
      *origAO->atomCounts = *countSimulationOutput.atomCounts;
    }
    args.additionalOutput = origAO;
  }

  return result;
}

namespace {
template <typename ReturnType, typename FuncType>
std::vector<std::unique_ptr<ReturnType>> mtgetFingerprints(
    FuncType func, const std::vector<const ROMol *> &mols, int numThreads) {
  std::vector<std::uint32_t> *fromAtoms = nullptr;
  std::vector<std::uint32_t> *ignoreAtoms = nullptr;
  std::vector<std::uint32_t> *customAtomInvariants = nullptr;
  std::vector<std::uint32_t> *customBondInvariants = nullptr;
  int confId = -1;
  AdditionalOutput *additionalOutput = nullptr;
  FingerprintFuncArguments args(fromAtoms, ignoreAtoms, confId,
                                additionalOutput, customAtomInvariants,
                                customBondInvariants);

  std::vector<std::unique_ptr<ReturnType>> result;
  auto numThreadsToUse = getNumThreadsToUse(numThreads);
  unsigned int nmols = mols.size();
  result.reserve(nmols);
  if (numThreadsToUse == 1) {
    for (auto i = 0u; i < nmols; ++i) {
      if (!mols[i]) {
        result.emplace_back(std::unique_ptr<ReturnType>());
      } else {
        result.emplace_back(std::move(func(*mols[i], args)));
      }
    }
  }
#ifdef RDK_BUILD_THREADSAFE_SSS
  else {
    std::vector<std::vector<std::unique_ptr<ReturnType>>> accum(
        numThreadsToUse);
    std::vector<std::thread> tg;
    for (auto ti = 0u; ti < numThreadsToUse; ++ti) {
      auto lfunc = [&](unsigned int tidx) {
        for (auto midx = tidx; midx < mols.size(); midx += numThreadsToUse) {
          if (!mols[midx]) {
            accum[tidx].emplace_back(std::unique_ptr<ReturnType>());
          } else {
            accum[tidx].emplace_back(std::move(func(*mols[midx], args)));
          }
        }
      };
      tg.emplace_back(std::thread(lfunc, ti));
    }
    for (auto &thread : tg) {
      if (thread.joinable()) {
        thread.join();
      }
    }
    for (auto midx = 0u; midx < mols.size(); ++midx) {
      auto tidx = midx % numThreadsToUse;
      auto jidx = midx / numThreadsToUse;
      result.emplace_back(std::move(accum[tidx][jidx]));
    }
  }
#endif
  return result;
}
}  // namespace

template <typename OutputType>
std::vector<std::unique_ptr<ExplicitBitVect>>
FingerprintGenerator<OutputType>::getFingerprints(
    const std::vector<const ROMol *> &mols, int numThreads) const {
  auto fpfunc = [&](const ROMol &mol, FingerprintFuncArguments &args) {
    return this->getFingerprint(mol, args);
  };
  return mtgetFingerprints<ExplicitBitVect, decltype(fpfunc)>(fpfunc, mols,
                                                              numThreads);
}

template <typename OutputType>
std::vector<std::unique_ptr<SparseBitVect>>
FingerprintGenerator<OutputType>::getSparseFingerprints(
    const std::vector<const ROMol *> &mols, int numThreads) const {
  auto fpfunc = [&](const ROMol &mol, FingerprintFuncArguments &args) {
    return this->getSparseFingerprint(mol, args);
  };
  return mtgetFingerprints<SparseBitVect, decltype(fpfunc)>(fpfunc, mols,
                                                            numThreads);
}

template <typename OutputType>
std::vector<std::unique_ptr<SparseIntVect<std::uint32_t>>>
FingerprintGenerator<OutputType>::getCountFingerprints(
    const std::vector<const ROMol *> &mols, int numThreads) const {
  auto fpfunc = [&](const ROMol &mol, FingerprintFuncArguments &args) {
    return this->getCountFingerprint(mol, args);
  };
  return mtgetFingerprints<SparseIntVect<std::uint32_t>, decltype(fpfunc)>(
      fpfunc, mols, numThreads);
}

template <typename OutputType>
std::vector<std::unique_ptr<SparseIntVect<OutputType>>>
FingerprintGenerator<OutputType>::getSparseCountFingerprints(
    const std::vector<const ROMol *> &mols, int numThreads) const {
  auto fpfunc = [&](const ROMol &mol, FingerprintFuncArguments &args) {
    return this->getSparseCountFingerprint(mol, args);
  };
  return mtgetFingerprints<SparseIntVect<OutputType>, decltype(fpfunc)>(
      fpfunc, mols, numThreads);
}

template RDKIT_FINGERPRINTS_EXPORT std::unique_ptr<SparseIntVect<std::uint32_t>>
FingerprintGenerator<std::uint32_t>::getSparseCountFingerprint(
    const ROMol &mol, FingerprintFuncArguments &args) const;

template RDKIT_FINGERPRINTS_EXPORT std::unique_ptr<SparseIntVect<std::uint64_t>>
FingerprintGenerator<std::uint64_t>::getSparseCountFingerprint(
    const ROMol &mol, FingerprintFuncArguments &args) const;

template RDKIT_FINGERPRINTS_EXPORT std::unique_ptr<SparseBitVect>
FingerprintGenerator<std::uint32_t>::getSparseFingerprint(
    const ROMol &mol, FingerprintFuncArguments &args) const;

template RDKIT_FINGERPRINTS_EXPORT std::unique_ptr<SparseBitVect>
FingerprintGenerator<std::uint64_t>::getSparseFingerprint(
    const ROMol &mol, FingerprintFuncArguments &args) const;

template RDKIT_FINGERPRINTS_EXPORT std::unique_ptr<SparseIntVect<std::uint32_t>>
FingerprintGenerator<std::uint32_t>::getCountFingerprint(
    const ROMol &mol, FingerprintFuncArguments &args) const;

template RDKIT_FINGERPRINTS_EXPORT std::unique_ptr<SparseIntVect<std::uint32_t>>
FingerprintGenerator<std::uint64_t>::getCountFingerprint(
    const ROMol &mol, FingerprintFuncArguments &args) const;

template RDKIT_FINGERPRINTS_EXPORT std::unique_ptr<ExplicitBitVect>
FingerprintGenerator<std::uint32_t>::getFingerprint(
    const ROMol &mol, FingerprintFuncArguments &args) const;

template RDKIT_FINGERPRINTS_EXPORT std::unique_ptr<ExplicitBitVect>
FingerprintGenerator<std::uint64_t>::getFingerprint(
    const ROMol &mol, FingerprintFuncArguments &args) const;

template RDKIT_FINGERPRINTS_EXPORT std::vector<std::unique_ptr<ExplicitBitVect>>
FingerprintGenerator<std::uint32_t>::getFingerprints(
    const std::vector<const ROMol *> &mols, int numThreads) const;

template RDKIT_FINGERPRINTS_EXPORT std::vector<std::unique_ptr<ExplicitBitVect>>
FingerprintGenerator<std::uint64_t>::getFingerprints(
    const std::vector<const ROMol *> &mols, int numThreads) const;

template RDKIT_FINGERPRINTS_EXPORT std::vector<std::unique_ptr<SparseBitVect>>
FingerprintGenerator<std::uint32_t>::getSparseFingerprints(
    const std::vector<const ROMol *> &mols, int numThreads) const;

template RDKIT_FINGERPRINTS_EXPORT std::vector<std::unique_ptr<SparseBitVect>>
FingerprintGenerator<std::uint64_t>::getSparseFingerprints(
    const std::vector<const ROMol *> &mols, int numThreads) const;

template RDKIT_FINGERPRINTS_EXPORT
    std::vector<std::unique_ptr<SparseIntVect<std::uint32_t>>>
    FingerprintGenerator<std::uint32_t>::getCountFingerprints(
        const std::vector<const ROMol *> &mols, int numThreads) const;

template RDKIT_FINGERPRINTS_EXPORT
    std::vector<std::unique_ptr<SparseIntVect<std::uint32_t>>>
    FingerprintGenerator<std::uint64_t>::getCountFingerprints(
        const std::vector<const ROMol *> &mols, int numThreads) const;

template RDKIT_FINGERPRINTS_EXPORT
    std::vector<std::unique_ptr<SparseIntVect<std::uint32_t>>>
    FingerprintGenerator<std::uint32_t>::getSparseCountFingerprints(
        const std::vector<const ROMol *> &mols, int numThreads) const;

template RDKIT_FINGERPRINTS_EXPORT
    std::vector<std::unique_ptr<SparseIntVect<std::uint64_t>>>
    FingerprintGenerator<std::uint64_t>::getSparseCountFingerprints(
        const std::vector<const ROMol *> &mols, int numThreads) const;

SparseIntVect<std::uint64_t> *getSparseCountFP(const ROMol &mol,
                                               FPType fPType) {
  std::vector<const ROMol *> tempVect(1, &mol);
  return (*getSparseCountFPBulk(tempVect, fPType))[0];
}

SparseBitVect *getSparseFP(const ROMol &mol, FPType fPType) {
  std::vector<const ROMol *> tempVect(1, &mol);
  return (*getSparseFPBulk(tempVect, fPType))[0];
}

SparseIntVect<std::uint32_t> *getCountFP(const ROMol &mol, FPType fPType) {
  std::vector<const ROMol *> tempVect(1, &mol);
  return (*getCountFPBulk(tempVect, fPType))[0];
}

ExplicitBitVect *getFP(const ROMol &mol, FPType fPType) {
  std::vector<const ROMol *> tempVect(1, &mol);
  return (*getFPBulk(tempVect, fPType))[0];
}

std::vector<SparseIntVect<std::uint64_t> *> *getSparseCountFPBulk(
    const std::vector<const ROMol *> molVector, FPType fPType) {
  FingerprintGenerator<std::uint64_t> *generator = nullptr;
  switch (fPType) {
    case FPType::AtomPairFP: {
      generator = AtomPair::getAtomPairGenerator<std::uint64_t>();
      break;
    }
    case FPType::MorganFP: {
      generator = MorganFingerprint::getMorganGenerator<std::uint64_t>(2);
      break;
    }
    case FPType::RDKitFP: {
      generator = RDKitFP::getRDKitFPGenerator<std::uint64_t>();
      break;
    }
    case FPType::TopologicalTorsionFP: {
      generator =
          TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>();
      break;
    }
    default: {
      throw UnimplementedFPException(
          "Fingerprint type not implemented for getSparseCountFP");
    }
  }
  auto *res = new std::vector<SparseIntVect<std::uint64_t> *>();

  for (const auto *mol : molVector) {
    res->push_back(generator->getSparseCountFingerprint(*mol));
  }

  delete generator;
  return res;
}

std::vector<SparseBitVect *> *getSparseFPBulk(
    const std::vector<const ROMol *> molVector, FPType fPType) {
  FingerprintGenerator<std::uint64_t> *generator = nullptr;
  switch (fPType) {
    case FPType::AtomPairFP: {
      generator = AtomPair::getAtomPairGenerator<std::uint64_t>();
      break;
    }
    case FPType::MorganFP: {
      generator = MorganFingerprint::getMorganGenerator<std::uint64_t>(2);
      break;
    }
    case FPType::RDKitFP: {
      generator = RDKitFP::getRDKitFPGenerator<std::uint64_t>();
      break;
    }
    case FPType::TopologicalTorsionFP: {
      generator =
          TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>();
      break;
    }
    default: {
      throw UnimplementedFPException(
          "Fingerprint type not implemented for getSparseFP");
    }
  }
  auto *res = new std::vector<SparseBitVect *>();

  for (const auto *mol : molVector) {
    res->push_back(generator->getSparseFingerprint(*mol));
  }

  delete generator;
  return res;
}

std::vector<SparseIntVect<std::uint32_t> *> *getCountFPBulk(
    const std::vector<const ROMol *> molVector, FPType fPType) {
  FingerprintGenerator<std::uint64_t> *generator = nullptr;
  switch (fPType) {
    case FPType::AtomPairFP: {
      generator = AtomPair::getAtomPairGenerator<std::uint64_t>();
      break;
    }
    case FPType::MorganFP: {
      generator = MorganFingerprint::getMorganGenerator<std::uint64_t>(2);
      break;
    }
    case FPType::RDKitFP: {
      generator = RDKitFP::getRDKitFPGenerator<std::uint64_t>();
      break;
    }
    case FPType::TopologicalTorsionFP: {
      generator =
          TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>();
      break;
    }
    default: {
      throw UnimplementedFPException(
          "Fingerprint type not implemented for getCountFP");
    }
  }
  auto *res = new std::vector<SparseIntVect<std::uint32_t> *>();

  for (const auto *mol : molVector) {
    res->push_back(generator->getCountFingerprint(*mol));
  }

  delete generator;
  return res;
}

std::vector<ExplicitBitVect *> *getFPBulk(
    const std::vector<const ROMol *> molVector, FPType fPType) {
  FingerprintGenerator<std::uint64_t> *generator = nullptr;
  switch (fPType) {
    case FPType::AtomPairFP: {
      generator = AtomPair::getAtomPairGenerator<std::uint64_t>();
      break;
    }
    case FPType::MorganFP: {
      generator = MorganFingerprint::getMorganGenerator<std::uint64_t>(2);
      break;
    }
    case FPType::RDKitFP: {
      generator = RDKitFP::getRDKitFPGenerator<std::uint64_t>();
      break;
    }
    case FPType::TopologicalTorsionFP: {
      generator =
          TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>();
      break;
    }
    default: {
      throw UnimplementedFPException(
          "Fingerprint type not implemented for getFP");
    }
  }
  auto *res = new std::vector<ExplicitBitVect *>();

  for (const auto *mol : molVector) {
    res->push_back(generator->getFingerprint(*mol));
  }

  delete generator;
  return res;
}

}  // namespace RDKit
