//
//  Copyright (C) 2018 Boran Adas, Google Summer of Code
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

namespace RDKit {

template <typename OutputType>
FingerprintArguments<OutputType>::FingerprintArguments(
    const bool countSimulation, const std::vector<std::uint32_t> countBounds,
    std::uint32_t fpSize, std::uint32_t numBitsPerFeature)
    : d_countSimulation(countSimulation),
      d_countBounds(countBounds),
      d_fpSize(fpSize),
      d_numBitsPerFeature(numBitsPerFeature) {
  PRECONDITION(!countSimulation || !countBounds.empty(),
               "bad count bounds provided");
  PRECONDITION(d_numBitsPerFeature > 0, "numBitsPerFeature must be >0");
}

template FingerprintArguments<std::uint32_t>::FingerprintArguments(
    bool countSimulation, const std::vector<std::uint32_t> countBounds,
    std::uint32_t fpSize, std::uint32_t numBitsPerFeature);

template FingerprintArguments<std::uint64_t>::FingerprintArguments(
    bool countSimulation, const std::vector<std::uint32_t> countBounds,
    std::uint32_t fpSize, std::uint32_t numBitsPerFeature);

template <typename OutputType>
std::string FingerprintArguments<OutputType>::commonArgumentsString() const {
  return "Common arguments : countSimulation=" +
         std::to_string(d_countSimulation) +
         " fpSize=" + std::to_string(d_fpSize) +
         " bitsPerFeature=" + std::to_string(d_numBitsPerFeature);
}

template <typename OutputType>
FingerprintGenerator<OutputType>::FingerprintGenerator(
    AtomEnvironmentGenerator<OutputType> *atomEnvironmentGenerator,
    FingerprintArguments<OutputType> *fingerprintArguments,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    BondInvariantsGenerator *bondInvariantsGenerator, bool ownsAtomInvGenerator,
    bool ownsBondInvGenerator)
    : df_ownsAtomInvGenerator(ownsAtomInvGenerator),
      df_ownsBondInvGenerator(ownsBondInvGenerator) {
  this->dp_atomEnvironmentGenerator = atomEnvironmentGenerator;
  this->dp_fingerprintArguments = fingerprintArguments;
  this->dp_atomInvariantsGenerator = atomInvariantsGenerator;
  this->dp_bondInvariantsGenerator = bondInvariantsGenerator;
}

template FingerprintGenerator<std::uint32_t>::FingerprintGenerator(
    AtomEnvironmentGenerator<std::uint32_t> *atomEnvironmentGenerator,
    FingerprintArguments<std::uint32_t> *fingerprintArguments,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    BondInvariantsGenerator *bondInvariantsGenerator, bool ownsAtomInvGenerator,
    bool ownsBondInvGenerator);

template FingerprintGenerator<std::uint64_t>::FingerprintGenerator(
    AtomEnvironmentGenerator<std::uint64_t> *atomEnvironmentGenerator,
    FingerprintArguments<std::uint64_t> *fingerprintArguments,
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
SparseIntVect<OutputType>
    *FingerprintGenerator<OutputType>::getFingerprintHelper(
        const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
        const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
        const AdditionalOutput *additionalOutput,
        const std::vector<std::uint32_t> *customAtomInvariants,
        const std::vector<std::uint32_t> *customBondInvariants,
        const std::uint64_t fpSize) const {
  if (additionalOutput) {
    if (additionalOutput->atomCounts) {
      additionalOutput->atomCounts->resize(mol.getNumAtoms());
      std::fill(additionalOutput->atomCounts->begin(),
                additionalOutput->atomCounts->end(), 0);
    }
    if (additionalOutput->atomToBits) {
      additionalOutput->atomToBits->resize(mol.getNumAtoms());
      std::fill(additionalOutput->atomToBits->begin(),
                additionalOutput->atomToBits->end(),
                std::vector<std::uint64_t>());
    }
    if (additionalOutput->bitInfoMap) {
      additionalOutput->bitInfoMap->clear();
    }
    if (additionalOutput->bitPaths) {
      additionalOutput->bitPaths->clear();
    }
  }
  bool hashResults = false;
  if (fpSize != 0) {
    hashResults = true;
  }

  std::unique_ptr<std::vector<std::uint32_t>> atomInvariants = nullptr;
  if (customAtomInvariants) {
    atomInvariants.reset(new std::vector<std::uint32_t>(*customAtomInvariants));
  } else if (dp_atomInvariantsGenerator) {
    atomInvariants.reset(dp_atomInvariantsGenerator->getAtomInvariants(mol));
  }

  std::unique_ptr<std::vector<std::uint32_t>> bondInvariants = nullptr;
  if (customBondInvariants) {
    bondInvariants.reset(new std::vector<std::uint32_t>(*customBondInvariants));
  } else if (dp_bondInvariantsGenerator) {
    bondInvariants.reset(dp_bondInvariantsGenerator->getBondInvariants(mol));
  }

  // create all atom environments that will generate the bit-ids that will make
  // up the fingerprint
  auto atomEnvironments = dp_atomEnvironmentGenerator->getEnvironments(
      mol, dp_fingerprintArguments, fromAtoms, ignoreAtoms, confId,
      additionalOutput, atomInvariants.get(), bondInvariants.get(),
      hashResults);

  // allocate the result
  SparseIntVect<OutputType> *res = nullptr;

  if (fpSize != 0) {
    res = new SparseIntVect<OutputType>(fpSize);
  } else {
    res =
        new SparseIntVect<OutputType>(dp_fingerprintArguments->getResultSize());
  }

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
                                    additionalOutput, hashResults, fpSize);

    auto bitId = seed;
    if (fpSize != 0) {
      bitId %= fpSize;
    }
    res->setVal(bitId, res->getVal(bitId) + 1);
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
      }
    }
  }

  return res;
}

template <typename OutputType>
SparseIntVect<OutputType>
    *FingerprintGenerator<OutputType>::getSparseCountFingerprint(
        const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
        const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
        const AdditionalOutput *additionalOutput,
        const std::vector<std::uint32_t> *customAtomInvariants,
        const std::vector<std::uint32_t> *customBondInvariants) const {
  return getFingerprintHelper(mol, fromAtoms, ignoreAtoms, confId,
                              additionalOutput, customAtomInvariants,
                              customBondInvariants);
}

// todo getSparseFingerprint does not completely produce the same output as
// getSparseCountFingerprint. Count simulation and potential 64 bit outputs
// makes size limiting necessary for getSparseFingerprint. This can be
// changed if there is another way to avoid the size limitation of SparseBitVect
template <typename OutputType>
SparseBitVect *FingerprintGenerator<OutputType>::getSparseFingerprint(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput,
    const std::vector<std::uint32_t> *customAtomInvariants,
    const std::vector<std::uint32_t> *customBondInvariants) const {
  // make sure the result will fit into SparseBitVect
  std::uint32_t resultSize =
      std::min((std::uint64_t)std::numeric_limits<std::uint32_t>::max(),
               (std::uint64_t)dp_fingerprintArguments->getResultSize());

  std::uint32_t effectiveSize = resultSize;
  if (dp_fingerprintArguments->d_countSimulation) {
    // effective size needs to be smaller than result size to compansate for
    // count simulation
    effectiveSize /= dp_fingerprintArguments->d_countBounds.size();
  }

  SparseIntVect<OutputType> *tempResult = getFingerprintHelper(
      mol, fromAtoms, ignoreAtoms, confId, additionalOutput,
      customAtomInvariants, customBondInvariants, effectiveSize);

  auto *result = new SparseBitVect(resultSize);

  for (auto val : tempResult->getNonzeroElements()) {
    if (dp_fingerprintArguments->d_countSimulation) {
      for (unsigned int i = 0;
           i < dp_fingerprintArguments->d_countBounds.size(); ++i) {
        // for every bound in the d_countBounds in dp_fingerprintArguments, set
        // a bit if the occurrence count is equal or higher than the bound for
        // that bit
        const auto &bounds_count = dp_fingerprintArguments->d_countBounds;
        if (val.second >= static_cast<int>(bounds_count[i])) {
          result->setBit(val.first * bounds_count.size() + i);
        }
      }
    } else {
      result->setBit(val.first);
    }
  }

  delete tempResult;
  return result;
}

template <typename OutputType>
SparseIntVect<std::uint32_t>
    *FingerprintGenerator<OutputType>::getCountFingerprint(
        const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
        const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
        const AdditionalOutput *additionalOutput,
        const std::vector<std::uint32_t> *customAtomInvariants,
        const std::vector<std::uint32_t> *customBondInvariants) const {
  SparseIntVect<OutputType> *tempResult = getFingerprintHelper(
      mol, fromAtoms, ignoreAtoms, confId, additionalOutput,
      customAtomInvariants, customBondInvariants,
      dp_fingerprintArguments->d_fpSize);

  auto *result =
      new SparseIntVect<std::uint32_t>(dp_fingerprintArguments->d_fpSize);
  for (auto val : tempResult->getNonzeroElements()) {
    result->setVal(val.first, val.second);
  }

  delete tempResult;
  return result;
}

template <typename OutputType>
ExplicitBitVect *FingerprintGenerator<OutputType>::getFingerprint(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput,
    const std::vector<std::uint32_t> *customAtomInvariants,
    const std::vector<std::uint32_t> *customBondInvariants) const {
  std::uint32_t effectiveSize = dp_fingerprintArguments->d_fpSize;
  if (dp_fingerprintArguments->d_countSimulation) {
    // effective size needs to be smaller than result size to compansate for
    // count simulation
    effectiveSize /= dp_fingerprintArguments->d_countBounds.size();
  }
  SparseIntVect<OutputType> *tempResult = getFingerprintHelper(
      mol, fromAtoms, ignoreAtoms, confId, additionalOutput,
      customAtomInvariants, customBondInvariants, effectiveSize);

  auto *result = new ExplicitBitVect(dp_fingerprintArguments->d_fpSize);
  for (auto val : tempResult->getNonzeroElements()) {
    if (dp_fingerprintArguments->d_countSimulation) {
      for (unsigned int i = 0;
           i < dp_fingerprintArguments->d_countBounds.size(); ++i) {
        // for every bound in the d_countBounds in dp_fingerprintArguments, set
        // a bit if the occurrence count is equal or higher than the bound for
        // that bit
        const auto &bounds_count = dp_fingerprintArguments->d_countBounds;
        if (val.second >= static_cast<int>(bounds_count[i])) {
          result->setBit(val.first * bounds_count.size() + i);
        }
      }
    } else {
      result->setBit(val.first);
    }
  }

  delete tempResult;
  return result;
}

template RDKIT_FINGERPRINTS_EXPORT SparseIntVect<std::uint32_t>
    *FingerprintGenerator<std::uint32_t>::getSparseCountFingerprint(
        const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
        const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
        const AdditionalOutput *additionalOutput,
        const std::vector<std::uint32_t> *customAtomInvariants,
        const std::vector<std::uint32_t> *customBondInvariants) const;

template RDKIT_FINGERPRINTS_EXPORT SparseIntVect<std::uint64_t>
    *FingerprintGenerator<std::uint64_t>::getSparseCountFingerprint(
        const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
        const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
        const AdditionalOutput *additionalOutput,
        const std::vector<std::uint32_t> *customAtomInvariants,
        const std::vector<std::uint32_t> *customBondInvariants) const;

template RDKIT_FINGERPRINTS_EXPORT SparseBitVect *
FingerprintGenerator<std::uint32_t>::getSparseFingerprint(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput,
    const std::vector<std::uint32_t> *customAtomInvariants,
    const std::vector<std::uint32_t> *customBondInvariants) const;

template RDKIT_FINGERPRINTS_EXPORT SparseBitVect *
FingerprintGenerator<std::uint64_t>::getSparseFingerprint(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput,
    const std::vector<std::uint32_t> *customAtomInvariants,
    const std::vector<std::uint32_t> *customBondInvariants) const;

template RDKIT_FINGERPRINTS_EXPORT SparseIntVect<std::uint32_t>
    *FingerprintGenerator<std::uint32_t>::getCountFingerprint(
        const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
        const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
        const AdditionalOutput *additionalOutput,
        const std::vector<std::uint32_t> *customAtomInvariants,
        const std::vector<std::uint32_t> *customBondInvariants) const;

template RDKIT_FINGERPRINTS_EXPORT SparseIntVect<std::uint32_t>
    *FingerprintGenerator<std::uint64_t>::getCountFingerprint(
        const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
        const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
        const AdditionalOutput *additionalOutput,
        const std::vector<std::uint32_t> *customAtomInvariants,
        const std::vector<std::uint32_t> *customBondInvariants) const;

template RDKIT_FINGERPRINTS_EXPORT ExplicitBitVect *
FingerprintGenerator<std::uint32_t>::getFingerprint(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput,
    const std::vector<std::uint32_t> *customAtomInvariants,
    const std::vector<std::uint32_t> *customBondInvariants) const;

template RDKIT_FINGERPRINTS_EXPORT ExplicitBitVect *
FingerprintGenerator<std::uint64_t>::getFingerprint(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput,
    const std::vector<std::uint32_t> *customAtomInvariants,
    const std::vector<std::uint32_t> *customBondInvariants) const;

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
