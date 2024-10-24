//
//  Copyright (C) 2018-2022 Boran Adas and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef RD_FINGERPRINTGEN_H_2018_05
#define RD_FINGERPRINTGEN_H_2018_05

#include <DataStructs/SparseIntVect.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/SparseBitVect.h>
#include <utility>
#include <vector>
#include <memory>
#include <cstdint>

namespace RDKit {
class ROMol;

struct RDKIT_FINGERPRINTS_EXPORT AdditionalOutput {
  using atomToBitsType = std::vector<std::vector<std::uint64_t>>;
  using bitInfoMapType =
      std::map<std::uint64_t,
               std::vector<std::pair<std::uint32_t, std::uint32_t>>>;
  using bitPathsType = std::map<std::uint64_t, std::vector<std::vector<int>>>;
  using atomCountsType = std::vector<unsigned int>;

  // numAtoms long
  atomToBitsType *atomToBits = nullptr;

  // bitId -> vector of (atomId, radius) for morgan
  // bitId -> (atom1, atom2) for atom pairs
  bitInfoMapType *bitInfoMap = nullptr;

  // rdkit fp
  // maps bitId -> vector of bond paths
  bitPathsType *bitPaths = nullptr;

  // number of paths that set bits for each atom, must have the same size as
  // atom count for molecule.
  atomCountsType *atomCounts = nullptr;

  void allocateAtomToBits() {
    atomToBitsHolder.reset(new atomToBitsType);
    atomToBits = atomToBitsHolder.get();
  }
  void allocateBitInfoMap() {
    bitInfoMapHolder.reset(new bitInfoMapType);
    bitInfoMap = bitInfoMapHolder.get();
  }
  void allocateBitPaths() {
    bitPathsHolder.reset(new bitPathsType);
    bitPaths = bitPathsHolder.get();
  }
  void allocateAtomCounts() {
    atomCountsHolder.reset(new atomCountsType);
    atomCounts = atomCountsHolder.get();
  }

 private:
  std::unique_ptr<atomToBitsType> atomToBitsHolder;
  std::unique_ptr<bitInfoMapType> bitInfoMapHolder;
  std::unique_ptr<bitPathsType> bitPathsHolder;
  std::unique_ptr<atomCountsType> atomCountsHolder;
};

/*!
  \brief Abstract base class that holds molecule independent arguments that are
  common amongst all fingerprint types and classes inherited from this would
  hold fingerprint type specific arguments

 */
class RDKIT_FINGERPRINTS_EXPORT FingerprintArguments {
 public:
  FingerprintArguments(bool countSimulation,
                       const std::vector<std::uint32_t> countBounds,
                       std::uint32_t fpSize,
                       std::uint32_t numBitsPerFeature = 1,
                       bool includeChirality = false);
  bool df_countSimulation = false;
  bool df_includeChirality = false;
  std::vector<std::uint32_t> d_countBounds;
  std::uint32_t d_fpSize = 2048;
  std::uint32_t d_numBitsPerFeature = 1;

  /**
   \brief method that returns information string about the fingerprint specific
   argument set and the arguments themselves

   \return std::string information string
   */
  virtual std::string infoString() const = 0;

  /**
   \brief method that returns information string about common fingerprinting
   arguments' values

   \return std::string information string
   */
  std::string commonArgumentsString() const;

  virtual ~FingerprintArguments() {}
  FingerprintArguments() = default;
};

/*!
  \brief abstract base class that holds atom-environments that will be hashed to
  generate the fingerprint

 */
template <typename OutputType>
class RDKIT_FINGERPRINTS_EXPORT AtomEnvironment : private boost::noncopyable {
 public:
  /*!
    \brief calculates and returns the bit id to be set for this atom-environment

    \param arguments         Fingerprinting type specific molecule independent
    arguments
    \param atomInvariants    Atom-invariants to be used during hashing
    \param bondInvariants    Bond-invariants to be used during hashing
    \param hashResults   if set results will be ready to be modded

    \return OutputType  calculated bit id for this environment
   */
  virtual OutputType getBitId(FingerprintArguments *arguments,
                              const std::vector<std::uint32_t> *atomInvariants,
                              const std::vector<std::uint32_t> *bondInvariants,
                              AdditionalOutput *AdditionalOutput,
                              const bool hashResults = false,
                              const std::uint64_t fpSize = 0) const = 0;
  virtual void updateAdditionalOutput(AdditionalOutput *AdditionalOutput,
                                      size_t bitId) const = 0;

  virtual ~AtomEnvironment() {}
};

/*!
  \brief abstract base class that generates atom-environments from a molecule

 */
template <typename OutputType>
class RDKIT_FINGERPRINTS_EXPORT AtomEnvironmentGenerator
    : private boost::noncopyable {
 public:
  /*!
    \brief generate and return all atom-envorinments from a molecule

    \param mol               molecule to generate the atom-environments from
    \param arguments         fingerprint type specific molecule independent
    arguments
    \param fromAtoms         atoms to be used during environment generation,
    usage of this parameter depends on the implementation of different
    fingerprint types
    \param ignoreAtoms      atoms to be ignored during environment generation,
    usage of this parameter depends on the implementation of different
    fingerprint types
    \param confId           which conformation to use during environment
    generation, needed for some fingerprint types
    \param additionalOutput contains pointers for additional outputs of
    fingerprinting operation, usage depends on implementation of the fingerprint
    type
    \param atomInvariants   atom invariants to be used during environment
    generation, in some cases some of the hashing can be done during environment
    generation so it is also passed here
    \param bondInvariants   bond invariants to be used during environment
    generation, same as atomInvariants it might be needed
    \param hashResults   if set results will be ready to be modded

    \return std::vector<AtomEnvironment *>  atom-environments generated from
    this molecule
   */
  virtual std::vector<AtomEnvironment<OutputType> *> getEnvironments(
      const ROMol &mol, FingerprintArguments *arguments,
      const std::vector<std::uint32_t> *fromAtoms = nullptr,
      const std::vector<std::uint32_t> *ignoreAtoms = nullptr,
      const int confId = -1, const AdditionalOutput *additionalOutput = nullptr,
      const std::vector<std::uint32_t> *atomInvariants = nullptr,
      const std::vector<std::uint32_t> *bondInvariants = nullptr,
      const bool hashResults = false) const = 0;

  /**
   \brief method that returns information about this /c AtomEnvironmentGenerator
   and its arguments if any

   \return std::string information string
   */
  virtual std::string infoString() const = 0;
  /*!
    \brief Returns the size of the fingerprint based on arguments

    \return OutputType size of the fingerprint
   */
  virtual OutputType getResultSize() const = 0;

  const FingerprintArguments *dp_fingerprintArguments;

  virtual ~AtomEnvironmentGenerator() {}
};

/*!
  \brief abstract base class for atom invariants generators

 */
class RDKIT_FINGERPRINTS_EXPORT AtomInvariantsGenerator
    : private boost::noncopyable {
 public:
  /*!
    \brief get atom invariants from a molecule

    \param mol                   molecule to generate the atom invariants for

    \return std::vector<std::uint32_t> atom invariants generated for the given
    molecule
   */
  virtual std::vector<std::uint32_t> *getAtomInvariants(
      const ROMol &mol) const = 0;

  /**
   \brief method that returns information about this /c AtomInvariantsGenerator
   and its arguments

   \return std::string information string
   */
  virtual std::string infoString() const = 0;

  virtual ~AtomInvariantsGenerator() {}
  virtual AtomInvariantsGenerator *clone() const = 0;
};

/*!
  \brief abstract base class for bond invariants generators

 */
class RDKIT_FINGERPRINTS_EXPORT BondInvariantsGenerator
    : private boost::noncopyable {
 public:
  /*!
    \brief get bond invariants from a molecule

    \param mol                   molecule to generate the bond invariants for

    \return std::vector<std::uint32_t> bond invariants generated for the given
    molecule
   */
  virtual std::vector<std::uint32_t> *getBondInvariants(
      const ROMol &mol) const = 0;

  /**
 \brief method that returns information about this /c BondInvariantsGenerator
 and its arguments

 \return std::string information string
 */
  virtual std::string infoString() const = 0;

  virtual ~BondInvariantsGenerator() {}
  virtual BondInvariantsGenerator *clone() const = 0;
};  // namespace RDKit

/*!
  \brief struct that makes calling the fingerprint generation functions easier

  FingerprintFuncArguments doesn't own any of the pointers it stores.

 */
struct FingerprintFuncArguments {
  const std::vector<std::uint32_t> *fromAtoms = nullptr;
  const std::vector<std::uint32_t> *ignoreAtoms = nullptr;
  int confId = -1;
  AdditionalOutput *additionalOutput = nullptr;
  const std::vector<std::uint32_t> *customAtomInvariants = nullptr;
  const std::vector<std::uint32_t> *customBondInvariants = nullptr;
  FingerprintFuncArguments() = default;
  FingerprintFuncArguments(
      const std::vector<std::uint32_t> *fromAtoms_arg,
      const std::vector<std::uint32_t> *ignoreAtoms_arg, int confId_arg,
      AdditionalOutput *additionalOutput_arg,
      const std::vector<std::uint32_t> *customAtomInvariants_arg,
      const std::vector<std::uint32_t> *customBondInvariants_arg)
      : fromAtoms(fromAtoms_arg),
        ignoreAtoms(ignoreAtoms_arg),
        confId(confId_arg),
        additionalOutput(additionalOutput_arg),
        customAtomInvariants(customAtomInvariants_arg),
        customBondInvariants(customBondInvariants_arg) {};
};

/*!
  \brief class that generates same fingerprint style for different output
  formats

 */
template <typename OutputType>
class RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator
    : private boost::noncopyable {
  FingerprintArguments *dp_fingerprintArguments;
  AtomEnvironmentGenerator<OutputType> *dp_atomEnvironmentGenerator;
  AtomInvariantsGenerator *dp_atomInvariantsGenerator;
  BondInvariantsGenerator *dp_bondInvariantsGenerator;
  const bool df_ownsAtomInvGenerator;
  const bool df_ownsBondInvGenerator;

  std::unique_ptr<SparseIntVect<OutputType>> getFingerprintHelper(
      const ROMol &mol, FingerprintFuncArguments &args,
      const std::uint64_t fpSize = 0) const;

 public:
  FingerprintGenerator(
      AtomEnvironmentGenerator<OutputType> *atomEnvironmentGenerator,
      FingerprintArguments *fingerprintArguments,
      AtomInvariantsGenerator *atomInvariantsGenerator = nullptr,
      BondInvariantsGenerator *bondInvariantsGenerator = nullptr,
      bool ownsAtomInvGenerator = false, bool ownsBondInvGenerator = false);

  ~FingerprintGenerator();

  FingerprintArguments *getOptions() { return dp_fingerprintArguments; };
  const FingerprintArguments *getOptions() const {
    return dp_fingerprintArguments;
  };

  std::unique_ptr<SparseIntVect<OutputType>> getSparseCountFingerprint(
      const ROMol &mol, FingerprintFuncArguments &args) const;

  std::unique_ptr<SparseBitVect> getSparseFingerprint(
      const ROMol &mol, FingerprintFuncArguments &args) const;

  std::unique_ptr<SparseIntVect<std::uint32_t>> getCountFingerprint(
      const ROMol &mol, FingerprintFuncArguments &args) const;

  std::unique_ptr<ExplicitBitVect> getFingerprint(
      const ROMol &mol, FingerprintFuncArguments &args) const;

  std::vector<std::unique_ptr<ExplicitBitVect>> getFingerprints(
      const std::vector<const ROMol *> &mols, int numThreads = 1) const;

  std::vector<std::unique_ptr<SparseBitVect>> getSparseFingerprints(
      const std::vector<const ROMol *> &mols, int numThreads = 1) const;

  std::vector<std::unique_ptr<SparseIntVect<std::uint32_t>>>
  getCountFingerprints(const std::vector<const ROMol *> &mols,
                       int numThreads = 1) const;

  std::vector<std::unique_ptr<SparseIntVect<OutputType>>>
  getSparseCountFingerprints(const std::vector<const ROMol *> &mols,
                             int numThreads = 1) const;

  SparseIntVect<OutputType> *getSparseCountFingerprint(
      const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms = nullptr,
      const std::vector<std::uint32_t> *ignoreAtoms = nullptr, int confId = -1,
      AdditionalOutput *additionalOutput = nullptr,
      const std::vector<std::uint32_t> *customAtomInvariants = nullptr,
      const std::vector<std::uint32_t> *customBondInvariants = nullptr) const {
    FingerprintFuncArguments ffa(fromAtoms, ignoreAtoms, confId,
                                 additionalOutput, customAtomInvariants,
                                 customBondInvariants);
    return getSparseCountFingerprint(mol, ffa).release();
  };

  SparseBitVect *getSparseFingerprint(
      const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms = nullptr,
      const std::vector<std::uint32_t> *ignoreAtoms = nullptr, int confId = -1,
      AdditionalOutput *additionalOutput = nullptr,
      const std::vector<std::uint32_t> *customAtomInvariants = nullptr,
      const std::vector<std::uint32_t> *customBondInvariants = nullptr) const {
    FingerprintFuncArguments ffa(fromAtoms, ignoreAtoms, confId,
                                 additionalOutput, customAtomInvariants,
                                 customBondInvariants);
    return getSparseFingerprint(mol, ffa).release();
  };

  SparseIntVect<std::uint32_t> *getCountFingerprint(
      const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms = nullptr,
      const std::vector<std::uint32_t> *ignoreAtoms = nullptr, int confId = -1,
      AdditionalOutput *additionalOutput = nullptr,
      const std::vector<std::uint32_t> *customAtomInvariants = nullptr,
      const std::vector<std::uint32_t> *customBondInvariants = nullptr) const {
    FingerprintFuncArguments ffa(fromAtoms, ignoreAtoms, confId,
                                 additionalOutput, customAtomInvariants,
                                 customBondInvariants);
    return getCountFingerprint(mol, ffa).release();
  };

  ExplicitBitVect *getFingerprint(
      const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms = nullptr,
      const std::vector<std::uint32_t> *ignoreAtoms = nullptr, int confId = -1,
      AdditionalOutput *additionalOutput = nullptr,
      const std::vector<std::uint32_t> *customAtomInvariants = nullptr,
      const std::vector<std::uint32_t> *customBondInvariants = nullptr) const {
    FingerprintFuncArguments ffa(fromAtoms, ignoreAtoms, confId,
                                 additionalOutput, customAtomInvariants,
                                 customBondInvariants);
    return getFingerprint(mol, ffa).release();
  };

  std::string infoString() const;
};

template RDKIT_FINGERPRINTS_EXPORT ExplicitBitVect *
FingerprintGenerator<std::uint32_t>::getFingerprint(
    const ROMol &, const std::vector<std::uint32_t> *,
    const std::vector<std::uint32_t> *, int, AdditionalOutput *,
    const std::vector<std::uint32_t> *,
    const std::vector<std::uint32_t> *) const;
template RDKIT_FINGERPRINTS_EXPORT SparseBitVect *
FingerprintGenerator<std::uint32_t>::getSparseFingerprint(
    const ROMol &, const std::vector<std::uint32_t> *,
    const std::vector<std::uint32_t> *, int, AdditionalOutput *,
    const std::vector<std::uint32_t> *,
    const std::vector<std::uint32_t> *) const;
template RDKIT_FINGERPRINTS_EXPORT SparseIntVect<std::uint32_t> *
FingerprintGenerator<std::uint32_t>::getCountFingerprint(
    const ROMol &, const std::vector<std::uint32_t> *,
    const std::vector<std::uint32_t> *, int, AdditionalOutput *,
    const std::vector<std::uint32_t> *,
    const std::vector<std::uint32_t> *) const;
template RDKIT_FINGERPRINTS_EXPORT SparseIntVect<std::uint32_t> *
FingerprintGenerator<std::uint32_t>::getSparseCountFingerprint(
    const ROMol &, const std::vector<std::uint32_t> *,
    const std::vector<std::uint32_t> *, int, AdditionalOutput *,
    const std::vector<std::uint32_t> *,
    const std::vector<std::uint32_t> *) const;
template RDKIT_FINGERPRINTS_EXPORT ExplicitBitVect *
FingerprintGenerator<std::uint64_t>::getFingerprint(
    const ROMol &, const std::vector<std::uint32_t> *,
    const std::vector<std::uint32_t> *, int, AdditionalOutput *,
    const std::vector<std::uint32_t> *,
    const std::vector<std::uint32_t> *) const;
template RDKIT_FINGERPRINTS_EXPORT SparseBitVect *
FingerprintGenerator<std::uint64_t>::getSparseFingerprint(
    const ROMol &, const std::vector<std::uint32_t> *,
    const std::vector<std::uint32_t> *, int, AdditionalOutput *,
    const std::vector<std::uint32_t> *,
    const std::vector<std::uint32_t> *) const;
template RDKIT_FINGERPRINTS_EXPORT SparseIntVect<std::uint32_t> *
FingerprintGenerator<std::uint64_t>::getCountFingerprint(
    const ROMol &, const std::vector<std::uint32_t> *,
    const std::vector<std::uint32_t> *, int, AdditionalOutput *,
    const std::vector<std::uint32_t> *,
    const std::vector<std::uint32_t> *) const;
template RDKIT_FINGERPRINTS_EXPORT SparseIntVect<std::uint64_t> *
FingerprintGenerator<std::uint64_t>::getSparseCountFingerprint(
    const ROMol &, const std::vector<std::uint32_t> *,
    const std::vector<std::uint32_t> *, int, AdditionalOutput *,
    const std::vector<std::uint32_t> *,
    const std::vector<std::uint32_t> *) const;

enum class FPType { AtomPairFP, MorganFP, RDKitFP, TopologicalTorsionFP };

//! used to indicate errors for unimplemented fp types in convenience
//! functions
class RDKIT_FINGERPRINTS_EXPORT UnimplementedFPException
    : public std::exception {
 public:
  //! construct with an error message
  UnimplementedFPException(const char *msg) : _msg(msg) {}
  //! construct with an error message
  UnimplementedFPException(std::string msg) : _msg(std::move(msg)) {}
  //! get the error message
  const char *what() const noexcept override { return _msg.c_str(); }
  ~UnimplementedFPException() noexcept override = default;

 private:
  std::string _msg;
};

// convenience functions, fingerprint generation with default values

RDKIT_FINGERPRINTS_EXPORT SparseIntVect<std::uint64_t> *getSparseCountFP(
    const ROMol &mol, FPType fPType);

RDKIT_FINGERPRINTS_EXPORT SparseBitVect *getSparseFP(const ROMol &mol,
                                                     FPType fPType);

RDKIT_FINGERPRINTS_EXPORT SparseIntVect<std::uint32_t> *getCountFP(
    const ROMol &mol, FPType fPType);

RDKIT_FINGERPRINTS_EXPORT ExplicitBitVect *getFP(const ROMol &mol,
                                                 FPType fPType);

RDKIT_FINGERPRINTS_EXPORT std::vector<SparseIntVect<std::uint64_t> *> *
getSparseCountFPBulk(const std::vector<const ROMol *> molVector, FPType fPType);

RDKIT_FINGERPRINTS_EXPORT std::vector<SparseBitVect *> *getSparseFPBulk(
    const std::vector<const ROMol *> molVector, FPType fPType);

RDKIT_FINGERPRINTS_EXPORT std::vector<SparseIntVect<std::uint32_t> *> *
getCountFPBulk(const std::vector<const ROMol *> molVector, FPType fPType);

RDKIT_FINGERPRINTS_EXPORT std::vector<ExplicitBitVect *> *getFPBulk(
    const std::vector<const ROMol *> molVector, FPType fPType);

}  // namespace RDKit

#endif
