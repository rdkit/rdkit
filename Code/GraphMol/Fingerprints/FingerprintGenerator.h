//
//  Copyright (C) 2018 Boran Adas, Google Summer of Code
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
#include <cstdint>

namespace RDKit {
class ROMol;

struct RDKIT_FINGERPRINTS_EXPORT AdditionalOutput {
  // will review this structure once more fingerprint types are implemented

  std::vector<std::vector<std::uint64_t>> *atomToBits;

  std::map<std::uint32_t, std::vector<std::pair<std::uint32_t, std::uint32_t>>>
      *bitInfoMap;
  // morgan fp
  // maps bitId -> vector of (atomId, radius)

  std::pair<std::vector<std::vector<std::uint32_t>>,
            std::map<std::uint32_t, std::vector<std::vector<int>>>> *bitInfo;
  // rdkit fp
  // first part, vector of bits set for each atom, must have the same size as
  // atom count for molecule
  // second part, maps bitId -> vector of paths

  std::vector<unsigned int> *atomCounts;
  // number of paths that set bits for each atom, must have the same size as
  // atom count for molecule
};

/*!
  \brief Abstract base class that holds molecule independent arguments that are
  common amongst all fingerprint types and classes inherited from this would
  hold fingerprint type specific arguments

 */
template <typename OutputType>
class RDKIT_FINGERPRINTS_EXPORT FingerprintArguments
    : private boost::noncopyable {
 public:
  FingerprintArguments(bool countSimulation,
                       const std::vector<std::uint32_t> countBounds,
                       std::uint32_t fpSize,
                       std::uint32_t numBitsPerFeature = 1);
  const bool d_countSimulation;
  const std::vector<std::uint32_t> d_countBounds;
  const std::uint32_t d_fpSize;
  const std::uint32_t d_numBitsPerFeature;

  /*!
    \brief Returns the size of the fingerprint based on arguments

    \return OutputType size of the fingerprint
   */
  virtual OutputType getResultSize() const = 0;

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

  virtual ~FingerprintArguments(){};
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
  virtual OutputType getBitId(FingerprintArguments<OutputType> *arguments,
                              const std::vector<std::uint32_t> *atomInvariants,
                              const std::vector<std::uint32_t> *bondInvariants,
                              const AdditionalOutput *AdditionalOutput,
                              const bool hashResults = false) const = 0;

  virtual ~AtomEnvironment(){};
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
      const ROMol &mol, FingerprintArguments<OutputType> *arguments,
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

  virtual ~AtomEnvironmentGenerator(){};
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

  virtual ~AtomInvariantsGenerator(){};
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

  virtual ~BondInvariantsGenerator(){};
  virtual BondInvariantsGenerator *clone() const = 0;
};  // namespace RDKit

/*!
  \brief class that generates same fingerprint style for different output
  formats

 */
template <typename OutputType>
class RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator
    : private boost::noncopyable {
  FingerprintArguments<OutputType> *dp_fingerprintArguments;
  AtomEnvironmentGenerator<OutputType> *dp_atomEnvironmentGenerator;
  AtomInvariantsGenerator *dp_atomInvariantsGenerator;
  BondInvariantsGenerator *dp_bondInvariantsGenerator;
  const bool df_ownsAtomInvGenerator;
  const bool df_ownsBondInvGenerator;

  SparseIntVect<OutputType> *getFingerprintHelper(
      const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms = nullptr,
      const std::vector<std::uint32_t> *ignoreAtoms = nullptr,
      const int confId = -1, const AdditionalOutput *additionalOutput = nullptr,
      const std::vector<std::uint32_t> *customAtomInvariants = nullptr,
      const std::vector<std::uint32_t> *customBondInvariants = nullptr,
      const std::uint64_t fpSize = 0) const;

 public:
  FingerprintGenerator(
      AtomEnvironmentGenerator<OutputType> *atomEnvironmentGenerator,
      FingerprintArguments<OutputType> *fingerprintArguments,
      AtomInvariantsGenerator *atomInvariantsGenerator = nullptr,
      BondInvariantsGenerator *bondInvariantsGenerator = nullptr,
      bool ownsAtomInvGenerator = false, bool ownsBondInvGenerator = false);

  ~FingerprintGenerator();

  SparseIntVect<OutputType> *getSparseCountFingerprint(
      const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms = nullptr,
      const std::vector<std::uint32_t> *ignoreAtoms = nullptr,
      const int confId = -1, const AdditionalOutput *additionalOutput = nullptr,
      const std::vector<std::uint32_t> *customAtomInvariants = nullptr,
      const std::vector<std::uint32_t> *customBondInvariants = nullptr) const;

  SparseBitVect *getSparseFingerprint(
      const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms = nullptr,
      const std::vector<std::uint32_t> *ignoreAtoms = nullptr,
      const int confId = -1, const AdditionalOutput *additionalOutput = nullptr,
      const std::vector<std::uint32_t> *customAtomInvariants = nullptr,
      const std::vector<std::uint32_t> *customBondInvariants = nullptr) const;

  SparseIntVect<std::uint32_t> *getCountFingerprint(
      const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms = nullptr,
      const std::vector<std::uint32_t> *ignoreAtoms = nullptr,
      const int confId = -1, const AdditionalOutput *additionalOutput = nullptr,
      const std::vector<std::uint32_t> *customAtomInvariants = nullptr,
      const std::vector<std::uint32_t> *customBondInvariants = nullptr) const;

  ExplicitBitVect *getFingerprint(
      const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms = nullptr,
      const std::vector<std::uint32_t> *ignoreAtoms = nullptr,
      const int confId = -1, const AdditionalOutput *additionalOutput = nullptr,
      const std::vector<std::uint32_t> *customAtomInvariants = nullptr,
      const std::vector<std::uint32_t> *customBondInvariants = nullptr) const;

  std::string infoString() const;
};

enum class FPType { AtomPairFP, MorganFP, RDKitFP, TopologicalTorsionFP };

//! used to indicate errors for unimplemented fp types in convenience functions
class RDKIT_FINGERPRINTS_EXPORT UnimplementedFPException
    : public std::exception {
 public:
  //! construct with an error message
  UnimplementedFPException(const char *msg) : _msg(msg){};
  //! construct with an error message
  UnimplementedFPException(const std::string &msg) : _msg(msg){};
  //! get the error message
  const char *what() const noexcept override { return _msg.c_str(); };
  ~UnimplementedFPException() noexcept {};

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

RDKIT_FINGERPRINTS_EXPORT std::vector<SparseIntVect<std::uint32_t> *>
    *getCountFPBulk(const std::vector<const ROMol *> molVector, FPType fPType);

RDKIT_FINGERPRINTS_EXPORT std::vector<ExplicitBitVect *> *getFPBulk(
    const std::vector<const ROMol *> molVector, FPType fPType);

}  // namespace RDKit

#endif
