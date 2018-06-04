#ifndef RD_FINGERPRINTGEN_H_2018_05
#define RD_FINGERPRINTGEN_H_2018_05

#include <DataStructs/SparseIntVect.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/SparseBitVect.h>

namespace RDKit {
class ROMol;

union AdditionalOutput {
  // will hold differently stuctured additonal output

  std::map<boost::uint32_t,
           std::vector<std::pair<boost::uint32_t, boost::uint32_t>>>
      *bitInfoMap;
  // morgan fp
  // maps bitId -> vector of (atomId, radius)

  std::pair<std::vector<std::vector<boost::uint32_t>>,
            std::map<boost::uint32_t, std::vector<std::vector<int>>>> *bitInfo;
  // rdkit fp
  // first part, vector of bits set for each atom, must have the same size as
  // atom count for molecule
  // second part, maps bitId -> vector of paths

  std::vector<unsigned int> *atomCounts;
  // number of paths that set bits for each atom, must have the same size as
  // atom count for molecule
};

/*!
  /brief Abstract base class that holds arguments that are common amongst all
  fingerprint types and classes inherited from this would hold fingerprint type
  specific arguments

 */
class FingerprintArguments {
 public:
  FingerprintArguments(const bool countSimulation);
  const bool d_countSimulation;

  /*!
    /brief Returns the size of the fingerprint based on arguments

    /return boost::uint64_t size of the fingerprint
   */
  virtual boost::uint64_t getResultSize() const = 0;

  virtual ~FingerprintArguments() = 0;
};

/*!
  /brief abstract base class that holds atom-environments that will be hashed to
  generate the fingerprint

 */
class AtomEnvironment {
 public:
  /*!
    /brief calculates and returns the bit id to be set for this atom-environment

    /param arguments         Fingerprinting type specific arguments
    /param atomInvariants    Atom-invariants to be used during hashing
    /param bondInvariants    Bond-invariants to be used during hashing

    /return boost::uint32_t  calculated bit id for this environment
   */
  virtual boost::uint32_t getBitId(
      FingerprintArguments *arguments,
      const std::vector<boost::uint32_t> *atomInvariants = nullptr,
      const std::vector<boost::uint32_t> *bondInvariants = nullptr) const = 0;

  virtual ~AtomEnvironment() = 0;
};

/*!
  /brief abstract base class that generates atom-environments from a molecule

 */
class AtomEnvironmentGenerator {
 public:
  /*!
    /brief generate and return all atom-envorinments from a molecule

    /param mol               molecule to generate the atom-environments from
    /param arguments         fingerprint type specific arguments
    /param fromAtoms         atoms to be used during environment generation,
    usage of this parameter depends on the implementation of different
    fingerprint types
    /param ignoreAtoms      atoms to be ignored during environment generation,
    usage of this parameter depends on the implementation of different
    fingerprint types
    /param confId           which conformation to use during environment
    generation, needed for some fingerprint types
    /param additionalOutput contains pointers for additional outputs of
    fingerprinting operation, usage depends on implementation of the fingerprint
    type
    /param atomInvariants   atom invariants to be used during environment
    generation, in some cases some of the hashing can be done during environment
    generation so it is also passed here
    /param bondInvariants   bond invariants to be used during environment
    generation, same as atomInvariants it might be needed

    /return std::vector<AtomEnvironment *>  atom-environments generated from
    this molecule
   */
  virtual std::vector<AtomEnvironment *> getEnvironments(
      const ROMol &mol, FingerprintArguments *arguments,
      const std::vector<boost::uint32_t> *fromAtoms = nullptr,
      const std::vector<boost::uint32_t> *ignoreAtoms = nullptr,
      const int confId = -1, const AdditionalOutput *additionalOutput = nullptr,
      const std::vector<boost::uint32_t> *atomInvariants = nullptr,
      const std::vector<boost::uint32_t> *bondInvariants = nullptr) const = 0;

  virtual void cleanUpEnvironments(
      std::vector<AtomEnvironment *> atomEnvironments) const = 0;

  virtual ~AtomEnvironmentGenerator() = 0;
};

/*!
  /brief abstract base class for atom invariants generators

 */
class AtomInvariantsGenerator {
 public:
  /*!
    /brief get atom invariants from a molecule

    /param mol                   molecule to generate the atom invariants for
    /param fingerprintArguments  fingerprinting type specific arguments

    /return std::vector<boost::uint32_t> atom invariants generated for the given
    molecule
   */
  virtual std::vector<boost::uint32_t> getAtomInvariants(
      const ROMol &mol,
      const FingerprintArguments *fingerprintArguments) const = 0;
};

/*!
  /brief abstract base class for bond invariants generators

 */
class BondInvariantsGenerator {
 public:
  /*!
    /brief get bond invariants from a molecule

    /param mol                   molecule to generate the bond invariants for
    /param fingerprintArguments  fingerprinting type specific arguments

    /return std::vector<boost::uint32_t> bond invariants generated for the given
    molecule
   */
  virtual std::vector<boost::uint32_t> getBondInvariants(
      const ROMol &mol,
      const FingerprintArguments *fingerprintArguments) const = 0;
};

/*!
  /brief class that generates same fingerprint style for different output formats
  
 */
class FingerprintGenerator {
  FingerprintArguments *dp_fingerprintArguments;
  AtomEnvironmentGenerator *dp_atomEnvironmentGenerator;
  AtomInvariantsGenerator *dp_atomInvariantsGenerator;
  BondInvariantsGenerator *dp_bondInvariantsGenerator;

 public:
  FingerprintGenerator(
      AtomEnvironmentGenerator *atomEnvironmentGenerator,
      FingerprintArguments *fingerprintArguments,
      AtomInvariantsGenerator *atomInvariantsGenerator = nullptr,
      BondInvariantsGenerator *bondInvariantsGenerator = nullptr);


  void cleanUpResources();

  SparseIntVect<boost::uint32_t> *getFingerprint(
      const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms = nullptr,
      const std::vector<boost::uint32_t> *ignoreAtoms = nullptr,
      const int confId = -1,
      const AdditionalOutput *additionalOutput = nullptr) const;

  SparseBitVect *getFingerprintAsBitVect(
      const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms = nullptr,
      const std::vector<boost::uint32_t> *ignoreAtoms = nullptr,
      const int confId = -1,
      const AdditionalOutput *additionalOutput = nullptr) const;

  SparseIntVect<boost::uint32_t> *getFoldedFingerprint(
      const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms = nullptr,
      const std::vector<boost::uint32_t> *ignoreAtoms = nullptr,
      const AdditionalOutput *additionalOutput = nullptr) const;

  ExplicitBitVect *getFoldedFingerprintAsBitVect(
      const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms = nullptr,
      const std::vector<boost::uint32_t> *ignoreAtoms = nullptr,
      const AdditionalOutput *additionalOutput = nullptr) const;
};
}  // namespace RDKit

#endif
