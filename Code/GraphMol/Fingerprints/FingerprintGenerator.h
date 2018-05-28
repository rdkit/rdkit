#ifndef RD_FINGERPRINTGEN_H_2018_05
#define RD_FINGERPRINTGEN_H_2018_05

#include <DataStructs/SparseIntVect.h>
#include <DataStructs/ExplicitBitVect.h>

namespace RDKit {
class ROMol;

class AtomEnvironmentGenerator {};

union AdditionalOutput {
  // will hold differently stuctured additonal output
}

class AtomInvariants {
  // arguments

 public:
  std::vector<boost::uint32_t> getAtomInvariants(const ROMol &mol);
};

class BondInvariants {
  // arguments

 public:
  std::vector<boost::uint32_t> getBondInvariants(const ROMol &mol);
};

class FingerprintGenerator {
  AtomEnvironmentGenerator atomEnvironmentGenerator;
  AtomInvariants atomInvariant;
  BondInvariants bondInvariant;
  // more data to hold arguments and fp type

 public:
  SparseIntVect<boost::uint32_t> *getFingerprint(
      const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms = 0,
      const std::vector<boost::uint32_t> *ignoreAtoms = 0,
      const int confId = -1,  // for atom pair fp
      const AdditionalOutput *additionalOutput = 0) const;

  ExplicitBitVect *getFingerprintAsBitVect(
      const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms = 0,
      const std::vector<boost::uint32_t> *ignoreAtoms = 0,
      const int confId = -1,
      const AdditionalOutput *additionalOutput = 0) const;

  SparseIntVect<boost::uint32_t> *getHashedFingerprint(
      const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms = 0,
      const std::vector<boost::uint32_t> *ignoreAtoms = 0,
      const AdditionalOutput *additionalOutput = 0) const;

  ExplicitBitVect *getHashedFingerprintAsBitVect(
      const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms = 0,
      const std::vector<boost::uint32_t> *ignoreAtoms = 0,
      const AdditionalOutput *additionalOutput = 0) const;

  SparseIntVect<boost::uint32_t> *getCountFingerprint(
      const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms = 0,
      const std::vector<boost::uint32_t> *ignoreAtoms = 0,
      const AdditionalOutput *additionalOutput = 0) const;

  ExplicitBitVect *getCountFingerprintAsBitVect(
      const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms = 0,
      const std::vector<boost::uint32_t> *ignoreAtoms = 0,
      const AdditionalOutput *additionalOutput = 0) const;
}
}  // namespace RDKit

#endif
