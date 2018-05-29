#ifndef RD_FINGERPRINTGEN_H_2018_05
#define RD_FINGERPRINTGEN_H_2018_05

#include <DataStructs/SparseIntVect.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/SparseBitVect.h>

namespace RDKit {
class ROMol;

class AtomEnvironmentGenerator {};

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

  SparseBitVect *getFingerprintAsBitVect(
      const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms = 0,
      const std::vector<boost::uint32_t> *ignoreAtoms = 0,
      const int confId = -1,
      const AdditionalOutput *additionalOutput = 0) const;

  SparseIntVect<boost::uint32_t> *getFoldedFingerprint(
      const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms = 0,
      const std::vector<boost::uint32_t> *ignoreAtoms = 0,
      const AdditionalOutput *additionalOutput = 0) const;

  ExplicitBitVect *getFoldedFingerprintAsBitVect(
      const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms = 0,
      const std::vector<boost::uint32_t> *ignoreAtoms = 0,
      const AdditionalOutput *additionalOutput = 0) const;
}
}  // namespace RDKit

#endif
