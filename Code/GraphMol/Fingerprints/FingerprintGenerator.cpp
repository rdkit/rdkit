
#include <DataStructs/SparseIntVect.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/SparseBitVect.h>
#include "FingerprintGenerator.h"

namespace RDKit {

template <class T1, class T2>
FingerprintGenerator<T1, T2>::FingerprintGenerator(
    T1 atomEnvironmentGenerator, T2 fingerprintArguments,
    AtomInvariants *atomInvariants, BondInvariants *bondInvariants) {
  this->atomEnvironmentGenerator = atomEnvironmentGenerator;
  this->fingerprintArguments = fingerprintArguments;
  this->asCommonArguments = &this->fingerprintArguments;
  this->atomInvariants = atomInvariants;
  this->bondInvariants = bondInvariants;
}

template <class T1, class T2>
SparseIntVect<boost::uint32_t> *FingerprintGenerator<T1, T2>::getFingerprint(
    const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms,
    const std::vector<boost::uint32_t> *ignoreAtoms,
    const int confId,  // for atom pair fp
    const AdditionalOutput *additionalOutput) const {
  return 0;
}

template <class T1, class T2>
SparseBitVect *FingerprintGenerator<T1, T2>::getFingerprintAsBitVect(
    const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms,
    const std::vector<boost::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput) const {
  return 0;
}

template <class T1, class T2>
SparseIntVect<boost::uint32_t>
    *FingerprintGenerator<T1, T2>::getFoldedFingerprint(
        const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms,
        const std::vector<boost::uint32_t> *ignoreAtoms,
        const AdditionalOutput *additionalOutput) const {
  return 0;
}

template <class T1, class T2>
ExplicitBitVect *FingerprintGenerator<T1, T2>::getFoldedFingerprintAsBitVect(
    const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms,
    const std::vector<boost::uint32_t> *ignoreAtoms,
    const AdditionalOutput *additionalOutput) const {
  return 0;
}

}  // namespace RDKit