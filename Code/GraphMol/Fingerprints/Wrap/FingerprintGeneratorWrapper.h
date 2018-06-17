#ifndef RD_FINGERPRINTWRAP_H_2018_06
#define RD_FINGERPRINTWRAP_H_2018_06

#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <boost/python.hpp>

namespace python = boost::python;

namespace RDKit {
namespace FingerprintWrapper {

class FingerprintGeneratorWrapper {
 public:
  FingerprintGenerator *dp_fingerprintGenerator;

  SparseIntVect<std::uint32_t> *getFingerprint(const ROMol &mol,
                                               python::object py_fromAtoms,
                                               python::object py_ignoreAtoms,
                                               const int confId) const;

  SparseBitVect *getFingerprintAsBitVect(const ROMol &mol,
                                         python::object py_fromAtoms,
                                         python::object py_ignoreAtoms,
                                         const int confId) const;

  SparseIntVect<std::uint32_t> *getFoldedFingerprint(
      const ROMol &mol, python::object py_fromAtoms,
      python::object py_ignoreAtoms, const int confId) const;

  ExplicitBitVect *getFoldedFingerprintAsBitVect(const ROMol &mol,
                                                 python::object py_fromAtoms,
                                                 python::object py_ignoreAtoms,
                                                 const int confId) const;

  FingerprintGeneratorWrapper();
  ~FingerprintGeneratorWrapper();
};

}  // namespace FingerprintWrapper
}  // namespace RDKit

#endif
