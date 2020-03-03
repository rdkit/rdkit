/*! \file MACCS.h

*/
#include <RDGeneral/export.h>
#ifndef __RD_MACCSFPS_H__
#define __RD_MACCSFPS_H__
#include <string>

class ExplicitBitVect;
namespace RDKit {
class ROMol;
namespace MACCSFingerprints {
const std::string maccsFingerprintVersion = "2.0.0";

//! returns the MACCS keys fingerprint for a molecule
/*!
  The result is a 167-bit vector. There are 166 public keys, but
  to maintain consistency with other software packages they are
  numbered from 1.

  \param mol:    the molecule to be fingerprinted

  \return a pointer to the fingerprint. The client is
  responsible for calling delete on this.

*/
RDKIT_FINGERPRINTS_EXPORT ExplicitBitVect *getFingerprintAsBitVect(
    const ROMol &mol);
}  // namespace MACCSFingerprints
}  // namespace RDKit

#endif
