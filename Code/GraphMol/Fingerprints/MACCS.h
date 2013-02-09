/*! \file MACCS.h

*/
#ifndef __RD_MACCSFPS_H__
#define __RD_MACCSFPS_H__

class ExplicitBitVect;
namespace RDKit {
  class ROMol;
  namespace MACCSFingerprints {
    ExplicitBitVect *getFingerprintAsBitVect(const ROMol &mol);
  }
}

#endif
