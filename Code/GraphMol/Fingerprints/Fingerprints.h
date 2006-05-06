//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef _RD_FINGERPRINTS_H_
#define _RD_FINGERPRINTS_H_

class ExplicitBitVect;
namespace RDKit{
  class ROMol;
  //! \brief Generates a topological (Daylight like) fingerprint for a molecule
  /*!

    \param mol:          the molecule to be fingerprinted
    \param minPath:      the minimum path length (in bonds) to be included
    \param maxPath:      the minimum path length (in bonds) to be included
    \param fpSize:       the size of the fingerprint
    \param nBitsPerHash: the number of bits to be set by each path
    \param useHs:        toggles inclusion of Hs in distinguishing paths from
                         each other.

    \return the molecular fingerprint, as an ExplicitBitVect

    <b>Notes:</b>
      - the caller is responsible for <tt>delete</tt>ing the result
    
  */
  ExplicitBitVect *DaylightFingerprintMol(const ROMol &mol,int minPath=1,int maxPath=7,
					  int fpSize=2048,int nBitsPerHash=4,
					  bool useHs=true);

}

#endif
