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
    \param tgtDensity:   
    \param minSize:    
                         

    \return the molecular fingerprint, as an ExplicitBitVect

    <b>Notes:</b>
      - the caller is responsible for <tt>delete</tt>ing the result
    
  */
  ExplicitBitVect *DaylightFingerprintMol(const ROMol &mol,
					  unsigned int minPath=1,unsigned int maxPath=7,
					  unsigned int fpSize=2048,unsigned int nBitsPerHash=4,
					  bool useHs=true,
					  double tgtDensity=0.0,unsigned int minSize=128);

}

#endif
