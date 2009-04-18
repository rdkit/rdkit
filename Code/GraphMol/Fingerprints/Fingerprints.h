//
//  Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef _RD_FINGERPRINTS_H_
#define _RD_FINGERPRINTS_H_

class ExplicitBitVect;
namespace RDKit{
  const std::string DaylightFingerprintMolVersion="2.0.0-deprecated";
  class ROMol;
  //! \brief Generates a topological (Daylight like) fingerprint for a molecule
  /*!
    \deprecated NOTE: This function is deprecated. Please use RDKFingerprintMol()
                instead.

    \param mol:          the molecule to be fingerprinted
    \param minPath:      the minimum path length (in bonds) to be included
    \param maxPath:      the minimum path length (in bonds) to be included
    \param fpSize:       the size of the fingerprint
    \param nBitsPerHash: the number of bits to be set by each path
    \param useHs:        toggles inclusion of Hs in distinguishing paths from
                         each other.
    \param tgtDensity:   if the generated fingerprint is below this density, it will
                         be folded until the density is reached.
    \param minSize:      the minimum size to which the fingerprint will be
                         folded

    \return the molecular fingerprint, as an ExplicitBitVect

    <b>Notes:</b>
      - the caller is responsible for <tt>delete</tt>ing the result
    
  */
  ExplicitBitVect *DaylightFingerprintMol(const ROMol &mol,
					  unsigned int minPath=1,unsigned int maxPath=7,
					  unsigned int fpSize=2048,unsigned int nBitsPerHash=4,
					  bool useHs=true,
					  double tgtDensity=0.0,unsigned int minSize=128);

  //! \brief Generates a topological (Daylight like) fingerprint for a molecule
  //!        using an alternate (faster) hashing algorithm  
  /*!

    \param mol:          the molecule to be fingerprinted
    \param minPath:      the minimum path length (in bonds) to be included
    \param maxPath:      the minimum path length (in bonds) to be included
    \param fpSize:       the size of the fingerprint
    \param nBitsPerHash: the number of bits to be set by each path
    \param useHs:        toggles inclusion of Hs in distinguishing paths from
                         each other.
    \param tgtDensity:   if the generated fingerprint is below this density, it will
                         be folded until the density is reached.
    \param minSize:      the minimum size to which the fingerprint will be
                         folded

    \return the molecular fingerprint, as an ExplicitBitVect

    <b>Notes:</b>
      - the caller is responsible for <tt>delete</tt>ing the result
    
  */
  const std::string RDKFingerprintMolVersion="1.0.0";
  ExplicitBitVect *RDKFingerprintMol(const ROMol &mol,
                                    unsigned int minPath=1,unsigned int maxPath=7,
                                    unsigned int fpSize=2048,unsigned int nBitsPerHash=4,
                                    bool useHs=true,
                                     double tgtDensity=0.0,unsigned int minSize=128);

  //! \brief Generates a topological (Daylight like) fingerprint for a molecule
  //!        using an alternate (faster) hashing algorithm  
  /*!

    <b>Experimental:</b> This function is experimental. The API or results may change from
    release to release.
    
    \param mol:          the molecule to be fingerprinted
    \param layers:       the layers to be included (see below)
    \param minPath:      the minimum path length (in bonds) to be included
    \param maxPath:      the minimum path length (in bonds) to be included
    \param fpSize:       the size of the fingerprint
    \param tgtDensity:   if the generated fingerprint is below this density, it will
                         be folded until the density is reached.
    \param minSize:      the minimum size to which the fingerprint will be
                         folded
    \param atomCounts:   if provided, this will be used to provide the count of the number
                         of paths that set bits each atom is involved in. The vector should
                         have at least as many entries as the molecule has atoms and is not
                         zeroed out here.
    \param setOnlyBits:  if provided, only bits that are set in this bit vector will be set
                         in the result. This is essentially the same as doing:
                            (*res) &= (*setOnlyBits);
                         but also has an impact on the atomCounts (if being used)

    \return the molecular fingerprint, as an ExplicitBitVect

    <b>Notes:</b>
      - the caller is responsible for <tt>delete</tt>ing the result

    <b>Layer definitions:</b>
       - 0x01: pure topology
       - 0x02: bond order
       - 0x04: atom types
       - 0x08: presence of rings
       - 0x10: ring sizes
  */
  const std::string LayeredFingerprintMolVersion="0.2.0";
  ExplicitBitVect *LayeredFingerprintMol(const ROMol &mol,
                                         unsigned int layerFlags=0xFFFFFFFF,
                                         unsigned int minPath=1,unsigned int maxPath=7,
                                         unsigned int fpSize=2048,
                                         double tgtDensity=0.0,unsigned int minSize=128,
                                         std::vector<unsigned int> *atomCounts=0,
                                         ExplicitBitVect *setOnlyBits=0);
}

#endif
