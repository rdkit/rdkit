//
//  Copyright (C) 2007-2013 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*! \file AtomPairs.h


  A few quick notes about fingerprint size and the way chirality is handled in these functions.

  By default the atom-pair and topologic-torsion fingerprints do not include any information about
  chirality; the atom invariants only include information about the atomic number,
  number of pi electrons, and degree.
  When chirality is included, two additional bits are added to the atom invariants to flag R/S/no
  chirality. These additional bits change the size of the atom invariants and either the size
  of the final fingerprint (atom pairs) or the maximum allowed path length (torsions). This means
  that even fingerprints for achiral molecules are different when includeChirality is true.
  
*/
#ifndef __RD_ATOMPAIRS_H__
#define __RD_ATOMPAIRS_H__

#include <DataStructs/SparseIntVect.h>
#include <DataStructs/BitVects.h>
#include <boost/cstdint.hpp>
namespace RDKit {
  class Atom;

  namespace AtomPairs {
    const std::string atomPairsVersion="1.1.0";
    const unsigned int numTypeBits=4;
    const unsigned int atomNumberTypes[1<<numTypeBits]={5,6,7,8,9,14,15,16,17,33,34,35,51,52,43};
    const unsigned int numPiBits=2;
    const unsigned int maxNumPi=(1<<numPiBits)-1;
    const unsigned int numBranchBits=3;
    const unsigned int maxNumBranches=(1<<numBranchBits)-1;
    const unsigned int numChiralBits=2;
    const unsigned int codeSize=numTypeBits+numPiBits+numBranchBits;
    const unsigned int numPathBits=5;
    const unsigned int maxPathLen=(1<<numPathBits)-1;
    const unsigned int numAtomPairFingerprintBits=numPathBits+2*codeSize; // note that this is only accurate if chirality is not included
    
    //! returns a numeric code for the atom (the atom's hash in the
    //! atom-pair scheme)
    /*!
      \param atom            the atom to be considered
      \param branchSubtract  (optional) a constant to subtract from
      the number of neighbors when the hash
      is calculated (used in the topological
      torsions code)
      \param includeChirality toggles the inclusions of bits indicating R/S chirality
    */
    boost::uint32_t getAtomCode(const Atom *atom,unsigned int branchSubtract=0,bool includeChirality=false);

    //! returns an atom pair hash based on two atom hashes and the
    //! distance between the atoms.
    /*!
      \param codeI  the hash for the first atom
      \param codeJ  the hash for the second atom
      \param dist   the distance (number of bonds) between the two
      atoms
      \param includeChirality toggles the inclusions of bits indicating R/S chirality
    */
    boost::uint32_t getAtomPairCode(boost::uint32_t codeI,boost::uint32_t codeJ,
                                    unsigned int dist,bool includeChirality=false);

    //! returns the atom-pair fingerprint for a molecule
    /*!
      The algorithm used is described here:
      R.E. Carhart, D.H. Smith, R. Venkataraghavan; "Atom Pairs as
      Molecular Features in Structure-Activity Studies: Definition
      and Applications" JCICS 25, 64-73 (1985).

      
      \param mol:   the molecule to be fingerprinted
      \param minLength:   minimum distance between atoms to be
                          considered in a pair. Default is 1 bond.
      \param maxLength:   maximum distance between atoms to be
                          considered in a pair.
                          Default is maxPathLen-1 bonds.
      \param fromAtoms:   if provided, only atom pairs that involve
                          the specified atoms will be included in the
                          fingerprint
      \param ignoreAtoms: if provided, any atom pairs that include 
                          the specified atoms will not be included in the
                          fingerprint
      \param atomInvariants: a list of invariants to use for the atom hashes
                             note: only the first \c codeSize bits of each
                             invariant are used.
      \param includeChirality: if set, chirality will be used in the atom invariants
                               (note: this is ignored if atomInvariants are provided)  
      \param use2D:       if set, the 2D (topological) distance matrix is used.
      

      \return a pointer to the fingerprint. The client is
      responsible for calling delete on this.

    */
    SparseIntVect<boost::int32_t> *
    getAtomPairFingerprint(const ROMol &mol,
                           unsigned int minLength,unsigned int maxLength,
                           const std::vector<boost::uint32_t> *fromAtoms=0,
                           const std::vector<boost::uint32_t> *ignoreAtoms=0,
                           const std::vector<boost::uint32_t> *atomInvariants=0,
                           bool includeChirality=false,
                           bool use2D=true);
    //! \overload
    SparseIntVect<boost::int32_t> *
    getAtomPairFingerprint(const ROMol &mol,
                           const std::vector<boost::uint32_t> *fromAtoms=0,
                           const std::vector<boost::uint32_t> *ignoreAtoms=0,
                           const std::vector<boost::uint32_t> *atomInvariants=0,
                           bool includeChirality=false,
                           bool use2D=true);


    //! returns the hashed atom-pair fingerprint for a molecule
    /*!
      \param mol:   the molecule to be fingerprinted
      \param nBits:   the length of the fingerprint to generate
      \param minLength:   minimum distance between atoms to be
                          considered in a pair. Default is 1 bond.
      \param maxLength:   maximum distance between atoms to be
                          considered in a pair.
                          Default is maxPathLen-1 bonds.
      \param fromAtoms:   if provided, only atom pairs that involve
                          the specified atoms will be included in the
                          fingerprint
      \param ignoreAtoms: if provided, any atom pairs that include 
                          the specified atoms will not be included in the
                          fingerprint
      \param atomInvariants: a list of invariants to use for the atom hashes
                             note: only the first \c codeSize bits of each
                             invariant are used.
      \param includeChirality: if set, chirality will be used in the atom invariants
                               (note: this is ignored if atomInvariants are provided)  
      \param use2D:       if set, the 2D (topological) distance matrix is used.

      \return a pointer to the fingerprint. The client is
      responsible for calling delete on this.

    */
    SparseIntVect<boost::int32_t> *
    getHashedAtomPairFingerprint(const ROMol &mol,
                                 unsigned int nBits=2048,
                                 unsigned int minLength=1,
                                 unsigned int maxLength=maxPathLen-1,
                                 const std::vector<boost::uint32_t> *fromAtoms=0,
                                 const std::vector<boost::uint32_t> *ignoreAtoms=0,
                                 const std::vector<boost::uint32_t> *atomInvariants=0,
                                 bool includeChirality=false,
                                 bool use2D=true);
    //! returns the hashed atom-pair fingerprint for a molecule as a bit vector
    /*!
      \param mol:   the molecule to be fingerprinted
      \param nBits:   the length of the fingerprint to generate
      \param minLength:   minimum distance between atoms to be
                          considered in a pair. Default is 1 bond.
      \param maxLength:   maximum distance between atoms to be
                          considered in a pair.
                          Default is maxPathLen-1 bonds.
      \param fromAtoms:   if provided, only atom pairs that involve
                          the specified atoms will be included in the
                          fingerprint
      \param ignoreAtoms: if provided, any atom pairs that include 
                          the specified atoms will not be included in the
                          fingerprint
      \param atomInvariants: a list of invariants to use for the atom hashes
                             note: only the first \c codeSize bits of each
                             invariant are used.
      \param nBitsPerEntry: number of bits to use in simulating counts
      \param includeChirality: if set, chirality will be used in the atom invariants
                               (note: this is ignored if atomInvariants are provided)  
      \param use2D:       if set, the 2D (topological) distance matrix is used.

      \return a pointer to the fingerprint. The client is
      responsible for calling delete on this.

    */
    ExplicitBitVect *
    getHashedAtomPairFingerprintAsBitVect(const ROMol &mol,
                                          unsigned int nBits=2048,
                                          unsigned int minLength=1,
                                          unsigned int maxLength=maxPathLen-1,
                                          const std::vector<boost::uint32_t> *fromAtoms=0,
                                          const std::vector<boost::uint32_t> *ignoreAtoms=0,
                                          const std::vector<boost::uint32_t> *atomInvariants=0,
                                          unsigned int nBitsPerEntry=4,
                                          bool includeChirality=false,
                                          bool use2D=true);
                                          


    //! returns an topological torsion hash based on the atom hashes
    //! passed in
    /*!
      \param atomCodes  the vector of atom hashes
    */
    boost::uint64_t getTopologicalTorsionCode(const std::vector<boost::uint32_t> &atomCodes,bool includeChirality=false);

    //! returns the topological-torsion fingerprint for a molecule
    /*!
      The algorithm used is described here:
      R. Nilakantan, N. Bauman, J. S. Dixon, R. Venkataraghavan;
      "Topological Torsion: A New Molecular Descriptor for SAR Applications.
      Comparison with Other Descriptors" JCICS 27, 82-85 (1987).

      \param mol:         the molecule to be fingerprinted
      \param targetSize:  the number of atoms to include in the "torsions"
      \param fromAtoms:   if provided, only torsions that start or end at
                          the specified atoms will be included in the
                          fingerprint
      \param ignoreAtoms: if provided, any torsions that include 
                          the specified atoms will not be included in the
                          fingerprint
      \param atomInvariants: a list of invariants to use for the atom hashes
                             note: only the first \c codeSize bits of each
                             invariant are used.
      \param includeChirality: if set, chirality will be used in the atom invariants
                               (note: this is ignored if atomInvariants are provided)  

      \return a pointer to the fingerprint. The client is
      responsible for calling delete on this.

    */
    SparseIntVect<boost::int64_t > *
    getTopologicalTorsionFingerprint(const ROMol &mol,
                                     unsigned int targetSize=4,
                                     const std::vector<boost::uint32_t> *fromAtoms=0,
                                     const std::vector<boost::uint32_t> *ignoreAtoms=0,
                                     const std::vector<boost::uint32_t> *atomInvariants=0,
                                     bool includeChirality=false
                                     );
    //! returns a hashed topological-torsion fingerprint for a molecule
    /*!
      The algorithm used is described here:
      R. Nilakantan, N. Bauman, J. S. Dixon, R. Venkataraghavan;
      "Topological Torsion: A New Molecular Descriptor for SAR Applications.
      Comparison with Other Descriptors" JCICS 27, 82-85 (1987).

      \param mol:         the molecule to be fingerprinted
      \param nBits:       number of bits to include in the fingerprint
      \param targetSize:  the number of atoms to include in the "torsions"
      \param fromAtoms:   if provided, only torsions that start or end at
                          the specified atoms will be included in the
                          fingerprint
      \param ignoreAtoms: if provided, any torsions that include 
                          the specified atoms will not be included in the
                          fingerprint
      \param atomInvariants: a list of invariants to use for the atom hashes
                             note: only the first \c codeSize bits of each
                             invariant are used.
      \param includeChirality: if set, chirality will be used in the atom invariants
                               (note: this is ignored if atomInvariants are provided)  

      \return a pointer to the fingerprint. The client is
      responsible for calling delete on this.

    */
    SparseIntVect<boost::int64_t > *
    getHashedTopologicalTorsionFingerprint(const ROMol &mol,
                                           unsigned int nBits=2048,
                                           unsigned int targetSize=4,
                                           const std::vector<boost::uint32_t> *fromAtoms=0,
                                           const std::vector<boost::uint32_t> *ignoreAtoms=0,
                                           const std::vector<boost::uint32_t> *atomInvariants=0,
                                           bool includeChirality=false);
    //! returns a hashed topological-torsion fingerprint for a molecule as a bit vector
    /*!
      \param mol:         the molecule to be fingerprinted
      \param nBits:       number of bits to include in the fingerprint
      \param targetSize:  the number of atoms to include in the "torsions"
      \param fromAtoms:   if provided, only torsions that start or end at
                          the specified atoms will be included in the
                          fingerprint
      \param ignoreAtoms: if provided, any torsions that include 
                          the specified atoms will not be included in the
                          fingerprint
      \param atomInvariants: a list of invariants to use for the atom hashes
                             note: only the first \c codeSize bits of each
                             invariant are used.
      \param nBitsPerEntry: number of bits to use in simulating counts
      \param includeChirality: if set, chirality will be used in the atom invariants
                               (note: this is ignored if atomInvariants are provided)  

      \return a pointer to the fingerprint. The client is
      responsible for calling delete on this.

    */
    ExplicitBitVect *
    getHashedTopologicalTorsionFingerprintAsBitVect(const ROMol &mol,
                                                    unsigned int nBits=2048,
                                                    unsigned int targetSize=4,
                                                    const std::vector<boost::uint32_t> *fromAtoms=0,
                                                    const std::vector<boost::uint32_t> *ignoreAtoms=0,
                                                    const std::vector<boost::uint32_t> *atomInvariants=0,
                                                    unsigned int nBitsPerEntry=4,
                                                    bool includeChirality=false);
  }    
}

#endif
