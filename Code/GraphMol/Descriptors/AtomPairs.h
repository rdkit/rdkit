//
//  Copyright (C) 2007 Greg Landrum
//
//   @@ All Rights Reserved  @@
//

/*! \file AtomPairs.h

  \brief Use MolDescriptors.h in client code.

*/
#ifndef __RD_ATOMPAIRS_H__
#define __RD_ATOMPAIRS_H__

#include <DataStructs/SparseIntVect.h>
namespace RDKit {
  class Atom;

  namespace Descriptors {
    namespace AtomPairs {
      const std::string atomPairsVersion="1.0.0";
      const unsigned int numTypeBits=4;
      const unsigned int atomNumberTypes[1<<numTypeBits]={5,6,7,8,9,14,15,16,17,33,34,35,51,52,43};
      const unsigned int numPiBits=2;
      const unsigned int maxNumPi=(1<<numPiBits)-1;
      const unsigned int numBranchBits=3;
      const unsigned int maxNumBranches=(1<<numBranchBits)-1;
      const unsigned int codeSize=numTypeBits+numPiBits+numBranchBits;
      const unsigned int numPathBits=5;
      const unsigned int maxPathLen=(1<<numPathBits)-1;
      const unsigned int numAtomPairFingerprintBits=numPathBits+2*codeSize;
    
      //! returns a numeric code for the atom (the atom's hash in the
      //! atom-pair scheme)
      /*!
	\param atom            the atom to be considered
	\param branchSubtract  (optional) a constant to subtract from
                               the number of neighbors when the hash
                               is calculated (used in the topological
                               torsions code)
      */
      unsigned int getAtomCode(const Atom *atom,unsigned int branchSubtract=0);

      //! returns an atom pair hash based on two atom hashes and the
      //! distance between the atoms.
      /*!
	\param codeI  the hash for the first atom
	\param codeJ  the hash for the second atom
	\param dist   the distance (number of bonds) between the two
                      atoms
       */
      unsigned int getAtomPairCode(unsigned int codeI,unsigned int codeJ,
				   unsigned int dist);

      //! returns the atom-pair fingerprint for a molecule
      /*!
        The algorithm used is described here:
        R.E. Carhart, D.H. Smith, R. Venkataraghavan; "Atom Pairs as
          Molecular Features in Structure-Activity Studies: Definition
	  and Applications" JCICS 25, 64-73 (1985).

      
        \param mol:   the molecule to be fingerprinted
	\return a pointer to the fingerprint. The client is
                responsible for calling delete on this.

       */
      SparseIntVect<int> *getAtomPairFingerprint(const ROMol &mol);

      //! returns an topological torsion hash based on the atom hashes
      //! passed in
      /*!
	\param atomCodes  the vector of atom hashes
       */
      unsigned long long int getTopologicalTorsionCode(const std::vector<unsigned int> &atomCodes);

      //! returns the topological-torsion fingerprint for a molecule
      /*!
        The algorithm used is described here:
        R. Nilakantan, N. Bauman, J. S. Dixon, R. Venkataraghavan;
        "Topological Torsion: A New Molecular Descriptor for SAR Applications.
        Comparison with Other Descriptors" JCICS 27, 82-85 (1987).

        \param mol:         the molecule to be fingerprinted
        \param targetSize:  the number of atoms to include in the torsions

	\return a pointer to the fingerprint. The client is
                responsible for calling delete on this.

       */
      SparseIntVect<long long int> *getTopologicalTorsionFingerprint(const ROMol &mol,
								     unsigned int targetSize=4);
    }    
  } // end of namespace Descriptors
}

#endif
