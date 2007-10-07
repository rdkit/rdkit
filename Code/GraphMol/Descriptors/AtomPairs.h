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
    
      unsigned int getAtomCode(const Atom *atom,unsigned int branchSubtract=0);

      unsigned int getAtomPairCode(unsigned int codeI,unsigned int codeJ,
				   unsigned int dist);
      SparseIntVect<int> *getAtomPairFingerprint(const ROMol &mol);
    }    
  } // end of namespace Descriptors
}

#endif
