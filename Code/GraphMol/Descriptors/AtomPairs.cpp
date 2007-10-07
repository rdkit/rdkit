// $Id$
//
//  Copyright (C) 2007 Greg Landrum
//
//   @@ All Rights Reserved  @@
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/AtomPairs.h>
#include <DataStructs/SparseIntVect.h>

namespace RDKit{
  namespace Descriptors {
    namespace AtomPairs {
      unsigned int numPiElectrons(const Atom *atom){
	PRECONDITION(atom,"no atom");
	unsigned int res=0;
	if(atom->getIsAromatic()){
	  res=1;
	} else if(atom->getHybridization() != Atom::SP3){
	  CHECK_INVARIANT(static_cast<unsigned int>(atom->getExplicitValence())>=atom->getDegree(),
			  "explicit valence exceeds atom degree");
	  res = atom->getExplicitValence()-atom->getDegree();
	}
	return res;
      }
    
      unsigned int getAtomCode(const Atom *atom,unsigned int branchSubtract){
	PRECONDITION(atom,"no atom");
	unsigned int code;

	unsigned int numBranches=0;
	if(atom->getDegree()>branchSubtract){
	  numBranches=atom->getDegree()-branchSubtract;
	}

	code=numBranches%maxNumBranches;
	unsigned int nPi=numPiElectrons(atom)%maxNumPi;
	code |= nPi<<numBranchBits;

	unsigned int typeIdx=0;
	unsigned int nTypes=1<<numTypeBits;
	while(typeIdx<nTypes){
	  if(atomNumberTypes[typeIdx]==static_cast<unsigned int>(atom->getAtomicNum())){
	    break;
	  } else if(atomNumberTypes[typeIdx]>static_cast<unsigned int>(atom->getAtomicNum())){
	    typeIdx=nTypes;
	    break;
	  }
	  ++typeIdx;
	}
	if(typeIdx==nTypes) --typeIdx;
	code |= typeIdx<<(numBranchBits+numPiBits);
	return code;
      };

      unsigned int getAtomPairCode(unsigned int codeI,unsigned int codeJ,
				   unsigned int dist){
	PRECONDITION(dist<maxPathLen,"dist too long");
	unsigned int res=dist;
	res |= std::min(codeI,codeJ) << numPathBits;
	res |= std::max(codeI,codeJ) << (numPathBits+codeSize);
	return res;
      }

      SparseIntVect<int> *
      getAtomPairFingerprint(const ROMol &mol){
	SparseIntVect<int> *res=new SparseIntVect<int>(1<<numAtomPairFingerprintBits);
	const double *dm = MolOps::getDistanceMat(mol);
	const unsigned int nAtoms=mol.getNumAtoms();

	std::vector<unsigned int> atomCodes;
	for(ROMol::ConstAtomIterator atomItI=mol.beginAtoms();
	    atomItI!=mol.endAtoms();++atomItI){
	  atomCodes.push_back(getAtomCode(*atomItI));
	}
	for(ROMol::ConstAtomIterator atomItI=mol.beginAtoms();
	    atomItI!=mol.endAtoms();++atomItI){
	  unsigned int i=(*atomItI)->getIdx();
	  for(ROMol::ConstAtomIterator atomItJ=atomItI+1;
	      atomItJ!=mol.endAtoms();++atomItJ){
	    unsigned int j=(*atomItJ)->getIdx();
	    unsigned int dist=static_cast<unsigned int>(floor(dm[i*nAtoms+j]));
	    if(dist<maxPathLen){
	      unsigned int bitId=getAtomPairCode(atomCodes[i],atomCodes[j],dist);
	      res->setVal(bitId,(*res)[bitId]+1);
	    }
	  }
	}
	return res;
      }

    } // end of namespace AtomPairs
  } // end of namespace Descriptors
} // end of namespace RDKit
