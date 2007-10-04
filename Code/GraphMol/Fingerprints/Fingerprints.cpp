// $Id$
//
//  Copyright (C) 2003-2007 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <GraphMol/RDKitBase.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>
#include "Fingerprints.h"
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <RDGeneral/Invariant.h>
#include <boost/random.hpp>
#include <limits.h>


namespace RDKit{
  // caller owns the result, it must be deleted
  ExplicitBitVect *DaylightFingerprintMol(const ROMol &mol,unsigned int minPath,
					  unsigned int maxPath,
					  unsigned int fpSize,unsigned int nBitsPerHash,
					  bool useHs,
					  double tgtDensity,unsigned int minSize){
    PRECONDITION(minPath!=0,"minPath==0");
    PRECONDITION(maxPath>=minPath,"maxPath<minPath");
    PRECONDITION(fpSize!=0,"fpSize==0");
    PRECONDITION(nBitsPerHash!=0,"nBitsPerHash==0");

    typedef boost::minstd_rand rng_type;
    typedef boost::uniform_int<> distrib_type;
    typedef boost::variate_generator<rng_type &,distrib_type> source_type;
    rng_type generator(42u);

    //
    // if we generate arbitrarily sized ints then mod them down to the
    // appropriate size, we can guarantee that a fingerprint of
    // size x has the same bits set as one of size 2x that's been folded
    // in half.  This is a nice guarantee to have.
    //
    distrib_type dist(0,INT_MAX);
    source_type randomSource(generator,dist);

    ExplicitBitVect *res = new ExplicitBitVect(fpSize);
    INT_PATH_LIST_MAP allPaths = findAllSubgraphsOfLengthsMtoN(mol,minPath,maxPath,
							       useHs);
    for(INT_PATH_LIST_MAP_CI paths=allPaths.begin();paths!=allPaths.end();paths++){
      for( PATH_LIST_CI pathIt=paths->second.begin();
	   pathIt!=paths->second.end();
	   pathIt++ ){
	const PATH_TYPE &path=*pathIt;
	float balabanJ = static_cast<float>(MolOps::computeBalabanJ(mol,true,true,
								    &path,false));
	
	unsigned long seed = *(unsigned long *)(&balabanJ);
	generator.seed(seed);

	for(int i=0;i<nBitsPerHash;i++){
	  unsigned int bit = randomSource();
	  bit %= fpSize;
	  res->SetBit(bit);
	}
      }
    }

    // EFF: this could be faster by folding by more than a factor
    // of 2 each time, but we're not going to be spending much
    // time here anyway
    while( static_cast<double>(res->GetNumOnBits())/res->GetNumBits() < tgtDensity &&
	   res->GetNumBits() >= 2*minSize ){
      ExplicitBitVect *tmpV=FoldFingerprint(*res,2);
      delete res;
      res = tmpV;
    }
    
    return res;
  }



}
