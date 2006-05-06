// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <GraphMol/RDKitBase.h>
#include <DataStructs/ExplicitBitVect.h>
#include "Fingerprints.h"
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <RDGeneral/Invariant.h>
#include <boost/random.hpp>
#include <limits.h>


namespace RDKit{
  // caller owns the result, it must be deleted
  ExplicitBitVect *DaylightFingerprintMol(const ROMol &mol,int minPath,int maxPath,
					  int fpSize,int nBitsPerHash,bool useHs){
    PRECONDITION(maxPath>=minPath,"maxPath<minPath");
    PRECONDITION(fpSize>0,"bad fingerprint length");
    PRECONDITION(nBitsPerHash>0,"bad nBitsPerHash");

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
	
	//std::copy(path.begin(),path.end(),std::ostream_iterator<int>(std::cout," "));
	//std::cout << std::endl;
	//std::cout << " " << balabanJ << std::endl;
	unsigned long seed = *(unsigned long *)(&balabanJ);
	//std::cout << MolToSmiles(*subMol) << " " << seed << std::endl;
	//std::cout << "\t" << seed << std::endl;
	generator.seed(seed);

	//std::cout << "\t\t";
	for(int i=0;i<nBitsPerHash;i++){
	  unsigned int bit = randomSource();
	  bit %= fpSize;
	  //std::cout << bit << " ";
	  res->SetBit(bit);
	}
      }
    }
    return res;
  }



}
