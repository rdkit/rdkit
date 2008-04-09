// $Id$
//
//  Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
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
#include <boost/functional/hash.hpp>
#include <algorithm>

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

	for(unsigned int i=0;i<nBitsPerHash;i++){
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


  // caller owns the result, it must be deleted
  ExplicitBitVect *RDKFingerprintMol2(const ROMol &mol,unsigned int minPath,
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
#ifdef VERBOSE_FINGERPRINTING        
        std::cerr<<"Path: ";
        std::copy(path.begin(),path.end(),std::ostream_iterator<int>(std::cerr,", "));
        std::cerr<<std::endl;
#endif
        // initialize the bond hashes to the number of neighbors the bond has in the path:
        std::vector<unsigned int> bondNbrs(path.size());
        std::fill(bondNbrs.begin(),bondNbrs.end(),0);
        std::set<unsigned int> atomsInPath;
        std::vector<unsigned int> bondHashes;
        bondHashes.reserve(path.size()+1);
        for(unsigned int i=0;i<path.size();++i){
          const Bond *bi = mol.getBondWithIdx(path[i]);
          atomsInPath.insert(bi->getBeginAtomIdx());
          atomsInPath.insert(bi->getEndAtomIdx());
          for(unsigned int j=i+1;j<path.size();++j){
            const Bond *bj = mol.getBondWithIdx(path[j]);
            if(bi->getBeginAtomIdx()==bj->getBeginAtomIdx() ||
               bi->getBeginAtomIdx()==bj->getEndAtomIdx() ||
               bi->getEndAtomIdx()==bj->getBeginAtomIdx() ||
               bi->getEndAtomIdx()==bj->getEndAtomIdx() ){
              ++bondNbrs[i];
              ++bondNbrs[j];
            }
          }
#ifdef VERBOSE_FINGERPRINTING        
          std::cerr<<"   bond("<<i<<"):"<<bondNbrs[i]<<std::endl;
#endif
          // we have the count of neighbors for bond bi, compute its hash:
          unsigned int a1Hash,a2Hash;
          a1Hash = (bi->getBeginAtom()->getAtomicNum()%128)<<1 | bi->getBeginAtom()->getIsAromatic();
          a2Hash = (bi->getEndAtom()->getAtomicNum()%128)<<1 | bi->getEndAtom()->getIsAromatic();
          if(a1Hash<a2Hash) std::swap(a1Hash,a2Hash);
          unsigned int bondHash;
          if(bi->getIsAromatic()){
            // makes sure aromatic bonds always hash the same:
            bondHash = Bond::AROMATIC;
          } else {
            bondHash = bi->getBondType();
          }
          unsigned int nBitsInHash=0;
          unsigned int ourHash=bondNbrs[i]%8; // 3 bits here
          nBitsInHash+=3;
          ourHash |= (bondHash%16)<<nBitsInHash; // 4 bits here
          nBitsInHash+=4;
          ourHash |= a1Hash<<nBitsInHash; // 8 bits
          nBitsInHash+=8;
          ourHash |= a2Hash<<nBitsInHash; // 8 bits
          bondHashes.insert(std::lower_bound(bondHashes.begin(),bondHashes.end(),ourHash),ourHash);
        }

        // finally, we will add the number of distinct atoms in the path at the end
        // of the vect. This allows us to distinguish C1CC1 from CC(C)C
        bondHashes.push_back(atomsInPath.size());
        
        // hash the path to generate a seed:
        boost::hash<std::vector<unsigned int> > vectHasher;
	unsigned long seed = vectHasher(bondHashes);

#ifdef VERBOSE_FINGERPRINTING        
        std::cerr<<" hash: "<<seed<<std::endl;
#endif
        // originally it seemed like a good idea to track hashes we've already
        // seen in order to avoid resetting them. In some benchmarking I did, that
        // seemed to actually result in a longer runtime (at least when using
        // an std::set to store the hashes)
        generator.seed(seed);
        for(unsigned int i=0;i<nBitsPerHash;i++){
          unsigned int bit = randomSource();
          bit %= fpSize;
          res->SetBit(bit);
        }
      }
    }

    // EFF: this could be faster by folding by more than a factor
    // of 2 each time, but we're not going to be spending much
    // time here anyway
    if(tgtDensity>0.0){
      while( static_cast<double>(res->GetNumOnBits())/res->GetNumBits() < tgtDensity &&
             res->GetNumBits() >= 2*minSize ){
        ExplicitBitVect *tmpV=FoldFingerprint(*res,2);
        delete res;
        res = tmpV;
      }
    }
    
    return res;
  }



}
