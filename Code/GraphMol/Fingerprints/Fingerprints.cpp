// $Id$
//
//  Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>
#include "Fingerprints.h"
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <RDGeneral/Invariant.h>
#include <boost/random.hpp>
#include <limits.h>
#include <boost/functional/hash.hpp>
#include <boost/cstdint.hpp>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>

namespace RDKit{
  // caller owns the result, it must be deleted
  ExplicitBitVect *DaylightFingerprintMol(const ROMol &mol,
                                          unsigned int minPath,
					  unsigned int maxPath,
					  unsigned int fpSize,
                                          unsigned int nBitsPerHash,
					  bool useHs,
					  double tgtDensity,
                                          unsigned int minSize){
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
	
        //boost::hash<float> floatHasher;
	//unsigned long seed = floatHasher(balabanJ);
        //unsigned long seed = *(unsigned long *)&balabanJ;
        boost::uint32_t seed = *(boost::uint32_t *)&balabanJ;
	generator.seed(static_cast<rng_type::result_type>(seed));
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
  ExplicitBitVect *RDKFingerprintMol(const ROMol &mol,unsigned int minPath,
                                      unsigned int maxPath,
                                      unsigned int fpSize,unsigned int nBitsPerHash,
                                      bool useHs,
                                      double tgtDensity,unsigned int minSize){
    PRECONDITION(minPath!=0,"minPath==0");
    PRECONDITION(maxPath>=minPath,"maxPath<minPath");
    PRECONDITION(fpSize!=0,"fpSize==0");
    PRECONDITION(nBitsPerHash!=0,"nBitsPerHash==0");

    typedef boost::mt19937 rng_type;
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
    boost::dynamic_bitset<> atomsInPath(mol.getNumAtoms());
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
        atomsInPath.reset();
        std::vector<unsigned int> bondHashes;
        bondHashes.reserve(path.size()+1);
        for(unsigned int i=0;i<path.size();++i){
          const Bond *bi = mol.getBondWithIdx(path[i]);
          atomsInPath.set(bi->getBeginAtomIdx());
          atomsInPath.set(bi->getEndAtomIdx());
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
          bondHashes.push_back(ourHash);
        }
        std::sort(bondHashes.begin(),bondHashes.end());

        // finally, we will add the number of distinct atoms in the path at the end
        // of the vect. This allows us to distinguish C1CC1 from CC(C)C
        bondHashes.push_back(atomsInPath.count());
        
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
        generator.seed(static_cast<rng_type::result_type>(seed));
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

  // caller owns the result, it must be deleted
  ExplicitBitVect *LayeredFingerprintMol(const ROMol &mol,
                                         unsigned int layerFlags,
                                         unsigned int minPath,
                                         unsigned int maxPath,
                                         unsigned int fpSize,
                                         double tgtDensity,unsigned int minSize){
    PRECONDITION(minPath!=0,"minPath==0");
    PRECONDITION(maxPath>=minPath,"maxPath<minPath");
    PRECONDITION(fpSize!=0,"fpSize==0");

    // if we generate arbitrarily sized ints then mod them down to the
    // appropriate size, we can guarantee that a fingerprint of
    // size x has the same bits set as one of size 2x that's been folded
    // in half.  This is a nice guarantee to have.
    //
    boost::hash<std::vector<unsigned int> > vectHasher;
    ExplicitBitVect *res = new ExplicitBitVect(fpSize);
    INT_PATH_LIST_MAP allPaths = findAllSubgraphsOfLengthsMtoN(mol,minPath,maxPath);
    boost::dynamic_bitset<> atomsInPath(mol.getNumAtoms());
    for(INT_PATH_LIST_MAP_CI paths=allPaths.begin();paths!=allPaths.end();paths++){
      for( PATH_LIST_CI pathIt=paths->second.begin();
	   pathIt!=paths->second.end();
	   pathIt++ ){
	const PATH_TYPE &path=*pathIt;

        std::vector< std::vector<unsigned int> > hashLayers(10);

        // calculate the number of neighbors each bond has in the path:
        std::vector<unsigned int> bondNbrs(path.size(),0);
        atomsInPath.reset();
        for(unsigned int i=0;i<path.size();++i){
          const Bond *bi = mol.getBondWithIdx(path[i]);
          atomsInPath.set(bi->getBeginAtomIdx());
          atomsInPath.set(bi->getEndAtomIdx());
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
          // we have the count of neighbors for bond bi, compute its hash layers:
          unsigned int ourHash=0;
          unsigned int nBitsInHash=0;

          if(layerFlags & 0x1){
            // layer 1: straight topology
            ourHash |= bondNbrs[i]%8; // 3 bits here
            nBitsInHash+=3;
            hashLayers[0].push_back(ourHash);
          }
          if(layerFlags & 0x2){
            // layer 2: include bond orders:
            unsigned int bondHash;
            if(bi->getIsAromatic()){
              // makes sure aromatic bonds always hash the same:
              bondHash = Bond::AROMATIC;
            } else {
              bondHash = bi->getBondType();
            }
            ourHash |= (bondHash%16)<<nBitsInHash; // 4 bits here
            nBitsInHash+=4;
            hashLayers[1].push_back(ourHash);
          }
          if(layerFlags & 0x4){
            // layer 3: include atom types:
            unsigned int a1Hash,a2Hash;
            a1Hash = (bi->getBeginAtom()->getAtomicNum()%128)<<1 | bi->getBeginAtom()->getIsAromatic();
            a2Hash = (bi->getEndAtom()->getAtomicNum()%128)<<1 | bi->getEndAtom()->getIsAromatic();
            if(a1Hash<a2Hash) std::swap(a1Hash,a2Hash);
            ourHash |= a1Hash<<nBitsInHash; // 8 bits
            nBitsInHash+=8;
            ourHash |= a2Hash<<nBitsInHash; // 8 bits
            nBitsInHash+=8;
            hashLayers[2].push_back(ourHash);
          }
          if(layerFlags & 0x8){
            // layer 4: include ring information
            ourHash |= queryIsBondInRing(bi)<<nBitsInHash; // 1 bit
            nBitsInHash++;
            hashLayers[3].push_back(ourHash);
          }
          if(layerFlags & 0x10){
            // layer 5: include ring size information
            ourHash |= (queryBondMinRingSize(bi)%8)<<nBitsInHash; // 3 bits here
            nBitsInHash+=3;
            hashLayers[4].push_back(ourHash);
          }
        }

        for(std::vector< std::vector<unsigned int> >::iterator layerIt=hashLayers.begin();
            layerIt!=hashLayers.end();++layerIt){
          if(!layerIt->size()) continue;
          // ----
          std::sort(layerIt->begin(),layerIt->end());
        
          // finally, we will add the number of distinct atoms in the path at the end
          // of the vect. This allows us to distinguish C1CC1 from CC(C)C
          layerIt->push_back(atomsInPath.count());

          // hash the path to generate a seed:
          unsigned long seed = vectHasher(*layerIt);
          res->SetBit(seed%fpSize);
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
    }
    return res;
  }



}
