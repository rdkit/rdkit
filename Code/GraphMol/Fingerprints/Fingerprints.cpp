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
  namespace {
    bool isComplexQuery(const Bond *b){
      if( !b->hasQuery()) return false;
      // negated things are always complex:
      if( b->getQuery()->getNegation()) return true;
      std::string descr=b->getQuery()->getDescription();
      if(descr=="BondOrder") return false;
      if(descr=="BondAnd" || descr=="BondOr" || descr=="BondXor") return true;
      
      return true;
    }
    bool isComplexQuery(const Atom *a){
      if( !a->hasQuery()) return false;
      // negated things are always complex:
      if( a->getQuery()->getNegation()) return true;
      std::string descr=a->getQuery()->getDescription();
      if(descr=="AtomAtomicNum") return false;
      if(descr=="AtomOr" || descr=="AtomXor") return true;
      if(descr=="AtomAnd"){
        Queries::Query<int,Atom const *,true>::CHILD_VECT_CI childIt=a->getQuery()->beginChildren();
        if( (*childIt)->getDescription()=="AtomAtomicNum" &&
            ((*(childIt+1))->getDescription()=="AtomIsAliphatic" ||
             (*(childIt+1))->getDescription()=="AtomIsAromatic") &&
            (childIt+2)==a->getQuery()->endChildren()){
          return false;
        }
        return true;
      }
      
      return true;
    }
  } // end of anonymous namespace

  // caller owns the result, it must be deleted
  // NOTE: This function is deprecated. Please use RDKFingerprintMol() instead.
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
                                         double tgtDensity,unsigned int minSize,
                                         std::vector<unsigned int> *atomCounts,
                                         ExplicitBitVect *setOnlyBits){
    PRECONDITION(minPath!=0,"minPath==0");
    PRECONDITION(maxPath>=minPath,"maxPath<minPath");
    PRECONDITION(fpSize!=0,"fpSize==0");
    PRECONDITION(!atomCounts || atomCounts->size()>=mol.getNumAtoms(),"bad atomCounts size");
    PRECONDITION(!setOnlyBits || setOnlyBits->GetNumBits()==fpSize,"bad setOnlyBits size");

    if(!mol.getRingInfo()->isInitialized()){
      MolOps::findSSSR(mol);
    }
    
#ifdef LAYEREDFP_USE_MT
    // create a mersenne twister with customized parameters. 
    // The standard parameters (used to create boost::mt19937) 
    // result in an RNG that's much too computationally intensive
    // to seed.
    typedef boost::random::mersenne_twister<boost::uint32_t,32,4,2,31,0x9908b0df,11,7,0x9d2c5680,15,0xefc60000,18, 3346425566U>  rng_type;
    
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
#endif

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

        // should we keep this path?
        bool keepPath=true;
        for(unsigned int i=0;i<path.size();++i){
          const Bond *bi = mol.getBondWithIdx(path[i]);
          if(isComplexQuery(bi) ||
             isComplexQuery(bi->getBeginAtom()) ||
             isComplexQuery(bi->getEndAtom())) {
            keepPath=false;
            break;
          }
        }

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
            ourHash = bondNbrs[i]%8; // 3 bits here
            hashLayers[0].push_back(ourHash);
          }
          nBitsInHash+=3;
          if(layerFlags & 0x2 && keepPath){
            // layer 2: include bond orders:
            unsigned int bondHash;
            if(bi->getIsAromatic()){
              // makes sure aromatic bonds always hash the same:
              bondHash = Bond::AROMATIC;
            } else if(bi->getBondType()==Bond::SINGLE &&
                      bi->getBeginAtom()->getIsAromatic() &&
                      bi->getEndAtom()->getIsAromatic() &&
                      queryIsBondInRing(bi)
                      ){
              // a special case that comes up if we're using these to filter
              // substructure matches:
              //   This molecule: 
              //     Cn1ccc2nn(C)c(=O)c-2c1C
              //   which has a non-aromatic bridging bond between aromatic
              //   atoms, matches the SMARTS query:
              //    Cc1ncccc1
              //   because unspecified bonds in SMARTS are aromatic or single
              // We need to make sure to capture this case.  
              bondHash = Bond::AROMATIC;
            } else {
              bondHash = bi->getBondType();
            }
            ourHash = (bondHash%16)<<nBitsInHash; // 4 bits here
            hashLayers[1].push_back(ourHash);
          }
          nBitsInHash+=4;
          if(layerFlags & 0x4 && keepPath){
            //std::cerr<<" consider: "<<bi->getBeginAtomIdx()<<" - " <<bi->getEndAtomIdx()<<std::endl;
            // layer 3: include atom types:
            unsigned int a1Hash,a2Hash;
            a1Hash = (bi->getBeginAtom()->getAtomicNum()%128)<<1 | bi->getBeginAtom()->getIsAromatic();
            a2Hash = (bi->getEndAtom()->getAtomicNum()%128)<<1 | bi->getEndAtom()->getIsAromatic();
            if(a1Hash<a2Hash) std::swap(a1Hash,a2Hash);
            ourHash = a1Hash<<nBitsInHash; // 8 bits
            ourHash |= a2Hash<<(nBitsInHash+8); // 8 bits
            hashLayers[2].push_back(ourHash);
          }
          nBitsInHash += 16;
          if(layerFlags & 0x8 && keepPath){
            // layer 4: include ring information
            ourHash = queryIsBondInRing(bi)<<nBitsInHash; // 1 bit
            hashLayers[3].push_back(ourHash);
          }
          nBitsInHash++;
          if(layerFlags & 0x10 && keepPath){
            // layer 5: include ring size information
            ourHash = (queryBondMinRingSize(bi)%8)<<nBitsInHash; // 3 bits here
            hashLayers[4].push_back(ourHash);
          }
          nBitsInHash+=3;
        }
        unsigned int l=0;
        bool flaggedPath=false;
        for(std::vector< std::vector<unsigned int> >::iterator layerIt=hashLayers.begin();
            layerIt!=hashLayers.end();++layerIt,++l){
          if(!layerIt->size()) continue;
          // ----
          std::sort(layerIt->begin(),layerIt->end());
        
          // finally, we will add the number of distinct atoms in the path at the end
          // of the vect. This allows us to distinguish C1CC1 from CC(C)C
          layerIt->push_back(atomsInPath.count());

          // hash the path to generate a seed:
          unsigned long seed = vectHasher(*layerIt);
#if 0
          if(!res->GetBit(seed%fpSize)) {
            std::cerr<<"seed "<<l<<": "<<seed<<"->"<<(seed%fpSize)<<std::endl;
          } else {
            std::cerr<<"         "<<l<<": "<<seed<<"->"<<(seed%fpSize)<<std::endl;
          }
#endif
          // NOTE: since we're only generating a single number here, it seems like it
          // might make sense to just use the hash itself. In early testing of these
          // fingerprints, this led to a bunch of collisions between layers 3 and 4,

          unsigned int bitId;
#ifdef LAYEREDFP_USE_MT
          // One solution to this problem is to go back to using a PRNG:
          generator.seed(static_cast<rng_type::result_type>(seed));
          bitId=randomSource()%fpSize;
#else
          // The other solution is to shift the seed so that we look at different
          // bits for the different layers:
          //std::cerr<<"layer: "<<l<<" seed: "<<seed<<std::endl;
          seed = seed>>l;
          //std::cerr<<"   >>> "<<seed<<std::endl;
          bitId=seed%fpSize;
          //std::cerr<<"   bit "<<bitId<<std::endl;
#endif
          if(!setOnlyBits || (*setOnlyBits)[bitId]){
            res->SetBit(bitId);
            if(atomCounts && !flaggedPath){
              for(unsigned int aIdx=0;aIdx<atomsInPath.size();++aIdx){
                if(atomsInPath[aIdx]){
                  (*atomCounts)[aIdx]+=1;
                }
              }
              flaggedPath=true;
            }
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
    }
    return res;
  }



}
