// $Id$
//
//  Copyright (C) 2007-2011 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <DataStructs/SparseIntVect.h>
#include <RDGeneral/hash/hash.hpp>
#include <boost/cstdint.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/foreach.hpp>

namespace RDKit{
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
    
    boost::uint32_t getAtomCode(const Atom *atom,unsigned int branchSubtract){
      PRECONDITION(atom,"no atom");
      boost::uint32_t code;

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

    boost::uint32_t getAtomPairCode(boost::uint32_t codeI,boost::uint32_t codeJ,
                                    unsigned int dist){
      PRECONDITION(dist<maxPathLen,"dist too long");
      boost::uint32_t res=dist;
      res |= std::min(codeI,codeJ) << numPathBits;
      res |= std::max(codeI,codeJ) << (numPathBits+codeSize);
      return res;
    }

    template <typename T1,typename T2>
    void updateElement(SparseIntVect<T1> &v,T2 elem){
      v.setVal(elem,v.getVal(elem)+1);
    }

    template <typename T1>
    void updateElement(ExplicitBitVect &v,T1 elem){
      v.setBit(elem%v.getNumBits());
    }

    template <typename T>
    void setAtomPairBit(boost::uint32_t i, boost::uint32_t j,boost::uint32_t nAtoms,
                        const std::vector<boost::uint32_t> &atomCodes,
                        const double *dm,T *bv,
                        unsigned int minLength,unsigned int maxLength){
      unsigned int dist=static_cast<unsigned int>(floor(dm[i*nAtoms+j]));
      if(dist>=minLength && dist<=maxLength){
        boost::uint32_t bitId=getAtomPairCode(atomCodes[i],atomCodes[j],dist);
        updateElement(*bv,static_cast<boost::uint32_t>(bitId));
      }
    }

    SparseIntVect<boost::int32_t> *getAtomPairFingerprint(const ROMol &mol,
                                                          const std::vector<boost::uint32_t> *fromAtoms,
                                                          const std::vector<boost::uint32_t> *ignoreAtoms
                                                          ){
      return getAtomPairFingerprint(mol,1,maxPathLen-1,fromAtoms);
    };

    SparseIntVect<boost::int32_t> *
    getAtomPairFingerprint(const ROMol &mol,unsigned int minLength,unsigned int maxLength,
                           const std::vector<boost::uint32_t> *fromAtoms,
                           const std::vector<boost::uint32_t> *ignoreAtoms
                           ){
      PRECONDITION(minLength<=maxLength,"bad lengths provided");
      SparseIntVect<boost::int32_t> *res=new SparseIntVect<boost::int32_t>(1<<numAtomPairFingerprintBits);
      const double *dm = MolOps::getDistanceMat(mol);
      const unsigned int nAtoms=mol.getNumAtoms();

      std::vector<boost::uint32_t> atomCodes;
      for(ROMol::ConstAtomIterator atomItI=mol.beginAtoms();
          atomItI!=mol.endAtoms();++atomItI){
        atomCodes.push_back(getAtomCode(*atomItI));
      }
      for(ROMol::ConstAtomIterator atomItI=mol.beginAtoms();
          atomItI!=mol.endAtoms();++atomItI){
        unsigned int i=(*atomItI)->getIdx();
        if(ignoreAtoms &&
           std::find(ignoreAtoms->begin(),ignoreAtoms->end(),i)!=ignoreAtoms->end()){
          continue;
        }
        if(!fromAtoms){
          for(ROMol::ConstAtomIterator atomItJ=atomItI+1;
              atomItJ!=mol.endAtoms();++atomItJ){
            unsigned int j=(*atomItJ)->getIdx();
            if(ignoreAtoms &&
               std::find(ignoreAtoms->begin(),ignoreAtoms->end(),j)!=ignoreAtoms->end()){
              continue;
            }
            setAtomPairBit(i,j,nAtoms,atomCodes,dm,res,minLength,maxLength);
          }
        } else {
          BOOST_FOREACH(boost::uint32_t j,*fromAtoms){
            if(j!=i){
              if(ignoreAtoms &&
                 std::find(ignoreAtoms->begin(),ignoreAtoms->end(),j)!=ignoreAtoms->end()){
                continue;
              }
              setAtomPairBit(i,j,nAtoms,atomCodes,dm,res,minLength,maxLength);
            }
          }
        }
      }
      return res;
    }

    SparseIntVect<boost::int32_t> *
    getHashedAtomPairFingerprint(const ROMol &mol,unsigned int nBits,
                                 unsigned int minLength,unsigned int maxLength,
                                 const std::vector<boost::uint32_t> *fromAtoms,
                                 const std::vector<boost::uint32_t> *ignoreAtoms
                                 ){
      PRECONDITION(minLength<=maxLength,"bad lengths provided");
      SparseIntVect<boost::int32_t> *res=new SparseIntVect<boost::int32_t>(nBits);
      const double *dm = MolOps::getDistanceMat(mol);
      const unsigned int nAtoms=mol.getNumAtoms();

      std::vector<boost::uint32_t> atomCodes;
      for(ROMol::ConstAtomIterator atomItI=mol.beginAtoms();
          atomItI!=mol.endAtoms();++atomItI){
        atomCodes.push_back(getAtomCode(*atomItI));
      }
      for(ROMol::ConstAtomIterator atomItI=mol.beginAtoms();
          atomItI!=mol.endAtoms();++atomItI){
        unsigned int i=(*atomItI)->getIdx();
        if(ignoreAtoms &&
           std::find(ignoreAtoms->begin(),ignoreAtoms->end(),i)!=ignoreAtoms->end()){
          continue;
        }
        if(!fromAtoms){
          for(ROMol::ConstAtomIterator atomItJ=atomItI+1;
              atomItJ!=mol.endAtoms();++atomItJ){
            unsigned int j=(*atomItJ)->getIdx();
            if(ignoreAtoms &&
               std::find(ignoreAtoms->begin(),ignoreAtoms->end(),j)!=ignoreAtoms->end()){
              continue;
            }
            unsigned int dist=static_cast<unsigned int>(floor(dm[i*nAtoms+j]));
            if(dist>=minLength && dist<=maxLength){
              boost::uint32_t bit=0;
              gboost::hash_combine(bit,std::min(atomCodes[i],atomCodes[j]));
              gboost::hash_combine(bit,dist);
              gboost::hash_combine(bit,std::max(atomCodes[i],atomCodes[j]));
              updateElement(*res,static_cast<boost::int32_t>(bit%nBits));
            }
          }
        } else {
          BOOST_FOREACH(boost::uint32_t j,*fromAtoms){
            if(j!=i){
              if(ignoreAtoms &&
                 std::find(ignoreAtoms->begin(),ignoreAtoms->end(),j)!=ignoreAtoms->end()){
                continue;
              }
              unsigned int dist=static_cast<unsigned int>(floor(dm[i*nAtoms+j]));
              if(dist>=minLength && dist<=maxLength){
                boost::uint32_t bit=0;
                gboost::hash_combine(bit,std::min(atomCodes[i],atomCodes[j]));
                gboost::hash_combine(bit,dist);
                gboost::hash_combine(bit,std::max(atomCodes[i],atomCodes[j]));
                updateElement(*res,static_cast<boost::int32_t>(bit%nBits));
              }
            }
          }
        }
      }
      return res;
    }

    ExplicitBitVect *
    getHashedAtomPairFingerprintAsBitVect(const ROMol &mol,unsigned int nBits,
                                          unsigned int minLength,unsigned int maxLength,
                                          const std::vector<boost::uint32_t> *fromAtoms,
                                           const std::vector<boost::uint32_t> *ignoreAtoms,
                                           unsigned int nBitsPerEntry
                                          ){
      PRECONDITION(minLength<=maxLength,"bad lengths provided");
      static int bounds[4] = {1,2,4,8};
      SparseIntVect<boost::int32_t> *sres=getHashedAtomPairFingerprint(mol,nBits,minLength,maxLength,
                                                                       fromAtoms,ignoreAtoms);
      ExplicitBitVect *res=new ExplicitBitVect(nBits*nBitsPerEntry);

      if(nBitsPerEntry!=4){
        BOOST_FOREACH(SparseIntVect<boost::int64_t>::StorageType::value_type val,sres->getNonzeroElements()){
          for(unsigned int i=0;i<nBitsPerEntry;++i){
            if(val.second>i) res->setBit(val.first*nBitsPerEntry+i);
          }        
        }
      } else {
        BOOST_FOREACH(SparseIntVect<boost::int64_t>::StorageType::value_type val,sres->getNonzeroElements()){
          for(unsigned int i=0;i<nBitsPerEntry;++i){
            if(val.second>=bounds[i]){
              res->setBit(val.first*nBitsPerEntry+i);
            }
          }        
        }
      }
      delete sres;
      return res;
    }      

    
    boost::uint64_t 
    getTopologicalTorsionCode(const std::vector<boost::uint32_t> &pathCodes){
      bool reverseIt=false;
      unsigned int i=0;
      unsigned int j=pathCodes.size()-1;
      while(i<j){
        if(pathCodes[i]>pathCodes[j]){
          reverseIt=true;
          break;
        } else if( pathCodes[i]<pathCodes[j]){
          break;
        }
        ++i;
        --j;
      }

      boost::uint64_t res=0;
      if(reverseIt){
        //std::cerr<<"r";
        for(unsigned int i=0;i<pathCodes.size();++i){
          res |= static_cast<boost::uint64_t>(pathCodes[pathCodes.size()-i-1])<<(codeSize*i);
        }
      }else{
        //std::cerr<<" ";
        for(unsigned int i=0;i<pathCodes.size();++i){
          res |= static_cast<boost::uint64_t>(pathCodes[i])<<(codeSize*i);
        }
      }
      //for(unsigned int i=0;i<pathCodes.size();++i){
      //  std::cerr<<atomCodes[i]<<" ";
      //}
      //std::cerr<<res<<std::endl;
        
      return res;
    }

    size_t
    getTopologicalTorsionHash(const std::vector<boost::uint32_t> &pathCodes){
      bool reverseIt=false;
      unsigned int i=0;
      unsigned int j=pathCodes.size()-1;
      while(i<j){
        if(pathCodes[i]>pathCodes[j]){
          reverseIt=true;
          break;
        } else if( pathCodes[i]<pathCodes[j]){
          break;
        }
        ++i;
        --j;
      }

      boost::uint32_t res=0;
      if(reverseIt){
        for(unsigned int i=0;i<pathCodes.size();++i){
          gboost::hash_combine(res,pathCodes[pathCodes.size()-i-1]);
        }
      }else{
        for(unsigned int i=0;i<pathCodes.size();++i){
          gboost::hash_combine(res,pathCodes[i]);
        }
      }
      return res;
    }

    SparseIntVect<boost::int64_t> *
    getTopologicalTorsionFingerprint(const ROMol &mol,unsigned int targetSize,
                                     const std::vector<boost::uint32_t> *fromAtoms,
                                     const std::vector<boost::uint32_t> *ignoreAtoms
                                     ){
      boost::uint64_t sz=1;
      sz=(sz<<(targetSize*codeSize));
      // NOTE: this -1 is incorrect but it's needed for backwards compatibility.
      //  hopefully we'll never have a case with a torsion that hits this.
      //
      //  mmm, bug compatible.
      sz-=1;
      SparseIntVect<boost::int64_t> *res=new SparseIntVect<boost::int64_t>(sz);

      std::vector<boost::uint32_t> atomCodes;
      for(ROMol::ConstAtomIterator atomItI=mol.beginAtoms();
          atomItI!=mol.endAtoms();++atomItI){
        atomCodes.push_back(getAtomCode(*atomItI));
      }

      boost::dynamic_bitset<> *fromAtomsBV=0;
      if(fromAtoms){
        fromAtomsBV = new boost::dynamic_bitset<>(mol.getNumAtoms());
        BOOST_FOREACH(boost::uint32_t fAt,*fromAtoms){
          fromAtomsBV->set(fAt);
        }
      }
      boost::dynamic_bitset<> *ignoreAtomsBV=0;
      if(ignoreAtoms){
        ignoreAtomsBV = new boost::dynamic_bitset<>(mol.getNumAtoms());
        BOOST_FOREACH(boost::uint32_t fAt,*ignoreAtoms){
          ignoreAtomsBV->set(fAt);
        }
      }

      PATH_LIST paths=findAllPathsOfLengthN(mol,targetSize,false);
      for(PATH_LIST::const_iterator pathIt=paths.begin();
          pathIt!=paths.end();++pathIt){
        bool keepIt=true;
        if(fromAtomsBV){
          keepIt=false;
        }
        std::vector<boost::uint32_t> pathCodes;
        const PATH_TYPE &path=*pathIt;
        if(fromAtomsBV){
          if(fromAtomsBV->test(static_cast<boost::uint32_t>(path.front())) ||
             fromAtomsBV->test(static_cast<boost::uint32_t>(path.back()))){
            keepIt=true;
          }
        }
        if(keepIt && ignoreAtomsBV){
          BOOST_FOREACH(int pElem,path){
            if(ignoreAtomsBV->test(pElem)){
              keepIt=false;
              break;
            }
          }
        }
        if(keepIt){
          for(unsigned int i=0;i<targetSize;++i){
            unsigned int code=atomCodes[path[i]];
            // subtract off the branching number:
            if(i==0 || i==targetSize-1){
              code-=1;
            } else {
              code-=2;
            }
            pathCodes.push_back(code);
          }
          boost::int64_t code=getTopologicalTorsionCode(pathCodes);
          updateElement(*res,code);
        }
      }
      if(fromAtomsBV) delete fromAtomsBV;
      if(ignoreAtomsBV) delete ignoreAtomsBV;

      return res;
    }

    namespace {
      template <typename T>
      void TorsionFpCalc(T *res,
                         const ROMol &mol,                    
                         unsigned int nBits,
                         unsigned int targetSize,
                         const std::vector<boost::uint32_t> *fromAtoms,
                         const std::vector<boost::uint32_t> *ignoreAtoms){

        std::vector<boost::uint32_t> atomCodes;
        for(ROMol::ConstAtomIterator atomItI=mol.beginAtoms();
            atomItI!=mol.endAtoms();++atomItI){
          atomCodes.push_back(getAtomCode(*atomItI));
        }

        boost::dynamic_bitset<> *fromAtomsBV=0;
        if(fromAtoms){
          fromAtomsBV = new boost::dynamic_bitset<>(mol.getNumAtoms());
          BOOST_FOREACH(boost::uint32_t fAt,*fromAtoms){
            fromAtomsBV->set(fAt);
          }
        }
        boost::dynamic_bitset<> *ignoreAtomsBV=0;
        if(ignoreAtoms){
          ignoreAtomsBV = new boost::dynamic_bitset<>(mol.getNumAtoms());
          BOOST_FOREACH(boost::uint32_t fAt,*ignoreAtoms){
            ignoreAtomsBV->set(fAt);
          }
        }
        
        PATH_LIST paths=findAllPathsOfLengthN(mol,targetSize,false);
        for(PATH_LIST::const_iterator pathIt=paths.begin();
            pathIt!=paths.end();++pathIt){
          bool keepIt=true;
          if(fromAtomsBV){
            keepIt=false;
          }
          const PATH_TYPE &path=*pathIt;
          if(fromAtomsBV){
            if(fromAtomsBV->test(static_cast<boost::uint32_t>(path.front())) ||
               fromAtomsBV->test(static_cast<boost::uint32_t>(path.back()))){
              keepIt=true;
            }
          }
          if(keepIt && ignoreAtomsBV){
            BOOST_FOREACH(int pElem,path){
              if(ignoreAtomsBV->test(pElem)){
                keepIt=false;
                break;
              }
            }
          }
          if(keepIt){
            std::vector<boost::uint32_t> pathCodes(targetSize);
            for(unsigned int i=0;i<targetSize;++i){
              unsigned int code=atomCodes[path[i]];
              // subtract off the branching number:
              if(i==0 || i==targetSize-1){
                code-=1;
              } else {
                code-=2;
              }
              pathCodes[i]=code;
            }
            size_t bit=getTopologicalTorsionHash(pathCodes);
            updateElement(*res,bit%nBits);
          }
        }
        if(fromAtomsBV) delete fromAtomsBV;
        if(ignoreAtomsBV) delete ignoreAtomsBV;
      }
    } // end of local namespace
    SparseIntVect<boost::int64_t> *
    getHashedTopologicalTorsionFingerprint(const ROMol &mol,
                                           unsigned int nBits,
                                           unsigned int targetSize,
                                           const std::vector<boost::uint32_t> *fromAtoms,
                                           const std::vector<boost::uint32_t> *ignoreAtoms){

      SparseIntVect<boost::int64_t> *res=new SparseIntVect<boost::int64_t>(nBits);
      TorsionFpCalc(res,mol,nBits,targetSize,fromAtoms,ignoreAtoms);
      return res;
    }

    ExplicitBitVect *
    getHashedTopologicalTorsionFingerprintAsBitVect(const ROMol &mol,
                                                    unsigned int nBits,
                                                    unsigned int targetSize,
                                                    const std::vector<boost::uint32_t> *fromAtoms,
                                                    const std::vector<boost::uint32_t> *ignoreAtoms,
                                                    unsigned int nBitsPerEntry){
      static int bounds[4] = {1,2,4,8};
      SparseIntVect<boost::int64_t> *sres=new SparseIntVect<boost::int64_t>(nBits);
      TorsionFpCalc(sres,mol,nBits,targetSize,fromAtoms,ignoreAtoms);
      ExplicitBitVect *res=new ExplicitBitVect(nBits*nBitsPerEntry);

      if(nBitsPerEntry!=4){
        BOOST_FOREACH(SparseIntVect<boost::int64_t>::StorageType::value_type val,sres->getNonzeroElements()){
          for(unsigned int i=0;i<nBitsPerEntry;++i){
            if(val.second>i) res->setBit(val.first*nBitsPerEntry+i);
          }        
        }
      } else {
        BOOST_FOREACH(SparseIntVect<boost::int64_t>::StorageType::value_type val,sres->getNonzeroElements()){
          //std::cerr<<" "<<val.first<<"("<<val.second<<"):";
          for(unsigned int i=0;i<nBitsPerEntry;++i){
            if(val.second>=bounds[i]){
              res->setBit(val.first*nBitsPerEntry+i);
              //std::cerr<<" "<<i;
            }
          }        
          //std::cerr<<std::endl;
        }
      }
      delete sres;
      return res;
    }
  } // end of namespace AtomPairs
} // end of namespace RDKit
