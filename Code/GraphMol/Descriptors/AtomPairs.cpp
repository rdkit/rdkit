// $Id$
//
//  Copyright (C) 2007-2008 Greg Landrum
//
//   @@ All Rights Reserved  @@
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/AtomPairs.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <DataStructs/SparseIntVect.h>
#include <boost/functional/hash.hpp>

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

      void setAtomPairBit(boost::uint32_t i, boost::uint32_t j,boost::uint32_t nAtoms,
                          const std::vector<boost::uint32_t> &atomCodes,
                          const double *dm,SparseIntVect<boost::int32_t> *bv,
                          unsigned int minLength,unsigned int maxLength){
        unsigned int dist=static_cast<unsigned int>(floor(dm[i*nAtoms+j]));
        if(dist>=minLength && dist<=maxLength){
          boost::uint32_t bitId=getAtomPairCode(atomCodes[i],atomCodes[j],dist);
          bv->setVal(bitId,(*bv)[bitId]+1);
        }
      }

      SparseIntVect<boost::int32_t> *getAtomPairFingerprint(const ROMol &mol,
                                                            const std::vector<boost::uint32_t> *fromAtoms){
        return getAtomPairFingerprint(mol,1,maxPathLen-1,fromAtoms);
      };

      SparseIntVect<boost::int32_t> *
      getAtomPairFingerprint(const ROMol &mol,unsigned int minLength,unsigned int maxLength,
                             const std::vector<boost::uint32_t> *fromAtoms){
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
          if(!fromAtoms){
            for(ROMol::ConstAtomIterator atomItJ=atomItI+1;
                atomItJ!=mol.endAtoms();++atomItJ){
              unsigned int j=(*atomItJ)->getIdx();
              setAtomPairBit(i,j,nAtoms,atomCodes,dm,res,minLength,maxLength);
            }
          } else {
            for(std::vector<boost::uint32_t>::const_iterator bvIt=fromAtoms->begin();
                bvIt!=fromAtoms->end();++bvIt){
              if(*bvIt!=i){
                setAtomPairBit(i,*bvIt,nAtoms,atomCodes,dm,res,minLength,maxLength);
              }
            }
          }
        }
        return res;
      }

      SparseIntVect<boost::int32_t> *
      getHashedAtomPairFingerprint(const ROMol &mol,unsigned int nBits,
                                   unsigned int minLength,unsigned int maxLength){
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
          for(ROMol::ConstAtomIterator atomItJ=atomItI+1;
              atomItJ!=mol.endAtoms();++atomItJ){
            unsigned int j=(*atomItJ)->getIdx();
            unsigned int dist=static_cast<unsigned int>(floor(dm[i*nAtoms+j]));
            if(dist>=minLength && dist<=maxLength){
              std::size_t bit=0;
              boost::hash_combine(bit,std::min(atomCodes[i],atomCodes[j]));
              boost::hash_combine(bit,dist);
              boost::hash_combine(bit,std::max(atomCodes[i],atomCodes[j]));
              res->setVal(bit%nBits,(*res)[bit%nBits]+1);
            }
          }
        }
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

        size_t res=0;
        if(reverseIt){
          for(unsigned int i=0;i<pathCodes.size();++i){
            boost::hash_combine(res,pathCodes[pathCodes.size()-i-1]);
          }
        }else{
          for(unsigned int i=0;i<pathCodes.size();++i){
            boost::hash_combine(res,pathCodes[i]);
          }
        }
        return res;
      }

      
      SparseIntVect<boost::int64_t> *
      getTopologicalTorsionFingerprint(const ROMol &mol,unsigned int targetSize,
                                       const std::vector<boost::uint32_t> *fromAtoms){
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

        PATH_LIST paths=findAllPathsOfLengthN(mol,targetSize,false);
        for(PATH_LIST::const_iterator pathIt=paths.begin();
            pathIt!=paths.end();++pathIt){
          bool keepIt=true;
          if(fromAtoms){
            keepIt=false;
          }
          std::vector<boost::uint32_t> pathCodes;
          const PATH_TYPE &path=*pathIt;
          if(fromAtoms){
            if(std::find(fromAtoms->begin(),fromAtoms->end(),
                         static_cast<boost::uint32_t>(path.front()))!=fromAtoms->end() ||
               std::find(fromAtoms->begin(),fromAtoms->end(),
                         static_cast<boost::uint32_t>(path.back()))!=fromAtoms->end()
               ){
              keepIt=true;
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
            res->setVal(code,res->getVal(code)+1);
          }
        }
        return res;
      }

      SparseIntVect<boost::int64_t> *
      getHashedTopologicalTorsionFingerprint(const ROMol &mol,
                                             unsigned int nBits,
                                             unsigned int targetSize,
                                             const std::vector<boost::uint32_t> *fromAtoms){
        SparseIntVect<boost::int64_t> *res=new SparseIntVect<boost::int64_t>(nBits);

        std::vector<boost::uint32_t> atomCodes;
        for(ROMol::ConstAtomIterator atomItI=mol.beginAtoms();
            atomItI!=mol.endAtoms();++atomItI){
          atomCodes.push_back(getAtomCode(*atomItI));
        }

        PATH_LIST paths=findAllPathsOfLengthN(mol,targetSize,false);
        for(PATH_LIST::const_iterator pathIt=paths.begin();
            pathIt!=paths.end();++pathIt){
          bool keepIt=true;
          if(fromAtoms){
            keepIt=false;
          }
          const PATH_TYPE &path=*pathIt;
          if(fromAtoms){
            if(std::find(fromAtoms->begin(),fromAtoms->end(),
                         static_cast<boost::uint32_t>(path.front()))!=fromAtoms->end() ||
               std::find(fromAtoms->begin(),fromAtoms->end(),
                         static_cast<boost::uint32_t>(path.back()))!=fromAtoms->end()
               ){
              keepIt=true;
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
            res->setVal(bit%nBits,res->getVal(bit%nBits)+1);
          }
        }
        return res;
      }

    } // end of namespace AtomPairs
  } // end of namespace Descriptors
} // end of namespace RDKit
