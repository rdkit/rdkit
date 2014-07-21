// $Id$
//
//  Copyright (C) 2013 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "MolFragmenter.h"
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDBoost/Exceptions.h>
#include <GraphMol/RDKitBase.h>
#include <boost/dynamic_bitset.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/cstdint.hpp>
#include <vector>
#include <algorithm>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/StreamOps.h>
#include <boost/flyweight.hpp>
#include <boost/flyweight/no_tracking.hpp>
#include <boost/functional/hash.hpp>
#include <sstream>

namespace RDKit{
  namespace MolFragmenter{
    std::size_t hash_value(const FragmenterBondType &fbt){
      size_t res=boost::hash<int>()((int)fbt.bondType);
      boost::hash_combine(res,fbt.atom1Label);
      boost::hash_combine(res,fbt.atom2Label);
      return res;
    }
    bool operator==(const FragmenterBondType &v1,const FragmenterBondType &v2){
      return (v1.atom1Label==v2.atom1Label)&&
        (v1.atom2Label==v2.atom2Label)&&
        (v1.bondType==v2.bondType);
    }
    void constructFragmenterAtomTypes(std::istream *inStream,std::map<unsigned int,std::string> &defs,
                                      std::string comment,bool validate,
                                      std::map<unsigned int,ROMOL_SPTR> *environs){
      PRECONDITION(inStream,"no stream");
      defs.clear();
      unsigned int line=0;
      while(!inStream->eof()){
        ++line;
        std::string tempStr=getLine(inStream);
        if(tempStr=="" || tempStr.find(comment)==0 ) continue;
        std::vector<std::string> tokens;
        boost::split(tokens,tempStr,boost::is_any_of(" \t"),boost::token_compress_on);
        if(tokens.size()<2){
          BOOST_LOG(rdWarningLog)<<"line "<<line<<" is too short"<<std::endl;
          continue;
        }
        unsigned int idx=boost::lexical_cast<unsigned int>(tokens[0]);
        if(defs.find(idx)!=defs.end()){
          BOOST_LOG(rdWarningLog)<<"definition #"<<idx<<" encountered more than once. Using the first occurance."<<std::endl;
          continue;
        }
        if(validate || environs){
          ROMol *p=SmartsToMol(tokens[1]);
          if(!p){
            BOOST_LOG(rdWarningLog)<<"cannot convert SMARTS "<<tokens[1]<<" to molecule at line "<<line<<std::endl;
            continue;
          }
          if(!environs){
            delete p;
          } else {
            (*environs)[idx] = ROMOL_SPTR(p);
          }
        }
        defs[idx]=tokens[1];
      }
    }
    void constructFragmenterAtomTypes(const std::string &str,std::map<unsigned int,std::string> &defs,
                                      std::string comment,bool validate,
                                      std::map<unsigned int,ROMOL_SPTR> *environs){
      std::stringstream istr(str);
      constructFragmenterAtomTypes(&istr,defs,comment,validate,environs);
    }
    void constructBRICSAtomTypes(std::map<unsigned int,std::string> &defs,
                                 std::map<unsigned int,ROMOL_SPTR> *environs){
      /* 
         After some discussion, the L2 definitions ("N.pl3" in the original
         paper) have been removed and incorporated into a (almost) general
         purpose amine definition in L5 ("N.sp3" in the paper).
        
         The problem is one of consistency.
            Based on the original definitions you should get the following
            fragmentations:
              C1CCCCC1NC(=O)C -> C1CCCCC1N[2*].[1*]C(=O)C
              c1ccccc1NC(=O)C -> c1ccccc1[16*].[2*]N[2*].[1*]C(=O)C
            This difference just didn't make sense to us. By switching to
            the unified definition we end up with:
              C1CCCCC1NC(=O)C -> C1CCCCC1[15*].[5*]N[5*].[1*]C(=O)C
              c1ccccc1NC(=O)C -> c1ccccc1[16*].[5*]N[5*].[1*]C(=O)C
      */
      const std::string BRICSdefs="1 [C;D3]([#0,#6,#7,#8])(=O)\n\
3 [O;D2]-;!@[#0,#6,#1]\n\
5 [N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]\n\
9 [n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]\n\
10 [N;R;$(N(@C(=O))@[C,N,O,S])]\n\
11 [S;D2](-;!@[#0,#6])\n\
12 [S;D4]([#6,#0])(=O)(=O)\n\
6 [C;D3;!R](=O)-;!@[#0,#6,#7,#8]\n\
13 [C;$(C(-;@[C,N,O,S])-;@[N,O,S])]\n\
14 [c;$(c(:[c,n,o,s]):[n,o,s])]\n\
15 [C;$(C(-;@C)-;@C)]\n\
4 [C;!D1;!$(C=*)]-;!@[#6]\n\
7 [C;D2,D3]-[#6]\n\
8 [C;!R;!D1;!$(C!-*)]\n\
16 [c;$(c(:c):c)]";
      constructFragmenterAtomTypes(BRICSdefs,defs,"//",true,environs);
    }


    void constructFragmenterBondTypes(std::istream *inStream,
                                      const std::map<unsigned int,std::string> &atomTypes,
                                      std::vector<FragmenterBondType> &defs,
                                      std::string comment,bool validate,bool labelByConnector){
      PRECONDITION(inStream,"no stream");
      defs.clear();
      defs.resize(0);

      unsigned int line=0;
      while(!inStream->eof()){
        ++line;
        std::string tempStr=getLine(inStream);
        if(tempStr=="" || tempStr.find(comment)==0 ) continue;
        std::vector<std::string> tokens;
        boost::split(tokens,tempStr,boost::is_any_of(" \t"),boost::token_compress_on);
        if(tokens.size()<3){
          BOOST_LOG(rdWarningLog)<<"line "<<line<<" is too short"<<std::endl;
          continue;
        }
        unsigned int idx1=boost::lexical_cast<unsigned int>(tokens[0]);
        if(atomTypes.find(idx1)==atomTypes.end()){
          BOOST_LOG(rdWarningLog)<<"atom type #"<<idx1<<" not recognized."<<std::endl;
          continue;
        }
        unsigned int idx2=boost::lexical_cast<unsigned int>(tokens[1]);
        if(atomTypes.find(idx2)==atomTypes.end()){
          BOOST_LOG(rdWarningLog)<<"atom type #"<<idx2<<" not recognized."<<std::endl;
          continue;
        }
        std::string sma1=atomTypes.find(idx1)->second;
        std::string sma2=atomTypes.find(idx2)->second;
        std::string smarts="[$("+ sma1 +")]"+tokens[2]+"[$("+ sma2 +")]";
        ROMol *p=SmartsToMol(smarts);
        if(validate){
          if(!p){
            BOOST_LOG(rdWarningLog)<<"cannot convert SMARTS "<<smarts<<" to molecule at line "<<line<<std::endl;
            continue;
          }
        }
        FragmenterBondType fbt;
        fbt.atom1Type=idx1;
        fbt.atom2Type=idx2;
        if(labelByConnector){
          fbt.atom1Label=idx1;
          fbt.atom2Label=idx2;
        } else {
          fbt.atom1Label=idx2;
          fbt.atom2Label=idx1;
        }
        if(p){
          // for the purposes of replacing the bond, we'll use just the first
          // character to set the bond type (if we recognize it):
          switch(tokens[2][0]){
          case '-':
            fbt.bondType = Bond::SINGLE;break;
          case '=':
            fbt.bondType = Bond::DOUBLE;break;
          case '#':
            fbt.bondType = Bond::TRIPLE;break;
          case ':':
            fbt.bondType = Bond::AROMATIC;break;
          default:
            fbt.bondType = p->getBondWithIdx(0)->getBondType();
          }
          fbt.query = ROMOL_SPTR(p);
        } else {
          fbt.bondType=Bond::UNSPECIFIED;
          fbt.query=ROMOL_SPTR();
        }
        defs.push_back(fbt);
      }
    }

    void constructFragmenterBondTypes(const std::string &str,
                                      const std::map<unsigned int,std::string> &atomTypes,
                                      std::vector<FragmenterBondType> &defs,
                                      std::string comment,bool validate,
                                      bool labelByConnector){
      std::stringstream istr(str);
      constructFragmenterBondTypes(&istr,atomTypes,defs,comment,validate,labelByConnector);
    }
    void constructBRICSBondTypes(std::vector<FragmenterBondType> &defs){
      const std::string BRICSdefs=
"// L1\n\
1 3 -;!@\n\
1 5 -;!@\n\
1 10 -;!@\n\
// L3 \n\
3 4 -;!@\n\
3 13 -;!@\n\
3 14 -;!@\n\
3 15 -;!@\n\
3 16 -;!@\n\
// L4\n\
4 5 -;!@\n\
4 11 -;!@\n\
// L5\n\
5 12 -;!@\n\
5 14 -;!@\n\
5 16 -;!@\n\
5 13 -;!@\n\
5 15 -;!@\n\
// L6\n\
6 13 -;!@\n\
6 14 -;!@\n\
6 15 -;!@\n\
6 16 -;!@\n\
// L7\n\
7 7 =;!@\n\
// L8\n\
8 9 -;!@\n\
8 10 -;!@\n\
8 13 -;!@\n\
8 14 -;!@\n\
8 15 -;!@\n\
8 16 -;!@\n\
// L9\n\
9 13 -;!@ // not in original paper\n\
9 14 -;!@ // not in original paper\n\
9 15 -;!@\n\
9 16 -;!@\n\
// L10\n\
10 13 -;!@\n\
10 14 -;!@\n\
10 15 -;!@\n\
10 16 -;!@\n\
// L11\n\
11 13 -;!@\n\
11 14 -;!@\n\
11 15 -;!@\n\
11 16 -;!@\n\
// L12\n\
// none left\n\
// L13\n\
13 14 -;!@\n\
13 15 -;!@\n\
13 16 -;!@\n\
// L14\n\
14 14 -;!@ // not in original paper\n\
14 15 -;!@\n\
14 16 -;!@\n\
// L15\n\
15 16 -;!@\n\
// L16\n\
16 16 -;!@ // not in original paper";
      std::map<unsigned int,std::string> atTypes;
      constructBRICSAtomTypes(atTypes);
      constructFragmenterBondTypes(BRICSdefs,atTypes,defs,"//",true,false);
    }

    namespace {
      boost::uint64_t nextBitCombo(boost::uint64_t v){
        // code from: http://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
        boost::uint64_t t = (v | (v - 1)) + 1;  
        return t | ((((t & -t) / (v & -v)) >> 1) - 1);  
      }
    }
    
    void fragmentOnSomeBonds(const ROMol &mol,const std::vector<unsigned int> &bondIndices,
                             std::vector<ROMOL_SPTR> &resMols,
                             unsigned int maxToCut,
                             bool addDummies,
                             const std::vector< std::pair<unsigned int,unsigned int> > *dummyLabels,
                             const std::vector< Bond::BondType > *bondTypes){
      PRECONDITION( ( !dummyLabels || dummyLabels->size() == bondIndices.size() ), "bad dummyLabel vector");
      PRECONDITION( ( !bondTypes || bondTypes->size() == bondIndices.size() ), "bad bondType vector");
      if(bondIndices.size()>63) throw ValueErrorException("currently can only fragment on up to 63 bonds");
      boost::uint64_t state=(0x1<<maxToCut)-1;
      boost::uint64_t stop=0x1<<bondIndices.size();
      std::vector<unsigned int> fragmentHere(maxToCut);
      std::vector< std::pair<unsigned int,unsigned int> > *dummyLabelsHere=NULL;
      if(dummyLabels){
        dummyLabelsHere=new std::vector< std::pair<unsigned int,unsigned int> >(maxToCut);
      }
      std::vector< Bond::BondType > *bondTypesHere=NULL;
      if(bondTypes){
        bondTypesHere=new std::vector< Bond::BondType >(maxToCut);
      }
      while(state<stop){
        unsigned int nSeen=0;
        for(unsigned int i=0;i<bondIndices.size() && nSeen<maxToCut;++i){
          if(state&(0x1<<i)){
            fragmentHere[nSeen]=bondIndices[i];
            if(dummyLabelsHere) (*dummyLabelsHere)[nSeen]=(*dummyLabels)[i];
            if(bondTypesHere) (*bondTypesHere)[nSeen]=(*bondTypes)[i];
            ++nSeen;
          }
        }
        ROMol *nm=fragmentOnBonds(mol,fragmentHere,addDummies,dummyLabelsHere,bondTypesHere);
        resMols.push_back(ROMOL_SPTR(nm));
        
        state = nextBitCombo(state);
      }
      delete dummyLabelsHere;
      delete bondTypesHere;
    }

    ROMol *fragmentOnBonds(const ROMol &mol,const std::vector<unsigned int> &bondIndices,
                           bool addDummies,
                           const std::vector< std::pair<unsigned int,unsigned int> > *dummyLabels,
                           const std::vector< Bond::BondType > *bondTypes){
      PRECONDITION( ( !dummyLabels || dummyLabels->size() == bondIndices.size() ), "bad dummyLabel vector");
      PRECONDITION( ( !bondTypes || bondTypes->size() == bondIndices.size() ), "bad bondType vector");
      RWMol *res=new RWMol(mol);
      std::vector<Bond *> bondsToRemove;
      bondsToRemove.reserve(bondIndices.size());
      BOOST_FOREACH(unsigned int bondIdx,bondIndices){
        bondsToRemove.push_back(res->getBondWithIdx(bondIdx));
      }
      for(unsigned int i=0;i<bondsToRemove.size();++i){
        const Bond *bond=bondsToRemove[i];
        unsigned int bidx=bond->getBeginAtomIdx();
        unsigned int eidx=bond->getEndAtomIdx();
	    Bond::BondType bT=bond->getBondType();
        res->removeBond(bond->getBeginAtomIdx(),bond->getEndAtomIdx());
        if(addDummies){
          Atom *at1,*at2;
          at1 = new Atom(0);
          at2 = new Atom(0);
          if(dummyLabels){
            at1->setIsotope((*dummyLabels)[i].first);
            at2->setIsotope((*dummyLabels)[i].second);
          } else {
            at1->setIsotope(bidx);
            at2->setIsotope(eidx);
          }
          unsigned int idx1=res->addAtom(at1,false,true);
          if(bondTypes) bT=(*bondTypes)[i];
          res->addBond(eidx,at1->getIdx(),bT);
          unsigned int idx2=res->addAtom(at2,false,true);
          res->addBond(bidx,at2->getIdx(),bT);

          for(ROMol::ConformerIterator confIt=res->beginConformers();
              confIt!=res->endConformers();++confIt){
            Conformer *conf=(*confIt).get();
            conf->setAtomPos(idx1,conf->getAtomPos(bidx));
            conf->setAtomPos(idx2,conf->getAtomPos(eidx));
          }
        }
      }
      res->clearComputedProps();
      return static_cast<ROMol *>(res);
    }

    ROMol *fragmentOnBonds(const ROMol &mol,const std::vector<FragmenterBondType> &bondPatterns,
                           const std::map<unsigned int,ROMOL_SPTR> *atomEnvirons){
      std::vector<unsigned int> bondIndices;
      std::vector< std::pair<unsigned int,unsigned int> > dummyLabels;
      std::vector<Bond::BondType> bondTypes;

      std::map<unsigned int,bool> environsMatch;
      if(atomEnvirons){
        for(std::map<unsigned int,ROMOL_SPTR>::const_iterator iter=atomEnvirons->begin();
            iter!=atomEnvirons->end();++iter){
          MatchVectType mv;
          environsMatch[iter->first]=SubstructMatch(mol,*(iter->second),mv);
        }
      }

      
      boost::dynamic_bitset<> bondsUsed(mol.getNumBonds(),0);
      // the bond definitions are organized (more or less) general -> specific, so loop
      // over them backwards
      BOOST_REVERSE_FOREACH(const FragmenterBondType &fbt,bondPatterns){
        if(fbt.query->getNumAtoms()!=2 || fbt.query->getNumBonds()!=1){
            BOOST_LOG(rdErrorLog)<<"fragmentation queries must have 2 atoms and 1 bond"<<std::endl;
            continue;
        }
        if(atomEnvirons && (!environsMatch[fbt.atom1Type] || !environsMatch[fbt.atom2Type] ) )
          continue;
        //std::cerr<<"  >>> "<<fbt.atom1Label<<" "<<fbt.atom2Label<<std::endl;
        std::vector<MatchVectType> bondMatches;
        SubstructMatch(mol,*fbt.query.get(),bondMatches);
        BOOST_FOREACH(const MatchVectType &mv,bondMatches){
          const Bond *bond=mol.getBondBetweenAtoms(mv[0].second,mv[1].second);
          //std::cerr<<"          "<<bond->getIdx()<<std::endl;
          TEST_ASSERT(bond);
          if(bondsUsed[bond->getIdx()]){
            //BOOST_LOG(rdWarningLog)<<"bond #"<<bond->getIdx()<<" matched multiple times in decomposition. Later matches ignored."<<std::endl;
            continue;
          }
          bondsUsed.set(bond->getIdx());
          bondIndices.push_back(bond->getIdx());
          if(bond->getBeginAtomIdx()==static_cast<unsigned int>(mv[0].second)){
            dummyLabels.push_back(std::make_pair(fbt.atom1Label,fbt.atom2Label));
          } else {
            dummyLabels.push_back(std::make_pair(fbt.atom2Label,fbt.atom1Label));
          }
          bondTypes.push_back(fbt.bondType);
        }
      }
      return fragmentOnBonds(mol,bondIndices,true,&dummyLabels,&bondTypes);
    }

    boost::flyweight<std::vector<FragmenterBondType>,boost::flyweights::no_tracking> bondPatterns;
    boost::flyweight<std::map<unsigned int,ROMOL_SPTR>,boost::flyweights::no_tracking> atomEnvs;
    ROMol *fragmentOnBRICSBonds(const ROMol &mol){
      if(bondPatterns.get().size()==0){
        std::map<unsigned int,std::string> adefs;
        std::map<unsigned int,ROMOL_SPTR> aenvs;
        constructBRICSAtomTypes(adefs,&aenvs);
        atomEnvs=aenvs;
        std::vector<FragmenterBondType> tbondPatterns ;
        constructBRICSBondTypes(tbondPatterns);
        bondPatterns=tbondPatterns;
      }
      return fragmentOnBonds(mol,bondPatterns,&(atomEnvs.get()));
    }

  } // end of namespace MolFragmenter
} // end of namespace RDKit
