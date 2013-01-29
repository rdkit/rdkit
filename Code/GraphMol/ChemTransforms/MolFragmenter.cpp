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

#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDBoost/Exceptions.h>
#include <GraphMol/RDKitBase.h>
#include <boost/dynamic_bitset.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <algorithm>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/StreamOps.h>
#include <sstream>

namespace RDKit{
  namespace MolFragmenter{
    void constructFragmenterAtomTypes(std::istream *inStream,std::map<unsigned int,ROMOL_SPTR> &defs,
                                      std::string comment){
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
        ROMol *p=SmartsToMol(tokens[1]);
        if(!p){
          BOOST_LOG(rdWarningLog)<<"cannot convert SMARTS "<<tokens[1]<<" to molecule at line "<<line<<std::endl;
          continue;
        }
        defs[idx]=p;
      }
    };
    void constructFragmenterAtomTypes(const std::string &str,std::map<unsigned int,ROMOL_SPTR> &defs,
                                      std::string comment){
      std::stringstream istr(str);
      constructFragmenterAtomTypes(&istr,defs,comment);
    };
    void constructBRICSAtomTypes(std::map<unsigned int,ROMOL_SPTR> &defs){
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
4 [C;!D1;!$(C=*)]-;!@[#6]\n\
5 [N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]\n\
6 [C;D3;!R](=O)-;!@[#0,#6,#7,#8]\n\
7 [C;D2,D3]-[#6]\n\
8 [C;!R;!D1;!$(C!-*)]\n\
9 [n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]\n\
10 [N;R;$(N(@C(=O))@[C,N,O,S])]\n\
11 [S;D2](-;!@[#0,#6])\n\
12 [S;D4]([#6,#0])(=O)(=O)\n\
13 [C;$(C(-;@[C,N,O,S])-;@[N,O,S])]\n\
14 [c;$(c(:[c,n,o,s]):[n,o,s])]\n\
15 [C;$(C(-;@C)-;@C)]\n\
16 [c;$(c(:c):c)]"
        constructFragmenterAtomTypes(BRICSdefs,defs);
    }
    void constructFragmenterBondTypes(std::istream *istr,
                                      const std::map<unsigned int,ROMOL_SPTR> &atomTypes,
                                      std::vector<FragmenterBondType> &defs,
                                      std::string comment){
    }
    void constructFragmenterBondTypes(const std::string &str,
                                      const std::map<unsigned int,ROMOL_SPTR> &atomTypes,
                                      std::vector<FragmenterBondType> &defs,
                                      std::string comment){
      std::stringstream istr(str);
      constructFragmenterBondTypes(&istr,atomTypes,defs,comment);
    }
    void constructBRICSBondTypes(std::vector<FragmenterBondType> &defs){
      const std::string BRICSdefs=
"// L1
1 3 -
1 5 -
1 10 -
// L3 
3 4 -
3 13 -
3 14 -
3 15 -
3 16 -
// L4
4 5 -
4 11 -
// L5
5 12 -
5 14 -
5 16 -
5 13 -
5 15 -
// L6
6 13 -
6 14 -
6 15 -
6 16 -
// L7
7 7 =
// L8
8 9 -
8 10 -
8 13 -
8 14 -
8 15 -
8 16 -
// L9
9 13 - // not in original paper
9 14 - // not in original paper
9 15 -
9 16 -
// L10
10 13 -
10 14 -
10 15 -
10 16 -
// L11
11 13 -
11 14 -
11 15 -
11 16 -
// L12
// none left
// L13
13 14 -
13 15 -
13 16 -
// L14
14 14 - // not in original paper
14 15 -
14 16 -
// L15
15 16 -
// L16
16 16 - // not in original paper
";
      std::map<unsigned int,ROMOL_SPTR> atTypes;
      constructBRICSAtomTypes(atTypes);
      constructFragmenterBondTypes(BRICSdefs,atTypes,defs);
    }



    ROMol *fragmentOnBonds(const ROMol &mol,const std::vector<unsigned int> &bondIndices,
                           bool addDummies,
                           const std::vector< std::pair<unsigned int,unsigned int> > *dummyLabels){
      PRECONDITION( ( !dummyLabels || dummyLabels->size() == bondIndices.size() ), "bad dummy label vector");
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
          res->addBond(eidx,at1->getIdx(),bond->getBondType());
          unsigned int idx2=res->addAtom(at2,false,true);
          res->addBond(bidx,at2->getIdx(),bond->getBondType());

          for(ROMol::ConformerIterator confIt=res->beginConformers();
              confIt!=res->endConformers();++confIt){
            Conformer *conf=(*confIt).get();
            conf->setAtomPos(idx1,conf->getAtomPos(bidx));
            conf->setAtomPos(idx2,conf->getAtomPos(eidx));
          }
        }
      }
      return static_cast<ROMol *>(res);
    }
  } // end of namespace MolFragmenter
} // end of namespace RDKit
