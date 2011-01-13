//
//  Copyright (C) 2001-2005 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// std bits
#include <iostream>

// RD bits
#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <SmilesParse/SmilesParse.h>
#include "SubstructMatch.h"
#include <RDGeneral/RDLog.h>

// vflib bits
#include <argedit.h>
#include <argraph.h>
#include <vf_sub_state.h>
#include <match.h>

using namespace RDKit;

int main(int argc,char *argv[])
{
  RDLog::InitLogs();
  if(argc < 3){
    BOOST_LOG(rdErrorLog)<< "USAGE: test1 molSmiles querySmiles" << std::endl;
    exit(-1);
  }
    
  std::string mSmi(argv[1]);
  std::string qSmi(argv[2]);
  Mol *mol = SmilesToMol(mSmi);
  Mol *query = SmilesToMol(qSmi);
    
#if 1
  MatchVectType matchV;
  if(SubstructMatch(mol,query,matchV)){
    BOOST_LOG(rdInfoLog)<< "Got a match: " << std::endl;
    MatchVectType::iterator i;
    for(i=matchV.begin();i!=matchV.end();i++){
      BOOST_LOG(rdInfoLog)<< "\t" << i->first << " -> " << i->second << std::endl;
    }
    
  } else {
    BOOST_LOG(rdInfoLog)<< "No match" << std::endl;
  }
#else
  std::vector< MatchVectType > matches;
  int n = SubstructMatch(mol,query,matches);
  if(n){
    BOOST_LOG(rdInfoLog)<< "Got " << n << " matches:" << std::endl;
    std::vector<MatchVectType>::iterator i;
    for(i=matches.begin();i!=matches.end();i++){
      MatchVectType::iterator j;
      for(j=i->begin();j!=i->end();j++){
	BOOST_LOG(rdInfoLog)<< "\t" << j->first << " -> " << j->second << std::endl;
      }
      BOOST_LOG(rdInfoLog)<< std::endl;
    }
    
  } else {
    BOOST_LOG(rdInfoLog)<< "No match" << std::endl;
  }

  

#endif
  return 0;
}

  
