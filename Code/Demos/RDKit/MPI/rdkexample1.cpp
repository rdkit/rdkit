// $Id$
//
//  Copyright (C) 2009 Greg Landrum
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/RDLog.h>
#include <vector>
#include <algorithm>
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <iostream>
#include <cstdlib>
#include <string>

namespace mpi = boost::mpi;


int main(int argc, char* argv[])
{
  mpi::environment env(argc, argv);
  mpi::communicator world;

  // construct the data:
  std::vector<std::string> data;
  if (world.rank() == 0) {
    for(unsigned int i=0;i<100;++i){
      std::string txt(i+1,'C');
      data.push_back(txt);
    }
  }

  // broadcast it:
  broadcast(world,data,0);

  // process it:
  std::vector<unsigned int> res;
  std::vector<std::vector<unsigned int> > allRes;
  // start by finding our chunk:
  unsigned int nProcs=world.size(); 
  unsigned int chunkSize=data.size() / nProcs;
  unsigned int extraBits=data.size() % nProcs;

  // handle extra bits on the root node:
  if( world.rank() == 0 ){
    for(unsigned int i=0;i<extraBits;++i){
      const std::string &elem=data[i];
      res.push_back(elem.length());
    }
  }
  
  unsigned int pos=extraBits+world.rank()*chunkSize;
  for(unsigned int i=0;i<chunkSize;++i){
    const std::string &elem=data[pos++];
    res.push_back(elem.length());
  }

  if( world.rank() == 0 ){
    gather(world,res,allRes,0);
  } else {
    gather(world,res,0);
  }

  // reporting:
  if(world.rank()==0){
    for(unsigned int i=0;i<static_cast<unsigned int>(world.size());++i){
      std::cout<<"results from process "<<i<<": ";
      std::copy(allRes[i].begin(),allRes[i].end(),std::ostream_iterator<int,char>(std::cout, " "));
      std::cout<<std::endl;
    }

  }
  return 0;
} 

