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
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/RDLog.h>
#include <vector>
#include <algorithm>
#include <boost/mpi.hpp>
#include <boost/foreach.hpp>
#include <iostream>
#include <cstdlib>
#include <string>

namespace mpi = boost::mpi;

void broadcastMols(mpi::communicator &world,std::vector<RDKit::ROMOL_SPTR> &mols){
  std::vector<std::string> data;
  if (world.rank() == 0) {
    data.reserve(mols.size());
    BOOST_FOREACH( const RDKit::ROMOL_SPTR &ptr, mols ) {
      std::string pickle;
      RDKit::MolPickler::pickleMol(*ptr,pickle);
      data.push_back(pickle);
    }
  }
  broadcast(world,data,0);
  if (world.rank() != 0) {
    mols.reserve(data.size());
    BOOST_FOREACH( const std::string &pickle, data ) {
      RDKit::ROMol *mol=new RDKit::ROMol;
      RDKit::MolPickler::molFromPickle(pickle,mol);
      mols.push_back(RDKit::ROMOL_SPTR(mol));
    }
  }
}

int main(int argc, char* argv[])
{
  mpi::environment env(argc, argv);
  mpi::communicator world;

  // construct the data:
  std::vector<RDKit::ROMOL_SPTR> data;
  if (world.rank() == 0) {
    for(unsigned int i=0;i<100;++i){
      std::string txt(i+1,'C');
      RDKit::ROMol *m=RDKit::SmilesToMol(txt);
      data.push_back(RDKit::ROMOL_SPTR(m));
    }
  }

  // broadcast it:
  broadcastMols(world,data);

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
      res.push_back(data[i]->getNumAtoms());
    }
  }
  
  unsigned int pos=extraBits+world.rank()*chunkSize;
  for(unsigned int i=0;i<chunkSize;++i){
    res.push_back(data[pos++]->getNumAtoms());
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

