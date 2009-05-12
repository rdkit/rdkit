// $Id$
//
//  Copyright (C) 2009 Greg Landrum
//   All Rights Reserved
//

#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/RDLog.h>
#include <vector>
#include <algorithm>
#include <boost/mpi.hpp>
#include <iostream>
#include <cstdlib>
#include <string>

namespace mpi = boost::mpi;


int main(int argc, char* argv[])
{
  mpi::environment env(argc, argv);
  mpi::communicator world;

  std::srand(time(0) + world.rank());
  std::string txt(std::rand()%10+1,'C');

  RDKit::ROMol *m=RDKit::SmilesToMol(txt);

  if (world.rank() == 0) {
    std::vector<unsigned int> all_vals;
    gather(world, m->getNumAtoms(), all_vals, 0);
    for (int proc = 0; proc < world.size(); ++proc)
      std::cout << "Process #" << proc << " thought of " 
                << all_vals[proc] << std::endl;
  } else {
    gather(world, m->getNumAtoms(), 0);
  }
  delete m;
  return 0;
} 

