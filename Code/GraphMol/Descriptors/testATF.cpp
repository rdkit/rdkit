//  Created by Guillaume GODIN
//  Copyright (C) 2012-2020 Greg Landrum
//   @@ All Rights Reserved @@
//
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/RDLog.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>

#include <chrono>  // for high_resolution_clock

#include <GraphMol/Descriptors/AtomFeat.h>

void test1() {
 
  
    RDKit::ROMOL_SPTR m( RDKit::SmilesToMol( "CO") );    
    TEST_ASSERT(m);
    
    std::vector<double> res;

    int atomid = 0;
    RDKit::Descriptors::AtomFeat(*m, res, atomid);

    std::vector <double > exp{ 0. , 1. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ,
       0. , 0. , 1. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 1. , 0. , 0. ,
       0. , 0. , 0. , 1. , 0. , 0. , 0. , 0. , 1. , 0. , 0. , 0. , 0. ,
       0. , 0. , 0. , 0. , 0. , 0. , 0. , 1. , 0. , 0.5};

    // 49 features
    TEST_ASSERT(res.size() == 49);

    for (std::size_t i = 0; i < res.size() ; i++) {
     	TEST_ASSERT(fabs( res[i]-exp[i])< 0.001);
      std::cout << res[i] << "," ;
    }
    std::cout << "\nnum features: " << res.size()  << "\n" ; 

    std::cout << "DONE\n"; 
    
}

int main() {
  RDLog::InitLogs();
  test1();
}