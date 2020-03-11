//  Created by Guillaume GODIN
//  Copyright (C) 2020 Greg Landrum
//   @@ All Rights Reserved @@
//
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <RDGeneral/RDLog.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <GraphMol/GraphMol.h>
#include <boost/tokenizer.hpp>
#include <GraphMol/Descriptors/Augmentation.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;

void testAugmentation1(){
    std::cout << "===================== Testing Augmentation 1 =======================\n";
 
    auto mol = "C1CC=C1[NH3+]"_smiles; 
    unsigned int naug = 3;

    std::vector<std::string > exp {"C1C(=CC1)[NH3+]","C1CCC=1[NH3+]","[NH3+]C1=CCC1","[NH3+]C1CCC=1"};

    std::vector<std::string > res;

    RDKit::Descriptors::AugmentationVect(*mol, res, naug);
    unsigned int i =0;
    for ( const auto &v : res ) {
      TEST_ASSERT(v == exp[i]);
      //std::cout << v << ":" <<  exp[i] << "\n";
      i+=1;
    }
}

void testAugmentation2(){
    std::cout << "===================== Testing Augmentation 2  =======================\n";
 
    auto mol = "C1CC=C1[NH3+]"_smiles; 
    unsigned int naug = 10;

    std::vector<std::string > exp {"C1=C(CC1)[NH3+]","C1=C([NH3+])CC1","C1C=C(C1)[NH3+]",
    "C1CC([NH3+])=C1","C1CC=C1[NH3+]","C1CCC=1[NH3+]","[NH3+]C1=CCC1","[NH3+]C1CCC=1"};

    std::vector<std::string > res;

    RDKit::Descriptors::AugmentationVect(*mol, res, naug, false);
    unsigned int i =0;
    for ( const auto &v : res ) {
        TEST_ASSERT(v == exp[i]);
        //std::cout << v << ":" <<  exp[i] << "\n";
        i+=1;
    }
}

int main() {
  RDLog::InitLogs();
  testAugmentation1();
  testAugmentation2();
}