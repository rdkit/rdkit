//  Created by Guillaume GODIN
//  Copyright (C) 2012-2018 Greg Landrum
//   @@ All Rights Reserved @@
//
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <RDGeneral/RDLog.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <boost/tokenizer.hpp>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/Descriptors/BagOfBonds.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>
#include <Eigen/Dense>

using namespace Eigen;

struct BagMatrix {
  VectorXd UCM;
  std::vector<std::string > BagTag;
};

std::vector<std::string> tokenize(const std::string &s) {
        boost::char_separator<char> sep(", \n\r\t");
        boost::tokenizer<boost::char_separator<char>> tok(s, sep);
        std::vector<std::string> tokens;
        std::copy(tok.begin(), tok.end(), std::back_inserter<std::vector<std::string> >(tokens));
  return tokens;
}

void testBagofBonds(){
    std::cout << "===================== Testing BagofBonds =======================\n";
 
    std::vector<std::string> strVec ={"CCCCCC",     "CCC(C)CC", "CC(C)CCCN", "OCC(C)C(C)CN",
                           "CC(C)(C)CC", "CCCCCO",   "CCC(O)CC", "CC(O)(C)CCOC",
                           "c1ccccc1O",  "CCCl",     "CCBr",     "CCI",  "CCCBr",
                           "OC(=O)c1ccncc1C(=O)O"};

    std::map<std::string, unsigned int> Global = RDKit::Descriptors::BagOfBondsMap(strVec);
  
    std::string pathName = getenv("RDBASE");

    std::string fName =
         pathName + "/Code/GraphMol/Descriptors/test_data/BoBMat.out";

    std::ifstream instrm(fName.c_str());
    std::string line;
    std::vector<std::string> tokens; 

    // need to size of the dictionary to resize the data output
    unsigned int resize = 0;
     std::cout << "===================== vis BoB Dictonary ========================\n";


    for (auto const& it :  Global) {
        std::getline(instrm, line);
        tokens = tokenize(line);
        TEST_ASSERT(tokens[0] == it.first);
        TEST_ASSERT(std::stoi(tokens[1]) == it.second);
        resize +=it.second;
    }
    std::cout << "BoB Mapping Done\n";

    // BoB tests
    std::string fNamebob =
         pathName + "/Code/GraphMol/Descriptors/test_data/BoB.out";

    std::ifstream instrmbob(fNamebob.c_str());
    std::string lineBoB;
    std::vector<std::string> tokensBoB; 

    for (auto const& s : strVec) {

          std::getline(instrmbob, lineBoB);
          tokensBoB = tokenize(lineBoB);

          TEST_ASSERT(tokensBoB[0] == s);

          //std::cout << s << "\n";

          RDKit::ROMOL_SPTR mol1( RDKit::SmilesToMol( s ) );
          RDKit::ROMOL_SPTR mol( RDKit::MolOps::addHs( *mol1 ) );

          std::vector<double>  res;

          RDKit::DGeomHelpers::EmbedParameters params( RDKit::DGeomHelpers::ETKDGv2 );
          params.randomSeed = 0xf00d;
          RDKit::DGeomHelpers::EmbedMolecule( *mol , params );
          //RDKit::UFF::UFFOptimizeMolecule(*mol);
          //std::cout << "===================== make BoB ========================\n";
          
          RDKit::Descriptors::BagOfBondsVector(*mol, res, -1, 1,  Global);

          //std::cout << "===================== vis BoB ========================\n";

          for (unsigned int i=0; i<res.size(); i++) {

              //std::cout << res[i] << ",";
              TEST_ASSERT(std::fabs(std::stof(tokensBoB[i+1]) - res[i])< 0.01);
          }
          //std::cout << "\n";
          res.clear();
          res.resize(resize);

    }
} 


int main() {
  RDLog::InitLogs();
  testBagofBonds();

}