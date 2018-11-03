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
#include <GraphMol/Descriptors/CoulombMat.h>
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
  
    std::string pathName = "/Users/tgg/Github/rdkit"; //getenv("RDBASE");
    //std::cout << pathName << ",";

    std::string fName =
         pathName + "/Code/GraphMol/Descriptors/test_data/BoBMat.out";

    std::ifstream instrm(fName.c_str());
    std::string line;
    std::vector<std::string> tokens; 

    // need to size of the dictionary to resize the data output
    unsigned int resize = 0;
       // std::cout << "===================== vis BoB Dictonary ========================\n";

    int i = 0;
    for (auto const& i :  Global) {
        std::getline(instrm, line);
        tokens = tokenize(line);

        TEST_ASSERT(tokens[0] == i.first);
        TEST_ASSERT(std::stoi(tokens[1]) == i.second);

        //std::cout << tokens[0] << ":" << tokens[1] << "\n";
        //std::cout << i.first << ":" << i.second << ",";
        resize +=i.second;
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

          RDKit::ROMol *mol1;

          TEST_ASSERT(tokensBoB[0] == s);

          //std::cout << s << "\n";
          mol1 = RDKit::SmilesToMol(s);
          RDKit::ROMol *mol = RDKit::MolOps::addHs(*mol1);
          std::vector<double>  res;

          RDKit::DGeomHelpers::EmbedMolecule(*mol, 0, 1234);
          RDKit::UFF::UFFOptimizeMolecule(*mol);
          
          RDKit::Descriptors::BagOfBondsVector(*mol, res, -1, 1,  Global);

          //std::cout << "===================== vis BoB ========================\n";

          for (unsigned int i=0; i<res.size(); i++) {
              //std::cout << res[i]  << "," << tokensBoB[i+1] << ";";
              TEST_ASSERT(std::fabs(std::stof(tokensBoB[i+1]) - res[i])< 0.01);
          }
          //std::cout << "\n";
          res.clear();
          res.resize(resize);
          delete mol;
          delete mol1;

    }
    std::cout << "BoB Vectors Done\n";
} 

void testCoulombMat1(){
    std::cout << "===================== Testing CoulombMat GLobal =======================\n";
 
    std::string pathName = "/Users/tgg/Github/rdkit"; //getenv("RDBASE");

    // CM test 1
    std::string fNameCM =
         pathName + "/Code/GraphMol/Descriptors/test_data/CM1.out";

    std::ifstream instrmCM(fNameCM.c_str());
    std::string line;
    std::vector<std::string> tokens; 



    RDKit::ROMOL_SPTR mol1( RDKit::SmilesToMol( "NCCCCCO" ) );
    RDKit::ROMOL_SPTR mol( RDKit::MolOps::addHs( *mol1 ) );


    RDKit::DGeomHelpers::EmbedParameters params( RDKit::DGeomHelpers::ETKDG );
    params.randomSeed = 0xf00d;
    RDKit::DGeomHelpers::EmbedMolecule( *mol , params );
    std::vector<std::vector<double>> Mres;
    int nbmats= 10;
    int confId = -1;
    int padding = 23; // padding is the size of the result matrix so >= size max(atoms)
    int seed = 0xf00d;
    double rcut = 0; // not used if local is false

    // definition switch between the local & global than decay/reduced/sorted/eigen also alpha coefficient!
    //CoulombMat(const ROMol &mol, std::vector<std::vector<double>> &res, int confId, int nbmats,
    //    int seed, int padding, double rcut, bool local, bool decaying, bool reduced,  bool sorted,
    //    bool eigenval, int alpha);
    RDKit::Descriptors::CoulombMat(*mol, Mres, confId, nbmats, seed, padding, rcut, false, false, false, false, false, 2);

    //std::cout << "===================== vis BoB ========================\n";

    for ( const auto &v : Mres ) {

      std::getline(instrmCM, line);
      tokens = tokenize(line);

      int ti = 0;
      for ( double x : v )  {
          //std::cout << x << ",";
          TEST_ASSERT(std::fabs(std::stof(tokens[ti]) - x)< 0.001);
          ti++;
      } 
      //std::cout << "\n";
    }
    //std::cout << "\n";
    Mres.clear();
    Mres.resize(nbmats);

    
    std::cout << "CM test 1 mat Done\n";
}


void testCoulombMat2(){
    std::cout << "===================== Testing CoulombMat Local 1 =======================\n";
 
    std::string pathName = "/Users/tgg/Github/rdkit"; //getenv("RDBASE");

    ///////////////////////
    // TEST 1 LOCAL
    std::string fNameCM2 =
         pathName + "/Code/GraphMol/Descriptors/test_data/CM2.out";

    std::ifstream instrmCM2(fNameCM2.c_str());

    std::string line;
    std::vector<std::string> tokens; 



    RDKit::ROMOL_SPTR mol1( RDKit::SmilesToMol( "NCCCCCO" ) );
    RDKit::ROMOL_SPTR mol( RDKit::MolOps::addHs( *mol1 ) );


    RDKit::DGeomHelpers::EmbedParameters params( RDKit::DGeomHelpers::ETKDG );
    params.randomSeed = 0xf00d;
    RDKit::DGeomHelpers::EmbedMolecule( *mol , params );
    std::vector<std::vector<double>> Mres;
    int nbmats= 10;
    int confId = -1;
    int padding = 23; // padding is the size of the result matrix so >= size max(atoms)
    int seed = 0xf00d; 

    double rcut = 2.0; 
    // so local is ture 
    RDKit::Descriptors::CoulombMat(*mol, Mres, confId, nbmats, seed, padding, rcut, true, false, false, false, false, 1);

    //std::cout << "===================== vis BoB ========================\n";

    for ( const auto &v : Mres ) {

      std::getline(instrmCM2, line);
      tokens = tokenize(line);

      int ti = 0;
      for ( double x : v )  {
          //std::cout << x << ",";
          TEST_ASSERT(std::fabs(std::stof(tokens[ti]) - x)< 0.001);
          ti++;
      } 
      //std::cout << "\n";
    }
    //std::cout << "\n";
    Mres.clear();
    Mres.resize(nbmats);

    std::cout << "CM test 2 mat Done\n";
}


void testCoulombMat3(){
    std::cout << "===================== Testing CoulombMat Local 2 =======================\n";
 
    std::string pathName = "/Users/tgg/Github/rdkit"; //getenv("RDBASE");

    //////////////
    // TEST 2 LOCAL decaying
    std::string fNameCM3 =
         pathName + "/Code/GraphMol/Descriptors/test_data/CM3.out";

    std::ifstream instrmCM3(fNameCM3.c_str());
    std::string line;
    std::vector<std::string> tokens; 


    RDKit::ROMOL_SPTR mol1( RDKit::SmilesToMol( "NCCCCCO" ) );
    RDKit::ROMOL_SPTR mol( RDKit::MolOps::addHs( *mol1 ) );


    RDKit::DGeomHelpers::EmbedParameters params( RDKit::DGeomHelpers::ETKDG );
    params.randomSeed = 0xf00d;
    RDKit::DGeomHelpers::EmbedMolecule( *mol , params );
    std::vector<std::vector<double>> Mres;
    int nbmats= 10;
    int confId = -1;
    int padding = 23; // padding is the size of the result matrix so >= size max(atoms)
    int seed = 0xf00d;
   
    double rcut = 2.0; 

    // so local is ture 
    RDKit::Descriptors::CoulombMat(*mol, Mres, confId, nbmats, seed, padding, rcut, true, true, false, false, false, 1);

    //std::cout << "===================== vis BoB ========================\n";

    for ( const auto &v : Mres ) {

      std::getline(instrmCM3, line);
      tokens = tokenize(line);

      int ti = 0;
      for ( double x : v )  {
          //std::cout << x << ",";
          TEST_ASSERT(std::fabs(std::stof(tokens[ti]) - x)< 0.001);
          ti++;
      } 
      //std::cout << "\n";
    }
    //std::cout << "\n";
    Mres.clear();
    Mres.resize(nbmats);

    std::cout << "CM test 3 mat Done\n";
}



void testCoulombMat4(){
    std::cout << "===================== Testing CoulombMat Local 3 =======================\n";
 
    std::string pathName = "/Users/tgg/Github/rdkit"; //getenv("RDBASE");

    //////////////
    // TEST 3 LOCAL decaying 
    std::string fNameCM4 =
         pathName + "/Code/GraphMol/Descriptors/test_data/CM4.out";

    std::ifstream instrmCM4(fNameCM4.c_str());
    std::string line;
    std::vector<std::string> tokens; 


    RDKit::ROMOL_SPTR mol1( RDKit::SmilesToMol( "NCCCCCO" ) );
    RDKit::ROMOL_SPTR mol( RDKit::MolOps::addHs( *mol1 ) );


    RDKit::DGeomHelpers::EmbedParameters params( RDKit::DGeomHelpers::ETKDG );
    params.randomSeed = 0xf00d;
    RDKit::DGeomHelpers::EmbedMolecule( *mol , params );
    std::vector<std::vector<double>> Mres;
    int nbmats= 10;
    int confId = -1;
    int padding = 23; // padding is the size of the result matrix so >= size max(atoms)
    int seed = 0xf00d;
    double rcut = 2.5;
    int alpha = 3; 

    // so local is ture 
    RDKit::Descriptors::CoulombMat(*mol, Mres, confId, nbmats, seed, padding, rcut, true, true, false, false, false, alpha);

    //std::cout << "===================== vis BoB ========================\n";

    for ( const auto &v : Mres ) {

      std::getline(instrmCM4, line);
      tokens = tokenize(line);

      int ti = 0;
      for ( double x : v )  {
          //std::cout << x << ",";
          TEST_ASSERT(std::fabs(std::stof(tokens[ti]) - x)< 0.001);
          ti++;
      } 
      //std::cout << "\n";
    }
    //std::cout << "\n";
    Mres.clear();
    Mres.resize(nbmats);

    std::cout << "CM test 4 mat Done\n";
}


void testCoulombMat5(){
    std::cout << "===================== Testing CoulombMat Local 4 =======================\n";
 
    std::string pathName = "/Users/tgg/Github/rdkit"; //getenv("RDBASE");

    //////////////
    // TEST 4 LOCAL decaying & reduced
    std::string fNameCM5 =
         pathName + "/Code/GraphMol/Descriptors/test_data/CM5.out";

    std::ifstream instrmCM5(fNameCM5.c_str());
    std::string line;
    std::vector<std::string> tokens; 


    RDKit::ROMOL_SPTR mol1( RDKit::SmilesToMol( "NCCCCCO" ) );
    RDKit::ROMOL_SPTR mol( RDKit::MolOps::addHs( *mol1 ) );


    RDKit::DGeomHelpers::EmbedParameters params( RDKit::DGeomHelpers::ETKDG );
    params.randomSeed = 0xf00d;
    RDKit::DGeomHelpers::EmbedMolecule( *mol , params );
    std::vector<std::vector<double>> Mres;
    int nbmats= 10;
    int confId = -1;
    int padding = 23; // padding is the size of the result matrix so >= size max(atoms)
    int seed = 0xf00d;
    double rcut = 2.5;
    int alpha = 3; 

    // so local is ture 
    RDKit::Descriptors::CoulombMat(*mol, Mres, confId, nbmats, seed, padding, rcut, true, true, true, false, false, alpha);

    //std::cout << "===================== vis BoB ========================\n";

    for ( const auto &v : Mres ) {

      std::getline(instrmCM5, line);
      tokens = tokenize(line);

      int ti = 0;
      for ( double x : v )  {
          //std::cout << x << ",";
          TEST_ASSERT(std::fabs(std::stof(tokens[ti]) - x)< 0.001);
          ti++;
      } 
      //std::cout << "\n";
    }
    //std::cout << "\n";
    Mres.clear();
    Mres.resize(nbmats);

    std::cout << "CM test 5 mat Done\n";

} 

int main() {
  RDLog::InitLogs();
  testBagofBonds();
  testCoulombMat1();
  testCoulombMat2();
  testCoulombMat3();
  testCoulombMat4();
  testCoulombMat5();

}