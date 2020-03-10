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
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>
#include <Eigen/Dense>

using namespace Eigen;

std::vector<std::string> tokenize(const std::string &s) {
        boost::char_separator<char> sep(", \n\r\t");
        boost::tokenizer<boost::char_separator<char>> tok(s, sep);
        std::vector<std::string> tokens;
        std::copy(tok.begin(), tok.end(), std::back_inserter<std::vector<std::string> >(tokens));
  return tokens;
}



void testCoulombMat1(){
    std::cout << "===================== Testing CoulombMat GLobal =======================\n";
 
    std::string pathName = getenv("RDBASE");

    // CM test 1
    std::string fNameCM =
         pathName + "/Code/GraphMol/Descriptors/test_data/CM1.out";

    std::string mol_file =
         pathName + "/Code/GraphMol/Descriptors/test_data/bobmol.sdf";

    std::ifstream instrmCM(fNameCM.c_str());
    std::string line;
    std::vector<std::string> tokens; 

    RDKit::ROMOL_SPTR mol( RDKit::MolFileToMol( mol_file , true, false) );

    std::vector<std::vector<double>> Mres;
    int nbmats= 10;
    int confId = -1;
    int padding = 23; // padding is the size of the result matrix so >= size max(atoms)
    int seed = 0xf00d;
    double rcut = 0; // not used if local is false

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
 
    std::string pathName = getenv("RDBASE");

    ///////////////////////
    // TEST 1 LOCAL
    std::string fNameCM2 =
         pathName + "/Code/GraphMol/Descriptors/test_data/CM2.out";

    std::string mol_file =
         pathName + "/Code/GraphMol/Descriptors/test_data/bobmol.sdf";

    std::ifstream instrmCM2(fNameCM2.c_str());

    std::string line;
    std::vector<std::string> tokens; 

    RDKit::ROMOL_SPTR mol( RDKit::MolFileToMol( mol_file , true, false) );

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
 
    std::string pathName = getenv("RDBASE");

    //////////////
    // TEST 2 LOCAL decaying
    std::string fNameCM3 =
         pathName + "/Code/GraphMol/Descriptors/test_data/CM3.out";

    std::string mol_file =
         pathName + "/Code/GraphMol/Descriptors/test_data/bobmol.sdf";
 
    std::ifstream instrmCM3(fNameCM3.c_str());
    std::string line;
    std::vector<std::string> tokens; 

    RDKit::ROMOL_SPTR mol( RDKit::MolFileToMol( mol_file , true, false) );

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
 
    std::string pathName = getenv("RDBASE");

    //////////////
    // TEST 3 LOCAL decaying 
    std::string fNameCM4 =
         pathName + "/Code/GraphMol/Descriptors/test_data/CM4.out";

    std::string mol_file =
         pathName + "/Code/GraphMol/Descriptors/test_data/bobmol.sdf";

    std::ifstream instrmCM4(fNameCM4.c_str());
    std::string line;
    std::vector<std::string> tokens; 

    RDKit::ROMOL_SPTR mol( RDKit::MolFileToMol( mol_file , true, false) );

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
 
    std::string pathName = getenv("RDBASE");

    //////////////
    // TEST 4 LOCAL decaying & reduced
    std::string fNameCM5 =
         pathName + "/Code/GraphMol/Descriptors/test_data/CM5.out";

    std::string mol_file =
         pathName + "/Code/GraphMol/Descriptors/test_data/bobmol.sdf";

    std::ifstream instrmCM5(fNameCM5.c_str());
    std::string line;
    std::vector<std::string> tokens; 

    RDKit::ROMOL_SPTR mol( RDKit::MolFileToMol( mol_file , true, false) );

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
  testCoulombMat1();
  testCoulombMat2();
  testCoulombMat3();
  testCoulombMat4();
  testCoulombMat5();

}