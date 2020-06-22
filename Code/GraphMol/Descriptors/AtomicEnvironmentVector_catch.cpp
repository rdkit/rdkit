//
//  Copyright (C) 2020 Manan Goel
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "RDGeneral/test.h"
#include "catch.hpp"

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
#include <GraphMol/Descriptors/AtomicEnvironmentVector.h>
#include <Eigen/Dense>

using namespace Eigen;

std::vector<std::string> tokenize(const std::string &s) {
  boost::char_separator<char> sep(", \n\r\t");
  boost::tokenizer<boost::char_separator<char>> tok(s, sep);
  std::vector<std::string> tokens;
  std::copy(tok.begin(), tok.end(),
            std::back_inserter<std::vector<std::string>>(tokens));
  return tokens;
}

std::vector<std::vector<double>> EigenMatToSTLVector(ArrayXXd aev) {
  std::vector<std::vector<double>> aevOutput;
  for (auto i = 0; i < aev.rows(); i++) {
    std::vector<double> row;
    for (auto j = 0; j < aev.cols(); j++) {
      row.push_back(aev(i, j));
    }
    aevOutput.push_back(row);
  }
  return aevOutput;
}

TEST_CASE("Symmetry Function Accuracy", "[Symmetry Function]") {
  SECTION("CH4") {
    std::string pathName = getenv("RDBASE");

    std::string molFile =
        pathName + "/Code/GraphMol/Descriptors/test_data/CH4.mol";

    std::string fNameSF =
        pathName + "/Code/GraphMol/Descriptors/test_data/CH4.out";

    std::ifstream instrmSF;
    instrmSF.open(fNameSF);

    std::string line;
    std::vector<std::string> tokens;
    std::vector<std::vector<double>> expectedOutput;

    while (!instrmSF.eof()) {
      std::getline(instrmSF, line);
      tokens = tokenize(line);
      std::vector<double> row;
      for (auto v : tokens) {
        std::istringstream os(v);
        double d;
        os >> d;
        row.push_back(d);
      }
      expectedOutput.push_back(row);
    }
    // Each row of the received vector is sorted because in angular terms
    // order of triplets may change because of unstable sorting
    // To check if the values of the terms are consistent each row of the output
    // is sorted

    RDKit::ROMOL_SPTR mol(RDKit::MolFileToMol(molFile, true, false));
    int confId = -1;

    auto aev = RDKit::Descriptors::ANI::AtomicEnvironmentVector(*mol, confId);

    CHECK(aev.rows() == mol->getNumAtoms());
    CHECK(aev.cols() == 384);

    auto aevOutput = EigenMatToSTLVector(aev);

    for (unsigned int i = 0; i < expectedOutput.size(); i++) {
      for (unsigned int j = 0; j < expectedOutput[i].size(); j++) {
        auto diff = std::fabs(expectedOutput[i][j] - aevOutput[i][j]);
        CHECK(diff < 0.2);
      }
    }
  }

  SECTION("NH3") {
    std::string pathName = getenv("RDBASE");

    std::string molFile =
        pathName + "/Code/GraphMol/Descriptors/test_data/NH3.mol";

    std::string fNameSF =
        pathName + "/Code/GraphMol/Descriptors/test_data/NH3.out";

    std::ifstream instrmSF;
    instrmSF.open(fNameSF);

    std::string line;
    std::vector<std::string> tokens;
    std::vector<std::vector<double>> expectedOutput;

    while (!instrmSF.eof()) {
      std::getline(instrmSF, line);
      tokens = tokenize(line);
      std::vector<double> row;
      for (auto v : tokens) {
        std::istringstream os(v);
        double d;
        os >> d;
        row.push_back(d);
      }
      expectedOutput.push_back(row);
    }

    RDKit::ROMOL_SPTR mol(RDKit::MolFileToMol(molFile, true, false));
    int confId = -1;

    auto aev = RDKit::Descriptors::ANI::AtomicEnvironmentVector(*mol, confId);

    CHECK(aev.rows() == mol->getNumAtoms());
    CHECK(aev.cols() == 384);

    auto aevOutput = EigenMatToSTLVector(aev);

    for (unsigned int i = 0; i < expectedOutput.size(); i++) {
      for (unsigned int j = 0; j < expectedOutput[i].size(); j++) {
        auto diff = std::fabs(expectedOutput[i][j] - aevOutput[i][j]);
        CHECK(diff < 0.2);
      }
    }
  }

  SECTION("SO2") {
    std::string pathName = getenv("RDBASE");

    std::string molFile =
        pathName + "/Code/GraphMol/Descriptors/test_data/SO2.mol";
    std::string fNameSF =
        pathName + "/Code/GraphMol/Descriptors/test_data/SO2.out";

    std::ifstream instrmSF;
    instrmSF.open(fNameSF);

    std::string line;
    std::vector<std::string> tokens;
    std::vector<std::vector<double>> expectedOutput;

    while (!instrmSF.eof()) {
      std::getline(instrmSF, line);
      tokens = tokenize(line);
      std::vector<double> row;
      for (auto v : tokens) {
        std::istringstream os(v);
        double d;
        os >> d;
        row.push_back(d);
      }
      expectedOutput.push_back(row);
    }

    RDKit::ROMOL_SPTR mol(RDKit::MolFileToMol(molFile, true, false));
    int confId = -1;

    auto aev = RDKit::Descriptors::ANI::AtomicEnvironmentVector(*mol, confId);

    CHECK(aev.rows() == mol->getNumAtoms());
    CHECK(aev.cols() == 384);

    auto aevOutput = EigenMatToSTLVector(aev);

    for (unsigned int i = 0; i < expectedOutput.size(); i++) {
      for (unsigned int j = 0; j < expectedOutput[i].size(); j++) {
        auto diff = std::fabs(expectedOutput[i][j] - aevOutput[i][j]);
        CHECK(diff < 0.2);
      }
    }
  }

  SECTION("Ethanol") {
    std::string pathName = getenv("RDBASE");

    std::string molFile =
        pathName + "/Code/GraphMol/Descriptors/test_data/ethanol.sdf";

    std::string fNameSF =
        pathName + "/Code/GraphMol/Descriptors/test_data/ethanol.out";

    std::ifstream instrmSF;
    instrmSF.open(fNameSF);

    std::string line;
    std::vector<std::string> tokens;
    std::vector<std::vector<double>> expectedOutput;

    while (!instrmSF.eof()) {
      std::getline(instrmSF, line);
      tokens = tokenize(line);
      std::vector<double> row;
      for (auto v : tokens) {
        std::istringstream os(v);
        double d;
        os >> d;
        row.push_back(d);
      }
      expectedOutput.push_back(row);
    }

    RDKit::ROMOL_SPTR mol(RDKit::MolFileToMol(molFile, false, false));
    int confId = -1;

    auto aev = RDKit::Descriptors::ANI::AtomicEnvironmentVector(*mol, confId);

    CHECK(aev.rows() == mol->getNumAtoms());
    CHECK(aev.cols() == 384);


    auto aevOutput = EigenMatToSTLVector(aev);


    for (unsigned int i = 0; i < expectedOutput.size(); i++) {
      for (unsigned int j = 0; j < expectedOutput[i].size(); j++) {
        auto diff = std::fabs(expectedOutput[i][j] - aevOutput[i][j]);
        CHECK(diff < 0.2);
      }
    }
  }
}

TEST_CASE("Species Vector Generation", "[Atomic Species Encoding]") {
  SECTION("CH4") {
    std::string pathName = getenv("RDBASE");

    std::string molFile =
        pathName + "/Code/GraphMol/Descriptors/test_data/CH4.mol";
    RDKit::ROMOL_SPTR mol(RDKit::MolFileToMol(molFile, true, false));
    auto speciesVec = RDKit::Descriptors::ANI::GenerateSpeciesVector(*mol);
    auto numAtoms = mol->getNumAtoms();
    CHECK(speciesVec.size() == numAtoms);
    VectorXi expected(5);
    expected << 0, 0, 0, 1, 0;
    CHECK(((expected - speciesVec).array() == 0).count() == numAtoms);
    int atomNums[] = {1, 1, 1, 6, 1};
    speciesVec = RDKit::Descriptors::ANI::GenerateSpeciesVector(atomNums, 5);
    CHECK(speciesVec.size() == 5);
    CHECK(((expected - speciesVec).array() == 0).count() == numAtoms);
  }
  SECTION("SO2") {
    std::string pathName = getenv("RDBASE");

    std::string molFile =
        pathName + "/Code/GraphMol/Descriptors/test_data/SO2.mol";
    RDKit::ROMOL_SPTR mol(RDKit::MolFileToMol(molFile, true, false));
    auto speciesVec = RDKit::Descriptors::ANI::GenerateSpeciesVector(*mol);
    auto numAtoms = mol->getNumAtoms();
    CHECK(speciesVec.size() == numAtoms);
    VectorXi expected(3);
    expected << -1, 3, 3;
    CHECK(((expected - speciesVec).array() == 0).count() == numAtoms);
    int atomNums[] = {16, 8, 8};
    speciesVec = RDKit::Descriptors::ANI::GenerateSpeciesVector(atomNums, 3);
    CHECK(speciesVec.size() == 3);
    CHECK(((expected - speciesVec).array() == 0).count() == numAtoms);
  }
  SECTION("NH3") {
    std::string pathName = getenv("RDBASE");

    std::string molFile =
        pathName + "/Code/GraphMol/Descriptors/test_data/NH3.mol";
    RDKit::ROMOL_SPTR mol(RDKit::MolFileToMol(molFile, true, false));
    auto speciesVec = RDKit::Descriptors::ANI::GenerateSpeciesVector(*mol);
    auto numAtoms = mol->getNumAtoms();
    CHECK(speciesVec.size() == numAtoms);
    VectorXi expected(4);
    expected << 0, 0, 2, 0;
    CHECK(((expected - speciesVec).array() == 0).count() == numAtoms);
    int atomNums[] = {1, 1, 7, 1};
    speciesVec = RDKit::Descriptors::ANI::GenerateSpeciesVector(atomNums, 4);
    CHECK(speciesVec.size() == 4);
    CHECK(((expected - speciesVec).array() == 0).count() == numAtoms);
  }
  SECTION("Ethanol") {
    std::string pathName = getenv("RDBASE");

    std::string molFile =
        pathName + "/Code/GraphMol/Descriptors/test_data/ethanol.sdf";
    RDKit::ROMOL_SPTR mol(RDKit::MolFileToMol(molFile, true, false));
    auto speciesVec = RDKit::Descriptors::ANI::GenerateSpeciesVector(*mol);
    auto numAtoms = mol->getNumAtoms();
    CHECK(speciesVec.size() == numAtoms);
    VectorXi expected(9);
    expected << 3, 1, 1, 0, 0, 0, 0, 0, 0;
    CHECK(((expected - speciesVec).array() == 0).count() == numAtoms);


    int atomNums[] = {8, 6, 6, 1, 1, 1, 1, 1, 1};
    speciesVec = RDKit::Descriptors::ANI::GenerateSpeciesVector(atomNums, 9);
    CHECK(speciesVec.size() == sizeof(atomNums)/sizeof(int));
    CHECK(((expected - speciesVec).array() == 0).count() == numAtoms);
  }
}