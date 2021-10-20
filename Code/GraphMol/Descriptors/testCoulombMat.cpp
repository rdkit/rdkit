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

std::vector<std::string> tokenize(const std::string &s) {
  boost::char_separator<char> sep(", \n\r\t");
  boost::tokenizer<boost::char_separator<char>> tok(s, sep);
  std::vector<std::string> tokens;
  std::copy(tok.begin(), tok.end(),
            std::back_inserter<std::vector<std::string>>(tokens));
  return tokens;
}

void testCoulombMat1() {
  std::cerr
      << "===================== Testing CoulombMat  =======================\n";

  std::string pathName = getenv("RDBASE");

  // CM test 1
  std::string fNameCM =
      pathName + "/Code/GraphMol/Descriptors/test_data/CM1.out";

  std::string mol_file =
      pathName + "/Code/GraphMol/Descriptors/test_data/bobmol.sdf";

  std::ifstream instrmCM(fNameCM.c_str());
  std::string line;
  std::vector<std::string> tokens;

  RDKit::ROMOL_SPTR mol(RDKit::MolFileToMol(mol_file, true, false));

  std::vector<std::vector<double>> Mres;
  int confId = -1;
  RDKit::Descriptors::CoulombMat(*mol, Mres, confId);

  for (const auto &v : Mres) {
    std::getline(instrmCM, line);
    tokens = tokenize(line);

    unsigned int ti = 0;
    for (double x : v) {
      // std::cout << x << ",";
      TEST_ASSERT(std::fabs(std::stof(tokens[ti]) - x) < 0.001);
      ti++;
    }
    // std::cout << "\n";
  }
  // std::cout << "\n";

  std::cerr << "CM test 1 mat Done\n";
}

int main() {
  RDLog::InitLogs();
  testCoulombMat1();
}