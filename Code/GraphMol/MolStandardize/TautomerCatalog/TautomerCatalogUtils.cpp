//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "TautomerCatalogUtils.h"
#include <RDGeneral/BadFileException.h>
#include <boost/tokenizer.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
#include <fstream>
#include <string>

namespace RDKit {
namespace {

MolStandardize::TautomerTransform* getTautomer(const std::string& tmpStr) {
  MolStandardize::TautomerTransform* transform = nullptr;
  if (tmpStr.length() == 0) {
    // empty line
    return transform;
  }
  if (tmpStr.substr(0, 2) == "//") {
    // comment line
    return transform;
  }
  boost::char_separator<char> tabSep("\t");
  tokenizer tokens(tmpStr, tabSep);
  std::vector<std::string> result(tokens.begin(), tokens.end());

  // tautomer information to collect from each line
  std::string name;
  std::string smarts;
  std::string bond_str;
  std::string charge_str;

  // line must have at least two tab separated values
  if (result.size() < 2) {
    std::cout << "Invalid line." << std::endl;
    return transform;
  }
  // line only has name and smarts
  if (result.size() == 2) {
    name = result[0];
    smarts = result[1];
  }
  // line has name, smarts, bonds
  if (result.size() == 3) {
    name = result[0];
    smarts = result[1];
    bond_str = result[2];
  }
  // line has name, smarts, bonds, charges
  if (result.size() == 4) {
    name = result[0];
    smarts = result[1];
    bond_str = result[2];
    charge_str = result[3];
  }

  boost::erase_all(smarts, " ");
  boost::erase_all(name, " ");
  boost::erase_all(bond_str, " ");
  boost::erase_all(charge_str, " ");

  std::vector<Bond::BondType> bond_types =
      MolStandardize::stringToBondType(bond_str);
  std::vector<int> charges = MolStandardize::stringToCharge(charge_str);

  ROMol* tautomer = SmartsToMol(smarts);
  CHECK_INVARIANT(tautomer, smarts);
  tautomer->setProp(common_properties::_Name, name);
  transform =
      new MolStandardize::TautomerTransform(tautomer, bond_types, charges);

  return transform;
}
}  // namespace

namespace MolStandardize {

std::vector<Bond::BondType> stringToBondType(std::string bond_str) {
  std::vector<Bond::BondType> bonds;
  for (const auto& c : bond_str) {
    switch (c) {
      case '-':
        bonds.push_back(Bond::SINGLE);
        break;
      case '=':
        bonds.push_back(Bond::DOUBLE);
        break;
      case '#':
        bonds.push_back(Bond::TRIPLE);
        break;
      case ':':
        bonds.push_back(Bond::AROMATIC);
        break;
    }
  }
  return bonds;
}

std::vector<int> stringToCharge(std::string charge_str) {
  std::vector<int> charges;
  for (const auto& c : charge_str) {
    switch (c) {
      case '+':
        charges.push_back(1);
        break;
      case '0':
        charges.push_back(0);
        break;
      case '-':
        charges.push_back(-1);
        break;
      default:
        throw ValueErrorException("Charge symbol not recognised.");
    }
  }
  return charges;
}

std::vector<TautomerTransform> readTautomers(std::string fileName) {
  std::ifstream inStream(fileName.c_str());
  if ((!inStream) || (inStream.bad())) {
    std::ostringstream errout;
    errout << "Bad input file " << fileName;
    throw BadFileException(errout.str());
  }
  std::vector<TautomerTransform> tautomers = readTautomers(inStream);
  return tautomers;
}

std::vector<TautomerTransform> readTautomers(std::istream& inStream,
                                             int nToRead) {
  if (inStream.bad()) {
    throw BadFileException("Bad stream contents.");
  }
  std::vector<TautomerTransform> tautomers;
  if (nToRead > 0) {
    tautomers.reserve(nToRead);
  }
  const int MAX_LINE_LEN = 512;
  char inLine[MAX_LINE_LEN];
  std::string tmpstr;
  int nRead = 0;
  while (!inStream.eof() && (nToRead < 0 || nRead < nToRead)) {
    inStream.getline(inLine, MAX_LINE_LEN, '\n');
    tmpstr = inLine;
    // parse the molpair on this line (if there is one)
    TautomerTransform* transform = getTautomer(tmpstr);
    if (transform != nullptr) {
      //			std::cout << MolToSmiles(*(transform->Mol) ) <<
      // std::endl;
      tautomers.emplace_back(*transform);
      delete transform;
      nRead++;
    }
  }

  return tautomers;
}

}  // namespace MolStandardize
}  // namespace RDKit
