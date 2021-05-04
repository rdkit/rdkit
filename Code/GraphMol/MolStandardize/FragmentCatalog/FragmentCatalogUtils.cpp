//
//  Copyright (C) 2018-2021 Susan H. Leung and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "FragmentCatalogUtils.h"
#include <RDGeneral/BadFileException.h>
#include <boost/tokenizer.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
#include <fstream>
#include <string>

namespace RDKit {
namespace {
ROMol *getMol(const std::string &name, const std::string &smarts) {
  auto mol = SmartsToMol(smarts);
  if (!mol) {
    throw ValueErrorException("Failed parsing fragment SMARTS: " + smarts);
  }
  mol->setProp(common_properties::_Name, name);
  mol->setProp(common_properties::_fragSMARTS, smarts);
  return mol;
}

ROMol *getMol(std::string &&tmpStr) {
  // Remove whitespace
  boost::trim(tmpStr);

  if (tmpStr.length() == 0 || tmpStr.substr(0, 2) == "//") {
    // empty or comment line
    return nullptr;
  }

  boost::char_separator<char> tabSep("\t");
  tokenizer tokens(tmpStr, tabSep);
  tokenizer::iterator token = tokens.begin();

  // name of the functional groups
  std::string name = *token;
  ++token;

  // There must be a SMARTS expression
  if (token == tokens.end()) {
    throw ValueErrorException("no SMARTS found in input: " + tmpStr);
  }

  // grab the smarts:
  std::string smarts = *token;
  boost::erase_all(smarts, " ");
  ++token;
  return getMol(name, smarts);
}
}  // namespace

namespace MolStandardize {

std::vector<std::shared_ptr<ROMol>> readFuncGroups(std::string fileName) {
  std::ifstream inStream(fileName.c_str());
  if ((!inStream) || (inStream.bad())) {
    std::ostringstream errout;
    errout << "Bad input file " << fileName;
    throw BadFileException(errout.str());
  }
  std::vector<std::shared_ptr<ROMol>> funcGroups;
  funcGroups = readFuncGroups(inStream);
  return funcGroups;
}

std::vector<std::shared_ptr<ROMol>> readFuncGroups(std::istream &inStream,
                                                   int nToRead) {
  if (inStream.bad()) {
    throw BadFileException("Bad stream contents.");
  }

  const int MAX_LINE_LEN = 512;
  char inLine[MAX_LINE_LEN];
  int nRead = 0;

  std::vector<std::shared_ptr<ROMol>> funcGroups;
  while (!inStream.eof() && !inStream.fail() &&
         (nToRead < 0 || nRead < nToRead)) {
    inStream.getline(inLine, MAX_LINE_LEN, '\n');
    std::string tmpstr(inLine);
    // parse the molecule on this line (if there is one)
    std::shared_ptr<ROMol> mol(getMol(std::move(tmpstr)));
    if (mol) {
      funcGroups.push_back(mol);
      nRead++;
    }
  }
  return funcGroups;
}

std::vector<std::shared_ptr<ROMol>> readFuncGroups(
    const std::vector<std::pair<std::string, std::string>> &data) {
  std::vector<std::shared_ptr<ROMol>> funcGroups;
  for (const auto &pr : data) {
    auto mol = getMol(pr.first, pr.second);
    if (mol) {
      funcGroups.push_back(std::shared_ptr<ROMol>(mol));
    }
  }
  return funcGroups;
}

}  // namespace MolStandardize
}  // namespace RDKit
