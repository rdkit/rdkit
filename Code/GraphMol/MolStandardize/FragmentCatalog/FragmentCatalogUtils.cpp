//
//  Copyright (C) 2018 Susan H. Leung
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
ROMol *getSmarts(const std::string &tmpStr) {
  ROMol *mol = nullptr;
  if (tmpStr.length() == 0) {
    // empty line
    return mol;
  }
  if (tmpStr.substr(0, 2) == "//") {
    // comment line
    return mol;
  }

  boost::char_separator<char> tabSep("\t");
  tokenizer tokens(tmpStr, tabSep);
  tokenizer::iterator token = tokens.begin();

  // name of the functional groups
  std::string name = *token;
  //  boost::erase_all(name, " ");
  ++token;

  // grab the smarts:
  std::string smarts = *token;
  boost::erase_all(smarts, " ");
  ++token;

  mol = SmartsToMol(smarts);
  CHECK_INVARIANT(mol, smarts);
  mol->setProp(common_properties::_Name, name);
  mol->setProp(common_properties::_fragSMARTS, smarts);
  return mol;
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
  std::vector<std::shared_ptr<ROMol>> funcGroups;
  funcGroups.clear();
  if (inStream.bad()) {
    throw BadFileException("Bad stream contents.");
  }
  const int MAX_LINE_LEN = 512;
  char inLine[MAX_LINE_LEN];
  std::string tmpstr;
  int nRead = 0;
  while (!inStream.eof() && (nToRead < 0 || nRead < nToRead)) {
    inStream.getline(inLine, MAX_LINE_LEN, '\n');
    tmpstr = inLine;
    // parse the molecule on this line (if there is one)
    std::shared_ptr<ROMol> mol(getSmarts(tmpstr));
    if (mol) {
      funcGroups.push_back(mol);
      nRead++;
    }
  }
  return funcGroups;
}

}  // namespace MolStandardize
}  // namespace RDKit
