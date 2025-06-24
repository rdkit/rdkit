//
//  Copyright (C) 2018-2021 Susan H. Leung and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "AcidBaseCatalogUtils.h"
#include <RDGeneral/BadFileException.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
#include <fstream>
#include <string>

namespace RDKit {
namespace {

std::pair<ROMol *, ROMol *> *getPair(const std::string &name,
                                     const std::string &acid_smarts,
                                     const std::string &base_smarts) {
  ROMol *acid(SmartsToMol(acid_smarts));
  if (!acid) {
    throw ValueErrorException("Failed parsing acid SMARTS: " + acid_smarts);
  }
  ROMol *base(SmartsToMol(base_smarts));
  if (!base) {
    delete acid;
    throw ValueErrorException("Failed parsing base SMARTS: " + base_smarts);
  }

  acid->setProp(common_properties::_Name, name);
  base->setProp(common_properties::_Name, name);

  return new std::pair<ROMol *, ROMol *>(acid, base);
}
std::pair<ROMol *, ROMol *> *getPair(std::string tmpStr) {
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
  boost::erase_all(name, " ");
  ++token;

  // grab the acid:
  std::string acid_smarts = *token;
  boost::erase_all(acid_smarts, " ");
  ++token;

  // grab the base:
  std::string base_smarts = *token;
  boost::erase_all(base_smarts, " ");
  ++token;

  return getPair(name, acid_smarts, base_smarts);
}
}  // namespace

namespace MolStandardize {

std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> readPairs(std::string fileName) {
  std::ifstream inStream(fileName.c_str());
  if ((!inStream) || (inStream.bad())) {
    std::ostringstream errout;
    errout << "Bad input file " << fileName;
    throw BadFileException(errout.str());
  }
  std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> mol_pairs =
      readPairs(inStream);
  return mol_pairs;
}

std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> readPairs(std::istream &inStream,
                                                         int nToRead) {
  std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> pairs;
  pairs.clear();
  if (inStream.bad()) {
    throw BadFileException("Bad stream contents.");
  }
  const int MAX_LINE_LEN = 512;
  char inLine[MAX_LINE_LEN];
  std::string tmpstr;
  int nRead = 0;
  while (!inStream.eof() && !inStream.fail() &&
         (nToRead < 0 || nRead < nToRead)) {
    inStream.getline(inLine, MAX_LINE_LEN, '\n');
    tmpstr = inLine;
    // parse the molpair on this line (if there is one)
    std::shared_ptr<std::pair<ROMol *, ROMol *>> mol_pair(getPair(tmpstr));
    if (mol_pair != nullptr) {
      pairs.emplace_back(ROMOL_SPTR(mol_pair->first),
                         ROMOL_SPTR(mol_pair->second));
      nRead++;
    }
  }
  return pairs;
}

std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> readPairs(
    const std::vector<std::tuple<std::string, std::string, std::string>>
        &data) {
  std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> pairs;
  for (const auto &tpl : data) {
    std::shared_ptr<std::pair<ROMol *, ROMol *>> mol_pair(
        getPair(std::get<0>(tpl), std::get<1>(tpl), std::get<2>(tpl)));
    pairs.emplace_back(ROMOL_SPTR(mol_pair->first),
                       ROMOL_SPTR(mol_pair->second));
  }
  return pairs;
}
}  // namespace MolStandardize
}  // namespace RDKit
