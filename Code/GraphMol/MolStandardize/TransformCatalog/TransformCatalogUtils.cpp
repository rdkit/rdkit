//
//  Copyright (C) 2018-2021 Susan H. Leung and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "TransformCatalogUtils.h"
#include <RDGeneral/BadFileException.h>
#include <boost/tokenizer.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
#include <fstream>
#include <string>

namespace RDKit {
namespace {
ChemicalReaction *getReaction(const std::string &name,
                              const std::string &smarts) {
  auto transformation = RxnSmartsToChemicalReaction(smarts);
  if (!transformation) {
    throw ValueErrorException("Failed parsing reaction SMARTS: " + smarts);
  }
  transformation->setProp(common_properties::_Name, name);
  return transformation;
}

ChemicalReaction *getReaction(const std::string &tmpStr) {
  if (tmpStr.length() == 0 || tmpStr.substr(0, 2) == "//") {
    // empty line or comment
    return nullptr;
  }

  boost::char_separator<char> tabSep("\t");
  tokenizer tokens(tmpStr, tabSep);
  tokenizer::iterator token = tokens.begin();

  // name of the functional groups
  auto name = *token;
  boost::erase_all(name, " ");
  ++token;

  // grab the smirks:
  auto smirks = *token;
  boost::erase_all(smirks, " ");
  ++token;
  return getReaction(name, smirks);
}
}  // namespace

namespace MolStandardize {

std::vector<std::shared_ptr<ChemicalReaction>> readTransformations(
    std::string fileName) {
  std::ifstream inStream(fileName.c_str());
  if ((!inStream) || (inStream.bad())) {
    std::ostringstream errout;
    errout << "Bad input file " << fileName;
    throw BadFileException(errout.str());
  }
  std::vector<std::shared_ptr<ChemicalReaction>> transformations;
  transformations = readTransformations(inStream);
  return transformations;
}

std::vector<std::shared_ptr<ChemicalReaction>> readTransformations(
    std::istream &inStream, int nToRead) {
  std::vector<std::shared_ptr<ChemicalReaction>> transformations;
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
    // parse the reaction on this line (if there is one)
    auto transformation = getReaction(tmpstr);
    if (transformation) {
      transformations.push_back(
          std::shared_ptr<ChemicalReaction>(transformation));
      nRead++;
    }
  }
  return transformations;
}

std::vector<std::shared_ptr<ChemicalReaction>> readTransformations(
    const std::vector<std::pair<std::string, std::string>> &data) {
  std::vector<std::shared_ptr<ChemicalReaction>> transformations;
  for (const auto &pr : data) {
    auto transformation = getReaction(pr.first, pr.second);
    if (transformation) {
      transformations.push_back(
          std::shared_ptr<ChemicalReaction>(transformation));
    }
  }
  return transformations;
}

}  // namespace MolStandardize
}  // namespace RDKit
