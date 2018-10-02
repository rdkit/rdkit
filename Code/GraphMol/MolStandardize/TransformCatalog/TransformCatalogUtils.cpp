//
//  Copyright (C) 2018 Susan H. Leung
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
ChemicalReaction *getSmirks(const std::string &tmpStr) {
  ChemicalReaction *transformation = nullptr;
  if (tmpStr.length() == 0) {
    // empty line
    return transformation;
  }
  if (tmpStr.substr(0, 2) == "//") {
    // comment line
    return transformation;
  }

  boost::char_separator<char> tabSep("\t");
  tokenizer tokens(tmpStr, tabSep);
  tokenizer::iterator token = tokens.begin();

  // name of the functional groups
  std::string name = *token;
  boost::erase_all(name, " ");
  ++token;

  // grab the smirks:
  std::string smirks = *token;
  boost::erase_all(smirks, " ");
  ++token;

  transformation = RxnSmartsToChemicalReaction(smirks);
  CHECK_INVARIANT(transformation, smirks);
  transformation->setProp(common_properties::_Name, name);
  //  transformation->setProp(common_properties::_SMIRKS, smirks); // TODO
  //  RDGeneral/types.h does not have a common property to use?...
  return transformation;
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
  transformations.clear();
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
    // parse the reaction on this line (if there is one)
    std::shared_ptr<ChemicalReaction> transformation(getSmirks(tmpstr));
    if (transformation) {
      transformations.push_back(transformation);
      nRead++;
    }
  }
  return transformations;
}

}  // namespace MolStandardize
}  // namespace RDKit
