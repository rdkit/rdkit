//
//  Copyright (C) 2022 Sreya Gogineni and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <boost/lexical_cast.hpp>

#include "FileParsers.h"
#include "FileParserUtils.h"
#include <RDGeneral/StreamOps.h>

#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <exception>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <vector>

namespace RDKit {

void ParseExtraLine(const std::string &extraLine) {
  std::string whitespace{" \t"};
  if (extraLine.find_first_not_of(whitespace) != std::string::npos) {
    std::ostringstream errout;
    errout << "More lines than expected" << std::endl;
    throw FileParseException(errout.str());
  }
}

Atom *ParseXYZFileAtomLine(const std::string &atomLine, RDGeom::Point3D &pos,
                           unsigned int line) {
  std::string whitespace{" \t"};
  size_t delims[8];
  size_t prev = 0;
  for (unsigned int i = 0; i < 7; i++) {
    if (i % 2 == 0) {
      delims[i] = atomLine.find_first_not_of(whitespace, prev);
    } else {
      delims[i] = atomLine.find_first_of(whitespace, prev);
    }
    if (delims[i] == std::string::npos) {
      std::ostringstream errout;
      errout << "Missing coordinates on line " << line << std::endl;
      throw FileParseException(errout.str());
    }
    prev = delims[i];
  }
  delims[7] = atomLine.find_last_not_of(whitespace) + 1;

  // set conformer
  try {
    pos.x = FileParserUtils::toDouble(
        atomLine.substr(delims[2], delims[3] - delims[2]), false);
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '"
           << atomLine.substr(delims[2], delims[3] - delims[2])
           << "' to double on line " << line << std::endl;
    throw FileParseException(errout.str());
  }

  try {
    pos.y = FileParserUtils::toDouble(
        atomLine.substr(delims[4], delims[5] - delims[4]), false);
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '"
           << atomLine.substr(delims[4], delims[5] - delims[4])
           << "' to double on line " << line << std::endl;
    throw FileParseException(errout.str());
  }

  try {
    pos.z = FileParserUtils::toDouble(
        atomLine.substr(delims[6], delims[7] - delims[6]), false);
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '"
           << atomLine.substr(delims[6], delims[7] - delims[6])
           << "' to double on line " << line << std::endl;
    throw FileParseException(errout.str());
  }

  std::string symb{atomLine.substr(delims[0], delims[1] - delims[0])};
  if (symb.size() == 2 && symb[1] >= 'A' && symb[1] <= 'Z') {
    symb[1] = static_cast<char>(tolower(symb[1]));
  }

  Atom *atom;
  try {
    atom = new Atom(PeriodicTable::getTable()->getAtomicNumber(symb));
  } catch (const Invar::Invariant &e) {
    throw FileParseException(e.what());
  }

  return atom;
}

RWMol *XYZDataStreamToMol(std::istream &inStream) {
  unsigned int numAtoms = 0;

  std::string num{getLine(inStream)};
  try {
    numAtoms = FileParserUtils::toUnsigned(num);
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Unable to recognize the number of atoms: cannot convert '" << num
           << "' to unsigned int on line 0" << std::endl;
    throw FileParseException(errout.str());
  }

  std::string comment{getLine(inStream)};

  RWMol *mol = new RWMol();
  if (numAtoms) {
    Conformer *conf = new Conformer(numAtoms);
    if (!comment.empty()) {
      mol->setProp("_FileComments", comment);
    }
    for (unsigned int i = 0; i < numAtoms; i++) {
      if (inStream.eof()) {
        throw FileParseException("EOF hit while reading atoms");
      }
      RDGeom::Point3D pos;
      std::string atomLine{getLine(inStream)};
      Atom *atom = ParseXYZFileAtomLine(atomLine, pos, i + 2);
      unsigned int idx = mol->addAtom(atom, false, true);
      conf->setAtomPos(idx, pos);
    }
    mol->addConformer(conf);
  }

  while (!inStream.eof()) {
    std::string extraLine{getLine(inStream)};
    ParseExtraLine(extraLine);
  }

  return mol;
}

RWMol *XYZBlockToMol(const std::string &xyzBlock) {
  std::istringstream xyz(xyzBlock);

  RWMol *mol = nullptr;
  xyz.peek();
  if (!xyz.eof()) {
    mol = XYZDataStreamToMol(xyz);
  }

  return mol;
}

RWMol *XYZFileToMol(const std::string &fName) {
  std::ifstream xyzFile(fName);
  if (!xyzFile || (xyzFile.bad())) {
    std::ostringstream errout;
    errout << "Bad input file " << fName;
    throw BadFileException(errout.str());
  }

  RWMol *mol = nullptr;
  xyzFile.peek();
  if (!xyzFile.eof()) {
    mol = XYZDataStreamToMol(xyzFile);
  }

  return mol;
}

}  // namespace RDKit
