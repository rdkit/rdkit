//  Copyright (C) 2019 Eisuke Kawashima
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "FileParsers.h"
#include <fstream>
#include <sstream>
#include <boost/format.hpp>
#include <RDGeneral/BadFileException.h>

namespace RDKit {

std::string MolToXYZBlock(const ROMol &mol, int confId,
                          unsigned int precision) {
  if (!mol.getNumConformers()) {
    BOOST_LOG(rdErrorLog)
        << "Cannot write molecules with no conformers to XYZ block\n";
    return "";
  }

  const auto &conf = mol.getConformer(confId);
  const unsigned int nAtoms = mol.getNumAtoms();

  std::stringstream ss;
  ss << nAtoms << '\n';

  unsigned fieldWidth = 5 + precision;
  std::stringstream formatString;
  formatString << "%-3s %" << fieldWidth << "." << precision << "f %"
               << fieldWidth << "." << precision << "f %" << fieldWidth << "."
               << precision << "f\n";

  std::string name;
  if (mol.getPropIfPresent(common_properties::_Name, name)) {
    ss << name.substr(0, name.find_first_of('\n'));
  }
  ss << '\n';

  for (unsigned int i = 0; i < nAtoms; i++) {
    const auto &symbol = mol.getAtomWithIdx(i)->getSymbol();
    const auto &pos = conf.getAtomPos(i);
    ss << boost::format{formatString.str()} % symbol % pos.x % pos.y % pos.z;
  }
  return ss.str();
}

void MolToXYZFile(const ROMol &mol, const std::string &fName, int confId,
                  unsigned int precision) {
  std::ofstream outStream(fName);
  if (!outStream) {
    std::ostringstream errout;
    errout << "Bad output file " << fName;
    throw BadFileException(errout.str());
  }
  outStream << MolToXYZBlock(mol, confId, precision);
}

}  // namespace RDKit
