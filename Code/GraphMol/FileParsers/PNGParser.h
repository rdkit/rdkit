//
//  Copyright (C) 2020 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_PNGPARSER_H
#define RD_PNGPARSER_H

#include <RDGeneral/types.h>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/RDKitBase.h>

#include <boost/format.hpp>

#include <string>
#include <fstream>
#include <map>

namespace RDKit {

RDKIT_FILEPARSERS_EXPORT std::map<std::string, std::string> PNGStreamToMetadata(
    std::istream &inStream);
RDKIT_FILEPARSERS_EXPORT std::map<std::string, std::string> PNGFileToMetadata(
    const std::string fname) {
  std::ifstream inStream(fname.c_str(), std::ios_base::binary);
  if (!inStream || (inStream.bad())) {
    throw BadFileException((boost::format("Bad input file %s") % fname).str());
  }
  return PNGStreamToMetadata(inStream);
};

}  // namespace RDKit

#endif
