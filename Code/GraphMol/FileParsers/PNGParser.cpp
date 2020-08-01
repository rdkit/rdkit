//
//  Copyright (C) 2020 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// details of how to handle the PNG file taken from OpenBabel's PNG handling
// code:
// https://github.com/openbabel/openbabel/blob/master/src/formats/pngformat.cpp

#include "PNGParser.h"
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/StreamOps.h>
#include <vector>

namespace RDKit {

RDKIT_FILEPARSERS_EXPORT std::map<std::string, std::string> PNGStreamToMetadata(
    std::istream &inStream) {
  // confirm that it's a PNG file:
  std::vector<unsigned char> header = {137, 80, 78, 71, 13, 10, 26, 10};
  for (auto byte : header) {
    unsigned char ibyte;
    inStream.read((char *)&ibyte, 1);
    if (ibyte != byte) {
      throw FileParseException("PNG header not recognized");
    }
  }

  // the file is organized in chunks. Read through them until we find the tEXt
  // block FIX: at some point we'll want to also include zEXt here, but that
  // requires zlib
  while (inStream) {
    std::uint32_t blockLen;
    inStream.read((char *)&blockLen, sizeof(blockLen));
    blockLen = EndianSwapBytes<BIG_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(blockLen);
    std::cerr << "  block len: " << blockLen << std::endl;
    char bytes[4];
    inStream.read(bytes, 4);
    auto beginBlock = inStream.tellg();
    std::cerr << "bytes: " << bytes << std::endl;
    if (bytes[0] == 'I' && bytes[1] == 'E' && bytes[2] == 'N' &&
        bytes[3] == 'D') {
      std::cerr << " done! " << std::endl;
      break;
    }
    if (bytes[0] == 't' && bytes[1] == 'E' && bytes[2] == 'X' &&
        bytes[3] == 't') {
      std::cerr << "tEXt!" << std::endl;
    }
    inStream.seekg(beginBlock);
    inStream.ignore(blockLen + 4);  // the extra 4 bytes are the CRC
  }

  std::map<std::string, std::string> res;
  return res;
};

}  // namespace RDKit
