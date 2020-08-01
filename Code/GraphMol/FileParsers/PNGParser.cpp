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
#include <boost/crc.hpp>

namespace RDKit {

namespace PNGData {
const std::string smilesTag = "rdkitSMILES";
const std::string molTag = "rdkitMOL";
const std::string jsonTag = "rdkitJSON";
}  // namespace PNGData

namespace {
std::vector<unsigned char> pngHeader = {137, 80, 78, 71, 13, 10, 26, 10};
bool checkPNGHeader(std::istream &inStream) {
  for (auto byte : pngHeader) {
    unsigned char ibyte;
    inStream.read((char *)&ibyte, 1);
    if (ibyte != byte) {
      return false;
    }
  }
  return true;
}
}  // namespace

std::map<std::string, std::string> PNGStreamToMetadata(std::istream &inStream) {
  // confirm that it's a PNG file:
  if (!checkPNGHeader(inStream)) {
    throw FileParseException("PNG header not recognized");
  }
  std::map<std::string, std::string> res;
  // the file is organized in chunks. Read through them until we find the tEXt
  // block FIX: at some point we'll want to also include zEXt here, but that
  // requires zlib
  while (inStream) {
    std::uint32_t blockLen;
    inStream.read((char *)&blockLen, sizeof(blockLen));
    // PNG is big endian, make sure we handle the order correctly
    blockLen = EndianSwapBytes<BIG_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(blockLen);
    char bytes[5];
    bytes[4] = 0;
    inStream.read(bytes, 4);
    auto beginBlock = inStream.tellg();
    if (bytes[0] == 'I' && bytes[1] == 'E' && bytes[2] == 'N' &&
        bytes[3] == 'D') {
      break;
    }
    if (blockLen > 0 && bytes[0] == 't' && bytes[1] == 'E' && bytes[2] == 'X' &&
        bytes[3] == 't') {
      // in a tEXt block, read the key:
      std::string key;
      std::getline(inStream, key, '\0');
      auto dataLen = blockLen - key.size() - 1;
      std::string value(dataLen, (char)0);
      inStream.read(&value.front(), dataLen);
      res[key] = value;
    }
    inStream.seekg(beginBlock);
    inStream.ignore(blockLen + 4);  // the extra 4 bytes are the CRC
  }

  return res;
};

std::string addMetadataToPNGStream(
    std::istream &inStream,
    const std::map<std::string, std::string> &metadata) {
  // confirm that it's a PNG file:
  if (!checkPNGHeader(inStream)) {
    throw FileParseException("PNG header not recognized");
  }
  std::stringstream res;
  // write the header
  for (auto byte : pngHeader) {
    res << byte;
  }

  // copy over everything up to IEND
  bool foundEnd = false;
  std::uint32_t finalCRC;
  while (inStream) {
    std::uint32_t blockLen;
    inStream.read((char *)&blockLen, sizeof(blockLen));
    char bytes[4];
    inStream.read(bytes, 4);
    if (bytes[0] == 'I' && bytes[1] == 'E' && bytes[2] == 'N' &&
        bytes[3] == 'D') {
      foundEnd = true;
      inStream.read((char *)&finalCRC, sizeof(finalCRC));
      break;
    }
    res.write((char *)&blockLen, sizeof(blockLen));
    res.write(bytes, 4);
    // PNG is big endian, make sure we handle the order correctly
    blockLen = EndianSwapBytes<BIG_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(blockLen);
    std::string block(blockLen + 4, 0);
    inStream.read((char *)&block.front(),
                  blockLen + 4);  // the extra 4 bytes are the CRC
    res.write(block.c_str(), blockLen + 4);
  }
  if (!foundEnd) {
    throw FileParseException("did not find IEND block in PNG");
  }

  // write out the metadata:
  for (const auto pr : metadata) {
    std::stringstream blk;
    blk.write("tEXt", 4);
    // write the name along with a zero
    blk.write(pr.first.c_str(), pr.first.size() + 1);
    blk.write(pr.second.c_str(), pr.second.size());
    auto blob = blk.str();
    std::uint32_t blksize =
        blob.size() - 4;  // we don't include the tag in the size;
    boost::crc_32_type crc;
    crc.process_bytes((void const *)blob.c_str(), blob.size());
    std::uint32_t crcVal = crc.checksum();
    // PNG is big endian, make sure we handle the order correctly
    blksize = EndianSwapBytes<HOST_ENDIAN_ORDER, BIG_ENDIAN_ORDER>(blksize);

    res.write((char *)&blksize, sizeof(blksize));
    res.write(blob.c_str(), blob.size());
    // PNG is big endian, make sure we handle the order correctly
    crcVal = EndianSwapBytes<HOST_ENDIAN_ORDER, BIG_ENDIAN_ORDER>(crcVal);
    res.write((char *)&crcVal, sizeof(crcVal));
  }

  // write out the IEND block
  std::uint32_t blksize = 0;
  res.write((char *)&blksize, sizeof(blksize));

  const char *endTag = "IEND";
  res.write(endTag, 4);
  res.write((char *)&finalCRC, sizeof(finalCRC));
  return res.str();
}

ROMol *PNGStreamToMol(std::istream &inStream,
                      const SmilesParserParams &params) {
  ROMol *res = nullptr;
  auto metadata = PNGStreamToMetadata(inStream);
  bool formatFound = false;
  for (const auto pr : metadata) {
    if (pr.first == PNGData::smilesTag) {
      res = SmilesToMol(pr.second, params);
      formatFound = true;
    }
  }
  if (!formatFound) {
    throw FileParseException("No suitable metadata found.");
  }
  return res;
}

}  // namespace RDKit
