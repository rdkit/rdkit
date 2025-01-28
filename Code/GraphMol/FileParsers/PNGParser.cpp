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
#include <GraphMol/MolPickler.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/StreamOps.h>
#include <vector>
#include <boost/crc.hpp>
#include <boost/algorithm/string.hpp>

#include "FileParsers.h"
#ifdef RDK_USE_BOOST_IOSTREAMS
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#endif

namespace RDKit {

namespace PNGData {
const std::string smilesTag = "SMILES";
const std::string molTag = "MOL";
const std::string pklTag = "rdkitPKL";
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

#ifdef RDK_USE_BOOST_IOSTREAMS
std::string uncompressString(const std::string &ztext) {
  std::stringstream compressed(ztext);
  std::stringstream uncompressed;
  boost::iostreams::filtering_streambuf<boost::iostreams::input> bioOutstream;
  bioOutstream.push(boost::iostreams::zlib_decompressor());
  bioOutstream.push(compressed);
  boost::iostreams::copy(bioOutstream, uncompressed);
  return uncompressed.str();
}
std::string compressString(const std::string &text) {
  std::stringstream uncompressed(text);
  std::stringstream compressed;
  boost::iostreams::filtering_streambuf<boost::iostreams::input> bioOutstream;
  bioOutstream.push(boost::iostreams::zlib_compressor());
  bioOutstream.push(uncompressed);
  boost::iostreams::copy(bioOutstream, compressed);
  return compressed.str();
}

#endif
}  // namespace

std::vector<std::pair<std::string, std::string>> PNGStreamToMetadata(
    std::istream &inStream) {
  // confirm that it's a PNG file:
  if (!checkPNGHeader(inStream)) {
    throw FileParseException("PNG header not recognized");
  }
  std::vector<std::pair<std::string, std::string>> res;
  // the file is organized in chunks. Read through them until we find the tEXt
  // block FIX: at some point we'll want to also include zEXt here, but that
  // requires zlib
  while (inStream) {
    std::uint32_t blockLen;
    inStream.read((char *)&blockLen, sizeof(blockLen));
    if (inStream.fail()) {
      throw FileParseException("error when reading from PNG");
    }
    // PNG is big endian, make sure we handle the order correctly
    blockLen = EndianSwapBytes<BIG_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(blockLen);
    char bytes[4];
    inStream.read(bytes, 4);
    if (inStream.fail()) {
      throw FileParseException("error when reading from PNG");
    }
    auto beginBlock = inStream.tellg();
    if (bytes[0] == 'I' && bytes[1] == 'E' && bytes[2] == 'N' &&
        bytes[3] == 'D') {
      break;
    }
#ifndef RDK_USE_BOOST_IOSTREAMS
    bool alreadyWarned = false;
#endif
    if (blockLen > 0 &&
        ((bytes[0] == 't' && bytes[1] == 'E') ||
         (bytes[0] == 'z' && bytes[1] == 'T')) &&
        bytes[2] == 'X' && bytes[3] == 't') {
      // in a tEXt block, read the key:
      std::string key;
      std::getline(inStream, key, '\0');
      if (inStream.fail()) {
        throw FileParseException("error when reading from PNG");
      }
      auto dataLen = blockLen - key.size() - 1;
      std::string value;
      if (bytes[0] == 't') {
        value.resize(dataLen);
        inStream.read(&value.front(), dataLen);
        if (inStream.fail()) {
          throw FileParseException("error when reading from PNG");
        }
      } else if (bytes[0] == 'z') {
#ifdef RDK_USE_BOOST_IOSTREAMS
        value.resize(dataLen);
        inStream.read(&value.front(), dataLen);
        if (inStream.fail()) {
          throw FileParseException("error when reading from PNG");
        }
        value = uncompressString(value.substr(1, dataLen - 1));
#else
        value = "";
        if (!alreadyWarned) {
          BOOST_LOG(rdWarningLog)
              << "compressed metadata found in PNG, but the RDKit was not "
                 "compiled with support for this. Skipping it."
              << std::endl;
          alreadyWarned = true;
        }
#endif
      } else {
        CHECK_INVARIANT(0, "impossible value");
      }
      if (!value.empty()) {
        res.push_back(std::make_pair(key, value));
      }
    }
    inStream.seekg(beginBlock);
    inStream.ignore(blockLen + 4);  // the extra 4 bytes are the CRC
  }

  return res;
};

std::string addMetadataToPNGStream(
    std::istream &inStream,
    const std::vector<std::pair<std::string, std::string>> &metadata,
    bool compressed) {
#ifndef RDK_USE_BOOST_IOSTREAMS
  compressed = false;
#endif
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
  for (const auto &pr : metadata) {
    std::stringstream blk;
    if (!compressed) {
      blk.write("tEXt", 4);
      // write the name along with a zero
      blk.write(pr.first.c_str(), pr.first.size() + 1);
      blk.write(pr.second.c_str(), pr.second.size());
    } else {
#ifdef RDK_USE_BOOST_IOSTREAMS
      blk.write("zTXt", 4);
      // write the name along with a zero
      blk.write(pr.first.c_str(), pr.first.size() + 1);
      // write the compressed data
      // first a zero for the "compression method":
      blk.write("\0", 1);
      auto dest = compressString(pr.second);
      blk.write((const char *)dest.c_str(), dest.size());
#else
      // we shouldn't get here since we disabled compressed at the beginning of
      // the function, but check to be sure
      CHECK_INVARIANT(0, "compression support not enabled");
#endif
    }
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

std::string addMolToPNGStream(const ROMol &mol, std::istream &iStream,
                              bool includePkl, bool includeSmiles,
                              bool includeMol) {
  std::vector<std::pair<std::string, std::string>> metadata;
  if (includePkl) {
    std::string pkl;
    MolPickler::pickleMol(mol, pkl);
    metadata.push_back(std::make_pair(augmentTagName(PNGData::pklTag), pkl));
  }
  if (includeSmiles) {
    std::string smi = MolToCXSmiles(mol);
    metadata.push_back(std::make_pair(augmentTagName(PNGData::smilesTag), smi));
  }
  if (includeMol) {
    bool includeStereo = true;
    int confId = -1;
    bool kekulize = false;
    std::string mb = MolToMolBlock(mol, includeStereo, confId, kekulize);
    metadata.push_back(std::make_pair(augmentTagName(PNGData::molTag), mb));
  }
  return addMetadataToPNGStream(iStream, metadata);
};

ROMol *PNGStreamToMol(std::istream &inStream,
                      const SmilesParserParams &params) {
  ROMol *res = nullptr;
  auto metadata = PNGStreamToMetadata(inStream);
  bool formatFound = false;
  for (const auto &pr : metadata) {
    if (boost::starts_with(pr.first, PNGData::pklTag)) {
      res = new ROMol(pr.second);
      formatFound = true;
    } else if (boost::starts_with(pr.first, PNGData::smilesTag)) {
      res = SmilesToMol(pr.second, params);
      formatFound = true;
    } else if (boost::starts_with(pr.first, PNGData::molTag)) {
      res = MolBlockToMol(pr.second, params.sanitize, params.removeHs);
      formatFound = true;
    }
    if (formatFound) {
      break;
    }
  }
  if (!formatFound) {
    throw FileParseException("No suitable metadata found.");
  }
  return res;
}

std::vector<std::unique_ptr<ROMol>> PNGStreamToMols(
    std::istream &inStream, const std::string &tagToUse,
    const SmilesParserParams &params) {
  std::vector<std::unique_ptr<ROMol>> res;
  auto metadata = PNGStreamToMetadata(inStream);
  for (const auto &pr : metadata) {
    if (!boost::starts_with(pr.first, tagToUse)) {
      continue;
    }
    if (boost::starts_with(pr.first, PNGData::pklTag)) {
      res.emplace_back(new ROMol(pr.second));
    } else if (boost::starts_with(pr.first, PNGData::smilesTag)) {
      res.emplace_back(SmilesToMol(pr.second, params));
    } else if (boost::starts_with(pr.first, PNGData::molTag)) {
      res.emplace_back(
          MolBlockToMol(pr.second, params.sanitize, params.removeHs));
    }
  }
  return res;
}

}  // namespace RDKit
