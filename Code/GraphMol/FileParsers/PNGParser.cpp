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
#include <GraphMol/Chirality.h>
#include <vector>

#include "FileParsers.h"
#include <RDGeneral/BoostStartInclude.h>
#include <boost/crc.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#ifdef RDK_USE_BOOST_IOSTREAMS
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#endif
#include <RDGeneral/BoostEndInclude.h>
#if !defined(RDK_USE_BOOST_IOSTREAMS) && defined(RDK_USE_STANDALONE_ZLIB)
#ifndef ZLIB_CONST
#define ZLIB_CONST
#endif
#include <zlib.h>
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

#if defined(RDK_USE_BOOST_IOSTREAMS)
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
#elif defined(RDK_USE_STANDALONE_ZLIB)
std::string zlibActOnString(const std::string &inText, bool compress) {
  static const char *deflatePrefix = "de";
  static const char *inflatePrefix = "in";
  const char *zlibActionPrefix;
  int (*zlibAction)(z_streamp, int);
  int (*zlibEnd)(z_streamp);
  int zRetCode;
  std::string res;
  constexpr uInt BUF_SIZE = 65536;
  std::vector<char> outBuf(BUF_SIZE);
  z_stream zs{};
  zs.next_in = reinterpret_cast<const Bytef *>(inText.c_str());
  zs.avail_in = static_cast<uInt>(inText.size());
  zs.next_out = reinterpret_cast<Bytef *>(outBuf.data());
  zs.avail_out = static_cast<uInt>(outBuf.size());
  if (compress) {
    zlibActionPrefix = deflatePrefix;
    zlibAction = deflate;
    zlibEnd = deflateEnd;
    zRetCode = deflateInit(&zs, Z_DEFAULT_COMPRESSION);
  } else {
    zlibActionPrefix = inflatePrefix;
    zlibAction = inflate;
    zlibEnd = inflateEnd;
    zRetCode = inflateInit(&zs);
  }
  if (zRetCode != Z_OK) {
    BOOST_LOG(rdWarningLog)
        << "Failed to initialize zlib stream (" << zRetCode << ")";
    if (zs.msg) {
      BOOST_LOG(rdWarningLog) << ": " << zs.msg;
    }
    BOOST_LOG(rdWarningLog) << "." << std::endl;
    zlibEnd(&zs);
    return "";
  }
  while (zRetCode == Z_OK) {
    zRetCode = zlibAction(&zs, zs.avail_in ? Z_NO_FLUSH : Z_FINISH);
    if (zRetCode != Z_OK && zRetCode != Z_STREAM_END) {
      BOOST_LOG(rdWarningLog) << "Failed to " << zlibActionPrefix
                              << "flate zlib stream (" << zRetCode << ")";
      if (zs.msg) {
        BOOST_LOG(rdWarningLog) << ": " << zs.msg;
      }
      BOOST_LOG(rdWarningLog) << "." << std::endl;
      zlibEnd(&zs);
      return "";
    }
    if (!zs.avail_out) {
      res += std::string(outBuf.data(), BUF_SIZE);
      zs.next_out = reinterpret_cast<Bytef *>(outBuf.data());
      zs.avail_out = BUF_SIZE;
    }
  }
  auto residual = zs.total_out - res.size();
  if (residual) {
    res += std::string(outBuf.data(), residual);
  }
  zlibEnd(&zs);
  return res;
}

std::string uncompressString(const std::string &ztext) {
  return zlibActOnString(ztext, false);
}
std::string compressString(const std::string &text) {
  return zlibActOnString(text, true);
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
#if !defined(RDK_USE_BOOST_IOSTREAMS) && !defined(RDK_USE_STANDALONE_ZLIB)
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
      } else {
#if defined(RDK_USE_BOOST_IOSTREAMS) || defined(RDK_USE_STANDALONE_ZLIB)
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
      }
      if (!value.empty()) {
        res.emplace_back(key, value);
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
#if !defined(RDK_USE_BOOST_IOSTREAMS) && !defined(RDK_USE_STANDALONE_ZLIB)
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
  for (const auto &[key, value] : metadata) {
    std::stringstream blk;
    if (!compressed) {
      blk.write("tEXt", 4);
      // write the name along with a zero
      blk.write(key.c_str(), key.size() + 1);
      blk.write(value.c_str(), value.size());
    } else {
#if defined(RDK_USE_BOOST_IOSTREAMS) || defined(RDK_USE_STANDALONE_ZLIB)
      blk.write("zTXt", 4);
      // write the name along with a zero
      blk.write(key.c_str(), key.size() + 1);
      // write the compressed data
      // first a zero for the "compression method":
      blk.write("\0", 1);
      auto dest = compressString(value);
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
                              const PNGMetadataParams &params) {
  std::vector<std::pair<std::string, std::string>> metadata;
  if (params.includePkl) {
    std::string pkl;
    MolPickler::pickleMol(mol, pkl, params.propertyFlags);
    metadata.emplace_back(augmentTagName(PNGData::pklTag), pkl);
  }
  if (params.includeSmiles) {
    std::string smi =
        MolToCXSmiles(mol, params.smilesWriteParams, params.cxSmilesFlags,
                      params.restoreBondDirs);
    metadata.emplace_back(augmentTagName(PNGData::smilesTag), smi);
  }
  if (params.includeMol) {
    std::unique_ptr<ROMol> molOrigWedging;
    const ROMol *molRef = &mol;
    bool includeStereo = true;
    int confId = -1;
    bool kekulize = false;
    if (params.restoreBondDirs == RestoreBondDirOptionTrue) {
      molOrigWedging.reset(new ROMol(mol));
      Chirality::reapplyMolBlockWedging(*molOrigWedging);
      molRef = molOrigWedging.get();
    }
    std::string mb = MolToMolBlock(*molRef, includeStereo, confId, kekulize);
    metadata.emplace_back(augmentTagName(PNGData::molTag), mb);
  }
  return addMetadataToPNGStream(iStream, metadata);
};

ROMol *PNGStreamToMol(std::istream &inStream,
                      const SmilesParserParams &params) {
  ROMol *res = nullptr;
  auto metadata = PNGStreamToMetadata(inStream);
  bool formatFound = false;
  for (const auto &[key, value] : metadata) {
    if (boost::starts_with(key, PNGData::pklTag)) {
      res = new ROMol(value);
      formatFound = true;
    } else if (boost::starts_with(key, PNGData::smilesTag)) {
      res = SmilesToMol(value, params);
      formatFound = true;
    } else if (boost::starts_with(key, PNGData::molTag)) {
      res = MolBlockToMol(value, params.sanitize, params.removeHs);
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
  for (const auto &[key, value] : metadata) {
    if (!boost::starts_with(key, tagToUse)) {
      continue;
    }
    if (boost::starts_with(key, PNGData::pklTag)) {
      res.emplace_back(new ROMol(value));
    } else if (boost::starts_with(key, PNGData::smilesTag)) {
      res.emplace_back(SmilesToMol(value, params));
    } else if (boost::starts_with(key, PNGData::molTag)) {
      res.emplace_back(MolBlockToMol(value, params.sanitize, params.removeHs));
    }
  }
  return res;
}

}  // namespace RDKit
