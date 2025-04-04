//
//  Copyright (C) 2020 Shrey Aryan
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef GENERAL_FILE_READER_H
#define GENERAL_FILE_READER_H
#include <RDGeneral/BadFileException.h>
#include <RDStreams/streams.h>

#include <boost/algorithm/string.hpp>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "MolSupplier.h"
#include "MultithreadedSDMolSupplier.h"
#include "MultithreadedSmilesMolSupplier.h"

namespace RDKit {
namespace FileParsers = v2::FileParsers;
namespace GeneralMolSupplier {
struct SupplierOptions {
  bool takeOwnership = true;
  bool sanitize = true;
  bool removeHs = true;
  bool strictParsing = true;

  std::string delimiter = "\t";
  int smilesColumn = 0;
  int nameColumn = 1;
  bool titleLine = true;

  std::string nameRecord = "";
  int confId2D = -1;
  int confId3D = 0;

  int numWriterThreads = 0;
};
//! current supported file formats
const std::vector<std::string> supportedFileFormats{
    "sdf", "mae", "maegz", "sdfgz", "smi", "csv", "txt", "tsv", "tdt"};
//! current supported compression formats
const std::vector<std::string> supportedCompressionFormats{"gz"};

//! given file path determines the file and compression format
//! returns true on success, otherwise false
//! Note: Error handeling is done in the getSupplier method

inline void determineFormat(const std::string path, std::string &fileFormat,
                            std::string &compressionFormat) {
  //! filename without compression format
  std::string basename;
  //! Special case maegz.
  //! NOTE: also supporting case-insensitive filesystems
  if (boost::algorithm::iends_with(path, ".maegz")) {
    fileFormat = "mae";
    compressionFormat = "gz";
    return;
  } else if (boost::algorithm::iends_with(path, ".sdfgz")) {
    fileFormat = "sdf";
    compressionFormat = "gz";
    return;
  } else if (boost::algorithm::iends_with(path, ".gz")) {
    compressionFormat = "gz";
    basename = path.substr(0, path.size() - 3);
  } else if (boost::algorithm::iends_with(path, ".zst") ||
             boost::algorithm::iends_with(path, ".bz2") ||
             boost::algorithm::iends_with(path, ".7z")) {
    throw BadFileException(
        "Unsupported compression extension (.zst, .bz2, .7z) given path: " +
        path);
  } else {
    basename = path;
    compressionFormat = "";
  }
  for (auto const &suffix : supportedFileFormats) {
    if (boost::algorithm::iends_with(basename, "." + suffix)) {
      fileFormat = suffix;
      return;
    }
  }
  throw BadFileException(
      "Unsupported structure or compression extension given path: " + path);
}

//! returns a new MolSupplier object based on the file name instantiated
//! with the relevant options provided in the SupplierOptions struct
/*!
    <b>Note:</b>
      - the caller owns the memory and therefore the pointer must be deleted
*/

inline std::unique_ptr<FileParsers::MolSupplier> getSupplier(
    const std::string &path, const struct SupplierOptions &opt) {
  std::string fileFormat = "";
  std::string compressionFormat = "";
  //! get the file and compression format form the path
  determineFormat(path, fileFormat, compressionFormat);

  std::istream *strm;
  if (compressionFormat.empty()) {
    strm = new std::ifstream(path.c_str(), std::ios::in | std::ios::binary);
  } else {
#ifdef RDK_USE_BOOST_IOSTREAMS
    strm = new gzstream(path);
#else
    throw BadFileException(
        "compressed files are only supported if the RDKit is built with boost::iostreams support");
#endif
  }

  if ((!(*strm)) || strm->bad()) {
    std::ostringstream errout;
    errout << "Bad input file " << path;
    delete strm;
    throw BadFileException(errout.str());
  }
  strm->peek();
  if (strm->bad() || strm->eof()) {
    std::ostringstream errout;
    errout << "Invalid input file " << path;
    delete strm;
    throw BadFileException(errout.str());
  }

#ifdef RDK_BUILD_THREADSAFE_SSS
  FileParsers::MultithreadedMolSupplier::Parameters params;
  params.numWriterThreads = getNumThreadsToUse(opt.numWriterThreads);
#endif
  //! Dispatch to the appropriate supplier
  if (fileFormat == "sdf") {
    FileParsers::MolFileParserParams parseParams;
    parseParams.sanitize = opt.sanitize;
    parseParams.removeHs = opt.removeHs;
    parseParams.strictParsing = opt.strictParsing;
#ifdef RDK_BUILD_THREADSAFE_SSS
    if (params.numWriterThreads > 1) {
      return std::make_unique<FileParsers::MultithreadedSDMolSupplier>(
          strm, true, params, parseParams);
    }
#endif
    return std::make_unique<FileParsers::ForwardSDMolSupplier>(strm, true,
                                                               parseParams);
  }

  else if (fileFormat == "smi" || fileFormat == "csv" || fileFormat == "txt" ||
           fileFormat == "tsv") {
    FileParsers::SmilesMolSupplierParams parseParams;
    parseParams.delimiter = opt.delimiter;
    parseParams.smilesColumn = opt.smilesColumn;
    parseParams.nameColumn = opt.nameColumn;
    parseParams.titleLine = opt.titleLine;
    parseParams.parseParameters.sanitize = opt.sanitize;
#ifdef RDK_BUILD_THREADSAFE_SSS
    if (params.numWriterThreads > 1) {
      return std::make_unique<FileParsers::MultithreadedSmilesMolSupplier>(
          strm, true, params, parseParams);
    }
#endif
    return std::make_unique<FileParsers::SmilesMolSupplier>(strm, true,
                                                            parseParams);
  }
#ifdef RDK_BUILD_MAEPARSER_SUPPORT
  else if (fileFormat == "mae") {
    FileParsers::MaeMolSupplierParams parseParams;
    parseParams.sanitize = opt.sanitize;
    parseParams.removeHs = opt.removeHs;
    return std::make_unique<FileParsers::MaeMolSupplier>(strm, true,
                                                         parseParams);
  }
#endif
  else if (fileFormat == "tdt") {
    FileParsers::TDTMolSupplierParams parseParams;
    parseParams.nameRecord = opt.nameRecord;
    parseParams.confId2D = opt.confId2D;
    parseParams.confId3D = opt.confId3D;
    parseParams.parseParameters.sanitize = opt.sanitize;
    return std::make_unique<FileParsers::TDTMolSupplier>(strm, true,
                                                         parseParams);
  }
  throw BadFileException("Unsupported file format: " + fileFormat);
}

}  // namespace GeneralMolSupplier
}  // namespace RDKit
#endif
