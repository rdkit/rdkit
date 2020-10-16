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

  unsigned int numWriterThreads = 0;
};
//! current supported file formats
const std::vector<std::string> supportedFileFormats{
    "sdf", "mae", "maegz", "sdfgz", "smi", "csv", "txt", "tsv", "tdt"};
//! current supported compression formats
const std::vector<std::string> supportedCompressionFormats{"gz"};

//! given file path determines the file and compression format
//! returns true on success, otherwise false
//! Note: Error handeling is done in the getSupplier method

void determineFormat(const std::string path, std::string& fileFormat,
                     std::string& compressionFormat) {
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
  for (auto const& suffix : supportedFileFormats) {
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

std::unique_ptr<MolSupplier> getSupplier(const std::string& path,
                                         const struct SupplierOptions& opt) {
  std::string fileFormat = "";
  std::string compressionFormat = "";
  //! get the file and compression format form the path
  determineFormat(path, fileFormat, compressionFormat);

  std::istream* strm;
  if (compressionFormat.empty()) {
    strm = new std::ifstream(path.c_str());
  } else {
    strm = new gzstream(path);
  }

  //! Dispatch to the appropriate supplier
  if (fileFormat == "sdf") {
#ifdef RDK_THREADSAFE_SSS
    if (opt.numWriterThreads > 0) {
      MultithreadedSDMolSupplier* sdsup = new MultithreadedSDMolSupplier(
          strm, true, opt.sanitize, opt.removeHs, opt.strictParsing,
          opt.numWriterThreads);
      std::unique_ptr<MolSupplier> p(sdsup);
      return p;
    }
#endif
    ForwardSDMolSupplier* sdsup = new ForwardSDMolSupplier(
        strm, true, opt.sanitize, opt.removeHs, opt.strictParsing);
    std::unique_ptr<MolSupplier> p(sdsup);
    return p;
  }

  else if (fileFormat == "smi" || fileFormat == "csv" || fileFormat == "txt" ||
           fileFormat == "tsv") {
#ifdef RDK_THREADSAFE_SSS
    if (opt.numWriterThreads > 0) {
      MultithreadedSmilesMolSupplier* smsup =
          new MultithreadedSmilesMolSupplier(
              strm, true, opt.delimiter, opt.smilesColumn, opt.nameColumn,
              opt.titleLine, opt.sanitize, opt.numWriterThreads);
      std::unique_ptr<MolSupplier> p(smsup);
      return p;
    }
#endif
    SmilesMolSupplier* smsup =
        new SmilesMolSupplier(strm, true, opt.delimiter, opt.smilesColumn,
                              opt.nameColumn, opt.titleLine, opt.sanitize);
    std::unique_ptr<MolSupplier> p(smsup);
    return p;
  }

  else if (fileFormat == "mae") {
    MaeMolSupplier* maesup =
        new MaeMolSupplier(strm, true, opt.sanitize, opt.removeHs);
    std::unique_ptr<MolSupplier> p(maesup);
    return p;
  }

  else if (fileFormat == "tdt") {
    TDTMolSupplier* tdtsup = new TDTMolSupplier(
        strm, true, opt.nameRecord, opt.confId2D, opt.confId3D, opt.sanitize);
    std::unique_ptr<MolSupplier> p(tdtsup);
    return p;
  }
  throw BadFileException("Unsupported fileFormat: " + fileFormat);
}

}  // namespace GeneralMolSupplier
}  // namespace RDKit
#endif
