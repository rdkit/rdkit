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

#include <iostream>
#include <string>
#include <vector>

#include "MolSupplier.h"

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
};
//! current supported file formats
std::vector<std::string> supportedFileFormats{"sdf", "mae", "maegz", "smi",
                                              "csv", "txt", "tsv",   "tdt"};
//! current supported compression formats
std::vector<std::string> supportedCompressionFormats{"gz"};

//! determines whether file and compression format are valid
bool validate(std::string& fileFormat, std::string& compressionFormat) {
  bool flag_fileFormat = true;
  bool flag_compressionFormat = true;

  //! Case 1: No compression format
  if (compressionFormat.empty()) {
    //! Unconventional format types
    if (fileFormat.compare("maegz") == 0) {
      fileFormat = "mae";
      compressionFormat = "gz";
      return true;
    }

    //! Usual file formats
    return std::find(supportedFileFormats.begin(), supportedFileFormats.end(),
                     fileFormat) != supportedFileFormats.end();
  }
  //! Case 2: With compression format
  else {
    flag_fileFormat =
        std::find(supportedFileFormats.begin(), supportedFileFormats.end(),
                  fileFormat) != supportedFileFormats.end();
    flag_compressionFormat =
        std::find(supportedCompressionFormats.begin(),
                  supportedCompressionFormats.end(),
                  compressionFormat) != supportedCompressionFormats.end();

    //! Case A: Valid file format and valid compression format
    if (flag_fileFormat && flag_compressionFormat) {
      return true;
    }

    //! Case B, C: Valid file format but not valid compression format
    //! 						or invalid file format
    //! and valid compression format
    else if ((flag_fileFormat && !flag_compressionFormat) ||
             (!flag_fileFormat && flag_compressionFormat)) {
      return false;
    }

    //! Case D: Invalid file format and invalid compression format
    else {
      //! it is possible that we have read a file with name *.2.txt (for
      //! example) thus fileformat = "2" and compressionformat = "txt", so we
      //! try to set the file format as the compression format and the
      //! compression format as an empty string. Then we check of validity.
      fileFormat = compressionFormat;
      compressionFormat = "";
      return std::find(supportedFileFormats.begin(), supportedFileFormats.end(),
                       fileFormat) != supportedFileFormats.end();
    }
  }
}

//! returns the file name givens the file path
std::string getFileName(const std::string path) {
  char delimiter = '/';
  char delimiter_win = '\\';
  auto n = path.length();
  std::string fname = "";

  auto slash1 = path.rfind(delimiter, n);
  auto slash2 = path.rfind(delimiter_win, n);
  if (slash1 == std::string::npos && slash2 != std::string::npos) {
    fname += path.substr(slash2 + 1, n - slash1);
  } else if (slash1 != std::string::npos && slash2 == std::string::npos) {
    fname += path.substr(slash1 + 1, n - slash1);
  } else if (slash1 != std::string::npos && slash2 != std::string::npos) {
    if (n - slash1 > n - slash2) {
      fname += path.substr(slash2 + 1, n - slash2);
    } else {
      fname += path.substr(slash1 + 1, n - slash1);
    }
  } else {
    throw std::invalid_argument(
        "Unable to determine filename from path: no back or forward slash "
        "found");
  }

  return fname;
}

//! given file path determines the file and compression format
void determineFormat(const std::string path, std::string& fileFormat,
                     std::string& compressionFormat) {
  std::string fileName = getFileName(path);
  int dots = std::count(fileName.begin(), fileName.end(), '.');

  if (dots == 0) {
    throw std::invalid_argument(
        "Unable to determine file format: no filename extension found");
  }

  else if (dots == 1) {
    //! there is a file format but no compression format
    int pos = fileName.rfind(".");
    fileFormat = fileName.substr(pos + 1);
    if (!validate(fileFormat, compressionFormat)) {
      throw std::invalid_argument(
          "Unable to determine file format: unsupported filename extension");
    }
  } else {
    //! there is a file and compression format
    int n = fileName.length();
    int p1 = fileName.rfind(".");
    int p2 = fileName.rfind(".", p1 - 1);
    fileFormat = fileName.substr(p2 + 1, (p1 - p2) - 1);
    //! possible compression format
    compressionFormat = fileName.substr(p1 + 1, (n - p1) + 1);
    if (!validate(fileFormat, compressionFormat)) {
      throw std::invalid_argument(
          "Unable to determine file format: unsupported extension");
    }
  }
}

//! returns a MolSupplier object based on the file name instantiated
//! with the relevant options provided in the SupplierOptions struct
/*!
    <b>Note:</b>
      - the caller owns the memory and therefore the pointer must be deleted
*/

MolSupplier* getSupplier(const std::string& path,
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

  //! Handle the case when there is no compression format
  if (fileFormat.compare("sdf") == 0) {
    ForwardSDMolSupplier* sdsup = new ForwardSDMolSupplier(
        strm, true, opt.sanitize, opt.removeHs, opt.strictParsing);
    return sdsup;
  }

  else if (fileFormat == "smi") {
    SmilesMolSupplier* smsup =
        new SmilesMolSupplier(strm, true, opt.delimiter, opt.smilesColumn,
                              opt.nameColumn, opt.titleLine, opt.sanitize);
    return smsup;
  }

  else if (fileFormat == "csv") {
    SmilesMolSupplier* smsup =
        new SmilesMolSupplier(strm, true, opt.delimiter, opt.smilesColumn,
                              opt.nameColumn, opt.titleLine, opt.sanitize);
    return smsup;
  }

  else if (fileFormat == "txt") {
    SmilesMolSupplier* smsup =
        new SmilesMolSupplier(strm, true, opt.delimiter, opt.smilesColumn,
                              opt.nameColumn, opt.titleLine, opt.sanitize);
    return smsup;
  }

  else if (fileFormat == "tsv") {
    SmilesMolSupplier* smsup =
        new SmilesMolSupplier(strm, true, opt.delimiter, opt.smilesColumn,
                              opt.nameColumn, opt.titleLine, opt.sanitize);
    return smsup;
  }

  else if (fileFormat == "mae") {
    MaeMolSupplier* maesup =
        new MaeMolSupplier(strm, true, opt.sanitize, opt.removeHs);
    return maesup;
  }

  else if (fileFormat == "tdt") {
    TDTMolSupplier* tdtsup = new TDTMolSupplier(
        strm, true, opt.nameRecord, opt.confId2D, opt.confId3D, opt.sanitize);
    return tdtsup;
  }
}

}  // namespace GeneralMolSupplier
}  // namespace RDKit
#endif
