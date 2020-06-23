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
#include <string>
#include <iostream>
#include <vector>
#include "MolSupplier.h"
#include <RDGeneral/BadFileException.h>
#include <RDStreams/streams.h>

namespace RDKit {
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

bool valid(std::string& fileFormat, std::string& compressionFormat) {
  //! current formats supported
  std::vector<std::string> fileFormats{"sdf", "mae", "smi", "csv",
                                       "txt", "tsv", "tdt"};
  std::vector<std::string> compressionFormats{"gz"};

  //! Case 1: Unconventional format types
  if (fileFormat.compare("maegz") == 0) {
    fileFormat = "mae";
    compressionFormat = "gz";
    return true;
  }

  //! Case 2: Either the filename has a file format or compression format or
  //! both
  bool flag_fileFormat = std::find(fileFormats.begin(), fileFormats.end(),
                                   fileFormat) != fileFormats.end();

  if (!compressionFormat.empty()) {
    bool flag_compressionFormat =
        std::find(compressionFormats.begin(), compressionFormats.end(),
                  compressionFormat) != compressionFormats.end();

    //! if the compression type is not valid then
    if (!flag_compressionFormat) {
      fileFormat = compressionFormat;
      compressionFormat = "";
      flag_fileFormat = std::find(fileFormats.begin(), fileFormats.end(),
                                  fileFormat) != fileFormats.end();
    }
  }
  return flag_fileFormat;
}

std::string getFileName(const std::string path) {
  char delimiter = '/';
  std::string fname = "";
  auto slash = path.rfind(delimiter, path.length());
  if (slash != std::string::npos) {
    fname += path.substr(slash + 1, path.length() - slash);
  }
  return fname;
}

void determineFormat(const std::string path, std::string& fileFormat,
                     std::string& compressionFormat) {
  std::string fileName = getFileName(path);
  int dots = std::count(fileName.begin(), fileName.end(), '.');

  if (dots == 0)
    throw std::invalid_argument(
        "Recieved Invalid File Format, no extension or compression");

  else if (dots == 1) {
    //! there is a file format but no compression format
    int pos = fileName.find(".");
    fileFormat = fileName.substr(pos + 1);
    if (!valid(fileFormat, compressionFormat))
      throw std::invalid_argument("Recieved Invalid File Format");
  } else {
    //! there is a file and compression format
    int n = fileName.length();
    int p1 = fileName.rfind(".");
    int p2 = fileName.rfind(".", p1 - 1);
    fileFormat = fileName.substr(p2 + 1, (p1 - p2) - 1);
    //! possible compression format
    compressionFormat = fileName.substr(p1 + 1, (n - p1) + 1);
    if (!valid(fileFormat, compressionFormat))
      throw std::invalid_argument(
          "Recieved Invalid File or Compression Format");
  }
}

MolSupplier* getSupplier(const std::string& path,
                         const struct SupplierOptions opt) {
  std::string fileFormat = "";
  std::string compressionFormat = "";
  determineFormat(
      path, fileFormat,
      compressionFormat);  //! get the file and compression format form the path

  if (compressionFormat.empty()) {
    std::ifstream* strm = new std::ifstream(path.c_str());
    if (fileFormat.compare("sdf") == 0) {
      ForwardSDMolSupplier* sdsup =
          new ForwardSDMolSupplier(strm, opt.takeOwnership, opt.sanitize,
                                   opt.removeHs, opt.strictParsing);
      return sdsup;
    }

    else if (fileFormat == "smi") {
      SmilesMolSupplier* smsup = new SmilesMolSupplier(
          strm, opt.takeOwnership, opt.delimiter, opt.smilesColumn,
          opt.nameColumn, opt.titleLine, opt.sanitize);
      return smsup;
    }

    else if (fileFormat == "csv") {
      SmilesMolSupplier* smsup = new SmilesMolSupplier(
          strm, opt.takeOwnership, opt.delimiter, opt.smilesColumn,
          opt.nameColumn, opt.titleLine, opt.sanitize);
      return smsup;
    }

    else if (fileFormat == "txt") {
      SmilesMolSupplier* smsup = new SmilesMolSupplier(
          strm, opt.takeOwnership, opt.delimiter, opt.smilesColumn,
          opt.nameColumn, opt.titleLine, opt.sanitize);
      return smsup;
    }

    else if (fileFormat == "tsv") {
      SmilesMolSupplier* smsup = new SmilesMolSupplier(
          strm, opt.takeOwnership, opt.delimiter, opt.smilesColumn,
          opt.nameColumn, opt.titleLine, opt.sanitize);
      return smsup;
    }

    else if (fileFormat == "mae") {
      MaeMolSupplier* maesup = new MaeMolSupplier(strm, opt.takeOwnership,
                                                  opt.sanitize, opt.removeHs);
      return maesup;
    }

    else if (fileFormat == "tdt") {
      TDTMolSupplier* tdtsup =
          new TDTMolSupplier(strm, opt.takeOwnership, opt.nameRecord,
                             opt.confId2D, opt.confId3D, opt.sanitize);
      return tdtsup;
    }
  } else {
    auto* strm = new gzstream(path);

    if (fileFormat == "sdf") {
      ForwardSDMolSupplier* sdsup =
          new ForwardSDMolSupplier(strm, opt.takeOwnership, opt.sanitize,
                                   opt.removeHs, opt.strictParsing);
      return sdsup;
    }

    else if (fileFormat == "smi") {
      SmilesMolSupplier* smsup = new SmilesMolSupplier(
          strm, opt.takeOwnership, opt.delimiter, opt.smilesColumn,
          opt.nameColumn, opt.titleLine, opt.sanitize);
      return smsup;
    }

    else if (fileFormat == "csv") {
      SmilesMolSupplier* smsup = new SmilesMolSupplier(
          strm, opt.takeOwnership, opt.delimiter, opt.smilesColumn,
          opt.nameColumn, opt.titleLine, opt.sanitize);
      return smsup;
    }

    else if (fileFormat == "txt") {
      SmilesMolSupplier* smsup = new SmilesMolSupplier(
          strm, opt.takeOwnership, opt.delimiter, opt.smilesColumn,
          opt.nameColumn, opt.titleLine, opt.sanitize);
      return smsup;
    }

    else if (fileFormat == "tsv") {
      SmilesMolSupplier* smsup = new SmilesMolSupplier(
          strm, opt.takeOwnership, opt.delimiter, opt.smilesColumn,
          opt.nameColumn, opt.titleLine, opt.sanitize);
      return smsup;
    }

    else if (fileFormat == "mae") {
      MaeMolSupplier* maesup = new MaeMolSupplier(strm, opt.takeOwnership,
                                                  opt.sanitize, opt.removeHs);
      return maesup;
    }

    else if (fileFormat == "tdt") {
      TDTMolSupplier* tdtsup =
          new TDTMolSupplier(strm, opt.takeOwnership, opt.nameRecord,
                             opt.confId2D, opt.confId3D, opt.sanitize);
      return tdtsup;
    }
  }
}
}  // namespace RDKit
#endif
