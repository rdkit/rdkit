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
bool validate(const std::string fileFormat,
              const std::string compressionFormat) {
  //! Case 1: No compression format
  if (compressionFormat.empty()) {
    return std::find(supportedFileFormats.begin(), supportedFileFormats.end(),
                     fileFormat) != supportedFileFormats.end();
  }
  //! Case 2: With compression format

  bool flagFileFormat =
      std::find(supportedFileFormats.begin(), supportedFileFormats.end(),
                fileFormat) != supportedFileFormats.end();
  bool flagCompressionFormat =
      std::find(supportedCompressionFormats.begin(),
                supportedCompressionFormats.end(),
                compressionFormat) != supportedCompressionFormats.end();

  return flagFileFormat && flagCompressionFormat;
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
    fname += path.substr(slash2 + 1);
  } else if (slash1 != std::string::npos && slash2 == std::string::npos) {
    fname += path.substr(slash1 + 1);
  } else if (slash1 != std::string::npos && slash2 != std::string::npos) {
    fname += path.substr(std::max(slash1, slash2) + 1);
  } else {
    //! in this case we assume that the path name is the filename
    fname += path;
  }

  return fname;
}

//! given file path determines the file and compression format
void determineFormat(const std::string path, std::string& fileFormat,
                     std::string& compressionFormat) {
	//! filename without compression format
	std::string basename;  
	//! Special case maegz.
	//! NOTE: also supporting case-insensitive filesystems
	if(boost::algorithm::ends_with(path, ".maegz")){
		fileFormat = "mae";
		compressionFormat = "gz";
		return;
	}
	else if(boost::algorithm::ends_with(path, ".gz")){
		compressionFormat = "gz";
		basename = path.substr(0, path.size() - 3);
	} else if( boost::algorithm::ends_with(path, ".zst") ||
      boost::algorithm::ends_with(path, ".bz2") ||
      boost::algorithm::ends_with(path, ".7z")){
		throw std::invalid_argument("Unsupported compression extension");	
	} else {
 		basename = path;
		compressionFormat = "";
	}
	for (auto const& suffix: supportedFileFormats) {
		 if (boost::algorithm::ends_with(basename, "." + suffix)) {
				 fileFormat = suffix;
				 return;
		 }
	}
	throw std::invalid_argument("Unsupported structure or compression extension");
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
