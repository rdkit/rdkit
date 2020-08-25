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

//! given file path determines the file and compression format
void determineFormat(const std::string path, std::string& fileFormat,
                     std::string& compressionFormat) {
	//! filename without compression format
	std::string basename;  
	//! Special case maegz.
	//! NOTE: also supporting case-insensitive filesystems
	if(boost::algorithm::iends_with(path, ".maegz")){
		fileFormat = "mae";
		compressionFormat = "gz";
		return;
	}
	else if(boost::algorithm::iends_with(path, ".gz")){
		compressionFormat = "gz";
		basename = path.substr(0, path.size() - 3);
	} else if( boost::algorithm::iends_with(path, ".zst") ||
      boost::algorithm::iends_with(path, ".bz2") ||
      boost::algorithm::iends_with(path, ".7z")){
		throw std::invalid_argument("Unsupported compression extension");	
	} else {
 		basename = path;
		compressionFormat = "";
	}
	for (auto const& suffix: supportedFileFormats) {
		 if (boost::algorithm::iends_with(basename, "." + suffix)) {
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

  //! Dispatch to the appropriate supplier
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
