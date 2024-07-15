//
//  Copyright (C) 2009-2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/SanitException.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/algorithm/string.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include "MolSupplier.h"
#include "FileParsers.h"
#include "FileParserUtils.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace RDKit {
std::string strip(const std::string &orig) {
  // FIX: this can be more efficient
  // strip the end of line, white spaces and tabs
  std::string res =
      boost::trim_right_copy_if(orig, boost::is_any_of(" \t\r\n"));
  res = boost::trim_left_copy_if(res, boost::is_any_of(" \t\r\n"));
  return res;
}

namespace v2 {
namespace FileParsers {

ForwardSDMolSupplier::ForwardSDMolSupplier(std::istream *inStream,
                                           bool takeOwnership,
                                           const MolFileParserParams &params) {
  PRECONDITION(inStream, "bad stream");
  init();
  dp_inStream = inStream;
  df_owner = takeOwnership;
  d_params = params;
  POSTCONDITION(dp_inStream, "bad instream");
}

void ForwardSDMolSupplier::init() {
  dp_inStream = nullptr;
  df_owner = false;
  df_end = false;
  d_line = 0;
  df_processPropertyLists = true;
}

void ForwardSDMolSupplier::reset() {
  UNDER_CONSTRUCTION("reset() not supported for ForwardSDMolSuppliers();");
}

void ForwardSDMolSupplier::readMolProps(ROMol &mol) {
  PRECONDITION(dp_inStream, "no stream");
  d_line++;
  bool hasProp = false;
  bool warningIssued = false;
  std::string dlabel = "";
  std::string inl;
  std::getline(*dp_inStream, inl);
  std::string_view tempStr = inl;

  // FIX: report files missing the $$$$ marker
  while (!dp_inStream->eof() && !dp_inStream->fail() && 
         (tempStr.empty() || tempStr.at(0) != '$' || tempStr.substr(0, 4) != "$$$$")) {
    tempStr = FileParserUtils::strip(tempStr);
    if (!tempStr.empty()) {
      if (tempStr.at(0) == '>') {  // data header line: start of a data item
        // ignore all other crap and seek for for a data label enclosed
        // by '<' and '>'
        // FIX: "CTfile.pdf" (page 51) says that the data header line does not
        // have to contain a data label (instead can have something line field
        // id into a MACCS db). But we do not currently know what to do in this
        // situation - so ignore such data items for now
        hasProp = true;
        warningIssued = false;
        tempStr = tempStr.substr(1);    // remove the first ">" sign
        size_t sl = tempStr.find("<");  // begin datalabel
        size_t se = tempStr.find(">");  // end datalabel
        if ((sl == std::string::npos) || (se == std::string::npos) ||
            (se == (sl + 1))) {
          // we either do not have a data label or the label is empty
          // no data label ignore until next data item
          // i.e. until we hit a blank line
          d_line++;
          std::getline(*dp_inStream, inl);
          tempStr = inl;
          auto stmp = FileParserUtils::strip(tempStr);
          while (stmp.length() != 0) {
            d_line++;
            std::getline(*dp_inStream, inl);
            tempStr = inl;
            if (dp_inStream->eof()) {
              throw FileParseException("End of data field name not found");
            }
          }
        } else {
          dlabel = tempStr.substr(sl + 1, se - sl - 1);
          // we know the label - now read in the relevant properties
          // until we hit a blank line
          d_line++;
          std::getline(*dp_inStream, inl);
          tempStr = inl;

          std::string prop = "";
          auto stmp = FileParserUtils::strip(tempStr);
          int nplines = 0;  // number of lines for this property
          while (!stmp.empty() || (!tempStr.empty() && (tempStr.at(0) == ' ' ||
							tempStr.at(0) == '\t'))) {
            nplines++;
            if (nplines > 1) {
              prop += "\n";
            }
            // take off \r if it's still in the property:
	    if (!tempStr.empty()) {
	      if (tempStr.back() == '\r') {
		tempStr = tempStr.substr(0, tempStr.size() - 1);
	      }
	      prop += tempStr;
	    }
            d_line++;
            // erase inl in case the file does not end with a carriage
            // return (we will end up in an infinite loop if we don't do
            // this and we do not check for EOF in this while loop body)
            inl.erase();
            std::getline(*dp_inStream, inl);
            tempStr = inl;
	    if (tempStr.empty()) {
	      stmp = tempStr;
	    } else {
	      stmp = FileParserUtils::strip(tempStr);
	    }
          }
          mol.setProp(dlabel, prop);
          if (df_processPropertyLists) {
            // apply this as an atom property list if that's appropriate
            FileParserUtils::processMolPropertyList(mol, dlabel);
          }
        }
      } else {
        if (d_params.strictParsing) {
          // at this point we should always be at a line starting with '>'
          // following a blank line. If this is not true and df_strictParsing
          // is true, then throw an exception, otherwise truncate the rest of
          // the data field following the blank line until the next '>' or EOF
          // and issue a warning
          // FIX: should we be deleting the molecule (which is probably fine)
          // because we couldn't read the data ???
          throw FileParseException("Problems encountered parsing data fields");
        } else {
          if (!warningIssued) {
            if (hasProp) {
              BOOST_LOG(rdWarningLog)
                  << "Property <" << dlabel << "> will be truncated after "
                  << "the first blank line" << std::endl;
            } else {
              BOOST_LOG(rdWarningLog)
                  << "Spurious data before the first property will be "
                     "ignored"
                  << std::endl;
            }
            warningIssued = true;
          }
        }
      }
    }
    d_line++;
    std::getline(*dp_inStream, inl);
    tempStr = inl;
  }
}

std::unique_ptr<RWMol> ForwardSDMolSupplier::next() {
  PRECONDITION(dp_inStream, "no stream");

  if (dp_inStream->eof()) {
    // FIX: we should probably be throwing an exception here
    df_end = true;
    return nullptr;
  }

  return _next();
}

std::unique_ptr<RWMol> ForwardSDMolSupplier::_next() {
  PRECONDITION(dp_inStream, "no stream");

  std::string tempStr;
  std::unique_ptr<RWMol> res;
  if (dp_inStream->eof()) {
    df_end = true;
    return res;
  }

  df_eofHitOnRead = false;
  unsigned int line = d_line;
  try {
    MolFromMolDataStream(*dp_inStream, line, d_params).swap(res);
    // there's a special case when trying to read an empty string that
    // we get an empty molecule after only reading a single line without any
    // additional error state.
    if (!res && dp_inStream->eof() && (line - d_line < 2)) {
      df_eofHitOnRead = true;
    }
    d_line = line;
    if (res) {
      this->readMolProps(*res);
    } else if (!dp_inStream->eof()) {
      // FIX: report files missing the $$$$ marker
      std::getline(*dp_inStream, tempStr);
      ++d_line;
      while (!dp_inStream->eof() && !dp_inStream->fail() &&
             (tempStr.at(0) != '$' || tempStr.substr(0, 4) != "$$$$")) {
        std::getline(*dp_inStream, tempStr);
        ++d_line;
      }
    }
  } catch (FileParseException &fe) {
    if (d_line < static_cast<int>(line)) {
      d_line = line;
    }
    // we couldn't read a mol block or the data for the molecule. In this case
    // advance forward in the stream until we hit the next record and then
    // rethrow
    // the exception. This should allow us to read the next molecule.
    BOOST_LOG(rdErrorLog) << "ERROR: " << fe.what() << std::endl;
    BOOST_LOG(rdErrorLog)
        << "ERROR: moving to the beginning of the next molecule\n";

    // FIX: report files missing the $$$$ marker
    d_line++;
    std::getline(*dp_inStream, tempStr);
    while (!dp_inStream->eof() && !dp_inStream->fail() &&
           (tempStr.empty() || tempStr.at(0) != '$' || tempStr.substr(0, 4) != "$$$$")) {
      d_line++;
      std::getline(*dp_inStream, tempStr);
    }
  } catch (MolSanitizeException &se) {
    if (d_line < static_cast<int>(line)) {
      d_line = line;
    }
    // We couldn't sanitize a molecule we got - write out an error message and
    // move to
    // the beginning of the next molecule
    BOOST_LOG(rdErrorLog)
        << "ERROR: Could not sanitize molecule ending on line " << d_line
        << std::endl;
    BOOST_LOG(rdErrorLog) << "ERROR: " << se.what() << "\n";

    d_line++;
    std::getline(*dp_inStream, tempStr);
    if (dp_inStream->eof()) {
      df_eofHitOnRead = true;
    }
    while (!dp_inStream->eof() && !dp_inStream->fail() &&
           (tempStr.empty() || tempStr.at(0) != '$' || tempStr.substr(0, 4) != "$$$$")) {
      d_line++;
      std::getline(*dp_inStream, tempStr);
    }
  } catch (...) {
    if (dp_inStream->eof()) {
      df_eofHitOnRead = true;
    }
    if (d_line < static_cast<int>(line)) {
      d_line = line;
    }

    BOOST_LOG(rdErrorLog) << "Unexpected error hit on line " << d_line
                          << std::endl;
    BOOST_LOG(rdErrorLog)
        << "ERROR: moving to the beginning of the next molecule\n";
    d_line++;
    std::getline(*dp_inStream, tempStr);
    if (dp_inStream->eof()) {
      df_eofHitOnRead = true;
    }
    while (!dp_inStream->eof() && !dp_inStream->fail() &&
           (tempStr.empty() || tempStr.at(0) != '$' || tempStr.substr(0, 4) != "$$$$")) {
      d_line++;
      std::getline(*dp_inStream, tempStr);
    }
  }
  if (dp_inStream->eof()) {
    // FIX: we should probably be throwing an exception here
    df_end = true;
  }
  return res;
}

void ForwardSDMolSupplier::checkForEnd() {
  PRECONDITION(dp_inStream, "no stream");
  // we will call it end of file if we have more than 4 contiguous empty lines
  // or we reach end of file in the meantime
  if (dp_inStream->eof()) {
    df_end = true;
    return;
  }
  // we are not at the end of file, check for blank lines
  unsigned int nempty = 0;
  std::string tempStr;
  for (unsigned int i = 0; i < 4; i++) {
    tempStr = getLine(dp_inStream);
    if (dp_inStream->eof()) {
      df_end = true;
      return;
    }
    if (tempStr.find_first_not_of(" \t\r\n") == std::string::npos) {
      ++nempty;
    }
  }
  if (nempty == 4) {
    df_end = true;
  }
}

bool ForwardSDMolSupplier::atEnd() {
  PRECONDITION(dp_inStream, "no stream");
  return df_end;
}
}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
