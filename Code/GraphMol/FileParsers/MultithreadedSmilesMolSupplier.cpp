#ifdef RDK_BUILD_THREADSAFE_SSS
//
//  Copyright (C) 2020 Shrey Aryan
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "MultithreadedSmilesMolSupplier.h"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/tokenizer.hpp>

#include <sstream>
#include <string>
#include <vector>

using tokenizer = boost::tokenizer<boost::char_separator<char>>;

namespace RDKit {

inline static bool lineIsEmptyOrComment(const std::string &line) {
  return line.empty() || line[0] == '#' ||
         line.find_first_not_of(" \t\r\n") == std::string::npos;
}

namespace v2 {
namespace FileParsers {
MultithreadedSmilesMolSupplier::MultithreadedSmilesMolSupplier(
    const std::string &fileName, const Parameters &params,
    const SmilesMolSupplierParams &parseParams) {
  dp_inStream = openAndCheckStream(fileName);
  CHECK_INVARIANT(dp_inStream, "bad instream");
  CHECK_INVARIANT(!(dp_inStream->eof()), "early EOF");
  // set df_takeOwnership = true
  initFromSettings(true, params, parseParams);
  skipComments();
  if (d_parseParams.titleLine) {
    processTitleLine();
  }
}

MultithreadedSmilesMolSupplier::MultithreadedSmilesMolSupplier(
    std::istream *inStream, bool takeOwnership, const Parameters &params,
    const SmilesMolSupplierParams &parseParams) {
  CHECK_INVARIANT(inStream, "bad instream");
  CHECK_INVARIANT(!(inStream->eof()), "early EOF");
  dp_inStream = inStream;
  initFromSettings(takeOwnership, params, parseParams);
  skipComments();
  if (d_parseParams.titleLine) {
    processTitleLine();
  }
}

MultithreadedSmilesMolSupplier::MultithreadedSmilesMolSupplier() {
  dp_inStream = nullptr;
  initFromSettings(true, d_params, d_parseParams);
}

void MultithreadedSmilesMolSupplier::initFromSettings(
    bool takeOwnership, const Parameters &params,
    const SmilesMolSupplierParams &parseParams) {
  MultithreadedMolSupplier::initFromSettings(takeOwnership, params);
  d_parseParams = parseParams;
  d_line = -1;
}

void MultithreadedSmilesMolSupplier::skipComments() {
  PRECONDITION(dp_inStream, "bad stream");
  if (this->atEnd()) {
    return;
  }

  // Peek at the first character of each line until we we find one
  // that is not a comment. We can't look at the full line and then
  // rewind, because the stream might not support it.
  std::string tempStr;
  while (!dp_inStream->eof() && !dp_inStream->fail()) {
    auto peek = dp_inStream->peek();
    if (dp_inStream->eof()) {
      // No next line, this is the end of the file
      df_readerDone = true;
      break;
    }
    if (peek != '#' && peek != '\r' && peek != '\n') {
      // Not a comment or empty line, leave
      break;
    }
    // Consume the comment line
    if (!std::getline(*dp_inStream, tempStr)) {
      // shouldn´t happen, since we were able to peek
      df_readerDone = true;
      break;
    }
    ++d_line;
  }
}

// --------------------------------------------------
//
//  Reads and processes the title line
//
void MultithreadedSmilesMolSupplier::processTitleLine() {
  PRECONDITION(dp_inStream, "bad stream");
  if (this->atEnd()) {
    // We expected some data, but didn't get any!
    df_eofHitOnRead = true;
    return;
  }

  // loop until we get a valid line. Since we already
  // called skipComments() before this, we shouldn't see
  // any comments
  std::string tempStr;
  while (!dp_inStream->eof() && !dp_inStream->fail() &&
         lineIsEmptyOrComment(tempStr)) {
    tempStr = getLine(dp_inStream);
    ++d_line;
  }
  if (tempStr.empty()) {
    // We expected some data, but didn't get any!
    df_eofHitOnRead = true;
    return;
  }

  boost::char_separator<char> sep(d_parseParams.delimiter.c_str(), "",
                                  boost::keep_empty_tokens);

  auto tokens = tokenizer(tempStr, sep);
  d_props.assign(tokens.begin(), tokens.end());
  for (auto &pname : d_props) {
    boost::trim_if(pname, boost::is_any_of(" \t\r\n"));
  }
}

bool MultithreadedSmilesMolSupplier::extractNextRecord(std::string &record,
                                                       unsigned int &lineNum,
                                                       unsigned int &index) {
  PRECONDITION(dp_inStream, "bad stream");
  if (dp_inStream->eof()) {
    df_readerDone = true;
    return false;
  }

  record.clear();
  std::string tempStr;
  while (!dp_inStream->eof() && !dp_inStream->fail() &&
         lineIsEmptyOrComment(tempStr)) {
    if (!std::getline(*dp_inStream, tempStr)) {
      // We expected some data, but didn't get any!
      break;
    };
    ++d_line;
  }

  // SmilesMolSupplier skips comments and blank lines and does not expose an
  // extra null record when none remain.
  if (lineIsEmptyOrComment(tempStr) &&
      (dp_inStream->eof() || dp_inStream->fail())) {
    if (d_lastReadRecordId == 0) {
      df_eofHitOnRead = true;
    }
    df_readerDone = true;
    return false;
  }

  record = tempStr;
  lineNum = d_line;
  ++d_lastReadRecordId;
  index = d_lastReadRecordId;
  return true;
}

std::unique_ptr<RWMol> MultithreadedSmilesMolSupplier::processMoleculeRecord(
    const std::string &record, unsigned int lineNum) {
  // -----------
  // tokenize the input line:
  // -----------
  boost::char_separator<char> sep(d_parseParams.delimiter.c_str(), "",
                                  boost::keep_empty_tokens);

  auto tokens = tokenizer(record, sep);
  std::vector<std::string> recs(tokens.begin(), tokens.end());

  if (recs.size() <= static_cast<unsigned int>(d_parseParams.smilesColumn)) {
    std::ostringstream errout;
    errout << "ERROR: line #" << lineNum << "does not contain enough tokens\n";
    throw FileParseException(errout.str());
  }

  for (auto &rec : recs) {
    boost::trim_if(rec, boost::is_any_of(" \t\r\n"));
  }

  // -----------
  // get the smiles and create a molecule
  // -----------
  auto res = MolFromSmiles(recs[d_parseParams.smilesColumn],
                           d_parseParams.parseParameters);
  if (!res) {
    std::stringstream errout;
    errout << "Cannot create molecule from : '"
           << recs[d_parseParams.smilesColumn] << "'";
    throw SmilesParseException(errout.str());
  }

  // -----------
  // get the name (if there's a name column)
  // -----------
  if (d_parseParams.nameColumn == -1) {
    // if no name defaults it to the line number we read it from string
    res->setProp(common_properties::_Name, std::to_string(lineNum));
  } else {
    if (d_parseParams.nameColumn >= static_cast<int>(recs.size())) {
      BOOST_LOG(rdWarningLog)
          << "WARNING: no name column found on line " << lineNum << std::endl;
    } else {
      res->setProp(common_properties::_Name, recs[d_parseParams.nameColumn]);
    }
  }

  // -----------
  // read in the properties
  // -----------
  for (unsigned int col = 0; col < recs.size(); col++) {
    if (static_cast<int>(col) == d_parseParams.smilesColumn ||
        static_cast<int>(col) == d_parseParams.nameColumn) {
      continue;
    }
    if (d_props.size() > col && !d_props[col].empty()) {
      res->setProp(d_props[col], recs[col]);
    } else {
      std::string pname = "Column_";
      pname += std::to_string(col);
      res->setProp(pname, recs[col]);
    }
  }
  return res;
}

}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
#endif
