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

namespace RDKit {
MultithreadedSmilesMolSupplier::MultithreadedSmilesMolSupplier(
    const std::string &fileName, const std::string &delimiter, int smilesColumn,
    int nameColumn, bool titleLine, bool sanitize,
    unsigned int numWriterThreads, size_t sizeInputQueue,
    size_t sizeOutputQueue) {
  dp_inStream = openAndCheckStream(fileName);
  CHECK_INVARIANT(dp_inStream, "bad instream");
  CHECK_INVARIANT(!(dp_inStream->eof()), "early EOF");
  // set df_takeOwnership = true
  initFromSettings(true, delimiter, smilesColumn, nameColumn, titleLine,
                   sanitize, numWriterThreads, sizeInputQueue, sizeOutputQueue);
  startThreads();
  POSTCONDITION(dp_inStream, "bad instream");
}

MultithreadedSmilesMolSupplier::MultithreadedSmilesMolSupplier(
    std::istream *inStream, bool takeOwnership, const std::string &delimiter,
    int smilesColumn, int nameColumn, bool titleLine, bool sanitize,
    unsigned int numWriterThreads, size_t sizeInputQueue,
    size_t sizeOutputQueue) {
  CHECK_INVARIANT(inStream, "bad instream");
  CHECK_INVARIANT(!(inStream->eof()), "early EOF");
  dp_inStream = inStream;
  initFromSettings(takeOwnership, delimiter, smilesColumn, nameColumn,
                   titleLine, sanitize, numWriterThreads, sizeInputQueue,
                   sizeOutputQueue);
  startThreads();
  POSTCONDITION(dp_inStream, "bad instream");
}

MultithreadedSmilesMolSupplier::MultithreadedSmilesMolSupplier() {
  dp_inStream = nullptr;
  initFromSettings(true, "", 0, 1, true, true, 1, 5, 5);
  startThreads();
}

MultithreadedSmilesMolSupplier::~MultithreadedSmilesMolSupplier() {
  if (df_owner && dp_inStream) {
    delete dp_inStream;
    df_owner = false;
    dp_inStream = nullptr;
  }
}

void MultithreadedSmilesMolSupplier::initFromSettings(
    bool takeOwnership, const std::string &delimiter, int smilesColumn,
    int nameColumn, bool titleLine, bool sanitize,
    unsigned int numWriterThreads, size_t sizeInputQueue,
    size_t sizeOutputQueue) {
  df_owner = takeOwnership;
  d_delim = delimiter;
  d_smi = smilesColumn;
  d_name = nameColumn;
  df_title = titleLine;
  df_sanitize = sanitize;
  d_numWriterThreads = getNumThreadsToUse(numWriterThreads);
  d_sizeInputQueue = sizeInputQueue;
  d_sizeOutputQueue = sizeOutputQueue;
  d_inputQueue =
      new ConcurrentQueue<std::tuple<std::string, unsigned int, unsigned int>>(
          d_sizeInputQueue);
  d_outputQueue =
      new ConcurrentQueue<std::tuple<ROMol *, std::string, unsigned int>>(
          d_sizeOutputQueue);
  df_end = false;
  d_line = -1;
}

bool MultithreadedSmilesMolSupplier::getEnd() const {
  PRECONDITION(dp_inStream, "no stream");
  return df_end;
}

// --------------------------------------------------
//
//  Reads and processes the title line
//
void MultithreadedSmilesMolSupplier::processTitleLine() {
  PRECONDITION(dp_inStream, "bad stream");
  std::string tempStr = getLine(dp_inStream);
  // loop until we get a valid line
  while (!dp_inStream->eof() && !dp_inStream->fail() &&
         ((tempStr[0] == '#') || (strip(tempStr).size() == 0))) {
    tempStr = getLine(dp_inStream);
  }
  boost::char_separator<char> sep(d_delim.c_str(), "",
                                  boost::keep_empty_tokens);
  tokenizer tokens(tempStr, sep);
  for (tokenizer::iterator tokIter = tokens.begin(); tokIter != tokens.end();
       ++tokIter) {
    std::string pname = strip(*tokIter);
    d_props.push_back(pname);
  }
}

bool MultithreadedSmilesMolSupplier::extractNextRecord(std::string &record,
                                                       unsigned int &lineNum,
                                                       unsigned int &index) {
  PRECONDITION(dp_inStream, "bad stream");
  if (dp_inStream->eof()) {
    df_end = true;
    return false;
  }

  // need to process title line
  // if we have not called next yet and the current record id = 1
  // then we are seeking the first record
  if (d_lastRecordId == 0 && d_currentRecordId == 1) {
    if (df_title) {
      this->processTitleLine();
    }
  }
  std::string tempStr = getLine(dp_inStream);
  record = "";
  while (!dp_inStream->eof() && !dp_inStream->fail() &&
         ((tempStr[0] == '#') || (strip(tempStr).size() == 0))) {
    tempStr = getLine(dp_inStream);
  }

  record = tempStr;
  lineNum = d_line;
  index = d_currentRecordId;
  ++d_currentRecordId;
  return true;
}

ROMol *MultithreadedSmilesMolSupplier::processMoleculeRecord(
    const std::string &record, unsigned int lineNum) {
  ROMol *res = nullptr;

  // -----------
  // tokenize the input line:
  // -----------
  boost::char_separator<char> sep(d_delim.c_str(), "",
                                  boost::keep_empty_tokens);
  tokenizer tokens(record, sep);
  STR_VECT recs;
  for (tokenizer::iterator tokIter = tokens.begin(); tokIter != tokens.end();
       ++tokIter) {
    std::string rec = strip(*tokIter);
    recs.push_back(rec);
  }
  if (recs.size() <= static_cast<unsigned int>(d_smi)) {
    std::ostringstream errout;
    errout << "ERROR: line #" << lineNum << "does not contain enough tokens\n";
    throw FileParseException(errout.str());
  }

  // -----------
  // get the smiles and create a molecule
  // -----------
  SmilesParserParams params;
  params.sanitize = df_sanitize;
  params.allowCXSMILES = false;
  params.parseName = false;
  res = SmilesToMol(recs[d_smi], params);
  if (!res) {
    std::stringstream errout;
    errout << "Cannot create molecule from : '" << recs[d_smi] << "'";
    throw SmilesParseException(errout.str());
  }

  // -----------
  // get the name (if there's a name column)
  // -----------
  if (d_name == -1) {
    // if no name defaults it to the line number we read it from string
    std::ostringstream tstr;
    tstr << lineNum;
    std::string mname = tstr.str();
    res->setProp(common_properties::_Name, mname);
  } else {
    if (d_name >= static_cast<int>(recs.size())) {
      BOOST_LOG(rdWarningLog)
          << "WARNING: no name column found on line " << lineNum << std::endl;
    } else {
      res->setProp(common_properties::_Name, recs[d_name]);
    }
  }

  // -----------
  // read in the properties
  // -----------
  for (unsigned int col = 0; col < recs.size(); col++) {
    if (static_cast<int>(col) == d_smi || static_cast<int>(col) == d_name) {
      continue;
    }
    std::string pname, pval;
    if (d_props.size() > col) {
      pname = d_props[col];
    } else {
      pname = "Column_";
      std::stringstream ss;
      ss << col;
      pname += ss.str();
    }

    pval = recs[col];
    res->setProp(pname, pval);
  }
  return res;
}

}  // namespace RDKit
#endif
