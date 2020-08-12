//
//  Copyright (C) 2020 Shrey Aryan
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "MultithreadedSDMolSupplier.h"

#include "FileParserUtils.h"

namespace RDKit {
MultithreadedSDMolSupplier::MultithreadedSDMolSupplier(
    const std::string &fileName, bool sanitize, bool removeHs,
    bool strictParsing, unsigned int numWriterThreads, size_t sizeInputQueue,
    size_t sizeOutputQueue) {
  dp_inStream = openAndCheckStream(fileName);
  _init(true, sanitize, removeHs, strictParsing, numWriterThreads,
        sizeInputQueue, sizeOutputQueue);
  POSTCONDITION(dp_inStream, "bad instream");
  startThreads();
}

MultithreadedSDMolSupplier::MultithreadedSDMolSupplier(
    std::istream *inStream, bool takeOwnership, bool sanitize, bool removeHs,
    bool strictParsing, unsigned int numWriterThreads, size_t sizeInputQueue,
    size_t sizeOutputQueue) {
  PRECONDITION(inStream, "bad stream");
  dp_inStream = inStream;
  _init(takeOwnership, sanitize, removeHs, strictParsing, numWriterThreads,
        sizeInputQueue, sizeOutputQueue);
  POSTCONDITION(dp_inStream, "bad instream");
  startThreads();
}

MultithreadedSDMolSupplier::MultithreadedSDMolSupplier() {
  dp_inStream = nullptr;
  _init(false, true, true, true, 2, 5, 5);
  startThreads();
}

void MultithreadedSDMolSupplier::_init(bool takeOwnership, bool sanitize,
                                       bool removeHs, bool strictParsing,
                                       unsigned int numWriterThreads,
                                       size_t sizeInputQueue,
                                       size_t sizeOutputQueue) {
  df_owner = takeOwnership;
  df_sanitize = sanitize;
  df_removeHs = removeHs;
  df_strictParsing = strictParsing;
  d_numWriterThreads = numWriterThreads;
  d_sizeInputQueue = sizeInputQueue;
  d_sizeOutputQueue = sizeOutputQueue;
  d_inputQueue =
      new ConcurrentQueue<std::tuple<std::string, unsigned int, unsigned int>>(
          d_sizeInputQueue);
  d_outputQueue =
      new ConcurrentQueue<std::tuple<ROMol *, std::string, unsigned int>>(
          d_sizeOutputQueue);

  df_end = false;
  d_line = 0;
  df_processPropertyLists = true;
}

MultithreadedSDMolSupplier::~MultithreadedSDMolSupplier() {
  if (df_owner && dp_inStream) {
    delete dp_inStream;
    df_owner = false;
    dp_inStream = nullptr;
  }
}

// ensures that there is a line available to be read
// from the file, implementation identical to the method in
// in ForwardSDMolSupplier
void MultithreadedSDMolSupplier::checkForEnd() {
  PRECONDITION(dp_inStream, "no stream");
  // we will call it end of file if we have more than 4 contiguous empty lines
  // or we reach end of file in the meantime
  if (dp_inStream->eof()) {
    df_end = true;
    return;
  }

  // we are not at the end of file, check for blank lines
  unsigned int numEmpty = 0;
  std::string tempStr;
  // in case df_end is not set then, reset file pointer
  std::streampos holder = dp_inStream->tellg();
  for (unsigned int i = 0; i < 4; i++) {
    tempStr = getLine(dp_inStream);
    if (dp_inStream->eof()) {
      df_end = true;
      break;
    }
    if (tempStr.find_first_not_of(" \t\r\n") == std::string::npos) {
      ++numEmpty;
    }
  }
  if (numEmpty == 4) {
    df_end = true;
  }
  // we need to reset the file pointer to read the next record
  if (!df_end) {
    dp_inStream->clear();
    dp_inStream->seekg(holder);
  }
}

bool MultithreadedSDMolSupplier::getEnd() const {
  PRECONDITION(dp_inStream, "no stream");
  return df_end;
}

bool MultithreadedSDMolSupplier::extractNextRecord(std::string &record,
                                                   unsigned int &lineNum,
                                                   unsigned int &index) {
  PRECONDITION(dp_inStream, "no stream");
  if (dp_inStream->eof()) {
    df_end = true;
    return false;
  }

  std::string currentStr, prevStr;
  record = "";
  lineNum = d_line;
  while (!dp_inStream->eof() && !dp_inStream->fail() &&
         (prevStr.find_first_not_of(" \t\r\n") != std::string::npos ||
          currentStr[0] != '$' || currentStr.substr(0, 4) != "$$$$")) {
    prevStr = currentStr;
    std::getline(*dp_inStream, currentStr);
    record += currentStr + "\n";
    ++d_line;
    if (prevStr.find_first_not_of(" \t\r\n") == std::string::npos &&
        currentStr[0] == '$' && currentStr.substr(0, 4) == "$$$$") {
      this->checkForEnd();
    }
  }
  //	PrintThread{} << "record: \n" << record;
  //	PrintThread{} << "------------------------------------------\n";
  index = d_currentRecordId;
  ++d_currentRecordId;
  return true;
}

ROMol *MultithreadedSDMolSupplier::processMoleculeRecord(
    const std::string &record, unsigned int lineNum) {
  PRECONDITION(dp_inStream, "no stream");
  std::istringstream inStream(record);
  auto res = MolDataStreamToMol(inStream, lineNum, df_sanitize, df_removeHs,
                                df_strictParsing);
  return res;
}

}  // namespace RDKit
