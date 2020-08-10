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

namespace RDKit {
struct PrintThread : public std::stringstream {
  static inline std::mutex cout_mutex;
  ~PrintThread() {
    std::lock_guard<std::mutex> l{cout_mutex};
    std::cout << rdbuf();
  }
};

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
  d_lastMolPos = dp_inStream->tellg();
}

MultithreadedSDMolSupplier::~MultithreadedSDMolSupplier() {
  if (df_owner && dp_inStream) {
    delete dp_inStream;
    df_owner = false;
    dp_inStream = nullptr;
  }
}

// ensures that there is a line available to be read
// from the file:
void MultithreadedSDMolSupplier::checkForEnd() {
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

bool MultithreadedSDMolSupplier::getEnd() const { return df_end; }

bool MultithreadedSDMolSupplier::extractNextRecord(std::string &record,
                                                   unsigned int &lineNum,
                                                   unsigned int &index) {
  PRECONDITION(dp_inStream, "no stream");
  if (dp_inStream->eof()) {
    df_end = true;
    return false;
  }
  record = "";
  std::string tempStr;
  dp_inStream->seekg(d_lastMolPos);
  lineNum = d_line;
  while (!dp_inStream->eof() && !dp_inStream->fail() &&
         (tempStr[0] != '$' || tempStr.substr(0, 4) != "$$$$")) {
    tempStr = getLine(dp_inStream);
    record += tempStr + "\n";
    if (tempStr[0] == '$' && tempStr.substr(0, 4) == "$$$$") {
      std::streampos posHold = dp_inStream->tellg();
      this->checkForEnd();
      if (!this->df_end) {
        d_lastMolPos = posHold;
      }
    }
    ++d_line;
  }

  index = d_currentRecordId;
  ++d_currentRecordId;
  return true;
}

ROMol *MultithreadedSDMolSupplier::processMoleculeRecord(
    const std::string &record, unsigned int lineNum) {
  PRECONDITION(dp_inStream, "no stream");
  auto res = MolBlockToMol(record, df_sanitize, df_removeHs, df_strictParsing);
  return res;
}

}  // namespace RDKit
