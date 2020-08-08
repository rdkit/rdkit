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
MultithreadedSDMolSupplier::MultithreadedSDMolSupplier(
    const std::string &fileName, bool sanitize, bool removeHs,
    bool strictParsing, unsigned int numWriterThreads, size_t sizeInputQueue,
    size_t sizeOutputQueue) {
  dp_inStream = openAndCheckStream(fileName);
  _init(true, sanitize, removeHs, strictParsing, numWriterThreads,
        sizeInputQueue, sizeOutputQueue);
  d_molpos.push_back(dp_inStream->tellg());
  this->checkForEnd();
  if (df_end) {
    // checkForEnd() sets d_len if we're at EOF. undo that (was GitHub issue
    // 19):
    d_len = 0;
  }
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
  d_molpos.push_back(dp_inStream->tellg());
  this->checkForEnd();
  if (df_end) {
    // checkForEnd() sets d_len if we're at EOF. undo that (was GitHub issue
    // 19):
    d_len = 0;
  }
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
  d_len = -1;
  d_last = 0;
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
    d_len = rdcast<int>(d_molpos.size());
    return;
  }
  // we are not at the end of file, check for blank lines
  unsigned int nempty = 0;
  std::string tempStr, stmp;
  for (unsigned int i = 0; i < 4; i++) {
    tempStr = getLine(dp_inStream);
    if (dp_inStream->eof()) {
      df_end = true;
      d_len = rdcast<int>(d_molpos.size());
      return;
    }
    if (tempStr.find_first_not_of(" \t\r\n") == std::string::npos) {
      ++nempty;
    }
  }
  if (nempty == 4) {
    df_end = true;
    d_len = rdcast<int>(d_molpos.size());
  }
}

bool MultithreadedSDMolSupplier::moveTo(unsigned int idx) {
  PRECONDITION(dp_inStream, "no stream");

  // dp_inStream->seekg() is called for all idx values
  // and earlier calls to next() may have put the stream into a bad state
  dp_inStream->clear();

  // move until we hit the desired idx
  if (idx < d_molpos.size()) {
    dp_inStream->seekg(d_molpos[idx]);
    d_last = idx;
  } else {
    std::string tempStr;
    dp_inStream->seekg(d_molpos.back());
    d_last = rdcast<int>(d_molpos.size()) - 1;
    while (d_last < static_cast<int>(idx) && !dp_inStream->eof() &&
           !dp_inStream->fail()) {
      d_line++;
      tempStr = getLine(dp_inStream);

      if (tempStr[0] == '$' && tempStr.substr(0, 4) == "$$$$") {
        std::streampos posHold = dp_inStream->tellg();
        this->checkForEnd();
        if (!this->df_end) {
          d_molpos.push_back(posHold);
          d_last++;
        }
      }
    }
    // if we reached end of file without reaching "idx" we have an index error
    if (dp_inStream->eof()) {
      d_len = rdcast<int>(d_molpos.size());
      return false;
    }
  }
  return true;
}

bool MultithreadedSDMolSupplier::getEnd() const { return df_end; }

bool MultithreadedSDMolSupplier::extractNextRecord(std::string &record,
                                                   unsigned int &lineNum,
                                                   unsigned int &index) {
  PRECONDITION(dp_inStream, "no stream");
  if (df_end && d_last >= d_len) {
    return false;
  }
  // set the stream to the current position
  dp_inStream->seekg(d_molpos[d_last]);

  // finally if we reached the end of the file set end to be true
  if (dp_inStream->eof()) {
    // FIX: we should probably be throwing an exception here
    df_end = true;
    d_len = rdcast<int>(d_molpos.size());
    return false;
  }

  // get the next record, using the getItemText method
  unsigned int idx = d_last;
  if (moveTo(idx)) {
    std::streampos begP = d_molpos[idx];
    std::streampos endP;
    if (moveTo(idx + 1)) {
      endP = d_molpos[idx + 1];
    } else {
      dp_inStream->clear();
      dp_inStream->seekg(0, std::ios_base::end);
      endP = dp_inStream->tellg();
    }
    d_last = idx;
    auto *buff = new char[endP - begP];
    dp_inStream->seekg(begP);
    dp_inStream->read(buff, endP - begP);
    std::string tempStr(buff, endP - begP);
    record = tempStr;
    index = d_currentRecordId;
    lineNum = d_line;
    ++d_currentRecordId;
    ++d_last;
    std::streampos posHold = dp_inStream->tellg();
    this->checkForEnd();
    if (!this->df_end && d_last >= static_cast<int>(d_molpos.size())) {
      d_molpos.push_back(posHold);
    }
    delete[] buff;
    return true;
  }

  return false;
}

void MultithreadedSDMolSupplier::readMolProps(ROMol *mol) {
  PRECONDITION(dp_inStream, "no stream");
  PRECONDITION(mol, "no molecule");
  d_line++;
  bool hasProp = false;
  bool warningIssued = false;
  std::string tempStr;
  std::string dlabel = "";
  std::getline(*dp_inStream, tempStr);

  // FIX: report files missing the $$$$ marker
  while (!dp_inStream->eof() && !dp_inStream->fail() &&
         (tempStr[0] != '$' || tempStr.substr(0, 4) != "$$$$")) {
    tempStr = strip(tempStr);
    if (tempStr != "") {
      if (tempStr[0] == '>') {  // data header line: start of a data item
        // ignore all other crap and seek for for a data label enclosed
        // by '<' and '>'
        // FIX: "CTfile.pdf" (page 51) says that the data header line does not
        // have to contain a data label (instead can have something line field
        // id into a MACCS db). But we do not currently know what to do in this
        // situation - so ignore such data items for now
        hasProp = true;
        warningIssued = false;
        tempStr.erase(0, 1);            // remove the first ">" sign
        size_t sl = tempStr.find("<");  // begin datalabel
        size_t se = tempStr.find(">");  // end datalabel
        if ((sl == std::string::npos) || (se == std::string::npos) ||
            (se == (sl + 1))) {
          // we either do not have a data label or the label is empty
          // no data label ignore until next data item
          // i.e. until we hit a blank line
          d_line++;
          std::getline(*dp_inStream, tempStr);
          std::string stmp = strip(tempStr);
          while (stmp.length() != 0) {
            d_line++;
            std::getline(*dp_inStream, tempStr);
            if (dp_inStream->eof()) {
              throw FileParseException("End of data field name not found");
            }
          }
        } else {
          dlabel = tempStr.substr(sl + 1, se - sl - 1);
          // we know the label - now read in the relevant properties
          // until we hit a blank line
          d_line++;
          std::getline(*dp_inStream, tempStr);

          std::string prop = "";
          std::string stmp = strip(tempStr);
          int nplines = 0;  // number of lines for this property
          while (stmp.length() != 0 || tempStr[0] == ' ' ||
                 tempStr[0] == '\t') {
            nplines++;
            if (nplines > 1) {
              prop += "\n";
            }
            // take off \r if it's still in the property:
            if (tempStr[tempStr.length() - 1] == '\r') {
              tempStr.erase(tempStr.length() - 1);
            }
            prop += tempStr;
            d_line++;
            // erase tempStr in case the file does not end with a carrier
            // return (we will end up in an infinite loop if we don't do
            // this and we do not check for EOF in this while loop body)
            tempStr.erase();
            std::getline(*dp_inStream, tempStr);
            stmp = strip(tempStr);
          }
          mol->setProp(dlabel, prop);
          if (df_processPropertyLists) {
            // apply this as an atom property list if that's appropriate
            FileParserUtils::processMolPropertyList(*mol, dlabel);
          }
        }
      } else {
        if (df_strictParsing) {
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
    std::getline(*dp_inStream, tempStr);
  }
}

ROMol *MultithreadedSDMolSupplier::processMoleculeRecord(
    const std::string &record, unsigned int lineNum) {
  PRECONDITION(dp_inStream, "no stream");
  auto res = MolBlockToMol(record, df_sanitize, df_removeHs, df_strictParsing);
  return res;
}

}  // namespace RDKit
