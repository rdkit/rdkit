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
#include "MultithreadedSDMolSupplier.h"

#include "FileParserUtils.h"

namespace RDKit {
namespace v2 {
namespace FileParsers {
MultithreadedSDMolSupplier::MultithreadedSDMolSupplier(
    const std::string &fileName, const Parameters &params,
    const MolFileParserParams &parseParams) {
  dp_inStream = openAndCheckStream(fileName);
  initFromSettings(true, params, parseParams);
  POSTCONDITION(dp_inStream, "bad instream");
}

MultithreadedSDMolSupplier::MultithreadedSDMolSupplier(
    std::istream *inStream, bool takeOwnership, const Parameters &params,
    const MolFileParserParams &parseParams) {
  PRECONDITION(inStream, "bad stream");
  dp_inStream = inStream;
  initFromSettings(takeOwnership, params, parseParams);
  POSTCONDITION(dp_inStream, "bad instream");
}

MultithreadedSDMolSupplier::MultithreadedSDMolSupplier() {
  dp_inStream = nullptr;
  initFromSettings(false, d_params, d_parseParams);
}

void MultithreadedSDMolSupplier::initFromSettings(
    bool takeOwnership, const Parameters &params,
    const MolFileParserParams &parseParams) {
  df_owner = takeOwnership;
  d_params = params;
  d_parseParams = parseParams;
  d_params.numWriterThreads = getNumThreadsToUse(params.numWriterThreads);
  d_inputQueue.reset(
      new ConcurrentQueue<std::tuple<std::string, unsigned int, unsigned int>>(
          d_params.sizeInputQueue));
  d_outputQueue.reset(
      new ConcurrentQueue<std::tuple<RWMol *, std::string, unsigned int>>(
          d_params.sizeOutputQueue));

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

  /*
    // we are not at the end of file, check for blank lines
    unsigned int numEmpty = 0;
    std::string tempStr;
    // in case df_end is not set then, reset file pointer
    std::streampos holder = dp_inStream->tellg();
          if(static_cast<long int>(holder) == -1){ std::cerr << "putan\n";
    return;} for (unsigned int i = 0; i < 4; i++) { tempStr =
    getLine(dp_inStream); if (dp_inStream->eof()) { df_end = true; break;
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
  */
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
         ((prevStr.find_first_not_of(" \t\r\n") != std::string::npos &&
           prevStr.find("M  END") != 0) ||
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
  index = d_currentRecordId;
  ++d_currentRecordId;
  return true;
}

void MultithreadedSDMolSupplier::readMolProps(RWMol &mol,
                                              std::istringstream &inStream) {
  PRECONDITION(inStream, "no stream");
  bool hasProp = false;
  bool warningIssued = false;
  std::string tempStr;
  std::string dlabel = "";
  std::getline(inStream, tempStr);

  // FIX: report files missing the $$$$ marker
  while (!inStream.eof() && !inStream.fail() &&
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
          std::getline(inStream, tempStr);
          std::string stmp = strip(tempStr);
          while (stmp.length() != 0) {
            std::getline(inStream, tempStr);
            if (inStream.eof()) {
              throw FileParseException("End of data field name not found");
            }
          }
        } else {
          dlabel = tempStr.substr(sl + 1, se - sl - 1);
          // we know the label - now read in the relevant properties
          // until we hit a blank line
          std::getline(inStream, tempStr);

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
            // erase tempStr in case the file does not end with a carrier
            // return (we will end up in an infinite loop if we don't do
            // this and we do not check for EOF in this while loop body)
            tempStr.erase();
            std::getline(inStream, tempStr);
            stmp = strip(tempStr);
          }
          mol.setProp(dlabel, prop);
          if (df_processPropertyLists) {
            // apply this as an atom property list if that's appropriate
            FileParserUtils::processMolPropertyList(mol, dlabel);
          }
        }
      } else {
        if (d_parseParams.strictParsing) {
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
                  << "Property <" << dlabel
                  << "> will be truncated after the first blank line\n";
            } else {
              BOOST_LOG(rdWarningLog) << "Spurious data before the first "
                                         "property will be ignored\n";
            }
            warningIssued = true;
          }
        }
      }
    }
    std::getline(inStream, tempStr);
  }
}

RWMol *MultithreadedSDMolSupplier::processMoleculeRecord(
    const std::string &record, unsigned int lineNum) {
  PRECONDITION(dp_inStream, "no stream");
  std::istringstream inStream(record);
  auto res =
      v2::FileParsers::MolFromMolDataStream(inStream, lineNum, d_parseParams);
  if (res) {
    this->readMolProps(*res, inStream);
    return res.release();
  }
  return nullptr;
}
}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
#endif
