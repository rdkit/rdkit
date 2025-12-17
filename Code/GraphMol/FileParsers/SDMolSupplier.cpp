//
//  Copyright (C) 2002-2022 Greg Landrum and other RDKit contributors
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

#include <boost/algorithm/string.hpp>
#include "MolSupplier.h"
#include "FileParsers.h"

#include <fstream>
#include <sstream>
#include <string>

namespace RDKit {
namespace v2 {
namespace FileParsers {

SDMolSupplier::SDMolSupplier(const std::string &fileName,
                             const MolFileParserParams &params) {
  init();
  dp_inStream = openAndCheckStream(fileName);
  df_owner = true;
  d_molpos.push_back(dp_inStream->tellg());
  d_params = params;
  this->checkForEnd();
  if (df_end) {
    // checkForEnd() sets d_len if we're at EOF. undo that (was GitHub issue
    // 19):
    d_len = 0;
  }
  POSTCONDITION(dp_inStream, "bad instream");
}

SDMolSupplier::SDMolSupplier(std::istream *inStream, bool takeOwnership,
                             const MolFileParserParams &params) {
  PRECONDITION(inStream, "bad stream");
  init();
  dp_inStream = inStream;
  df_owner = takeOwnership;
  d_molpos.push_back(dp_inStream->tellg());
  d_params = params;
  this->checkForEnd();
  if (df_end) {
    // checkForEnd() sets d_len if we're at EOF. undo that (was GitHub issue
    // 19):
    d_len = 0;
  }
  POSTCONDITION(dp_inStream, "bad instream");
}

void SDMolSupplier::init() {
  ForwardSDMolSupplier::init();
  d_len = -1;
  d_last = 0;
}

void SDMolSupplier::setData(const std::string &text) {
  if (dp_inStream && df_owner) {
    delete dp_inStream;
  }
  init();
  std::istream *tmpStream = nullptr;
  tmpStream = static_cast<std::istream *>(
      new std::istringstream(text, std::ios_base::binary));
  dp_inStream = tmpStream;
  df_owner = true;
  d_molpos.push_back(dp_inStream->tellg());
  this->checkForEnd();
  if (df_end) {
    // checkForEnd() sets d_len if we're at EOF. undo that (was GitHub issue
    // 19):
    d_len = 0;
  }
  POSTCONDITION(dp_inStream, "bad instream");
}

void SDMolSupplier::setData(const std::string &text,
                            const MolFileParserParams &params) {
  d_params = params;
  setData(text);
}

void SDMolSupplier::checkForEnd() {
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

void SDMolSupplier::peekCheckForEnd(char* bufPtr, char* bufEnd, std::streampos molStartPos) {
  PRECONDITION(dp_inStream, "no stream");
  int emptyLines = 0;
  char* p = bufPtr;
  
  while (p < bufEnd) {
      if (!std::isspace(*p)) {
          return; 
      }
      if (*p == '\n') {
          emptyLines++;
          if (emptyLines >= 4) { // the 4th empty line found
              this->df_end = true;
              this->d_len = rdcast<int>(this->d_molpos.size());
              return;
          }
      }
      p++;
  }

  // buffer was exhausted without finding 4 empty lines or data. Need to check the stream.
  std::streampos saveBlockPos = dp_inStream->tellg();
  
  dp_inStream->clear();
  dp_inStream->seekg(molStartPos);
  
  // run the standard (slow) logic
  this->checkForEnd();
  
  // restore the stream position
  dp_inStream->clear();
  dp_inStream->seekg(saveBlockPos);
}

void SDMolSupplier::reset() {
  PRECONDITION(dp_inStream, "no stream");
  dp_inStream->clear();
  dp_inStream->seekg(0, std::ios::beg);
  df_end = false;
  d_last = 0;
  d_line = 0;
}

std::unique_ptr<RWMol> SDMolSupplier::next() {
  PRECONDITION(dp_inStream, "no stream");
  if (df_end && d_last >= d_len) {
    throw FileParseException("EOF hit.");
  }

  // set the stream to the current position
  dp_inStream->seekg(d_molpos[d_last]);

  std::string tempStr;
  // finally if we reached the end of the file set end to be true
  if (dp_inStream->eof()) {
    // FIX: we should probably be throwing an exception here
    df_end = true;
    d_len = rdcast<int>(d_molpos.size());
    return nullptr;
  }

  auto res = _next();

  ++d_last;
  std::streampos posHold = dp_inStream->tellg();
  this->checkForEnd();
  if (!this->df_end && d_last >= static_cast<int>(d_molpos.size())) {
    d_molpos.push_back(posHold);
  }

  return res;
}

std::string SDMolSupplier::getItemText(unsigned int idx) {
  PRECONDITION(dp_inStream, "no stream");
  unsigned int holder = d_last;
  moveTo(idx);
  std::streampos begP = d_molpos[idx];
  std::streampos endP;
  try {
    moveTo(idx + 1);
    endP = d_molpos[idx + 1];
  } catch (FileParseException &) {
    dp_inStream->clear();
    dp_inStream->seekg(0, std::ios_base::end);
    endP = dp_inStream->tellg();
  }
  d_last = holder;
  auto *buff = new char[endP - begP];
  dp_inStream->seekg(begP);
  dp_inStream->read(buff, endP - begP);
  std::string res(buff, endP - begP);
  delete[] buff;
  return res;
}

void SDMolSupplier::moveTo(unsigned int idx) {
  PRECONDITION(dp_inStream, "no stream");

  // dp_inStream->seekg() is called for all idx values
  // and earlier calls to next() may have put the stream into a bad state
  dp_inStream->clear();

  // move until we hit the desired idx
  if (idx < d_molpos.size()) {
    dp_inStream->seekg(d_molpos[idx]);
    d_last = idx;
  } 
  // actually scan with buffering
  else {
    dp_inStream->seekg(d_molpos.back());
    d_last = rdcast<int>(d_molpos.size()) - 1;

    const size_t CHUNK_SIZE = 65536;
    const size_t OVERLAP = 4;
    std::vector<char> buffer(CHUNK_SIZE + OVERLAP);
    std::fill(buffer.begin(), buffer.begin() + OVERLAP, '\n'); //safe init
    std::streampos currentStreamPos = dp_inStream->tellg();
    bool foundTarget = false;

    while (dp_inStream->good() && !foundTarget) {
      dp_inStream->read(&buffer[OVERLAP], CHUNK_SIZE);
      std::streamsize bytesRead = dp_inStream->gcount();
      if (bytesRead == 0) break;

      char *bufStart = &buffer[0];
      char *bufEnd = bufStart + OVERLAP + bytesRead;
      char *ptr = bufStart + 1;

      while (true) {
        static const char dollarSigns[] = "$$$$";
        auto match = std::search(ptr, bufEnd, dollarSigns, dollarSigns + 4);
        if (match == bufEnd) break;
        if (*(match - 1) == '\n') {
          char *nlPos = match + 4;
          while (nlPos < bufEnd && *nlPos != '\n') ++nlPos;
          if (nlPos < bufEnd) ++nlPos;
          std::streampos posHold = currentStreamPos + std::streamoff(nlPos - bufStart - OVERLAP);
          bool atTrueEOF = (bytesRead < (std::streamsize)CHUNK_SIZE) && (nlPos >= bufEnd);
          if (!atTrueEOF) {
            this->peekCheckForEnd(nlPos, bufEnd, posHold);//the optimized peek version
            if (!this->df_end) {
              d_molpos.push_back(posHold);
              d_last++;
              if (static_cast<unsigned int>(d_last) == idx) { //not really needed but this way we only index as much as needed
                  foundTarget = true;
                  break; 
              }
            }
          }
        }
        ptr = match + 4;
      }
      
      if (foundTarget) break;

      if (bytesRead >= OVERLAP) std::memcpy(&buffer[0], bufEnd - OVERLAP, OVERLAP);
      currentStreamPos += bytesRead;
    }

    if (foundTarget) {
        dp_inStream->clear();
        dp_inStream->seekg(d_molpos[idx]);
        d_last = idx;
    } else {
        // if we reached end of file without reaching "idx" we have an index error
        d_len = rdcast<int>(d_molpos.size());
        std::ostringstream errout;
        errout << "ERROR: Index error (idx = " << idx << ") : "
               << " we do not have enough mol blocks";
        throw FileParseException(errout.str());
    }
  }
}

std::unique_ptr<RWMol> SDMolSupplier::operator[](unsigned int idx) {
  PRECONDITION(dp_inStream, "no stream");
  // get the molecule with index idx
  moveTo(idx);
  return next();
}

unsigned int SDMolSupplier::length() {
  PRECONDITION(dp_inStream, "no stream");
  // return the number of mol blocks in the sdfile
  if (d_len > 0 || (df_end && d_len == 0)) {
    return d_len;
  } else {
    d_len = rdcast<int>(d_molpos.size());
    dp_inStream->seekg(d_molpos.back());
    const size_t CHUNK_SIZE = 65536;
    const size_t OVERLAP = 4; // to catch "$$$$" at chunk boundaries ("...\n$$ <new chunk> $$...")
    std::vector<char> buffer(CHUNK_SIZE + OVERLAP);
    std::fill(buffer.begin(), buffer.begin() + OVERLAP, '\n'); //safe init
    std::streampos currentStreamPos = dp_inStream->tellg();

    while (dp_inStream->good()) {
      dp_inStream->read(&buffer[OVERLAP], CHUNK_SIZE);
      std::streamsize bytesRead = dp_inStream->gcount();
      if (bytesRead == 0) break; // EOF
      char *bufStart = &buffer[0];
      char *bufEnd = bufStart + OVERLAP + bytesRead;
      char *ptr = bufStart + 1;

      while (true) {
        static const char dollarSigns[] = "$$$$";
        auto match = std::search(ptr, bufEnd, dollarSigns, dollarSigns + 4);
        if (match == bufEnd) break;
        if (*(match - 1) == '\n') { //ensure $$$$ is at start of line
          char *nlPos = match + 4;
          while (nlPos < bufEnd && *nlPos != '\n') ++nlPos;
          if (nlPos < bufEnd) ++nlPos; //move past newline
          std::streampos posHold = currentStreamPos + std::streamoff(nlPos - bufStart - OVERLAP);
          bool atTrueEOF = (bytesRead < (std::streamsize)CHUNK_SIZE) && (nlPos >= bufEnd);
          if (!atTrueEOF) {
            this->peekCheckForEnd(nlPos, bufEnd, posHold);
            if (!this->df_end) {
              d_molpos.push_back(posHold);
              ++d_len;
            }
          }
        }
        ptr = match + 4;
      }
      if (bytesRead >= OVERLAP) std::memcpy(&buffer[0], bufEnd - OVERLAP, OVERLAP);
      currentStreamPos += bytesRead;
    }
    // now remember to set the stream to the last position we want to read
    dp_inStream->clear();
    dp_inStream->seekg(d_molpos[d_last]);
    df_end = false;
    return d_len;
  }
}

bool SDMolSupplier::atEnd() {
  PRECONDITION(dp_inStream, "no stream");
  return df_end;
}

void SDMolSupplier::setStreamIndices(const std::vector<std::streampos> &locs) {
  d_molpos.clear();
  d_molpos.resize(locs.size());
  std::copy(locs.begin(), locs.end(), d_molpos.begin());
  this->reset();
  d_len = rdcast<int>(d_molpos.size());
}
}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
