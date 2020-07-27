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

#include <iomanip>
#include <mutex>
#include <sstream>
namespace RDKit {
//! method for thread-safe printing, only for debugging
struct PrintThread : public std::stringstream {
  static inline std::mutex cout_mutex;
  ~PrintThread() {
    std::lock_guard<std::mutex> l{cout_mutex};
    std::cout << rdbuf();
  }
};

MultithreadedSmilesMolSupplier::MultithreadedSmilesMolSupplier(
    const std::string &fileName, const std::string &delimiter, int smilesColumn,
    int nameColumn, bool titleLine, bool sanitize, int numWriterThreads,
    size_t sizeInputQueue, size_t sizeOutputQueue) {
  init();
  dp_inStream = openAndCheckStream(fileName);
  CHECK_INVARIANT(dp_inStream, "bad instream");
  CHECK_INVARIANT(!(dp_inStream->eof()), "early EOF");

  df_owner = false;
  if (numWriterThreads == -1) {
    d_numWriterThreads = (int)getNumThreadsToUse(numWriterThreads);
  } else {
    d_numWriterThreads =
        std::min(numWriterThreads, (int)getNumThreadsToUse(numWriterThreads));
  }

  d_sizeInputQueue = sizeInputQueue;
  d_sizeOutputQueue = sizeOutputQueue;
  d_inputQueue = new ConcurrentQueue<std::tuple<std::string, unsigned int>>(
      d_sizeInputQueue);
  d_outputQueue = new ConcurrentQueue<ROMol *>(d_sizeOutputQueue);

  d_delim = delimiter;
  df_sanitize = sanitize;
  df_title = titleLine;
  d_smi = smilesColumn;
  d_name = nameColumn;
  df_end = false;

  this->checkForEnd();
  POSTCONDITION(dp_inStream, "bad instream");
}

MultithreadedSmilesMolSupplier::MultithreadedSmilesMolSupplier(
    std::istream *inStream, bool takeOwnership, const std::string &delimiter,
    int smilesColumn, int nameColumn, bool titleLine, bool sanitize,
    int numWriterThreads, size_t sizeInputQueue, size_t sizeOutputQueue) {
  init();
  CHECK_INVARIANT(inStream, "bad instream");
  CHECK_INVARIANT(!(inStream->eof()), "early EOF");
  dp_inStream = inStream;
  df_owner = takeOwnership;
  if (numWriterThreads == -1) {
    d_numWriterThreads = (int)getNumThreadsToUse(numWriterThreads);
  } else {
    d_numWriterThreads =
        std::min(numWriterThreads, (int)getNumThreadsToUse(numWriterThreads));
  }
  d_sizeInputQueue = sizeInputQueue;
  d_sizeOutputQueue = sizeOutputQueue;
  d_inputQueue = new ConcurrentQueue<std::tuple<std::string, unsigned int>>(
      d_sizeInputQueue);
  d_outputQueue = new ConcurrentQueue<ROMol *>(d_sizeOutputQueue);

  d_delim = delimiter;
  df_sanitize = sanitize;
  df_title = titleLine;
  d_smi = smilesColumn;
  d_name = nameColumn;
  df_end = false;
  this->checkForEnd();
  POSTCONDITION(dp_inStream, "bad instream");
}

MultithreadedSmilesMolSupplier::MultithreadedSmilesMolSupplier() { init(); }

MultithreadedSmilesMolSupplier::~MultithreadedSmilesMolSupplier() {
  if (df_owner && dp_inStream) {
    delete dp_inStream;
  }
}

void MultithreadedSmilesMolSupplier::init() {
  dp_inStream = nullptr;
  df_owner = true;
  d_numWriterThreads = 1;
  d_sizeInputQueue = 10;
  d_sizeOutputQueue = 10;
  d_inputQueue = new ConcurrentQueue<std::tuple<std::string, unsigned int>>(
      d_sizeInputQueue);
  d_outputQueue = new ConcurrentQueue<ROMol *>(d_sizeOutputQueue);

  dp_inStream = nullptr;
  df_owner = true;
  df_end = false;
  d_next = -1;
  d_line = -1;
}

// --------------------------------------------------
//
//  Returns the position of the beginning of the next
//  non-comment line in the input stream. -1 is returned if
//  no line could be read;
//
//  Side-effects:
//    - If EOF is hit without finding a valid line, the df_end
//      flag will be set.
//    - Our d_line counter is incremented for each line read
//
long int MultithreadedSmilesMolSupplier::skipComments() {
  PRECONDITION(dp_inStream, "bad stream");
  if (this->getEnd()) {
    return -1;
  }

  std::streampos prev = dp_inStream->tellg();
  std::string tempStr = this->nextLine();
  if (!df_end) {
    // if we didn't immediately hit EOF, loop until we get a valid line:
    while ((tempStr[0] == '#') || (strip(tempStr).size() == 0)) {
      prev = dp_inStream->tellg();
      tempStr = this->nextLine();
      if (this->getEnd()) {
        break;
      }
    }
  }
  // if we hit EOF without getting a proper line, return -1:
  if (tempStr.empty() || (tempStr[0] == '#') || (strip(tempStr).size() == 0)) {
    return -1;
  }
  return static_cast<long int>(prev);
}

// ensures that there is a line available to be read
// from the file:
void MultithreadedSmilesMolSupplier::checkForEnd() {
  PRECONDITION(dp_inStream, "no stream");
  int pos = this->skipComments();
  if (pos != -1) {
    d_line = -1;
    dp_inStream->seekg(0);
    df_end = false;
  }
}

bool MultithreadedSmilesMolSupplier::getEnd() { return df_end; }

// --------------------------------------------------
//
//  Returns true if EOF is hit else false.
//
//  Side-effects:
//    - If EOF is hit without reading anything, the df_end
//      flag will be set.
//    - If a real line is read, our d_line counter is
//      incremented
//
// --------------------------------------------------
std::string MultithreadedSmilesMolSupplier::nextLine() {
  PRECONDITION(dp_inStream, "bad stream");
  if (df_end) {
    return "";
  }
  std::string record = getLine(dp_inStream);
  if (record == "") {
    // got an empty string, check to see if we hit EOF:
    if (dp_inStream->eof() || dp_inStream->bad()) {
      // yes, set our flag:
      df_end = true;
    }
  } else if (dp_inStream->eof()) {
    // we got some data before hitting EOF. So clear the
    // flag on inStream
    dp_inStream->clear();
  }
  d_line++;
  return record;
}

bool MultithreadedSmilesMolSupplier::extractNextRecord(std::string &record,
                                                       unsigned int &lineNum) {
  PRECONDITION(dp_inStream, "bad stream");
  if (this->getEnd()) {
    return -1;
  }

  std::streampos prev = dp_inStream->tellg();
  std::string tempStr = this->nextLine();
  if (!df_end) {
    // if we didn't immediately hit EOF, loop until we get a valid line:
    while ((tempStr[0] == '#') || (strip(tempStr).size() == 0)) {
      prev = dp_inStream->tellg();
      tempStr = this->nextLine();
      if (this->getEnd()) {
        break;
      }
    }
  }
  record = tempStr;
  lineNum = d_line;
  // if we hit EOF without getting a proper line, return -1:
  if (tempStr.empty() || (tempStr[0] == '#') || (strip(tempStr).size() == 0)) {
    return false;
  }
  return true;
}

ROMol *MultithreadedSmilesMolSupplier::processMoleculeRecord(
    const std::string &record, unsigned int lineNum) {
  ROMol *res = nullptr;
  try {
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
      errout << "ERROR: line #" << lineNum
             << "does not contain enough tokens\n";
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
    unsigned int iprop = 0;
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
      iprop++;
    }

  } catch (const SmilesParseException &pe) {
    // Couldn't parse the passed in smiles
    // Simply print out a message
    BOOST_LOG(rdErrorLog) << "ERROR: Smiles parse error on line " << lineNum
                          << "\n";
    BOOST_LOG(rdErrorLog) << "ERROR: " << pe.what() << "\n";
    res = nullptr;
  } catch (const MolSanitizeException &se) {
    // We couldn't sanitize the molecule
    //  write out an error message
    BOOST_LOG(rdErrorLog) << "ERROR: Could not sanitize molecule on line "
                          << lineNum << std::endl;
    BOOST_LOG(rdErrorLog) << "ERROR: " << se.what() << "\n";
    res = nullptr;
  } catch (...) {
    //  write out an error message
    BOOST_LOG(rdErrorLog) << "ERROR: Could not process molecule on line "
                          << lineNum << std::endl;
    res = nullptr;
  }

  return res;
}

}  // namespace RDKit
