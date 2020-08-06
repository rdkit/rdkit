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
  _init(true, delimiter, smilesColumn, nameColumn, titleLine, sanitize,
        numWriterThreads, sizeInputQueue, sizeOutputQueue);
  this->checkForEnd();
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
  _init(takeOwnership, delimiter, smilesColumn, nameColumn, titleLine, sanitize,
        numWriterThreads, sizeInputQueue, sizeOutputQueue);
  this->checkForEnd();
  startThreads();
  POSTCONDITION(dp_inStream, "bad instream");
}

MultithreadedSmilesMolSupplier::MultithreadedSmilesMolSupplier() {
  dp_inStream = nullptr;
  _init(true, "", 0, 1, true, true, 2, 5, 5);
  startThreads();
}

MultithreadedSmilesMolSupplier::~MultithreadedSmilesMolSupplier() {
  if (df_owner && dp_inStream) {
    delete dp_inStream;
  }
}

void MultithreadedSmilesMolSupplier::_init(bool takeOwnership,
                                           const std::string &delimiter,
                                           int smilesColumn, int nameColumn,
                                           bool titleLine, bool sanitize,
                                           unsigned int numWriterThreads,
                                           size_t sizeInputQueue,
                                           size_t sizeOutputQueue) {
  df_owner = takeOwnership;
  d_delim = delimiter;
  d_smi = smilesColumn;
  d_name = nameColumn;
  df_title = titleLine;
  df_sanitize = sanitize;
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
  d_len = -1;
  d_next = -1;
  d_line = -1;
  d_molpos.clear();
  d_lineNums.clear();
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
// from the file
void MultithreadedSmilesMolSupplier::checkForEnd() {
  PRECONDITION(dp_inStream, "no stream");
  int pos = this->skipComments();
  if (pos != -1) {
    d_line = -1;
    dp_inStream->seekg(0);
    df_end = false;
  }
}

bool MultithreadedSmilesMolSupplier::getEnd() const { return df_end; }

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

// --------------------------------------------------
//
//  Reads and processes the title line
//
void MultithreadedSmilesMolSupplier::processTitleLine() {
  PRECONDITION(dp_inStream, "bad stream");
  int pos = this->skipComments();
  if (pos >= 0) {
    dp_inStream->seekg(pos);

    std::string tempStr = getLine(dp_inStream);
    boost::char_separator<char> sep(d_delim.c_str(), "",
                                    boost::keep_empty_tokens);
    tokenizer tokens(tempStr, sep);
    for (tokenizer::iterator tokIter = tokens.begin(); tokIter != tokens.end();
         ++tokIter) {
      std::string pname = strip(*tokIter);
      d_props.push_back(pname);
    }
  }
}

// --------------------------------------------------
//
//  Moves to the position of a particular entry in the
//  stream.
//
//  If insufficient entries are present, the method returns false
//	instead of throwing an exception
//
bool MultithreadedSmilesMolSupplier::moveTo(unsigned int idx) {
  PRECONDITION(dp_inStream, "bad instream");
  // get the easy situations (boundary conditions) out of the
  // way first:
  if (d_len > -1 && idx >= static_cast<unsigned int>(d_len)) {
    df_end = true;
    return false;
  }

  // dp_inStream->seekg() is called for all idx values
  // and earlier calls to next() may have put the stream into a bad state
  dp_inStream->clear();

  // -----------
  // Case 1: we have already read the particular entry:
  //
  // Set the stream position and return
  // -----------
  if (!d_molpos.empty() && d_molpos.size() > idx) {
    dp_inStream->clear();  // clear the EOF tag if it has been set
    df_end = false;
    dp_inStream->seekg(d_molpos[idx]);
    d_next = idx;
    d_line = d_lineNums[idx];
    return true;
  }

  // -----------
  // Case 2: we haven't read the entry, so move forward until
  //   we've gone far enough.
  // -----------
  if (d_molpos.empty()) {
    // if we are just starting out, process the title line
    dp_inStream->seekg(0);
    if (df_title) {
      this->processTitleLine();
    }
  } else {
    // move to the last position we've seen:
    dp_inStream->seekg(d_molpos.back());
    // read that line:
    std::string tmp = getLine(dp_inStream);
  }

  // the stream pointer is now at the last thing we read in
  while (d_molpos.size() <= idx) {
    int nextP = this->skipComments();
    if (nextP < 0) {
      return false;
    } else {
      d_molpos.emplace_back(nextP);
      d_lineNums.push_back(d_line);
      if (d_molpos.size() == idx + 1 && df_end) {
        // boundary condition: we could read the point we were looking for
        // but not the next one.
        // indicate that we've reached EOF:
        dp_inStream->clear();
        dp_inStream->seekg(0, std::ios_base::end);
        d_len = d_molpos.size();
        break;
      }
    }
  }

  POSTCONDITION(d_molpos.size() > idx, "not enough lines");
  dp_inStream->seekg(d_molpos[idx]);
  d_next = idx;
  return true;
}

bool MultithreadedSmilesMolSupplier::extractNextRecord(std::string &record,
                                                       unsigned int &lineNum,
                                                       unsigned int &index) {
  PRECONDITION(dp_inStream, "bad stream");
  if (d_next < 0) {
    d_next = 0;
  }

  if (moveTo(d_next)) {
    // bad index length
    CHECK_INVARIANT(static_cast<int>(d_molpos.size()) > d_next,
                    "bad index length");
    // ---------
    // if we get here we can just build the molecule:
    // ---------
    // set the stream to the relevant position:
    dp_inStream->clear();  // clear the EOF tag if it has been set
    dp_inStream->seekg(d_molpos[d_next]);
    d_line = d_lineNums[d_next];
    // grab the line:
    std::string inLine = getLine(dp_inStream);

    if (d_len < 0 && this->skipComments() < 0) {
      d_len = d_molpos.size();
    }

    // make sure the line number is correct:
    if (d_next < static_cast<int>(d_lineNums.size())) {
      d_line = d_lineNums[d_next];
    }

    ++d_next;
    // if we just hit the last one, simulate EOF:
    if (d_len > 0 && d_next == d_len) {
      df_end = true;
    }

    record = inLine;
    lineNum = d_line;
    index = d_currentRecordId;
    ++d_currentRecordId;
    return true;
  }

  return false;
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
