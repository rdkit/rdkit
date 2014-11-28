// $Id$
//
//  Copyright (C) 2002-2011 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/RDLog.h>
#include "MolSupplier.h"
#include "FileParsers.h"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>

namespace RDKit {
  SmilesMolSupplier::SmilesMolSupplier(){
    init();
  }


  SmilesMolSupplier::SmilesMolSupplier(const std::string &fileName, 
                                       const std::string &delimiter,
                                       int smilesColumn,
                                       int nameColumn, 
                                       bool titleLine,
                                       bool sanitize) {
    init();

    // FIX: this binary mode of opening file is here because of a bug in VC++ 6.0
    // the function "tellg" does not work correctly if we do not open it this way
    // Need to check if this has been fixed in VC++ 7.0
    std::ifstream *tmpStream = new std::ifstream(fileName.c_str(),
                                                 std::ios_base::binary);
    
    if (!tmpStream || (!(*tmpStream)) || (tmpStream->bad()) ) {
      std::ostringstream errout;
      errout << "Bad input file " << fileName;
      throw BadFileException(errout.str());
    }
    dp_inStream = static_cast<std::istream *>(tmpStream);
    CHECK_INVARIANT(dp_inStream,"bad instream");
    CHECK_INVARIANT(!(dp_inStream->eof()),"early EOF");

    d_delim = delimiter;
    df_sanitize = sanitize;
    df_title = titleLine;
    d_smi = smilesColumn;
    d_name = nameColumn;
    df_end=false;

    //if(d_title) processTitleLine();
    this->checkForEnd();
    POSTCONDITION(dp_inStream,"bad instream");
  }


  SmilesMolSupplier::SmilesMolSupplier(std::istream *inStream, bool takeOwnership,
                                       const std::string &delimiter,
                                       int smilesColumn,
                                       int nameColumn, 
                                       bool titleLine,
                                       bool sanitize) {
    CHECK_INVARIANT(inStream,"bad instream");
    CHECK_INVARIANT(!(inStream->eof()),"early EOF");

    init();
    dp_inStream = inStream;
    df_owner=takeOwnership;
    d_delim = delimiter;
    df_sanitize = sanitize;
    df_title = titleLine;
    d_smi = smilesColumn;
    d_name = nameColumn;
    df_end=false;
    this->checkForEnd();
    POSTCONDITION(dp_inStream,"bad instream");
  }



  SmilesMolSupplier::~SmilesMolSupplier() {
    if (df_owner && dp_inStream) {
      delete dp_inStream;
    }
  }

  void SmilesMolSupplier::init(){
    dp_inStream=0;
    df_owner = true;
    df_end = false;

    d_len = -1;
    d_next = -1;
    d_line = -1;
    d_molpos.clear();
    d_lineNums.clear();
  }

  void SmilesMolSupplier::setData(const std::string &text,
                                  const std::string &delimiter,
                                  int smilesColumn,
                                  int nameColumn, 
                                  bool titleLine,
                                  bool sanitize) {
    if(dp_inStream && df_owner) delete dp_inStream;
    init();

    dp_inStream=new std::stringstream(text);

    d_delim = delimiter;
    df_sanitize = sanitize;
    df_title = titleLine;
    d_smi = smilesColumn;
    d_name = nameColumn;
    df_end=false;

    this->checkForEnd();
    POSTCONDITION(dp_inStream,"bad instream");
  }


  // ensures that there is a line available to be read
  // from the file:
  void SmilesMolSupplier::checkForEnd() {
    PRECONDITION(dp_inStream,"no stream");
    int pos=this->skipComments();
    if(pos != -1){
      d_line = -1;
      dp_inStream->seekg(0);
      df_end=false;
    }
  }

  void SmilesMolSupplier::reset() {
    PRECONDITION(dp_inStream,"no stream");
    dp_inStream->clear();

    df_end=0;
    if (d_molpos.size() > 0) {
      dp_inStream->seekg(d_molpos.front());
      d_next = 0;
      d_line = 0;
    }
    else {
      dp_inStream->seekg(0);
      d_next = -1;
      d_line = -1;
    }
  }

  ROMol *SmilesMolSupplier::processLine(std::string inLine) {
    ROMol *res = NULL;

    try{
      // -----------
      // tokenize the input line:
      // -----------
      boost::char_separator<char> sep(d_delim.c_str(),"",boost::keep_empty_tokens);
      tokenizer tokens(inLine,sep);
      STR_VECT recs;
      for(tokenizer::iterator tokIter=tokens.begin();
          tokIter!=tokens.end();++tokIter){
        std::string rec = strip(*tokIter);
        recs.push_back(rec);
      }
      if(recs.size()<=static_cast<unsigned int>(d_smi)){
        std::ostringstream errout;
        errout << "ERROR: line #" << d_line  << "does not contain enough tokens\n";
        throw FileParseException(errout.str());
      }

      // -----------
      // get the smiles and create a molecule
      // -----------
      res = SmilesToMol(recs[d_smi], 0, df_sanitize);
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
        tstr << d_line;
        std::string mname = tstr.str();
        res->setProp("_Name", mname);
      }
      else {
        if(d_name>=static_cast<int>(recs.size())){
          BOOST_LOG(rdWarningLog)<<"WARNING: no name column found on line "<<d_line<<std::endl;
        } else {
          res->setProp("_Name", recs[d_name]);
        }
      }

      // -----------
      // read in the properties 
      // -----------
      unsigned int iprop = 0;
      for (unsigned int col = 0; col < recs.size(); col++) {
        if(static_cast<int>(col)==d_smi || static_cast<int>(col)==d_name) continue;
        std::string pname, pval;
        if (d_props.size() > col) {
          pname = d_props[col];
        }
        else {
          pname = "Column_";
          std::stringstream ss;
          ss << col;
          pname += ss.str();
        }

        pval = recs[col];
        res->setProp(pname, pval);
        iprop++;
      }
            
    }
    catch(const SmilesParseException &pe) {
      // Couldn't parse the passed in smiles
      // Simply print out a message
      BOOST_LOG(rdErrorLog) << "ERROR: Smiles parse error on line " << d_line << "\n";
      BOOST_LOG(rdErrorLog) << "ERROR: " << pe.message() << "\n";
      res = NULL;
    }
    catch(const MolSanitizeException &se) {
      // We couldn't sanitize the molecule
      //  write out an error message
      BOOST_LOG(rdErrorLog) << "ERROR: Could not sanitize molecule on line " << d_line << std::endl;
      BOOST_LOG(rdErrorLog) << "ERROR: " << se.message() << "\n";
      res = NULL;
    }
    catch(...) {
      //  write out an error message
      BOOST_LOG(rdErrorLog) << "ERROR: Could not process molecule on line " << d_line << std::endl;
      res = NULL;
    }
    
    return res;
  }
  

  // --------------------------------------------------
  //
  //  Returns the next available line in the input stream.
  //
  //  Side-effects:  
  //    - If EOF is hit without reading anything, the df_end
  //      flag will be set.
  //    - If a real line is read, our d_line counter is
  //      incremented
  //
  // --------------------------------------------------
  std::string  SmilesMolSupplier::nextLine(){
    PRECONDITION(dp_inStream,"bad stream");
    if(df_end) return "";
    std::string tempStr = getLine(dp_inStream);

    if( tempStr == "" ) {
      // got an empty string, check to see if we hit EOF:
      if(dp_inStream->eof()) {
        // yes, set our flag:
        df_end = true;
      }
    } else if(dp_inStream->eof()) {
      // we got some data before hitting EOF. So clear the
      // flag on inStream
      dp_inStream->clear();
    }
    d_line++;
    return tempStr;
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
  long int SmilesMolSupplier::skipComments(){
    PRECONDITION(dp_inStream,"bad stream");
    if(this->atEnd()) return -1;

    std::streampos prev = dp_inStream->tellg();
    std::string tempStr = this->nextLine();
    if(!df_end){
      // if we didn't immediately hit EOF, loop until we get a valid line:
      while((tempStr[0] == '#') || (strip(tempStr).size() == 0)) {
        prev = dp_inStream->tellg();
        tempStr = this->nextLine();
        if(this->atEnd()) break;
      }
    }
    // if we hit EOF without getting a proper line, return -1:
    if( tempStr.empty() ||
        (tempStr[0] == '#') ||
        (strip(tempStr).size() == 0) ) {
      return -1;
    }
    return static_cast<long int>(prev);
  }

  // --------------------------------------------------
  //
  //  Reads and processes the title line
  //
  void SmilesMolSupplier::processTitleLine() {
    PRECONDITION(dp_inStream,"bad stream");
    int pos = this->skipComments();
    if(pos>=0){
      dp_inStream->seekg(pos);

      std::string tempStr=getLine(dp_inStream);
      boost::char_separator<char> sep(d_delim.c_str(),"",boost::keep_empty_tokens);
      tokenizer tokens(tempStr,sep);
      for(tokenizer::iterator tokIter=tokens.begin();
          tokIter!=tokens.end();++tokIter){
        std::string pname = strip(*tokIter);
        d_props.push_back(pname);
      }
    }
  }
  
  
  std::string SmilesMolSupplier::getItemText(unsigned int idx){
    PRECONDITION(dp_inStream,"no stream");
    unsigned int holder=d_next;
    bool endHolder=df_end;
    // this throws the relevant exception if we go too far:
    moveTo(idx);
    std::string res=getLine(dp_inStream);
    d_next=holder;
    df_end=endHolder;
    return res;
  }

  // --------------------------------------------------
  //
  //  Moves to the position of a particular entry in the
  //  stream.
  //
  //  If insufficient entries are present, a FileParseException
  //    will be thrown
  //
  void SmilesMolSupplier::moveTo(unsigned int idx) {
    PRECONDITION(dp_inStream,"bad instream");
    // get the easy situations (boundary conditions) out of the
    // way first:
    if( d_len>-1 && idx>=static_cast<unsigned int>(d_len)) {
      df_end=true;
      std::ostringstream errout;
      errout << "ERROR: Index error (idx = " << idx  << "): " << "ran out of lines\n";
      throw FileParseException(errout.str());
    }

    // dp_inStream->seekg() is called for all idx values
    // and earlier calls to next() may have put the stream into a bad state
    dp_inStream->clear();

    // -----------
    // Case 1: we have already read the particular entry:
    //
    // Set the stream position and return
    // -----------
    if(!d_molpos.empty() && d_molpos.size()>idx){
      dp_inStream->clear(); // clear the EOF tag if it has been set
      df_end = false;
      dp_inStream->seekg(d_molpos[idx]);
      d_next = idx;
      d_line = d_lineNums[idx];
      return;
    }

    // -----------
    // Case 2: we haven't read the entry, so move forward until
    //   we've gone far enough.
    // -----------
    if(d_molpos.empty()){
      // if we are just starting out, process the title line
      dp_inStream->seekg(0);
      if(df_title) this->processTitleLine();
    } else {
      // move to the last position we've seen:
      dp_inStream->seekg(d_molpos.back());
      // read that line:
      std::string tmp=getLine(dp_inStream);
    }

    // the stream pointer is now at the last thing we read in
    while(d_molpos.size()<=idx){
      int nextP = this->skipComments();
      if(nextP<0){
        std::ostringstream errout;
        errout << "ERROR: Index error (idx = " << idx  << "): " << "ran out of lines\n";
        throw FileParseException(errout.str());
      } else {
        d_molpos.push_back(nextP);
        d_lineNums.push_back(d_line);
        if(d_molpos.size()==idx+1 && df_end) {
          // boundary condition: we could read the point we were looking for
          // but not the next one.
          // indicate that we've reached EOF:
          dp_inStream->clear();
          dp_inStream->seekg(0,std::ios_base::end);
          d_len=d_molpos.size();
          break;
        }
      }
    }

    POSTCONDITION(d_molpos.size()>idx,"not enough lines");
    dp_inStream->seekg(d_molpos[idx]);
    d_next = idx;
    return;
  }

  // ----------------------------------------------------------------------
  //
  //  Grabs and returns the next molecule from the input stream.
  //  After processing the line, the file is advanced to the next
  //  position in the file (skipping blank and comment lines).
  //
  //  Throws a FileParseException if EOF has already been hit.
  //
  ROMol *SmilesMolSupplier::next() {
    PRECONDITION(dp_inStream,"no stream");
    ROMol *res=NULL;

    if(d_next<0){
      d_next=0;
    }

    // This throws an exception if it fails:
    moveTo(d_next);
    CHECK_INVARIANT(static_cast<int>(d_molpos.size())>d_next,"bad index length");

    // ---------
    // if we get here we can just build the molecule:
    // ---------
    // set the stream to the relevant position:
    dp_inStream->clear(); // clear the EOF tag if it has been set
    dp_inStream->seekg(d_molpos[d_next]);
    d_line = d_lineNums[d_next];
    // grab the line:
    std::string inLine=getLine(dp_inStream);
    // and process it:
    res=this->processLine(inLine);

    // if we don't already know the length of the supplier,
    // check if we can read another line:
    if(d_len<0&&this->skipComments()<0){
      d_len=d_molpos.size();
    }

    // make sure the line number is correct:
    if(d_next < static_cast<int>(d_lineNums.size()) ){
      d_line = d_lineNums[d_next];
    }

    ++d_next;
    // if we just hit the last one, simulate EOF:
    if(d_len>0 && d_next==d_len ) {
      df_end=true;
    }
    return res;
  }


  // ----------------------------------------------------------------------
  //
  //  Grabs and returns a particular molecule from the input stream.
  //
  //  Raises a FileParseException on failure.
  //
  ROMol *SmilesMolSupplier::operator[](unsigned int idx) {
    PRECONDITION(dp_inStream,"no stream");

    // ---------
    // move to the appropriate location in the file:
    // ---------
    moveTo(idx);

    // ---------
    // and then pull the molecule:
    // ---------
    ROMol *res = next();

    return res;
  }

  // ----------------------------------------------------------------------
  //
  //  Returns the number of entries in the input stream
  //
  unsigned int SmilesMolSupplier::length() {
    PRECONDITION(dp_inStream,"no stream");
    // return the number of molecule lines in the file
    if (d_len > 0) {
      return d_len;
    }
    else {
      std::streampos oPos = dp_inStream->tellg();
      if(d_molpos.size()){
        // we've already read some molecules, go to the last
        // one and read it in to initialize our location:
        dp_inStream->seekg(d_molpos.back());
        // skip that line and then continue:
        this->skipComments();
      } else {
        // process the title line if need be:
        if(df_title) this->processTitleLine();
      }
      int pos = this->skipComments();
      while(pos>=0){
        d_molpos.push_back(pos);
        d_lineNums.push_back(d_line);
        pos = this->skipComments();
      }
      // now remember to set the stream to its original position:
      dp_inStream->seekg(oPos);
      d_len = d_molpos.size();
      return d_len;
    }
  }

  bool SmilesMolSupplier::atEnd() {
    return df_end;
  }
}

          
