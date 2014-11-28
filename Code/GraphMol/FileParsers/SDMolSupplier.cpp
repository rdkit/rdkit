// $Id$
//
//  Copyright (C) 2002-2012 Greg Landrum and Rational Discovery LLC
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
#include <iostream>
#include <sstream>
#include <string>

namespace RDKit {
  

  SDMolSupplier::SDMolSupplier(const std::string &fileName, bool sanitize, bool removeHs,
                               bool strictParsing){
    init();
    // FIX: this binary mode of opening file is here because of a bug in VC++ 6.0
    // the function "tellg" does not work correctly if we do not open it this way
    //   Jan 2009: Confirmed that this is still the case in visual studio 2008
    std::istream *tmpStream=0;
    tmpStream = static_cast<std::istream *>(new std::ifstream(fileName.c_str(), std::ios_base::binary));
    if (!tmpStream || (!(*tmpStream)) || (tmpStream->bad()) ) {
      std::ostringstream errout;
      errout << "Bad input file " << fileName;
      throw BadFileException(errout.str());
    }

    //dp_inStream = static_cast<std::istream *>(tmpStream);
    dp_inStream = tmpStream;
    df_owner=true;
    d_molpos.push_back(dp_inStream->tellg());
    df_sanitize = sanitize;
    df_removeHs = removeHs;
    df_strictParsing = strictParsing;
    this->checkForEnd();
    if(df_end){
      // checkForEnd() sets d_len if we're at EOF. undo that (was GitHub issue 19):
      d_len=0;
    }
    POSTCONDITION(dp_inStream,"bad instream");
  }


  SDMolSupplier::SDMolSupplier(std::istream *inStream, bool takeOwnership,
                               bool sanitize, bool removeHs,bool strictParsing){
    PRECONDITION(inStream,"bad stream");
    init();
    dp_inStream = inStream;
    df_owner=takeOwnership;
    d_molpos.push_back(dp_inStream->tellg());
    df_sanitize = sanitize;
    df_removeHs = removeHs;
    df_strictParsing = strictParsing;
    this->checkForEnd();
    if(df_end){
      // checkForEnd() sets d_len if we're at EOF. undo that (was GitHub issue 19):
      d_len=0;
    }
    POSTCONDITION(dp_inStream,"bad instream");
  }

  void SDMolSupplier::init(){
    ForwardSDMolSupplier::init();
    d_len = -1;
    d_last = 0;
  }

  void SDMolSupplier::setData(const std::string &text,
                              bool sanitize, bool removeHs){
    if(dp_inStream && df_owner) delete dp_inStream;
    init();
    std::istream *tmpStream=0;
    tmpStream = static_cast<std::istream *>(new std::istringstream(text, std::ios_base::binary));
    dp_inStream = tmpStream;
    df_owner=true;
    d_molpos.push_back(dp_inStream->tellg());
    df_sanitize = sanitize;
    df_removeHs=removeHs;
    this->checkForEnd();
    if(df_end){
      // checkForEnd() sets d_len if we're at EOF. undo that (was GitHub issue 19):
      d_len=0;
    }
    POSTCONDITION(dp_inStream,"bad instream");
  }

  
  void SDMolSupplier::checkForEnd() {
    PRECONDITION(dp_inStream,"no stream");
    // we will call it end of file if we have more than 4 contiguous empty lines
    // or we reach end of file in the meantime
    if (dp_inStream->eof()) {
      df_end = true;
      d_len = d_molpos.size();
      return;
    }
    // we are not at the end of file, check for blank lines
    unsigned int nempty = 0;
    std::string tempStr, stmp;
    for (unsigned int i = 0; i < 4 ; i++) {
      tempStr = getLine(dp_inStream);
      if (dp_inStream->eof()) {
        df_end = true;
        d_len = d_molpos.size();
        return;
      }
      if(tempStr.find_first_not_of(" \t\r\n")==std::string::npos){
        ++nempty;
      }
    }
    if (nempty == 4) {
      df_end = true;
      d_len = d_molpos.size();
    }
  }

  void SDMolSupplier::reset() {
    PRECONDITION(dp_inStream,"no stream");
    dp_inStream->clear();
    dp_inStream->seekg(0, std::ios::beg);
    df_end = false;
    d_last = 0;
    d_line = 0;
  }
    
  ROMol *SDMolSupplier::next() {
    PRECONDITION(dp_inStream,"no stream");
    if(df_end && d_last>=d_len){
      throw FileParseException("EOF hit.");
    }

    // set the stream to the current position
    dp_inStream->seekg(d_molpos[d_last]);

    std::string tempStr;
    ROMol *res = NULL;
    // finally if we reached the end of the file set end to be true
    if (dp_inStream->eof()) {
      //FIX: we should probably be throwing an exception here
      df_end = true;
      d_len = d_molpos.size();
      return res;
    }

    res = _next();

    ++d_last;
    std::streampos posHold=dp_inStream->tellg();
    this->checkForEnd();
    if (!this->df_end && d_last >= static_cast<int>(d_molpos.size())) {
      d_molpos.push_back(posHold);
    }

    return res;
  }

  std::string SDMolSupplier::getItemText(unsigned int idx){
    PRECONDITION(dp_inStream,"no stream");
    unsigned int holder=d_last;
    moveTo(idx);
    std::streampos begP=d_molpos[idx];
    std::streampos endP;
    try {
      moveTo(idx+1);
      endP=d_molpos[idx+1];
    } catch (FileParseException &) {
      dp_inStream->clear();
      dp_inStream->seekg(0,std::ios_base::end);
      endP=dp_inStream->tellg();
    }
    d_last=holder;
    char *buff=new char[endP-begP];
    dp_inStream->seekg(begP);
    dp_inStream->read(buff,endP-begP);
    std::string res(buff,endP-begP);
    delete [] buff;
    return res;
  }

  void SDMolSupplier::moveTo(unsigned int idx) {
    PRECONDITION(dp_inStream,"no stream");

    // dp_inStream->seekg() is called for all idx values
    // and earlier calls to next() may have put the stream into a bad state
    dp_inStream->clear();

    // move until we hit the desired idx
    if (idx < d_molpos.size() ) {
      dp_inStream->seekg(d_molpos[idx]);
      d_last = idx;
    } else {
      std::string tempStr;
      dp_inStream->seekg(d_molpos.back());
      d_last = d_molpos.size() - 1;
      while ((d_last < static_cast<int>(idx)) && (!dp_inStream->eof()) ) {
        d_line++;
        tempStr = getLine(dp_inStream);
        
        if (tempStr[0]=='$' && tempStr.substr(0,4)=="$$$$") {
          std::streampos posHold=dp_inStream->tellg();
          this->checkForEnd();
          if (!this->df_end){
            d_molpos.push_back(posHold);
            d_last++;
          }
        }
      }
      // if we reached end of file without reaching "idx" we have an index error
      if (dp_inStream->eof()) {
        d_len = d_molpos.size();
        std::ostringstream errout;
        errout << "ERROR: Index error (idx = " << idx  << ") : " << " we do no have enough mol blocks";
        throw FileParseException(errout.str());
      }
    }
  }

  ROMol *SDMolSupplier::operator[](unsigned int idx) {
    PRECONDITION(dp_inStream,"no stream");
    // get the molecule with index idx
    moveTo(idx);
    return next();
    
  }

  unsigned int SDMolSupplier::length() {
    PRECONDITION(dp_inStream,"no stream");
    // return the number of mol blocks in the sdfile
    if (d_len > 0 || (df_end && d_len==0) ) {
      return d_len;
    } else {
      std::string tempStr;
      d_len = d_molpos.size();
      dp_inStream->seekg(d_molpos.back());
      while (!dp_inStream->eof()) {
        std::getline(*dp_inStream,tempStr);
        if (tempStr.length()>=4 && tempStr[0]=='$' && tempStr[1]=='$' && tempStr[2]=='$' && tempStr[3]=='$'){
          std::streampos posHold=dp_inStream->tellg();
          // don't worry about the last molecule:
          this->checkForEnd();
          if (!this->df_end){
            d_molpos.push_back(posHold);
            ++d_len;
          }
        }
      }
      // now remember to set the stream to the last postion we want to read
      dp_inStream->clear();
      dp_inStream->seekg(d_molpos[d_last]);
      return d_len;
    }
  }

  bool SDMolSupplier::atEnd() {
    PRECONDITION(dp_inStream,"no stream");
    return df_end;
  }

  void SDMolSupplier::setStreamIndices(const std::vector<std::streampos> &locs){
    d_molpos.clear();
    d_molpos.resize(locs.size());
    std::copy(locs.begin(),locs.end(),d_molpos.begin());
    this->reset();
    d_len = d_molpos.size();
  }
}
 




