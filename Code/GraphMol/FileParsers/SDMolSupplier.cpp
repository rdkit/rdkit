// $Id$
//
//  Copyright (C) 2002-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <RDGeneral/FileParseException.h>
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
  

  std::string strip(const std::string &orig) {
    // FIX: this can be more efficeint
    // strip the end of line, white spaces and tabs
    std::string res = boost::trim_right_copy_if(orig,boost::is_any_of(" \t\r\n"));
    res = boost::trim_left_copy_if(res,boost::is_any_of(" \t\r\n"));
    return res;
  }
    
    
  SDMolSupplier::SDMolSupplier(const std::string &fileName, bool sanitize, bool removeHs){
    init();
    // FIX: this binary mode of opening file is here because of a bug in VC++ 6.0
    // the function "tellg" does not work correctly if we do not open it this way
    // Need to check if this has been fixed in VC++ 7.0
    std::istream *tmpStream=0;
    tmpStream = static_cast<std::istream *>(new std::ifstream(fileName.c_str(), std::ios_base::binary));
    if (!tmpStream || (!(*tmpStream)) || (tmpStream->bad()) ) {
      std::ostringstream errout;
      errout << "Bad input file " << fileName;
      throw FileParseException(errout.str());
    }

    //dp_inStream = static_cast<std::istream *>(tmpStream);
    dp_inStream = tmpStream;
    df_owner=true;
    d_molpos.push_back(dp_inStream->tellg());
    df_sanitize = sanitize;
    df_removeHs = removeHs;
    this->checkForEnd();
    POSTCONDITION(dp_inStream,"bad instream");
  }

  void SDMolSupplier::init(){
    dp_inStream=0;
    df_owner = false;
    df_end = false;
    d_len = -1;
    d_last = 0;
    d_line = 0;
  }

  SDMolSupplier::~SDMolSupplier() {
    if (df_owner) {
      delete dp_inStream;
    }
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
      stmp = strip(tempStr);
      if (stmp.length() == 0) {
        nempty++;
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
    
  void SDMolSupplier::readMolProps(ROMol *mol){
    PRECONDITION(dp_inStream,"no stream");
    PRECONDITION(mol,"no molecule");
    d_line++;
    std::string tempStr = getLine(dp_inStream);
    // FIX: report files missing the $$$$ marker
    while(!(dp_inStream->eof()) && (tempStr[0]!='$'||tempStr.substr(0,4)!="$$$$") ){
      tempStr = strip(tempStr);
      if(tempStr!=""){
        if (tempStr[0] == '>') { // data header line: start of a data item
            // ignore all other crap and seek for for a data label enclosed
            // by '<' and '>' 
            // FIX: "CTfile.pdf" (page 51) says that the a data header line does not
            // have to contain a data label (instead can have something line field 
            // id into a MACCS db). But we do not currently know what to do in this 
            // situation - so ignore such data items for now
            tempStr.erase(0,1); // remove the first ">" sign
            int sl = tempStr.find("<"); // begin datalabel
            int se = tempStr.find(">"); // end datalabel
            if ((sl == -1) || (se == -1) || (se == (sl+1)) ) {
              // we either do not have a data label or the label is emtpy
              // no data label ignore until next data item
              // i.e. until we hit a blank line
              d_line++;
              tempStr = getLine(dp_inStream);
              std::string stmp = strip(tempStr);
              while (stmp.length() != 0) {
                d_line++;
                tempStr = getLine(dp_inStream);
                if(dp_inStream->eof()) throw FileParseException("End of data field name not found");
              }
            } else {
              std::string dlabel = tempStr.substr(sl+1, se-sl-1);
              // we know the label - now read in the relevant properties
              // until we hit a blank line
              d_line++;
              tempStr = getLine(dp_inStream);
            
              std::string prop;
              std::string stmp = strip(tempStr);
              int nplines = 0; // number of lines for this property
              while (stmp.length() != 0) {
                nplines++;
                if (nplines > 1) {
                  prop += "\n";
                }
                prop += tempStr;
                d_line++;
                tempStr = getLine(dp_inStream);
                stmp = strip(tempStr);
              }
              mol->setProp(dlabel, prop);
            }
        } else {
          // at this point we should always be at a line starting with '>'
          // following a blank line. If this is not true throw an exception
          // FIX: should we be deleting the molecule (which is probably fine)
          // because we couldn't read the data ???
          throw FileParseException("Problems encountered parsing data fields");
        }
      }
      d_line++;
      tempStr = getLine(dp_inStream);
    }  
  }
  
  ROMol *SDMolSupplier::next() {
    PRECONDITION(dp_inStream,"no stream");
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

    unsigned int line=d_line;
    try {
      res = MolDataStreamToMol(dp_inStream, line, df_sanitize, df_removeHs);
      d_line=line;
      this->readMolProps(res);  
    }
    catch (FileParseException &fe) {
      if(d_line<static_cast<int>(line)) d_line=line;
      // we couldn't read a mol block or the data for the molecule. In this case
      // advance forward in the stream until we hit the next record and then rethrow
      // the exception. This should allow us to read the next molecule.
      BOOST_LOG(rdErrorLog) << "ERROR: on line " << d_line << " " << fe.message() << std::endl;
      BOOST_LOG(rdErrorLog) << "ERROR: moving to the begining of the next molecule\n";
      
      // FIX: report files missing the $$$$ marker
      while(!(dp_inStream->eof()) && (tempStr[0]!='$'||tempStr.substr(0,4)!="$$$$") ){
        d_line++;
        tempStr = getLine(dp_inStream);
      }
    }
    catch (MolSanitizeException &se) {
      if(d_line<static_cast<int>(line)) d_line=line;
      // We couldn't sanitize a molecule we got - write out an error message and move to
      // the beginning of the next molecule
      BOOST_LOG(rdErrorLog) << "ERROR: Could not sanitize molecule ending on line " << d_line << std::endl;
      BOOST_LOG(rdErrorLog) << "ERROR: " << se.message() << "\n";
      
      while(!(dp_inStream->eof()) && (tempStr[0]!='$'||tempStr.substr(0,4)!="$$$$") ){
        d_line++;
        tempStr = getLine(dp_inStream);
      }
    } catch (...) {
      if(d_line<static_cast<int>(line)) d_line=line;
      
      BOOST_LOG(rdErrorLog) << "Unexpected error hit on line " << d_line << std::endl;
      BOOST_LOG(rdErrorLog) << "ERROR: moving to the begining of the next molecule\n";
    }
    d_last++;
    unsigned int posHold=dp_inStream->tellg();
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
    unsigned int begP=d_molpos[idx];
    unsigned int endP;
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

    // move until we hit the desired idx
    if (idx < d_molpos.size() ) {
      dp_inStream->clear();
      dp_inStream->seekg(d_molpos[idx]);
      d_last = idx;
    } else {
      std::string tempStr;
      d_last = d_molpos.size() - 1;
      dp_inStream->seekg(d_molpos.back());
      while ((d_last < static_cast<int>(idx)) && (!dp_inStream->eof()) ) {
        d_line++;
        tempStr = getLine(dp_inStream);
        
        if (tempStr[0]=='$' && tempStr.substr(0,4)=="$$$$") {
          unsigned int posHold=dp_inStream->tellg();
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
    if (d_len > 0) {
      return d_len;
    }
    else {
      std::string tempStr;
      d_len = d_molpos.size();
      dp_inStream->seekg(d_molpos.back());
      while (!dp_inStream->eof()) {
        tempStr = getLine(dp_inStream);
        
        if (tempStr[0]=='$' && tempStr[1]=='$' && tempStr[2]=='$' && tempStr[3]=='$'){
          unsigned int posHold=dp_inStream->tellg();
          // don't worry about the last molecule:
          this->checkForEnd();
          if (!this->df_end){
            d_molpos.push_back(posHold);
            d_len++;
            dp_inStream->seekg(posHold);

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
 




