// $Id$
//
//  Copyright (C) 2009-2012 Greg Landrum
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
  std::string strip(const std::string &orig) {
    // FIX: this can be more efficeint
    // strip the end of line, white spaces and tabs
    std::string res = boost::trim_right_copy_if(orig,boost::is_any_of(" \t\r\n"));
    res = boost::trim_left_copy_if(res,boost::is_any_of(" \t\r\n"));
    return res;
  }
    
  ForwardSDMolSupplier::ForwardSDMolSupplier(std::istream *inStream, bool takeOwnership,
                                             bool sanitize, bool removeHs,bool strictParsing){
    PRECONDITION(inStream,"bad stream");
    init();
    dp_inStream = inStream;
    df_owner=takeOwnership;
    df_sanitize = sanitize;
    df_removeHs = removeHs;
    df_strictParsing = strictParsing;
    POSTCONDITION(dp_inStream,"bad instream");
  }

  void ForwardSDMolSupplier::init(){
    dp_inStream=0;
    df_owner = false;
    df_end = false;
    d_line = 0;
  }

  void ForwardSDMolSupplier::reset() {
    UNDER_CONSTRUCTION("reset() not supported for ForwardSDMolSuppliers();");
  }
    
  void ForwardSDMolSupplier::readMolProps(ROMol *mol){
    PRECONDITION(dp_inStream,"no stream");
    PRECONDITION(mol,"no molecule");
    d_line++;
    std::string tempStr;
    std::getline(*dp_inStream,tempStr);

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
              std::getline(*dp_inStream,tempStr);
              std::string stmp = strip(tempStr);
              while (stmp.length() != 0) {
                d_line++;
                std::getline(*dp_inStream,tempStr);
                if(dp_inStream->eof()) throw FileParseException("End of data field name not found");
              }
            } else {
              std::string dlabel = tempStr.substr(sl+1, se-sl-1);
              // we know the label - now read in the relevant properties
              // until we hit a blank line
              d_line++;
              std::getline(*dp_inStream,tempStr);

              std::string prop="";
              std::string stmp = strip(tempStr);
              int nplines = 0; // number of lines for this property
              while (stmp.length() != 0 || tempStr[0]==' ' || tempStr[0]=='\t') {
                nplines++;
                if (nplines > 1) {
                  prop += "\n";
                }
                // take off \r if it's still in the property:
                if (tempStr[tempStr.length()-1]=='\r'){
                  tempStr.erase(tempStr.length()-1);
                }
                prop += tempStr;
                d_line++;
                // erase tempStr in case the file does not end with a carrier
                // return (we will end up in an infinite loop if we don't do
                // this and we do not check for EOF in this while loop body)
                tempStr.erase();
                std::getline(*dp_inStream,tempStr);
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
      std::getline(*dp_inStream,tempStr);
    }  
  }
  
  ROMol *ForwardSDMolSupplier::next() {
    PRECONDITION(dp_inStream,"no stream");
    ROMol *res = NULL;

    if (dp_inStream->eof()) {
      //FIX: we should probably be throwing an exception here
      df_end = true;
      return res;
    }

    res=_next();
    return res;
  }

  ROMol *ForwardSDMolSupplier::_next() {
    PRECONDITION(dp_inStream,"no stream");

    std::string tempStr;
    ROMol *res = NULL;
    if (dp_inStream->eof()) {
      df_end = true;
      return res;
    }

    unsigned int line=d_line;
    try {
      res = MolDataStreamToMol(dp_inStream, line, df_sanitize, df_removeHs,df_strictParsing);
      d_line=line;
      if(res){
        this->readMolProps(res);
      } else if(!dp_inStream->eof()) {
        // FIX: report files missing the $$$$ marker
        std::getline(*dp_inStream,tempStr);
        ++d_line;
        while(!(dp_inStream->eof()) && (tempStr[0]!='$'||tempStr.substr(0,4)!="$$$$") ){
          std::getline(*dp_inStream,tempStr);
          ++d_line;
        }
      }
    }
    catch (FileParseException &fe) {
      if(d_line<static_cast<int>(line)) d_line=line;
      // we couldn't read a mol block or the data for the molecule. In this case
      // advance forward in the stream until we hit the next record and then rethrow
      // the exception. This should allow us to read the next molecule.
      BOOST_LOG(rdErrorLog) << "ERROR: " << fe.message() << std::endl;
      BOOST_LOG(rdErrorLog) << "ERROR: moving to the begining of the next molecule\n";
      
      // FIX: report files missing the $$$$ marker
      while(!(dp_inStream->eof()) && (tempStr[0]!='$'||tempStr.substr(0,4)!="$$$$") ){
        d_line++;
        std::getline(*dp_inStream,tempStr);
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
        std::getline(*dp_inStream,tempStr);
      }
    } catch (...) {
      if(d_line<static_cast<int>(line)) d_line=line;
      
      BOOST_LOG(rdErrorLog) << "Unexpected error hit on line " << d_line << std::endl;
      BOOST_LOG(rdErrorLog) << "ERROR: moving to the begining of the next molecule\n";
      while(!(dp_inStream->eof()) && (tempStr[0]!='$'||tempStr.substr(0,4)!="$$$$") ){
        d_line++;
        std::getline(*dp_inStream,tempStr);
      }
    }
    if (dp_inStream->eof()) {
      //FIX: we should probably be throwing an exception here
      df_end = true;
    }
    return res;
  }


  void ForwardSDMolSupplier::checkForEnd() {
    PRECONDITION(dp_inStream,"no stream");
    // we will call it end of file if we have more than 4 contiguous empty lines
    // or we reach end of file in the meantime
    if (dp_inStream->eof()) {
      df_end = true;
      return;
    }
    // we are not at the end of file, check for blank lines
    unsigned int nempty = 0;
    std::string tempStr;
    for (unsigned int i = 0; i < 4 ; i++) {
      tempStr = getLine(dp_inStream);
      if (dp_inStream->eof()) {
        df_end = true;
        return;
      }
      if(tempStr.find_first_not_of(" \t\r\n")==std::string::npos){
        ++nempty;
      }
    }
    if (nempty == 4) {
      df_end = true;
    }
  }


  bool ForwardSDMolSupplier::atEnd() {
    PRECONDITION(dp_inStream,"no stream");
    return df_end;
  }
}
 




