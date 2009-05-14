// $Id$
//
//  Copyright (C) 2009 Greg Landrum
//
//   @@ All Rights Reserved  @@
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
                                             bool sanitize, bool removeHs){
    PRECONDITION(inStream,"bad stream");
    init();
    dp_inStream = inStream;
    df_owner=takeOwnership;
    df_sanitize = sanitize;
    df_removeHs = removeHs;
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
              while (stmp.length() != 0) {
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

    std::string tempStr;
    ROMol *res = NULL;
    // finally if we reached the end of the file set end to be true
    if (dp_inStream->eof()) {
      //FIX: we should probably be throwing an exception here
      df_end = true;
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

    return res;
  }

  bool ForwardSDMolSupplier::atEnd() {
    PRECONDITION(dp_inStream,"no stream");
    return df_end;
  }
}
 




