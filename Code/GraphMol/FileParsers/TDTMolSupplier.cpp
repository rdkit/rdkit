// $Id$
//
//  Copyright (C) 2005-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/RDLog.h>
#include "MolSupplier.h"
#include "FileParsers.h"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/LocaleSwitcher.h>


#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>


#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace RDKit {
  namespace TDTParseUtils {
    typedef boost::tokenizer<boost::escaped_list_separator<char> > CommaTokenizer;

    /*
     * if inStream is valid, we'll allow the numbers to be broken across multiple
     * lines.
     *
     * This will throw a boost::bad_lexical_cast exception if it hits a bogus number
     *
     */
    template <typename T>
    void ParseNumberList(std::string inLine,
                         std::vector<T> &res,                    
                         std::istream *inStream=0){
      bool foundEnd=false;
      while(!foundEnd) {
        CommaTokenizer commaTok(inLine);
        for(CommaTokenizer::const_iterator commaTokIt=commaTok.begin();
            commaTokIt!=commaTok.end();
            commaTokIt++){
          std::string number=*commaTokIt;
          bool atEnd= number.find(";>")!=std::string::npos;
          boost::trim_if(number,boost::is_any_of(" \r\n\t;>"));
          if(number!="" && !atEnd ){
            res.push_back(boost::lexical_cast<T>(number));
          } else if(atEnd){
            // that's it, we're done:
            foundEnd=true;
            break;
          }
        }
        if(foundEnd || !inStream || inStream->eof()){
          break;
        } else {
          std::getline(*inStream,inLine);
        }
      }
      if(!foundEnd){
        throw FileParseException("no end tag found for numeric list");
      }
    }
    

  } // end of namespace TDTParseUtils

  TDTMolSupplier::TDTMolSupplier(){
    init();
  }

  TDTMolSupplier::TDTMolSupplier(const std::string &fileName,
                                 const std::string &nameRecord,
                                 int confId2D,
                                 int confId3D,
                                 bool sanitize){
    init();
    d_confId2D=confId2D;
    d_confId3D=confId3D;
    d_nameProp=nameRecord;
    // FIX: this binary moe of opening file is here because of a bug in VC++ 6.0
    // the function "tellg" does not work correctly if we do not open it this way
    // Need to check if this has been fixed in VC++ 7.0
    std::istream *tmpStream=0;
    tmpStream = static_cast<std::istream *>(new std::ifstream(fileName.c_str(), std::ios_base::binary));
    if (!tmpStream || (!(*tmpStream)) || (tmpStream->bad()) ) {
      std::ostringstream errout;
      errout << "Bad input file " << fileName;
      throw BadFileException(errout.str());
    }

    dp_inStream = tmpStream;
    df_owner = true;
    this->advanceToNextRecord();
    d_molpos.push_back(dp_inStream->tellg());
    df_sanitize = sanitize;
    this->checkForEnd();
  }


  TDTMolSupplier::TDTMolSupplier(std::istream *inStream, bool takeOwnership,
                                 const std::string &nameRecord,
                                 int confId2D,
                                 int confId3D,
                                 bool sanitize){
    CHECK_INVARIANT(inStream,"bad instream");
    CHECK_INVARIANT(!(inStream->eof()),"early EOF");
    init();
    dp_inStream = inStream;
    df_owner = takeOwnership;
    d_confId2D=confId2D;
    d_confId3D=confId3D;
    d_nameProp=nameRecord;
    this->advanceToNextRecord();
    d_molpos.push_back(dp_inStream->tellg());
    df_sanitize = sanitize;
    this->checkForEnd();
  }


  void TDTMolSupplier::init(){
    dp_inStream=0;
    df_owner = false;
    df_end = false;
    d_len = -1;
    d_last = 0;
    d_line = 0;
  }
  TDTMolSupplier::~TDTMolSupplier() {
    if (df_owner && dp_inStream) {
      delete dp_inStream;
    }
  }


  void TDTMolSupplier::setData(const std::string &text,
                               const std::string &nameRecord,
                               int confId2D,
                               int confId3D,
                               bool sanitize){
    if(dp_inStream && df_owner) delete dp_inStream;
    init();
    d_confId2D=confId2D;
    d_confId3D=confId3D;
    d_nameProp=nameRecord;
    std::istream *tmpStream=0;
    tmpStream = static_cast<std::istream *>(new std::istringstream(text, std::ios_base::binary));
    dp_inStream = tmpStream;
    df_owner=true;
    this->advanceToNextRecord();
    d_molpos.push_back(dp_inStream->tellg());
    df_sanitize = sanitize;
    this->checkForEnd();
    POSTCONDITION(dp_inStream,"bad instream");
  }



  
  bool TDTMolSupplier::advanceToNextRecord(){
    PRECONDITION(dp_inStream,"no stream");
    std::streampos pos;
    bool res=false;
    while(1){
      if(dp_inStream->eof()) return false;
      pos = dp_inStream->tellg();
      std::string inL;
      std::getline(*dp_inStream,inL);
      if(inL.find("$SMI<")==0){
        res=true;
        break;
      }
    }
    dp_inStream->clear();
    dp_inStream->seekg(pos);
    return res;
  }
  
  void TDTMolSupplier::checkForEnd() {
    PRECONDITION(dp_inStream,"no stream");
    if (dp_inStream->eof()) {
      df_end = true;
      // the -1 here is because by the time we get here we've already pushed on the
      // position of the next line:
      d_len = d_molpos.size()-1;
      return;
    }

    // we are not at the end of file, but check for blank lines:
    std::string tempStr;
    std::getline(*dp_inStream,tempStr);
    boost::trim_left_if(tempStr,boost::is_any_of(" \t\r\n"));
    if(tempStr.length()==0){
      df_end = true;
      // the -1 here is because by the time we get here we've already pushed on the
      // position of the next line:
      d_len = d_molpos.size()-1;
    }
    return;
  }

  void TDTMolSupplier::reset() {
    PRECONDITION(dp_inStream,"no stream");
    dp_inStream->clear();
    
    dp_inStream->seekg(0, std::ios::beg);
    df_end = false;
    d_last = 0;
    d_line = 0;
  }

  ROMol *TDTMolSupplier::parseMol(std::string inLine){
    PRECONDITION(dp_inStream,"no stream");
    Utils::LocaleSwitcher ls;
    std::size_t startP=inLine.find("<");
    std::size_t endP=inLine.find_last_of(">");
    std::string smiles = inLine.substr(startP+1,endP-startP-1);
    ROMol *res = SmilesToMol(smiles,0,df_sanitize);

    if(res && res->getNumAtoms()>0){
      // -----------
      //   Process the properties:
      d_line++;
      std::getline(*dp_inStream,inLine);
      while(!dp_inStream->eof() && inLine.find("|")!=0){
        endP=inLine.find("<");
        std::string propName = inLine.substr(0,endP);
        boost::trim_if(propName,boost::is_any_of(" \t"));
        startP = endP+1;

        if(propName=="2D" && d_confId2D>=0){
          std::string rest=inLine.substr(startP,inLine.size()-startP);
          std::vector<double> coords;
          TDTParseUtils::ParseNumberList(rest,coords,dp_inStream);
          Conformer *conf=new Conformer(res->getNumAtoms());
          conf->setId(d_confId2D);
          conf->set3D(false);
          for(unsigned int atIdx=0;atIdx<res->getNumAtoms();atIdx++){
            if(2*atIdx+1 < coords.size()){
              conf->setAtomPos(atIdx,RDGeom::Point3D(coords[2*atIdx],coords[2*atIdx+1],0.0));
            } else {
              // we're going to let this slide... but maybe we should do something else?
            }
          }
          res->addConformer(conf,false);
        } else if(propName=="3D" && d_confId3D>=0){
          std::string rest=inLine.substr(startP,inLine.size()-startP);
          std::vector<double> coords;
          TDTParseUtils::ParseNumberList(rest,coords,dp_inStream);
          Conformer *conf=new Conformer(res->getNumAtoms());
          conf->setId(d_confId3D);
          conf->set3D(true);
          for(unsigned int atIdx=0;atIdx<res->getNumAtoms();atIdx++){
            if(3*atIdx+2 < coords.size()){
              conf->setAtomPos(atIdx,RDGeom::Point3D(coords[3*atIdx],
                                                     coords[3*atIdx+1],
                                                     coords[3*atIdx+2]));
            } else {
              // we're going to let this slide... but maybe we should do something else?
            }
          }
          res->addConformer(conf,false);
        } else {
          endP=inLine.find_last_of(">");
          if(endP==std::string::npos){
            std::ostringstream errout;
            errout << "no end tag found for property" << propName;
            throw FileParseException(errout.str());
          } else {
            std::string propVal = inLine.substr(startP,endP-startP);
            res->setProp(propName,propVal);
            if(propName==d_nameProp) res->setProp("_Name",propVal);
          }
        }
        std::getline(*dp_inStream,inLine);
      }
    }    
    
    return res;
  }
  
  ROMol *TDTMolSupplier::next() {
    PRECONDITION(dp_inStream,"no stream");
    // set the stream to the appropriate position
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

    // start by finding the $SMI element (we're assuming that this starts the block)
    std::string tempp;
    d_line++;
    std::getline(*dp_inStream,tempp);
    while(tempp.find("$SMI<")!=0 && !dp_inStream->eof()){
      d_line++;
      std::getline(*dp_inStream,tempp);
    }
    if(tempp.find("$SMI<")==0) {
      try {
        res = parseMol(tempp);
      }
      catch (MolSanitizeException &se) {
        // We couldn't sanitize a molecule we got - write out an error message and move to
        BOOST_LOG(rdErrorLog) << "ERROR: Could not sanitize molecule ending on line " << d_line << std::endl;
        BOOST_LOG(rdErrorLog) << "ERROR: " << se.message() << "\n";
        while(!(dp_inStream->eof()) && tempStr.find("|") != 0){
          d_line++;
          std::getline(*dp_inStream,tempStr);
        }
      }
    }
    d_last++;
    if (d_last >= static_cast<int>(d_molpos.size())) {
      d_molpos.push_back(dp_inStream->tellg());
    }
    this->checkForEnd();
    return res;
  }

  std::string TDTMolSupplier::getItemText(unsigned int idx){
    PRECONDITION(dp_inStream,"no stream");
    unsigned int holder=d_last;
    moveTo(idx);
    std::streampos begP=d_molpos[idx];
    bool endHolder=df_end;
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
    df_end=endHolder;
    char *buff=new char[endP-begP];
    dp_inStream->seekg(begP);
    dp_inStream->read(buff,endP-begP);
    std::string res(buff,endP-begP);
    delete [] buff;
    return res;
  }

  void TDTMolSupplier::moveTo(unsigned int idx) {
    PRECONDITION(dp_inStream,"no stream");
    CHECK_INVARIANT(idx >= 0, "");

    // dp_inStream->seekg() is called for all idx values
    // and earlier calls to next() may have put the stream into a bad state
    dp_inStream->clear();

    // move until we hit the desired idx
    if (idx < d_molpos.size() ) {
      dp_inStream->seekg(d_molpos[idx]);
      d_last = idx;
    }
    else {
      std::string tempStr;
      d_last = d_molpos.size() - 1;
      dp_inStream->seekg(d_molpos.back());
      while ((d_last < static_cast<int>(idx)) && (!dp_inStream->eof()) ) {
        d_line++;
        std::getline(*dp_inStream,tempStr);
        
        if (tempStr.find("|") == 0) {
          d_molpos.push_back(dp_inStream->tellg());
          d_last++;
        }
      }
      // if we reached end of file without reaching "idx" we have an index error
      if (dp_inStream->eof()) {
        d_len = d_molpos.size();
        std::ostringstream errout;
        errout << "ERROR: Index error (idx = " << idx  << ") : " << " we do no have enough molecule blocks";
        throw FileParseException(errout.str());
      }
    }
  }

  ROMol *TDTMolSupplier::operator[](unsigned int idx) {
    PRECONDITION(dp_inStream,"no stream");
    // get the molecule with index idx
    moveTo(idx);
    return next();
    
  }

  unsigned int TDTMolSupplier::length() {
    PRECONDITION(dp_inStream,"no stream");
    // return the number of mol blocks in the sdfile
    if (d_len > 0) {
      return d_len;
    }
    else {
      std::string tempStr;
      d_len = d_molpos.size();
      dp_inStream->seekg(d_molpos.back());
      std::string inL;
      std::getline(*dp_inStream,inL);
      while(this->advanceToNextRecord()){
        d_molpos.push_back(dp_inStream->tellg());
        d_len++;
        std::getline(*dp_inStream,inL);
      }
      // now remember to set the stream to the last postion we want to read
      dp_inStream->clear();
      dp_inStream->seekg(d_molpos[d_last]);
      return d_len;
    }
  }

  bool TDTMolSupplier::atEnd() {
    PRECONDITION(dp_inStream,"no stream");
    return df_end;
  }

  
}

          
