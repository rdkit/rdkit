// $Id$
//
//  Copyright (C) 2004-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "FeatureParser.h"
#include "MolChemicalFeatureDef.h"
#include <RDGeneral/StreamOps.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <boost/shared_ptr.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

#include <fstream>
#include <sstream>
#include <map>

namespace RDKit {
  namespace Local {
    typedef boost::tokenizer<boost::escaped_list_separator<char> > CommaTokenizer;

    void getNextLine(std::istream &inStream,std::string &line,unsigned int &lineNo){
      if(inStream.eof()) return;
      line = "";
      bool continuationLine=false;
      while(!inStream.eof()){
        std::string tmpLine;
        std::getline(inStream,tmpLine);
        lineNo++;
        //std::cerr << ">> " << lineNo << " " << tmpLine << std::endl;
        if(tmpLine=="") continue;
        if(tmpLine[0]!='#'){
          // strip space at the end to check for a continuation line:
          std::string stripLine=boost::trim_right_copy_if(tmpLine,
                                                          boost::is_any_of(" \t\r\n"));
          if(stripLine=="") continue;
          if(stripLine[stripLine.size()-1]!='\\'){
            if(continuationLine){
              // if it's a continuation line, strip any whitespace:
              boost::trim_if(tmpLine,boost::is_any_of(" \t\r\n"));
            }
            line += tmpLine;
            return;
          } else {
            continuationLine=true;
            boost::trim_if(tmpLine,boost::is_any_of(" \t\r\n"));
            line += tmpLine.substr(0,tmpLine.size()-1);
          }
        }
      } 
    }
    
    // ------------------------------------------------------
    bool expandAndTestSmarts(std::string &smarts,
                             const std::map<std::string,std::string> &atomTypeDefs) {

      for(std::map<std::string,std::string>::const_iterator mapIt=atomTypeDefs.begin();
          mapIt!=atomTypeDefs.end();mapIt++){
        std::string atomName=mapIt->first;
        std::string atomSma=mapIt->second;
        boost::replace_all(smarts,atomName,atomSma);
      }

      RWMol *mol=0;
      try{
        mol=SmartsToMol(smarts);
      } catch (SmilesParseException &) {
        return false;
      }
      if(mol){
        delete mol;
      } else {
        return false;
      }
      return true;
    }

    // ------------------------------------------------------
    void parseAtomType(const std::string &inLine,
                       std::map<std::string,std::string> &atomTypeDefs,
                       const unsigned int &lineNo) {
      boost::char_separator<char> sep(" \t");
      boost::tokenizer<boost::char_separator<char> > tok(inLine,sep);
      boost::tokenizer<boost::char_separator<char> >::iterator tokIt=tok.begin();
      if(tokIt==tok.end()){
        throw FeatureFileParseException(lineNo,inLine,
                                        "empty input line for AtomType");         
      }
      std::string keyword=boost::to_upper_copy(*tokIt);
      if(keyword!="ATOMTYPE"){
        throw FeatureFileParseException(lineNo,inLine,
                                        "bad input line for AtomType");   
      }


      tokIt++;
      if(tokIt==tok.end()){
        throw FeatureFileParseException(lineNo,inLine,
                                        "bad AtomType line, missing label");      
      }
      std::string atomType=*tokIt;
      bool negater=false;
      if(atomType[0]=='!'){
        atomType.erase(0,1);
        negater=true;
      }
      atomType="{" + atomType + "}";
      tokIt++;
      if(tokIt==tok.end()){
        throw FeatureFileParseException(lineNo,inLine,
                                        "bad AtomType line, missing definition");         
      }
      std::string sma;
      if(atomTypeDefs.count(atomType)){
        std::string base=atomTypeDefs[atomType];
        sma="$(" + *tokIt + ")";
        if(negater){
          std::string toAdd="[!"+sma+";";
          boost::replace_first(base,"[",toAdd);
        } else {
          std::string toAdd=","+sma+"]";
          boost::replace_last(base,"]",toAdd);
        }
        sma = base;
      } else {
        sma="$(" + *tokIt + ")";
      }
      // make it a valid smarts definition for an atom:
      sma = "["+sma+"]";
      // make sure we get sensible SMARTS:
      if(!expandAndTestSmarts(sma,atomTypeDefs)){
        std::string msg="invalid SMARTS in AtomType (" + atomType + "): "+sma;
        throw FeatureFileParseException(lineNo,inLine,msg);
      }
      // now cut the brackets back off:
      sma = sma.substr(1,sma.size()-2);
      atomTypeDefs[atomType] = sma;
    }

    // ------------------------------------------------------
    MolChemicalFeatureDef *
    parseFeatureDef(std::istream &inStream,
                    const std::string &inLine,
                    unsigned int &lineNo,
                    const std::map<std::string,std::string> &atomTypeDefs){
      std::string nextLine=inLine;
      MolChemicalFeatureDef *res=0;

      // handle a blank or comment first line:
      boost::trim_if(nextLine,boost::is_any_of(" \t\r\n"));
      while(nextLine=="" || nextLine[0]=='#'){
        Local::getNextLine(inStream,nextLine,lineNo);
        // need to check for EOS before we strip:
        if(nextLine=="") {
          // we hit EOS:
          throw FeatureFileParseException(lineNo,inLine,
                                          "EOF hit parsing feature definition");
        }
        boost::trim_if(nextLine,boost::is_any_of(" \t\r\n"));
      }

      boost::char_separator<char> sep(" \t");
      boost::tokenizer<boost::char_separator<char> > tok(nextLine,sep);
      if(tok.begin()==tok.end()){
        throw FeatureFileParseException(lineNo,inLine,
                                        "bad DefineFeature line, no tokens found");
      }
      boost::tokenizer<boost::char_separator<char> >::iterator tokIt=tok.begin();
      tokIt++;
      if(tokIt==tok.end()){
        throw FeatureFileParseException(lineNo,inLine,
                                        "bad DefineFeature line, missing subtype");
      }
      std::string subType=*tokIt;
      tokIt++;
      if(tokIt==tok.end()){
        throw FeatureFileParseException(lineNo,inLine,
                                        "bad DefineFeature line, missing pattern");
      }
      std::string pattern=*tokIt;

      //---------------
      // make sure we get sensible SMARTS:
      //
      if(!expandAndTestSmarts(pattern,atomTypeDefs)){
        std::string msg="invalid SMARTS in DefineFeature for type "+subType+": "+pattern;
        throw FeatureFileParseException(lineNo,inLine,msg);
      }


      //---------------
      // read out the rest of the definition
      //
      std::vector<double> weights;
      std::string family="";
      bool foundEnd=false;
      Local::getNextLine(inStream,nextLine,lineNo);
      //std::getline(inStream,nextLine);
      while(nextLine != ""){
        boost::trim_if(nextLine,boost::is_any_of(" \t\r\n"));
        if(nextLine != "" && nextLine[0]!='#'){
          tok.assign(nextLine,sep);
          tokIt = tok.begin();
          std::string token=boost::to_upper_copy(*tokIt);
          if(token=="ENDFEATURE"){
            foundEnd=true;
            break;
          } else if(token=="FAMILY"){
            tokIt++;
            if(tokIt==tok.end()){
              std::string msg="bad Type line for feature: "+subType;
              throw FeatureFileParseException(lineNo,inLine,msg);
            }
            family = *tokIt;
          } else if(token=="WEIGHTS"){
            tokIt++;
            if(tokIt==tok.end()){
              std::string msg="bad Weights line for feature: "+subType;
              throw FeatureFileParseException(lineNo,inLine,msg);
            }
            CommaTokenizer commaTok(*tokIt);
            for(CommaTokenizer::const_iterator commaTokIt=commaTok.begin();
                commaTokIt!=commaTok.end();commaTokIt++){
              std::string number=*commaTokIt;
              try {
                weights.push_back(boost::lexical_cast<double>(number));
              } catch (boost::bad_lexical_cast &){
                std::string msg="bad weight value (" + number + ") for feature: "+subType;
                throw FeatureFileParseException(lineNo,inLine,msg);
              }
            }
          } else {
            std::string msg="bad input line for feature: "+subType;
            throw FeatureFileParseException(lineNo,inLine,msg);
          }
        }
        Local::getNextLine(inStream,nextLine,lineNo);
        //std::getline(inStream,nextLine);
      }
      if(!foundEnd){
        std::string msg="could not find EndFeature line for feature: "+subType;
        throw FeatureFileParseException(lineNo,inLine,msg);
      }
      if(family==""){
        std::string msg="did not find Family definition for feature: "+subType;
        throw FeatureFileParseException(lineNo,inLine,msg);
      }

      //---------------
      // Build the feature definition
      //
      res = new MolChemicalFeatureDef(pattern,family,subType);
      if(weights.size()){
        res->setWeights(weights);
        res->normalizeWeights();
      }
      
      
      return res;
    }
  } // end of namespace Local

  // ------------------------------------------------------
  int parseFeatureData(const std::string &defnText,
                       MolChemicalFeatureDef::CollectionType &featDefs) {
    std::stringstream ss(defnText);
    return parseFeatureData(ss,featDefs);
  }

  
  // ------------------------------------------------------
  int parseFeatureData(std::istream &inStream,
                       MolChemicalFeatureDef::CollectionType &res){
    unsigned int lineNo=0;
    std::string inLine;
    Local::getNextLine(inStream,inLine,lineNo);
    std::map<std::string,std::string> atomTypeDefs;
    while(!inStream.eof()) {
      // clean any whitespace off the line:
      boost::trim_if(inLine,boost::is_any_of(" \t\r\n"));
      if(inLine != "" && inLine[0]!='#' && inLine[0]!='\n'){
        boost::tokenizer<> tok(inLine);
        boost::tokenizer<>::iterator tokIt=tok.begin();
        std::string token=boost::to_upper_copy(*tokIt);
        if(token=="ATOMTYPE"){
          Local::parseAtomType(inLine,atomTypeDefs,lineNo);
        } else if(token=="DEFINEFEATURE"){
          MolChemicalFeatureDef *fDef=Local::parseFeatureDef(inStream,inLine,
                                                          lineNo,atomTypeDefs);
          if(fDef) res.push_back(boost::shared_ptr<MolChemicalFeatureDef>(fDef));
        } else {
          throw FeatureFileParseException(lineNo,inLine,"bad or missing keyword");
        }
      }
      //std::getline(inStream,inLine);
      Local::getNextLine(inStream,inLine,lineNo);
    }
    return 0;
  }

  // ------------------------------------------------------
  int parseFeatureFile(const std::string &fileName,
                       MolChemicalFeatureDef::CollectionType &res){
    std::ifstream inStream(fileName.c_str());
    if( !inStream || inStream.eof() ){
      return -1;
    }
    return parseFeatureData(inStream,res);
  }


}
