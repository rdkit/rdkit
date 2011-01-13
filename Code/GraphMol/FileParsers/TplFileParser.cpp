// $Id$
//
//  Copyright (C) 2007-2010 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "FileParsers.h"
#include "FileParserUtils.h"
#include "MolFileStereochem.h"
#include <RDGeneral/StreamOps.h>

#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <typeinfo>

namespace RDKit{
  void ParseTPLAtomLine(std::string text,unsigned int lineNum,RWMol *mol,
                        Conformer *conf){
    PRECONDITION(mol,"no molecule");
    PRECONDITION(conf,"no conformer");
    std::vector<std::string> splitLine;
    boost::split(splitLine,text,boost::is_any_of(" \t"),boost::token_compress_on);
    if(splitLine.size()<8){
      std::ostringstream errout;
      errout << "Atom line " << lineNum << " only has " << splitLine.size() << " tokens. 8 are required."<<std::endl;
      throw FileParseException(errout.str()) ;
    }
    Atom *atom=new Atom(splitLine[1]);
    unsigned int atomId;
    atomId=mol->addAtom(atom,false,true);
    
    atom->setFormalCharge(FileParserUtils::stripSpacesAndCast<int>(splitLine[2]));
    double partialChg=FileParserUtils::stripSpacesAndCast<double>(splitLine[3]);
    atom->setProp("TPLCharge",partialChg);
    double xp=FileParserUtils::stripSpacesAndCast<double>(splitLine[4]);
    double yp=FileParserUtils::stripSpacesAndCast<double>(splitLine[5]);
    double zp=FileParserUtils::stripSpacesAndCast<double>(splitLine[6]);
    // coords in TPL files are in picometers, adjust:
    xp/=100.;
    yp/=100.;
    zp/=100.;
    conf->setAtomPos(atomId,RDGeom::Point3D(xp,yp,zp));

    unsigned int nBonds=FileParserUtils::stripSpacesAndCast<unsigned int>(splitLine[7]);
    // the only remaining info we care about is stereochem, and then only if
    // the number of bonds is 4:
    if(nBonds==4 && splitLine.size()>8+nBonds){
      std::string stereoChem=splitLine[8+nBonds];
      atom->setProp("TPLStereoFlag",stereoChem);
    }
  }

  void ParseTPLBondLine(std::string text,unsigned int lineNum,RWMol *mol){
    PRECONDITION(mol,"no molecule");

    std::vector<std::string> splitLine;
    boost::split(splitLine,text,boost::is_any_of(" \t"),boost::token_compress_on);
    if(splitLine.size()<5){
      std::ostringstream errout;
      errout << "Bond line " << lineNum << " only has " << splitLine.size() << " tokens. 5 are required."<<std::endl;
      throw FileParseException(errout.str()) ;
    }

    std::string tplOrder=boost::trim_copy(splitLine[1]).substr(0,3);
    Bond::BondType bondOrder;
    if(tplOrder=="1.5"){
      bondOrder=Bond::AROMATIC;
    } else if(tplOrder=="1.0"){
      bondOrder=Bond::SINGLE;
    } else if(tplOrder=="2.0"){
      bondOrder=Bond::DOUBLE;
    } else if(tplOrder=="3.0"){
      bondOrder=Bond::TRIPLE;
    } else {
      std::ostringstream errout;
      errout << "Bond line " << lineNum << " has unknown order: " << tplOrder <<std::endl;
      throw FileParseException(errout.str()) ;
    }
    unsigned int idx1,idx2;
    idx1 = FileParserUtils::stripSpacesAndCast<unsigned int>(splitLine[2])-1;
    idx2 = FileParserUtils::stripSpacesAndCast<unsigned int>(splitLine[3])-1;
    
    unsigned int bondIdx=mol->addBond(idx1,idx2,bondOrder)-1;
    std::string stereoFlag1="";
    std::string stereoFlag2="";
    if(splitLine.size()>4){
      stereoFlag1=splitLine[4];
      if(splitLine.size()>5){
        stereoFlag2=splitLine[5];
      }
    }
    mol->getBondWithIdx(bondIdx)->setProp("TPLBondDir1",stereoFlag1);
    mol->getBondWithIdx(bondIdx)->setProp("TPLBondDir2",stereoFlag2);
  }

  Conformer *ParseConfData(std::istream *inStream,unsigned int &line,RWMol *mol,
                           unsigned int confId){
    PRECONDITION(inStream,"no stream");
    PRECONDITION(mol,"no mol");

    std::string tempStr;
    std::vector<std::string> splitLine;

    line++;
    tempStr = getLine(inStream);
    boost::split(splitLine,tempStr,boost::is_any_of(" \t"),boost::token_compress_on);
    if(splitLine[0]!="NAME"){
      std::ostringstream errout;
      errout << "Did not find NAME tag on line " << line << " while reading conformer  " << confId << std::endl;
      throw FileParseException(errout.str()) ;
    }
    std::ostringstream propName;
    propName << "Conf_" <<mol->getNumConformers()<<"_Name";
    mol->setProp(propName.str(),boost::trim_copy(tempStr.substr(4,tempStr.size()-4)));
    
    Conformer *conf=new Conformer(mol->getNumAtoms());
    for(unsigned int i=0;i<mol->getNumAtoms();++i){
      line++;
      tempStr = getLine(inStream);
      if(inStream->eof()){
        std::ostringstream errout;
        errout << "EOF hit while reading conformer  " << confId << std::endl;
        throw FileParseException(errout.str()) ;
      }
      boost::trim(tempStr);
      boost::split(splitLine,tempStr,
                   boost::is_any_of(" \t"),boost::token_compress_on);
      if(splitLine.size()<3){
        std::ostringstream errout;
        errout << "Did not find enough fields on line " << line << " while reading conformer  " << confId << std::endl;
        throw FileParseException(errout.str()) ;
      }
      double xp=FileParserUtils::stripSpacesAndCast<double>(splitLine[0]);
      double yp=FileParserUtils::stripSpacesAndCast<double>(splitLine[1]);
      double zp=FileParserUtils::stripSpacesAndCast<double>(splitLine[2]);
      // coords in TPL files are in picometers, adjust:
      xp/=100.;
      yp/=100.;
      zp/=100.;
      conf->setAtomPos(i,RDGeom::Point3D(xp,yp,zp));
    }
    return conf;
  }
  
  
  //*************************************
  //
  // Every effort has been made to adhere to the BioCad tpl definition 
  //  
  //*************************************
  RWMol *TPLDataStreamToMol(std::istream *inStream, unsigned int &line,
                            bool sanitize,bool skipFirstConf){
    PRECONDITION(inStream,"no stream");
    std::string tempStr;
    std::vector<std::string> splitText;

    // format line:
    line++;
    tempStr = getLine(inStream);
    if(inStream->eof()){
      return NULL;
    }
    // comment line:
    line++;
    tempStr = getLine(inStream);
    if(inStream->eof()){
      return NULL;
    }
    // optional name line:
    line++;
    tempStr = getLine(inStream);
    if(inStream->eof()){
      return NULL;
    }
    RWMol *res = new RWMol();
    if(tempStr.size()>=4 && tempStr.substr(0,4)=="NAME"){
      tempStr = boost::trim_copy(tempStr.substr(4,tempStr.size()-4));
      res->setProp("_Name",tempStr);
      line++;
      tempStr = getLine(inStream);
      if(inStream->eof()){
        return res;
      }
    }
    if(tempStr.size()>=4 && tempStr.substr(0,4)=="PROP"){
      line++;
      tempStr = getLine(inStream);
      if(inStream->eof()){
        return res;
      }
    }
    
    // we're at the counts line:
    boost::split(splitText,tempStr,boost::is_any_of(" \t"),boost::token_compress_on);
    unsigned int nAtoms,nBonds;
    nAtoms = FileParserUtils::stripSpacesAndCast<unsigned int>(splitText[0]);
    nBonds = FileParserUtils::stripSpacesAndCast<unsigned int>(splitText[1]);

    Conformer *conf = new Conformer(nAtoms);
    conf->setId(0);
    for(unsigned int i=0;i<nAtoms;++i){
      line++;
      tempStr = getLine(inStream);
      if(inStream->eof()){
        throw FileParseException("EOF hit while reading atoms.") ;
      }
      ParseTPLAtomLine(tempStr,line,res,conf);
    }
    res->addConformer(conf, true);

    for(unsigned int i=0;i<nBonds;++i){
      line++;
      tempStr = getLine(inStream);
      if(inStream->eof()){
        throw FileParseException("EOF hit while reading bonds.") ;
      }
      ParseTPLBondLine(tempStr,line,res);
    }

    line++;
    tempStr = getLine(inStream);
    if(inStream->eof()){
      return res;
    }
    unsigned int nConfs=0;
    if(tempStr.size()>=5 && tempStr.substr(0,5)=="CONFS"){
      boost::split(splitText,tempStr,boost::is_any_of(" \t"),boost::token_compress_on);
      nConfs = FileParserUtils::stripSpacesAndCast<unsigned int>(splitText[1]);
    }
    for(unsigned int i=0;i<nConfs;++i){
      Conformer *conf=ParseConfData(inStream,line,res,i+1);
      if(i>0 || !skipFirstConf){
        conf->setId(i+1);
        res->addConformer(conf, true);
      } else {
        delete conf;
      }
      // there should be a blank line:
      line++;
      tempStr = getLine(inStream);
      boost::trim(tempStr);
      if(!inStream->eof() && tempStr != ""){
        throw FileParseException("Found a non-blank line between conformers.") ;
      }
    }
    if(sanitize && res){
      MolOps::sanitizeMol(*res);
    }
    
    return res;
  }

  //------------------------------------------------
  //
  //  Read a molecule from a file
  //
  //------------------------------------------------
  RWMol *TPLFileToMol(std::string fName, bool sanitize,bool skipFirstConf){
    std::ifstream inStream(fName.c_str());
    if (!inStream || (inStream.bad()) ) {
      std::ostringstream errout;
      errout << "Bad input file " << fName;
      throw BadFileException(errout.str());
    }
    RWMol *res=NULL;
    if(!inStream.eof()){
      unsigned int line = 0;
      res=TPLDataStreamToMol(&inStream, line, sanitize,skipFirstConf);
    }
    return res;
  }    

#if 0  

  RWMol *MolDataStreamToMol(std::istream &inStream, unsigned int &line,
                            bool sanitize){
    return MolDataStreamToMol(&inStream,line,sanitize);
  };
  //------------------------------------------------
  //
  //  Read a molecule from a string
  //
  //------------------------------------------------
  RWMol *MolBlockToMol(const std::string &molBlock, bool sanitize){
    std::istringstream inStream(molBlock);
    RWMol *res=NULL;
    unsigned int line = 0;
    return MolDataStreamToMol(inStream, line, sanitize);
  }    

#endif

}

