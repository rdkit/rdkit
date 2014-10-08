// $Id$
//
//  Copyright (C) 2002-2014 Greg Landrum and Rational Discovery LLC
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
#include <GraphMol/RDKitQueries.h>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/RDLog.h>

#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/LocaleSwitcher.h>
#include <typeinfo>
#include <exception>
#include <sstream>
#include <locale>
#include <stdlib.h>


namespace RDKit{
  class MolFileUnhandledFeatureException : public std::exception {
  public:
    //! construct with an error message
    explicit MolFileUnhandledFeatureException(const char *msg) : _msg(msg) {};
    //! construct with an error message
    explicit MolFileUnhandledFeatureException(const std::string msg) : _msg(msg) {};
    //! get the error message
    const char *message () const { return _msg.c_str(); };
    ~MolFileUnhandledFeatureException () throw () {};
  private:
    std::string _msg;
  };

  namespace FileParserUtils {
    int toInt(const std::string &input,bool acceptSpaces){
      int res=0;
      // don't need to worry about locale stuff here because
      // we're not going to have delimiters
      res=strtol(input.c_str(),NULL,10);
      if(!res && !acceptSpaces && input[0]==' '){
	std::string trimmed=boost::trim_copy(input);
	if(trimmed.length()==0) throw boost::bad_lexical_cast();
      }
      return res;
    }

    double toDouble(const std::string &input,bool acceptSpaces){
      double res=atof(input.c_str());
      if(res==0.0 && !acceptSpaces && input[0]==' '){
	std::string trimmed=boost::trim_copy(input);
	if(trimmed.length()==0) throw boost::bad_lexical_cast();
      }
      return res;
    }

    std::string getV3000Line(std::istream *inStream,unsigned int &line){
      // FIX: technically V3K blocks are case-insensitive. We should really be
      // up-casing everything here.
      PRECONDITION(inStream,"bad stream");
      std::string res,tempStr;

      ++line;
      tempStr = getLine(inStream);
      if(tempStr.size()<7 || tempStr.substr(0,7) != "M  V30 "){
        std::ostringstream errout;
        errout << "Line "<<line<<" does not start with 'M  V30 '"<<std::endl;
        throw FileParseException(errout.str()) ;
      }
      // FIX: do we need to handle trailing whitespace after a -?
      while(tempStr[tempStr.length()-1]=='-'){
        // continuation character, append what we read:
        res += tempStr.substr(7,tempStr.length()-8);
        // and then read another line: 
        ++line;
        tempStr = getLine(inStream);
        if(tempStr.size()<7 || tempStr.substr(0,7) != "M  V30 "){
          std::ostringstream errout;
          errout << "Line "<<line<<" does not start with 'M  V30 '"<<std::endl;
          throw FileParseException(errout.str()) ;
        }
      }
      res += tempStr.substr(7,tempStr.length()-7);
     
      return res;
    }

    Atom *replaceAtomWithQueryAtom(RWMol *mol,Atom *atom){
      PRECONDITION(mol,"bad molecule");
      PRECONDITION(atom,"bad atom");
      if(atom->hasQuery()) return atom;

      QueryAtom qa(*atom);
      unsigned int idx=atom->getIdx();

      if(atom->getFormalCharge()!=0){
	qa.expandQuery(makeAtomFormalChargeQuery(atom->getFormalCharge()));
      }
      if(atom->hasProp("_hasMassQuery")){
	qa.expandQuery(makeAtomMassQuery(static_cast<int>(atom->getMass())));
      }
      mol->replaceAtom(idx,&qa);
      return mol->getAtomWithIdx(idx);
    }

  }
  using RDKit::FileParserUtils::getV3000Line;

  namespace {
    void completeQueryAndChildren(ATOM_EQUALS_QUERY *query,Atom *tgt,int magicVal){
      PRECONDITION(query,"no query");
      PRECONDITION(tgt,"no atom");
      if(query->getVal()==magicVal){
        int tgtVal=query->getDataFunc()(tgt);
        query->setVal(tgtVal);
      }
      QueryAtom::QUERYATOM_QUERY::CHILD_VECT_CI childIt;
      for(childIt=query->beginChildren();childIt!=query->endChildren();++childIt){
        completeQueryAndChildren((ATOM_EQUALS_QUERY *)(childIt->get()),tgt,magicVal);
      }
    }
    void CompleteMolQueries(RWMol *mol,int magicVal=-0xDEADBEEF){
      for (ROMol::AtomIterator ai=mol->beginAtoms();
           ai != mol->endAtoms(); ++ai){
        if((*ai)->hasQuery()){
          ATOM_EQUALS_QUERY *query=static_cast<ATOM_EQUALS_QUERY *>((*ai)->getQuery());
          completeQueryAndChildren(query,*ai,magicVal);
        }
      }
    }
 

    //*************************************
    //
    // Every effort has been made to adhere to MDL's standard
    // for mol files
    //  
    //*************************************

    void ParseOldAtomList(RWMol *mol,const std::string &text,unsigned int line){
      PRECONDITION(mol,"bad mol");
      unsigned int idx;
      try {
        idx = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(0,3))-1;
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(0,3) << " to int on line "<<line;
        throw FileParseException(errout.str()) ;
      }

      RANGE_CHECK(0,idx,mol->getNumAtoms()-1);
      QueryAtom a(*(mol->getAtomWithIdx(idx)));

      ATOM_OR_QUERY *q = new ATOM_OR_QUERY;
      q->setDescription("AtomOr");
    
      switch(text[4]){
      case 'T':
        q->setNegation(true);
        break;
      case 'F':
        q->setNegation(false);
        break;
      default:
        std::ostringstream errout;
        errout << "Unrecognized atom-list query modifier: " << text[14]<<" on line "<<line;
        throw FileParseException(errout.str()) ;
      }          
    
      int nQueries;
      try {
        nQueries = FileParserUtils::toInt(text.substr(9,1));
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(9,1) << " to int on line "<<line;
        throw FileParseException(errout.str()) ;
      }

      RANGE_CHECK(0,nQueries,5);
      for(int i=0;i<nQueries;i++){
        int pos = 11+i*4;
        int atNum;
        try {
          atNum = FileParserUtils::toInt(text.substr(pos,3));
        }
        catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << text.substr(pos,3) << " to int on line "<<line;
          throw FileParseException(errout.str()) ;
        }
        RANGE_CHECK(0,atNum,200);  // goofy!
        q->addChild(QueryAtom::QUERYATOM_QUERY::CHILD_TYPE(makeAtomNumQuery(atNum)));
        if(!i) a.setAtomicNum(atNum);
      }
    
      a.setQuery(q);
      mol->replaceAtom(idx,&a); 
    };
  
    void ParseChargeLine(RWMol *mol, const std::string &text,bool firstCall,unsigned int line) {
      PRECONDITION(mol,"bad mol");
      PRECONDITION(text.substr(0,6)==std::string("M  CHG"),"bad charge line");
    

      // if this line is specified all the atom other than those specified
      // here should carry a charge of 0; but we should only do this once:
      if(firstCall){
        for (ROMol::AtomIterator ai = mol->beginAtoms();
             ai != mol->endAtoms(); ++ai) {
          (*ai)->setFormalCharge(0);
        }
      }

      int ie, nent;
      try {
        nent = FileParserUtils::toInt(text.substr(6,3));
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(6,3) << " to int on line "<<line;
        throw FileParseException(errout.str()) ;
      }
      int spos = 9;
      for (ie = 0; ie < nent; ie++) {
        int aid, chg;
        try {
          aid = FileParserUtils::toInt(text.substr(spos,4));
          spos += 4;
          chg = FileParserUtils::toInt(text.substr(spos,4));
          spos += 4;
          mol->getAtomWithIdx(aid-1)->setFormalCharge(chg);
        }
        catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << text.substr(spos,4) << " to int on line "<<line;
          throw FileParseException(errout.str()) ;
        }
      }
    }

    bool SGroupOK(std::string typ){
      const char *cfailTyps[11] = {
        // polymer sgroups:
        "SRU","MON","COP","CRO","GRA","MOD","MER","ANY",
        // formulations/mixtures:
        "COM","MIX","FOR"
      };
      std::vector<std::string> failTyps(cfailTyps,cfailTyps+11);
      return std::find(failTyps.begin(),failTyps.end(),typ)==failTyps.end();
    }
    
    void ParseSGroup2000STYLine(RWMol *mol, const std::string &text,unsigned int line){
      PRECONDITION(mol,"bad mol");
      PRECONDITION(text.substr(0,6)==std::string("M  STY"),"bad STY line");

      int nent;
      try {
        nent = FileParserUtils::toInt(text.substr(6,3));
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(6,3) << " to int on line "<<line;
        throw FileParseException(errout.str()) ;
      }

      unsigned int spos = 9;
      for (int ie = 0; ie < nent; ie++) {
        if(text.size()<spos+8){
          std::ostringstream errout;
          errout << "SGroup line too short: '" << text<<"' on line "<<line;
          throw FileParseException(errout.str()) ;
        }
#if 0        
        int nbr;
        try {
          nbr = FileParserUtils::toInt(text.substr(spos,4));
        } 
        catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << text.substr(spos,3) << " to int on line "<<line;
          throw FileParseException(errout.str()) ;
        }
#endif
        spos += 4;
        std::string typ = text.substr(spos+1,3);
        if(!SGroupOK(typ)){
          std::ostringstream errout;
          errout << "S group "<<typ;
          throw MolFileUnhandledFeatureException(errout.str()) ;
        } else {
          BOOST_LOG(rdWarningLog) << " S group " << typ <<" ignored on line "<<line<<std::endl;
        }
        spos += 4;
      }
    }

    void ParseRadicalLine(RWMol *mol, const std::string &text,bool firstCall,unsigned int line) {
      PRECONDITION(mol,"bad mol");
      PRECONDITION(text.substr(0,6)==std::string("M  RAD"),"bad charge line");

      // if this line is specified all the atom other than those specified
      // here should carry a charge of 0; but we should only do this once:
      if(firstCall){
        for (ROMol::AtomIterator ai = mol->beginAtoms();
             ai != mol->endAtoms(); ++ai) {
          (*ai)->setFormalCharge(0);
        }
      }

      int ie, nent;
      try {
        nent = FileParserUtils::toInt(text.substr(6,3));
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(6,3) << " to int on line "<<line;
        throw FileParseException(errout.str()) ;
      }
      int spos = 9;
      for (ie = 0; ie < nent; ie++) {
        int aid, rad;
        std::ostringstream errout;
      
        try {
          aid = FileParserUtils::toInt(text.substr(spos,4));
          spos += 4;
          rad = FileParserUtils::toInt(text.substr(spos,4));
          spos += 4;

          switch(rad) {
          case 1:
            mol->getAtomWithIdx(aid-1)->setNumRadicalElectrons(2);
            break;
          case 2:
            mol->getAtomWithIdx(aid-1)->setNumRadicalElectrons(1);
            break;
          case 3:
            mol->getAtomWithIdx(aid-1)->setNumRadicalElectrons(2);
            break;
          default:
            errout << "Unrecognized radical value " << rad << " for atom "<< aid-1 << " on line "<<line<<std::endl;
            throw FileParseException(errout.str()) ;
          }
        }
        catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << text.substr(spos,4) << " to int on line "<<line;
          throw FileParseException(errout.str()) ;
        }
      }
    }

    void ParseIsotopeLine(RWMol *mol, const std::string &text,unsigned int line){
      PRECONDITION(mol,"bad mol");
      PRECONDITION(text.substr(0,6)==std::string("M  ISO"),"bad isotope line");
    
      unsigned int nent;
      try {
        nent = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(6,3));
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(6,3) << " to int on line "<<line;
        throw FileParseException(errout.str()) ;
      }
      unsigned int spos = 9;
      for (unsigned int ie = 0; ie < nent; ie++) {
        unsigned int aid;
        try {
          aid = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(spos,4));
          spos += 4;
          Atom *atom=mol->getAtomWithIdx(aid-1); 
          if(text.size()>=spos+4 && text.substr(spos,4)!="    "){
            int isotope = FileParserUtils::toInt(text.substr(spos,4));
            if(isotope<0){
              BOOST_LOG(rdErrorLog) << " atom "<<aid<<" has a negative isotope value. line:  "<<line<<std::endl;
            } else {
              atom->setIsotope(isotope);
              spos += 4;
            }
          } else {
            atom->setMass(PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum()));
          }
        }
        catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << text.substr(spos,4) << " to int on line "<<line;
          throw FileParseException(errout.str()) ;
        }
      }
    
    }

    void ParseSubstitutionCountLine(RWMol *mol, const std::string &text,unsigned int line){
      PRECONDITION(mol,"bad mol");
      PRECONDITION(text.substr(0,6)==std::string("M  SUB"),"bad SUB line");
    
      unsigned int nent;
      try {
        nent = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(6,3));
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(6,3) << " to int on line "<<line;
        throw FileParseException(errout.str()) ;
      }
      unsigned int spos = 9;
      for (unsigned int ie = 0; ie < nent; ie++) {
        unsigned int aid;
        int count;
        try {
          aid = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(spos,4));
          spos += 4;
          Atom *atom=mol->getAtomWithIdx(aid-1); 
          if(text.size()>=spos+4 && text.substr(spos,4)!="    "){
            count = FileParserUtils::toInt(text.substr(spos,4));
            if(count==0) continue;
            ATOM_EQUALS_QUERY *q=makeAtomExplicitDegreeQuery(0);
            switch(count){
            case -1:
              q->setVal(0);break;
            case -2:
              q->setVal(atom->getDegree());break;
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
              q->setVal(count);break;
            case 6:
              BOOST_LOG(rdWarningLog) << " atom degree query with value 6 found. This will not match degree >6. The MDL spec says it should.  line: "<<line;
              q->setVal(6);break;
            default:
              std::ostringstream errout;
              errout << "Value " << count << " is not supported as a degree query. line: "<<line;
              throw FileParseException(errout.str()) ;
            }
            if(!atom->hasQuery()){
	      atom=FileParserUtils::replaceAtomWithQueryAtom(mol,atom);
            }
            atom->expandQuery(q,Queries::COMPOSITE_AND);
            spos += 4;
          }
        }
        catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << text.substr(spos,4) << " to int on line "<<line;
          throw FileParseException(errout.str()) ;
        }
      }
    }

    void ParseUnsaturationLine(RWMol *mol, const std::string &text,unsigned int line){
      PRECONDITION(mol,"bad mol");
      PRECONDITION(text.substr(0,6)==std::string("M  UNS"),"bad UNS line");
    
      unsigned int nent;
      try {
        nent = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(6,3));
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(6,3) << " to int on line "<<line;
        throw FileParseException(errout.str()) ;
      }
      unsigned int spos = 9;
      for (unsigned int ie = 0; ie < nent; ie++) {
        unsigned int aid;
        int count;
        try {
          aid = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(spos,4));
          spos += 4;
          Atom *atom=mol->getAtomWithIdx(aid-1); 
          if(text.size()>=spos+4 && text.substr(spos,4)!="    "){
            count = FileParserUtils::toInt(text.substr(spos,4));
            if(count==0){
              continue;
            } else if(count==1){
              ATOM_EQUALS_QUERY *q=makeAtomUnsaturatedQuery();
              if(!atom->hasQuery()){
		atom=FileParserUtils::replaceAtomWithQueryAtom(mol,atom);
              }
              atom->expandQuery(q,Queries::COMPOSITE_AND);
            } else {
              std::ostringstream errout;
              errout << "Value " << count << " is not supported as an unsaturation query (only 0 and 1 are allowed). line: "<<line;
              throw FileParseException(errout.str()) ;
            }
          }
        }catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << text.substr(spos,4) << " to int on line "<<line;
          throw FileParseException(errout.str()) ;
        }

      }
    }

    void ParseRingBondCountLine(RWMol *mol, const std::string &text,unsigned int line){
      PRECONDITION(mol,"bad mol");
      PRECONDITION(text.substr(0,6)==std::string("M  RBC"),"bad RBC line");
    
      unsigned int nent;
      try {
        nent = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(6,3));
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(6,3) << " to int on line "<<line;
        throw FileParseException(errout.str()) ;
      }
      unsigned int spos = 9;
      for (unsigned int ie = 0; ie < nent; ie++) {
        unsigned int aid;
        int count;
        try {
          aid = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(spos,4));
          spos += 4;
          Atom *atom=mol->getAtomWithIdx(aid-1); 
          if(text.size()>=spos+4 && text.substr(spos,4)!="    "){
            count = FileParserUtils::toInt(text.substr(spos,4));
            if(count==0) continue;
            ATOM_EQUALS_QUERY *q=makeAtomRingBondCountQuery(0);
            switch(count){
            case -1:
              q->setVal(0);break;
            case -2:
              q->setVal(-0xDEADBEEF);
              mol->setProp("_NeedsQueryScan",1);
              break;
            case 1:
            case 2:
            case 3:
              q->setVal(count);break;
            case 4:
              delete q;
              q = static_cast<ATOM_EQUALS_QUERY *>(new ATOM_LESSEQUAL_QUERY);
              q->setVal(4);
              q->setDescription("AtomRingBondCount");
              q->setDataFunc(queryAtomRingBondCount);
              break;
            default:
              std::ostringstream errout;
              errout << "Value " << count << " is not supported as a ring-bond count query. line: "<<line;
              throw FileParseException(errout.str()) ;
            }
            if(!atom->hasQuery()){
	      atom=FileParserUtils::replaceAtomWithQueryAtom(mol,atom);
            }
            atom->expandQuery(q,Queries::COMPOSITE_AND);
            spos += 4;
          }
        }
        catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << text.substr(spos,4) << " to int on line "<<line;
          throw FileParseException(errout.str()) ;
        }
      }
    }

    void ParseZCHLine(RWMol *mol, const std::string &text,unsigned int line){
      // part of Alex Clark's ZBO proposal
      // from JCIM 51:3149-57 (2011)
      PRECONDITION(mol,"bad mol");
      PRECONDITION(text.substr(0,6)==std::string("M  ZCH"),"bad ZCH line");
    
      unsigned int nent;
      try {
        nent = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(6,3));
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(6,3) << " to int on line "<<line;
        throw FileParseException(errout.str()) ;
      }
      unsigned int spos = 9;
      for (unsigned int ie = 0; ie < nent; ie++) {
        unsigned int aid=0;
        int val=0;
        try {
          aid = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(spos,4));
          spos += 4;
          if(text.size()>=spos+4 && text.substr(spos,4)!="    "){
            val = FileParserUtils::stripSpacesAndCast<int>(text.substr(spos,4));
          }
          if(!aid || aid>mol->getNumAtoms() || aid==0 ){
            std::ostringstream errout;
            errout << "Bad ZCH specification on line "<<line;
            throw FileParseException(errout.str()) ;
          }
          spos += 4;
          --aid;
          Atom *atom=mol->getAtomWithIdx(aid);
          if(!atom){
            std::ostringstream errout;
            errout << "Atom "<<aid<<" from ZCH specification on line "<<line<<" not found";
            throw FileParseException(errout.str()) ;
          } else {
            atom->setFormalCharge(val);
          }
        }
        catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << text.substr(spos,4) << " to int on line "<<line;
          throw FileParseException(errout.str()) ;
        }
      }
    }
    void ParseHYDLine(RWMol *mol, const std::string &text,unsigned int line){
      // part of Alex Clark's ZBO proposal
      // from JCIM 51:3149-57 (2011)
      PRECONDITION(mol,"bad mol");
      PRECONDITION(text.substr(0,6)==std::string("M  HYD"),"bad HYD line");
    
      unsigned int nent;
      try {
        nent = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(6,3));
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(6,3) << " to int on line "<<line;
        throw FileParseException(errout.str()) ;
      }
      unsigned int spos = 9;
      for (unsigned int ie = 0; ie < nent; ie++) {
        unsigned int aid=0;
        int val=-1;
        try {
          aid = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(spos,4));
          spos += 4;
          if(text.size()>=spos+4 && text.substr(spos,4)!="    "){
            val = FileParserUtils::stripSpacesAndCast<int>(text.substr(spos,4));
          }
          if(!aid || aid>mol->getNumAtoms() || aid==0 ){
            std::ostringstream errout;
            errout << "Bad HYD specification on line "<<line;
            throw FileParseException(errout.str()) ;
          }
          spos += 4;
          --aid;
          Atom *atom=mol->getAtomWithIdx(aid);
          if(!atom){
            std::ostringstream errout;
            errout << "Atom "<<aid<<" from HYD specification on line "<<line<<" not found";
            throw FileParseException(errout.str()) ;
          } else {
            if(val >=0 ){
              atom->setProp("_ZBO_H",true);
              atom->setNumExplicitHs(val);
            }
          }
        }
        catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << text.substr(spos,4) << " to int on line "<<line;
          throw FileParseException(errout.str()) ;
        }
      }
    }
    void ParseZBOLine(RWMol *mol, const std::string &text,unsigned int line){
      // part of Alex Clark's ZBO proposal
      // from JCIM 51:3149-57 (2011)
      PRECONDITION(mol,"bad mol");
      PRECONDITION(text.substr(0,6)==std::string("M  ZBO"),"bad ZBO line");
    
      unsigned int nent;
      try {
        nent = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(6,3));
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(6,3) << " to int on line "<<line;
        throw FileParseException(errout.str()) ;
      }
      unsigned int spos = 9;
      for (unsigned int ie = 0; ie < nent; ie++) {
        unsigned int bid=0;
        unsigned int order=0;
        try {
          bid = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(spos,4));
          spos += 4;
          if(text.size()>=spos+4 && text.substr(spos,4)!="    "){
            order = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(spos,4));
          }
          if(!bid || bid>mol->getNumBonds() || bid==0){
            std::ostringstream errout;
            errout << "Bad ZBO specification on line "<<line;
            throw FileParseException(errout.str()) ;
          }
          spos += 4;
          --bid;
          Bond *bnd=mol->getBondWithIdx(bid);
          if(!bnd){
            std::ostringstream errout;
            errout << "Bond "<<bid<<" from ZBO specification on line "<<line<<" not found";
            throw FileParseException(errout.str()) ;
          } else {
            if(order==0){
              bnd->setBondType(Bond::ZERO);
            } else {
              bnd->setBondType(static_cast<Bond::BondType>(order));
            }
          }
        }
        catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << text.substr(spos,4) << " to int on line "<<line;
          throw FileParseException(errout.str()) ;
        }
      }
    }

    void ParseNewAtomList(RWMol *mol,const std::string &text,unsigned int line){
      if(text.size()<15){
        std::ostringstream errout;
        errout << "Atom list line too short: '"<<text<<"'";
        throw FileParseException(errout.str()) ;
      }
      PRECONDITION(mol,"bad mol");
      PRECONDITION(text.substr(0,6)==std::string("M  ALS"),"bad atom list line");
    
      unsigned int idx;
      try {
        idx = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(7,3))-1;
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(7,3) << " to int on line "<<line;
        throw FileParseException(errout.str()) ;
      }
      RANGE_CHECK(0,idx,mol->getNumAtoms()-1);
      QueryAtom *a=0;
    
      int nQueries;
      try {
        nQueries = FileParserUtils::toInt(text.substr(10,3));
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(10,3) << " to int on line "<<line;
        throw FileParseException(errout.str()) ;
      }

      ASSERT_INVARIANT(nQueries>0,"no queries provided");
      for(unsigned int i=0;i<static_cast<unsigned int>(nQueries);i++){
        unsigned int pos = 16+i*4;
        if(text.size()<pos+4){
          std::ostringstream errout;
          errout << "Atom list line too short: '"<<text<<"' on line "<<line;
          throw FileParseException(errout.str()) ;
        }

        std::string atSymb = text.substr(pos,4);
        atSymb.erase(atSymb.find(" "),atSymb.size());
        int atNum = PeriodicTable::getTable()->getAtomicNumber(atSymb);
        if(!i){
          a = new QueryAtom(*(mol->getAtomWithIdx(idx)));
          // replace the query:
          Atom::QUERYATOM_QUERY *oq=a->getQuery();
          a->setAtomicNum(atNum);
          a->setQuery(makeAtomNumQuery(atNum));
          delete oq;
        } else {
          a->expandQuery(makeAtomNumQuery(atNum),Queries::COMPOSITE_OR,true);
        }
      }
      ASSERT_INVARIANT(a,"no atom built");
      
      switch(text[14]){
      case 'T':
        a->getQuery()->setNegation(true);
        break;
      case 'F':
        a->getQuery()->setNegation(false);
        break;
      default:
        std::ostringstream errout;
        errout << "Unrecognized atom-list query modifier: " << text[14]<<" on line "<<line;
        throw FileParseException(errout.str()) ;
      }          

      mol->replaceAtom(idx,a); 
    };
  
    void ParseV3000RGroups(RWMol *mol,Atom *&atom,const std::string &text,unsigned int line){
      PRECONDITION(mol,"bad mol");
      PRECONDITION(atom,"bad atom");
      if(text[0]!='('||text[text.size()-1]!=')'){
        std::ostringstream errout;
        errout << "Bad RGROUPS specification " << text << " on line "<<line<<". Missing parens.";
        throw FileParseException(errout.str()) ;
      }
      std::vector<std::string> splitToken;
      std::string resid=text.substr(1,text.size()-2);
      boost::split(splitToken,resid,boost::is_any_of(" "));
      if(splitToken.size()<1){
        std::ostringstream errout;
        errout << "Bad RGROUPS specification " << text << " on line "<<line<<". Missing values.";
        throw FileParseException(errout.str()) ;
      }
      unsigned int nRs;
      try {
        nRs = FileParserUtils::stripSpacesAndCast<unsigned int>(splitToken[0]);
      } catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << splitToken[0] << " to int on line"<<line;
        throw FileParseException(errout.str()) ;
      }
      if(splitToken.size()<nRs+1){
        std::ostringstream errout;
        errout << "Bad RGROUPS specification " << text << " on line "<<line<<". Not enough values.";
        throw FileParseException(errout.str()) ;
      }
      for(unsigned int i=0;i<nRs;++i){
        unsigned int rLabel;
        try {
          rLabel = FileParserUtils::stripSpacesAndCast<unsigned int>(splitToken[i+1]);
        } catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << splitToken[i+1] << " to int on line"<<line;
          throw FileParseException(errout.str()) ;
        }
        atom=FileParserUtils::replaceAtomWithQueryAtom(mol,atom);
        atom->setProp("_MolFileRLabel",rLabel);
        std::string dLabel="R"+boost::lexical_cast<std::string>(rLabel);
        atom->setProp("dummyLabel",dLabel);
        atom->setIsotope(rLabel);
        atom->setQuery(makeAtomNullQuery());
      }
    }
    
    void ParseRGroupLabels(RWMol *mol,const std::string &text,unsigned int line){
      PRECONDITION(mol,"bad mol");
      PRECONDITION(text.substr(0,6)==std::string("M  RGP"),"bad R group label line");
    
      int nLabels;
      try {
        nLabels = FileParserUtils::toInt(text.substr(6,3));
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(6,3) << " to int on line "<<line;
        throw FileParseException(errout.str()) ;
      }

      for(int i=0;i<nLabels;i++){
        int pos = 10+i*8;
        unsigned int atIdx;
        try {
          atIdx = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(pos,3));
        }
        catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << text.substr(pos,3) << " to int on line "<<line;
          throw FileParseException(errout.str()) ;
        }
        unsigned int rLabel;
        try {
          rLabel = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(pos+4,3));
        }
        catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << text.substr(pos+4,3) << " to int on line "<<line;
          throw FileParseException(errout.str()) ;
        }
        atIdx-=1;
        if(atIdx>mol->getNumAtoms()){
          std::ostringstream errout;
          errout << "Attempt to set R group label on nonexistent atom " << atIdx<< " on line "<<line;
          throw FileParseException(errout.str()) ;
        }
        QueryAtom qatom(*(mol->getAtomWithIdx(atIdx)));
        qatom.setProp("_MolFileRLabel",rLabel);

        // set the dummy label so that this is shown correctly
        // in other pieces of the code :
        // (this was sf.net issue 3316600)
        std::string dLabel="R"+boost::lexical_cast<std::string>(rLabel);
        qatom.setProp("dummyLabel",dLabel);
        
        // the CTFile spec (June 2005 version) technically only allows
        // R labels up to 32. Since there are three digits, we'll accept
        // anything: so long as it's positive and less than 1000:
        if(rLabel>0 && rLabel<999){
          qatom.setIsotope(rLabel);
        }
        qatom.setQuery(makeAtomNullQuery());
        mol->replaceAtom(atIdx,&qatom); 
      }
    };
  
    void ParseAtomAlias(RWMol *mol,std::string text,const std::string &nextLine,unsigned int line){
      PRECONDITION(mol,"bad mol");
      PRECONDITION(text.substr(0,2)==std::string("A "),"bad atom alias line");
      
      unsigned int idx;
      try {
        idx = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(3,3))-1;
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(3,3) << " to int on line "<<line;
        throw FileParseException(errout.str()) ;
      }
      RANGE_CHECK(0,idx,mol->getNumAtoms()-1);
      Atom *at = mol->getAtomWithIdx(idx);
      at->setProp("molFileAlias",nextLine);
    };
  
    void ParseAtomValue(RWMol *mol,std::string text,unsigned int line){
      PRECONDITION(mol,"bad mol");
      PRECONDITION(text.substr(0,2)==std::string("V "),"bad atom value line");
      
      unsigned int idx;
      try {
        idx = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(3,3))-1;
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(3,3) << " to int on line"<<line;
        throw FileParseException(errout.str()) ;
      }
      RANGE_CHECK(0,idx,mol->getNumAtoms()-1);
      Atom *at = mol->getAtomWithIdx(idx);
      at->setProp("molFileValue",text.substr(7,text.length()-7));
    };

    Atom *ParseMolFileAtomLine(const std::string text, RDGeom::Point3D &pos,unsigned int line) {
      Atom *res = new Atom;
      std::string symb;
      int massDiff,chg,hCount;

      if(text.size()<34){
        std::ostringstream errout;
        errout << "Atom line too short: '"<<text<<"' on line "<<line;
        throw FileParseException(errout.str()) ;
      }

      try {
        pos.x = FileParserUtils::toDouble(text.substr(0,10));
        pos.y = FileParserUtils::toDouble(text.substr(10,10));
        pos.z = FileParserUtils::toDouble(text.substr(20,10));
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot process coordinates on line "<<line;
        throw FileParseException(errout.str()) ;
      }
      symb = text.substr(31,3);
      symb = symb.substr(0,symb.find(' '));
    
      // REVIEW: should we handle missing fields at the end of the line?
      massDiff=0;
      if(text.size()>=36 && text.substr(34,2)!=" 0"){
        try {
          massDiff = FileParserUtils::toInt(text.substr(34,2),true);
        }
        catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << text.substr(34,2) << " to into on line "<<line;
          throw FileParseException(errout.str()) ;
        }
      }    
      chg=0;
      if(text.size()>=39 && text.substr(36,3)!="  0"){
        try {
          chg = FileParserUtils::toInt(text.substr(36,3),true);
        }
        catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << text.substr(36,3) << " to int on line "<<line;
          throw FileParseException(errout.str()) ;
        }
      }
      hCount = 0;
      if(text.size()>=45 && text.substr(42,3)!="  0"){
        try {
          hCount = FileParserUtils::toInt(text.substr(42,3),true);
        }
        catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << text.substr(42,3) << " to int on line "<<line;
          throw FileParseException(errout.str()) ;
        }
      }
      if(symb=="L" || symb=="A" || symb=="Q" || symb=="*" || symb=="LP"
         || symb=="R" || symb=="R#" || (symb[0]=='R' && symb>="R0" && symb<="R99") ){
        if(symb=="A"||symb=="Q"||symb=="*"){
          QueryAtom *query=new QueryAtom(0);
          if(symb=="*"){
            // according to the MDL spec, these match anything
            query->setQuery(makeAtomNullQuery());
          } else if(symb=="Q"){
            ATOM_OR_QUERY *q = new ATOM_OR_QUERY;
            q->setDescription("AtomOr");
            q->setNegation(true);
            q->addChild(QueryAtom::QUERYATOM_QUERY::CHILD_TYPE(makeAtomNumQuery(6)));
            q->addChild(QueryAtom::QUERYATOM_QUERY::CHILD_TYPE(makeAtomNumQuery(1)));
            query->setQuery(q);
          } else if(symb=="A"){
            query->setQuery(makeAtomNumQuery(1));
            query->getQuery()->setNegation(true);
          }
          delete res;
          res=query;  
          // queries have no implicit Hs:
          res->setNoImplicit(true);
        } else {
          res->setAtomicNum(0);
        }
        if(massDiff==0&&symb[0]=='R'){
          if(symb.length()>1){
            std::string rlabel="";
            rlabel = symb.substr(1,symb.length()-1);
            int rnumber;
            try {
              rnumber = boost::lexical_cast<int>(rlabel);
            } catch (boost::bad_lexical_cast &) {
              rnumber=-1;
            }
            if(rnumber>=0) res->setIsotope(rnumber);
          }
        }
      } else if( symb=="D" ){  // mol blocks support "D" and "T" as shorthand... handle that.
        res->setAtomicNum(1); 
        res->setIsotope(2);
      } else if( symb=="T" ){  // mol blocks support "D" and "T" as shorthand... handle that.
        res->setAtomicNum(1);
        res->setIsotope(3);
      } else {
        res->setAtomicNum(PeriodicTable::getTable()->getAtomicNumber(symb));
        res->setMass(PeriodicTable::getTable()->getAtomicWeight(res->getAtomicNum()));
      }
    
      //res->setPos(pX,pY,pZ);
      if(chg!=0) res->setFormalCharge(4-chg);

      // FIX: this does not appear to be correct
      if(hCount==1){
        res->setNoImplicit(true);
      }
    
      if(massDiff!=0) {
        int defIso=PeriodicTable::getTable()->getMostCommonIsotope(res->getAtomicNum());
        int dIso=defIso+massDiff;
        if(dIso<0){
          BOOST_LOG(rdWarningLog) << " atom "<<res->getIdx()<<" has a negative isotope offset. line:  "<<line<<std::endl;
        }
        res->setIsotope(dIso);
        res->setMass(PeriodicTable::getTable()->getMassForIsotope(res->getAtomicNum(),
                                                                  dIso));
	res->setProp("_hasMassQuery",true);
      }
    
      if(text.size()>=42 && text.substr(39,3)!="  0"){
        int parity=0;
        try {
          parity = FileParserUtils::toInt(text.substr(39,3),true);
        }
        catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << text.substr(39,3) << " to int on line "<<line;
          throw FileParseException(errout.str()) ;
        }
        res->setProp("molParity",parity);
      }

      if(text.size()>=48 && text.substr(45,3)!="  0"){
        int stereoCare=0;
        try {
          stereoCare = FileParserUtils::toInt(text.substr(45,3),true);
        }
        catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << text.substr(45,3) << " to int on line "<<line;
          throw FileParseException(errout.str()) ;
        }
        res->setProp("molStereoCare",stereoCare);
      }
      if(text.size()>=51 && text.substr(48,3)!="  0"){
        int totValence=0;
        try {
          totValence= FileParserUtils::toInt(text.substr(48,3),true);
        }
        catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << text.substr(48,3) << " to int on line "<<line;
          throw FileParseException(errout.str()) ;
        }
        if(totValence!=0){
          // only set if it's a non-default value
          res->setProp("molTotValence",totValence);
        }
      }
      if(text.size()>=57 && text.substr(54,3)!="  0"){
        int rxnRole=0;
        try {
          rxnRole= FileParserUtils::toInt(text.substr(54,3),true);
        }
        catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << text.substr(54,3) << " to int on line "<<line;
          throw FileParseException(errout.str()) ;
        }
        if(rxnRole!=0){
          // only set if it's a non-default value
          res->setProp("molRxnRole",rxnRole);
        }
      }
      if(text.size()>=60 && text.substr(57,3)!="  0"){
        int rxnComponent=0;
        try {
          rxnComponent= FileParserUtils::toInt(text.substr(57,3),true);
        }
        catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << text.substr(57,3) << " to int on line "<<line;
          throw FileParseException(errout.str()) ;
        }
        if(rxnComponent!=0){
          // only set if it's a non-default value
          res->setProp("molRxnComponent",rxnComponent);
        }
      }
      if(text.size()>=63 && text.substr(60,3)!="  0"){
        int atomMapNumber=0;
        try {
          atomMapNumber = FileParserUtils::toInt(text.substr(60,3),true);
        }
        catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << text.substr(60,3) << " to int on line "<<line;
          throw FileParseException(errout.str()) ;
        }
        res->setProp("molAtomMapNumber",atomMapNumber);
      }
      if(text.size()>=66 && text.substr(63,3)!="  0"){
        int inversionFlag=0;
        try {
          inversionFlag= FileParserUtils::toInt(text.substr(63,3),true);
        }
        catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << text.substr(63,3) << " to int on line "<<line;
          throw FileParseException(errout.str()) ;
        }
        res->setProp("molInversionFlag",inversionFlag);
      }
      if(text.size()>=69 && text.substr(66,3)!="  0"){
        int exactChangeFlag=0;
        try {
          exactChangeFlag = FileParserUtils::toInt(text.substr(66,3),true);
        }
        catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Cannot convert " << text.substr(66,3) << " to int on line "<<line;
          throw FileParseException(errout.str()) ;
        }
        res->setProp("molExactChangeFlag",exactChangeFlag);
      }
      return res;
    };
  
    Bond *ParseMolFileBondLine(const std::string &text,unsigned int line){
      int idx1,idx2,bType,stereo;
      int spos = 0;

      if(text.size()<9){
        std::ostringstream errout;
        errout << "Bond line too short: '"<<text<<"' on line "<<line;
        throw FileParseException(errout.str()) ;
      }

      try {
        idx1 = FileParserUtils::toInt(text.substr(spos,3));
        spos += 3;
        idx2 = FileParserUtils::toInt(text.substr(spos,3));
        spos += 3;
        bType = FileParserUtils::toInt(text.substr(spos,3));  
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(spos,3) << " to int on line "<<line;
        throw FileParseException(errout.str()) ;
      }
    
      // adjust the numbering
      idx1--;idx2--;

      Bond::BondType type;
      Bond *res=0;  
      switch(bType){
      case 1: type = Bond::SINGLE;res = new Bond;break;
      case 2: type = Bond::DOUBLE;res = new Bond;break;
      case 3: type = Bond::TRIPLE;res = new Bond;break;
      case 4: type = Bond::AROMATIC;res = new Bond;break;
      case 0:
        type = Bond::UNSPECIFIED;
        res = new Bond;
        BOOST_LOG(rdWarningLog) << "bond with order 0 found on line "<<line<<". This is not part of the MDL specification."<<std::endl;
        break;
      default:
        type = Bond::UNSPECIFIED;
        // it's a query bond of some type
        res = new QueryBond;
        if(bType == 8){
          BOND_NULL_QUERY *q;
          q = makeBondNullQuery();
          res->setQuery(q);
        } else if (bType==5 || bType==6 || bType==7 ){
          BOND_OR_QUERY *q;
          q = new BOND_OR_QUERY;
          if(bType == 5){
            // single or double
            q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(makeBondOrderEqualsQuery(Bond::SINGLE)));
            q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(makeBondOrderEqualsQuery(Bond::DOUBLE)));
            q->setDescription("BondOr");
          } else if(bType == 6){
            // single or aromatic
            q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(makeBondOrderEqualsQuery(Bond::SINGLE)));
            q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(makeBondOrderEqualsQuery(Bond::AROMATIC)));      
            q->setDescription("BondOr");
          } else if(bType == 7){
            // double or aromatic
            q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(makeBondOrderEqualsQuery(Bond::DOUBLE)));
            q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(makeBondOrderEqualsQuery(Bond::AROMATIC)));
            q->setDescription("BondOr");
          }
          res->setQuery(q);
        } else {
          BOND_NULL_QUERY *q;
          q = makeBondNullQuery();
          res->setQuery(q);
          BOOST_LOG(rdWarningLog) << "unrecognized query bond type, " << bType <<", found on line "<<line<<". Using an \"any\" query."<<std::endl;          
        }
        break;
      }
      res->setBeginAtomIdx(idx1);
      res->setEndAtomIdx(idx2);
      res->setBondType(type);

      if( text.size() >= 12 && text.substr(9,3)!="  0"){
        try {
          stereo = FileParserUtils::toInt(text.substr(9,3));
          switch(stereo){
          case 0:
            res->setBondDir(Bond::NONE);
            break;
          case 1:
            res->setBondDir(Bond::BEGINWEDGE);
            break;
          case 6:
            res->setBondDir(Bond::BEGINDASH);
            break;
          case 3: // "either" double bond
            res->setBondDir(Bond::EITHERDOUBLE);
	    res->setStereo(Bond::STEREOANY);
	    break;
          case 4: // "either" single bond
            res->setBondDir(Bond::UNKNOWN);
            break;
          }
        } catch (boost::bad_lexical_cast) {
          ;
        }
      }
      if( text.size() >= 18 && text.substr(15,3)!="  0"){
        try {
          int topology = FileParserUtils::toInt(text.substr(15,3));
          if(topology){
            if(!res->hasQuery()){
              QueryBond *qBond=new QueryBond(*res);
              delete res;
              res = qBond;
            }
            BOND_EQUALS_QUERY *q=makeBondIsInRingQuery();
            switch(topology){
            case 1:
              break;
            case 2:
              q->setNegation(true);
              break;
            default:
              std::ostringstream errout;
              errout << "Unrecognized bond topology specifier: " << topology<<" on line "<<line;
              throw FileParseException(errout.str()) ;
            }
            res->expandQuery(q);          
          }
        } catch (boost::bad_lexical_cast) {
          ;
        }
      }
      if( text.size() >= 21 && text.substr(18,3)!="  0"){
        try {
          int reactStatus = FileParserUtils::toInt(text.substr(18,3));
          res->setProp("molReactStatus",reactStatus);
        } catch (boost::bad_lexical_cast) {
          ;
        }
      }
      return res;
    };  

    void ParseMolBlockAtoms(std::istream *inStream,unsigned int &line,
                           unsigned int nAtoms,RWMol *mol,Conformer *conf){
      PRECONDITION(inStream,"bad stream");
      PRECONDITION(mol,"bad molecule");
      PRECONDITION(conf,"bad conformer");
      for(unsigned int i=0;i<nAtoms;++i){
        ++line;
        std::string tempStr = getLine(inStream);
        if(inStream->eof()){
          throw FileParseException("EOF hit while reading atoms");
        }
        RDGeom::Point3D pos;
        Atom *atom = ParseMolFileAtomLine(tempStr, pos, line);
        unsigned int aid = mol->addAtom(atom,false,true);
        conf->setAtomPos(aid, pos);
      }
    }

    // returns whether or not any sign of chirality was detected
    void ParseMolBlockBonds(std::istream *inStream,unsigned int &line,
			    unsigned int nBonds,RWMol *mol,bool &chiralityPossible){
      PRECONDITION(inStream,"bad stream");
      PRECONDITION(mol,"bad molecule");
      for(unsigned int i=0;i<nBonds;++i){
        ++line;
        std::string tempStr = getLine(inStream);
        if(inStream->eof()){
          throw FileParseException("EOF hit while reading bonds");
        }
        Bond *bond = ParseMolFileBondLine(tempStr,line);
        // if we got an aromatic bond set the flag on the bond and the connected atoms
        if (bond->getBondType() == Bond::AROMATIC) {
          bond->setIsAromatic(true);
          mol->getAtomWithIdx(bond->getBeginAtomIdx())->setIsAromatic(true);
          mol->getAtomWithIdx(bond->getEndAtomIdx())->setIsAromatic(true);
        }
        // if the bond might have chirality info associated with it, set a flag:
        if(bond->getBondDir() != Bond::NONE && bond->getBondDir() != Bond::UNKNOWN){
          chiralityPossible=true;
        }
        mol->addBond(bond,true);
      }
    }

    bool ParseMolBlockProperties(std::istream *inStream,unsigned int &line,
                                 RWMol *mol){
      PRECONDITION(inStream,"bad stream");
      PRECONDITION(mol,"bad molecule");
      // older mol files can have an atom list block here
      std::string tempStr = getLine(inStream);
      ++line;
      if( tempStr[0] != 'M' && tempStr[0] != 'A'
          && tempStr[0] != 'V' && tempStr[0] != 'G' && tempStr[0] != 'S'){
        ParseOldAtomList(mol,tempStr,line);
      }

      bool fileComplete=false;
      bool firstChargeLine=true;
      std::string lineBeg=tempStr.substr(0,6);
      while(!inStream->eof() && lineBeg!="M  END" && tempStr.substr(0,4)!="$$$$"){
        if(tempStr[0]=='A'){
          line++;
          std::string nextLine = getLine(inStream);
          if(lineBeg!="M  END"){
            ParseAtomAlias(mol,tempStr,nextLine,line);
          }
        } else if(tempStr[0]=='G'){
          BOOST_LOG(rdWarningLog)<<" deprecated group abbreviation ignored on line "<<line<<std::endl;
          // we need to skip the next line, which holds the abbreviation:
          line++;
          tempStr = getLine(inStream);
        } else if(tempStr[0]=='V'){
          ParseAtomValue(mol,tempStr,line);
        } else if(lineBeg=="S  SKP") {
          int nToSkip=FileParserUtils::toInt(tempStr.substr(6,3));
          if(nToSkip<0) {
            std::ostringstream errout;
            errout << "negative skip value "<<nToSkip<<" on line "<<line;
            throw FileParseException(errout.str()) ;
          }
          for(unsigned int i=0;i<static_cast<unsigned int>(nToSkip);++i){
            ++line;
            tempStr=getLine(inStream);
          }
        }
        else if(lineBeg=="M  ALS") ParseNewAtomList(mol,tempStr,line);
        else if(lineBeg=="M  ISO") ParseIsotopeLine(mol,tempStr,line);
        else if(lineBeg=="M  RGP") ParseRGroupLabels(mol,tempStr,line);
        else if(lineBeg=="M  RBC") ParseRingBondCountLine(mol,tempStr,line);
        else if(lineBeg=="M  SUB") ParseSubstitutionCountLine(mol,tempStr,line);
        else if(lineBeg=="M  UNS") ParseUnsaturationLine(mol,tempStr,line);
        else if(lineBeg=="M  CHG") {
          ParseChargeLine(mol, tempStr,firstChargeLine,line);
          firstChargeLine=false;
        }
        else if(lineBeg=="M  RAD") {
          ParseRadicalLine(mol, tempStr,firstChargeLine,line);
          firstChargeLine=false;
        }
        else if(lineBeg=="M  STY") {
          ParseSGroup2000STYLine(mol, tempStr,line);
        }
        else if(lineBeg=="M  ZBO") ParseZBOLine(mol,tempStr,line);
        else if(lineBeg=="M  ZCH") ParseZCHLine(mol,tempStr,line);
        else if(lineBeg=="M  HYD") ParseHYDLine(mol,tempStr,line);
        line++;
        tempStr = getLine(inStream);
        lineBeg=tempStr.substr(0,6);
      }
      if(tempStr[0]=='M'&&tempStr.substr(0,6)=="M  END"){
        fileComplete=true;
      }
      return fileComplete;
    }

    Atom *ParseV3000AtomSymbol(std::string token,unsigned int &line){
      bool negate=false;
      boost::trim(token);
      std::string cpy=token;
      boost::to_upper(cpy);
      if(cpy.size()>3 && cpy.substr(0,3)=="NOT"){
	negate=true;
	token = token.substr(3,token.size()-3);
	boost::trim(token);
      }

      Atom *res=0;
      if(token[0]=='['){
        // atom list:
        if(token[token.length()-1]!=']'){
          std::ostringstream errout;
          errout << "Bad atom token '"<<token<<"' on line: "<<line;
          throw FileParseException(errout.str()) ;
        }
        token = token.substr(1,token.size()-2);

        std::vector<std::string> splitToken;
        boost::split(splitToken,token,boost::is_any_of(","));

        for(std::vector<std::string>::const_iterator stIt=splitToken.begin();
            stIt!=splitToken.end();++stIt){
          std::string atSymb=boost::trim_copy(*stIt);
          if(atSymb=="") continue;
          int atNum = PeriodicTable::getTable()->getAtomicNumber(atSymb);
          if(!res){
            res = new QueryAtom(atNum);
          } else {
            res->expandQuery(makeAtomNumQuery(atNum),Queries::COMPOSITE_OR,true);
          }
        }
        res->getQuery()->setNegation(negate);
      } else {
        if(negate) {
          std::ostringstream errout;
          errout << "NOT tokens only supported for atom lists. line "<<line;
          throw FileParseException(errout.str()) ;
        }
        // it's a normal CTAB atom symbol:
        // NOTE: "R" and "R0"-"R99" are not in the v3K CTAB spec, but we're going to support them anyway
        if(token=="R" || 
           (token[0]=='R' && token>="R0" && token<="R99") ||
           token=="R#" || token=="A" || token=="Q" || token=="*"){
          if(token=="A"||token=="Q"||token=="*"){
            res=new QueryAtom(0);
            if(token=="*"){
              // according to the MDL spec, these match anything
              res->setQuery(makeAtomNullQuery());
            } else if(token=="Q"){
              ATOM_OR_QUERY *q = new ATOM_OR_QUERY;
              q->setDescription("AtomOr");
              q->setNegation(true);
              q->addChild(QueryAtom::QUERYATOM_QUERY::CHILD_TYPE(makeAtomNumQuery(6)));
              q->addChild(QueryAtom::QUERYATOM_QUERY::CHILD_TYPE(makeAtomNumQuery(1)));
              res->setQuery(q);
            } else if(token=="A"){
              res->setQuery(makeAtomNumQuery(1));
              res->getQuery()->setNegation(true);
            }
            // queries have no implicit Hs:
            res->setNoImplicit(true);
          } else {
            res = new Atom(1);
            res->setAtomicNum(0);
          }
          if(token[0]=='R' && token>="R0" && token<="R99"){
            std::string rlabel="";
            rlabel = token.substr(1,token.length()-1);
            int rnumber;
            try {
              rnumber = boost::lexical_cast<int>(rlabel);
            } catch (boost::bad_lexical_cast &) {
              rnumber=-1;
            }
            if(rnumber>=0) res->setIsotope(rnumber);
          }
        } else if( token=="D" ){  // mol blocks support "D" and "T" as shorthand... handle that.
          res = new Atom(1);
          res->setIsotope(2);
        } else if( token=="T" ){  // mol blocks support "D" and "T" as shorthand... handle that.
          res = new Atom(1);
          res->setIsotope(3);
        } else {
          res = new Atom(PeriodicTable::getTable()->getAtomicNumber(token));
          res->setMass(PeriodicTable::getTable()->getAtomicWeight(res->getAtomicNum()));
        }
      }
      
      
      POSTCONDITION(res,"no atom built");
      return res;
    }

    bool splitAssignToken(const std::string &token,std::string &prop,std::string &val){
      std::vector<std::string> splitToken;
      boost::split(splitToken,token,
                   boost::is_any_of("="));
      if(splitToken.size()!=2){
        return false;
      }
      prop = splitToken[0];
      boost::to_upper(prop);
      val = splitToken[1];
      return true;
    }

    template <class T>
    void ParseV3000AtomProps(RWMol *mol,Atom *& atom,
			     typename T::iterator &token,const T &tokens,
                             unsigned int &line){
      PRECONDITION(mol,"bad molecule");
      PRECONDITION(atom,"bad atom");
      std::ostringstream errout;
      while(token!=tokens.end()){
        std::string prop,val;
        if(!splitAssignToken(*token,prop,val)){
          errout << "Invalid atom property: " << *token << " for atom "<< atom->getIdx()+1 <<" on line "<<line<< std::endl;
          throw FileParseException(errout.str()) ;
        }

        if(prop=="CHG"){
          int charge=FileParserUtils::toInt(val);
          if(!atom->hasQuery()) {
            atom->setFormalCharge(charge);
          } else {
            atom->expandQuery(makeAtomFormalChargeQuery(charge));
          }
        } else if(prop=="RAD"){
          // FIX handle queries here
          switch( FileParserUtils::toInt(val) ){
          case 0: break;
          case 1:
            atom->setNumRadicalElectrons(2);break;
          case 2:
            atom->setNumRadicalElectrons(1);break;
          case 3:
            atom->setNumRadicalElectrons(2);break;
          default:
            errout << "Unrecognized RAD value " << val << " for atom "<< atom->getIdx()+1 <<" on line "<<line<< std::endl;
            throw FileParseException(errout.str()) ;
          }
        } else if(prop=="MASS"){
          // the documentation for V3000 CTABs says that this should contain the "absolute atomic weight" (whatever that means).
          // Online examples seem to have integer (isotope) values and Marvin won't even read something that has a float.
          // We'll go with the int
          int v;
          double dv;
          try{
            v=FileParserUtils::toInt(val);
          } catch (boost::bad_lexical_cast &) {
            try{
              dv=FileParserUtils::toDouble(val);
              v = static_cast<int>(floor(dv));
            } catch (boost::bad_lexical_cast &){
              v=-1;
            }
          }
          if(v<0){
            errout << "Bad value for MASS :" << val << " for atom "<< atom->getIdx()+1 <<" on line "<<line << std::endl;
            throw FileParseException(errout.str()) ;
          } else {
	    if(!atom->hasQuery()) {
	      atom->setIsotope(v);
	    } else {
	      atom->expandQuery(makeAtomIsotopeQuery(v));
	    }
	  }
        } else if(prop=="CFG"){
          int cfg=FileParserUtils::toInt(val);
          switch(cfg){
          case 0: break;
          case 1:
          case 2:
          case 3:
            atom->setProp("molParity",cfg);
            break;
          default:
            errout << "Unrecognized CFG value : " << val << " for atom "<< atom->getIdx()+1 <<" on line " << line<< std::endl;
            throw FileParseException(errout.str()) ;
          }
        } else if(prop=="HCOUNT"){
	  if(val!="0"){
	    int hcount=FileParserUtils::toInt(val);
	    if(!atom->hasQuery()) {
	      atom=FileParserUtils::replaceAtomWithQueryAtom(mol,atom);
	    }
	    if(hcount==-1) hcount=0;
	    atom->expandQuery(makeAtomHCountQuery(hcount));
	  }
        } else if(prop=="UNSAT"){
	  if(val=="1"){
	    if(!atom->hasQuery()) {
	      atom=FileParserUtils::replaceAtomWithQueryAtom(mol,atom);
	    } 
	    atom->expandQuery(makeAtomUnsaturatedQuery());
	  }
        } else if(prop=="RBCNT"){
	  if(val!="0"){
	    int rbcount=FileParserUtils::toInt(val);
	    if(!atom->hasQuery()) {
	      atom=FileParserUtils::replaceAtomWithQueryAtom(mol,atom);
	    }
	    if(rbcount==-1) rbcount=0;
	    atom->expandQuery(makeAtomRingBondCountQuery(rbcount));
	  }
        } else if(prop=="VAL"){
	  if(val!="0"){
	    int totval=FileParserUtils::toInt(val);
	    atom->setProp("molTotValence",totval);
	  }
        } else if(prop=="RGROUPS"){
          ParseV3000RGroups(mol,atom,val,line);
          // FIX
        }
        ++token;
      }
    }

    void tokenizeV3000Line(std::string line,std::vector<std::string> &tokens){
      bool inQuotes=false,inParens=false;
      unsigned int start=0;
      unsigned int pos=0;
      while(pos<line.size()){
        if(line[pos]==' ' || line[pos]=='\t'){
          if(start == pos){
            ++start;
            ++pos;
          } else if( !inQuotes && !inParens){
            tokens.push_back(line.substr(start,pos-start));
            ++pos;
            start=pos;
          } else {
            ++pos;
          }
        } else if(line[pos]==')' && inParens){
          tokens.push_back(line.substr(start,pos-start+1));
          inParens=false;
          ++pos;
          start=pos;
        } else if(line[pos]=='(' && !inQuotes){
          inParens=true;
          ++pos;
        } else if(line[pos]=='"' && !inParens){
          if(pos+1<line.size() && line[pos+1]=='"'){
            pos+=2;
          } else if(inQuotes){
            // don't push on the quotes themselves
            tokens.push_back(line.substr(start+1,pos-start-1));
            ++pos;
            start=pos;
            inQuotes=false;
          } else {
            ++pos;
            inQuotes=true;
          }
        } else {
          ++pos;
        }
      }
      if(start!=pos){
        tokens.push_back(line.substr(start,line.size()-start));
      }
#if 0
      std::cerr<<"tokens: ";
      std::copy(tokens.begin(),tokens.end(),std::ostream_iterator<std::string>(std::cerr,"|"));
      std::cerr<<std::endl;
#endif
    }

    void ParseV3000AtomBlock(std::istream *inStream,unsigned int &line,
                             unsigned int nAtoms,RWMol *mol, Conformer *conf){
      PRECONDITION(inStream,"bad stream");
      PRECONDITION(nAtoms>0,"bad atom count");
      PRECONDITION(mol,"bad molecule");
      PRECONDITION(conf,"bad conformer");
      std::string tempStr;
      std::vector<std::string> splitLine;

      tempStr = getV3000Line(inStream,line);
      if(tempStr.length()<10 || tempStr.substr(0,10) != "BEGIN ATOM"){
        std::ostringstream errout;
        errout<<"BEGIN ATOM line not found on line "<<line;
        throw FileParseException(errout.str()) ;
      }
      for(unsigned int i=0;i<nAtoms;++i){

        tempStr = getV3000Line(inStream,line);
        std::string trimmed=boost::trim_copy(tempStr);

        std::vector<std::string> tokens;
        std::vector<std::string>::iterator token;

        tokenizeV3000Line(trimmed,tokens);
        token=tokens.begin();

        if(token==tokens.end()) {
          std::ostringstream errout;
          errout << "Bad atom line : '"<<tempStr<<"' on line"<<line;
          throw FileParseException(errout.str()) ;
        }
        unsigned int molIdx=atoi(token->c_str());

        // start with the symbol:
        ++token;
        if(token==tokens.end()) {
          std::ostringstream errout;
          errout << "Bad atom line : '"<<tempStr<<"' on line "<<line;
          throw FileParseException(errout.str()) ;
        }
        Atom *atom=ParseV3000AtomSymbol(*token,line);

        // now the position;
        RDGeom::Point3D pos;
        ++token;
        if(token==tokens.end()) {
          std::ostringstream errout;
          errout << "Bad atom line : '"<<tempStr<<"' on line "<<line;
          throw FileParseException(errout.str()) ;
        }
        pos.x=atof(token->c_str());
        ++token;
        if(token==tokens.end()) {
          std::ostringstream errout;
          errout << "Bad atom line : '"<<tempStr<<"' on line "<<line;
          throw FileParseException(errout.str()) ;
        }
        pos.y=atof(token->c_str());
        ++token;
        if(token==tokens.end()) {
          std::ostringstream errout;
          errout << "Bad atom line : '"<<tempStr<<"' on line "<<line;
          throw FileParseException(errout.str()) ;
        }
        pos.z=atof(token->c_str());
        // the map number:
        ++token;
        if(token==tokens.end()) {
          std::ostringstream errout;
          errout << "Bad atom line : '"<<tempStr<<"' on line "<<line;
          throw FileParseException(errout.str()) ;
        }
        int mapNum=atoi(token->c_str());
	if(mapNum>0){
	  atom->setProp("molAtomMapNumber",mapNum);
	}
        ++token;
        
        unsigned int aid=mol->addAtom(atom,false,true);

        // additional properties this may change the atom,
	// so be careful with it:
        ParseV3000AtomProps(mol,atom,token,tokens,line);

        mol->setAtomBookmark(atom,molIdx);
        conf->setAtomPos(aid,pos);
      }
      tempStr = getV3000Line(inStream,line);
      if(tempStr.length()<8 || tempStr.substr(0,8) != "END ATOM"){
        std::ostringstream errout;
        errout<<"END ATOM line not found on line "<<line;
        throw FileParseException(errout.str()) ;
      }

      if(mol->hasProp("_2DConf")){
        conf->set3D(false);
        mol->clearProp("_2DConf");
      } else if(mol->hasProp("_3DConf")){
        conf->set3D(true);
        mol->clearProp("_3DConf");
      }
    }
    void ParseV3000BondBlock(std::istream *inStream,unsigned int &line,
                             unsigned int nBonds,RWMol *mol,
			     bool &chiralityPossible){
      PRECONDITION(inStream,"bad stream");
      PRECONDITION(nBonds>0,"bad bond count");
      PRECONDITION(mol,"bad molecule");

      std::string tempStr;
      std::vector<std::string> splitLine;

      tempStr = getV3000Line(inStream,line);
      if(tempStr.length()<10 || tempStr.substr(0,10) != "BEGIN BOND"){
        throw FileParseException("BEGIN BOND line not found") ;
      }
      for(unsigned int i=0;i<nBonds;++i){
        tempStr = boost::trim_copy(getV3000Line(inStream,line));
        boost::split(splitLine,tempStr,
                     boost::is_any_of(" \t"),boost::token_compress_on);
        if(splitLine.size()<4){
          std::ostringstream errout;
          errout << "bond line "<<line<<" is too short";
          throw FileParseException(errout.str()) ;
        }
        Bond *bond;
        unsigned int bondIdx=atoi(splitLine[0].c_str());
        unsigned int bType=atoi(splitLine[1].c_str());
        unsigned int a1Idx=atoi(splitLine[2].c_str());
        unsigned int a2Idx=atoi(splitLine[3].c_str());

        switch(bType){
        case 1: bond = new Bond(Bond::SINGLE);break;
        case 2: bond = new Bond(Bond::DOUBLE);break;
        case 3: bond = new Bond(Bond::TRIPLE);break;
        case 4: bond = new Bond(Bond::AROMATIC);bond->setIsAromatic(true);break;
        case 0:
          bond = new Bond(Bond::UNSPECIFIED);
          BOOST_LOG(rdWarningLog) << "bond with order 0 found on line "<<line<<". This is not part of the MDL specification."<<std::endl;
          break;
        default:
          // it's a query bond of some type
          bond = new QueryBond;
          if(bType == 8){
            BOND_NULL_QUERY *q;
            q = makeBondNullQuery();
            bond->setQuery(q);
          } else if (bType==5 || bType==6 || bType==7 ){
            BOND_OR_QUERY *q;
            q = new BOND_OR_QUERY;
            if(bType == 5){
              // single or double
              q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(makeBondOrderEqualsQuery(Bond::SINGLE)));
              q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(makeBondOrderEqualsQuery(Bond::DOUBLE)));
              q->setDescription("BondOr");
            } else if(bType == 6){
              // single or aromatic
              q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(makeBondOrderEqualsQuery(Bond::SINGLE)));
              q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(makeBondOrderEqualsQuery(Bond::AROMATIC)));      
              q->setDescription("BondOr");
            } else if(bType == 7){
              // double or aromatic
              q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(makeBondOrderEqualsQuery(Bond::DOUBLE)));
              q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(makeBondOrderEqualsQuery(Bond::AROMATIC)));
              q->setDescription("BondOr");
            }
            bond->setQuery(q);
          } else {
            BOND_NULL_QUERY *q;
            q = makeBondNullQuery();
            bond->setQuery(q);
            BOOST_LOG(rdWarningLog) << "unrecognized query bond type, " << bType <<", found on line"<<line<<". Using an \"any\" query."<<std::endl;          
          }
          break;
        }

        // additional bond properties:
        unsigned int lPos=4;
        std::ostringstream errout;
        while(lPos<splitLine.size()){
          std::string prop,val;
          if(!splitAssignToken(splitLine[lPos],prop,val)){
            errout << "bad bond property '"<<splitLine[lPos]<<"' on line "<<line;
            throw FileParseException(errout.str()) ;
          }
          if(prop=="CFG"){
            unsigned int cfg=atoi(val.c_str());
            switch(cfg){
            case 0: break;
            case 1:
              bond->setBondDir(Bond::BEGINWEDGE);
	      chiralityPossible=true;
              break;
            case 2:
              if(bType==1) bond->setBondDir(Bond::UNKNOWN);
              else if(bType==2){
		bond->setBondDir(Bond::EITHERDOUBLE);
		bond->setStereo(Bond::STEREOANY);
	      }
              break;
            case 3:
              bond->setBondDir(Bond::BEGINDASH);
	      chiralityPossible=true;
              break;
            default:
              errout << "bad bond CFG "<<val<<"' on line "<<line;
              throw FileParseException(errout.str()) ;
            }
          } else if(prop=="TOPO"){
            if(val!="0"){
              if(!bond->hasQuery()){
                QueryBond *qBond=new QueryBond(*bond);
                delete bond;
                bond=qBond;
              }
              BOND_EQUALS_QUERY *q=makeBondIsInRingQuery();
              if(val=="1"){
                // nothing
              } else if(val=="2"){
                q->setNegation(true);
              } else {
                errout << "bad bond TOPO "<<val<<"' on line "<<line;
                throw FileParseException(errout.str()) ;
              }
              bond->expandQuery(q);          
            }
          } else if(prop=="RXCTR"){
            int reactStatus = FileParserUtils::toInt(val);
            bond->setProp("molReactStatus",reactStatus);
          } else if(prop=="STBOX"){
          }
          ++lPos;
        }

        bond->setBeginAtomIdx(mol->getAtomWithBookmark(a1Idx)->getIdx());
        bond->setEndAtomIdx(mol->getAtomWithBookmark(a2Idx)->getIdx());
        mol->addBond(bond,true);
        if(bond->getIsAromatic()){
          mol->getAtomWithIdx(bond->getBeginAtomIdx())->setIsAromatic(true);
          mol->getAtomWithIdx(bond->getEndAtomIdx())->setIsAromatic(true);
        }
        mol->setBondBookmark(bond,bondIdx);
      }
      tempStr = getV3000Line(inStream,line);
      if(tempStr.length()<8 || tempStr.substr(0,8) != "END BOND"){
        std::ostringstream errout;
        errout<<"END BOND line not found at line "<<line;
        throw FileParseException(errout.str());
      }
    }

    void ProcessMolProps(RWMol *mol){
      PRECONDITION(mol,"no molecule");
      for(RWMol::AtomIterator atomIt=mol->beginAtoms();
          atomIt!=mol->endAtoms();
          ++atomIt) {
        Atom *atom=*atomIt;
        if(atom->hasProp("molTotValence") && !atom->hasProp("_ZBO_H")){
          int totV;
          atom->getProp("molTotValence",totV);
          if(totV==0) continue;
          atom->setNoImplicit(true);
          if(totV==15 // V2000
             || totV==-1 // v3000
             ){
            atom->setNumExplicitHs(0);
          } else {
            if(atom->getExplicitValence()>totV){
              BOOST_LOG(rdWarningLog) << "atom " << atom->getIdx() <<" has specified valence ("<<totV<<") smaller than the drawn valence "<<atom->getExplicitValence()<<"."<<std::endl;
              atom->setNumExplicitHs(0);
            } else {
              atom->setNumExplicitHs(totV-atom->getExplicitValence());
            }
          }
        }
      }
    }
    
  } // end of local namespace
  namespace FileParserUtils {
    bool ParseV3000CTAB(std::istream *inStream,unsigned int &line,
			RWMol *mol, Conformer *&conf,
			bool &chiralityPossible,unsigned int &nAtoms,
			unsigned int &nBonds,
                        bool strictParsing,
                        bool expectMEND){
      PRECONDITION(inStream,"bad stream");
      PRECONDITION(mol,"bad molecule");

      std::string tempStr;
      std::vector<std::string> splitLine;

      bool fileComplete=false;

      tempStr = getV3000Line(inStream,line);
      boost::to_upper(tempStr);
      if(tempStr.length()<10 || tempStr.substr(0,10) != "BEGIN CTAB"){
        std::ostringstream errout;
        errout << "BEGIN CTAB line not found on line "<<line;
        throw FileParseException(errout.str()) ;
      }
      
      tempStr = getV3000Line(inStream,line);
      boost::to_upper(tempStr);
      if(tempStr.size()<8 || tempStr.substr(0,7)!="COUNTS "){
        std::ostringstream errout;
        errout << "Bad counts line : '"<<tempStr<<"' on line "<<line;
        throw FileParseException(errout.str()) ;
      }
      std::string trimmed=boost::trim_copy(tempStr.substr(7,tempStr.length()-7));
      boost::split(splitLine,trimmed,boost::is_any_of(" \t"),boost::token_compress_on);
      if(splitLine.size()<2){
        std::ostringstream errout;
        errout << "Bad counts line : '"<<tempStr<<"' on line "<<line;
        throw FileParseException(errout.str()) ;
      }

      nAtoms=FileParserUtils::toInt(splitLine[0]);
      nBonds=FileParserUtils::toInt(splitLine[1]);
      if(!nAtoms){
        throw FileParseException("molecule has no atoms");
      }
      conf = new Conformer(nAtoms);
      
      unsigned int nSgroups=0,n3DConstraints=0,chiralFlag=0;
      if(splitLine.size()>2) nSgroups = FileParserUtils::toInt(splitLine[2]);
      if(splitLine.size()>3) n3DConstraints = FileParserUtils::toInt(splitLine[3]);
      if(splitLine.size()>4) chiralFlag = FileParserUtils::toInt(splitLine[4]);

      ParseV3000AtomBlock(inStream,line,nAtoms,mol,conf);
      if(nBonds){
	ParseV3000BondBlock(inStream,line,nBonds,mol,chiralityPossible);
      }

      if(nSgroups){
        tempStr = getV3000Line(inStream,line);
        boost::to_upper(tempStr);
        if(tempStr.length()<12 || tempStr.substr(0,12) != "BEGIN SGROUP"){
          std::ostringstream errout;
          errout<<"BEGIN SGROUP line not found on line "<<line;
          throw FileParseException(errout.str());
        }
        for(unsigned int si=0;si<nSgroups;++si){
          tempStr = getV3000Line(inStream,line);
	  boost::to_upper(tempStr);
          std::vector<std::string> localSplitLine;
          boost::split(localSplitLine,tempStr,boost::is_any_of(" \t"),boost::token_compress_on);
          std::string typ = localSplitLine[1];
          if(strictParsing && !SGroupOK(typ)){
            std::ostringstream errout;
            errout << "S group "<<typ<<" on line "<<line;
            throw MolFileUnhandledFeatureException(errout.str()) ;
          } else {
            BOOST_LOG(rdWarningLog) << " S group " << typ <<" ignored on line "<<line<<"."<<std::endl;
          }
        }
        tempStr = getV3000Line(inStream,line);
        boost::to_upper(tempStr);
        if(tempStr.length()<10 || tempStr.substr(0,10) != "END SGROUP"){
          std::ostringstream errout;
          errout<<"END SGROUP line not found on line "<<line;
          throw FileParseException(errout.str()) ;
        }
      }
      if(n3DConstraints){
        BOOST_LOG(rdWarningLog)<<"3d constraint information in mol block igored at line "<<line<<std::endl;
        tempStr = getV3000Line(inStream,line);
	boost::to_upper(tempStr);
        if(tempStr.length()<11 || tempStr.substr(0,11) != "BEGIN OBJ3D"){
          std::ostringstream errout;
          errout << "BEGIN OBJ3D line not found on line "<<line;
          throw FileParseException(errout.str()) ;
        }
        for(unsigned int i=0;i<n3DConstraints;++i) tempStr = getV3000Line(inStream,line);
        tempStr = getV3000Line(inStream,line);
	boost::to_upper(tempStr);
        if(tempStr.length()<9 || tempStr.substr(0,9) != "END OBJ3D"){
          std::ostringstream errout;
          errout << "END OBJ3D line not found on line "<<line;
          throw FileParseException(errout.str()) ;
        }
      }
      
      tempStr = getV3000Line(inStream,line);
      // do link nodes:
      boost::to_upper(tempStr);
      while(tempStr.length()>8 && tempStr.substr(0,8)=="LINKNODE"){
        tempStr = getV3000Line(inStream,line);
	boost::to_upper(tempStr);
      }

      while(tempStr.length()>5 && tempStr.substr(0,5)=="BEGIN"){
        // skip blocks we don't know how to read
        BOOST_LOG(rdWarningLog)<<"skipping block at line "<<line<<": "<<tempStr<<std::endl;
        tempStr = getV3000Line(inStream,line);
        
        while(tempStr.length()<3 || tempStr.substr(0,3)!="END"){
          tempStr = getV3000Line(inStream,line);
        }
        tempStr = getV3000Line(inStream,line);
      }

      boost::to_upper(tempStr);
      if(tempStr.length()<8 || tempStr.substr(0,8) != "END CTAB"){
        throw FileParseException("END CTAB line not found") ;
      }

      if(expectMEND){
        tempStr = getLine(inStream);
        ++line;
        if(tempStr[0]=='M'&&tempStr.substr(0,6)=="M  END"){
          fileComplete=true;
        }
      } else {
        fileComplete = true;
      }

      mol->addConformer(conf, true);
      conf=0;

      return fileComplete;
    }

    bool ParseV2000CTAB(std::istream *inStream,unsigned int &line,
			RWMol *mol, Conformer *&conf,
			bool &chiralityPossible,unsigned int &nAtoms,
			unsigned int &nBonds,bool strictParsing){
        conf = new Conformer(nAtoms);
        if(nAtoms==0){
          conf->set3D(false);
        } else {
          ParseMolBlockAtoms(inStream,line,nAtoms,mol,conf);

          if(mol->hasProp("_2DConf")){
            conf->set3D(false);
            mol->clearProp("_2DConf");
          } else if(mol->hasProp("_3DConf")){
            conf->set3D(true);
            mol->clearProp("_3DConf");
          }
        }
        mol->addConformer(conf, true);
        conf=0;
        
        ParseMolBlockBonds(inStream,line,nBonds,mol,chiralityPossible);
      
        bool fileComplete=ParseMolBlockProperties(inStream,line,mol);
	return fileComplete;
    }

  }  // end of FileParserUtils namespace

  //------------------------------------------------
  //
  //  Read a molecule from a stream
  //
  //------------------------------------------------
  RWMol *MolDataStreamToMol(std::istream *inStream, unsigned int &line, bool sanitize,
                            bool removeHs,bool strictParsing){
    PRECONDITION(inStream,"no stream");
    std::string tempStr;
    bool fileComplete=false;
    bool chiralityPossible = false;
    Utils::LocaleSwitcher ls;
    // mol name
    line++;
    tempStr = getLine(inStream);
    if(inStream->eof()){
      return NULL;
    }
    RWMol *res = new RWMol();
    res->setProp("_Name", tempStr);

    // info
    line++;
    tempStr = getLine(inStream);
    res->setProp("_MolFileInfo", tempStr);
    if(tempStr.length()>=22){
      std::string dimLabel=tempStr.substr(20,2);
      if(dimLabel=="2d"||dimLabel=="2D"){
        res->setProp("_2DConf",1);
      } else if(dimLabel=="3d"||dimLabel=="3D"){
        res->setProp("_3DConf",1);
      }
    }
    // comments
    line++;
    tempStr = getLine(inStream);
    res->setProp("_MolFileComments", tempStr);
        
    unsigned int nAtoms=0,nBonds=0,nLists=0,chiralFlag=0,nsText=0,nRxnComponents=0;
    int nReactants=0,nProducts=0,nIntermediates=0;
    // counts line, this is where we really get started
    line++;
    tempStr = getLine(inStream);

    if(tempStr.size()<6){
      if(res){
        delete res;
        res = NULL;
      }
      std::ostringstream errout;
      errout << "Counts line too short: '"<<tempStr<<"' on line"<<line;
      throw FileParseException(errout.str()) ;
    }

    unsigned int spos = 0;
    // this needs to go into a try block because if the lexical_cast throws an
    // exception we want to catch and delete mol before leaving this function
    try {
      nAtoms = FileParserUtils::toInt(tempStr.substr(spos,3));
      spos = 3;
      nBonds = FileParserUtils::toInt(tempStr.substr(spos,3));
      spos = 6;
    } catch (boost::bad_lexical_cast &) {
      if(res){
        delete res;
        res = NULL;
      }
      std::ostringstream errout;
      errout << "Cannot convert " << tempStr.substr(spos,3) << " to int on line "<<line;
      throw FileParseException(errout.str()) ;
    }
    try {
      spos = 6;
      if(tempStr.size()>=9)
        nLists = FileParserUtils::toInt(tempStr.substr(spos,3));

      spos = 12;
      if(tempStr.size()>=spos+3)
        chiralFlag = FileParserUtils::toInt(tempStr.substr(spos,3));

      spos = 15;
      if(tempStr.size()>=spos+3)
        nsText = FileParserUtils::toInt(tempStr.substr(spos,3));

      spos = 18;
      if(tempStr.size()>=spos+3)
        nRxnComponents = FileParserUtils::toInt(tempStr.substr(spos,3));

      spos = 21;
      if(tempStr.size()>=spos+3)
        nReactants   = FileParserUtils::toInt(tempStr.substr(spos,3));

      spos = 24;
      if(tempStr.size()>=spos+3)
        nProducts   = FileParserUtils::toInt(tempStr.substr(spos,3));

      spos = 27;
      if(tempStr.size()>=spos+3)
        nIntermediates = FileParserUtils::toInt(tempStr.substr(spos,3));

    } catch (boost::bad_lexical_cast &) {
      // some SD files (such as some from NCI) lack all the extra information
      // on the header line, so ignore problems parsing there.
    }

    unsigned int ctabVersion=2000;
    if(tempStr.size()>35){
      if(tempStr.size()<39 || tempStr[34]!='V' ){
        std::ostringstream errout;
        errout<<"CTAB version string invalid at line "<<line;
        if(strictParsing){
          delete res;
          res=NULL;
          throw FileParseException(errout.str());
        } else {
          BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
        }
      } else if(tempStr.substr(34,5)=="V3000"){
        ctabVersion=3000;
      } else if(tempStr.substr(34,5)!="V2000"){
        std::ostringstream errout;
        errout << "Unsupported CTAB version: '"<< tempStr.substr(34,5) << "' at line " << line;
        if(strictParsing){
          delete res;
          res = NULL;
          throw FileParseException(errout.str()) ;
        } else {
          BOOST_LOG(rdWarningLog) << errout.str() <<std::endl;
        }
      }
    }

    if(chiralFlag){
      res->setProp("_MolFileChiralFlag",chiralFlag);
    }

    Conformer *conf=0;
    try {
      if(ctabVersion==2000){
        fileComplete=FileParserUtils::ParseV2000CTAB(inStream,line,
						     res,conf,chiralityPossible,
						     nAtoms,nBonds,strictParsing);
      } else {
        if(nAtoms!=0 || nBonds!=0){
          std::ostringstream errout;
          errout << "V3000 mol blocks should have 0s in the initial counts line. (line: "<<line<<")";
          if(strictParsing){
            delete res;
            res = NULL;
            throw FileParseException(errout.str()) ;
          } else {
            BOOST_LOG(rdWarningLog)<<errout.str()<<std::endl ;
          }
        }
        fileComplete=FileParserUtils::ParseV3000CTAB(inStream,line,
						     res,conf,chiralityPossible,
						     nAtoms,nBonds,strictParsing);
      }
    } catch (MolFileUnhandledFeatureException &e) { 
      // unhandled mol file feature, just delete the result 
      delete res;
      delete conf;
      res=NULL;
      conf=NULL;
      BOOST_LOG(rdErrorLog) << " Unhandled CTAB feature: " << e.message() <<" on line: "<<line<<". Molecule skipped."<<std::endl;

      if(!inStream->eof()) tempStr = getLine(inStream);
      ++line;
      while(!inStream->eof() && tempStr.substr(0,6)!="M  END" && tempStr.substr(0,4)!="$$$$"){
        tempStr = getLine(inStream);
        ++line;
      }
      if(!inStream->eof() || tempStr.substr(0,6)=="M  END" || tempStr.substr(0,4)=="$$$$")
      fileComplete=true;
      else
        fileComplete=false;        
    } catch (FileParseException &e) { 
      // catch our exceptions and throw them back after cleanup
      delete res;
      delete conf;
      res=NULL;
      conf=NULL;
      throw e;
    }

    if(!fileComplete){
      delete res;
      delete conf;
      res=NULL;
      conf=NULL;
      std::ostringstream errout;
      errout << "Problems encountered parsing Mol data, M  END missing around line "<< line;
      throw FileParseException(errout.str());
    }

    if (res ) {
      // calculate explicit valence on each atom:
      for(RWMol::AtomIterator atomIt=res->beginAtoms();
          atomIt!=res->endAtoms();
          ++atomIt) {
        (*atomIt)->calcExplicitValence(false);
      }

      // postprocess mol file flags
      ProcessMolProps(res);

      // update the chirality and stereo-chemistry
      //
      // NOTE: we detect the stereochemistry before sanitizing/removing
      // hydrogens because the removal of H atoms may actually remove
      // the wedged bond from the molecule.  This wipes out the only
      // sign that chirality ever existed and makes us sad... so first
      // perceive chirality, then remove the Hs and sanitize.
      //
      // One exception to this (of course, there's always an exception):
      // DetectAtomStereoChemistry() needs to check the number of
      // implicit hydrogens on atoms to detect if things can be
      // chiral. However, if we ask for the number of implicit Hs before
      // we've called MolOps::cleanUp() on the molecule, we'll get
      // exceptions for common "weird" cases like a nitro group
      // mis-represented as -N(=O)=O.  *SO*... we need to call
      // cleanUp(), then detect the stereochemistry.
      // (this was Issue 148)
      //
      const Conformer &conf = res->getConformer();
      if(chiralityPossible){
        MolOps::cleanUp(*res);
        DetectAtomStereoChemistry(*res, &conf);
      }

      if ( sanitize ) {
        try {
          if(removeHs){
            MolOps::removeHs(*res,false,false);
          } else {
            MolOps::sanitizeMol(*res);
          }

          // now that atom stereochem has been perceived, the wedging
          // information is no longer needed, so we clear
          // single bond dir flags:
          ClearSingleBondDirFlags(*res);
      
          // unlike DetectAtomStereoChemistry we call DetectBondStereoChemistry 
          // here after sanitization because we need the ring information:
          DetectBondStereoChemistry(*res, &conf);
        }
        catch (...){
          delete res;
          res=NULL;
          throw;
        }
        MolOps::assignStereochemistry(*res,true,true,true);
      } else {
        // we still need to do something about double bond stereochemistry
        // (was github issue 337)
        DetectBondStereoChemistry(*res, &conf);
      }

      if(res->hasProp("_NeedsQueryScan")){
        res->clearProp("_NeedsQueryScan");
        CompleteMolQueries(res);
      }
    }
    return res;
  };
  

  RWMol *MolDataStreamToMol(std::istream &inStream, unsigned int &line,
                            bool sanitize, bool removeHs,bool strictParsing){
    return MolDataStreamToMol(&inStream,line,sanitize,removeHs,strictParsing);
  };
  //------------------------------------------------
  //
  //  Read a molecule from a string
  //
  //------------------------------------------------
  RWMol *MolBlockToMol(const std::string &molBlock, bool sanitize, bool removeHs,bool strictParsing){
    std::istringstream inStream(molBlock);
    unsigned int line = 0;
    return MolDataStreamToMol(inStream, line, sanitize, removeHs,strictParsing);
  }    


  //------------------------------------------------
  //
  //  Read a molecule from a file
  //
  //------------------------------------------------
  RWMol *MolFileToMol(std::string fName, bool sanitize, bool removeHs,bool strictParsing){
    std::ifstream inStream(fName.c_str());
    if (!inStream || (inStream.bad()) ) {
      std::ostringstream errout;
      errout << "Bad input file " << fName;
      throw BadFileException(errout.str());
    }
    RWMol *res=NULL;
    if(!inStream.eof()){
      unsigned int line = 0;
      res=MolDataStreamToMol(inStream, line, sanitize, removeHs,strictParsing);
    }
    return res;
  }    
}

