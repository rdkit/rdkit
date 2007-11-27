// $Id$
//
//  Copyright (C) 2002-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#include "FileParsers.h"
#include "MolFileStereochem.h"
#include <GraphMol/RDKitQueries.h>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/RDLog.h>

#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <RDGeneral/FileParseException.h>
#include <typeinfo>

namespace RDKit{
  
  
  // it's kind of stinky that we have to do this, but as of g++3.2 and
  // boost 1.30, on linux calls to lexical_cast<int>(std::string)
  // crash if the string starts with spaces.
  template <typename T>
  T stripSpacesAndCast(const std::string &input,bool acceptSpaces=false){
    std::string trimmed=boost::trim_copy(input);
    if(acceptSpaces && trimmed==""){
      return 0;
    } else {
      return boost::lexical_cast<T>(trimmed);
    }
  }

  //*************************************
  //
  // Every effort has been made to adhere to MDL's standard
  // for mol files
  //  
  //*************************************
  
  void ParseOldAtomList(RWMol *mol,std::string text){
    unsigned int idx;
    try {
      idx = stripSpacesAndCast<unsigned int>(text.substr(0,3))-1;
    }
    catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert " << text.substr(0,3) << " to int";
      throw FileParseException(errout.str()) ;
    }

    RANGE_CHECK(0,idx,mol->getNumAtoms()-1);
    QueryAtom a;
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
      errout << "Unrecognized atom-list query modifier: " << text[14];
      throw FileParseException(errout.str()) ;
    }          
    
    int nQueries;
    try {
      nQueries = stripSpacesAndCast<int>(text.substr(9,1));
    }
    catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert " << text.substr(9,1) << " to int";
      throw FileParseException(errout.str()) ;
    }

    RANGE_CHECK(0,nQueries,5);
    for(int i=0;i<nQueries;i++){
      int pos = 11+i*4;
      int atNum;
      try {
        atNum = stripSpacesAndCast<int>(text.substr(pos,3));
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(pos,3) << " to int";
        throw FileParseException(errout.str()) ;
      }
      RANGE_CHECK(0,atNum,200);  // goofy!
      q->addChild(QueryAtom::QUERYATOM_QUERY::CHILD_TYPE(makeAtomNumEqualsQuery(atNum)));
    }
    
    a.setQuery(q);
    mol->replaceAtom(idx,&a); 
  };
  
  void ParseChargeLine(RWMol *mol, std::string text,bool firstCall) {
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
      nent = stripSpacesAndCast<int>(text.substr(6,3));
    }
    catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert " << text.substr(6,3) << " to int";
      throw FileParseException(errout.str()) ;
    }
    int spos = 9;
    for (ie = 0; ie < nent; ie++) {
      int aid, chg;
      try {
        aid = stripSpacesAndCast<int>(text.substr(spos,4));
        spos += 4;
        chg = stripSpacesAndCast<int>(text.substr(spos,4));
        spos += 4;
        mol->getAtomWithIdx(aid-1)->setFormalCharge(chg);
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(spos,4) << " to int";
        throw FileParseException(errout.str()) ;
      }
    }
  }

  void ParseIsotopeLine(RWMol *mol, std::string text){
    PRECONDITION(text.substr(0,6)==std::string("M  ISO"),"bad isotope line");
    
    unsigned int nent;
    try {
      nent = stripSpacesAndCast<unsigned int>(text.substr(6,3));
    }
    catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert " << text.substr(6,3) << " to int";
      throw FileParseException(errout.str()) ;
    }
    unsigned int spos = 9;
    for (unsigned int ie = 0; ie < nent; ie++) {
      unsigned int aid;
      int mass;
      try {
        aid = stripSpacesAndCast<unsigned int>(text.substr(spos,4));
        spos += 4;
        Atom *atom=mol->getAtomWithIdx(aid-1); 
        if(text.size()>=spos+4 && text.substr(spos,4)!="    "){
          mass = stripSpacesAndCast<int>(text.substr(spos,4));
          atom->setMass(static_cast<double>(mass));
          spos += 4;
        } else {
          atom->setMass(PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum()));
        }
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(spos,4) << " to int";
        throw FileParseException(errout.str()) ;
      }
    }
    
  }

  void ParseSubstitutionCountLine(RWMol *mol, std::string text){
    PRECONDITION(text.substr(0,6)==std::string("M  SUB"),"bad SUB line");
    
    unsigned int nent;
    try {
      nent = stripSpacesAndCast<unsigned int>(text.substr(6,3));
    }
    catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert " << text.substr(6,3) << " to int";
      throw FileParseException(errout.str()) ;
    }
    unsigned int spos = 9;
    for (unsigned int ie = 0; ie < nent; ie++) {
      unsigned int aid;
      int count;
      try {
        aid = stripSpacesAndCast<unsigned int>(text.substr(spos,4));
        spos += 4;
        Atom *atom=mol->getAtomWithIdx(aid-1); 
        if(text.size()>=spos+4 && text.substr(spos,4)!="    "){
          count = stripSpacesAndCast<int>(text.substr(spos,4));
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
              BOOST_LOG(rdWarningLog) << " atom degree query with value 6 found. This will not match degree >6. The MDL spec says it should.";
              q->setVal(6);break;
            default:
              std::ostringstream errout;
              errout << "Value " << count << " is not supported as a degree query.";
              throw FileParseException(errout.str()) ;
          }
          if(!atom->hasQuery()){
            QueryAtom a(atom->getAtomicNum());
            mol->replaceAtom(aid-1,&a);           
            atom = mol->getAtomWithIdx(aid-1);
          }
          atom->expandQuery(q,Queries::COMPOSITE_AND);
          spos += 4;
        }
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(spos,4) << " to int";
        throw FileParseException(errout.str()) ;
      }
    }
  }

  void ParseNewAtomList(RWMol *mol,std::string text){
    PRECONDITION(text.substr(0,6)==std::string("M  ALS"),"bad atom list line");
    
    unsigned int idx;
    try {
      idx = stripSpacesAndCast<unsigned int>(text.substr(9,3))-1;
    }
    catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert " << text.substr(9,3) << " to int";
      throw FileParseException(errout.str()) ;
    }
    RANGE_CHECK(0,idx,mol->getNumAtoms()-1);
    QueryAtom a;
    ATOM_OR_QUERY *q = new ATOM_OR_QUERY;
    q->setDescription("AtomOr");
    
    switch(text[14]){
    case 'T':
      q->setNegation(true);
      break;
    case 'F':
      q->setNegation(false);
      break;
    default:
      std::ostringstream errout;
      errout << "Unrecognized atom-list query modifier: " << text[14];
      throw FileParseException(errout.str()) ;
    }          
    
    int nQueries;
    try {
      nQueries = stripSpacesAndCast<int>(text.substr(10,3));
    }
    catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert " << text.substr(10,3) << " to int";
      throw FileParseException(errout.str()) ;
    }

    for(int i=0;i<nQueries;i++){
      int pos = 16+i*4;
      std::string atSymb = text.substr(pos,4);
      atSymb.erase(atSymb.find(" "),atSymb.size());
      int atNum = PeriodicTable::getTable()->getAtomicNumber(atSymb);
      q->addChild(QueryAtom::QUERYATOM_QUERY::CHILD_TYPE(makeAtomNumEqualsQuery(atNum)));
    }
    a.setQuery(q);
    mol->replaceAtom(idx,&a); 
  };
  
  void ParseRGroupLabels(RWMol *mol,std::string text){
    PRECONDITION(text.substr(0,6)==std::string("M  RGP"),"bad R group label line");
    
    int nLabels;
    try {
      nLabels = stripSpacesAndCast<int>(text.substr(6,3));
    }
    catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert " << text.substr(6,3) << " to int";
      throw FileParseException(errout.str()) ;
    }

    for(int i=0;i<nLabels;i++){
      int pos = 10+i*8;
      unsigned int atIdx;
      try {
        atIdx = stripSpacesAndCast<unsigned int>(text.substr(pos,3));
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(pos,3) << " to int";
        throw FileParseException(errout.str()) ;
      }
      unsigned int rLabel;
      try {
        rLabel = stripSpacesAndCast<unsigned int>(text.substr(pos+4,3));
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(pos+4,3) << " to int";
        throw FileParseException(errout.str()) ;
      }
      atIdx-=1;
      if(atIdx>mol->getNumAtoms()){
        std::ostringstream errout;
        errout << "Attempt to set R group label on nonexistent atom " << atIdx;
        throw FileParseException(errout.str()) ;
      }
      QueryAtom atom;
      atom.setProp("_MolFileRLabel",rLabel);
      std::string tmpLabel="X";
      if(rLabel>0 && rLabel<10){
        switch(rLabel){
        case 1:
          tmpLabel="Xa";
          break;
        case 2:
          tmpLabel="Xb";
          break;
        case 3:
          tmpLabel="Xc";
          break;
        case 4:
          tmpLabel="Xd";
          break;
        case 5:
          tmpLabel="Xf";
          break;
        case 6:
          tmpLabel="Xg";
          break;
        case 7:
          tmpLabel="Xh";
          break;
        case 8:
          tmpLabel="Xi";
          break;
        case 9:
          tmpLabel="Xj";
          break;
        }
      }
      atom.setProp("dummyLabel",tmpLabel.c_str());
      atom.setQuery(makeAtomNullQuery());
      mol->replaceAtom(atIdx,&atom); 
    }
  };
  
  void ParseAtomAlias(RWMol *mol,std::string text,std::string &nextLine){
    PRECONDITION(text.substr(0,2)==std::string("A "),"bad atom alias line");
      
    unsigned int idx;
    try {
      idx = stripSpacesAndCast<unsigned int>(text.substr(3,3))-1;
    }
    catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert " << text.substr(3,3) << " to int";
      throw FileParseException(errout.str()) ;
    }
    RANGE_CHECK(0,idx,mol->getNumAtoms()-1);
    Atom *at = mol->getAtomWithIdx(idx);
    at->setProp("molFileAlias",nextLine);
  };
  
  Atom *ParseMolFileAtomLine(const std::string text, RDGeom::Point3D &pos) {
    Atom *res = new Atom;
    //double pX,pY,pZ;
    std::string symb;
    int massDiff,chg,hCount;
    int stereoCare,totValence;
    int atomMapNumber,inversionFlag,exactChangeFlag;
    
    try {
      pos.x = stripSpacesAndCast<double>(text.substr(0,10));
      pos.y = stripSpacesAndCast<double>(text.substr(10,10));
      pos.z = stripSpacesAndCast<double>(text.substr(20,10));
    }
    catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot process coordinates.";
      throw FileParseException(errout.str()) ;
    }
    symb = text.substr(31,3);
    symb = symb.substr(0,symb.find(' '));
    
    // REVIEW: should we handle missing fields at the end of the line?
    massDiff=0;
    if(text.size()>=36 && text.substr(34,2)!=" 0"){
      try {
        massDiff = stripSpacesAndCast<int>(text.substr(34,2),true);
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(34,2) << " to int";
        throw FileParseException(errout.str()) ;
      }
    }    
    chg=0;
    if(text.size()>=39 && text.substr(36,3)!="  0"){
      try {
        chg = stripSpacesAndCast<int>(text.substr(36,3),true);
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(36,3) << " to int";
        throw FileParseException(errout.str()) ;
      }
    }
    hCount = 0;
    if(text.size()>=45 && text.substr(42,3)!="  0"){
      try {
        hCount = stripSpacesAndCast<int>(text.substr(42,3),true);
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(42,3) << " to int";
        throw FileParseException(errout.str()) ;
      }
    }
    if(symb=="L" || symb=="A" || symb=="Q" || symb=="*" || symb=="LP"
       || symb=="R" || symb=="R#" || (symb>="R0" && symb<="R9") ){

      if(symb=="A"||symb=="Q"||symb=="*"){
        // according to the MDL spec, these match anything 
        QueryAtom *query=new QueryAtom(0);
        query->setQuery(makeAtomNullQuery());
        delete res;
        res=query;  
      } else {
        res->setAtomicNum(0);
      }
      if(symb=="*") res->setProp("dummyLabel",std::string("X"));
      else if(symb=="R0") res->setProp("dummyLabel",std::string("Xa"));
      else if(symb=="R1") res->setProp("dummyLabel",std::string("Xb"));
      else if(symb=="R2") res->setProp("dummyLabel",std::string("Xc"));
      else if(symb=="R3") res->setProp("dummyLabel",std::string("Xd"));
      else if(symb=="R4") res->setProp("dummyLabel",std::string("Xf"));
      else if(symb=="R5") res->setProp("dummyLabel",std::string("Xg"));
      else if(symb=="R6") res->setProp("dummyLabel",std::string("Xh"));
      else if(symb=="R7") res->setProp("dummyLabel",std::string("Xi"));
      else if(symb=="R8") res->setProp("dummyLabel",std::string("Xj"));
      else if(symb=="R9") res->setProp("dummyLabel",std::string("Xk"));
      else if(symb=="R#") res->setProp("dummyLabel",std::string("X"));
      else res->setProp("dummyLabel",symb);
      
      
    } else if( symb=="D" ){  // mol blocks support "D" and "T" as shorthand... handle that.
      res->setAtomicNum(1); 
      res->setMass(2.014);
    } else if( symb=="T" ){  // mol blocks support "D" and "T" as shorthand... handle that.
      res->setAtomicNum(1);
      res->setMass(3.016);
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
      // FIX: this isn't precisely correct because we should be doing the difference w.r.t. most abundant species.
      res->setMass(res->getMass()+massDiff);
    }
    
#if 1
    stereoCare=0;
    if(text.size()>=48 && text.substr(45,3)!="  0"){
      try {
        stereoCare = stripSpacesAndCast<int>(text.substr(45,3),true);
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(45,3) << " to int";
        throw FileParseException(errout.str()) ;
      }
      res->setProp("molStereoCare",stereoCare);
    }
    totValence=0;
    if(text.size()>=51 && text.substr(48,3)!="  0"){
      try {
        totValence= stripSpacesAndCast<int>(text.substr(48,3),true);
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(48,3) << " to int";
        throw FileParseException(errout.str()) ;
      }
      res->setProp("molTotValence",totValence);
    }
    atomMapNumber=0;
    if(text.size()>=63 && text.substr(60,3)!="  0"){
      try {
        atomMapNumber = stripSpacesAndCast<int>(text.substr(60,3),true);
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(60,3) << " to int";
        throw FileParseException(errout.str()) ;
      }
      res->setProp("molAtomMapNumber",atomMapNumber);
    }
    inversionFlag=0;
    if(text.size()>=66 && text.substr(63,3)!="  0"){
      try {
        inversionFlag= stripSpacesAndCast<int>(text.substr(63,3),true);
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(63,3) << " to int";
        throw FileParseException(errout.str()) ;
      }
      res->setProp("molInversionFlag",inversionFlag);
    }
    exactChangeFlag=0;
    if(text.size()>=69 && text.substr(66,3)!="  0"){
      try {
        exactChangeFlag = stripSpacesAndCast<int>(text.substr(66,3),true);
      }
      catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Cannot convert " << text.substr(66,3) << " to int";
        throw FileParseException(errout.str()) ;
      }
      res->setProp("molExactChangeFlag",exactChangeFlag);
    }
#endif    
    return res;
  };
  
  Bond *ParseMolFileBondLine(const std::string text){
    int idx1,idx2,bType,stereo;
    int spos = 0;
    try {
      idx1 = stripSpacesAndCast<int>(text.substr(spos,3));
      spos += 3;
      idx2 = stripSpacesAndCast<int>(text.substr(spos,3));
      spos += 3;
      bType = stripSpacesAndCast<int>(text.substr(spos,3));  
    }
    catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert " << text.substr(spos,3) << " to int";
      throw FileParseException(errout.str()) ;
    }
    
    
    // adjust the numbering
    idx1--;idx2--;

    Bond::BondType type;
    Bond *res;  
    switch(bType){
    case 1: type = Bond::SINGLE;res = new Bond;break;
    case 2: type = Bond::DOUBLE;res = new Bond;break;
    case 3: type = Bond::TRIPLE;res = new Bond;break;
    case 4: type = Bond::AROMATIC;res = new Bond;break;
    default:
      type = Bond::UNSPECIFIED;
      // it's a query bond of some type
      res = new QueryBond;
      if(bType == 8){
        BOND_NULL_QUERY *q;
        q = makeBondNullQuery();
        res->setQuery(q);
      } else {
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
      } 
      break;
    }
    res->setBeginAtomIdx(idx1);
    res->setEndAtomIdx(idx2);
    res->setBondType(type);

    if( text.size() >= 12 && text.substr(9,3)!="  0")
      try {
        stereo = stripSpacesAndCast<int>(text.substr(9,3));
        //res->setProp("stereo",stereo);
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
          break;
        case 4: // "either" single bond
          res->setBondDir(Bond::UNKNOWN);
          break;
        }
      } catch (boost::bad_lexical_cast) {
        ;
      }
#if 1
    if( text.size() >= 18 && text.substr(15,3)!="  0")
      try {
        int topology = stripSpacesAndCast<int>(text.substr(15,3));
        QueryBond *qBond=new QueryBond(*res);
        BOND_EQUALS_QUERY *q=makeBondIsInRingQuery();
        switch(topology){
        case 1:
          break;
        case 2:
          q->setNegation(true);
          break;
        default:
          std::ostringstream errout;
          errout << "Unrecognized bond topology specifier: " << topology;
          throw FileParseException(errout.str()) ;
        }
        qBond->expandQuery(q);          
        delete res;
        res = qBond;
      } catch (boost::bad_lexical_cast) {
        ;
      }
    
    if( text.size() >= 21 && text.substr(18,3)!="  0")
      try {
        int reactStatus = stripSpacesAndCast<int>(text.substr(18,3));
        res->setProp("molReactStatus",reactStatus);
      } catch (boost::bad_lexical_cast) {
        ;
      }
#endif    
    return res;
  };  
  
  //------------------------------------------------
  //
  //  Read a molecule from a stream
  //
  //------------------------------------------------
  RWMol *MolDataStreamToMol(std::istream *inStream, unsigned int &line, bool sanitize,
                            bool removeHs){
    PRECONDITION(inStream,"no stream");
    std::string tempStr;
    bool fileComplete=false;
    bool chiralityPossible = false;

    // mol name
    line++;
    tempStr = getLine(inStream);
    if(inStream->eof()){
      return NULL;
    }
    //std::cerr <<"\tnew" <<std::endl;
    RWMol *res = new RWMol();
    res->setProp("_Name", tempStr);
    bool firstChargeLine=true;

    // info
    line++;
    tempStr = getLine(inStream);
    res->setProp("_MolFileInfo", tempStr);
    // comments
    line++;
    tempStr = getLine(inStream);
    res->setProp("_MolFileComments", tempStr);

    // FIX: if name is NULL, copy in _MolFileInfo or something
        
    int nAtoms=0,nBonds=0,nLists=0,chiralFlag=0,nsText=0,nRxnComponents=0;
    int nReactants=0,nProducts=0,nIntermediates=0;
    // counts line, this is where we really get started
    line++;
    tempStr = getLine(inStream);
    
    // this needs to go into a try block because if the lexical_cast throws an
    // exception we want to catch and delete mol before leaving this function
    
    //std::cerr <<"\tcounts" <<std::endl;
    unsigned int spos = 0;
    try {
      // it *sucks* that the lexical_cast stuff above doesn't work on linux        
      nAtoms = stripSpacesAndCast<int>(tempStr.substr(0,3));
      spos = 3;
      nBonds = stripSpacesAndCast<int>(tempStr.substr(3,3));
      spos = 6;
    } catch (boost::bad_lexical_cast &) {
      if (res) {
        delete res;
      }
      
      std::ostringstream errout;
      errout << "Cannot convert " << tempStr.substr(spos,3) << " to int";
      throw FileParseException(errout.str()) ;
    }
    try {
      spos = 6;
      if(tempStr.size()>=9)
        nLists = stripSpacesAndCast<int>(tempStr.substr(spos,3));

      spos = 12;
      if(tempStr.size()>=spos+3)
        chiralFlag = stripSpacesAndCast<int>(tempStr.substr(spos,3));

      spos = 15;
      if(tempStr.size()>=spos+3)
        nsText = stripSpacesAndCast<int>(tempStr.substr(spos,3));

      spos = 18;
      if(tempStr.size()>=spos+3)
        nRxnComponents = stripSpacesAndCast<int>(tempStr.substr(spos,3));

      spos = 21;
      if(tempStr.size()>=spos+3)
        nReactants   = stripSpacesAndCast<int>(tempStr.substr(spos,3));

      spos = 24;
      if(tempStr.size()>=spos+3)
        nProducts   = stripSpacesAndCast<int>(tempStr.substr(spos,3));

      spos = 27;
      if(tempStr.size()>=spos+3)
        nIntermediates = stripSpacesAndCast<int>(tempStr.substr(spos,3));

    } catch (boost::bad_lexical_cast &) {
      // some SD files (such as some from NCI) lack all the extra information
      // on the header line, so ignore problems parsing there.
    }

    //std::cerr <<"\tnAtoms: " << nAtoms << " " << nBonds << std::endl;
    try {
      if(nAtoms<=0){
        throw FileParseException("molecule has no atoms");
      }
      int i;
      Conformer *conf = new Conformer(nAtoms);
      //std::cerr <<"\trAts" << std::endl;
      for(i=0;i<nAtoms;i++){
        line++;
        tempStr = getLine(inStream);
        CHECK_INVARIANT(!inStream->eof(),"premature EOF");
        RDGeom::Point3D pos;
        Atom *atom = ParseMolFileAtomLine(tempStr, pos);
        unsigned int aid = res->addAtom(atom,false,true);
        conf->setAtomPos(aid, pos);
      }
      //std::cerr <<"\taddConf" << std::endl;
      res->addConformer(conf, true);

      //std::cerr <<"\tbonds" << std::endl;
      for(i=0;i<nBonds;i++){
        line++;
        tempStr = getLine(inStream);
        CHECK_INVARIANT(!inStream->eof(),"premature EOF");
        Bond *bond = ParseMolFileBondLine(tempStr);
        // if we got an aromatic bond set the flag on the bond and the connected atoms
        if (bond->getBondType() == Bond::AROMATIC) {
          bond->setIsAromatic(true);
          res->getAtomWithIdx(bond->getBeginAtomIdx())->setIsAromatic(true);
          res->getAtomWithIdx(bond->getEndAtomIdx())->setIsAromatic(true);
        }
        // if the bond might have chirality info associated with it, set a flag:
        if(bond->getBondDir() != Bond::NONE && bond->getBondDir() != Bond::UNKNOWN){
          chiralityPossible=true;
        }
        res->addBond(bond,true);
      }
      //std::cerr <<"\tok" << std::endl;
      
      // older mol files can have an atom list block here
      line++;
      tempStr = getLine(inStream);
      if( tempStr[0] != 'M' && tempStr[0] != 'A'
          && tempStr[0] != 'V' && tempStr[0] != 'G'){
        ParseOldAtomList(res,tempStr);
      }
      
      while(!inStream->eof() && tempStr[0] != 'M'){
        if(tempStr[0]=='A'){
          line++;
          std::string nextLine = getLine(inStream);
          if(tempStr.substr(0,6)!="M  END"){
            ParseAtomAlias(res,tempStr,nextLine);
          }
        }
        line++;
        tempStr = getLine(inStream);
      }

      //tempStr = inLine;
      std::string lineBeg=tempStr.substr(0,6);
      while(!inStream->eof() && lineBeg!="M  END" && tempStr.substr(0,4)!="$$$$"){
        if(lineBeg=="M  ALS") ParseNewAtomList(res,tempStr);
        else if(lineBeg=="M  ISO") ParseIsotopeLine(res,tempStr);
        else if(lineBeg=="M  RGP") ParseRGroupLabels(res,tempStr);
        else if(lineBeg=="M  SUB") ParseSubstitutionCountLine(res,tempStr);
        else if(lineBeg=="M  CHG") {
          ParseChargeLine(res, tempStr,firstChargeLine);
          firstChargeLine=false;
        }
        line++;
        tempStr = getLine(inStream);
        lineBeg=tempStr.substr(0,6);
      }
      if(tempStr[0]=='M'&&tempStr.substr(0,6)=="M  END"){
        fileComplete=true;
      }
    }
    catch (FileParseException &e) { // catch any exception because of lexical_casting etc
      // and throw them back after cleanup
      if (res) {
        delete res;
      }
      throw e;
    }

    if(!fileComplete){
      delete res;
      throw FileParseException("Problems encountered parsing Mol data, M  END ");
    }


    // calculate explicit valence on each atom:
    for(RWMol::AtomIterator atomIt=res->beginAtoms();
        atomIt!=res->endAtoms();
        atomIt++) {
      (*atomIt)->calcExplicitValence();
    }
    //std::cerr << "bloop "<<line << std::endl;

    if (res && sanitize ) {
      // update the chirality and stereo-chemistry and stuff:
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
      if(chiralityPossible){
        MolOps::cleanUp(*res);
        const Conformer &conf = res->getConformer();
        DetectAtomStereoChemistry(*res, &conf);
      }

      try {
        if(removeHs){
          ROMol *tmp=MolOps::removeHs(*res,false,false);
          delete res;
          res = static_cast<RWMol *>(tmp);
        } else {
          MolOps::sanitizeMol(*res);
        }
        // unlike DetectAtomStereoChemistry we call DetectBondStereoChemistry here after
        // sanitization because the rings should have been perceived by now, in order to
        // correctly recognize double bonds that may be cis/trans type
        const Conformer &conf = res->getConformer();
        DetectBondStereoChemistry(*res, &conf);
      }
      catch (MolSanitizeException &se){
        delete res;
        throw se;
      }
    }

    
    return res;
  };
  

  RWMol *MolDataStreamToMol(std::istream &inStream, unsigned int &line,
                            bool sanitize, bool removeHs){
    return MolDataStreamToMol(&inStream,line,sanitize,removeHs);
  };
  //------------------------------------------------
  //
  //  Read a molecule from a string
  //
  //------------------------------------------------
  RWMol *MolBlockToMol(const std::string &molBlock, bool sanitize, bool removeHs){
    std::istringstream inStream(molBlock);
    RWMol *res=NULL;
    unsigned int line = 0;
    return MolDataStreamToMol(inStream, line, sanitize, removeHs);
  }    


  //------------------------------------------------
  //
  //  Read a molecule from a file
  //
  //------------------------------------------------
  RWMol *MolFileToMol(std::string fName, bool sanitize, bool removeHs){
    std::ifstream inStream(fName.c_str());
    if(!inStream){
      return NULL;
    }
    RWMol *res=NULL;
    if(!inStream.eof()){
      unsigned int line = 0;
      res=MolDataStreamToMol(inStream, line, sanitize, removeHs);
    }
    return res;
  }    
}

