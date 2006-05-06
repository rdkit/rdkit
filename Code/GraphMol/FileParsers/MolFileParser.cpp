// $Id: MolFileParser.cpp 4964 2006-02-18 00:22:34Z glandrum $
//
//  Copyright (C) 2002-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#include "FileParsers.h"
#include "MolFileStereochem.h"
#include <GraphMol/RDKitQueries.h>
#include <RDGeneral/StreamOps.h>

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
  T stripSpacesAndCast(const std::string &input){
    std::string trimmed=boost::trim_copy(input);
    T res;
    return boost::lexical_cast<T>(trimmed);
    try {
      res = boost::lexical_cast<T>(trimmed);
    }
    catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert " << input << " to " << typeid(T).name();
      
      throw FileParseException(errout.str()) ;
    }
    return res;
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
    
    // FIX: couldn't we support NOT lists using the setNegation on the query?
    CHECK_INVARIANT(text[4] != 'T',"[NOT] lists not currently supported");
    
    
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
  
  void ParseChargeLine(RWMol *mol, std::string text) {
    PRECONDITION(text.substr(0,6)==std::string("M  CHG"),"bad atom list line");
    
    // if this line is specified all the atom other than those specified
    // here should carry a charge of 0
    ROMol::AtomIterator ai; 
    for (ai = mol->beginAtoms(); ai != mol->endAtoms(); ai++) {
      (*ai)->setFormalCharge(0);
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
    
    CHECK_INVARIANT(text[14] != 'T',"[NOT] lists not currently supported");
    
    
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
  
  
  void ParseAtomAlias(RWMol *mol,std::string text,std::string &nextLine){
    PRECONDITION(text.substr(0,2)==std::string("A "),"bad atom list line");
      
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
    //int rxnComponentType,rxnComponentNumber,atomMapNumber,inversionFlag,exactChangeFlag;
    
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
    if(text.size()>=36){
      try {
	massDiff = stripSpacesAndCast<int>(text.substr(34,2));
      }
      catch (boost::bad_lexical_cast &) {
	std::ostringstream errout;
	errout << "Cannot convert " << text.substr(34,2) << " to int";
	throw FileParseException(errout.str()) ;
      }
    }    
    chg=0;
    if(text.size()>=39){
      try {
	chg = stripSpacesAndCast<int>(text.substr(36,3));
      }
      catch (boost::bad_lexical_cast &) {
	std::ostringstream errout;
	errout << "Cannot convert " << text.substr(36,3) << " to int";
	throw FileParseException(errout.str()) ;
      }
    }
    hCount = 0;
    if(text.size()>=45){
      try {
	// FIX: go ahead and at least parse the parity field
	hCount = stripSpacesAndCast<int>(text.substr(42,3));
      }
      catch (boost::bad_lexical_cast &) {
	std::ostringstream errout;
	errout << "Cannot convert " << text.substr(42,3) << " to int";
	throw FileParseException(errout.str()) ;
      }
    }
    if(symb=="L" || symb=="A" || symb=="Q" || symb=="*" || symb=="LP"
       || (symb>="R0" && symb<="R9") ){
      res->setAtomicNum(0);
    } else {
      res->setAtomicNum(PeriodicTable::getTable()->getAtomicNumber(symb));
      res->setMass(PeriodicTable::getTable()->getAtomicWeight(res->getAtomicNum()));
    }
    
    //res->setPos(pX,pY,pZ);
    if(chg!=0) res->setFormalCharge(4-chg);
    // FIX: this does not appear to be correct
    if(hCount==1) res->setNoImplicit(true);
    
    if(massDiff!=0) {
      res->setMass(res->getMass()+massDiff);
    }
    
#if 0
    stereoCare=0;
    if(text.size()>=48){
      try {
	stereoCare = stripSpacesAndCast<int>(text.substr(45,3));
      }
      catch (boost::bad_lexical_cast &) {
	std::ostringstream errout;
	errout << "Cannot convert " << text.substr(45,3) << " to int";
	throw FileParseException(errout.str()) ;
      }
    }
    totValence=0;
    if(text.size()>=51){
      try {
	totValence= stripSpacesAndCast<int>(text.substr(48,3));
      }
      catch (boost::bad_lexical_cast &) {
	std::ostringstream errout;
	errout << "Cannot convert " << text.substr(48,3) << " to int";
	throw FileParseException(errout.str()) ;
      }
    }
    rxnComponentType=0;
    if(text.size()>=57){
      try {
	rxnComponentType= stripSpacesAndCast<int>(text.substr(54,3));
      }
      catch (boost::bad_lexical_cast &) {
	std::ostringstream errout;
	errout << "Cannot convert " << text.substr(54,3) << " to int";
	throw FileParseException(errout.str()) ;
      }
    }
    rxnComponentNumber=0;
    if(text.size()>=60){
      try {
	rxnComponentNumber= stripSpacesAndCast<int>(text.substr(57,3));
      }
      catch (boost::bad_lexical_cast &) {
	std::ostringstream errout;
	errout << "Cannot convert " << text.substr(57,3) << " to int";
	throw FileParseException(errout.str()) ;
      }
    }
    atomMapNumber=0;
    if(text.size()>=63){
      try {
	atomMapNumber = stripSpacesAndCast<int>(text.substr(60,3));
      }
      catch (boost::bad_lexical_cast &) {
	std::ostringstream errout;
	errout << "Cannot convert " << text.substr(60,3) << " to int";
	throw FileParseException(errout.str()) ;
      }
    }
    inversionFlag=0;
    if(text.size()>=66){
      try {
	inversionFlag= stripSpacesAndCast<int>(text.substr(63,3));
      }
      catch (boost::bad_lexical_cast &) {
	std::ostringstream errout;
	errout << "Cannot convert " << text.substr(63,3) << " to int";
	throw FileParseException(errout.str()) ;
      }
    }
    exactChangeFlag=0;
    if(text.size()>=69){
      try {
	exactChangeFlag = stripSpacesAndCast<int>(text.substr(66,3));
      }
      catch (boost::bad_lexical_cast &) {
	std::ostringstream errout;
	errout << "Cannot convert " << text.substr(66,3) << " to int";
	throw FileParseException(errout.str()) ;
      }
    }
    

    // save it for later
    res->setProp("stereoCare",stereoCare);
    res->setProp("totValence",totValence);
    res->setProp("rxnComponentType",rxnComponentType);
    res->setProp("rxnComponentNumber",rxnComponentNumber);
    res->setProp("atomMapNumber",atomMapNumber);
    res->setProp("inversionFlag",inversionFlag);
    res->setProp("exactChangeFlag",exactChangeFlag);
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
	} else if(bType == 6){
	  // single or aromatic
	  q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(makeBondOrderEqualsQuery(Bond::SINGLE)));
	  q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(makeBondOrderEqualsQuery(Bond::AROMATIC)));      
	} else if(bType == 7){
	  // double or aromatic
	  q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(makeBondOrderEqualsQuery(Bond::DOUBLE)));
	  q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(makeBondOrderEqualsQuery(Bond::AROMATIC)));
	}
	res->setQuery(q);
      } 
      break;
      
    }
    res->setBeginAtomIdx(idx1);
    res->setEndAtomIdx(idx2);
    res->setBondType(type);

    if( text.size() >= 12)
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
#if 0
    if( text.size() >= 18 )
      try {
	topology = stripSpacesAndCast<int>(text.substr(15,3));
	res->setProp("topology",topology);
      } catch (boost::bad_lexical_cast) {
	;
      }
    
    if( text.size() >= 21 )
      try {
	reactStatus = stripSpacesAndCast<int>(text.substr(18,3));
	res->setProp("reactStatus",reactStatus);
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
  RWMol *MolDataStreamToMol(std::istream *inStream, unsigned int &line, bool sanitize){
    //char inLine[MOLFILE_MAXLINE];
    std::string tempStr;
    bool fileComplete=false;
    bool chiralityPossible = false;

    // mol name
    line++;
    tempStr = getLine(inStream);
    if(inStream->eof()){
      return NULL;
    }
    
    RWMol *res = new RWMol();
    std::string mname = tempStr;
    res->setProp("_Name", mname);
    
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

    try {
      int i;
      Conformer *conf = new Conformer(nAtoms);
      for(i=0;i<nAtoms;i++){
        line++;
	tempStr = getLine(inStream);
	CHECK_INVARIANT(!inStream->eof(),"premature EOF");
        RDGeom::Point3D pos;
	Atom *atom = ParseMolFileAtomLine(tempStr, pos);
        unsigned int aid = res->addAtom(atom,false,true);
        conf->setAtomPos(aid, pos);
      }
      res->addConformer(conf, true);

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
      
      // older mol files can have an atom list block here
      line++;
      tempStr = getLine(inStream);
      if( tempStr[0] != 'M' && tempStr[0] != 'A'
          && tempStr[0] != 'V' && tempStr[0] != 'G'){
        ParseOldAtomList(res,tempStr);
      }
      
      while(!inStream->eof() && tempStr[0] != 'M'){
        if(tempStr.find("A") == 0){
          line++;
	  std::string nextLine = getLine(inStream);
	  if(tempStr.find("M  END") != 0){
            ParseAtomAlias(res,tempStr,nextLine);
          }
        }
        line++;
	tempStr = getLine(inStream);
      }
      //tempStr = inLine;
      while(!inStream->eof() && tempStr.find("M  END") != 0 && tempStr.find("$$$$") != 0){
        if(tempStr.find("M  ALS") == 0) ParseNewAtomList(res,tempStr);
        if(tempStr.find("M  CHG") == 0) ParseChargeLine(res, tempStr);
        line++;
        tempStr = getLine(inStream);
      }
      if(tempStr.find("M  END")==0){
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
    if(res && chiralityPossible){
      MolOps::cleanUp(*res);
      const Conformer &conf = res->getConformer();
      DetectAtomStereoChemistry(*res, &conf);
    }

    if (res && sanitize) {
      try {
	ROMol *tmp=MolOps::removeHs(*res,false,false);
        // unlike DetectAtomStereoChemistry we call DetectBondStereoChemistry here after
        // sanitization because the rings should have been perceived by now, in order to
        // correctly recognize double bonds that may be cis/trans type
        const Conformer &conf = tmp->getConformer();
        DetectBondStereoChemistry(*tmp, &conf);
	delete res;
	res = static_cast<RWMol *>(tmp);
      }
      catch (MolSanitizeException &se){
        delete res;
        throw se;
      }
    }

    
    return res;
  };
  

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


  //------------------------------------------------
  //
  //  Read a molecule from a file
  //
  //------------------------------------------------
  RWMol *MolFileToMol(std::string fName, bool sanitize){
    std::ifstream inStream(fName.c_str());
    if(!inStream){
      return NULL;
    }
    RWMol *res=NULL;
    if(!inStream.eof()){
      unsigned int line = 0;
      res=MolDataStreamToMol(inStream, line, sanitize);
    }
    return res;
  }    
}

