// $Id$
//
//  Copyright (C) 2003-2011 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "FileParsers.h"
#include "MolFileStereochem.h"
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitQueries.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <boost/format.hpp>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>


namespace RDKit{
  
  //*************************************
  //
  // Every effort has been made to adhere to MDL's standard
  // for mol files
  //  
  //*************************************

  const std::string GetMolFileChargeInfo(const RWMol &mol){
    std::stringstream res;
    std::stringstream chgss;
    std::stringstream radss;
    std::stringstream massdiffss;
    unsigned int nChgs=0;
    unsigned int nRads=0;
    unsigned int nMassDiffs=0;
    for(ROMol::ConstAtomIterator atomIt=mol.beginAtoms();
        atomIt!=mol.endAtoms();++atomIt){
      const Atom *atom=*atomIt;
      if(atom->getFormalCharge()!=0){
        ++nChgs;
        chgss << boost::format(" %3d %3d") % (atom->getIdx()+1) % atom->getFormalCharge();
        if(nChgs==8){
          res << boost::format("M  CHG%3d")%nChgs << chgss.str()<<std::endl;
          chgss.str("");
          nChgs=0;
        }
      }
      unsigned int nRadEs=atom->getNumRadicalElectrons();
      if(nRadEs!=0 && atom->getTotalDegree()!=0){
        ++nRads;
        if(nRadEs%2){
          nRadEs=2;
        } else {
          nRadEs=3; // we use triplets, not singlets:
        }
        radss << boost::format(" %3d %3d") % (atom->getIdx()+1) % nRadEs;
        if(nRads==8){
          res << boost::format("M  RAD%3d")%nRads << radss.str()<<std::endl;
          radss.str("");
          nRads=0;
        }
      }
      if(!atom->hasQuery()){
        int isotope=atom->getIsotope();
        if(isotope!=0){
          ++nMassDiffs;
          massdiffss << boost::format(" %3d %3d") % (atom->getIdx()+1) % isotope;
          if(nMassDiffs==8){
            res << boost::format("M  ISO%3d")%nMassDiffs << massdiffss.str()<<std::endl;
            massdiffss.str("");
            nMassDiffs=0;
          }
        }
      }
    }
    if(nChgs){
      res << boost::format("M  CHG%3d")%nChgs << chgss.str()<<std::endl;
    }
    if(nRads){
      res << boost::format("M  RAD%3d")%nRads << radss.str()<<std::endl;
    }
    if(nMassDiffs){
      res << boost::format("M  ISO%3d")%nMassDiffs << massdiffss.str()<<std::endl;
    }
    return res.str();
  }

  bool hasComplexQuery(const Atom *atom){
    PRECONDITION(atom,"bad atom");
    bool res=false;
    if(atom->hasQuery()){
      res=true;
      // counter examples:
      //  1) atomic number
      //  2) the smarts parser inserts AtomAnd queries
      //     for "C" or "c":
      //
      std::string descr=atom->getQuery()->getDescription();
      if(descr=="AtomAtomicNum"){
	res=false;
      } else if(descr=="AtomAnd"){
	if((*atom->getQuery()->beginChildren())->getDescription()=="AtomAtomicNum"){
	  res=false;
	}
      }
    }
    return res;
  }

  bool isListQuery(const Atom::QUERYATOM_QUERY *q){
    // list queries are series of nested ors of AtomAtomicNum queries
    PRECONDITION(q,"bad query");
    bool res=false;
    std::string descr=q->getDescription();
    if(descr=="AtomOr"){
      res=true;
      for(Atom::QUERYATOM_QUERY::CHILD_VECT_CI cIt=q->beginChildren();
          cIt!=q->endChildren() && res;++cIt){
        std::string descr=(*cIt)->getDescription();
        // we don't allow negation of any children of the query:
        if((*cIt)->getNegation()){
          res=false;
        } else if(descr=="AtomOr"){
          res = isListQuery((*cIt).get());
        } else if(descr!="AtomAtomicNum"){
          res=false;
        }
      }
    }
    return res;
  }

  void getListQueryVals(const Atom::QUERYATOM_QUERY *q,INT_VECT &vals){
    // list queries are series of nested ors of AtomAtomicNum queries
    PRECONDITION(q,"bad query");
    std::string descr=q->getDescription();
    PRECONDITION(descr=="AtomOr","bad query");
    if(descr=="AtomOr"){
      for(Atom::QUERYATOM_QUERY::CHILD_VECT_CI cIt=q->beginChildren();
          cIt!=q->endChildren();++cIt){
        std::string descr=(*cIt)->getDescription();
        CHECK_INVARIANT((descr=="AtomOr"||descr=="AtomAtomicNum"),"bad query");
        // we don't allow negation of any children of the query:
        if(descr=="AtomOr"){
          getListQueryVals((*cIt).get(),vals);
        } else if(descr=="AtomAtomicNum"){
          vals.push_back(static_cast<ATOM_EQUALS_QUERY *>((*cIt).get())->getVal());
        }
      }
    }
  }

  bool hasListQuery(const Atom *atom){
    PRECONDITION(atom,"bad atom");
    bool res=false;
    if(atom->hasQuery()){
      res=isListQuery(atom->getQuery());
    }
    return res;
  }

  const std::string GetMolFileQueryInfo(const RWMol &mol){
    std::stringstream ss;
    boost::dynamic_bitset<> listQs(mol.getNumAtoms());
    for(ROMol::ConstAtomIterator atomIt=mol.beginAtoms();
	atomIt!=mol.endAtoms();++atomIt){
      if(hasListQuery(*atomIt)) listQs.set((*atomIt)->getIdx());
    }    
    for(ROMol::ConstAtomIterator atomIt=mol.beginAtoms();
	atomIt!=mol.endAtoms();++atomIt){
      if(!listQs[(*atomIt)->getIdx()] && hasComplexQuery(*atomIt)){
	std::string sma=SmartsWrite::GetAtomSmarts(static_cast<const QueryAtom *>(*atomIt));
	ss<< "V  "<<std::setw(3)<<(*atomIt)->getIdx()+1<<" "<<sma<<std::endl;
      }
    }

    for(ROMol::ConstAtomIterator atomIt=mol.beginAtoms();
	atomIt!=mol.endAtoms();++atomIt){
      if(listQs[(*atomIt)->getIdx()]){
        INT_VECT vals;
        getListQueryVals((*atomIt)->getQuery(),vals);
        ss<<"M  ALS "<<std::setw(3)<<(*atomIt)->getIdx()+1<<" ";
        ss<<std::setw(2)<<vals.size();
        if((*atomIt)->getQuery()->getNegation()){
          ss<<" T";
        } else {
          ss<<" F";
        }
        BOOST_FOREACH(int val,vals){
          ss<<" "<<std::setw(3)<<std::left<<(PeriodicTable::getTable()->getElementSymbol(val));
        }
        ss<<"\n";
      }

    }
    return ss.str();
  }

  const std::string GetMolFileRGroupInfo(const RWMol &mol){
    std::stringstream ss;
    unsigned int nEntries=0;
    for(ROMol::ConstAtomIterator atomIt=mol.beginAtoms();
	atomIt!=mol.endAtoms();++atomIt){
      if((*atomIt)->hasProp("_MolFileRLabel")){
        unsigned int lbl;
        (*atomIt)->getProp("_MolFileRLabel",lbl);
        ss<<" "<<std::setw(3)<<(*atomIt)->getIdx()+1<<" "<<std::setw(3)<<lbl;
        ++nEntries;
      }
    }
    std::stringstream ss2;
    if(nEntries) ss2<<"M  RGP"<<std::setw(3)<<nEntries<<ss.str()<<std::endl;
    return ss2.str();
  }


  const std::string GetMolFileAliasInfo(const RWMol &mol){
    std::stringstream ss;
    for(ROMol::ConstAtomIterator atomIt=mol.beginAtoms();
	atomIt!=mol.endAtoms();++atomIt){
      if((*atomIt)->hasProp("molFileAlias")){
        std::string lbl;
        (*atomIt)->getProp("molFileAlias",lbl);
        if (!lbl.empty())
          ss<<"A  "<<std::setw(3)<<(*atomIt)->getIdx()+1<<"\n"<<lbl<<"\n";
      }
    }
    return ss.str();
  }

  
  const std::string AtomGetMolFileSymbol(const Atom *atom){
    PRECONDITION(atom,"");

    std::string res;
    if(atom->hasProp("_MolFileRLabel")){
      res="R#";
    } else if(hasComplexQuery(atom)){
      if(hasListQuery(atom)){
        res="L";
      } else {
        res="*";
      }
    } else if(atom->getAtomicNum()){
      res=atom->getSymbol();
    } else {
      if(!atom->hasProp("dummyLabel")){
      	res = "R";
      } else {
        std::string symb;
        atom->getProp("dummyLabel",symb);
        if(symb=="*") res="R";
        else if(symb=="X") res="R";
        else if(symb=="Xa") res="R1";
        else if(symb=="Xb") res="R2";
      	else if(symb=="Xc") res="R3";
      	else if(symb=="Xd") res="R4";
      	else if(symb=="Xf") res="R5";
      	else if(symb=="Xg") res="R6";
      	else if(symb=="Xh") res="R7";
      	else if(symb=="Xi") res="R8";
      	else if(symb=="Xj") res="R9";
      	else res=symb;
      }
    }
    // pad the end with spaces
    while(res.size()<3) res += " ";
    return res;
  }

  namespace {
    bool compPair(const std::pair<unsigned int,RDGeom::Point3D> &v1,
                  const std::pair<unsigned int,RDGeom::Point3D> &v2) {
      return (v1.first < v2.first);
    }

    unsigned int getAtomParityFlag(const Atom *atom, const Conformer *conf){
      PRECONDITION(atom,"bad atom");
      PRECONDITION(conf,"bad conformer");
      if(!conf->is3D() ||
         !(atom->getDegree()>=3 && atom->getTotalDegree()==4)) return 0;

      const ROMol &mol=atom->getOwningMol();
      RDGeom::Point3D pos=conf->getAtomPos(atom->getIdx());
      std::vector< std::pair<unsigned int,RDGeom::Point3D> > vs;
      ROMol::ADJ_ITER nbrIdx,endNbrs;
      boost::tie(nbrIdx,endNbrs) = mol.getAtomNeighbors(atom);
      while(nbrIdx!=endNbrs){
        const Atom *at=mol.getAtomWithIdx(*nbrIdx);
        unsigned int idx=at->getIdx();
        RDGeom::Point3D v = conf->getAtomPos(idx);
        v -= pos;
        if(at->getAtomicNum()==1){
          idx += mol.getNumAtoms();
        }
        vs.push_back(std::make_pair(idx,v));
        ++nbrIdx;
      }
      std::sort(vs.begin(),vs.end(),compPair);
      double vol;
      if(vs.size()==4) {
        vol = vs[0].second.crossProduct(vs[1].second).dotProduct(vs[3].second);
      }  else {
        vol = -vs[0].second.crossProduct(vs[1].second).dotProduct(vs[2].second);
      }
      if(vol<0){
        return 2;
      } else if(vol>0) {
        return 1;
      } 
      return 0;
    }
  }
  const std::string GetMolFileAtomLine(const Atom *atom, const Conformer *conf=0){
    PRECONDITION(atom,"");
    std::string res;
    int massDiff,chg,stereoCare,hCount,totValence,rxnComponentType;
    int rxnComponentNumber,atomMapNumber,inversionFlag,exactChangeFlag;
    massDiff=0;
    chg=0;
    stereoCare=0;
    hCount=0;
    totValence=0;
    rxnComponentType=0;
    rxnComponentNumber=0;
    atomMapNumber=0;
    inversionFlag=0;
    exactChangeFlag=0;

    if(atom->hasProp("molAtomMapNumber")){
      atom->getProp("molAtomMapNumber",atomMapNumber);
    }
    
    unsigned int parityFlag=0;
    double x, y, z;
    x = y = z = 0.0;
    if (conf) {
      const RDGeom::Point3D pos = conf->getAtomPos(atom->getIdx());
      x = pos.x; y = pos.y; z = pos.z;
      if(conf->is3D() &&
         atom->getChiralTag()!=Atom::CHI_UNSPECIFIED &&
         atom->getChiralTag()!=Atom::CHI_OTHER
         && atom->getDegree()>=3 &&
         atom->getTotalDegree()==4 ){
        parityFlag=getAtomParityFlag(atom,conf);
      }
    } 
    if (atom->getNumRadicalElectrons()!=0){
      if(atom->getTotalDegree()==0){
        // Specify zero valence for elements/metals without neighbors
        // or hydrogens (degree 0) instead of writing them as radicals.
        totValence = 15;
      } else {
        // write the total valence for other radicals so that we have a chance of
        // reconstructing what was there.
        totValence = atom->getTotalValence();
      }
    }
    std::string symbol = AtomGetMolFileSymbol(atom);
    std::stringstream ss;
    ss << boost::format("%10.4f%10.4f%10.4f %3s%2d%3d%3d%3d%3d%3d  0%3d%3d%3d%3d%3d") % x % y % z % symbol.c_str() %
      massDiff%chg%parityFlag%hCount%stereoCare%totValence%rxnComponentType%
      rxnComponentNumber%atomMapNumber%inversionFlag%exactChangeFlag;
    res += ss.str();
    return res;
  };
  
  const std::string BondGetMolFileSymbol(const Bond *bond){
    PRECONDITION(bond,"");
    // FIX: should eventually recognize queries
    std::string res;
    switch(bond->getBondType()){
    case Bond::SINGLE:
      if(bond->getIsAromatic()){
        res="  4";
      } else {
        res="  1";
      }
      break;
    case Bond::DOUBLE: 
      if(bond->getIsAromatic()){
        res="  4";
      } else {
        res="  2";
      }
      break;
    case Bond::TRIPLE: res="  3";break;
    case Bond::AROMATIC: res="  4";break;
    default: res="  0";break;
    }
    return res;
    //return res.c_str();
  }

  // only valid for single bonds
  int BondGetDirCode(const Bond::BondDir dir){
    int res=0;
    switch(dir){
    case Bond::NONE: res=0;break;
    case Bond::BEGINWEDGE: res=1;break;
    case Bond::BEGINDASH: res=6;break;
    case Bond::UNKNOWN: res=4;break;
    default:
      break;
    }
    return res;
  }

  const std::string GetMolFileBondLine(const Bond *bond, const INT_MAP_INT &wedgeBonds,
                                 const Conformer *conf){
    PRECONDITION(bond,"");
    std::string symbol = BondGetMolFileSymbol(bond);
    int dirCode=0;
    
    Bond::BondDir dir=Bond::NONE;
    bool reverse = false;
    if(bond->getBondType()==Bond::SINGLE){
      // single bond stereo chemistry
      dir = DetermineBondWedgeState(bond, wedgeBonds, conf);
      dirCode = BondGetDirCode(dir);
      // if this bond needs to be wedged it is possible that this
      // wedging was determined by a chiral atom at the end of the
      // bond (instead of at the beginning). In this case we need to
      // reverse the begin and end atoms for the bond when we write
      // the mol file
      if ((dirCode == 1) || (dirCode == 6)) {
        INT_MAP_INT_CI wbi = wedgeBonds.find(bond->getIdx());
        if (static_cast<unsigned int>(wbi->second) != bond->getBeginAtomIdx()) {
          reverse = true;
        }
      }
    } else if (bond->getBondType()==Bond::DOUBLE) {
      // double bond stereochemistry - 
      // if the bond isn't specified, then it should go in the mol block
      // as "any", this was sf.net issue 2963522.
      // two caveats to this:
      // 1) if it's a ring bond, we'll only put the "any"
      //    in the mol block if the user specifically asked for it. 
      //    Constantly seeing crossed bonds in rings, though maybe 
      //    technically correct, is irritating.
      // 2) if it's a terminal bond (where there's no chance of
      //    stereochemistry anyway), we also skip the any.
      //    this was sf.net issue 3009756
      if (bond->getStereo() <= Bond::STEREOANY){
        if(bond->getStereo()==Bond::STEREOANY){
	  dirCode = 3;
	} else if(!(bond->getOwningMol().getRingInfo()->numBondRings(bond->getIdx())) &&
		  bond->getBeginAtom()->getDegree()>1 && bond->getEndAtom()->getDegree()>1){
          dirCode = 3;
        }
      }
    }

    std::stringstream ss;
    if (reverse) {
      // switch the begin and end atoms on the bond line
      ss << std::setw(3) << bond->getEndAtomIdx()+1;
      ss << std::setw(3) << bond->getBeginAtomIdx()+1;
    } else {
      ss << std::setw(3) << bond->getBeginAtomIdx()+1;
      ss << std::setw(3) << bond->getEndAtomIdx()+1;
    }
    ss << symbol;
    ss << " " << std::setw(2) << dirCode;
    return ss.str();
  }
    
  //------------------------------------------------
  //
  //  gets a mol block as a string
  //
  //------------------------------------------------
  std::string MolToMolBlock(const ROMol &mol,bool includeStereo, int confId, bool kekulize){
    ROMol tromol(mol);
    RWMol &trwmol = static_cast<RWMol &>(tromol);
    // NOTE: kekulize the molecule before writing it out
    // because of the way mol files handle aromaticity
    if(kekulize) MolOps::Kekulize(trwmol);

#if 0
    if(includeStereo){
      // assign "any" status to any stereo bonds that are not 
      // marked with "E" or "Z" code - these bonds need to be explictly written
      // out to the mol file
      MolOps::findPotentialStereoBonds(trwmol);
      // now assign stereo code if any have been specified by the directions on
      // single bonds
      MolOps::assignStereochemistry(trwmol);
    }
#endif
    const RWMol &tmol = const_cast<RWMol &>(trwmol);

    std::string res;

    unsigned int nAtoms,nBonds,nLists,chiralFlag,nsText,nRxnComponents;
    unsigned int nReactants,nProducts,nIntermediates;
    nAtoms = tmol.getNumAtoms();
    nBonds = tmol.getNumBonds();
    nLists = 0;
    chiralFlag = 0;
    nsText=0;
    nRxnComponents=0;
    nReactants=0;
    nProducts=0;
    nIntermediates=0;

    if(mol.hasProp("_MolFileChiralFlag")){
      mol.getProp("_MolFileChiralFlag",chiralFlag);
    }
    
    const Conformer *conf;
    if(confId<0 && tmol.getNumConformers()==0){
      conf=0;
    } else {
      conf = &(tmol.getConformer(confId));
    }

    if(tmol.hasProp("_Name")){
      std::string name;
      tmol.getProp("_Name",name);
      res += name;
    }
    res += "\n";

    // info
    if(tmol.hasProp("MolFileInfo")){
      std::string info;
      tmol.getProp("MolFileInfo",info);
      res += info;
    } else {
      std::stringstream ss;
      ss<<"  "<<std::setw(8)<<"RDKit";
      ss<<std::setw(10)<<"";
      if(conf){
        if(conf->is3D()){
          ss<<"3D";
        } else {
          ss<<"2D";
        }
      }
      res += ss.str();
    }
    res += "\n";
    // comments
    if(tmol.hasProp("MolFileComments")){
      std::string info;
      tmol.getProp("MolFileComments",info);
      res += info;
    }
    res += "\n";

    std::stringstream ss;
    ss<<std::setw(3)<<nAtoms;
    ss<<std::setw(3)<<nBonds;
    ss<<std::setw(3)<<nLists;
    ss<<std::setw(3)<<0;
    ss<<std::setw(3)<<chiralFlag;
    ss<<std::setw(3)<<nsText;
    ss<<std::setw(3)<<nRxnComponents;
    ss<<std::setw(3)<<nReactants;
    ss<<std::setw(3)<<nProducts;
    ss<<std::setw(3)<<nIntermediates;
    ss<<"999 V2000\n";
    res += ss.str();
    for(ROMol::ConstAtomIterator atomIt=tmol.beginAtoms();
	atomIt!=tmol.endAtoms();++atomIt){
      res += GetMolFileAtomLine(*atomIt, conf);
      res += "\n";
    }

    INT_MAP_INT wedgeBonds = pickBondsToWedge(tmol);
    for(ROMol::ConstBondIterator bondIt=tmol.beginBonds();
	bondIt!=tmol.endBonds();++bondIt){
      res += GetMolFileBondLine(*bondIt, wedgeBonds, conf);
      res += "\n";
    }

    res += GetMolFileChargeInfo(tmol);
    res += GetMolFileRGroupInfo(tmol);
    res += GetMolFileQueryInfo(tmol);
    res += GetMolFileAliasInfo(tmol);

    // FIX: R-group logic, SGroups and 3D features etc.
    res += "M  END\n";
    return res;
  }
  
  //------------------------------------------------
  //
  //  Dump a molecule to a file
  //
  //------------------------------------------------
  void MolToMolFile(const ROMol &mol,std::string fName,bool includeStereo, int confId, bool kekulize){
    std::ofstream *outStream = new std::ofstream(fName.c_str());
    if (!outStream || !(*outStream) || outStream->bad() ) {
      std::ostringstream errout;
      errout << "Bad output file " << fName;
      throw BadFileException(errout.str());
    }
    std::string outString = MolToMolBlock(mol,includeStereo,confId,kekulize);
    *outStream  << outString;
    delete outStream;
  }    
}

