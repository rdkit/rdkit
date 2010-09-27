// $Id$
//
//  Copyright (C) 2003-2009 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "FragCatalogEntry.h"

#include <RDGeneral/types.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/StreamOps.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <boost/cstdint.hpp>
#include <RDGeneral/hash/hash.hpp>

namespace RDKit {
  
  FragCatalogEntry::FragCatalogEntry(const ROMol *omol, 
				     const PATH_TYPE &path, 
				     const MatchVectType &aidToFid) {
    PRECONDITION(omol,"bad mol");
    // start with the assumption that this entry is not participating in
    //  any find of fingerprinting
    d_aToFmap.clear();
    setBitId(-1); 
    INT_MAP_INT aIdxMap; // a map from atom id in omol to the new atoms id in mol
    dp_mol = Subgraphs::pathToSubmol(*omol, path, false, aIdxMap); // Using Subgraphs functionality
    d_order = path.size();
    
    // using aIdxMap initialize the location (and their IDs) of the 
    // functional groups on dp_mol
    for(MatchVectType::const_iterator mvtci=aidToFid.begin();
	mvtci!=aidToFid.end();
	mvtci++){
      int oldAid = mvtci->first;
      if(aIdxMap.find(oldAid)!=aIdxMap.end()){
	int newAid = aIdxMap[oldAid];
	if(d_aToFmap.find(newAid)!=d_aToFmap.end()){
	  d_aToFmap[newAid].push_back(mvtci->second);
	}else {
	  INT_VECT tmpVect;
	  tmpVect.clear();
	  tmpVect.push_back(mvtci->second);
	  d_aToFmap[newAid] = tmpVect;
	}
      }
    }
    dp_props = new Dict();

    d_descrip="";
  }

  FragCatalogEntry::FragCatalogEntry(const std::string &pickle){
    d_aToFmap.clear();
    dp_props = new Dict();
    this->initFromString(pickle);
  }


  void FragCatalogEntry::setDescription(const FragCatParams *params){
    PRECONDITION(params,"");
    INT_INT_VECT_MAP::const_iterator fMapIt;
    for(fMapIt=d_aToFmap.begin();fMapIt!=d_aToFmap.end();fMapIt++){
      int atIdx = fMapIt->first;
      INT_VECT fGroups = fMapIt->second;
      std::string label="",temp;

      INT_VECT::const_iterator fGroupIdx=fGroups.begin();
      const ROMol *fGroup;
      for(unsigned int i=0;i<fGroups.size()-1;i++){
	fGroup=params->getFuncGroup(*fGroupIdx);
	fGroup->getProp("_Name",temp);
	label += "(<" + temp + ">)";
	fGroupIdx++;
      }
      fGroup=params->getFuncGroup(*fGroupIdx);
      fGroup->getProp("_Name",temp);
      label += "<"+temp+">";
      dp_mol->getAtomWithIdx(atIdx)->setProp("_supplementalSmilesLabel",label);
    }
    std::string smi = MolToSmiles(*dp_mol);
    //std::cout << "----" << smi << "----" << std::endl;
    d_descrip = smi;
  };


  bool FragCatalogEntry::match(const FragCatalogEntry *other, double tol) const {
    PRECONDITION(other,"bad fragment to compare");
    //std::cerr << " MATCH: "<<d_order<<" " << other->getOrder()<<std::endl;
    if (d_order != other->getOrder()) {
      return false;
    }
    // now check if both the entries have the same number of functional groups
    const INT_INT_VECT_MAP &oFgpMap = other->getFuncGroupMap();
    //std::cerr << "     "<<oFgpMap.size() <<" " <<d_aToFmap.size()<<std::endl;
    if (oFgpMap.size() != d_aToFmap.size()) {
      return false;
    }

    // now check if the IDs are the same 
    INT_INT_VECT_MAP_CI tfi, ofi;
    for (tfi = d_aToFmap.begin(); tfi != d_aToFmap.end(); tfi++) {
      bool found = false;
      //std::cerr << "     "<< (tfi->second[0]) << ":";
      for (ofi = oFgpMap.begin(); ofi != oFgpMap.end(); ofi++) {
        //std::cerr << " "<< (ofi->second[0]);
	if (tfi->second == ofi->second) {
	  found = true;
	  break;
	}
      }
      //std::cerr<<std::endl;
      if (!found) {
	return false;
      }
    }

    // FIX: if might be better if we just do the balaban first and then 
    // move onto eigen values
    Subgraphs::DiscrimTuple tdiscs, odiscs;
    odiscs = other->getDiscrims();
    
    
    //double x1 = boost::tuples::get<0>(odiscs);
    //std::cout << x1 << "\n";
    tdiscs = this->getDiscrims();
#if 0
    std::cout << "DISCRIMS: " << d_descrip  << " ";
    std::cout << tdiscs.get<0>() << " " << tdiscs.get<1>() << " " << tdiscs.get<2>();
    std::cout << "  -- "<<odiscs.get<0>() << " " << odiscs.get<1>() << " " << odiscs.get<2>();
    std::cout << std::endl;
#endif
    // REVIEW: need an overload of feq that handles tuples in MolOps, or wherever
    // DiscrimTuple is defined
    if (!(feq(boost::tuples::get<0>(tdiscs), boost::tuples::get<0>(odiscs),tol)) ||
	!(feq(boost::tuples::get<1>(tdiscs), boost::tuples::get<1>(odiscs),tol)) || 
	!(feq(boost::tuples::get<2>(tdiscs), boost::tuples::get<2>(odiscs),tol)) ) {
      return false;
    }
    
    // FIX: this may not be enough
    // we may have to do teh actual isomorphism mapping
    return true;
  }

  Subgraphs::DiscrimTuple FragCatalogEntry::getDiscrims() const {
    Subgraphs::DiscrimTuple res;
    if (this->hasProp("Discrims")) {
      this->getProp("Discrims", res);
    }
    else {
      PATH_TYPE path;
      for(unsigned int i=0;i<dp_mol->getNumBonds();++i) path.push_back(i);

      // create invariant additions to reflect the functional groups attached to the atoms
      std::vector<boost::uint32_t> funcGpInvars;
      gboost::hash<INT_VECT > vectHasher;
      for(ROMol::AtomIterator atomIt=dp_mol->beginAtoms();
	  atomIt!=dp_mol->endAtoms();
	  ++atomIt){
	unsigned int aid = (*atomIt)->getIdx();
        boost::uint32_t invar=0;
	INT_INT_VECT_MAP_CI mapPos = d_aToFmap.find(aid);
	if(mapPos!= d_aToFmap.end()){
	  INT_VECT fGroups = mapPos->second;
          std::sort(fGroups.begin(),fGroups.end());
          invar=vectHasher(fGroups);
	}
        funcGpInvars.push_back(invar);
      }
      res = Subgraphs::calcPathDiscriminators(*dp_mol,path,true,&funcGpInvars);
      this->setProp("Discrims", res);
    }

    //std::cout << "DISCRIMS: " << d_descrip  << " ";
    //std::cout << res.get<0>() << " " << res.get<1>() << " " << res.get<2>();
    //std::cout << std::endl;
    return res;
  }


  void FragCatalogEntry::toStream(std::ostream &ss) const {
    MolPickler::pickleMol(*dp_mol,ss);

    boost::int32_t tmpInt;
    tmpInt = getBitId();
    streamWrite(ss,tmpInt);
      
    tmpInt = d_descrip.size();
    streamWrite(ss,tmpInt);
    ss.write(d_descrip.c_str(),tmpInt*sizeof(char));

    tmpInt=d_order;
    streamWrite(ss,tmpInt);

    tmpInt = d_aToFmap.size();
    streamWrite(ss,tmpInt);
    for(INT_INT_VECT_MAP::const_iterator iivmci=d_aToFmap.begin();
	iivmci!=d_aToFmap.end();
	iivmci++){
      tmpInt=iivmci->first;
      streamWrite(ss,tmpInt);
      INT_VECT tmpVect=iivmci->second;
      tmpInt=tmpVect.size();
      streamWrite(ss,tmpInt);
      for(INT_VECT_CI ivci=tmpVect.begin();ivci!=tmpVect.end();ivci++){
        tmpInt=*ivci;
	streamWrite(ss,tmpInt);
      }
    }
  }
    
  std::string FragCatalogEntry::Serialize() const {
    std::stringstream ss(std::ios_base::binary|std::ios_base::out|std::ios_base::in);
    toStream(ss);
    return ss.str();
  }

  void FragCatalogEntry::initFromStream(std::istream &ss){
    // the molecule:
    dp_mol = new ROMol();
    MolPickler::molFromPickle(ss,*dp_mol);

    boost::int32_t tmpInt;
    // the bitId:
    streamRead(ss,tmpInt);
    setBitId(tmpInt);

    // the description:
    streamRead(ss,tmpInt);
    char *tmpText=new char[tmpInt+1];
    ss.read(tmpText,tmpInt*sizeof(char));
    tmpText[tmpInt]=0;
    d_descrip = tmpText;
    delete [] tmpText;

    streamRead(ss,tmpInt);
    d_order=tmpInt;

    // now the map:
    streamRead(ss,tmpInt);
    for(int i=0;i<tmpInt;i++){
      boost::int32_t key,value,size;
      streamRead(ss,key);
      streamRead(ss,size);
      INT_VECT tmpVect;
      tmpVect.clear();
      for(int j=0;j<size;j++){
	streamRead(ss,value);
	tmpVect.push_back(value);
      }
      d_aToFmap[key] = tmpVect;
    }
  }

  void FragCatalogEntry::initFromString(const std::string &text){
    std::stringstream ss(std::ios_base::binary|std::ios_base::out|std::ios_base::in);
    // initialize the stream:
    ss.write(text.c_str(),text.length());
    // now start reading out values:
    initFromStream(ss);
  }
	
}
