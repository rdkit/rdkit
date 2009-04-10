// $Id$
//
//  Copyright (C) 2003-2009 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
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
    dp_mol = Subgraphs::PathToSubmol(*omol, path, false, aIdxMap); // Using Subgraphs functionality
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
      int nAts = dp_mol->getNumAtoms();
      int nBds = dp_mol->getNumBonds();
      double *dMat;
      // Get the distance matrix without touching the diagonal elements
      //  (we'll modify those ourselves below):
      dMat = MolOps::getDistanceMat(*dp_mol, true, false);

      //  Our discriminators (based upon eigenvalues of the distance matrix
      //  and Balaban indices) are good, but is a case they don't properly
      //  handle by default.  These two fragments:
      //    C-ccc and c-ccc
      //  give the same discriminators because their distance matrices are
      //  identical.  We'll work around this by modifying the diagonal
      //  elements of the distance matrix corresponding to aromatic atoms


      // now modify the diagonal element of distance matrix to reflect:
      //   1) the identity of the atom
      //   2) whether or not the atom is aromatic
      //   3) the functional groups attached to the atom
      int aid;
      for(ROMol::AtomIterator atomIt=dp_mol->beginAtoms();
	  atomIt!=dp_mol->endAtoms();
	  atomIt++){
	double diagElem;
	aid = (*atomIt)->getIdx();
	// start with atom identity and aromaticity:
	if((*atomIt)->getIsAromatic()){
	  diagElem = 6.0 / ((*atomIt)->getAtomicNum() + 0.5);
	} else {
	  diagElem = 6.0 / (*atomIt)->getAtomicNum();
	}

	INT_INT_VECT_MAP_CI mapPos = d_aToFmap.find(aid);
	if(mapPos!= d_aToFmap.end()){
	  // figure out the product of the primes for the functional group
	  // indices on this atom:
	  double prod = 1.0;
	  const INT_VECT &fGroups = mapPos->second;
	  for(INT_VECT_CI fGroupIdx=fGroups.begin();
	      fGroupIdx != fGroups.end();
	      fGroupIdx++){
	    prod *= firstThousandPrimes[*fGroupIdx];
	  }
	  // divide that until it's less than 10:
	  while(prod>10.) prod /= 10.;
	  // and add it to the diagonal:
	  diagElem += prod;
	}

	dMat[aid*nAts+aid] = diagElem;
      }
      // now do the discriminators
      res = Subgraphs::computeDiscriminators(dMat, nBds, nAts);

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
