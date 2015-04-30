// $Id$
//
//  Copyright (C) 2006-2008 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "MolCatalogEntry.h"

#include <RDGeneral/types.h>
#include <RDGeneral/StreamOps.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <iostream>
#include <sstream>
#include <boost/cstdint.hpp>


namespace RDKit {
  
  MolCatalogEntry::MolCatalogEntry(const ROMol *omol){
    PRECONDITION(omol,"bad mol");
    setBitId(-1); 
    dp_props = new Dict();
    d_descrip="";
    dp_mol = omol;
  }

  MolCatalogEntry::MolCatalogEntry(const MolCatalogEntry &other){
    setBitId(other.getBitId()); 
    d_descrip=other.d_descrip;
    dp_props=0;
    dp_mol=0;
    if(other.dp_props){
      dp_props = new Dict(*other.dp_props);
    }
    if(other.dp_mol){
      dp_mol = new ROMol(*other.dp_mol);
    }
  }

  MolCatalogEntry::~MolCatalogEntry() { 
    //std::cerr << "mce: " << dp_mol <<" " <<dp_props << std::endl;
    delete dp_mol;
    dp_mol=NULL;

    delete dp_props;
    dp_props=NULL;

    //std::cerr << "<< done" << std::endl;
  }

  void MolCatalogEntry::setMol(const ROMol *omol){
    PRECONDITION(omol,"bad mol");
    delete dp_mol;
    dp_mol = omol;
  }

  void MolCatalogEntry::toStream(std::ostream &ss) const {
    PRECONDITION(dp_mol,"bad mol");
    MolPickler::pickleMol(*dp_mol,ss);

    boost::int32_t tmpInt;
    tmpInt = getBitId();
    streamWrite(ss,tmpInt);
      
    tmpInt = getOrder();
    streamWrite(ss,tmpInt);
      
    tmpInt = d_descrip.size();
    streamWrite(ss,tmpInt);
    ss.write(d_descrip.c_str(),tmpInt*sizeof(char));

  }    

  std::string MolCatalogEntry::Serialize() const {
    std::stringstream ss(std::ios_base::binary|std::ios_base::out|std::ios_base::in);
    toStream(ss);
    return ss.str();
  }

  void MolCatalogEntry::initFromStream(std::istream &ss){
    delete dp_mol;
    dp_mol=NULL;
    delete dp_props;
    dp_props=NULL;

    // the molecule:
    dp_mol = new ROMol();
    MolPickler::molFromPickle(ss,*const_cast<ROMol *>(dp_mol));

    dp_props = new Dict();
    
    boost::int32_t tmpInt;
    // the bitId:
    streamRead(ss,tmpInt);
    setBitId(tmpInt);

    // the order:
    streamRead(ss,tmpInt);
    setOrder(tmpInt);

    // the description:
    streamRead(ss,tmpInt);
    char *tmpText=new char[tmpInt+1];
    ss.read(tmpText,tmpInt*sizeof(char));
    tmpText[tmpInt]=0;
    d_descrip = tmpText;
    delete [] tmpText;
  }

  void MolCatalogEntry::initFromString(const std::string &text){
    std::stringstream ss(std::ios_base::binary|std::ios_base::out|std::ios_base::in);
    // initialize the stream:
    ss.write(text.c_str(),text.length());
    // now start reading out values:
    initFromStream(ss);
  }
	
}
