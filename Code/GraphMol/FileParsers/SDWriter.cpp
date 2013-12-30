// $Id$
//
//  Copyright (C) 2003-2010 Greg Landrum and Rational Discovery LLC
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

#include "MolWriters.h"
#include "FileParsers.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <boost/any.hpp>

namespace RDKit {
  SDWriter::SDWriter(std::string fileName) {
    if(fileName!= "-"){
      std::ofstream *tmpStream = new std::ofstream(fileName.c_str());
      df_owner=true;
      if ( !tmpStream || !(*tmpStream) || (tmpStream->bad()) ) {
        std::ostringstream errout;
        errout << "Bad output file " << fileName;
        throw BadFileException(errout.str());
      }
      dp_ostream = static_cast<std::ostream *>(tmpStream);
    } else {
      dp_ostream = static_cast<std::ostream *>(&std::cout);
      df_owner=false;
    }
    d_molid = 0;
    df_kekulize=true;
    df_forceV3000=false;
  }

  SDWriter::SDWriter(std::ostream *outStream,bool takeOwnership) {
    PRECONDITION(outStream,"null stream");
    if (outStream->bad()){
      throw FileParseException("Bad output stream");
    }
    dp_ostream = outStream;
    df_owner = takeOwnership;
    d_molid = 0;
    df_kekulize=true;
    df_forceV3000=false;
  }

  SDWriter::~SDWriter() {
    // close the writer if it's still open:
    if(dp_ostream!=NULL) close();
  }

  void SDWriter::setProps(const STR_VECT &propNames) {
    if (d_molid > 0) {
      BOOST_LOG(rdWarningLog) << "WARNING: Setting property list after a few molecules have been written\n";
    }
    
    d_props = propNames;
  }

  void SDWriter::write(const ROMol &mol, int confId) {
    PRECONDITION(dp_ostream,"no output stream");

    // write the molecule 
    (*dp_ostream) << MolToMolBlock(mol, true, confId, df_kekulize, df_forceV3000);

    // now write the properties
    STR_VECT_CI pi;
    if (d_props.size() > 0) {
      // check if we have any properties the user specified to write out
      // in which loop over them and write them out
      for (pi = d_props.begin(); pi != d_props.end(); pi++) {
	if (mol.hasProp(*pi)) {
	  writeProperty(mol, (*pi));
	}
      }
    }
    else {
      // if use did not specify any properties, write all non computed properties
      // out to the file
      STR_VECT properties = mol.getPropList();
      STR_VECT compLst;
      if (mol.hasProp(detail::computedPropName)) {
	mol.getProp(detail::computedPropName, compLst);
      }
      STR_VECT_CI pi;
      for (pi = properties.begin(); pi != properties.end(); pi++) {

	// ignore any of the following properties
	if ( ((*pi) == detail::computedPropName) || 
	     ((*pi) == "_Name") ||
	     ((*pi) == "_MolFileInfo") ||
	     ((*pi) == "_MolFileComments") ||
             ((*pi) == "_MolFileChiralFlag")) {
	  continue;
	}

	// check if this property is not computed
	if (std::find(compLst.begin(), compLst.end(), (*pi)) == compLst.end()) {
	  writeProperty(mol, (*pi));
	}
      }
    }
    // add the $$$$ that marks the end of a molecule
    (*dp_ostream) << "$$$$\n";

    ++d_molid;
  }

  void SDWriter::writeProperty(const ROMol &mol, std::string name) {
    PRECONDITION(dp_ostream,"no output stream");

    // write the property value
    // FIX: we will assume for now that the desired property value is 
    // catable to a string 
    std::string pval;
    try {
      mol.getProp(name, pval);
    } catch (boost::bad_any_cast &){
      return;
    }

    // write the property header line
    (*dp_ostream) << ">  <" << name << ">  " << "(" << d_molid+1 << ") " << "\n";
    
    (*dp_ostream) << pval << "\n";
    
    // empty line after the property
    (*dp_ostream) << "\n";
  }
    
}

      
