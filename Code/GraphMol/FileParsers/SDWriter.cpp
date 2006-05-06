// $Id: SDWriter.cpp 4964 2006-02-18 00:22:34Z glandrum $
//
//  Copyright (C) 2003-2006  Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/RDLog.h>

#include "MolWriters.h"
#include "FileParsers.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace RDKit {
  SDWriter::SDWriter(std::string fileName) {
    std::ofstream *tmpStream = new std::ofstream(fileName.c_str());
    
    if ((!(*tmpStream)) || (tmpStream->bad()) ) {
      std::ostringstream errout;
      errout << "Bad output file " << fileName;
      throw FileParseException(errout.str());
    }

    dp_ostream = static_cast<std::ostream *>(tmpStream);

    d_owner = true;
    d_molid = 0;
  }

  SDWriter::~SDWriter() {
    // if we've written any mols, finish with a "$$$$" line
    if (d_molid > 0) {
      (*dp_ostream) << "$$$$\n";
    }

    if (d_owner) {
      delete dp_ostream;
    }
  }

  void SDWriter::setProps(const STR_VECT &propNames) {
    if (d_molid > 0) {
      BOOST_LOG(rdWarningLog) << "WARNING: Setting property list after a few molecules have been written\n";
    }
    
    d_props = propNames;
  }

  void SDWriter::write(ROMol &mol, int confId) {
    //start by writing a "$$$$" line unless this is the first line
    if (d_molid > 0) {
      (*dp_ostream) << "$$$$\n";
    }

    // write the molecule 
    (*dp_ostream) << MolToMolBlock(mol, true, confId);

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
      if (mol.hasProp("computedProps")) {
	mol.getProp("computedProps", compLst);
      }
      STR_VECT_CI pi;
      for (pi = properties.begin(); pi != properties.end(); pi++) {

	// ignore any of the following properties
	if ( ((*pi) == "computedProps") || 
	     ((*pi) == "_Name") ||
	     ((*pi) == "_MolFileInfo") ||
	     ((*pi) == "_MolFileComments")) {
	  continue;
	}

	// check if this property is not computed
	if (std::find(compLst.begin(), compLst.end(), (*pi)) == compLst.end()) {
	  writeProperty(mol, (*pi));
	}
      }
    }

    d_molid++;
  }

  void SDWriter::writeProperty(const ROMol &mol, std::string name) {
    // write the property header line
    (*dp_ostream) << ">  <" << name << ">  " << "(" << d_molid+1 << ") " << "\n";
    
    // write the property value
    // FIX: we will assume for now that the desired property value is 
    // catable to a string 
    std::string pval;
    mol.getProp(name, pval);
    (*dp_ostream) << pval << "\n";
    
    // empty line after the property
    (*dp_ostream) << "\n";
  }
    
}

      
