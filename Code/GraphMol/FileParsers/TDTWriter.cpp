// $Id$
//
//  Copyright (C) 2005-2006  Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include "MolWriters.h"
#include "FileParsers.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <boost/algorithm/string.hpp>

namespace RDKit {
  TDTWriter::TDTWriter(std::string fileName) {
    std::ofstream *tmpStream = new std::ofstream(fileName.c_str());
    
    if ((!(*tmpStream)) || (tmpStream->bad()) ) {
      std::ostringstream errout;
      errout << "Bad output file " << fileName;
      throw FileParseException(errout.str());
    }

    dp_ostream = static_cast<std::ostream *>(tmpStream);

    d_owner = true;
    d_molid = 0;
    d_numDigits=4;
    df_write2D=false;
    df_writeNames=true;
  }

  TDTWriter::TDTWriter(std::ostream *outStream, bool takeOwnership) {
    PRECONDITION(outStream,"null stream");
    if (outStream->bad()) {
      throw FileParseException("Bad output stream");
    }
    dp_ostream = outStream;
    d_owner = takeOwnership;
    d_molid = 0;
    d_numDigits=4;
    df_write2D=false;
    df_writeNames=true;
  }

  TDTWriter::~TDTWriter() {
    // if we've written any mols, finish with a "|" line
    if (d_molid > 0) {
      CHECK_INVARIANT(dp_ostream,"null outstream even though molecules were written");
      (*dp_ostream) << "|\n";
    }

    if (d_owner) {
      delete dp_ostream;
    }
  }

  void TDTWriter::setProps(const STR_VECT &propNames) {
    if (d_molid > 0) {
      BOOST_LOG(rdWarningLog) << "WARNING: Setting property list after a few molecules have been written\n";
    }
    
    d_props = propNames;
  }

  void TDTWriter::write(ROMol &mol, int confId) {
    CHECK_INVARIANT(dp_ostream,"no output stream");
    //start by writing a "|" line unless this is the first line
    if (d_molid > 0) {
      (*dp_ostream) << "|\n";
    }

    // write the molecule 
    (*dp_ostream) << "$SMI<" << MolToSmiles(mol) << ">\n";

    if(df_writeNames && mol.hasProp("_Name")){
      std::string name;
      mol.getProp("_Name",name);
      (*dp_ostream) << "NAME<" << name << ">\n";
    }
    
    // do we need to write coordinates?
    if(mol.getNumConformers()){
      const Conformer &conf = mol.getConformer(confId);
      if(df_write2D){
	(*dp_ostream) << "2D<";
      } else {
	(*dp_ostream) << "3D<";
      }
      const RDGeom::POINT3D_VECT &coords=conf.getPositions();
      int nAts=conf.getNumAtoms();
      for(int i=0;i<nAts;i++){
	(*dp_ostream) << std::setprecision(d_numDigits) << coords[i].x << ",";
	(*dp_ostream) << std::setprecision(d_numDigits) << coords[i].y;
	if(!df_write2D){
	  (*dp_ostream) << "," << std::setprecision(d_numDigits) << coords[i].z;
	}
	if(i!=nAts-1) (*dp_ostream) << ",";
      }
      (*dp_ostream) << ";>\n";

    }
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

  void TDTWriter::writeProperty(const ROMol &mol, std::string name) {
    PRECONDITION(dp_ostream,"no output stream");
    (*dp_ostream) << name << "<";
    
    // write the property value
    // FIX: we will assume for now that the desired property value is 
    // catable to a string 
    std::string pval;

    // we need to remove any line breaks in the output, replace them with spaces
    mol.getProp(name, pval);
    boost::replace_all(pval,"\n"," ");
    (*dp_ostream) << pval << ">\n";
  }
    

}

      
