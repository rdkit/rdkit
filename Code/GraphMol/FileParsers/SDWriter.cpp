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
SDWriter::SDWriter(const std::string &fileName) {
  if (fileName != "-") {
    auto *tmpStream = new std::ofstream(fileName.c_str());
    df_owner = true;
    if (!(*tmpStream) || (tmpStream->bad())) {
      delete tmpStream;
      std::ostringstream errout;
      errout << "Bad output file " << fileName;
      throw BadFileException(errout.str());
    }
    dp_ostream = static_cast<std::ostream *>(tmpStream);
  } else {
    dp_ostream = static_cast<std::ostream *>(&std::cout);
    df_owner = false;
  }
  d_molid = 0;
  df_kekulize = true;
  df_forceV3000 = false;
}

SDWriter::SDWriter(std::ostream *outStream, bool takeOwnership) {
  PRECONDITION(outStream, "null stream");
  if (outStream->bad()) {
    throw FileParseException("Bad output stream");
  }
  dp_ostream = outStream;
  df_owner = takeOwnership;
  d_molid = 0;
  df_kekulize = true;
  df_forceV3000 = false;
}

SDWriter::~SDWriter() {
  // close the writer if it's still open:
  if (dp_ostream != nullptr) {
    close();
  }
}

void SDWriter::setProps(const STR_VECT &propNames) {
  if (d_molid > 0) {
    BOOST_LOG(rdWarningLog) << "WARNING: Setting property list after a few "
                               "molecules have been written\n";
  }

  d_props = propNames;
}

namespace {
void _writePropToStream(std::ostream *dp_ostream, const ROMol &mol,
                        const std::string &name, int d_molid) {
  PRECONDITION(dp_ostream, "no output stream");

  // write the property value
  // FIX: we will assume for now that the desired property value is
  // catable to a string
  std::string pval;
  try {
    mol.getProp(name, pval);
  } catch (boost::bad_any_cast &) {
    return;
  }

  // warn and skip if we include a new line
  if (name.find("\n") != std::string::npos) {
    BOOST_LOG(rdWarningLog)
        << "WARNING: Skipping property " << name
        << " because the name includes a newline" << std::endl;
    return;
  }
  if (pval.find("\r\n\r\n") != std::string::npos ||
      pval.find("\n\n") != std::string::npos) {
    BOOST_LOG(rdWarningLog)
        << "WARNING: Skipping property " << name
        << " because the value includes an illegal blank line" << std::endl;
    return;
  }
  // write the property header line
  (*dp_ostream) << ">  <" << name << ">  ";
  if (d_molid >= 0) {
    (*dp_ostream) << "(" << d_molid + 1 << ") ";
  }
  (*dp_ostream) << "\n";

  (*dp_ostream) << pval << "\n";

  // empty line after the property
  (*dp_ostream) << "\n";
}
void _MolToSDStream(std::ostream *dp_ostream, const ROMol &mol, int confId,
                    bool df_kekulize, bool df_forceV3000, int d_molid,
                    STR_VECT *props) {
  PRECONDITION(dp_ostream, "no output stream");

  // write the molecule
  (*dp_ostream) << MolToMolBlock(mol, true, confId, df_kekulize, df_forceV3000);

  // now write the properties
  STR_VECT_CI pi;
  if (props && props->size() > 0) {
    // check if we have any properties the user specified to write out
    // in which loop over them and write them out
    for (pi = props->begin(); pi != props->end(); pi++) {
      if (mol.hasProp(*pi)) {
        _writePropToStream(dp_ostream, mol, (*pi), d_molid);
      }
    }
  } else {
    // if use did not specify any properties, write all non computed properties
    // out to the file
    STR_VECT properties = mol.getPropList();
    STR_VECT compLst;
    mol.getPropIfPresent(RDKit::detail::computedPropName, compLst);

    STR_VECT_CI pi;
    for (pi = properties.begin(); pi != properties.end(); pi++) {
      // ignore any of the following properties
      if (((*pi) == RDKit::detail::computedPropName) ||
          ((*pi) == common_properties::_Name) || ((*pi) == "_MolFileInfo") ||
          ((*pi) == "_MolFileComments") ||
          ((*pi) == common_properties::_MolFileChiralFlag)) {
        continue;
      }

      // check if this property is not computed
      if (std::find(compLst.begin(), compLst.end(), (*pi)) == compLst.end()) {
        _writePropToStream(dp_ostream, mol, (*pi), d_molid);
      }
    }
  }
  // add the $$$$ that marks the end of a molecule
  (*dp_ostream) << "$$$$\n";
}
}  // namespace

std::string SDWriter::getText(const ROMol &mol, int confId, bool kekulize,
                              bool forceV3000, int molid, STR_VECT *propNames) {
  std::stringstream sstr;
  _MolToSDStream(&sstr, mol, confId, kekulize, forceV3000, molid, propNames);
  return sstr.str();
};

void SDWriter::write(const ROMol &mol, int confId) {
  PRECONDITION(dp_ostream, "no output stream");
  _MolToSDStream(dp_ostream, mol, confId, df_kekulize, df_forceV3000, d_molid,
                 &d_props);
  ++d_molid;
}

void SDWriter::writeProperty(const ROMol &mol, const std::string &name) {
  PRECONDITION(dp_ostream, "no output stream");

  _writePropToStream(dp_ostream, mol, name, d_molid);
}
}  // namespace RDKit
