// $Id$
//
//  Copyright (C) 2003-2008  Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/RDLog.h>

#ifdef RDK_BUILD_MAEPARSER_SUPPORT
#undef RDK_BUILD_MAEPARSER_SUPPORT
#endif
#include "MolWriters.h"
#include "FileParsers.h"

namespace RDKit {

SmilesWriter::SmilesWriter(const std::string &fileName,
                           const std::string &delimiter,
                           const std::string &nameHeader, bool includeHeader,
                           bool isomericSmiles, bool kekuleSmiles) {
  if (fileName != "-") {
    auto *tmpStream = new std::ofstream(fileName.c_str());
    if (!(*tmpStream) || (tmpStream->bad())) {
      delete tmpStream;
      std::ostringstream errout;
      errout << "Bad output file " << fileName;
      throw BadFileException(errout.str());
    }
    dp_ostream = static_cast<std::ostream *>(tmpStream);
    df_owner = true;
  } else {
    dp_ostream = static_cast<std::ostream *>(&std::cout);
    df_owner = false;
  }
  this->init(delimiter, nameHeader, includeHeader, isomericSmiles,
             kekuleSmiles);
}

SmilesWriter::SmilesWriter(std::ostream *outStream, std::string delimiter,
                           std::string nameHeader, bool includeHeader,
                           bool takeOwnership, bool isomericSmiles,
                           bool kekuleSmiles) {
  PRECONDITION(outStream, "null stream");
  if (outStream->bad()) {
    throw FileParseException("Bad output stream.");
  }

  dp_ostream = outStream;
  df_owner = takeOwnership;
  this->init(delimiter, nameHeader, includeHeader, isomericSmiles,
             kekuleSmiles);
}
void SmilesWriter::init(const std::string &delimiter,
                        const std::string &nameHeader, bool includeHeader,
                        bool isomericSmiles, bool kekuleSmiles) {
  d_molid = 0;
  d_delim = delimiter;
  d_nameHeader = nameHeader;
  df_includeHeader = includeHeader;
  df_isomericSmiles = isomericSmiles;
  df_kekuleSmiles = kekuleSmiles;
  // these are set by the setProps function as required
  d_props.clear();
}

void SmilesWriter::setProps(const STR_VECT &propNames) {
  if (d_molid > 0) {
    BOOST_LOG(rdErrorLog)
        << "ERROR: Atleast one molecule has already been written\n";
    BOOST_LOG(rdErrorLog)
        << "ERROR: Cannot set properties now - ignoring setProps\n";
    return;
  }
  d_props = propNames;
}

void SmilesWriter::dumpHeader() const {
  CHECK_INVARIANT(dp_ostream, "no output stream");
  if (df_includeHeader) {
    (*dp_ostream) << "SMILES" << d_delim;
    if (d_nameHeader != "") {
      (*dp_ostream) << d_nameHeader << d_delim;
    }

    if (d_props.size() > 0) {
      auto pi = d_props.begin();
      (*dp_ostream) << (*pi);
      pi++;
      while (pi != d_props.end()) {
        (*dp_ostream) << d_delim << (*pi);
        pi++;
      }
    }
    (*dp_ostream) << "\n";
  }
}

SmilesWriter::~SmilesWriter() {
  // close the writer if it's still open:
  if (dp_ostream != nullptr) {
    close();
  }
}

void SmilesWriter::write(const ROMol &mol, int) {
  CHECK_INVARIANT(dp_ostream, "no output stream");
  if (d_molid <= 0 && df_includeHeader) {
    dumpHeader();
  }

  std::string name = "";
  std::string smi = MolToSmiles(mol, df_isomericSmiles, df_kekuleSmiles);
  (*dp_ostream) << smi;
  if (d_nameHeader != "") {
    if (!mol.getPropIfPresent(common_properties::_Name, name) ||
        name.size() == 0) {
      std::stringstream tstream;
      tstream << d_molid;
      name = tstream.str();
    }

    (*dp_ostream) << d_delim << name;
  }

  STR_VECT_CI pi;
  for (pi = d_props.begin(); pi != d_props.end(); pi++) {
    std::string pval;
    // FIX: we will assume that any property that the user requests is castable
    // to
    // a std::string
    if (mol.getPropIfPresent(*pi, pval)) {
      (*dp_ostream) << d_delim << pval;
    } else {
      (*dp_ostream) << d_delim << "";
    }
  }
  (*dp_ostream) << "\n";
  d_molid++;
}
}  // namespace RDKit
